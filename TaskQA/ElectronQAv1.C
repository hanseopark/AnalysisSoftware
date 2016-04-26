/****************************************************************************************************************************
 ******         provided by Gamma Conversion Group, PWG4,                                                                                                               *****
 ******         Ana Marin, marin@physi.uni-heidelberg.de                                                                                                        *****
 ******         Kathrin Koch, kkoch@physi.uni-heidelberg.de                                                                                                     *****
 ******         Friederike Bock, friederike.bock@cern.ch                                                                                                        *****
 *****************************************************************************************************************************/

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
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TMinuit.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TPolyLine.h";
#include "TGraph2D.h"
#include "TMatrixDSym.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "ProduceFinalResults.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "THnSparse.h"
#include "ElectronQAv1.h"
//#include "THnSparse.h"






void ElectronQAv1(TString outputData="",TString outputMC="", TString fCutSelection="8000011002093603007200000000_9047540023310262371_01031035009000",TString Suffix="eps",TString energy="pPb_5.023TeV", TString period="",Int_t mode = 9){

   
  
	
	
	TString outputDir = Form("%s/%s/%s/ElectronQA",fCutSelection.Data(),energy.Data(),Suffix.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	
		
	
	
	StyleSettingsThesis(Suffix);
	SetPlotStyle();
	
	gStyle->SetOptTitle(kFALSE); 
	
	
	TFile fileData(outputData.Data());
	TH1::AddDirectory(kFALSE);
	fileData.Print();
	
	cout<<"Entro pedro"<<endl;
	
	
	
	
	
	TString nameOutputDir = "";
	
	if( mode == 9 || mode == 1 ){
	  
	     nameOutputDir = "GammaConvDalitzV1";
	}
	
	
	else if( mode == 6 || mode == 7 ){
	  
	      nameOutputDir = "GammaCaloDalitz";
	  
	}
	
	
	
	
	TList* GammaConvDalitzV1Data = (TList*)fileData.Get(nameOutputDir.Data());

	        
        
        if( ! GammaConvDalitzV1Data ) {
	  
                    Int_t iTrain = 1;

                    while( !GammaConvDalitzV1Data && iTrain <= 100 ) {
                            
                          GammaConvDalitzV1Data = (TList*)fileData.Get(Form("%s_%d",nameOutputDir.Data(),iTrain));
                           
                          iTrain++;
                    
                    }
                     
                    if( ! GammaConvDalitzV1Data ) {
		      
		        GammaConvDalitzV1Data = (TList*)fileData.Get("GammaConvDalitzV1QA");
			
			if( ! GammaConvDalitzV1Data ) {		
			    cout<<"ERROR Data: GammaConvDalitzV1 is not found in the file"<<endl;
			    return;
			}
                   }            
        }
        


        
	if ( mode == 9  ){
	  
		ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection,kTRUE);
		fEventCutSelection = fGammaCutSelection(0,7);
		fGammaCutSelection = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
		cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
		
		
	} else {
	  
	
		ReturnSeparatedCutNumberAdvanced(fCutSelection, fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
		cout<<"Testing 0: " << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() <<"\t"<<fClusterCutSelection.Data()<<"\t"<< fElectronCutSelection.Data()<<"\t"<<fMesonCutSelection.Data()<<endl;
	}    
	
	
	
	
	fCutSelectionMC         =  fCutSelection;
	fEventCutSelectionMC    =  fEventCutSelection;
	fGammaCutSelectionMC    =  fGammaCutSelection;
	fClusterCutSelectionMC  =  fClusterCutSelection;
	fElectronCutSelectionMC =  fElectronCutSelection;
	fMesonCutSelectionMC    =  fMesonCutSelection;
	
	if( energy.CompareTo("pPb_5.023TeV") == 0 && outputMC.Contains("68") ){
	  
	  fCutSelectionMC         =  "80000113_00200009360300007200004000_20405400233002223710_0263103500900000";
	  fElectronCutSelectionMC =  "20405400233002223710";

	  
	}
	
	


	cout<<"Data\t"<<"EventSelection: "<< fEventCutSelection.Data()   <<", PhotonCut: " << fGammaCutSelection.Data()   <<", ElectronCut: " << fElectronCutSelection.Data()  <<", MesonCut: " << fMesonCutSelection.Data() << endl;
	cout<<"MC\t"  <<"EventSelection: "<< fEventCutSelectionMC.Data() <<", PhotonCut: " << fGammaCutSelectionMC.Data() <<", ElectronCut: " << fElectronCutSelectionMC.Data()<<", MesonCut: "  <<fMesonCutSelectionMC.Data()<<endl;

	
	
	GammaConvDalitzV1Data->Print();

        TList* fHistosGammaConversionDalitzData = (TList*) GammaConvDalitzV1Data->FindObject(Form("Cut Number %s",fCutSelection.Data()) );
	

        if( ! fHistosGammaConversionDalitzData ){
	  
                cout<<Form("Data: Cut Number %s",fCutSelection.Data())<<" is not found in the file"<<endl;
                return;
		
        }
											   //ElectronCuts_30105400000003300000
											   
        TString fElectronCutsPreSelection = "30105400000003300000";							      ;
        TList* fElectronCutsPreData  = (TList*) GammaConvDalitzV1Data->FindObject(Form("ElectronCuts_%s",fElectronCutsPreSelection.Data()));
        
        
	
	if ( !fElectronCutsPreData ) {
	  
	    cout<<"El error esta aqui:     "<<Form("ElectronCuts_%s",fElectronCutsPreSelection.Data())<<endl;
	    //return;
	  
	}
	
	
        TString fConvCutsPreSelection = "30105400000003300000"								      ;
        TList* fConvCutsPreData  = (TList*) GammaConvDalitzV1Data->FindObject(Form("ConvCuts_%s",fConvCutsPreSelection.Data()));
      
        
	/*if ( !fConvCutsPreData ) {
	  
	    cout<<"El error esta aqui:     "<<Form("ElectronCuts_%s",fConvCutsPreSelection.Data())<<endl;
	    return;
	  
	}*/
        
       
        TList *fElectronCutsData   = (TList*)fHistosGammaConversionDalitzData->FindObject(Form("ElectronCuts_%s",fCutSelection.Data() ) );
       
        if(  !fElectronCutsData  ) {
	  
	    fElectronCutsData   = (TList*)fHistosGammaConversionDalitzData->FindObject(Form("ElectronCuts_%s",fElectronCutSelection.Data() ) );
	    
	    if( ! fElectronCutsData ) {
	      cout<<Form("Data: ElectronCuts_%s",fElectronCutSelection.Data())<<" is not found in the file"<<endl;
	      return;
	    }
	  
        }
       
        TList *fConvCutsData   = (TList*)fHistosGammaConversionDalitzData->FindObject(Form("ConvCuts_%s",fGammaCutSelection.Data() ) );
       
        if( ! fConvCutsData ) {
	 
	 cout<<Form("ConvCuts_%s",fGammaCutSelection.Data())<<" is not found in the file"<<endl;
         return;
	  
        }
       
	TList *fESDContainerData   = (TList*)fHistosGammaConversionDalitzData->FindObject(Form("%s ESD histograms",fCutSelection.Data()));
       
	if(  !fESDContainerData  ) {
	  
	     cout<<Form("%s ESD histograms",fCutSelection.Data())<<" is not found in the file"<<endl;
	     return;
	 
       }
       
       
	TList  *fQAFolderData       = (TList*)fHistosGammaConversionDalitzData->FindObject(Form("%s QA histograms",fCutSelection.Data()));
       
       if( !fQAFolderData ) {
	 
	    cout<<Form("%s QA histograms",fCutSelection.Data())<<" is not found in the file"<<endl;
	 
       }
       
       
       
       /////////////////////////////////////////////////////////////////////////////////////////////
       
       
       
        
        fEventQualityData = (TH1F*)fESDContainerData->FindObject("NEvents");
	
        if ( energy.CompareTo("PbPb_2.76TeV") == 0 || energy.CompareTo("pPb_5.023TeV") == 0 ){ 
	
		nEventsData = fEventQualityData->GetBinContent(1);
	   
	 } else {	   
	   
		nEventsData =  GetNEvents(fEventQualityData);
		
	 }
       
   
   
      
	
	fESD_NumberOfGoodESDTracksData = (TH1F*)fESDContainerData->FindObject("GoodESDTracks");
		
      
        mean_Data = fESD_NumberOfGoodESDTracksData->GetMean(); 
	
	NormFactorData = (1./nEventsData)*(1./mean_Data);
	//NormFactorData = 1./mean_Data;
	
	cout<<"Mean Data: "<<mean_Data<<endl;
	cout<<"Normalization Factor: "<<NormFactorData<<endl;
	
	
	
       
       
	
	
       Float_t  rebinEta = 1;
       Float_t  rebinPhi = 1;
       Double_t YMinEta  = 1;
       Double_t YMaxEta  = 1;
       
        if ( energy.CompareTo("pPb_5.023TeV") == 0 ){
	  
	      rebinEta = 4;
	      rebinPhi = 4;
	      YMinEta  = 1e-3;
	      YMaxEta  = 3e-2;
	      
	      
	} else if ( energy.CompareTo("7TeV") == 0 ){
	  
	      rebinEta = 4;
	      rebinPhi = 4;
	      YMinEta  = 1e-3;
	      YMaxEta  = 3e-2;
	      
	} else if ( energy.CompareTo("2.76TeV") == 0 ){
	  
	      rebinEta = 8;
	      rebinPhi = 4;
	      
	      YMinEta  = 1e-3;
	      YMaxEta  = 3e-2;
	}
       
       
       //fElectronCutSelection.
       

       
       
       
       
       Float_t floatLocationLeft[4] = {0.2,0.95,0.035, 0.02};
       
       TString collisionSystem  = ReturnFullCollisionsSystem(energy);
       
       TLatex *processSystem = new TLatex(0.63, 0.92,Form(collisionSystem.Data()));
        processSystem->SetNDC();
        processSystem->SetTextColor(1);
        processSystem->SetTextSize(0.04);
	processSystem->SetTextFont(42);
       
        hESDConvGammaEta_data                            	= (TH1F*)fESDContainerData->FindObject("ESD_ConvGamma_Eta");
	
	if(  fQAFolderData  ) {
	
       	hESDDalitzElectronAfterPt_data 				= (TH1F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_Pt");
	hESDDalitzPositronAfterPt_data  			= (TH1F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_Pt");
	hESDDalitzElectronAfterEta_data 			= (TH1F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_Eta");
	hESDDalitzPositronAfterEta_data 			= (TH1F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_Eta");
	hESDDalitzElectronAfterEtaPCut_data               	= (TH1F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_Eta_PCut");
	hESDDalitzPositronAfterEtaPCut_data               	= (TH1F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_Eta_PCut");
	hESDDalitzElectronAfterPhi_data 			= (TH1F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_Phi");
	hESDDalitzPositronAfterPhi_data 			= (TH1F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_Phi");
        hESDDalitzElectronAfterNClsITS_data			= (TH1F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_NClsITS");
	hESDDalitzPositronAfterNClsITS_data			= (TH1F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_NClsITS");
      
        hESDDalitzElectronAfterNFindClsTPC_data			= (TH2F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_NFindClsTPC");
	hESDDalitzPositronAfterNFindClsTPC_data			= (TH2F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_NFindClsTPC");
	hESDDalitzElectronAfterNFindClsTPCPCut_data		= (TH1F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_NFindClsTPC_PCut");
	hESDDalitzPositronAfterNFindClsTPCPCut_data		= (TH1F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_NFindClsTPC_PCut");


	hESDDalitzElectronAfterNClsTPC_data			= (TH2F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_NClsTPC");
	hESDDalitzPositronAfterNClsTPC_data			= (TH2F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_NClsTPC");
		
	hESDDalitzElectronAfterNClsTPCPCut_data			= (TH1F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_NClsTPC_PCut");
	hESDDalitzPositronAfterNClsTPCPCut_data			= (TH1F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_NClsTPC_PCut");
	
	
	hESDDalitzElectronAfterNCrossedRowsTPC_data		= (TH2F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_NCrossedRowsTPC");
	hESDDalitzPositronAfterNCrossedRowsTPC_data		= (TH2F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_NCrossedRowsTPC");
	
	hESDDalitzElectronAfterNCrossedRowsTPCPCut_data         = (TH1F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_NCrossedRowsTPC_PCut");
	hESDDalitzPositronAfterNCrossedRowsTPCPCut_data         = (TH1F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_NCrossedRowsTPC_PCut");
	
		
	
	hESDDalitzPosEleAfterDCAxy_data				= (TH2F*)fQAFolderData->FindObject("ESD_DalitzPosEle_After_DCAxy");
	hESDDalitzPosEleAfterDCAz_data  			= (TH2F*)fQAFolderData->FindObject("ESD_DalitzPosEle_After_DCAz");
      
	
	hESDDalitzElectronAfterTPCdEdxVsP_data			= (TH2F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_TPCdEdxVsP");
        hESDDalitzPositronAfterTPCdEdxVsP_data			= (TH2F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_TPCdEdxVsP");
	hESDDalitzElectronAfterTPCdEdxSignalVsP_data		= (TH2F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_TPCdEdxSignalVsP");
        hESDDalitzPositronAfterTPCdEdxSignalVsP_data		= (TH2F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_TPCdEdxSignalVsP");
      
	hESDDalitzElectronAfterTPCdEdxVsEta_data		= (TH2F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_TPCdEdxVsEta");
        hESDDalitzPositronAfterTPCdEdxVsEta_data		= (TH2F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_TPCdEdxVsEta");
      
	hESDDalitzElectronAfterTPCdEdxVsPhi_data		= (TH2F*)fQAFolderData->FindObject("ESD_DalitzElectron_After_TPCdEdxVsPhi");
        hESDDalitzPositronAfterTPCdEdxVsPhi_data		= (TH2F*)fQAFolderData->FindObject("ESD_DalitzPositron_After_TPCdEdxVsPhi");
	
        hESDMotherPhi_data					= (TH1F*)fQAFolderData->FindObject("ESD_DalitzMother_Phi");
        hESDEposEnegPsiPairDPhi_data				= (TH2F*)fQAFolderData->FindObject("ESD_EposEneg_PsiPair_DPhi");
        hESDEposEnegInvMassPt_data				= (TH2F*)fQAFolderData->FindObject("ESD_EposEneg_InvMassPt");
        hESDEposEnegLikeSignBackInvMassPt_data			= (TH2F*)fQAFolderData->FindObject("ESD_EposEneg_LikeSignBack_InvMassPt");
	
	sESDConvGammaZR_data                              	= (THnSparseF*)fQAFolderData->FindObject("ESD_ConvGamma_ZR");
	sESDConvGammaXY_data 					= (THnSparseF*)fQAFolderData->FindObject("ESD_ConvGamma_XY");
	sESDConvGammaXY_data 					= (THnSparseF*)fQAFolderData->FindObject("ESD_ConvGamma_XY");
	
	
	
      
	}
								  
	 
	if(  fConvCutsData ) {
      
	  hTPCdEdxSigafter_data 					= (TH2F*)fConvCutsData->FindObject(Form("Gamma_dEdxSig_after %s",fGammaCutSelection.Data()));
	  hTPCdEdxafter_data    					= (TH2F*)fConvCutsData->FindObject(Form("Gamma_dEdx_after %s",   fGammaCutSelection.Data()));
	
	}
	
	if( fElectronCutsPreData ){ //Just for thesis {
	  
	  
           hTPCdEdxSigbefore_data 					= (TH2F*)fElectronCutsPreData->FindObject(Form("Electron_dEdxSignal_before %s",fElectronCutsPreSelection.Data()));
	   hTPCdEdxbefore_data						= (TH2F*)fElectronCutsPreData->FindObject(Form("Electron_dEdx_before %s",fElectronCutsPreSelection.Data()));
	   
	   
	   hTPCdEdxSigbefore_data->Print();
  
	   cout<<"Hola entro a electron preselection "<<endl;
	  
	}
	
	if( fElectronCutsData ) {
	  
	   hTPCEledEdxSigafter_data 					= (TH2F*)fElectronCutsData->FindObject(Form("Electron_dEdxSignal_after %s",fElectronCutSelection.Data() ));
	   
	   //hTPCdEdxafter->Print();
	   
	   hTPCEledEdxafter_data 					= (TH2F*)fElectronCutsData->FindObject(Form("Electron_dEdx_after %s",fElectronCutSelection.Data() ));
  
	   
	   cout<<"Hola entro a electron after cuts "<<endl;
	
	  
	}
	
	
	
	//////////////////////////////////////////MonteCarlo///////////////////////////////////////////////////////////////
	
	
	
	TFile fileMC(outputMC.Data());
	
	TList* GammaConvDalitzV1MC = (TList*)fileMC.Get(nameOutputDir.Data());

	        
        
        if( ! GammaConvDalitzV1MC ) {

                    Int_t iTrain = 1;

                    while( !GammaConvDalitzV1MC && iTrain <= 100 ) {
                            
                          GammaConvDalitzV1MC = (TList*)fileMC.Get(Form("%s_%d",nameOutputDir.Data(),iTrain));
                           
                          iTrain++;
                    
                     }
                    if( ! GammaConvDalitzV1MC ) {
		      
		           GammaConvDalitzV1MC = (TList*)fileMC.Get("GammaConvDalitzV1QA");
			
			if( ! GammaConvDalitzV1MC ) {
			
			    cout<<"ERROR MC: GammaConvDalitzV1 is not found in the file"<<endl;
			    return;
			}

                        
                    }
            
        }

       
        TList* fHistosGammaConversionDalitzMC = (TList*) GammaConvDalitzV1MC->FindObject(Form("Cut Number %s",fCutSelectionMC.Data()) );
	

        if( ! fHistosGammaConversionDalitzMC ){
                cout<<Form("MC: Cut Number %s",fCutSelectionMC.Data())<<" is not found in the file"<<endl;
                return;
        }
        
	fHistosGammaConversionDalitzMC->Print();
      
       
       TList *fConvCutsMC   = (TList*)fHistosGammaConversionDalitzMC->FindObject(Form("ConvCuts_%s",fGammaCutSelectionMC.Data() ) );
       
       if( ! fConvCutsMC ) {
	 
	 cout<<Form("MC: ConvCuts_%s",fGammaCutSelection.Data())<<" is not found in the file"<<endl;
         return;
	  
       }
       
    
        TList* fElectronCutsPreMC  = (TList*) GammaConvDalitzV1MC->FindObject(Form("ElectronCuts_%s",fElectronCutsPreSelection.Data()));
   	if ( !fElectronCutsPreMC ) {
	  
	    cout<<"El error esta aqui:     "<<Form("ElectronCuts_%s",fElectronCutsPreSelection.Data())<<endl;
	    
	  
	}

       
       if( fElectronCutsPreMC ){ //Just for thesis {
	  
	  
           hTPCdEdxSigbefore_mc 					= (TH2F*)fElectronCutsPreMC->FindObject(Form("Electron_dEdxSignal_before %s",fElectronCutsPreSelection.Data()));
	   hTPCdEdxbefore_mc						= (TH2F*)fElectronCutsPreMC->FindObject(Form("Electron_dEdx_before %s",fElectronCutsPreSelection.Data()));
	   
	   
	   hTPCdEdxSigbefore_mc->Print();
  
	   cout<<"Hola entro a electron preselection "<<endl;
	  
	}
       
       
       
       TList *fElectronCutsMC   = (TList*)fHistosGammaConversionDalitzMC->FindObject(Form("ElectronCuts_%s",fCutSelectionMC.Data() ) );
       
       if( ! fElectronCutsMC ) {
	 
	  fElectronCutsMC   = (TList*)fHistosGammaConversionDalitzMC->FindObject( Form("ElectronCuts_%s",fElectronCutSelectionMC.Data() ) );
	  
	  if( ! fElectronCutsMC ) {
	    cout<<Form("MC: ElectronCuts_%s",fElectronCutSelectionMC.Data())<<" is not found in the file"<<endl;
	    return;
	  }
	  
       }
       
       
       TList *fESDContainerMC   = (TList*)fHistosGammaConversionDalitzMC->FindObject(Form("%s ESD histograms",fCutSelectionMC.Data()));
       
       if(  !fESDContainerMC  ) {
	  
	     cout<<Form("%s ESD histograms",fCutSelectionMC.Data())<<" is not found in the file"<<endl;
	     return;
	 
       }
       
      TList *fTrueHistograms = (TList*)fHistosGammaConversionDalitzMC->FindObject(Form("%s True histograms",fCutSelectionMC.Data()) );
      if(  !fTrueHistograms  ) {
	  
	     cout<<Form("%s True histograms",fCutSelectionMC.Data())<<" was not found in the file"<<endl;
	     return;
	 
       }  
  
      TList *fMCContainer    = (TList*)fHistosGammaConversionDalitzMC->FindObject(Form("%s MC histograms",fCutSelectionMC.Data()) );
      
      if(  !fMCContainer  ) {
	  
	     cout<<Form("%s MC histograms",fCutSelectionMC.Data())<<" was not found in the file"<<endl;
	     return;
	 
       }  
       
       TList  *fQAFolderMC          = (TList*)fHistosGammaConversionDalitzMC->FindObject(Form("%s QA histograms",fCutSelectionMC.Data()));
       
       if( !fQAFolderMC ) {
	 
	    cout<<Form("%s QA histograms",fCutSelectionMC.Data())<<" was not found in the file"<<endl;
	    //return;
	 
       }
       TList *fClusterOutputList = (TList*)fHistosGammaConversionDalitzMC->FindObject(Form("%s Cluster Output",fCutSelectionMC.Data()));
       
       if( !fClusterOutputList ) {
	 
	    cout<<Form("%s Cluster Output",fCutSelectionMC.Data())<<" was not found in the file"<<endl;
	 
       }
       
        
       
       
       
       fEventQualityMC = (TH1F*)fESDContainerMC->FindObject("NEvents");
				
		
		
	if ( energy.CompareTo("PbPb_2.76TeV") == 0 || energy.CompareTo("pPb_5.023TeV") == 0 ){ 
	
		nEventsMC = fEventQualityMC->GetBinContent(1);
	   
	 } else {
	   
		nEventsMC =  GetNEvents(fEventQualityMC);
	 }	
		
		
       
        fESD_NumberOfGoodESDTracksMC = (TH1F*)fESDContainerMC->FindObject("GoodESDTracks");
		
      
        mean_MC = fESD_NumberOfGoodESDTracksMC->GetMean();  
	
	NormFactorMC = 1./nEventsMC*1./mean_MC;
	//NormFactorMC = 1./mean_MC;
	
	
	
	cout<<"Normalization Factor: "<<NormFactorMC<<endl;
	
	cout<<"Mean MC: "<<mean_MC<<endl;
	
	
	
	
	hESDConvGammaEta_mc                             = (TH1F*)fESDContainerMC->FindObject("ESD_ConvGamma_Eta");
	
	if( fMCContainer ) {
	
	hMCPi0EposEnegInvMassPt_mc			= (TH2F*)fMCContainer->FindObject("MC_Pi0EposEneg_InvMassPt");
	hMCConvGammaR	                                = (TH1F*)fMCContainer->FindObject("MC_ConvGamma_R");
	hMCConvGammaPt                                  = (TH1F*)fMCContainer->FindObject("MC_ConvGamma_Pt");
	hMCAllGammaPt					= (TH1F*)fMCContainer->FindObject("MC_AllGamma_Pt");
	hMCConvGammaPtR					= (TH2F*)fMCContainer->FindObject("MC_ConvGamma_Pt_R");
	hMCAllPositronsPt				= (TH1F*)fMCContainer->FindObject("MC_AllPositrons_Pt");
	hMCAllElectronsPt				= (TH1F*)fMCContainer->FindObject("MC_AllElectrons_Pt");
	hMCDecayPositronPi0Pt				= (TH1F*)fMCContainer->FindObject("MC_DecayPositronPi0_Pt");
	hMCDecayElectronPi0Pt				= (TH1F*)fMCContainer->FindObject("MC_DecayElectronPi0_Pt");
	hMCConvGammaPi0Pt				= (TH1F*)fMCContainer->FindObject("MC_ConvGammaPi0_Pt");
	hMCAllGammaPi0Pt				= (TH1F*)fMCContainer->FindObject("MC_AllGammaPi0_Pt");
	
	
	
	fHistoMCMesonPt					= (TH1D*) fMCContainer->FindObject("MC_Pi0_Pt");   
        fHistoMCMesonGGPt               		= (TH1D*) fMCContainer->FindObject("MC_Pi0_GG_Pt");
      
	if( fHistoMCMesonPt )
	fHistoMCMesonDalitzPt           		= (TH1D*) fHistoMCMesonPt->Clone("MC_Pi0_Dalitz_Pt");
	
	}
						
	if( fQAFolderMC ) {
	       
	hESDDalitzElectronAfterPt_mc 			= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_Pt");
	hESDDalitzPositronAfterPt_mc  			= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_Pt");
	hESDDalitzElectronAfterEta_mc 			= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_Eta");
	hESDDalitzPositronAfterEta_mc 			= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_Eta");
	
	hESDDalitzElectronAfterEtaPCut_mc              	= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_Eta_PCut");
	hESDDalitzPositronAfterEtaPCut_mc              	= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_Eta_PCut");
	
	hESDDalitzElectronAfterPhi_mc 			= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_Phi");
	hESDDalitzPositronAfterPhi_mc 			= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_Phi");
      
	hESDDalitzElectronAfterNClsITS_mc		= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_NClsITS");
	hESDDalitzPositronAfterNClsITS_mc		= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_NClsITS");
      
      
	hESDDalitzElectronAfterNFindClsTPC_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_NFindClsTPC");
	hESDDalitzPositronAfterNFindClsTPC_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_NFindClsTPC");
	
	hESDDalitzElectronAfterNFindClsTPCPCut_mc	= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_NFindClsTPC_PCut");
	hESDDalitzPositronAfterNFindClsTPCPCut_mc	= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_NFindClsTPC_PCut");

	hESDDalitzElectronAfterNClsTPC_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_NClsTPC");
	hESDDalitzPositronAfterNClsTPC_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_NClsTPC");
	
	hESDDalitzElectronAfterNClsTPCPCut_mc		= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_NClsTPC_PCut");
	hESDDalitzPositronAfterNClsTPCPCut_mc		= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_NClsTPC_PCut");
	
	hESDDalitzElectronAfterNCrossedRowsTPC_mc	= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_NCrossedRowsTPC");
	hESDDalitzPositronAfterNCrossedRowsTPC_mc	= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_NCrossedRowsTPC");
	
	hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc	= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_NCrossedRowsTPC_PCut");
	hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc     = (TH1F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_NCrossedRowsTPC_PCut");
	
	hESDDalitzPosEleAfterDCAxy_mc			= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzPosEle_After_DCAxy");
	hESDDalitzPosEleAfterDCAz_mc  			= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzPosEle_After_DCAz");
      
	hESDDalitzElectronAfterTPCdEdxVsP_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_TPCdEdxVsP");
        hESDDalitzPositronAfterTPCdEdxVsP_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_TPCdEdxVsP");
	hESDDalitzElectronAfterTPCdEdxSignalVsP_mc	= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_TPCdEdxSignalVsP");
        hESDDalitzPositronAfterTPCdEdxSignalVsP_mc	= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_TPCdEdxSignalVsP");
      
	hESDDalitzElectronAfterTPCdEdxVsEta_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_TPCdEdxVsEta");
        hESDDalitzPositronAfterTPCdEdxVsEta_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_TPCdEdxVsEta");
      
	hESDDalitzElectronAfterTPCdEdxVsPhi_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzElectron_After_TPCdEdxVsPhi");
        hESDDalitzPositronAfterTPCdEdxVsPhi_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_DalitzPositron_After_TPCdEdxVsPhi");
	
        hESDMotherPhi_mc				= (TH1F*)fQAFolderMC->FindObject("ESD_DalitzMother_Phi");
        hESDEposEnegPsiPairDPhi_mc			= (TH2F*)fQAFolderMC->FindObject("ESD_EposEneg_PsiPair_DPhi");
        hESDEposEnegInvMassPt_mc			= (TH2F*)fQAFolderMC->FindObject("ESD_EposEneg_InvMassPt");
        hESDEposEnegLikeSignBackInvMassPt_mc		= (TH2F*)fQAFolderMC->FindObject("ESD_EposEneg_LikeSignBack_InvMassPt");
	
	hESDEposEnegInvMassPi0Pt_mc 			= (TH2F*)fQAFolderMC->FindObject("ESD_EposEneg_InvMassPi0Pt");
	
	sESDConvGammaZR_mc                              = (THnSparseF*)fQAFolderMC->FindObject("ESD_ConvGamma_ZR");
	sESDConvGammaXY_mc 				= (THnSparseF*)fQAFolderMC->FindObject("ESD_ConvGamma_XY");
	
	sESDMotherDalitzPlot_mc				= (THnSparseF*)fQAFolderMC->FindObject("ESD_Mother_DalitzPlot");
	sESDTruePi0DalitzDalitzPlot_mc			= (THnSparseF*)fQAFolderMC->FindObject("ESD_TruePi0Dalitz_DalitzPlot");
	sMCPi0DalitzDalitzPlot_mc			= (THnSparseF*)fQAFolderMC->FindObject("MC_Pi0Dalitz_DalitzPlot");
	
		
	
	hESDEposEnegInvMassPt				= (TH2F*)fQAFolderMC->FindObject("ESD_EposEneg_InvMassPt");
	
			
	
	cout<<"Entro al QA pedrito"<<endl;
	
	}
	
	if( fConvCutsMC ) {
	
	hTPCdEdxSigafter_mc 				= (TH2F*)fConvCutsMC->FindObject(Form("Gamma_dEdxSig_after %s",fGammaCutSelectionMC.Data()));
        hTPCdEdxafter_mc    				= (TH2F*)fConvCutsMC->FindObject(Form("Gamma_dEdx_after %s",fGammaCutSelectionMC.Data()));
	
	}
	
	if( fElectronCutsMC ) {
	  
	   hTPCEledEdxSigafter_mc 					= (TH2F*)fElectronCutsMC->FindObject(Form("Electron_dEdxSignal_after %s",fElectronCutSelection.Data() ));
	   hTPCEledEdxafter_mc 						= (TH2F*)fElectronCutsMC->FindObject(Form("Electron_dEdx_after %s",fElectronCutSelection.Data() ));
	   
	   cout<<Form("Electron_dEdx_after %s",fElectronCutSelection.Data())<<endl;
	   fElectronCutsMC->Print();
	   
	   if (!hTPCEledEdxafter_mc) return ;
    
	}
	
	
	
	///////////////////////////////////////////////True histograms /////////////////////////////////////////////////////////////
       
       gammaCutEffiLegend = "|#eta|<0.9, p_{T,min}=0.05 GeV/c, R_{#gamma conv,min}=5.0cm}";
       electronCutEffiLegend = "|#eta|<0.9";
       
       
       if( fTrueHistograms ) {
       
       hESDEposEnegTruePi0DalitzInvMassPt   		= (TH2F*)fTrueHistograms->FindObject("ESD_EposEneg_TruePi0Dalitz_InvMassPt");
       hESDEposEnegTrueEtaDalitzInvMassPt               = (TH2F*)fTrueHistograms->FindObject("ESD_EposEneg_TrueEtaDalitz_InvMassPt");
       hESDEposEnegTruePhotonInvMassPt			= (TH2F*)fTrueHistograms->FindObject("ESD_EposEneg_TruePhoton_InvMassPt");
       hESDEposEnegTrueInvMassPt			= (TH2F*)fTrueHistograms->FindObject("ESD_EposEneg_True_InvMassPt");
       hESDTrueMotherEposEnegInvMassPt			= (TH2F*)fTrueHistograms->FindObject("ESD_EposEneg_TrueMother_InvMassPt");
				
       		
       			
       
       hESDTrueConvGammaRMC 		  		= (TH1F*)fTrueHistograms->FindObject("ESD_TrueConvGamma_R_MC");
       hESDTrueConvGammaR 				= (TH1F*)fTrueHistograms->FindObject("ESD_TrueConvGamma_R");
       hESDTrueConvGammaPt				= (TH1F*)fTrueHistograms->FindObject("ESD_TrueConvGamma_Pt");
       hESDTruePositronPt				= (TH1F*)fTrueHistograms->FindObject("ESD_TruePositron_Pt");
       hESDTrueElectronPt				= (TH1F*)fTrueHistograms->FindObject("ESD_TrueElectron_Pt");
       hESDEposEnegTruePi0DalitzPsiPairDPhi		= (TH2F*)fTrueHistograms->FindObject("ESD_EposEneg_TruePi0Dalitz_PsiPair_DPhi");
       hESDEposEnegTrueEtaDalitzPsiPairDPhi		= (TH2F*)fTrueHistograms->FindObject("ESD_EposEneg_TrueEtaDalitz_PsiPair_DPhi");
       hESDEposEnegTruePhotonPsiPairDPhi		= (TH2F*)fTrueHistograms->FindObject("ESD_EposEneg_TruePhoton_PsiPair_DPhi");
       hESDTruePi0DalitzElectronPt			= (TH1F*)fTrueHistograms->FindObject("ESD_TruePi0DalitzElectron_Pt");
       hESDTruePi0DalitzPositronPt			= (TH1F*)fTrueHistograms->FindObject("ESD_TruePi0DalitzPositron_Pt");
       hESDTruePi0DalitzConvGammaPt			= (TH1F*)fTrueHistograms->FindObject("ESD_TruePi0DalitzConvGamma_Pt");
       hHistoTruePi0DalitzClusGammaPt			= (TH1F*)fTrueHistograms->FindObject("TruePi0DalitzClusGamma_Pt");
	
			
       }
       
       if( fClusterOutputList ) {
	 
	 hHistoTruePi0DalitzClusGammaPt			= (TH1F*)fClusterOutputList->FindObject("TruePi0DalitzClusGamma_Pt");
	 
	
       }
       
       
       
       
       Double_t Pi0DalitzBR = 0.0;
       
       if( fHistoMCMesonDalitzPt && fHistoMCMesonGGPt ){
	 
	  Pi0DalitzBR =  fHistoMCMesonDalitzPt->GetEntries() / ( fHistoMCMesonDalitzPt->GetEntries() + fHistoMCMesonGGPt->GetEntries() );
       
       }
		  
       
       
       
       if(  hMCConvGammaPtR && hMCAllGammaPt ){
	 
	 
	 
	
	 Double_t RMin[nBins] = {  0.0,   5.0,   7.5,  10.0,  12.5,  20.0,  35.0};
	 Double_t RMax[nBins] = {180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0};
	 
	 TH1F*  hMCConvGammaPtR_Bins[nBins];
	 TH1F*  hMCAllGammaPtRebin;
	
	 
	 hMCAllGammaPtRebin = (TH1F*) hMCAllGammaPt->Rebin(kNbins-1,"hMCConvGammaPtRebin01",xbins);
	 hMCAllGammaPtRebin->Sumw2();
	
	
	 for(UInt_t iBin = 0; iBin < nBins; iBin++){
	 
	    Int_t binMin, binMax;
	 
	    binMin = hMCConvGammaPtR->GetYaxis()->FindBin(RMin[iBin] + 0.05);
	    binMax = hMCConvGammaPtR->GetYaxis()->FindBin(RMax[iBin] - 0.05);
	  
	    TH1F*  temp = (TH1F*)hMCConvGammaPtR->ProjectionX("temp", binMin, binMax);
	    
	    hMCConvGammaEffiVsPtProb[iBin] =(TH1F*)temp->Rebin(kNbins-1,Form("hMCConvGammaEffiVsPtProb_%04d_%04d",(Int_t)(RMin[iBin]*10),(Int_t)(RMax[iBin]*10)),xbins);
	    hMCConvGammaEffiVsPtProb[iBin]->Divide(hMCConvGammaEffiVsPtProb[iBin],hMCAllGammaPtRebin,1,1,"B");
	    
	    cout<<"Pedro gonzalez"<<endl;
	    new TCanvas;
	    
	    hMCConvGammaEffiVsPtProb[iBin]->Draw();
	 
	 }
	
	 
       }
       
       
      
       if( hESDEposEnegPsiPairDPhi_mc && hESDEposEnegTruePi0DalitzPsiPairDPhi && hESDEposEnegTrueEtaDalitzPsiPairDPhi && hESDEposEnegTruePhotonPsiPairDPhi ) {
       
       
       TH1F*  	 hESDEposEnegPsiPair_mc 		= (TH1F*) hESDEposEnegPsiPairDPhi_mc->ProjectionY("ESD_EposEneg_PsiPair_mc");
       TH1F*     hESDEposEnegTruePi0DalitzPsiPair	= (TH1F*) hESDEposEnegTruePi0DalitzPsiPairDPhi->ProjectionY("ESD_EposEneg_TruePi0Dalitz_PsiPair");
       TH1F* 	 hESDEposEnegTrueEtaDalitzPsiPair	= (TH1F*) hESDEposEnegTrueEtaDalitzPsiPairDPhi->ProjectionY("ESD_EposEneg_TrueEtaDalitz_PsiPair");
       TH1F*     hESDEposEnegTruePhotonPsiPair          = (TH1F*) hESDEposEnegTruePhotonPsiPairDPhi->ProjectionY("ESD_EposEneg_TruePhoton_PsiPair");
       
       TH1F*  	 hESDEposEnegDPhi_mc 			= (TH1F*) hESDEposEnegPsiPairDPhi_mc->ProjectionX("ESD_EposEneg_PsiPair_mc");
       TH1F*     hESDEposEnegTruePi0DalitzDPhi		= (TH1F*) hESDEposEnegTruePi0DalitzPsiPairDPhi->ProjectionX("ESD_EposEneg_TruePi0Dalitz_PsiPair");
       TH1F* 	 hESDEposEnegTrueEtaDalitzDPhi		= (TH1F*) hESDEposEnegTrueEtaDalitzPsiPairDPhi->ProjectionX("ESD_EposEneg_TrueEtaDalitz_PsiPair");
       TH1F*     hESDEposEnegTruePhotonDPhi          	= (TH1F*) hESDEposEnegTruePhotonPsiPairDPhi->ProjectionX("ESD_EposEneg_TruePhoton_PsiPair");
       
	
	 
	 
       hESDEposEnegTruePi0DalitzPsiPairDPhi->SetMarkerColor(kGreen);
       hESDEposEnegTruePi0DalitzPsiPairDPhi->SetMarkerSize(0.1);
       //hESDEposEnegTruePi0DalitzPsiPairDPhi->SetMarkerColorAlpha(kGreen,0.55);
       hESDEposEnegTruePhotonPsiPairDPhi->SetMarkerColor(kRed);
       hESDEposEnegTruePhotonPsiPairDPhi->SetMarkerSize(0.1);
       //
       hESDEposEnegTruePhotonPsiPairDPhi->SetMarkerColorAlpha(kRed,0.55);
       
       hESDEposEnegTruePi0DalitzPsiPairDPhi->SetMarkerStyle(20);
       hESDEposEnegTruePhotonPsiPairDPhi->SetMarkerStyle(20);
       
        
       hESDEposEnegTruePi0DalitzPsiPairDPhi->SetLineColor(kGreen);
       hESDEposEnegTruePhotonPsiPairDPhi->SetLineColor(kRed);
       
       hESDEposEnegTruePi0DalitzPsiPairDPhi->SetTitle("");
       
       hESDEposEnegTruePi0DalitzPsiPairDPhi->SetXTitle("#Delta#phi");
       hESDEposEnegTruePi0DalitzPsiPairDPhi->SetYTitle("#psi-pair");
       hESDEposEnegTruePhotonPsiPairDPhi->SetTitle("");
       
       
       TCanvas * canvasEposENegPsiPairDPhi = new TCanvas("canvasEposENegPsiPairDPhi","",10,10,500,500);  // gives the page size
       canvasEposENegPsiPairDPhi->SetLeftMargin(0.12);
       canvasEposENegPsiPairDPhi->SetRightMargin(0.02);
       canvasEposENegPsiPairDPhi->SetBottomMargin(0.14);
       canvasEposENegPsiPairDPhi->SetTopMargin(0.02);
       
       
       
       canvasEposENegPsiPairDPhi->cd();
       
       TH2F* plot2DCanvas01 = new TH2F("plot2DCanvas01","plot2DCanvas01",400,-2.0,2.0,400,-2.0,2.0);
       SetStyleHistoTH2ForGraphs(plot2DCanvas01,"#Delta#phi","#psi_{pair}", 0.032,0.04, 0.032,0.04, 1.0,1.2);
	
	
       hESDEposEnegTruePi0DalitzPsiPairDPhi->SetStats(kFALSE);
       hESDEposEnegTruePhotonPsiPairDPhi->SetStats(kFALSE);
       
	
	
	
	TGraph2D* graph2DEposEnegTruePi0DalitzPsiPairDPhi = new TGraph2D(hESDEposEnegTruePi0DalitzPsiPairDPhi);
	TGraph2D* graph2DEposEnegTruePhotonPsiPairDPhi    = new TGraph2D(hESDEposEnegTruePhotonPsiPairDPhi);
	Double_t x[4]={.0,.0,.12,.0};
	Double_t y[4]={.6,-.6,.0,.6};
	TPolyLine *pline = new TPolyLine(4,x,y);
	pline->SetLineColor(kBlue);
	pline->SetLineWidth(4);
	
	TArrow * arrowPsi00 = new TArrow(-0.7,0.02,-0.03,0.57,0.02,"|>");
	arrowPsi00->SetAngle(60);
	arrowPsi00->SetFillColor(42);
	arrowPsi00->SetLineWidth(3);
	
	TArrow * arrowPsi01 = new TArrow(-0.7,0.02,-0.03,-0.57,0.02,"|>");
	arrowPsi01->SetAngle(60);
	arrowPsi01->SetFillColor(42);
	arrowPsi01->SetLineWidth(3);
	
	TArrow * arrowPhi00 = new TArrow(0.0,-0.8,-0.00,-0.62,0.02,"|>");
	arrowPhi00->SetAngle(60);
	arrowPhi00->SetFillColor(42);
	arrowPhi00->SetLineWidth(3);
	
	TArrow * arrowPhi01 = new TArrow(0.12,-0.8,0.12,-0.02,0.02,"|>");
	arrowPhi01->SetAngle(60);
	arrowPhi01->SetFillColor(42);
	arrowPhi01->SetLineWidth(3);
	
	
	TLatex *process00 = new TLatex(0.63, 0.92,Form("MC, %s",collisionSystem.Data()));
        process00->SetNDC();
        process00->SetTextColor(1);
        process00->SetTextSize(0.04);
	process00->SetTextFont(42);
	
	TLatex *Psi00 = new TLatex(0.20, 0.56,"|#Psi_{0}|");
	Psi00->SetNDC();
	Psi00->SetTextFont(42);
        Psi00->SetTextColor(1);
        Psi00->SetTextSize(0.04);
	TLatex *DeltaPhi00 = new TLatex(0.54, 0.17,"#phi_{0}");
	DeltaPhi00->SetNDC();
	DeltaPhi00->SetTextFont(42);
        DeltaPhi00->SetTextColor(1);
        DeltaPhi00->SetTextSize(0.04);
	TLatex *DeltaPhi01 = new TLatex(0.59, 0.17,"#phi_{1}");
	DeltaPhi01->SetNDC();
	DeltaPhi01->SetTextFont(42);
        DeltaPhi01->SetTextColor(1);
        DeltaPhi01->SetTextSize(0.04);
	
	
	
	
	
	graph2DEposEnegTruePi0DalitzPsiPairDPhi->SetNpx(20);
	graph2DEposEnegTruePi0DalitzPsiPairDPhi->SetNpy(20);
	graph2DEposEnegTruePhotonPsiPairDPhi->SetNpx(20);
	graph2DEposEnegTruePhotonPsiPairDPhi->SetNpy(20);
	
	plot2DCanvas01->GetXaxis()->SetRangeUser(-1.1,1.1);
	plot2DCanvas01->GetYaxis()->SetRangeUser(-1.0,1.0);
	plot2DCanvas01->Draw();
	
	
	graph2DEposEnegTruePi0DalitzPsiPairDPhi->Draw("same");
	graph2DEposEnegTruePhotonPsiPairDPhi->Draw("same");
	
	pline->Draw("sames");
	
	process00->Draw();

	Psi00->Draw("sames");
	DeltaPhi00->Draw("sames");
	DeltaPhi01->Draw("sames");
	arrowPsi00->Draw();
	arrowPsi01->Draw();
	arrowPhi00->Draw();
	arrowPhi01->Draw();
	
	
	
	TLegend* leg2 = new TLegend(0.15,0.86,0.32,0.95);
	leg2->SetTextSize(0.04);
	leg2->SetTextFont(42);			
	leg2->SetFillColor(0);
	leg2->SetFillStyle(0);
	leg2->SetLineColor(0);
	
	leg2->AddEntry(graph2DEposEnegTruePi0DalitzPsiPairDPhi,"#pi^{0} Dalitz","pf");
	leg2->AddEntry(graph2DEposEnegTruePhotonPsiPairDPhi,"#gamma_{conv}","pf");
	leg2->Draw("sames");
	
	
	canvasEposENegPsiPairDPhi->Update();
	canvasEposENegPsiPairDPhi->SaveAs(Form("%s/ESD_EposEneg_TruePi0_Dalitz_PsiPair_DPhi%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasEposENegPsiPairDPhi;	
	
	
	
	
	
	
       TCanvas * canvasEposENegPsiPair = new TCanvas("canvasEposENegPsiPair","",10,10,500,500);  // gives the page size
       canvasEposENegPsiPair->SetLeftMargin(0.10);
       canvasEposENegPsiPair->SetRightMargin(0.02);
       canvasEposENegPsiPair->SetBottomMargin(0.14);
       canvasEposENegPsiPair->SetTopMargin(0.02);
       canvasEposENegPsiPair->SetLogy(1);
       
       canvasEposENegPsiPair->cd();
       
       TH1F* plot1DCanvas01 = new TH1F("plot1DCanvas01","plot1DCanvas01",200,-2.0,2.0);
     
       SetStyleHistoTH1ForGraphs(plot1DCanvas01,"#Psi_{pair} RAD","counts", 0.032,0.04, 0.032,0.04, 1.0,1.0);
	
	
	
	plot1DCanvas01->GetYaxis()->SetRangeUser(1,6.5e+7);
	plot1DCanvas01->GetXaxis()->SetRangeUser(-0.81,0.81);
	plot1DCanvas01->SetLineWidth(2);
	
	TH1D * hESDEposEnegPsiPairClone_mc = (TH1D*)hESDEposEnegPsiPair_mc->Clone("hESDEposEnegPsiPairClone_mc");
	Double_t nEntries = hESDEposEnegPsiPairClone_mc->GetEffectiveEntries();
	
	//hESDEposEnegPsiPairClone_mc->Scale(1./nEntries);
	hESDEposEnegPsiPairClone_mc->SetFillStyle(3002);
	hESDEposEnegPsiPairClone_mc->SetFillColor(18);
	hESDEposEnegPsiPairClone_mc->SetLineColor(kBlack);
	
	TH1D* hESDEposEnegTruePi0DalitzPsiPairClone = (TH1D*)hESDEposEnegTruePi0DalitzPsiPair->Clone("hESDEposEnegTruePi0DalitzPsiPairClone");
	
	//hESDEposEnegTruePi0DalitzPsiPairClone->Scale(1./nEntries);
        hESDEposEnegTruePi0DalitzPsiPairClone->SetLineWidth(2);
	hESDEposEnegTruePi0DalitzPsiPairClone->SetFillStyle(3002);
	hESDEposEnegTruePi0DalitzPsiPairClone->SetFillColor(kGreen); 
        hESDEposEnegTruePi0DalitzPsiPairClone->SetLineColor(kGreen);
	
	TH1D* hESDEposEnegTruePhotonPsiPairClone = (TH1D*)hESDEposEnegTruePhotonPsiPair->Clone("hESDEposEnegTruePhotonPsiPairClone");
	
	//hESDEposEnegTruePhotonPsiPairClone->Scale(1./nEntries);
       	hESDEposEnegTruePhotonPsiPairClone->SetFillStyle(3002);
	hESDEposEnegTruePhotonPsiPairClone->SetFillColor(kRed); 
        hESDEposEnegTruePhotonPsiPairClone->SetLineColor(kRed);
       
     
	plot1DCanvas01->Draw();
	hESDEposEnegPsiPairClone_mc->Draw("same");
        hESDEposEnegTruePi0DalitzPsiPairClone->Draw("same");
	hESDEposEnegTruePhotonPsiPairClone->Draw("same");
	
	
        TLatex *process01 = new TLatex(0.63, 0.92,Form("MC, %s",collisionSystem.Data()));
        process01->SetNDC();
        process01->SetTextColor(1);
        process01->SetTextSize(0.04);
	process01->SetTextFont(42);
	
	
	
	TLegend*  legendEposENegPsiPair = new TLegend(0.14,0.80,0.5,0.93);
	legendEposENegPsiPair->SetTextSize(0.04);
	legendEposENegPsiPair->SetTextFont(42);
	legendEposENegPsiPair->SetFillColor(0);
	legendEposENegPsiPair->SetLineColor(0);
	legendEposENegPsiPair->SetFillStyle(0);
	legendEposENegPsiPair->AddEntry(hESDEposEnegPsiPairClone_mc,"All candidates","f");
	legendEposENegPsiPair->AddEntry(hESDEposEnegTruePi0DalitzPsiPairClone,"#pi^{0} Dalitz","f");
	legendEposENegPsiPair->AddEntry(hESDEposEnegTruePhotonPsiPairClone,"#gamma_{conv}","f");
	legendEposENegPsiPair->Draw();
	process01->Draw();
	
	
	canvasEposENegPsiPair->Update();
	canvasEposENegPsiPair->SaveAs(Form("%s/ESD_EposEneg_TruePi0_Dalitz_PsiPair%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasEposENegPsiPair;	
	
	
	
	
       TCanvas * canvasEposENegDPhi = new TCanvas("canvasEposENegDPhi","",10,10,500,500);  // gives the page size
       canvasEposENegDPhi->SetLeftMargin(0.10);
       canvasEposENegDPhi->SetRightMargin(0.02);
       canvasEposENegDPhi->SetBottomMargin(0.14);
       canvasEposENegDPhi->SetTopMargin(0.02);
       canvasEposENegDPhi->SetLogy(1);
       
       canvasEposENegDPhi->cd();
       
       TH1F* plot1DCanvas02 = new TH1F("plot1DCanvas02","plot1DCanvas02",200,-2.0,2.0);
     
       SetStyleHistoTH1ForGraphs(plot1DCanvas02,"#Delta#Phi RAD","counts", 0.032,0.04, 0.032,0.04, 1.0,1.0);
	
       plot1DCanvas02->GetYaxis()->SetRangeUser(1,6.5e+7);
       plot1DCanvas02->GetXaxis()->SetRangeUser(-0.81,0.81);
       plot1DCanvas02->SetLineWidth(2);
	
        
	
	hESDEposEnegDPhi_mc->SetFillStyle(3002);
	hESDEposEnegDPhi_mc->SetFillColor(18);
	hESDEposEnegDPhi_mc->SetLineColor(kBlack);
	
	
	hESDEposEnegTruePi0DalitzDPhi->SetFillStyle(3002);
	hESDEposEnegTruePi0DalitzDPhi->SetFillColor(kGreen); 
        hESDEposEnegTruePi0DalitzDPhi->SetLineColor(kGreen);
	
	
       	hESDEposEnegTruePhotonDPhi->SetFillStyle(3002);
	hESDEposEnegTruePhotonDPhi->SetFillColor(kRed); 
        hESDEposEnegTruePhotonDPhi->SetLineColor(kRed);
       
     
	
	hESDEposEnegTruePi0DalitzDPhi->SetLineWidth(2);
	hESDEposEnegTruePhotonDPhi->SetLineWidth(2);
		
	
	plot1DCanvas02->Draw();
	hESDEposEnegDPhi_mc->Draw("same");
	hESDEposEnegTruePhotonDPhi->Draw("same");
	hESDEposEnegTruePi0DalitzDPhi->Draw("same");
	
	
	
		
	TLegend*  legendEposENegDPhi = new TLegend(0.15,0.80,0.5,0.93);
	legendEposENegDPhi->SetTextSize(0.04);
	legendEposENegDPhi->SetTextFont(42);
	legendEposENegDPhi->SetFillColor(0);
	legendEposENegDPhi->SetLineColor(0);
	legendEposENegDPhi->SetFillStyle(0);
	legendEposENegDPhi->AddEntry(hESDEposEnegDPhi_mc,"All candidates","f");
	legendEposENegDPhi->AddEntry(hESDEposEnegTruePi0DalitzDPhi,"#pi^{0} Dalitz","f");
	legendEposENegDPhi->AddEntry(hESDEposEnegTruePhotonDPhi,"#gamma_{conv}","f");
	legendEposENegDPhi->Draw();
	process01->Draw();
	
	
	canvasEposENegDPhi->Update();
	canvasEposENegDPhi->SaveAs(Form("%s/ESD_EposEneg_TruePi0_Dalitz_DPhi%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasEposENegDPhi;	
	
	
	
       }	
	
	
	if( hESDTrueConvGammaRMC && hMCConvGammaR  ){
	  
	/*hESDTrueConvGammaPtvsEffi = (TH1F*)hESDTrueConvGammaPt->Rebin(kNbins-1,"hESDTrueConvGammaPtvsEffi",xbins);
	hESDTrueConvGammaPtvsEffi->Sumw2();
	
	TH1F* hMCConvGammaPtRebin = (TH1F*) hMCConvGammaPt->Rebin(kNbins-1,"hMCConvGammaPtRebin",xbins);
	hMCConvGammaPtRebin->Sumw2();*/
	  
	  
	hESDTrueConvGammaRvsEffi = (TH1F*)hESDTrueConvGammaRMC->Rebin(kNRbins-1,"hESDTrueConvGammaRvsEffi",xRbins);
	hESDTrueConvGammaRvsEffi->Sumw2();
	
	TH1F* hMCConvGammaRRebin = (TH1F*) hMCConvGammaR->Rebin(kNRbins-1,"hMCConvGammaRRebin",xRbins);
	hMCConvGammaR->Sumw2();
	
	hESDTrueConvGammaRvsEffi->Divide(hESDTrueConvGammaRvsEffi,hMCConvGammaRRebin,1,1,"B");
	  
	  
	TCanvas * canvasConvGammaRvsEffi = new TCanvas("canvasConvGammaRvsEffi","",10,10,600,470);  // gives the page size
	
	canvasConvGammaRvsEffi->SetLeftMargin(0.12);
	canvasConvGammaRvsEffi->SetRightMargin(0.03);
	canvasConvGammaRvsEffi->SetTopMargin(0.03);
	canvasConvGammaRvsEffi->SetBottomMargin(0.1);
	
	TLatex *latexConvGammaRvsEffi = new TLatex(0.35,0.30,gammaCutEffiLegend.Data());
        SetStyleTLatex( latexConvGammaRvsEffi, 0.038,4);
	
	canvasConvGammaRvsEffi->cd();
	
	hESDTrueConvGammaRvsEffi->SetStats(kFALSE);
	
        hESDTrueConvGammaRvsEffi->SetXTitle("R_{#gammaconv}");
	hESDTrueConvGammaRvsEffi->SetYTitle("Efficiency");
	
	
	hESDTrueConvGammaRvsEffi->GetYaxis()->SetRangeUser(0.01,1.0);
	hESDTrueConvGammaRvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hESDTrueConvGammaRvsEffi->GetYaxis()->SetTitleOffset(1.5);
	
	
	hESDTrueConvGammaRvsEffi->SetLineColor(kRed);
	hESDTrueConvGammaRvsEffi->SetMarkerStyle(7);
	hESDTrueConvGammaRvsEffi->SetMarkerColor(kRed);
	hESDTrueConvGammaRvsEffi->Draw();
	latexConvGammaRvsEffi->Draw();
	
	
	
	canvasConvGammaRvsEffi->Update();
	canvasConvGammaRvsEffi->SaveAs(Form("%s/GammaConv_Effi_Vs_Radius%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasConvGammaRvsEffi;		
	
	  
		  
	  
	}
	
	
	if( hESDTrueConvGammaPt && hMCConvGammaPt  ){
	  
	  
				    
	hESDTrueConvGammaPtvsEffi = (TH1F*)hESDTrueConvGammaPt->Rebin(kNbins-1,"hESDTrueConvGammaPtvsEffi",xbins);
	hESDTrueConvGammaPtvsEffi->Sumw2();
	
	TH1F* hMCConvGammaPtRebin = (TH1F*) hMCConvGammaPt->Rebin(kNbins-1,"hMCConvGammaPtRebin",xbins);
	hMCConvGammaPtRebin->Sumw2();
	
	hESDTrueConvGammaPtvsEffi->Divide(hESDTrueConvGammaPtvsEffi,hMCConvGammaPtRebin,1,1,"B");
	  
	  
	TCanvas * canvasConvGammaPtvsEffi = new TCanvas("canvasConvGammaPtvsEffi","",10,10,600,470);  // gives the page size
	canvasConvGammaPtvsEffi->SetLeftMargin(0.12);
	canvasConvGammaPtvsEffi->SetRightMargin(0.03);
	canvasConvGammaPtvsEffi->SetTopMargin(0.03);
	canvasConvGammaPtvsEffi->SetBottomMargin(0.1);
	
	canvasConvGammaPtvsEffi->cd();
	
	TLatex *latexConvGammaPtvsEffi = new TLatex(0.30,0.4,gammaCutEffiLegend.Data());
        SetStyleTLatex( latexConvGammaPtvsEffi, 0.038,4);
	
	hESDTrueConvGammaPtvsEffi->SetStats(kFALSE);
	
        hESDTrueConvGammaPtvsEffi->SetXTitle("p_{T} (GeV/c)");
	hESDTrueConvGammaPtvsEffi->SetYTitle("#gamma_{converted}  Efficiency");
	
	
	hESDTrueConvGammaPtvsEffi->GetYaxis()->SetRangeUser(0.01,1.0);
	hESDTrueConvGammaPtvsEffi->GetXaxis()->SetRangeUser(0.0,10.0);
	hESDTrueConvGammaPtvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hESDTrueConvGammaPtvsEffi->GetYaxis()->SetTitleOffset(1.5);
	
	hESDTrueConvGammaPtvsEffi->SetLineColor(kRed);
	hESDTrueConvGammaPtvsEffi->SetMarkerStyle(22);
	hESDTrueConvGammaPtvsEffi->SetMarkerColor(kRed);
	hESDTrueConvGammaPtvsEffi->Draw();
	latexConvGammaPtvsEffi->Draw();
	
	
	canvasConvGammaPtvsEffi->Update();
	canvasConvGammaPtvsEffi->SaveAs(Form("%s/GammaConv_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasConvGammaPtvsEffi;		
	
	  
		  
	  
	}
	
	
	
	if( hESDTruePi0DalitzConvGammaPt && hMCConvGammaPi0Pt  ){
	  
	  				    
	hESDTruePi0DalitzConvGammaPtvsEffi = (TH1F*)hESDTruePi0DalitzConvGammaPt->Rebin(kNbins-1,"hESDTruePi0DalitzConvGammaPtvsEffi",xbins);
	hESDTruePi0DalitzConvGammaPtvsEffi->Sumw2();
	
	TH1F* hMCConvGammaPi0PtRebin = (TH1F*) hMCConvGammaPi0Pt->Rebin(kNbins-1,"hMCConvGammaPi0PtRebin",xbins);
	hMCConvGammaPi0PtRebin->Sumw2();
	
	hESDTruePi0DalitzConvGammaPtvsEffi->Divide(hESDTruePi0DalitzConvGammaPtvsEffi,hMCConvGammaPi0PtRebin,1,1,"B");
	  
	  
	TCanvas * canvasTruePi0DalitzConvGammaPtvsEffi = new TCanvas("canvasTruePi0DalitzConvGammaPtvsEffi","",10,10,600,470);  // gives the page size
	canvasTruePi0DalitzConvGammaPtvsEffi->SetLeftMargin(0.12);
	canvasTruePi0DalitzConvGammaPtvsEffi->SetRightMargin(0.03);
	canvasTruePi0DalitzConvGammaPtvsEffi->SetTopMargin(0.03);
	canvasTruePi0DalitzConvGammaPtvsEffi->SetBottomMargin(0.1);
	
	canvasTruePi0DalitzConvGammaPtvsEffi->cd();
	
	TLatex *latexTruePi0DalitzConvGammaPtvsEffi = new TLatex(0.30,0.4,gammaCutEffiLegend.Data());
        SetStyleTLatex( latexTruePi0DalitzConvGammaPtvsEffi, 0.038,4);
	
	hESDTruePi0DalitzConvGammaPtvsEffi->SetStats(kFALSE);
	
        hESDTruePi0DalitzConvGammaPtvsEffi->SetXTitle("p_{T} (GeV/c)");
	hESDTruePi0DalitzConvGammaPtvsEffi->SetYTitle("#pi^{0}_#gamma_{converted}  Efficiency");
	
	
	hESDTruePi0DalitzConvGammaPtvsEffi->GetYaxis()->SetRangeUser(0.01,1.0);
	hESDTruePi0DalitzConvGammaPtvsEffi->GetXaxis()->SetRangeUser(0.0,10.0);
	hESDTruePi0DalitzConvGammaPtvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hESDTruePi0DalitzConvGammaPtvsEffi->GetYaxis()->SetTitleOffset(1.5);
	
	hESDTruePi0DalitzConvGammaPtvsEffi->SetLineColor(kRed);
	hESDTruePi0DalitzConvGammaPtvsEffi->SetMarkerStyle(22);
	hESDTruePi0DalitzConvGammaPtvsEffi->SetMarkerColor(kRed);
	hESDTruePi0DalitzConvGammaPtvsEffi->Draw();
	latexTruePi0DalitzConvGammaPtvsEffi->Draw();
	
	
	canvasTruePi0DalitzConvGammaPtvsEffi->Update();
	canvasTruePi0DalitzConvGammaPtvsEffi->SaveAs(Form("%s/TruePi0DalitzGammaConv_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasTruePi0DalitzConvGammaPtvsEffi;		
	
	  
		  
	  
	}
	
	
	
	if( hESDTruePi0DalitzConvGammaPt && hMCAllGammaPi0Pt  ){
	  
	  				    
	hESDTruePi0DalitzGammaPtvsEffi = (TH1F*)hESDTruePi0DalitzConvGammaPt->Rebin(kNbins-1,"hESDTruePi0DalitzGammaPtvsEffi",xbins);
	hESDTruePi0DalitzGammaPtvsEffi->Sumw2();
	
	TH1F* hMCAllGammaPi0PtRebin = (TH1F*) hMCAllGammaPi0Pt->Rebin(kNbins-1,"hMCAllGammaPi0PtRebin",xbins);
	hMCAllGammaPi0PtRebin->Sumw2();
	
	hESDTruePi0DalitzGammaPtvsEffi->Divide(hESDTruePi0DalitzGammaPtvsEffi,hMCAllGammaPi0PtRebin,1,1,"B");
	  
	  
	TCanvas * canvasTruePi0DalitzGammaPtvsEffi = new TCanvas("canvasTruePi0DalitzGammaPtvsEffi","",10,10,600,470);  // gives the page size
	canvasTruePi0DalitzGammaPtvsEffi->SetLeftMargin(0.12);
	canvasTruePi0DalitzGammaPtvsEffi->SetRightMargin(0.03);
	canvasTruePi0DalitzGammaPtvsEffi->SetTopMargin(0.03);
	canvasTruePi0DalitzGammaPtvsEffi->SetBottomMargin(0.1);
	
	canvasTruePi0DalitzGammaPtvsEffi->cd();
	
	TLatex *latexTruePi0DalitzGammaPtvsEffi = new TLatex(0.30,0.4,gammaCutEffiLegend.Data());
        SetStyleTLatex( latexTruePi0DalitzGammaPtvsEffi, 0.038,4);
	
	hESDTruePi0DalitzGammaPtvsEffi->SetStats(kFALSE);
	
        hESDTruePi0DalitzGammaPtvsEffi->SetXTitle("p_{T} (GeV/c)");
	hESDTruePi0DalitzGammaPtvsEffi->SetYTitle("#pi^{0}_#gamma  Efficiency");
	
	
	hESDTruePi0DalitzGammaPtvsEffi->GetYaxis()->SetRangeUser(0.01,1.0);
	hESDTruePi0DalitzGammaPtvsEffi->GetXaxis()->SetRangeUser(0.0,10.0);
	hESDTruePi0DalitzGammaPtvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hESDTruePi0DalitzGammaPtvsEffi->GetYaxis()->SetTitleOffset(1.5);
	
	hESDTruePi0DalitzGammaPtvsEffi->SetLineColor(kRed);
	hESDTruePi0DalitzGammaPtvsEffi->SetMarkerStyle(22);
	hESDTruePi0DalitzGammaPtvsEffi->SetMarkerColor(kRed);
	hESDTruePi0DalitzGammaPtvsEffi->Draw();
	latexTruePi0DalitzGammaPtvsEffi->Draw();
	
	
	canvasTruePi0DalitzGammaPtvsEffi->Update();
	canvasTruePi0DalitzGammaPtvsEffi->SaveAs(Form("%s/TruePi0DalitzGamma_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasTruePi0DalitzGammaPtvsEffi;		
	
	  
		  
	  
	}
	
	
	
	
	
	if( hMCConvGammaPt && hMCAllGammaPt ){
	
	
	
	  
	  
	  
	hMCGammaPtvsConvProb = (TH1F*) hMCConvGammaPt->Rebin(kNbins-1,"hMCGammaPtvsConvProb",xbins);
	hMCGammaPtvsConvProb->Sumw2();
		
	TH1F* hMCAllGammaPtRebin = (TH1F*) hMCAllGammaPt->Rebin(kNbins-1,"hMCAllGammaPtRebin",xbins);	
	hMCAllGammaPtRebin->Sumw2();
	
	hMCGammaPtvsConvProb->Divide(hMCGammaPtvsConvProb,hMCAllGammaPtRebin,1,1,"B");
	  
	  
	TCanvas * canvasGammaPtvsConvProb = new TCanvas("canvasGammaPtvsConvProb","",10,10,600,470);  // gives the page size
	
	canvasGammaPtvsConvProb->SetLeftMargin(0.12);
	canvasGammaPtvsConvProb->SetRightMargin(0.03);
	canvasGammaPtvsConvProb->SetTopMargin(0.03);
	canvasGammaPtvsConvProb->SetBottomMargin(0.1);
	
	TLatex *latexGammaPtvsConvProb = new TLatex(0.30,0.45,gammaCutEffiLegend.Data());
        SetStyleTLatex( latexGammaPtvsConvProb, 0.038,4);
	
	
	
   
	
	
	
	canvasGammaPtvsConvProb->cd();
	
	hMCGammaPtvsConvProb->SetStats(kFALSE);
	
        hMCGammaPtvsConvProb->SetXTitle("p_{T} (GeV/c)");
	hMCGammaPtvsConvProb->SetYTitle("Conversion probability");
	

	
	hMCGammaPtvsConvProb->GetYaxis()->SetRangeUser(0.001,0.1);
	hMCGammaPtvsConvProb->GetXaxis()->SetRangeUser(0.0,10.0);
	
	hMCGammaPtvsConvProb->GetXaxis()->SetTitleOffset(1.3);
	hMCGammaPtvsConvProb->GetYaxis()->SetTitleOffset(1.5);
	
	
	hMCGammaPtvsConvProb->SetLineColor(kRed);
	hMCGammaPtvsConvProb->SetMarkerStyle(22);
	hMCGammaPtvsConvProb->SetMarkerColor(kRed);
	hMCGammaPtvsConvProb->Draw();
	
	latexGammaPtvsConvProb->Draw();
	
	
	canvasGammaPtvsConvProb->Update();
	canvasGammaPtvsConvProb->SaveAs(Form("%s/GammaConv_ConvProb_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasGammaPtvsConvProb;		
	  
	    
	  
	}
	
	
	
	if( hESDTrueConvGammaPt && hMCConvGammaPt ){
	  
	  
				    
	hESDTrueConvGammaPtvsEffi = (TH1F*)hESDTrueConvGammaPt->Rebin(kNbins-1,"hESDTrueConvGammaPtvsEffi",xbins);
	hESDTrueConvGammaPtvsEffi->Sumw2();
	
	TH1F* hMCConvGammaPtRebin = (TH1F*) hMCConvGammaPt->Rebin(kNbins-1,"hMCConvGammaPtRebin",xbins);
	hMCConvGammaPtRebin->Sumw2();
	
	hESDTrueConvGammaPtvsEffi->Divide(hESDTrueConvGammaPtvsEffi,hMCConvGammaPtRebin,1,1,"B");
	  
	  
	TCanvas * canvasConvGammaPtvsEffi = new TCanvas("canvasConvGammaPtvsEffi","",10,10,600,470);  // gives the page size
	canvasConvGammaPtvsEffi->SetLeftMargin(0.12);
	canvasConvGammaPtvsEffi->SetRightMargin(0.03);
	canvasConvGammaPtvsEffi->SetTopMargin(0.03);
	canvasConvGammaPtvsEffi->SetBottomMargin(0.1);
	
	canvasConvGammaPtvsEffi->cd();
	
	TLatex *latexConvGammaPtvsEffi = new TLatex(0.30,0.4,gammaCutEffiLegend.Data());
        SetStyleTLatex( latexConvGammaPtvsEffi, 0.038,4);
	
	hESDTrueConvGammaPtvsEffi->SetStats(kFALSE);
	
        hESDTrueConvGammaPtvsEffi->SetXTitle("p_{T} (GeV/c)");
	hESDTrueConvGammaPtvsEffi->SetYTitle("#gamma_{converted}  Efficiency");
	
	
	hESDTrueConvGammaPtvsEffi->GetYaxis()->SetRangeUser(0.01,1.0);
	hESDTrueConvGammaPtvsEffi->GetXaxis()->SetRangeUser(0.0,10.0);
	hESDTrueConvGammaPtvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hESDTrueConvGammaPtvsEffi->GetYaxis()->SetTitleOffset(1.5);
	
	hESDTrueConvGammaPtvsEffi->SetLineColor(kRed);
	hESDTrueConvGammaPtvsEffi->SetMarkerStyle(22);
	hESDTrueConvGammaPtvsEffi->SetMarkerColor(kRed);
	hESDTrueConvGammaPtvsEffi->Draw();
	latexConvGammaPtvsEffi->Draw();
	
	
	canvasConvGammaPtvsEffi->Update();
	canvasConvGammaPtvsEffi->SaveAs(Form("%s/GammaConv_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasConvGammaPtvsEffi;		
	
	  
		  
	  
	}
	
	/////////////////////////////////////////////////////////////////////////////
	
	if( hHistoTruePi0DalitzClusGammaPt && hMCAllGammaPi0Pt ){
	  
	  
	hHistoTruePi0DalitzClusGammaPtvsEffi = (TH1F*)hHistoTruePi0DalitzClusGammaPt->Rebin(kNbins-1,"hHistoTruePi0DalitzClusGammaPtvsEffi",xbins);
	hHistoTruePi0DalitzClusGammaPtvsEffi->Sumw2();
	
	TH1F* hMCAllGammaPi0PtRebin = (TH1F*) hMCAllGammaPi0Pt->Rebin(kNbins-1,"hMCAllGammaPi0PtRebin",xbins);
	hMCAllGammaPi0PtRebin->Sumw2();
	
	hHistoTruePi0DalitzClusGammaPtvsEffi->Divide(hHistoTruePi0DalitzClusGammaPtvsEffi,hMCAllGammaPi0PtRebin,1,1,"B");
	  
	  
	TCanvas * canvasTruePi0DalitzClusGammaPtvsEffi = new TCanvas("canvasTruePi0DalitzClusGammaPtvsEffi","",10,10,600,470);  // gives the page size
	canvasTruePi0DalitzClusGammaPtvsEffi->SetLeftMargin(0.12);
	canvasTruePi0DalitzClusGammaPtvsEffi->SetRightMargin(0.03);
	canvasTruePi0DalitzClusGammaPtvsEffi->SetTopMargin(0.03);
	canvasTruePi0DalitzClusGammaPtvsEffi->SetBottomMargin(0.1);
	
	canvasTruePi0DalitzClusGammaPtvsEffi->cd();
	
	TLatex *latexTruePi0DalitzClusGammaPtvsEffi = new TLatex(0.30,0.4,gammaCutEffiLegend.Data());
        SetStyleTLatex( latexTruePi0DalitzClusGammaPtvsEffi, 0.038,4);
	
	hHistoTruePi0DalitzClusGammaPtvsEffi->SetStats(kFALSE);
	
        hHistoTruePi0DalitzClusGammaPtvsEffi->SetXTitle("p_{T} (GeV/c)");
	hHistoTruePi0DalitzClusGammaPtvsEffi->SetYTitle("#gamma_{converted}  Efficiency");
	
	
	hHistoTruePi0DalitzClusGammaPtvsEffi->GetYaxis()->SetRangeUser(0.01,1.0);
	hHistoTruePi0DalitzClusGammaPtvsEffi->GetXaxis()->SetRangeUser(0.0,10.0);
	hHistoTruePi0DalitzClusGammaPtvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hHistoTruePi0DalitzClusGammaPtvsEffi->GetYaxis()->SetTitleOffset(1.5);
	
	hHistoTruePi0DalitzClusGammaPtvsEffi->SetLineColor(kRed);
	hHistoTruePi0DalitzClusGammaPtvsEffi->SetMarkerStyle(22);
	hHistoTruePi0DalitzClusGammaPtvsEffi->SetMarkerColor(kRed);
	hHistoTruePi0DalitzClusGammaPtvsEffi->Draw();
	hHistoTruePi0DalitzClusGammaPtvsEffi->Draw();
	
	
	canvasTruePi0DalitzClusGammaPtvsEffi->Update();
	canvasTruePi0DalitzClusGammaPtvsEffi->SaveAs(Form("%s/TruePi0DalitzClusGamma_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasTruePi0DalitzClusGammaPtvsEffi;		
	
	  
		  
	  
	}
	
	
	
	
	
	
	
	
	if( hMCAllPositronsPt && hESDTruePositronPt ){
	
	
	
	hMCAllPositronsPtScaled    =  (TH1F*) hMCAllPositronsPt->Clone("hMCAllPositronsPtScaled");
	hMCAllPositronsPtScaledBR  =  (TH1F*) hMCAllPositronsPt->Clone("hMCAllPositronsPtScaledBR");
	
	hMCAllPositronsPtScaled->Sumw2();
	hMCAllPositronsPtScaled->Scale(1.0/nEventsMC);
	
	hMCAllPositronsPtScaledBR->Sumw2();
	hMCAllPositronsPtScaledBR->Scale(1.0/nEventsMC);
	hMCAllPositronsPtScaledBR->Scale(fPi0DalitzBRDPG/Pi0DalitzBR);
	
	
	hESDTruePositronPtScaled 	=  (TH1F*) hESDTruePositronPt->Clone("hESDTruePositronPtScaled");
	hESDTruePositronPtScaledBR	=  (TH1F*) hESDTruePositronPt->Clone("hESDTruePositronPtScaledBR");
	
	hESDTruePositronPtScaled->Sumw2();
	hESDTruePositronPtScaled->Scale(1.0/nEventsMC);
	
	hESDTruePositronPtScaledBR->Sumw2();
	hESDTruePositronPtScaledBR->Scale(1.0/nEventsMC);
	hESDTruePositronPtScaledBR->Scale(fPi0DalitzBRDPG/Pi0DalitzBR);
	
	
	  
	  
	hESDTruePositronPtvsEffi = (TH1F*) hESDTruePositronPt->Rebin(kNbins-1,"hESDTruePositronPtvsEffi",xbins);
	
	
	hESDTruePositronPtvsEffi->Sumw2();
	
	
		
	TH1F* hMCAllPositronsPtRebin = (TH1F*) hMCAllPositronsPt->Rebin(kNbins-1,"hMCAllPositronsPtRebin",xbins);
	
	
	hMCAllPositronsPtRebin->Sumw2();
	
	hESDTruePositronPtvsEffi->Divide(hESDTruePositronPtvsEffi,hMCAllPositronsPtRebin,1,1,"B");
	  
	  
	TCanvas * canvasPositronEffivsPt = new TCanvas("canvasPositronEffivsPt","",10,10,600,470);  // gives the page size
	
	canvasPositronEffivsPt->SetLeftMargin(0.08);
	canvasPositronEffivsPt->SetRightMargin(0.03);
	canvasPositronEffivsPt->SetTopMargin(0.03);
	canvasPositronEffivsPt->SetBottomMargin(0.1);
	canvasPositronEffivsPt->SetLogx(1);
	
	
	canvasPositronEffivsPt->cd();
	
	
	//TLatex *latexPositronEffivsPt = new TLatex(0.45,0.81,electronCutEffiLegend);
        //SetStyleTLatex( latexPositronEffivsPt, 0.038,4);
	
	
	hESDTruePositronPtvsEffi->SetStats(kFALSE);
	
        hESDTruePositronPtvsEffi->SetXTitle("p_{T} (GeV/c)");
	hESDTruePositronPtvsEffi->SetYTitle("#epsilon");
	

	
	hESDTruePositronPtvsEffi->GetYaxis()->SetRangeUser(0.001,1.0);
	hESDTruePositronPtvsEffi->GetXaxis()->SetRangeUser(0.1,10.0);
	
	hESDTruePositronPtvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hESDTruePositronPtvsEffi->GetYaxis()->SetTitleOffset(1.0);
	
	
	hESDTruePositronPtvsEffi->SetLineColor(kBlue);
	hESDTruePositronPtvsEffi->SetMarkerStyle(29);
	hESDTruePositronPtvsEffi->SetMarkerColor(kBlue);
	hESDTruePositronPtvsEffi->DrawCopy();
	
	//latexPositronEffivsPt->Draw();
	
	
	canvasPositronEffivsPt->Update();
	canvasPositronEffivsPt->SaveAs(Form("%s/Positron_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronEffivsPt;		
	  
	    
	  
	}
	
		
	if( hMCDecayElectronPi0Pt && hESDTruePi0DalitzElectronPt ){
	
	
		 
	hESDTruePi0DalitzElectronPtvsEffi = (TH1F*) hESDTruePi0DalitzElectronPt->Rebin(kNbins-1,"hESDTruePi0DalitzElectronPtvsEffi",xbins);
	
	
	hESDTruePi0DalitzElectronPtvsEffi->Sumw2();
	
	
		
	TH1F* hMCDecayElectronPi0PtRebin = (TH1F*) hMCDecayElectronPi0Pt->Rebin(kNbins-1,"hMCDecayElectronPi0PtRebin",xbins);
	
	
	hMCDecayElectronPi0PtRebin->Sumw2();
	
	hESDTruePi0DalitzElectronPtvsEffi->Divide(hESDTruePi0DalitzElectronPtvsEffi,hMCDecayElectronPi0PtRebin,1,1,"B");
	  
	  
	TCanvas * canvasPi0DalitzElectronEffivsPt = new TCanvas("canvasPi0DalitzElectronEffivsPt","",10,10,600,470);  // gives the page size
	
	canvasPi0DalitzElectronEffivsPt->SetLeftMargin(0.12);
	canvasPi0DalitzElectronEffivsPt->SetRightMargin(0.03);
	canvasPi0DalitzElectronEffivsPt->SetTopMargin(0.03);
	canvasPi0DalitzElectronEffivsPt->SetBottomMargin(0.1);
	
	TLatex *latexPi0DalitzElectronEffivsPt = new TLatex(0.45,0.81,electronCutEffiLegend);
        SetStyleTLatex( latexPi0DalitzElectronEffivsPt, 0.038,4);
	
	
	canvasPi0DalitzElectronEffivsPt->cd();
	
	hESDTruePi0DalitzElectronPtvsEffi->SetStats(kFALSE);
	
        hESDTruePi0DalitzElectronPtvsEffi->SetXTitle("p_{T} (GeV/c)");
	hESDTruePi0DalitzElectronPtvsEffi->SetYTitle("Electron Efficiency");
	

	
	hESDTruePi0DalitzElectronPtvsEffi->GetYaxis()->SetRangeUser(0.0,1.0);
	hESDTruePi0DalitzElectronPtvsEffi->GetXaxis()->SetRangeUser(0.08,10.0);
	
	hESDTruePi0DalitzElectronPtvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hESDTruePi0DalitzElectronPtvsEffi->GetYaxis()->SetTitleOffset(1.5);
	
	
	hESDTruePi0DalitzElectronPtvsEffi->SetLineColor(kRed);
	hESDTruePi0DalitzElectronPtvsEffi->SetMarkerStyle(22);
	hESDTruePi0DalitzElectronPtvsEffi->SetMarkerColor(kRed);
	hESDTruePi0DalitzElectronPtvsEffi->DrawCopy();
	latexPi0DalitzElectronEffivsPt->Draw();
	
	
	canvasPi0DalitzElectronEffivsPt->Update();
	canvasPi0DalitzElectronEffivsPt->SaveAs(Form("%s/Pi0DalitzElectron_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPi0DalitzElectronEffivsPt;		
	  
	    
	  
	}
	
	
	if( hMCDecayPositronPi0Pt && hESDTruePi0DalitzPositronPt ){
	
	
		 
	hESDTruePi0DalitzPositronPtvsEffi = (TH1F*) hESDTruePi0DalitzPositronPt->Rebin(kNbins-1,"hESDTruePi0DalitzPositronPtvsEffi",xbins);
	
	
	hESDTruePi0DalitzPositronPtvsEffi->Sumw2();
	
	
		
	TH1F* hMCDecayPositronPi0PtRebin = (TH1F*) hMCDecayPositronPi0Pt->Rebin(kNbins-1,"hMCDecayPositronPi0PtRebin",xbins);
	
	
	hMCDecayPositronPi0PtRebin->Sumw2();
	
	hESDTruePi0DalitzPositronPtvsEffi->Divide(hESDTruePi0DalitzPositronPtvsEffi,hMCDecayPositronPi0PtRebin,1,1,"B");
	  
	  
	TCanvas * canvasPi0DalitzPositronEffivsPt = new TCanvas("canvasPi0DalitzPositronEffivsPt","",10,10,600,470);  // gives the page size
	
	canvasPi0DalitzPositronEffivsPt->SetLeftMargin(0.12);
	canvasPi0DalitzPositronEffivsPt->SetRightMargin(0.03);
	canvasPi0DalitzPositronEffivsPt->SetTopMargin(0.03);
	canvasPi0DalitzPositronEffivsPt->SetBottomMargin(0.1);
	
	TLatex *latexPi0DalitzPositronEffivsPt = new TLatex(0.45,0.81,electronCutEffiLegend);
        SetStyleTLatex( latexPi0DalitzPositronEffivsPt, 0.038,4);
	
	
	canvasPi0DalitzPositronEffivsPt->cd();
	
	hESDTruePi0DalitzPositronPtvsEffi->SetStats(kFALSE);
	
        hESDTruePi0DalitzPositronPtvsEffi->SetXTitle("p_{T} (GeV/c)");
	hESDTruePi0DalitzPositronPtvsEffi->SetYTitle("Positron Efficiency");
	

	
	hESDTruePi0DalitzPositronPtvsEffi->GetYaxis()->SetRangeUser(0.0,1.0);
	hESDTruePi0DalitzPositronPtvsEffi->GetXaxis()->SetRangeUser(0.08,10.0);
	
	hESDTruePi0DalitzPositronPtvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hESDTruePi0DalitzPositronPtvsEffi->GetYaxis()->SetTitleOffset(1.5);
	
	
	hESDTruePi0DalitzPositronPtvsEffi->SetLineColor(kRed);
	hESDTruePi0DalitzPositronPtvsEffi->SetMarkerStyle(22);
	hESDTruePi0DalitzPositronPtvsEffi->SetMarkerColor(kRed);
	hESDTruePi0DalitzPositronPtvsEffi->DrawCopy();
	latexPi0DalitzPositronEffivsPt->Draw();
	
	
	canvasPi0DalitzPositronEffivsPt->Update();
	canvasPi0DalitzPositronEffivsPt->SaveAs(Form("%s/Pi0DalitzPositron_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPi0DalitzPositronEffivsPt;		
	  
	    
	  
	}
	
	
	if( hMCAllElectronsPt && hESDTrueElectronPt ){
	
	
	  
	hMCAllElectronsPtScaled  	=  (TH1F*) hMCAllElectronsPt->Clone("hMCAllElectronsPtScaled");
	hMCAllElectronsPtScaledBR  	=  (TH1F*) hMCAllElectronsPt->Clone("hMCAllElectronsPtScaledBR");
	
	hMCAllElectronsPtScaled->Sumw2();
	hMCAllElectronsPtScaled->Scale(1.0/nEventsMC);  
	
	hMCAllElectronsPtScaledBR->Sumw2();
	hMCAllElectronsPtScaledBR->Scale(1.0/nEventsMC);  
	hMCAllElectronsPtScaledBR->Scale(fPi0DalitzBRDPG/Pi0DalitzBR);
	
	
	
	
	hESDTrueElectronPtScaled 	=  (TH1F*) hESDTrueElectronPt->Clone("hESDTrueElectronPtScaled");
	hESDTrueElectronPtScaledBR  	=  (TH1F*) hESDTrueElectronPt->Clone("hESDTrueElectronPtScaledBR");
	
	
	hESDTrueElectronPtScaled->Sumw2();
	hESDTrueElectronPtScaled->Scale(1.0/nEventsMC);
	
	hESDTrueElectronPtScaledBR->Sumw2();
	hESDTrueElectronPtScaledBR->Scale(1.0/nEventsMC);
	hESDTrueElectronPtScaledBR->Scale(fPi0DalitzBRDPG/Pi0DalitzBR);
	
	
	
	
  
	hESDTrueElectronPtvsEffi = (TH1F*) hESDTrueElectronPt->Rebin(kNbins-1,"hESDTrueElectronPtvsEffi",xbins);
	
	
	hESDTrueElectronPtvsEffi->Sumw2();
	
	
		
	TH1F* hMCAllElectronsPtRebin = (TH1F*) hMCAllElectronsPt->Rebin(kNbins-1,"hMCAllElectronsPtRebin",xbins);
	
	
	hMCAllElectronsPtRebin->Sumw2();
	
	hESDTrueElectronPtvsEffi->Divide(hESDTrueElectronPtvsEffi,hMCAllElectronsPtRebin,1,1,"B");
	  
	  
	TCanvas * canvasElectronEffivsPt = new TCanvas("canvasElectronEffivsPt","",10,10,600,470);  // gives the page size
	
	canvasElectronEffivsPt->SetLeftMargin(0.08);
	canvasElectronEffivsPt->SetRightMargin(0.03);
	canvasElectronEffivsPt->SetTopMargin(0.03);
	canvasElectronEffivsPt->SetBottomMargin(0.1);
	
	//canvasElectronEffivsPt->SetLogy(1);
	canvasElectronEffivsPt->SetLogx(1);
	
	//TLatex *latexElectronEffivsPt = new TLatex(0.45,0.81,electronCutEffiLegend);
        //SetStyleTLatex( latexElectronEffivsPt, 0.038,4);
	
	
	canvasElectronEffivsPt->cd();
	
	hESDTrueElectronPtvsEffi->SetStats(kFALSE);
	
        hESDTrueElectronPtvsEffi->SetXTitle("p_{T} (GeV/c)");
	hESDTrueElectronPtvsEffi->SetYTitle("#epsilon^{e^{-}}");
	

	
	hESDTrueElectronPtvsEffi->GetYaxis()->SetRangeUser(0.001,1.0);
	hESDTrueElectronPtvsEffi->GetXaxis()->SetRangeUser(0.1,10.0);
	
	hESDTrueElectronPtvsEffi->GetXaxis()->SetTitleOffset(1.3);
	hESDTrueElectronPtvsEffi->GetYaxis()->SetTitleOffset(1.0);
	
	
	hESDTrueElectronPtvsEffi->SetLineColor(kRed);
	hESDTrueElectronPtvsEffi->SetMarkerStyle(22);
	hESDTrueElectronPtvsEffi->SetMarkerColor(kRed);
	hESDTrueElectronPtvsEffi->DrawCopy();
	//latexElectronEffivsPt->Draw();
	
	
	canvasElectronEffivsPt->Update();
	canvasElectronEffivsPt->SaveAs(Form("%s/Electron_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronEffivsPt;		
	  
	    
	  
	}
	
	
	if(  hESDTrueElectronPtvsEffi  &&  hESDTruePositronPtvsEffi  ){
	
	TCanvas * canvasElectronPositronEffivsPt = new TCanvas("canvasElectronPositronEffivsPt","",10,10,600,470);  // gives the page size
	
	canvasElectronPositronEffivsPt->SetLeftMargin(0.12);
	canvasElectronPositronEffivsPt->SetRightMargin(0.03);
	canvasElectronPositronEffivsPt->SetTopMargin(0.03);
	canvasElectronPositronEffivsPt->SetBottomMargin(0.1);
	
	//TLatex *latexElectronPositronEffivsPt = new TLatex(0.45,0.81,electronCutEffiLegend);
        //SetStyleTLatex( latexElectronPositronEffivsPt, 0.038,4);
	
	canvasElectronPositronEffivsPt->cd();
	
	hESDTrueElectronPtvsEffi->SetMarkerColor(kRed);
	hESDTrueElectronPtvsEffi->SetMarkerStyle(26);
	hESDTruePositronPtvsEffi->SetMarkerColor(kBlue);
	hESDTruePositronPtvsEffi->SetLineColor(kBlue);
	hESDTruePositronPtvsEffi->SetMarkerStyle(29);
	
	
	hESDTrueElectronPtvsEffi->DrawCopy();
	hESDTruePositronPtvsEffi->DrawCopy("sames");
	//latexElectronPositronEffivsPt->Draw();
	
	
	
	canvasElectronPositronEffivsPt->Update();
	canvasElectronPositronEffivsPt->SaveAs(Form("%s/ElectronPositron_Effi_Vs_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronEffivsPt;		

	    
	  
	}
	
    
	if( hESDDalitzElectronAfterPt_data && hESDDalitzElectronAfterPt_mc ){
	
	TCanvas * canvasElectronPt = new TCanvas("canvasElectronPt","",10,10,500,500);  // gives the page size
	canvasElectronPt->SetLeftMargin(0.16);
	canvasElectronPt->SetBottomMargin(0.1);
	canvasElectronPt->cd();
	
	
	TH1F* hESDDalitzElectronAfterPtRebin_data  = (TH1F*) hESDDalitzElectronAfterPt_data->Clone("hESDDalitzElectronAfterPtRebin_data");
	TH1F* hESDDalitzElectronAfterPtRebin_mc  = (TH1F*) hESDDalitzElectronAfterPt_mc->Clone("hESDDalitzElectronAfterPtRebin_mc");
	
	hESDDalitzElectronAfterPtRebin_data->SetStats(kFALSE);
        hESDDalitzElectronAfterPtRebin_mc->SetStats(kFALSE);
	
	hESDDalitzElectronAfterPtRebin_data->Scale(1./hESDDalitzElectronAfterPtRebin_data->Integral());

	
	
	
	hESDDalitzElectronAfterPtRebin_mc->Scale(1./hESDDalitzElectronAfterPtRebin_mc->Integral());
	
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterPtRebin_data,
						hESDDalitzElectronAfterPtRebin_mc, 
						"",pTLabelN.Data(),textYAxisPtDistE,
						kTRUE, 5.,1e-3,
						kFALSE,1e-5 ,3e-1, 
						kTRUE, 0.1,10.,1.2);	
	
	
	
	processSystem->Draw();
	
	canvasElectronPt->SetLogy(1);
	canvasElectronPt->SetLogx(1);
	canvasElectronPt->Update();
	
	canvasElectronPt->SaveAs(Form("%s/Electron_After_MassCut_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPt;		
	
		
	}
	
	if( hESDDalitzPositronAfterPt_data && hESDDalitzPositronAfterPt_mc ){
	
	TCanvas * canvasPositronPt = new TCanvas("canvasPositronPt","",10,10,500,500);  // gives the page size
	canvasPositronPt->SetLeftMargin(0.16);
	canvasPositronPt->SetBottomMargin(0.1);
	canvasPositronPt->cd();
	
	
	TH1F* hESDDalitzPositronAfterPtRebin_data  = (TH1F*) hESDDalitzPositronAfterPt_data->Clone("hESDDalitzPositronAfterPtRebin_data");
	TH1F* hESDDalitzPositronAfterPtRebin_mc  = (TH1F*) hESDDalitzPositronAfterPt_mc->Clone("hESDDalitzPositronAfterPtRebin_mc");
	
	hESDDalitzPositronAfterPtRebin_data->SetStats(kFALSE);
        hESDDalitzPositronAfterPtRebin_mc->SetStats(kFALSE);
	
	hESDDalitzPositronAfterPtRebin_data->Scale(1./hESDDalitzPositronAfterPtRebin_data->Integral());

	
	
	
	hESDDalitzPositronAfterPtRebin_mc->Scale(1./hESDDalitzPositronAfterPtRebin_mc->Integral());
	
	
	
	
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterPtRebin_data,
						hESDDalitzPositronAfterPtRebin_mc, 
						"",pTLabelP.Data(),textYAxisPtDistP,
						kTRUE,5.,1e-3,
						kFALSE,1e-5 ,6e-1, 
						kTRUE, 0.1,10.,1.2);	
	
	//TLatex *process01 = new TLatex(0.63, 0.92,Form(collisionSystem.Data()));
        //process01->SetNDC();
        //process01->SetTextColor(1);
        //process01->SetTextSize(0.04);
	//process01->SetTextFont(42);
	//process01->Draw();
	processSystem->Draw();
	
	canvasPositronPt->SetLogy(1);
	canvasPositronPt->SetLogx(1);
	canvasPositronPt->Update();
	
	canvasPositronPt->SaveAs(Form("%s/Positron_After_MassCut_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronPt;
	}
	
	if( hESDDalitzPositronAfterPt_data && hESDDalitzPositronAfterPt_mc && hESDDalitzElectronAfterPt_data && hESDDalitzElectronAfterPt_mc ){
	  
	TCanvas * canvasElectronPositronPt = new TCanvas("canvasElectronPositronPt","",10,10,500,500);  // gives the page size
	canvasElectronPositronPt->SetLeftMargin(0.17);
	canvasElectronPositronPt->SetBottomMargin(0.13);
	canvasElectronPositronPt->cd();
	
	
	TH1F* hESDDalitzElectronPositronAfterPt_data = (TH1F*)hESDDalitzPositronAfterPt_data->Clone("ESD_DalitzElectronPositron_After_Pt_data");
	hESDDalitzElectronPositronAfterPt_data->Add(hESDDalitzElectronAfterPt_data,1.0);
	
	TH1F* hESDDalitzElectronPositronAfterPt_mc = (TH1F*)hESDDalitzPositronAfterPt_mc->Clone("ESD_DalitzElectronPositron_After_Pt_mc");
	hESDDalitzElectronPositronAfterPt_mc->Add(hESDDalitzElectronAfterPt_mc,1.0);
	
	
	hESDDalitzElectronPositronAfterPt_data->SetStats(kFALSE);
        hESDDalitzElectronPositronAfterPt_mc->SetStats(kFALSE);
	
	
	
	
	
	hESDDalitzElectronPositronAfterPt_data->Scale(1./hESDDalitzElectronPositronAfterPt_data->Integral());
	
	hESDDalitzElectronPositronAfterPt_mc->Scale(1./hESDDalitzElectronPositronAfterPt_mc->Integral());
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronPositronAfterPt_data,
						hESDDalitzElectronPositronAfterPt_mc, 
						"", "p_{T} e^{#pm} (GeV/c) ",textYAxisPtDistEP,
						kFALSE, 5.,1e-12,
						kTRUE,1e-5 ,3e-1, 
						kTRUE, 0.1,20.,1.3);	
	canvasElectronPositronPt->SetLogy(1);
	canvasElectronPositronPt->SetLogx(1);
	canvasElectronPositronPt->Update();
	canvasElectronPositronPt->SaveAs(Form("%s/ElectronPositron_After_MassCut_Pt%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronPt;		
	  
	  	  	  
	}
	
	
	
	
	if( hESDDalitzElectronAfterEta_data && hESDDalitzElectronAfterEta_mc ) {
	  
	  
	 
	
	TCanvas * canvasElectronEta = new TCanvas("canvasElectronEta","",10,10,500,500);  // gives the page size
	canvasElectronEta->SetLogy(1);
	canvasElectronEta->SetLeftMargin(0.16);
	canvasElectronEta->cd();
	
	TH1F* hESDDalitzElectronAfterEtaRebin_data = (TH1F*)hESDDalitzElectronAfterEta_data->Clone("hESDDalitzElectronAfterEtaRebin_data");
	TH1F* hESDDalitzElectronAfterEtaRebin_mc   = (TH1F*)hESDDalitzElectronAfterEta_mc->Clone("hESDDalitzElectronAfterEtaRebin_data");
	
	
	hESDDalitzElectronAfterEtaRebin_data->SetStats(kFALSE);
        hESDDalitzElectronAfterEtaRebin_mc->SetStats(kFALSE);
	
	hESDDalitzElectronAfterEtaRebin_data->Rebin(rebinEta);
	hESDDalitzElectronAfterEtaRebin_mc->Rebin(rebinEta);
	
	hESDDalitzElectronAfterEtaRebin_data->Scale(1./hESDDalitzElectronAfterEtaRebin_data->GetEntries());

	
	
	
	hESDDalitzElectronAfterEtaRebin_mc->Scale(1./hESDDalitzElectronAfterEtaRebin_mc->GetEntries());
	
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterEtaRebin_data,
						hESDDalitzElectronAfterEtaRebin_mc, 
						"", "#eta",textYAxisEtaDistE,
						kTRUE, 5.,1e-4,
						kFALSE,YMinEta ,YMaxEta, 
						kTRUE, -1.5,1.5);	
	
	processSystem->Draw();
	
	canvasElectronEta->Update();
	canvasElectronEta->SaveAs(Form("%s/Electron_After_MassCut_Eta%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronEta;		
	
	
	}
	
	if( hESDDalitzElectronAfterEtaPCut_data != NULL && hESDDalitzElectronAfterEtaPCut_mc != NULL ) {
	  
	TCanvas * canvasElectronEtaPCut = new TCanvas("canvasElectronEtaPCut","",10,10,500,500);  // gives the page size
	canvasElectronEtaPCut->SetLogy(1);
	canvasElectronEtaPCut->SetLeftMargin(0.16);
	canvasElectronEtaPCut->cd();
	
	hESDDalitzElectronAfterEtaPCut_data->SetStats(kFALSE);
        hESDDalitzElectronAfterEtaPCut_mc->SetStats(kFALSE);
	//hESDDalitzElectronAfterEtaPCut_data->Sumw2();
	//hESDDalitzElectronAfterEtaPCut_mc->Sumw2();
	
	hESDDalitzElectronAfterEtaPCut_data->Rebin(rebinEta);
	hESDDalitzElectronAfterEtaPCut_mc->Rebin(rebinEta);
	
	hESDDalitzElectronAfterEtaPCut_data->Scale(1./hESDDalitzElectronAfterEtaPCut_data->GetEntries());

	
	
	
	hESDDalitzElectronAfterEtaPCut_mc->Scale(1./hESDDalitzElectronAfterEtaPCut_mc->GetEntries());
	
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterEtaPCut_data,
						hESDDalitzElectronAfterEtaPCut_mc, 
						"", "#eta",textYAxisEtaDistE,
						kTRUE, 1.2,1e-3,
						kFALSE,0. ,0., 
						kTRUE, -60.,60.);	
	canvasElectronEtaPCut->Update();
	canvasElectronEtaPCut->SaveAs(Form("%s/Electron_After_MassCut_Eta_PCut%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronEtaPCut;		
	delete hESDDalitzElectronAfterEtaPCut_data;
	delete hESDDalitzElectronAfterEtaPCut_mc;
	
	}
	
	if( hESDDalitzPositronAfterEta_data && hESDDalitzPositronAfterEta_mc ) {
	TCanvas * canvasPositronEta = new TCanvas("canvasPositronEta","",10,10,500,500);  // gives the page size
	canvasPositronEta->SetLogy(1);
	canvasPositronEta->SetLeftMargin(0.16);
	canvasPositronEta->cd();
	
	TH1F* hESDDalitzPositronAfterEtaRebin_data = (TH1F*)hESDDalitzPositronAfterEta_data->Clone("hESDDalitzPositronAfterEtaRebin_data");
	TH1F* hESDDalitzPositronAfterEtaRebin_mc   = (TH1F*)hESDDalitzPositronAfterEta_mc->Clone("hESDDalitzPositronAfterEtaRebin_data");
	
	
	hESDDalitzPositronAfterEtaRebin_data->SetStats(kFALSE);
        hESDDalitzPositronAfterEtaRebin_mc->SetStats(kFALSE);
	//hESDDalitzPositronAfterEtaRebin_data->Sumw2();
	//hESDDalitzPositronAfterEtaRebin_mc->Sumw2();
	
	hESDDalitzPositronAfterEtaRebin_data->Rebin(rebinEta);
	hESDDalitzPositronAfterEtaRebin_mc->Rebin(rebinEta);
	
	hESDDalitzPositronAfterEtaRebin_data->Scale(1./hESDDalitzPositronAfterEtaRebin_data->GetEntries());
        hESDDalitzPositronAfterEtaRebin_mc->Scale(1./hESDDalitzPositronAfterEtaRebin_mc->GetEntries());
	
	
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterEtaRebin_data,
						hESDDalitzPositronAfterEtaRebin_mc, 
						"", "#eta",textYAxisEtaDistP,
						kTRUE, 5.,1e-4,
						kFALSE,YMinEta ,YMaxEta, 
						kTRUE, -1.5,1.5);	
	
	processSystem->Draw();
	
	canvasPositronEta->Update();
	canvasPositronEta->SaveAs(Form("%s/Positron_After_MassCut_Eta%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronEta;
	
	}
		
	if( hESDDalitzPositronAfterEtaPCut_data != NULL && hESDDalitzPositronAfterEtaPCut_mc != NULL ) {
	  
	TCanvas * canvasPositronEtaPCut = new TCanvas("canvasPositronEtaPCut","",10,10,500,500);  // gives the page size
	canvasPositronEtaPCut->SetLogy(1);
	canvasPositronEtaPCut->SetLeftMargin(0.16);
	canvasPositronEtaPCut->cd();
	
	hESDDalitzPositronAfterEtaPCut_data->SetStats(kFALSE);
        hESDDalitzPositronAfterEtaPCut_mc->SetStats(kFALSE);
	//hESDDalitzPositronAfterEtaPCut_data->Sumw2();
	//hESDDalitzPositronAfterEtaPCut_mc->Sumw2();
	
	hESDDalitzPositronAfterEtaPCut_data->Rebin(rebinEta);
	hESDDalitzPositronAfterEtaPCut_mc->Rebin(rebinEta);
	
	hESDDalitzPositronAfterEtaPCut_data->Scale(1./hESDDalitzPositronAfterEtaPCut_data->GetEntries());

	
	
	
	hESDDalitzPositronAfterEtaPCut_mc->Scale(1./hESDDalitzPositronAfterEtaPCut_mc->GetEntries());
	
	
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterEtaPCut_data,
						hESDDalitzPositronAfterEtaPCut_mc, 
						"", "#eta",textYAxisEtaDistP,
						kTRUE, 1.2,1e-3,
						kFALSE,0. ,0., 
						kTRUE, -60.,60.);	
	canvasPositronEtaPCut->Update();
	canvasPositronEtaPCut->SaveAs(Form("%s/Positron_After_MassCut_Eta_PCut%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronEtaPCut;		
	delete hESDDalitzPositronAfterEtaPCut_data;
	delete hESDDalitzPositronAfterEtaPCut_mc;
	
	}
	
	
	if( hESDDalitzElectronAfterEta_data && hESDDalitzElectronAfterEta_mc && hESDDalitzPositronAfterEta_data && hESDDalitzPositronAfterEta_mc ) {
	  
	  
	 
	
	TCanvas * canvasElectronEta = new TCanvas("canvasElectronEta","",10,10,500,500);  // gives the page size
	canvasElectronEta->SetLogy(1);
	canvasElectronEta->SetLeftMargin(0.16);
	canvasElectronEta->cd();
	
	TH1F* hESDDalitzElectronPositronAfterEta_data = (TH1F*)hESDDalitzElectronAfterEta_data->Clone("hESDDalitzElectronPositronAfterEta_data");
	hESDDalitzElectronPositronAfterEta_data->Add(hESDDalitzPositronAfterEta_data,1.0);
	
	TH1F* hESDDalitzElectronPositronAfterEta_mc = (TH1F*)hESDDalitzElectronAfterEta_mc->Clone("hESDDalitzElectronPositronAfterEta_mc");
	hESDDalitzElectronPositronAfterEta_mc->Add(hESDDalitzPositronAfterEta_mc,1.0);
	
	
	hESDDalitzElectronPositronAfterEta_data->SetStats(kFALSE);
        hESDDalitzElectronPositronAfterEta_mc->SetStats(kFALSE);
	
	//hESDDalitzElectronPositronAfterEta_data->Sumw2();
	//hESDDalitzElectronPositronAfterEta_mc->Sumw2();
	
	hESDDalitzElectronPositronAfterEta_data->Rebin(rebinEta);
	hESDDalitzElectronPositronAfterEta_mc->Rebin(rebinEta);
	
	hESDDalitzElectronPositronAfterEta_data->Scale(1./hESDDalitzElectronPositronAfterEta_data->Integral());

	
	
	
	hESDDalitzElectronPositronAfterEta_mc->Scale(1./hESDDalitzElectronPositronAfterEta_mc->Integral());
	
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronPositronAfterEta_data,
						hESDDalitzElectronPositronAfterEta_mc, 
						"", "#eta",textYAxisEtaDistEP,
						kFALSE, 1.,1e-3,
						kTRUE,YMinEta ,YMaxEta, 
						kTRUE, -60.,60.);	
	canvasElectronEta->Update();
	canvasElectronEta->SaveAs(Form("%s/ElectronPositron_After_MassCut_Eta%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronEta;		
	
	
	}
	
	
	
	
	
	if( hESDConvGammaEta_data && hESDConvGammaEta_mc ) {
	
	TCanvas * canvasGammaEta = new TCanvas("canvasGammaEta","",10,10,500,500);  // gives the page size
	canvasGammaEta->SetLogy(1);
	canvasGammaEta->SetLeftMargin(0.14);
	canvasGammaEta->cd();
	
	hESDConvGammaEta_data->SetStats(kFALSE);
        hESDConvGammaEta_mc->SetStats(kFALSE);
	//hESDConvGammaEta_data->Sumw2();
	//hESDConvGammaEta_mc->Sumw2();
	
	hESDConvGammaEta_data->Scale(1./hESDConvGammaEta_data->GetEntries());

	
	
	
	hESDConvGammaEta_mc->Scale(1./hESDConvGammaEta_mc->GetEntries());
	
	
	
	DrawAutoGammaHistosTemp( hESDConvGammaEta_data,
						hESDConvGammaEta_mc, 
						"", "#eta",textYAxisEtaDistGamma,
						kTRUE, 5.,1e-3,
						kFALSE,0. ,0., 
						kTRUE, -60.,60.);	
	canvasGammaEta->Update();
	canvasGammaEta->SaveAs(Form("%s/Gamma_After_MassCut_Eta%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasGammaEta;		
	delete hESDConvGammaEta_data;
	delete hESDConvGammaEta_mc;
	
	}
	
	if( hESDDalitzElectronAfterPhi_data && hESDDalitzElectronAfterPhi_mc ) {
	
	  TCanvas * canvasElectronPhi = new TCanvas("canvasElectronPhi","",10,10,500,500);  // gives the page size
	  canvasElectronPhi->SetLeftMargin(0.14);
	  canvasElectronPhi->cd();
      	  canvasElectronPhi->SetLogy(1);

	  TH1F* hESDDalitzElectronAfterPhiRebin_data = (TH1F*)hESDDalitzElectronAfterPhi_data->Clone("hESDDalitzElectronAfterPhiRebin_data");
	  TH1F* hESDDalitzElectronAfterPhiRebin_mc   = (TH1F*)hESDDalitzElectronAfterPhi_mc->Clone("hESDDalitzElectronAfterPhiRebin_mc");
	  
	  hESDDalitzElectronAfterPhiRebin_data->SetStats(kFALSE);
	  hESDDalitzElectronAfterPhiRebin_mc->SetStats(kFALSE);
	  
	  hESDDalitzElectronAfterPhiRebin_data->Rebin(rebinPhi);
	  hESDDalitzElectronAfterPhiRebin_mc->Rebin(rebinPhi);
      
	  hESDDalitzElectronAfterPhiRebin_data->Scale(1.0/hESDDalitzElectronAfterPhiRebin_data->Integral());
	  hESDDalitzElectronAfterPhiRebin_mc->Scale(1.0/hESDDalitzElectronAfterPhiRebin_mc->Integral());
     
      
	  DrawAutoGammaHistosTemp( hESDDalitzElectronAfterPhiRebin_data,
						hESDDalitzElectronAfterPhiRebin_mc, 
						"", "#phi (rad) ",textYAxisPhiE,
						kTRUE,5.0,1e-4,
						kFALSE,1e-5 ,3e-2, 
						kTRUE, 0.,2*TMath::Pi());
	  
	  processSystem->Draw();
	  canvasElectronPhi->Update();
						
	  canvasElectronPhi->SaveAs(Form("%s/Electron_After_MassCut_Phi%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	  delete canvasElectronPhi;		
	  
	
	}

	if( hESDDalitzPositronAfterPhi_data && hESDDalitzPositronAfterPhi_mc ){
	  
	  
	  TCanvas * canvasPositronPhi = new TCanvas("canvasPositronPhi","",10,10,500,500);  // gives the page size
	  canvasPositronPhi->SetLeftMargin(0.14);
	  canvasPositronPhi->cd(1);
      	  canvasPositronPhi->SetLogy(1);
	  
	  
	

	  TH1F* hESDDalitzPositronAfterPhiRebin_data = (TH1F*)hESDDalitzPositronAfterPhi_data->Clone("hESDDalitzPositronAfterPhiRebin_data");
	  TH1F* hESDDalitzPositronAfterPhiRebin_mc   = (TH1F*)hESDDalitzPositronAfterPhi_mc->Clone("hESDDalitzPositronAfterPhiRebin_mc");
	  
	  hESDDalitzPositronAfterPhiRebin_data->SetStats(kFALSE);
	  hESDDalitzPositronAfterPhiRebin_mc->SetStats(kFALSE);
	
	  hESDDalitzPositronAfterPhiRebin_data->Rebin(rebinPhi);
	  hESDDalitzPositronAfterPhiRebin_mc->Rebin(rebinPhi);
      
	  hESDDalitzPositronAfterPhiRebin_data->Scale(1.0/hESDDalitzPositronAfterPhiRebin_data->Integral());
	  hESDDalitzPositronAfterPhiRebin_mc->Scale(1.0/hESDDalitzPositronAfterPhiRebin_mc->Integral());
     
      
	  DrawAutoGammaHistosTemp( hESDDalitzPositronAfterPhiRebin_data,
						hESDDalitzPositronAfterPhiRebin_mc, 
						"", "#phi (rad) ",textYAxisPhiP,
						kTRUE,5.0,1e-4,
						kFALSE,1e-5 ,3e-2, 
						kTRUE, 0.,2*TMath::Pi());
	  
	  processSystem->Draw();
	  canvasPositronPhi->Update();
						
	  canvasPositronPhi->SaveAs(Form("%s/Positron_After_MassCut_Phi%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	  delete canvasPositronPhi;		
	  
	
	}
	
	if( hESDDalitzPositronAfterPhi_data && hESDDalitzPositronAfterPhi_mc  && hESDDalitzElectronAfterPhi_data && hESDDalitzElectronAfterPhi_mc){
	  
	TCanvas * canvasElectronPositronPhi = new TCanvas("canvasElectronPositronPhi","",10,10,500,500);  // gives the page size
	canvasElectronPositronPhi->cd(1);
      	canvasElectronPositronPhi->SetLogy(1);
	
	canvasElectronPositronPhi->SetLeftMargin(0.14);
	canvasElectronPositronPhi->SetRightMargin(0.02);
	
	
	

        TH1F* hESDDalitzElectronPositronAfterPhi_data = (TH1F*) hESDDalitzPositronAfterPhi_data->Clone("ESD_DalitzElectronPositron_After_Phi_data");
	hESDDalitzElectronPositronAfterPhi_data->Add(hESDDalitzElectronAfterPhi_data,1.0);
	
	TH1F* hESDDalitzElectronPositronAfterPhi_mc = (TH1F*) hESDDalitzPositronAfterPhi_mc->Clone("ESD_DalitzElectronPositron_After_Phi_mc");
	hESDDalitzElectronPositronAfterPhi_mc->Add(hESDDalitzElectronAfterPhi_mc,1.0);
	
	
	hESDDalitzElectronPositronAfterPhi_data->SetStats(kFALSE);
	hESDDalitzElectronPositronAfterPhi_mc->SetStats(kFALSE);
	
	hESDDalitzElectronPositronAfterPhi_data->Rebin(rebinPhi);
	hESDDalitzElectronPositronAfterPhi_mc->Rebin(rebinPhi);
      
	hESDDalitzElectronPositronAfterPhi_data->Scale(1.0/hESDDalitzElectronPositronAfterPhi_data->Integral());
	hESDDalitzElectronPositronAfterPhi_mc->Scale(1.0/hESDDalitzElectronPositronAfterPhi_mc->Integral());
       
      
	DrawAutoGammaHistosTemp( hESDDalitzElectronPositronAfterPhi_data,
						hESDDalitzElectronPositronAfterPhi_mc, 
						"", "#phi (rad)",textYAxisPhiEP,
						kFALSE,1.6,1e-4,
						kTRUE,1e-5 ,2e-1, 
						kTRUE, 0.,2*TMath::Pi());	
	canvasElectronPositronPhi->Update();
	canvasElectronPositronPhi->SaveAs(Form("%s/ElectronPositron_After_MassCut_Phi%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronPhi;		
	
	
	}
	
		
	
	
	if( hESDDalitzElectronAfterNFindClsTPC_data ) {
	  
	TCanvas * canvasElectronNFindCls_data = new TCanvas("canvasElectronNFindCls_data","",10,10,500,500);  // gives the page size
	canvasElectronNFindCls_data->SetLeftMargin(0.13);
	canvasElectronNFindCls_data->SetRightMargin(0.12);
	canvasElectronNFindCls_data->SetTopMargin(0.05);
	canvasElectronNFindCls_data->SetBottomMargin(0.1);
	
	
	
	
	canvasElectronNFindCls_data->cd();
	canvasElectronNFindCls_data->SetLogz(1);
	
	hESDDalitzElectronAfterNFindClsTPC_data->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterNFindClsTPC_data,
						 "", "p (GeV/c)", "d#it{E}/d#it{x}","Data",
						 kTRUE,0., 15.,
						 kTRUE, 0.,1.5);
	
	
	canvasElectronNFindCls_data->Update();
	canvasElectronNFindCls_data->SaveAs(Form("%s/Electron_After_MassCut_NFindCls_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNFindCls_data;
	
	}
	
	if( hESDDalitzElectronAfterNFindClsTPC_mc ) {
	
	TCanvas * canvasElectronNFindCls_mc = new TCanvas("canvasElectronNFindCls_mc","",10,10,500,500);  // gives the page size
	canvasElectronNFindCls_mc->SetLeftMargin(0.13);
	canvasElectronNFindCls_mc->SetRightMargin(0.12);
	canvasElectronNFindCls_mc->SetTopMargin(0.05);
	canvasElectronNFindCls_mc->SetBottomMargin(0.1);
	
	
	
	
	canvasElectronNFindCls_mc->cd();
	canvasElectronNFindCls_mc->SetLogz(1);
	
	hESDDalitzElectronAfterNFindClsTPC_mc->SetStats(kFALSE);
       
	
	DrawAutoGammaHisto2D(	hESDDalitzElectronAfterNFindClsTPC_mc,
						 "", "p (GeV/c)", "d#it{E}/d#it{x}","MC",
						 kTRUE,0., 15.,
						 kTRUE, 0.,1.5);
	
	
	canvasElectronNFindCls_mc->Update();
	canvasElectronNFindCls_mc->SaveAs(Form("%s/Electron_After_MassCut_NFindCls_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNFindCls_mc;	
	
	}
	
	if( hESDDalitzPositronAfterNFindClsTPC_data ){
	  
	TCanvas * canvasPositronNFindCls_data = new TCanvas("canvasPositronNFindCls_data","",10,10,500,500);  // gives the page size
	canvasPositronNFindCls_data->SetLeftMargin(0.13);
	canvasPositronNFindCls_data->SetRightMargin(0.12);
	canvasPositronNFindCls_data->SetTopMargin(0.05);
	canvasPositronNFindCls_data->SetBottomMargin(0.1);
	
	
	
	
	canvasPositronNFindCls_data->cd();
	canvasPositronNFindCls_data->SetLogz(1);
	
	hESDDalitzPositronAfterNFindClsTPC_data->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterNFindClsTPC_data,
						 "", "p (GeV/c)", "d#it{E}/d#it{x}","Data",
						 kTRUE,0., 15.,
						 kTRUE, 0.,1.5);
	
	
	canvasPositronNFindCls_data->Update();
	canvasPositronNFindCls_data->SaveAs(Form("%s/Positron_After_MassCut_NFindCls_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNFindCls_data;
	}
	
	if( hESDDalitzPositronAfterNFindClsTPC_mc ) {
	
	TCanvas * canvasPositronNFindCls_mc = new TCanvas("canvasPositronNFindCls_mc","",10,10,500,500);  // gives the page size
	canvasPositronNFindCls_mc->SetLeftMargin(0.13);
	canvasPositronNFindCls_mc->SetRightMargin(0.12);
	canvasPositronNFindCls_mc->SetTopMargin(0.05);
	canvasPositronNFindCls_mc->SetBottomMargin(0.1);
	canvasPositronNFindCls_mc->cd();
	canvasPositronNFindCls_mc->SetLogz(1);
	
	hESDDalitzPositronAfterNFindClsTPC_mc->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterNFindClsTPC_mc,
						 "", "p (GeV/c)", "d#it{E}/d#it{x}","MC",
						 kTRUE,0., 15.,
						 kTRUE, 0.,1.5);
	
	
	canvasPositronNFindCls_mc->Update();
	canvasPositronNFindCls_mc->SaveAs(Form("%s/Positron_After_MassCut_NFindCls_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNFindCls_mc;		

	}
	
	if( hESDDalitzElectronAfterNFindClsTPC_mc && hESDDalitzElectronAfterNFindClsTPC_data ) {
		
	TCanvas * canvasElectronNFindClsTPC = new TCanvas("canvasElectronNFindClsTPC","",10,10,500,500);  // gives the page size
	
	canvasElectronNFindClsTPC->SetLeftMargin(0.18);
	canvasElectronNFindClsTPC->SetRightMargin(0.02);
	canvasElectronNFindClsTPC->SetTopMargin(0.03);
	canvasElectronNFindClsTPC->SetBottomMargin(0.13);
	canvasElectronNFindClsTPC->SetLogy(1);

	
	canvasElectronNFindClsTPC->cd();
	
	TH1F* hESDDalitzElectronAfterNFindClsTPC_Proj_mc   = (TH1F*)hESDDalitzElectronAfterNFindClsTPC_mc->ProjectionX("hESDDalitzElectronAfterNFindClsTPC_Proj_mc");
        TH1F* hESDDalitzElectronAfterNFindClsTPC_Proj_data = (TH1F*)hESDDalitzElectronAfterNFindClsTPC_data->ProjectionX("hESDDalitzElectronAfterNFindClsTPC_Proj_data");
	
	hESDDalitzElectronAfterNFindClsTPC_Proj_mc->SetStats(kFALSE);
	hESDDalitzElectronAfterNFindClsTPC_Proj_data->SetStats(kFALSE);
	
	
	
	hESDDalitzElectronAfterNFindClsTPC_Proj_mc->Scale(1./hESDDalitzElectronAfterNFindClsTPC_Proj_mc->GetEntries());
	hESDDalitzElectronAfterNFindClsTPC_Proj_data->Scale(1./hESDDalitzElectronAfterNFindClsTPC_Proj_data->GetEntries());
       
	
	hESDDalitzElectronAfterNFindClsTPC_Proj_mc->GetXaxis()->SetTitleOffset(1.3);
	hESDDalitzElectronAfterNFindClsTPC_Proj_data->GetXaxis()->SetTitleOffset(1.3);
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterNFindClsTPC_Proj_data,
						hESDDalitzElectronAfterNFindClsTPC_Proj_mc, 
						"", "N^{found}_{TPC clusters} / N^{Findable}_{TPC clusters}",textYAxisPtDistE,
						kTRUE,3.1,1e-5,
						kFALSE,1e-3 ,15, 
						kTRUE, 0.,1.5);	
	
	
	canvasElectronNFindClsTPC->Update();
	canvasElectronNFindClsTPC->SaveAs(Form("%s/Electron_After_MassCut_NFindClsTPC%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNFindClsTPC;		
	delete hESDDalitzElectronAfterNFindClsTPC_Proj_mc;
	delete hESDDalitzElectronAfterNFindClsTPC_Proj_data;
	
	}
	
	if( hESDDalitzElectronAfterNFindClsTPCPCut_mc && hESDDalitzElectronAfterNFindClsTPCPCut_data ) {
	
	
	TCanvas * canvasElectronNFindClsTPCPCut = new TCanvas("canvasElectronNFindClsTPCPCut","",10,10,500,500);  // gives the page size
	
	canvasElectronNFindClsTPCPCut->SetLeftMargin(0.18);
	canvasElectronNFindClsTPCPCut->SetRightMargin(0.02);
	canvasElectronNFindClsTPCPCut->SetTopMargin(0.03);
	canvasElectronNFindClsTPCPCut->SetBottomMargin(0.13);
	canvasElectronNFindClsTPCPCut->SetLogy(1);

	
	canvasElectronNFindClsTPCPCut->cd();
	
	
	hESDDalitzElectronAfterNFindClsTPCPCut_mc->SetStats(kFALSE);
	hESDDalitzElectronAfterNFindClsTPCPCut_data->SetStats(kFALSE);
	
	
		
	hESDDalitzElectronAfterNFindClsTPCPCut_mc->Scale(1./hESDDalitzElectronAfterNFindClsTPCPCut_mc->GetEntries());
	hESDDalitzElectronAfterNFindClsTPCPCut_data->Scale(1./hESDDalitzElectronAfterNFindClsTPCPCut_data->GetEntries());
	
	hESDDalitzElectronAfterNFindClsTPCPCut_mc->GetXaxis()->SetTitleOffset(1.3);
	hESDDalitzElectronAfterNFindClsTPCPCut_data->GetXaxis()->SetTitleOffset(1.3);
	
       
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterNFindClsTPCPCut_data,
						hESDDalitzElectronAfterNFindClsTPCPCut_mc, 
						"", "N^{found}_{TPC clusters} / N^{Findable}_{TPC clusters}",textYAxisPtDistE,
						kTRUE,3.1,1e-5,
						kFALSE,1e-3 ,15, 
						kTRUE, 0.,1.5);	
	
	
	canvasElectronNFindClsTPCPCut->Update();
	canvasElectronNFindClsTPCPCut->SaveAs(Form("%s/Electron_After_MassCut_NFindClsTPCPCut%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNFindClsTPCPCut;		
	delete hESDDalitzElectronAfterNFindClsTPCPCut_mc;
	delete hESDDalitzElectronAfterNFindClsTPCPCut_data;
	
	}
		
	if( hESDDalitzPositronAfterNFindClsTPC_mc && hESDDalitzPositronAfterNFindClsTPC_data ) {
	
	TCanvas * canvasPositronNFindClsTPC = new TCanvas("canvasPositronNFindClsTPC","",10,10,500,500);  // gives the page size
	
	canvasPositronNFindClsTPC->SetLeftMargin(0.18);
	canvasPositronNFindClsTPC->SetRightMargin(0.02);
	canvasPositronNFindClsTPC->SetTopMargin(0.03);
	canvasPositronNFindClsTPC->SetBottomMargin(0.13);
	canvasPositronNFindClsTPC->SetLogy(1);

	canvasPositronNFindClsTPC->cd();
	
	
	TH1F* hESDDalitzPositronAfterNFindClsTPC_Proj_mc   = (TH1F*)hESDDalitzPositronAfterNFindClsTPC_mc->ProjectionX("hESDDalitzPositronAfterNFindClsTPCPCut_Proj_mc");
        TH1F* hESDDalitzPositronAfterNFindClsTPC_Proj_data = (TH1F*)hESDDalitzPositronAfterNFindClsTPC_data->ProjectionX("hESDDalitzPositronAfterNFindClsTPCPCut_Proj_data");
	
	
	hESDDalitzPositronAfterNFindClsTPC_Proj_mc->SetStats(kFALSE);
	hESDDalitzPositronAfterNFindClsTPC_Proj_data->SetStats(kFALSE);
	
	
	
	hESDDalitzPositronAfterNFindClsTPC_Proj_mc->Scale(1./hESDDalitzPositronAfterNFindClsTPC_Proj_mc->GetEntries());
	hESDDalitzPositronAfterNFindClsTPC_Proj_data->Scale(1./hESDDalitzPositronAfterNFindClsTPC_Proj_data->GetEntries());
	
	hESDDalitzPositronAfterNFindClsTPC_Proj_mc->GetXaxis()->SetTitleOffset(1.3);
	hESDDalitzPositronAfterNFindClsTPC_Proj_data->GetXaxis()->SetTitleOffset(1.3);
	
       
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterNFindClsTPC_Proj_data,
						hESDDalitzPositronAfterNFindClsTPC_Proj_mc, 
						"", "N^{found}_{TPC clusters} / N^{Findable}_{TPC clusters}",textYAxisPtDistP,
						kTRUE, 3.1,1e-5,
						kFALSE,1e-3 ,15, 
						kTRUE, 0.,1.5);	
	
	processSystem->Draw();
	
	
	canvasPositronNFindClsTPC->Update();
	canvasPositronNFindClsTPC->SaveAs(Form("%s/Positron_After_MassCut_NFindClsTPC%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNFindClsTPC;		
	delete hESDDalitzPositronAfterNFindClsTPC_Proj_mc;
	delete hESDDalitzPositronAfterNFindClsTPC_Proj_data;
	
	}
	
	
	if( hESDDalitzPositronAfterNFindClsTPCPCut_mc && hESDDalitzPositronAfterNFindClsTPCPCut_data ) {  
	
	TCanvas * canvasPositronNFindClsTPCPCut = new TCanvas("canvasPositronNFindClsTPCPCut","",10,10,500,500);  // gives the page size
	
	canvasPositronNFindClsTPCPCut->SetLeftMargin(0.18);
	canvasPositronNFindClsTPCPCut->SetRightMargin(0.02);
	canvasPositronNFindClsTPCPCut->SetTopMargin(0.03);
	canvasPositronNFindClsTPCPCut->SetBottomMargin(0.13);
	canvasPositronNFindClsTPCPCut->SetLogy(1);

	canvasPositronNFindClsTPCPCut->cd();

	canvasPositronNFindClsTPCPCut->SetLeftMargin(0.16);
	canvasPositronNFindClsTPCPCut->cd();
	
	
	hESDDalitzPositronAfterNFindClsTPCPCut_mc->SetStats(kFALSE);
	hESDDalitzPositronAfterNFindClsTPCPCut_data->SetStats(kFALSE);
	
	
	
	hESDDalitzPositronAfterNFindClsTPCPCut_mc->Scale(1./hESDDalitzPositronAfterNFindClsTPCPCut_mc->GetEntries());
	hESDDalitzPositronAfterNFindClsTPCPCut_data->Scale(1./hESDDalitzPositronAfterNFindClsTPCPCut_data->GetEntries());
	
	hESDDalitzPositronAfterNFindClsTPCPCut_mc->GetXaxis()->SetTitleOffset(1.3);
	hESDDalitzPositronAfterNFindClsTPCPCut_data->GetXaxis()->SetTitleOffset(1.3);
	
       
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterNFindClsTPCPCut_data,
						hESDDalitzPositronAfterNFindClsTPCPCut_mc, 
						"", "N^{found}_{TPC clusters} / N^{Findable}_{TPC clusters}",textYAxisPtDistP,
						kTRUE, 3.1,1e-5,
						kFALSE,1e-3 ,15, 
						kTRUE, 0.,1.5);	
	
	
	canvasPositronNFindClsTPCPCut->Update();
	canvasPositronNFindClsTPCPCut->SaveAs(Form("%s/Positron_After_MassCut_NFindClsTPCPCut%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNFindClsTPCPCut;		
	delete hESDDalitzPositronAfterNFindClsTPCPCut_mc;
	delete hESDDalitzPositronAfterNFindClsTPCPCut_data;
	
	
	}
	
	if( hESDDalitzElectronAfterNClsTPC_data ) {
	
	TCanvas * canvasElectronNClsTPC_data = new TCanvas("canvasElectronNClsTPC_data","",10,10,500,500);  // gives the page size
	canvasElectronNClsTPC_data->SetLeftMargin(0.16);
	canvasElectronNClsTPC_data->cd();
	canvasElectronNClsTPC_data->SetLogy(1);
	
	hESDDalitzElectronAfterNClsTPC_data->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterNClsTPC_data,
						 "", "p (GeV/c)", "dE/dx","Data",
						 kTRUE,0., 15.,
						 kTRUE, 0.,200);
	
	
	canvasElectronNClsTPC_data->Update();
	canvasElectronNClsTPC_data->SaveAs(Form("%s/Electron_After_MassCut_NClsTPC_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNClsTPC_data;
	
	}
	
	if( hESDDalitzElectronAfterNClsTPC_mc ) {
	
	TCanvas * canvasElectronNClsTPC_mc = new TCanvas("canvasElectronNClsTPC_mc","",10,10,500,500);  // gives the page size
	canvasElectronNClsTPC_mc->SetLeftMargin(0.16);
	canvasElectronNClsTPC_mc->cd();
	canvasElectronNClsTPC_mc->SetLogy(1);
	hESDDalitzElectronAfterNClsTPC_mc->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterNClsTPC_mc,
						 "", "p (GeV/c)", "dE/dx","Data",
						 kTRUE,0., 15.,
						 kTRUE, 0.,200);
	
	
	canvasElectronNClsTPC_mc->Update();
	canvasElectronNClsTPC_mc->SaveAs(Form("%s/Electron_After_MassCut_NClsTPC_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNClsTPC_mc;
	
	}
	
	if( hESDDalitzPositronAfterNClsTPC_data ) {
			
	TCanvas * canvasPositronNClsTPC_data = new TCanvas("canvasPositronNClsTPC_data","",10,10,500,500);  // gives the page size
	canvasPositronNClsTPC_data->SetLeftMargin(0.16);
	canvasPositronNClsTPC_data->cd();
	canvasPositronNClsTPC_data->SetLogy(1);
	
	hESDDalitzPositronAfterNClsTPC_data->SetStats(kFALSE);
       
	
	DrawAutoGammaHisto2D(	hESDDalitzPositronAfterNClsTPC_data,
						 "", "p (GeV/c)", "dE/dx","Data",
						 kTRUE,0., 15.,
						 kTRUE, 0.,200);
	
	
	canvasPositronNClsTPC_data->Update();
	canvasPositronNClsTPC_data->SaveAs(Form("%s/Positron_After_MassCut_NClsTPC_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNClsTPC_data;		
	
	}
	
	if( hESDDalitzPositronAfterNClsTPC_mc ) {
	
	TCanvas * canvasPositronNClsTPC_mc = new TCanvas("canvasPositronNClsTPC_mc","",10,10,500,500);  // gives the page size
	canvasPositronNClsTPC_mc->SetLeftMargin(0.16);
	canvasPositronNClsTPC_mc->cd();
	canvasPositronNClsTPC_mc->SetLogy(1);
	hESDDalitzPositronAfterNClsTPC_mc->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterNClsTPC_mc,
						 "", "p (GeV/c)", "#mathrm{d}E/#mathrm{d}x","Data",
						 kTRUE,0., 15.,
						 kTRUE, 0.,200);
	
	
	canvasPositronNClsTPC_mc->Update();
	canvasPositronNClsTPC_mc->SaveAs(Form("%s/Positron_After_MassCut_NClsTPC_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNClsTPC_mc;
	
	}
		
	if( hESDDalitzElectronAfterNClsTPC_mc && hESDDalitzElectronAfterNClsTPC_data ){
	
	TCanvas * canvasElectronNClsTPC = new TCanvas("canvasElectronNClsTPC","",10,10,500,500);  // gives the page size
	canvasElectronNClsTPC->SetLogy(1);

	canvasElectronNClsTPC->SetLeftMargin(0.16);
	canvasElectronNClsTPC->cd();
	
	TH1F* hESDDalitzElectronAfterNClsTPC_Proj_mc   = (TH1F*)hESDDalitzElectronAfterNClsTPC_mc->ProjectionX("hESDDalitzElectronAfterNFindClsTPC_Proj_mc");
        TH1F* hESDDalitzElectronAfterNClsTPC_Proj_data = (TH1F*)hESDDalitzElectronAfterNClsTPC_data->ProjectionX("hESDDalitzElectronAfterNFindClsTPC_Proj_data");
	
	hESDDalitzElectronAfterNClsTPC_Proj_mc->SetStats(kFALSE);
	hESDDalitzElectronAfterNClsTPC_Proj_data->SetStats(kFALSE);
	
	
	  
	hESDDalitzElectronAfterNClsTPC_Proj_mc->Scale(1./hESDDalitzElectronAfterNClsTPC_Proj_mc->Integral());
	hESDDalitzElectronAfterNClsTPC_Proj_data->Scale(1./hESDDalitzElectronAfterNClsTPC_Proj_data->Integral());
	
	hESDDalitzElectronAfterNClsTPC_Proj_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzElectronAfterNClsTPC_Proj_data->GetXaxis()->SetTitleOffset(1.0);
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterNClsTPC_Proj_data,
						hESDDalitzElectronAfterNClsTPC_Proj_mc, 
						"",textXAxisTPCclsE.Data(),textYAxisTPCclsE,
						kTRUE, 2.5,1e-4,
						kFALSE,1e-3 ,15, 
						kTRUE, 0.,200.);	
	
	
	canvasElectronNClsTPC->Update();
	canvasElectronNClsTPC->SaveAs(Form("%s/Electron_After_MassCut_NClsTPC%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNClsTPC;		
	delete hESDDalitzElectronAfterNClsTPC_Proj_mc;
	delete hESDDalitzElectronAfterNClsTPC_Proj_data;
	
	}
	
	if( hESDDalitzElectronAfterNClsTPCPCut_mc && hESDDalitzElectronAfterNClsTPCPCut_data ) {
	
	
	TCanvas * canvasElectronNClsTPCPCut = new TCanvas("canvasElectronNClsTPCPCut","",10,10,500,500);  // gives the page size
	canvasElectronNClsTPCPCut->SetLogy(1);

	canvasElectronNClsTPCPCut->SetLeftMargin(0.16);
	canvasElectronNClsTPCPCut->cd();
	
	
	hESDDalitzElectronAfterNClsTPCPCut_mc->SetStats(kFALSE);
	hESDDalitzElectronAfterNClsTPCPCut_data->SetStats(kFALSE);
	
	
	//hESDDalitzElectronAfterNClsTPCPCut_mc->Sumw2();
	//hESDDalitzElectronAfterNClsTPCPCut_data->Sumw2();
       
	hESDDalitzElectronAfterNClsTPCPCut_mc->Scale(1./hESDDalitzElectronAfterNClsTPCPCut_mc->Integral());
	hESDDalitzElectronAfterNClsTPCPCut_data->Scale(1./hESDDalitzElectronAfterNClsTPCPCut_data->Integral());
	
	hESDDalitzElectronAfterNClsTPCPCut_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzElectronAfterNClsTPCPCut_data->GetXaxis()->SetTitleOffset(1.0);
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterNClsTPCPCut_data,
						hESDDalitzElectronAfterNClsTPCPCut_mc, 
						"", "N_{TPC_{cls}} e^{-}",textYAxisTPCclsE,
						kTRUE,2.5,1e-4,
						kFALSE,1e-3 ,15, 
						kTRUE, 0.,200.);	
	
	
	canvasElectronNClsTPCPCut->Update();
	canvasElectronNClsTPCPCut->SaveAs(Form("%s/Electron_After_MassCut_NClsTPCPCut%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNClsTPCPCut;		
	delete hESDDalitzElectronAfterNClsTPCPCut_mc;
	delete hESDDalitzElectronAfterNClsTPCPCut_data;
	
	}
	
	if( hESDDalitzPositronAfterNClsTPC_mc && hESDDalitzPositronAfterNClsTPC_data ) {
	
	TCanvas * canvasPositronNClsTPC = new TCanvas("canvasPositronNClsTPC","",10,10,500,500);  // gives the page size
	canvasPositronNClsTPC->SetLogy(1);

	canvasPositronNClsTPC->SetLeftMargin(0.16);
	canvasPositronNClsTPC->cd();
	
	TH1F* hESDDalitzPositronAfterNClsTPC_Proj_mc   = (TH1F*)hESDDalitzPositronAfterNClsTPC_mc->ProjectionX("hESDDalitzPositronAfterNClsTPC_Proj_mc");
        TH1F* hESDDalitzPositronAfterNClsTPC_Proj_data = (TH1F*)hESDDalitzPositronAfterNClsTPC_data->ProjectionX("hESDDalitzPositronAfterNClsTPC_Proj_data");
	
	hESDDalitzPositronAfterNClsTPC_Proj_mc->SetStats(kFALSE);
	hESDDalitzPositronAfterNClsTPC_Proj_data->SetStats(kFALSE);
	
	
	//hESDDalitzPositronAfterNClsTPC_Proj_mc->Sumw2();
	//hESDDalitzPositronAfterNClsTPC_Proj_data->Sumw2();
       
	hESDDalitzPositronAfterNClsTPC_Proj_mc->Scale(1./hESDDalitzPositronAfterNClsTPC_Proj_mc->Integral());
	hESDDalitzPositronAfterNClsTPC_Proj_data->Scale(1./hESDDalitzPositronAfterNClsTPC_Proj_data->Integral());
	
	hESDDalitzPositronAfterNClsTPC_Proj_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzPositronAfterNClsTPC_Proj_data->GetXaxis()->SetTitleOffset(1.0);
	
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterNClsTPC_Proj_data,
						hESDDalitzPositronAfterNClsTPC_Proj_mc, 
						"",textXAxisTPCclsP,textYAxisTPCclsP,
						kTRUE,2.5,1e-4,
						kFALSE,1e-3 ,15, 
						kFALSE, 0.,200.);	
	
	
	canvasPositronNClsTPC->Update();
	canvasPositronNClsTPC->SaveAs(Form("%s/Positron_After_MassCut_NClsTPC%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNClsTPC;		
	delete hESDDalitzPositronAfterNClsTPC_Proj_mc;
	delete hESDDalitzPositronAfterNClsTPC_Proj_data;
	
	}
	
	if( hESDDalitzPositronAfterNClsTPCPCut_mc && hESDDalitzPositronAfterNClsTPCPCut_data ){
	
	TCanvas * canvasPositronNClsTPCPCut = new TCanvas("canvasPositronNClsTPCPCut","",10,10,500,500);  // gives the page size
	canvasPositronNClsTPCPCut->SetLogy(1);

	canvasPositronNClsTPCPCut->SetLeftMargin(0.16);
	canvasPositronNClsTPCPCut->SetTopMargin(0.05);
	canvasPositronNClsTPCPCut->cd();
	
	
	hESDDalitzPositronAfterNClsTPCPCut_mc->SetStats(kFALSE);
	hESDDalitzPositronAfterNClsTPCPCut_data->SetStats(kFALSE);
	
	
	//hESDDalitzPositronAfterNClsTPCPCut_mc->Sumw2();
	//hESDDalitzPositronAfterNClsTPCPCut_data->Sumw2();
       
	hESDDalitzPositronAfterNClsTPCPCut_mc->Scale(1./hESDDalitzPositronAfterNClsTPCPCut_mc->GetEntries());
	hESDDalitzPositronAfterNClsTPCPCut_data->Scale(1./hESDDalitzPositronAfterNClsTPCPCut_data->GetEntries());
	
	hESDDalitzPositronAfterNClsTPCPCut_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzPositronAfterNClsTPCPCut_mc->GetXaxis()->SetTitleOffset(1.0);
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterNClsTPCPCut_data,
						hESDDalitzPositronAfterNClsTPCPCut_mc, 
						"", "Number of TPC clusters e^{+}",textYAxisPtDistP,
						kTRUE,2.5,1e-4,
						kFALSE,1e-3 ,15, 
						kTRUE, 0.,200.);	
	
	
	canvasPositronNClsTPCPCut->Update();
	canvasPositronNClsTPCPCut->SaveAs(Form("%s/Positron_After_MassCut_NClsTPCPCut%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNClsTPCPCut;		
	delete hESDDalitzPositronAfterNClsTPCPCut_mc;
	delete hESDDalitzPositronAfterNClsTPCPCut_data;
	
	}
	
	
	if( hESDDalitzPositronAfterNClsTPC_mc && hESDDalitzPositronAfterNClsTPC_data && hESDDalitzElectronAfterNClsTPC_mc && hESDDalitzElectronAfterNClsTPC_data) {
	
	TCanvas * canvasElectronPositronNClsTPC = new TCanvas("canvasElectronPositronNClsTPC","",10,10,500,500);  // gives the page size
	canvasElectronPositronNClsTPC->SetLogy(1);
	canvasElectronPositronNClsTPC->cd(1);
	
	canvasElectronPositronNClsTPC->SetLeftMargin(0.17);
	canvasElectronPositronNClsTPC->SetRightMargin(0.02);
	
	TH1F* hESDDalitzElectronPositronAfterNClsTPC_Proj_mc   = (TH1F*)hESDDalitzPositronAfterNClsTPC_mc->ProjectionX("hESDDalitzElectronPositronAfterNClsTPC_Proj_mc");
	TH1F* hESDDalitzElectronAfterNClsTPC_Proj_mc           = (TH1F*)hESDDalitzElectronAfterNClsTPC_mc->ProjectionX("hESDDalitzElectronAfterNClsTPC_Proj_mc");
	
	hESDDalitzElectronPositronAfterNClsTPC_Proj_mc->Add(hESDDalitzElectronAfterNClsTPC_Proj_mc,1.0);
	
       
        TH1F* hESDDalitzElectronPositronAfterNClsTPC_Proj_data   = (TH1F*)hESDDalitzPositronAfterNClsTPC_data->ProjectionX("hESDDalitzElectronPositronAfterNClsTPC_Proj_data");
	TH1F* hESDDalitzElectronAfterNClsTPC_Proj_data           = (TH1F*)hESDDalitzElectronAfterNClsTPC_data->ProjectionX("hESDDalitzElectronAfterNClsTPC_Proj_data");
	
	hESDDalitzElectronPositronAfterNClsTPC_Proj_data->Add(hESDDalitzElectronAfterNClsTPC_Proj_data,1.0);
	
	
	hESDDalitzElectronPositronAfterNClsTPC_Proj_mc->SetStats(kFALSE);
	hESDDalitzElectronPositronAfterNClsTPC_Proj_data->SetStats(kFALSE);
	
	   
	hESDDalitzElectronPositronAfterNClsTPC_Proj_mc->Scale(1./hESDDalitzElectronPositronAfterNClsTPC_Proj_mc->Integral());
	hESDDalitzElectronPositronAfterNClsTPC_Proj_data->Scale(1./hESDDalitzElectronPositronAfterNClsTPC_Proj_data->Integral());
	
	hESDDalitzElectronPositronAfterNClsTPC_Proj_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzElectronPositronAfterNClsTPC_Proj_data->GetXaxis()->SetTitleOffset(1.0);
	
	
	
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronPositronAfterNClsTPC_Proj_data,
						hESDDalitzElectronPositronAfterNClsTPC_Proj_mc, 
						"", "N_{TPCcls} e^{#pm}",textYAxisTPCclsEP,
						kTRUE,1.5,1e-4,
						kFALSE,1e-3 ,12, 
						kFALSE, 0.,200.);	
	
	
	canvasElectronPositronNClsTPC->Update();
	canvasElectronPositronNClsTPC->SaveAs(Form("%s/ElectronPositron_After_MassCut_NClsTPC%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronNClsTPC;		
	delete hESDDalitzElectronPositronAfterNClsTPC_Proj_mc;
	delete hESDDalitzElectronPositronAfterNClsTPC_Proj_data;
	
	}
	
	
		
	
	
	
	if( hESDDalitzPositronAfterNCrossedRowsTPC_mc ) {
		
	TCanvas * canvasPositronNCrossedRowsTPC_mc = new TCanvas("canvasPositronNCrossedRowsTPC_mc","",10,10,500,500);  // gives the page size
	canvasPositronNCrossedRowsTPC_mc->SetLeftMargin(0.12);
	canvasPositronNCrossedRowsTPC_mc->SetRightMargin(0.13);
	canvasPositronNCrossedRowsTPC_mc->SetTopMargin(0.05);
	canvasPositronNCrossedRowsTPC_mc->SetBottomMargin(0.12);
	
	
	hESDDalitzPositronAfterNCrossedRowsTPC_mc->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterNCrossedRowsTPC_mc,
						 "", "Number of crossed TPC rows e^{+}", "p_{T} (GeV/c)","MC",
						 kTRUE,0.,5.,
						 kTRUE, 0.,200,1.2,1.0);
	
	
	canvasPositronNCrossedRowsTPC_mc->Update();
	canvasPositronNCrossedRowsTPC_mc->SaveAs(Form("%s/Positron_After_MassCut_NCrossedRowsTPC_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNCrossedRowsTPC_mc;
	}
	
	
	if( hESDDalitzPositronAfterNCrossedRowsTPC_data ) {
	
	TCanvas * canvasPositronNCrossedRowsTPC_data = new TCanvas("canvasPositronNCrossedRowsTPC_data","",10,10,500,500);  // gives the page size
	canvasPositronNCrossedRowsTPC_data->SetLeftMargin(0.12);
	canvasPositronNCrossedRowsTPC_data->SetRightMargin(0.13);
	canvasPositronNCrossedRowsTPC_data->SetTopMargin(0.05);
	canvasPositronNCrossedRowsTPC_data->SetBottomMargin(0.12);
	
	
	hESDDalitzPositronAfterNCrossedRowsTPC_data->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterNCrossedRowsTPC_data,
						 "", "Number of crossed TPC rows e^{+}", "p_{T} (GeV/c)","Data",
						 kTRUE,0.,5.,
						 kTRUE, 0.,200,1.2,1.0);
	
	
	canvasPositronNCrossedRowsTPC_data->Update();
	canvasPositronNCrossedRowsTPC_data->SaveAs(Form("%s/Positron_After_MassCut_NCrossedRowsTPC_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNCrossedRowsTPC_data;
	
	}

	if( hESDDalitzElectronAfterNCrossedRowsTPC_mc ) {
	
	TCanvas * canvasElectronNCrossedRowsTPC_mc = new TCanvas("canvasElectronNCrossedRowsTPC_mc","",10,10,500,500);  // gives the page size
	canvasElectronNCrossedRowsTPC_mc->SetLeftMargin(0.12);
	canvasElectronNCrossedRowsTPC_mc->SetRightMargin(0.13);
	canvasElectronNCrossedRowsTPC_mc->SetTopMargin(0.05);
	canvasElectronNCrossedRowsTPC_mc->SetBottomMargin(0.12);
	
	
	hESDDalitzElectronAfterNCrossedRowsTPC_mc->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterNCrossedRowsTPC_mc,
						 "", "Number of crossed TPC rows e^{-}", "p_{T} (GeV/c)","MC",
						 kTRUE,0.,5.,
						 kTRUE, 0.,200,1.2,1.0);
	
	
	canvasElectronNCrossedRowsTPC_mc->Update();
	canvasElectronNCrossedRowsTPC_mc->SaveAs(Form("%s/Electron_After_MassCut_NCrossedRowsTPC_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNCrossedRowsTPC_mc;
	
	}
	
	
	if( hESDDalitzElectronAfterNCrossedRowsTPC_data ) {
	
	TCanvas * canvasElectronNCrossedRowsTPC_data = new TCanvas("canvasElectronNCrossedRowsTPC_data","",10,10,500,500);  // gives the page size
	canvasElectronNCrossedRowsTPC_data->SetLeftMargin(0.12);
	canvasElectronNCrossedRowsTPC_data->SetRightMargin(0.13);
	canvasElectronNCrossedRowsTPC_data->SetTopMargin(0.05);
	canvasElectronNCrossedRowsTPC_data->SetBottomMargin(0.12);
	
	
	hESDDalitzElectronAfterNCrossedRowsTPC_data->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterNCrossedRowsTPC_data,
						 "", "N_{TPC_{CrossedRows}}", "p_{T} (GeV/c)","Data",
						 kTRUE,0.,5.,
						 kTRUE, 0.,200,1.2,1.0);
	
	
	canvasElectronNCrossedRowsTPC_data->Update();
	canvasElectronNCrossedRowsTPC_data->SaveAs(Form("%s/Electron_After_MassCut_NCrossedRowsTPC_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNCrossedRowsTPC_data;
	
	}
	
	
	if( hESDDalitzPositronAfterNCrossedRowsTPC_mc && hESDDalitzPositronAfterNCrossedRowsTPC_data ){
	
	
	TCanvas * canvasPositronNCrossedRowsTPC = new TCanvas("canvasPositronNCrossedRowsTPC","",10,10,500,500);  // gives the page size
	canvasPositronNCrossedRowsTPC->SetLogy(1);

	canvasPositronNCrossedRowsTPC->SetLeftMargin(0.16);
	//canvasPositronNCrossedRowsTPC->SetTopMargin(0.05);
	canvasPositronNCrossedRowsTPC->cd();
	
	TH1F* hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc   = (TH1F*)hESDDalitzPositronAfterNCrossedRowsTPC_mc->ProjectionX("hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc");
        TH1F* hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data = (TH1F*)hESDDalitzPositronAfterNCrossedRowsTPC_data->ProjectionX("hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data");
	
	hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc->SetStats(kFALSE);
	hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data->SetStats(kFALSE);
	
	
	//hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc->Sumw2();
	//hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data->Sumw2();
       
	hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc->Scale(1./hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc->GetEntries());
	hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data->Scale(1./hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data->GetEntries());
	
	hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data->GetXaxis()->SetTitleOffset(1.0);
	
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data,
						hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc, 
						"", "N_{TPC_{CrossedRows}} e^{+}",textYAxisTPCcrossedRowsP,
						kFALSE,2.5,1e-4,
						kTRUE,1e-4 ,9e-1, 
						kTRUE, 0.,200.);	
	
	
	canvasPositronNCrossedRowsTPC->Update();
	canvasPositronNCrossedRowsTPC->SaveAs(Form("%s/Positron_After_MassCut_NCrossedRowsTPC%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNCrossedRowsTPC;		
	delete hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc;
	delete hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data;
	
	}
	
	if( hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc && hESDDalitzPositronAfterNCrossedRowsTPCPCut_data ) {
	
	TCanvas * canvasPositronNCrossedRowsTPCPCut = new TCanvas("canvasPositronNCrossedRowsTPCPCut","",10,10,500,500);  // gives the page size
	canvasPositronNCrossedRowsTPCPCut->SetLogy(1);

	canvasPositronNCrossedRowsTPCPCut->SetLeftMargin(0.16);
	//canvasPositronNCrossedRowsTPCPCut->SetTopMargin(0.05);
	canvasPositronNCrossedRowsTPCPCut->cd();
	
		
	hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc->SetStats(kFALSE);
	hESDDalitzPositronAfterNCrossedRowsTPCPCut_data->SetStats(kFALSE);
	
	
	//hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc->Sumw2();
	//hESDDalitzPositronAfterNCrossedRowsTPCPCut_data->Sumw2();
       
	hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc->Scale(1./hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc->GetEntries());
	hESDDalitzPositronAfterNCrossedRowsTPCPCut_data->Scale(1./hESDDalitzPositronAfterNCrossedRowsTPCPCut_data->GetEntries());
	
	hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzPositronAfterNCrossedRowsTPCPCut_data->GetXaxis()->SetTitleOffset(1.0);
	
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterNCrossedRowsTPCPCut_data,
						hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc, 
						"", "N_{TPC_{CrossedRows}} e^{+}",textYAxisTPCcrossedRowsP,
						kFALSE,2.5,1e-4,
						kTRUE,1e-4,9e-1, 
						kTRUE, 0.,200.);	
	
	
	canvasPositronNCrossedRowsTPCPCut->Update();
	canvasPositronNCrossedRowsTPCPCut->SaveAs(Form("%s/Positron_After_MassCut_NCrossedRowsTPCPCut%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNCrossedRowsTPCPCut;		
	delete hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc;
	delete hESDDalitzPositronAfterNCrossedRowsTPCPCut_data;
	
	}
	
	
	
	if( hESDDalitzElectronAfterNCrossedRowsTPC_mc && hESDDalitzElectronAfterNCrossedRowsTPC_data ) {
	
	TCanvas * canvasElectronNCrossedRowsTPC = new TCanvas("canvasElectronNCrossedRowsTPC","",10,10,500,500);  // gives the page size
	canvasElectronNCrossedRowsTPC->SetLogy(1);

	canvasElectronNCrossedRowsTPC->SetLeftMargin(0.16);
	canvasElectronNCrossedRowsTPC->cd();
	
	TH1F* hESDDalitzElectronAfterNCrossedRowsTPC_Proj_mc   = (TH1F*)hESDDalitzElectronAfterNCrossedRowsTPC_mc->ProjectionX("hESDDalitzElectronAfterNCrossedRowsTPC_Proj_mc");
        TH1F* hESDDalitzElectronAfterNCrossedRowsTPC_Proj_data = (TH1F*)hESDDalitzElectronAfterNCrossedRowsTPC_data->ProjectionX("hESDDalitzElectronAfterNCrossedRowsTPC_Proj_data");
	
	hESDDalitzElectronAfterNCrossedRowsTPC_Proj_mc->SetStats(kFALSE);
	hESDDalitzElectronAfterNCrossedRowsTPC_Proj_data->SetStats(kFALSE);
	
	
	//hESDDalitzElectronAfterNCrossedRowsTPC_Proj_mc->Sumw2();
	//hESDDalitzElectronAfterNCrossedRowsTPC_Proj_data->Sumw2();
       
	hESDDalitzElectronAfterNCrossedRowsTPC_Proj_mc->Scale(1./hESDDalitzElectronAfterNCrossedRowsTPC_Proj_mc->Integral());
	hESDDalitzElectronAfterNCrossedRowsTPC_Proj_data->Scale(1./hESDDalitzElectronAfterNCrossedRowsTPC_Proj_data->Integral());
	
	hESDDalitzElectronAfterNCrossedRowsTPC_Proj_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzElectronAfterNCrossedRowsTPC_Proj_data->GetXaxis()->SetTitleOffset(1.0);
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterNCrossedRowsTPC_Proj_data,
						hESDDalitzElectronAfterNCrossedRowsTPC_Proj_mc, 
						"", "N_{TPC_{CrossedRows}} e^{-}",textYAxisTPCcrossedRowsE,
						kFALSE,2.5,1e-4,
						kTRUE,1e-4 ,9e-1, 
						kTRUE, 0.,200.);	
	
	
	canvasElectronNCrossedRowsTPC->Update();
	canvasElectronNCrossedRowsTPC->SaveAs(Form("%s/Electron_After_MassCut_NCrossedRowsTPC%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNCrossedRowsTPC;		
	delete hESDDalitzElectronAfterNCrossedRowsTPC_Proj_mc;
	delete hESDDalitzElectronAfterNCrossedRowsTPC_Proj_data;
	
	}
	
	if( hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc && hESDDalitzElectronAfterNCrossedRowsTPCPCut_data ){
	
	
	TCanvas * canvasElectronNCrossedRowsTPCPCut = new TCanvas("canvasElectronNCrossedRowsTPCPCut","",10,10,500,500);  // gives the page size
	canvasElectronNCrossedRowsTPCPCut->SetLogy(1);

	canvasElectronNCrossedRowsTPCPCut->SetLeftMargin(0.16);
	canvasElectronNCrossedRowsTPCPCut->SetTopMargin(0.05);
	canvasElectronNCrossedRowsTPCPCut->cd();
	
		
	hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc->SetStats(kFALSE);
	hESDDalitzElectronAfterNCrossedRowsTPCPCut_data->SetStats(kFALSE);
	
	
	//hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc->Sumw2();
	//hESDDalitzElectronAfterNCrossedRowsTPCPCut_data->Sumw2();
       
	hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc->Scale(1./hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc->GetEntries());
	hESDDalitzElectronAfterNCrossedRowsTPCPCut_data->Scale(1./hESDDalitzElectronAfterNCrossedRowsTPCPCut_data->GetEntries());
	
	hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzElectronAfterNCrossedRowsTPCPCut_data->GetXaxis()->SetTitleOffset(1.0);
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterNCrossedRowsTPCPCut_data,
						hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc, 
						"", "Number of TPC clusters e^{-} ",textYAxisPtDistE,
						kTRUE,2.5,1e-4,
						kFALSE,1e-3 ,15, 
						kTRUE, 0.,200.);	
	
	
	canvasElectronNCrossedRowsTPCPCut->Update();
	canvasElectronNCrossedRowsTPCPCut->SaveAs(Form("%s/Electron_After_MassCut_NCrossedRowsTPCPCut%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNCrossedRowsTPCPCut;		
	delete hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc;
	delete hESDDalitzElectronAfterNCrossedRowsTPCPCut_data;
	
	}
	
	
	
	if( hESDDalitzElectronAfterNCrossedRowsTPC_mc && hESDDalitzElectronAfterNCrossedRowsTPC_data && hESDDalitzPositronAfterNCrossedRowsTPC_mc && hESDDalitzPositronAfterNCrossedRowsTPC_data ) {
	
	TCanvas * canvasElectronNCrossedRowsTPC = new TCanvas("canvasElectronNCrossedRowsTPC","",10,10,500,500);  // gives the page size
	canvasElectronNCrossedRowsTPC->SetLogy(1);

	canvasElectronNCrossedRowsTPC->SetLeftMargin(0.16);
	canvasElectronNCrossedRowsTPC->cd();
	
	TH1F* hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc   = (TH1F*)hESDDalitzElectronAfterNCrossedRowsTPC_mc->ProjectionX("hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc");
	TH1F* hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc           = (TH1F*)hESDDalitzPositronAfterNCrossedRowsTPC_mc->ProjectionX("hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc");
       
	hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc->Add(hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc,1.0);
	
	
        TH1F* hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_data   = (TH1F*)hESDDalitzElectronAfterNCrossedRowsTPC_data->ProjectionX("hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc");
	TH1F* hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data           = (TH1F*)hESDDalitzPositronAfterNCrossedRowsTPC_data->ProjectionX("hESDDalitzPositronAfterNCrossedRowsTPC_Proj_mc");
       
	hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_data->Add(hESDDalitzPositronAfterNCrossedRowsTPC_Proj_data,1.0);
	
	
	
	hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc->SetStats(kFALSE);
	hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_data->SetStats(kFALSE);
	
	
	//hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc->Sumw2();
	//hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_data->Sumw2();
       
	hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc->Scale(1./hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc->Integral());
	hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_data->Scale(1./hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_data->Integral());
	
	hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc->GetXaxis()->SetTitleOffset(1.0);
	hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_data->GetXaxis()->SetTitleOffset(1.0);
	
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_data,
						hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc, 
						"", "N_{TPCcrossedRows} e^{#pm}",textYAxisTPCcrossedRowsEP,
						kTRUE,1.5,1e-4,
						kFALSE,1e-4 ,9e-1, 
						kTRUE, 0.,200.);	
	
	
	canvasElectronNCrossedRowsTPC->Update();
	canvasElectronNCrossedRowsTPC->SaveAs(Form("%s/ElectronPositron_After_MassCut_NCrossedRowsTPC%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNCrossedRowsTPC;		
	delete hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_mc;
	delete hESDDalitzElectronPositronAfterNCrossedRowsTPC_Proj_data;
	
	}
	
	if( hESDDalitzElectronAfterNClsITS_mc && hESDDalitzElectronAfterNClsITS_data ) {
	  
	    Double_t YMax = 0;
	    Double_t YMin = 0;
	    Double_t XMin = 0;
	    Double_t XMax = 0;
	  
	 if( energy.CompareTo("7TeV") == 0 ){
	      
	    YMin = 1e-4;   YMax = 2;
	    XMin = -0.5;   XMax = 7.0; 
	
	  } else if( energy.CompareTo("2.76TeV") == 0 ){
	    
	    YMin = 1e-4;   YMax = 2;
	    XMin = -0.5;   XMax = 7.0; 
	    
	  } else if( energy.CompareTo("pPb_5.023TeV") == 0 ){
	    
	    YMin = 1e-4;   YMax = 2;
	    XMin = -0.5;   XMax = 7.0; 
	    
	  }
   
   
	TCanvas * canvasElectronNClsITS = new TCanvas("canvasElectronNClsITS","",10,10,500,500);  // gives the page size
	canvasElectronNClsITS->SetLogy(1);

	canvasElectronNClsITS->SetLeftMargin(0.16);
	canvasElectronNClsITS->cd();
	
	
	
	TH1F*  hESDDalitzElectronAfterNClsITSRebin_mc   = (TH1F*) hESDDalitzElectronAfterNClsITS_mc->Clone("hESDDalitzElectronAfterNClsITSRebin_mc");
	TH1F*  hESDDalitzElectronAfterNClsITSRebin_data = (TH1F*) hESDDalitzElectronAfterNClsITS_data->Clone("hESDDalitzElectronAfterNClsITSRebin_data");
	
	hESDDalitzElectronAfterNClsITSRebin_mc->SetStats(kFALSE);
	hESDDalitzElectronAfterNClsITSRebin_data->SetStats(kFALSE);
	
	
	//hESDDalitzElectronAfterNClsITSRebin_mc->Sumw2();
	//hESDDalitzElectronAfterNClsITSRebin_data->Sumw2();
       
	hESDDalitzElectronAfterNClsITSRebin_mc->Scale(1./hESDDalitzElectronAfterNClsITSRebin_mc->Integral());
	hESDDalitzElectronAfterNClsITSRebin_data->Scale(1./hESDDalitzElectronAfterNClsITSRebin_data->Integral());
	
	DrawAutoGammaHistosTemp( hESDDalitzElectronAfterNClsITSRebin_data,
						hESDDalitzElectronAfterNClsITSRebin_mc, 
						"", "N_{ITS_{cls}} e^{-}",textYAxisITSclsE,
						kFALSE, 15,1e-2,
						kTRUE,  YMin,YMax, 
						kTRUE,  XMin,XMax,1.,1.6);	
	
	
	canvasElectronNClsITS->Update();
	canvasElectronNClsITS->SaveAs(Form("%s/Electron_After_MassCut_NClsITS%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronNClsITS;		
	
	
	}
	
	
	if( hESDDalitzPositronAfterNClsITS_mc && hESDDalitzPositronAfterNClsITS_data ) {
	  
	    Double_t YMax = 0;
	    Double_t YMin = 0;
	    Double_t XMin = 0;
	    Double_t XMax = 0;
	  
	   if( energy.CompareTo("7TeV") == 0 ){
	      
	    YMin = 1e-4;   YMax = 2;
	    XMin = -0.5;   XMax = 7.0; 
	
	  } else if( energy.CompareTo("2.76TeV") == 0 ){
	    
	    YMin = 1e-4;   YMax = 2;
	    XMin = -0.5;   XMax = 7.0; 
	    
	  } else if( energy.CompareTo("pPb_5.023TeV") == 0 ){
	    
	    YMin = 1e-4;   YMax = 2;
	    XMin = -0.5;   XMax = 7.0; 
	    
	  }
	
	TCanvas * canvasPositronNClsITS = new TCanvas("canvasPositronNClsITS","",10,10,500,500);  // gives the page size
	canvasPositronNClsITS->SetLogy(1);

	canvasPositronNClsITS->SetLeftMargin(0.16);
	canvasPositronNClsITS->cd();
	
	
	
	TH1F*  hESDDalitzPositronAfterNClsITSRebin_mc   = (TH1F*) hESDDalitzPositronAfterNClsITS_mc->Clone("hESDDalitzPositronAfterNClsITSRebin_mc");
	TH1F*  hESDDalitzPositronAfterNClsITSRebin_data = (TH1F*) hESDDalitzPositronAfterNClsITS_data->Clone("hESDDalitzPositronAfterNClsITSRebin_data");
	
	hESDDalitzPositronAfterNClsITSRebin_mc->SetStats(kFALSE);
	hESDDalitzPositronAfterNClsITSRebin_data->SetStats(kFALSE);
	
	
	//hESDDalitzPositronAfterNClsITSRebin_mc->Sumw2();
	//hESDDalitzPositronAfterNClsITSRebin_data->Sumw2();
       
	hESDDalitzPositronAfterNClsITSRebin_mc->Scale(1./hESDDalitzPositronAfterNClsITSRebin_mc->Integral());
	hESDDalitzPositronAfterNClsITSRebin_data->Scale(1./hESDDalitzPositronAfterNClsITSRebin_data->Integral());
	
	DrawAutoGammaHistosTemp( hESDDalitzPositronAfterNClsITSRebin_data,
						hESDDalitzPositronAfterNClsITSRebin_mc, 
						"", "N_{ITS_{cls}} e^{+}",textYAxisITSclsP,
						kFALSE, 15,1e-2,
						kTRUE,  YMin,YMax, 
						kTRUE,  XMin,XMax,1.,1.6);	
	
	
	canvasPositronNClsITS->Update();
	canvasPositronNClsITS->SaveAs(Form("%s/Positron_After_MassCut_NClsITS%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronNClsITS;		
	
	}
	
	if( hESDDalitzElectronAfterNClsITS_mc && hESDDalitzElectronAfterNClsITS_data  &&  hESDDalitzPositronAfterNClsITS_mc && hESDDalitzPositronAfterNClsITS_data ){
	
	    Double_t YMax = 0;
	    Double_t YMin = 0;
	    Double_t XMin = 0;
	    Double_t XMax = 0;
	   if( energy.CompareTo("7TeV") == 0 ){
	      
	    YMin = 1e-4;   YMax = 2;
	    XMin = -0.5;   XMax = 7.0; 
	
	  } else if( energy.CompareTo("2.76TeV") == 0 ){
	    
	    YMin = 1e-4;   YMax = 2;
	    XMin = -0.5;   XMax = 7.0; 
	    
	  } else if( energy.CompareTo("pPb_5.023TeV") == 0 ){
	    
	    YMin = 1e-4;   YMax = 2;
	    XMin = -0.5;   XMax = 7.0; 
	    
	  }
	  
	   hESDDalitzElectronPositronAfterClsITS_mc = (TH1F*)hESDDalitzElectronAfterNClsITS_mc->Clone("ESD_DalitzElectronPositron_After_ClsITS_mc");
	   hESDDalitzElectronPositronAfterClsITS_mc->Add(hESDDalitzPositronAfterNClsITS_mc,1.0);
	   
	   hESDDalitzElectronPositronAfterClsITS_data = (TH1F*)hESDDalitzElectronAfterNClsITS_data->Clone("ESD_DalitzElectronPositron_After_ClsITS_data");
	   hESDDalitzElectronPositronAfterClsITS_data->Add(hESDDalitzPositronAfterNClsITS_data,1.0);
	   
	   
	   TCanvas * canvasElectronPositronNClsITS = new TCanvas("canvasElectronPositronNClsITS","",10,10,500,500);  // gives the page size
	   canvasElectronPositronNClsITS->SetLogy(1);

	   canvasElectronPositronNClsITS->SetLeftMargin(0.16);
	   canvasElectronPositronNClsITS->cd();
	
	
	   hESDDalitzElectronPositronAfterClsITS_mc->SetStats(kFALSE);
	   hESDDalitzElectronPositronAfterClsITS_data->SetStats(kFALSE);
	
	
	   //hESDDalitzElectronPositronAfterClsITS_mc->Sumw2();
	   //hESDDalitzElectronPositronAfterClsITS_data->Sumw2();
       
	   hESDDalitzElectronPositronAfterClsITS_mc->Scale(1./hESDDalitzElectronPositronAfterClsITS_mc->Integral());
	   hESDDalitzElectronPositronAfterClsITS_data->Scale(1./hESDDalitzElectronPositronAfterClsITS_data->Integral());
	
	   DrawAutoGammaHistosTemp( hESDDalitzElectronPositronAfterClsITS_data,
						hESDDalitzElectronPositronAfterClsITS_mc, 
						"", "N_{ITS_{cls}} e^{#pm}",textYAxisITSclsEP,
						kFALSE, 15,1e-3,
						kTRUE,YMin,YMax, 
						kTRUE,XMin,XMax,1.,1.6);	
	
	   
	   
	
	canvasElectronPositronNClsITS->Update();
	canvasElectronPositronNClsITS->SaveAs(Form("%s/ElectronPositron_After_MassCut_NClsITS%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronNClsITS;		
	  
	  
	}
	
	
	if( hESDDalitzPosEleAfterDCAxy_mc && hESDDalitzPosEleAfterDCAxy_data ) {
	
	TCanvas * canvasPosEleAfterDCAxy = new TCanvas("canvasPosEleAfterDCAxy","",10,10,500,500);  // gives the page size
	canvasPosEleAfterDCAxy->cd(1);
	canvasPosEleAfterDCAxy->SetLogy(1);

	canvasPosEleAfterDCAxy->SetLeftMargin(0.17);
	canvasPosEleAfterDCAxy->SetRightMargin(0.02);
	
	TH1F* hESDDalitzPosEleAfterDCAxy_Proj_mc   = (TH1F*)hESDDalitzPosEleAfterDCAxy_mc->ProjectionX("hESDDalitzPosEleAfterDCAxy_Proj_mc");
        TH1F* hESDDalitzPosEleAfterDCAxy_Proj_data = (TH1F*)hESDDalitzPosEleAfterDCAxy_data->ProjectionX("hESDDalitzPosEleAfterDCAxy_Proj_data");
	
		
	hESDDalitzPosEleAfterDCAxy_Proj_mc->Scale(1./hESDDalitzPosEleAfterDCAxy_Proj_mc->GetEffectiveEntries());
	hESDDalitzPosEleAfterDCAxy_Proj_data->Scale(1./hESDDalitzPosEleAfterDCAxy_Proj_data->GetEffectiveEntries());
	
	
	
	DrawAutoGammaHistosTemp( hESDDalitzPosEleAfterDCAxy_Proj_mc,
						hESDDalitzPosEleAfterDCAxy_Proj_data, 
						"", "DCA_{xy} (cm)",textYAxisDCAxy,
						kTRUE,5.0,1e-5,
						kFALSE,1e-3 ,5e-1, 
						kTRUE, -0.4,0.4,1.0,1.9);	
		
	
	processSystem->Draw();
	
	canvasPosEleAfterDCAxy->Update();
	canvasPosEleAfterDCAxy->SaveAs(Form("%s/ElectronPositron_After_MassCut_DCAxy%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPosEleAfterDCAxy;		
	delete hESDDalitzPosEleAfterDCAxy_Proj_mc;
	delete hESDDalitzPosEleAfterDCAxy_Proj_data;
	
	}
      
	
	if( hESDDalitzPosEleAfterDCAz_mc && hESDDalitzPosEleAfterDCAz_data ) {
	
	TCanvas * canvasPosEleAfterDCAz = new TCanvas("canvasPosEleAfterDCAz","",10,10,500,500);  // gives the page size
	canvasPosEleAfterDCAz->SetLogy(1);
	canvasPosEleAfterDCAz->cd(1);
	

	canvasPosEleAfterDCAz->SetLeftMargin(0.17);
	canvasPosEleAfterDCAz->SetRightMargin(0.02);
	
	
	
	TH1F* hESDDalitzPosEleAfterDCAz_Proj_mc   = (TH1F*)hESDDalitzPosEleAfterDCAz_mc->ProjectionX("hESDDalitzPosEleAfterDCAz_Proj_mc");
        TH1F* hESDDalitzPosEleAfterDCAz_Proj_data = (TH1F*)hESDDalitzPosEleAfterDCAz_data->ProjectionX("hESDDalitzPosEleAfterDCAz_Proj_data");
	
	hESDDalitzPosEleAfterDCAz_Proj_mc->SetStats(kFALSE);
	hESDDalitzPosEleAfterDCAz_Proj_data->SetStats(kFALSE);
	
	
	hESDDalitzPosEleAfterDCAz_Proj_mc->Scale(1./hESDDalitzPosEleAfterDCAz_Proj_mc->GetEntries());
	hESDDalitzPosEleAfterDCAz_Proj_data->Scale(1./hESDDalitzPosEleAfterDCAz_Proj_data->GetEntries());
       
	DrawAutoGammaHistosTemp( hESDDalitzPosEleAfterDCAz_Proj_mc,
						hESDDalitzPosEleAfterDCAz_Proj_data, 
						"", "DCA_{z} (cm)",textYAxisDCAz,
						kTRUE,5.0,1e-5,
						kFALSE,1e-3 ,5e-1, 
						kTRUE, -0.9,0.9,1.0,1.9);
	
	processSystem->Draw();
		
	canvasPosEleAfterDCAz->Update();
	canvasPosEleAfterDCAz->SaveAs(Form("%s/ElectronPositron_After_MassCut_DCAz%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPosEleAfterDCAz;		
	delete hESDDalitzPosEleAfterDCAz_Proj_mc;
	delete hESDDalitzPosEleAfterDCAz_Proj_data;
	
	}
	
	  
	if( hESDDalitzElectronAfterTPCdEdxVsP_data ) {
	
	TCanvas * canvasElectronTPCdEdxVsP_data = new TCanvas("canvasElectronTPCdEdxVsP_data","",10,10,500,500);  // gives the page size
	canvasElectronTPCdEdxVsP_data->SetLeftMargin(0.13);
	canvasElectronTPCdEdxVsP_data->SetRightMargin(0.12);
	canvasElectronTPCdEdxVsP_data->SetTopMargin(0.05);
	canvasElectronTPCdEdxVsP_data->SetBottomMargin(0.1);
	
	canvasElectronTPCdEdxVsP_data->cd();
	canvasElectronTPCdEdxVsP_data->SetLogx(1);
	canvasElectronTPCdEdxVsP_data->SetLogz(1);
	hESDDalitzElectronAfterTPCdEdxVsP_data->SetStats(kFALSE);
	hESDDalitzElectronAfterTPCdEdxVsP_data->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hESDDalitzElectronAfterTPCdEdxVsP_data,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e^{-}} (#sigma)","Data",
						 kTRUE,-4.1, 5.1,
						 kTRUE, 0.,200,0.9);
	
	
	canvasElectronTPCdEdxVsP_data->Update();
	canvasElectronTPCdEdxVsP_data->SaveAs(Form("%s/Electron_After_MassCut_TPCdEdxVsP_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronTPCdEdxVsP_data;
	
	}
	
	
	if( hESDDalitzElectronAfterTPCdEdxVsP_mc ) {
	
	TCanvas * canvasElectronTPCdEdxVsP_mc = new TCanvas("canvasElectronTPCdEdxVsP_mc","",10,10,500,500);  // gives the page size
	canvasElectronTPCdEdxVsP_mc->SetLeftMargin(0.13);
	canvasElectronTPCdEdxVsP_mc->SetRightMargin(0.12);
	canvasElectronTPCdEdxVsP_mc->SetTopMargin(0.05);
	canvasElectronTPCdEdxVsP_mc->SetBottomMargin(0.1);
	
	canvasElectronTPCdEdxVsP_mc->cd();
	canvasElectronTPCdEdxVsP_mc->SetLogz(1);
	canvasElectronTPCdEdxVsP_mc->SetLogx(1);
	hESDDalitzElectronAfterTPCdEdxVsP_mc->SetStats(kFALSE);
	hESDDalitzElectronAfterTPCdEdxVsP_mc->SetTitle("");
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterTPCdEdxVsP_mc,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e^{-}} (#sigma)","MC",
						 kTRUE,-4.1,5.1,
						 kTRUE, 0.,200,0.9);
	
	
	canvasElectronTPCdEdxVsP_mc->Update();
	canvasElectronTPCdEdxVsP_mc->SaveAs(Form("%s/Electron_After_MassCut_TPCdEdxVsP_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronTPCdEdxVsP_mc;		
	
	}
	
	
	if( hESDDalitzPositronAfterTPCdEdxVsP_data ) {
	
	TCanvas * canvasPositronTPCdEdxVsP_data = new TCanvas("canvasPositronTPCdEdxVsP_data","",10,10,500,500);  // gives the page size
	canvasPositronTPCdEdxVsP_data->SetLeftMargin(0.13);
	canvasPositronTPCdEdxVsP_data->SetRightMargin(0.12);
	canvasPositronTPCdEdxVsP_data->SetTopMargin(0.05);
	canvasPositronTPCdEdxVsP_data->SetBottomMargin(0.1);
	
	canvasPositronTPCdEdxVsP_data->cd();
	canvasPositronTPCdEdxVsP_data->SetLogz(1);
	canvasPositronTPCdEdxVsP_data->SetLogx(1);
	hESDDalitzPositronAfterTPCdEdxVsP_data->SetStats(kFALSE);
	hESDDalitzPositronAfterTPCdEdxVsP_data->SetTitle("");
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterTPCdEdxVsP_data,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e^{+}} (#sigma)","Data",
						 kTRUE,-4.1, 5.1,
						 kTRUE, 0.,200,0.9);
	
	
	canvasPositronTPCdEdxVsP_data->Update();
	canvasPositronTPCdEdxVsP_data->SaveAs(Form("%s/Positron_After_MassCut_TPCdEdxVsP_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronTPCdEdxVsP_data;
	
	}
	
	if( hESDDalitzPositronAfterTPCdEdxVsP_mc ) {
	
	TCanvas * canvasPositronTPCdEdxVsP_mc = new TCanvas("canvasPositronTPCdEdxVsP_mc","",10,10,500,500);  // gives the page size
	canvasPositronTPCdEdxVsP_mc->SetLeftMargin(0.13);
	canvasPositronTPCdEdxVsP_mc->SetRightMargin(0.12);
	canvasPositronTPCdEdxVsP_mc->SetTopMargin(0.05);
	canvasPositronTPCdEdxVsP_mc->SetBottomMargin(0.1);
	
	canvasPositronTPCdEdxVsP_mc->cd();
	canvasPositronTPCdEdxVsP_mc->SetLogz(1);
	canvasPositronTPCdEdxVsP_mc->SetLogx(1);
	hESDDalitzPositronAfterTPCdEdxVsP_mc->SetStats(kFALSE);
	hESDDalitzPositronAfterTPCdEdxVsP_mc->SetTitle("");
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterTPCdEdxVsP_mc,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e^{+}} (#sigma)","MC",
						 kTRUE,-4.1, 5.1,
						 kTRUE, 0.,200,0.9);
	
	
	canvasPositronTPCdEdxVsP_mc->Update();
	canvasPositronTPCdEdxVsP_mc->SaveAs(Form("%s/Positron_After_MassCut_TPCdEdxVsP_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronTPCdEdxVsP_mc;
	
	}
	
	
	if( hESDDalitzElectronAfterTPCdEdxVsP_data ) {
	  
	
	
	TH2F* hESDDalitzElectronPositronAfterTPCdEdxVsP_data = (TH2F*)hESDDalitzElectronAfterTPCdEdxVsP_data->Clone("hESDDalitzElectronPositronAfterTPCdEdxVsP_data");
	
	hESDDalitzElectronPositronAfterTPCdEdxVsP_data->Add(hESDDalitzPositronAfterTPCdEdxVsP_data,1.0);
	
	TCanvas * canvasElectronPositronTPCdEdxVsP_data = new TCanvas("canvasElectronPositronTPCdEdxVsP_data","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxVsP_data->SetLeftMargin(0.1);
	canvasElectronPositronTPCdEdxVsP_data->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxVsP_data->SetTopMargin(0.03);
	canvasElectronPositronTPCdEdxVsP_data->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxVsP_data->SetTickx();
        canvasElectronPositronTPCdEdxVsP_data->SetTicky();
      
	canvasElectronPositronTPCdEdxVsP_data->cd();
	canvasElectronPositronTPCdEdxVsP_data->SetLogx(1);
	canvasElectronPositronTPCdEdxVsP_data->SetLogz(1);
	hESDDalitzElectronPositronAfterTPCdEdxVsP_data->SetStats(kFALSE);
	hESDDalitzElectronPositronAfterTPCdEdxVsP_data->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hESDDalitzElectronPositronAfterTPCdEdxVsP_data,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e} (#sigma)","",
						 kTRUE,-4.1, 5.1,
						 kTRUE, 0.1,10.,1.11,0.9);
	
	 
	
	 
	
	TLatex *process03 = new TLatex(0.50,0.905,collisionSystem.Data());
        process03->SetNDC();
        process03->SetTextColor(1);
        process03->SetTextSize(0.040);
	process03->SetTextFont(42);
        process03->Draw();
	
	TLatex *DateTex03 = new TLatex(0.50,0.855,"Data");
	DateTex03->SetNDC();
	DateTex03->SetTextColor(1);
	DateTex03->SetTextSize(0.04);
	DateTex03->SetTextFont(42);
	DateTex03->Draw();
	
	
	canvasElectronPositronTPCdEdxVsP_data->Update();
	canvasElectronPositronTPCdEdxVsP_data->SaveAs(Form("%s/ElectronPositron_After_MassCut_TPCdEdxVsP_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxVsP_data;
	
	
	}
	
	
	if( hESDDalitzElectronAfterTPCdEdxVsP_mc ) {
	  
	
	
	TH2F* hESDDalitzElectronPositronAfterTPCdEdxVsP_mc = (TH2F*)hESDDalitzElectronAfterTPCdEdxVsP_mc->Clone("hESDDalitzElectronPositronAfterTPCdEdxVsP_mc");
	
	hESDDalitzElectronPositronAfterTPCdEdxVsP_mc->Add(hESDDalitzPositronAfterTPCdEdxVsP_mc,1.0);
	
	TCanvas * canvasElectronPositronTPCdEdxVsP_mc = new TCanvas("canvasElectronPositronTPCdEdxVsP_mc","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxVsP_mc->SetLeftMargin(0.1);
	canvasElectronPositronTPCdEdxVsP_mc->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxVsP_mc->SetTopMargin(0.03);
	canvasElectronPositronTPCdEdxVsP_mc->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxVsP_mc->SetTickx();
        canvasElectronPositronTPCdEdxVsP_mc->SetTicky();
      
	canvasElectronPositronTPCdEdxVsP_mc->cd();
	canvasElectronPositronTPCdEdxVsP_mc->SetLogx(1);
	canvasElectronPositronTPCdEdxVsP_mc->SetLogz(1);
	hESDDalitzElectronPositronAfterTPCdEdxVsP_mc->SetStats(kFALSE);
	hESDDalitzElectronPositronAfterTPCdEdxVsP_mc->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hESDDalitzElectronPositronAfterTPCdEdxVsP_mc,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e} (#sigma)","",
						 kTRUE,-4.1, 5.1,
						 kTRUE, 0.1,10.,1.11,0.9);
	
	
	 
	
	TLatex *process03 = new TLatex(0.50,0.905,collisionSystem.Data());
        process03->SetNDC();
        process03->SetTextColor(1);
        process03->SetTextSize(0.040);
	process03->SetTextFont(42);
        process03->Draw();
	
	TLatex *DateTex03 = new TLatex(0.50,0.855,"MC");
	DateTex03->SetNDC();
	DateTex03->SetTextColor(1);
	DateTex03->SetTextSize(0.04);
	DateTex03->SetTextFont(42);
	DateTex03->Draw();
	
	
	canvasElectronPositronTPCdEdxVsP_mc->Update();
	canvasElectronPositronTPCdEdxVsP_mc->SaveAs(Form("%s/ElectronPositron_After_MassCut_TPCdEdxVsP_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxVsP_mc;
	
	
	}
	
	
	
	
	if( hESDDalitzElectronAfterTPCdEdxSignalVsP_data ) {
	
	TCanvas * canvasElectronTPCdEdxSignalVsP_data = new TCanvas("canvasElectronTPCdEdxSignalVsP_data","",10,10,500,500);  // gives the page size
	canvasElectronTPCdEdxSignalVsP_data->SetLeftMargin(0.13);
	canvasElectronTPCdEdxSignalVsP_data->SetRightMargin(0.12);
	canvasElectronTPCdEdxSignalVsP_data->SetTopMargin(0.05);
	canvasElectronTPCdEdxSignalVsP_data->SetBottomMargin(0.1);
	
	canvasElectronTPCdEdxSignalVsP_data->cd();
	canvasElectronTPCdEdxSignalVsP_data->SetLogz(1);
	canvasElectronTPCdEdxSignalVsP_data->SetLogx(1);
	hESDDalitzElectronAfterTPCdEdxSignalVsP_data->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterTPCdEdxSignalVsP_data,
						 "", "p (GeV/c)", "dE/dx signal in TPC (arb. units)","Data",
						 kTRUE, 40.,140.,
						 kTRUE, 40.,140,0.90);
	 
	
	
	canvasElectronTPCdEdxSignalVsP_data->Update();
	canvasElectronTPCdEdxSignalVsP_data->SaveAs(Form("%s/Electron_After_MassCut_TPCdEdxSignalVsP_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronTPCdEdxSignalVsP_data;
	
	
	}
	
	if( hESDDalitzElectronAfterTPCdEdxSignalVsP_mc ) {
	
	
	TCanvas * canvasElectronTPCdEdxSignalVsP_mc = new TCanvas("canvasElectronTPCdEdxSignalVsP_mc","",10,10,500,500);  // gives the page size
	canvasElectronTPCdEdxSignalVsP_mc->SetLeftMargin(0.13);
	canvasElectronTPCdEdxSignalVsP_mc->SetRightMargin(0.12);
	canvasElectronTPCdEdxSignalVsP_mc->SetTopMargin(0.05);
	canvasElectronTPCdEdxSignalVsP_mc->SetBottomMargin(0.1);
	canvasElectronTPCdEdxSignalVsP_mc->cd();
	canvasElectronTPCdEdxSignalVsP_mc->SetLogz(1);
	canvasElectronTPCdEdxSignalVsP_mc->SetLogx(1);
	hESDDalitzElectronAfterTPCdEdxSignalVsP_mc->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterTPCdEdxSignalVsP_mc,
						 "", "p (GeV/c)", "dE/dx signal in TPC (arb. units)","MC",
						 kTRUE, 40., 140.,
						 kTRUE, 40., 140,0.9);
	 
	
	
	canvasElectronTPCdEdxSignalVsP_mc->Update();
	canvasElectronTPCdEdxSignalVsP_mc->SaveAs(Form("%s/Electron_After_MassCut_TPCdEdxSignalVsP_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronTPCdEdxSignalVsP_mc;
	
	}
	
	
	if( hESDDalitzPositronAfterTPCdEdxSignalVsP_data ) {
	
	TCanvas * canvasPositronTPCdEdxSignalVsP_data = new TCanvas("canvasPositronTPCdEdxSignalVsP_data","",10,10,500,500);  // gives the page size
	canvasPositronTPCdEdxSignalVsP_data->SetLeftMargin(0.13);
	canvasPositronTPCdEdxSignalVsP_data->SetRightMargin(0.12);
	canvasPositronTPCdEdxSignalVsP_data->SetTopMargin(0.05);
	canvasPositronTPCdEdxSignalVsP_data->SetBottomMargin(0.1);
	canvasPositronTPCdEdxSignalVsP_data->cd();
	canvasPositronTPCdEdxSignalVsP_data->SetLogz(1);
	canvasPositronTPCdEdxSignalVsP_data->SetLogx(1);
	hESDDalitzPositronAfterTPCdEdxSignalVsP_data->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterTPCdEdxSignalVsP_data,
						 "", "p (GeV/c)", "dE/dx signal in TPC (arb. units)","Data",
						 kTRUE, 40.,140.,
						 kTRUE, 40.,140,0.9);
	
	
	canvasPositronTPCdEdxSignalVsP_data->Update();
	canvasPositronTPCdEdxSignalVsP_data->SaveAs(Form("%s/Positron_After_MassCut_TPCdEdxSignalVsP_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronTPCdEdxSignalVsP_data;
	
	}
	
	
	if( hESDDalitzPositronAfterTPCdEdxSignalVsP_mc ) {
	
	
	TCanvas * canvasPositronTPCdEdxSignalVsP_mc = new TCanvas("canvasPositronTPCdEdxSignalVsP_mc","",10,10,500,500);  // gives the page size
	canvasPositronTPCdEdxSignalVsP_mc->SetLeftMargin(0.13);
	canvasPositronTPCdEdxSignalVsP_mc->SetRightMargin(0.12);
	canvasPositronTPCdEdxSignalVsP_mc->SetTopMargin(0.05);
	canvasPositronTPCdEdxSignalVsP_mc->SetBottomMargin(0.1);
	canvasPositronTPCdEdxSignalVsP_mc->cd();
	canvasPositronTPCdEdxSignalVsP_mc->SetLogz(1);
	canvasPositronTPCdEdxSignalVsP_mc->SetLogx(1);
	hESDDalitzPositronAfterTPCdEdxSignalVsP_mc->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterTPCdEdxSignalVsP_mc,
						 "", "p (GeV/c)", "dE/dx signal in TPC (arb. units)","MC",
						 kTRUE, 40.,140.,
						 kTRUE, 40.,140,0.9);
	 
	
	
	canvasPositronTPCdEdxSignalVsP_mc->Update();
	canvasPositronTPCdEdxSignalVsP_mc->SaveAs(Form("%s/Positron_After_MassCut_TPCdEdxSignalVsP_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronTPCdEdxSignalVsP_mc;
	
	}
	
	
	if( hESDDalitzElectronAfterTPCdEdxSignalVsP_data ) {
	
	
	TH2F* hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_data = (TH2F*)hESDDalitzElectronAfterTPCdEdxSignalVsP_data->Clone("hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_data");
	
	hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_data->Add(hESDDalitzPositronAfterTPCdEdxSignalVsP_data,1.0);
	
	TCanvas * canvasElectronPositronTPCdEdxSignalVsP_data = new TCanvas("canvasElectronPositronTPCdEdxSignalVsP_data","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxSignalVsP_data->SetLeftMargin(0.13);
	canvasElectronPositronTPCdEdxSignalVsP_data->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxSignalVsP_data->SetTopMargin(0.02);
	canvasElectronPositronTPCdEdxSignalVsP_data->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxSignalVsP_data->cd();
	canvasElectronPositronTPCdEdxSignalVsP_data->SetLogz(1);
	canvasElectronPositronTPCdEdxSignalVsP_data->SetLogx(1);
	hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_data->SetStats(kFALSE);
	hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_data->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_data,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x}","",
						 kTRUE,0.0, 180.0,
						 kTRUE, 0.,10.);
	 
	 
	TLatex *process03 = new TLatex(0.50,0.905,collisionSystem.Data());
        process03->SetNDC();
        process03->SetTextColor(1);
        process03->SetTextSize(0.040);
	process03->SetTextFont(42);
        process03->Draw();
	
	TLatex *DateTex03 = new TLatex(0.50,0.855,"Data");
	DateTex03->SetNDC();
	DateTex03->SetTextColor(1);
	DateTex03->SetTextSize(0.04);
	DateTex03->SetTextFont(42);
	DateTex03->Draw();
	
	TLatex *DateTexElectron = new TLatex(0.30,0.47,"e");
	DateTexElectron->SetNDC();
	DateTexElectron->SetTextColor(1);
	DateTexElectron->SetTextSize(0.04);
	DateTexElectron->SetTextFont(42);
	DateTexElectron->Draw();
	
	 
	
	
	
	canvasElectronPositronTPCdEdxSignalVsP_data->Update();
	canvasElectronPositronTPCdEdxSignalVsP_data->SaveAs(Form("%s/ElectronPositron_After_MassCut_TPCdEdxSignalVsP_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxSignalVsP_data;
	
	}
	
	
	if( hTPCdEdxbefore_data ) {
	  
	TCanvas * canvasElectronPositronTPCdEdxVsPBefore_data = new TCanvas("canvasElectronPositronTPCdEdxVsPBefore_data","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxVsPBefore_data->SetLeftMargin(0.13);
	canvasElectronPositronTPCdEdxVsPBefore_data->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxVsPBefore_data->SetTopMargin(0.02);
	canvasElectronPositronTPCdEdxVsPBefore_data->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxVsPBefore_data->cd();
	canvasElectronPositronTPCdEdxVsPBefore_data->SetLogz(1);
	canvasElectronPositronTPCdEdxVsPBefore_data->SetLogx(1);
	hTPCdEdxbefore_data->SetStats(kFALSE);
	hTPCdEdxbefore_data->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hTPCdEdxbefore_data,
						 "",pLabel.Data(), "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e} (#sigma)","",
						 kTRUE,-10.0, 10.0,
						 kTRUE, 0.1,10.);
	 

	 TPaveText *paveTextElectronPositronTPCdEdxVsPBefore = new TPaveText(0.55,0.90,0.89,0.98,"b1NDC");
	 paveTextElectronPositronTPCdEdxVsPBefore->SetBorderSize(1);
	 paveTextElectronPositronTPCdEdxVsPBefore->SetFillColor(kWhite);
	 paveTextElectronPositronTPCdEdxVsPBefore->AddText(collisionSystem.Data());
	 paveTextElectronPositronTPCdEdxVsPBefore->AddText("This thesis");
	 
	 paveTextElectronPositronTPCdEdxVsPBefore->Draw();
	 
		
	
	
	canvasElectronPositronTPCdEdxVsPBefore_data->Update();
	canvasElectronPositronTPCdEdxVsPBefore_data->SaveAs(Form("%s/ElectronPositron_TPCdEdxVsP_Before_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxVsPBefore_data;
	
	cout<<"Pedrito lleno el material"<<endl;
	  
	  
	  
	}
	
	if( hTPCdEdxbefore_mc ) {
	  
	TCanvas * canvasElectronPositronTPCdEdxVsPBefore_mc = new TCanvas("canvasElectronPositronTPCdEdxVsPBefore_mc","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxVsPBefore_mc->SetLeftMargin(0.13);
	canvasElectronPositronTPCdEdxVsPBefore_mc->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxVsPBefore_mc->SetTopMargin(0.02);
	canvasElectronPositronTPCdEdxVsPBefore_mc->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxVsPBefore_mc->cd();
	canvasElectronPositronTPCdEdxVsPBefore_mc->SetLogz(1);
	canvasElectronPositronTPCdEdxVsPBefore_mc->SetLogx(1);
	hTPCdEdxbefore_mc->SetStats(kFALSE);
	hTPCdEdxbefore_mc->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hTPCdEdxbefore_mc,
						 "",pLabel.Data(), "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e} (#sigma)","",
						 kTRUE,-10.0, 10.0,
						 kTRUE, 0.1,10.);
	 

	
	
	 TPaveText *paveTextElectronPositronTPCdEdxVsPBefore = new TPaveText(0.55,0.90,0.89,0.98,"b1NDC");
	 paveTextElectronPositronTPCdEdxVsPBefore->SetBorderSize(1);
	 paveTextElectronPositronTPCdEdxVsPBefore->SetFillColor(kWhite);
	 paveTextElectronPositronTPCdEdxVsPBefore->AddText(collisionSystem.Data());
	 paveTextElectronPositronTPCdEdxVsPBefore->AddText("MC");
	 
	 
	 paveTextElectronPositronTPCdEdxVsPBefore->Draw();
	  
		
	
	
	canvasElectronPositronTPCdEdxVsPBefore_mc->Update();
	canvasElectronPositronTPCdEdxVsPBefore_mc->SaveAs(Form("%s/ElectronPositron_TPCdEdxVsP_Before_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxVsPBefore_mc;
	
	cout<<"Pedrito lleno el material"<<endl;
	  
	  
	  
	}
	
	
	if( hTPCdEdxSigbefore_data ) {
	
	
	cout<<"Llego a preparar el plot"<<endl;  
	
	
	
	cout<<"No entiendo que pasa aqui"<<endl;
	
	
	TCanvas * canvasElectronPositronTPCdEdxSignalVsPBefore_data = new TCanvas("canvasElectronPositronTPCdEdxSignalVsPBefore_data","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxSignalVsPBefore_data->SetLeftMargin(0.13);
	canvasElectronPositronTPCdEdxSignalVsPBefore_data->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxSignalVsPBefore_data->SetTopMargin(0.02);
	canvasElectronPositronTPCdEdxSignalVsPBefore_data->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxSignalVsPBefore_data->cd();
	canvasElectronPositronTPCdEdxSignalVsPBefore_data->SetLogz(1);
	canvasElectronPositronTPCdEdxSignalVsPBefore_data->SetLogx(1);
	hTPCdEdxSigbefore_data->SetStats(kFALSE);
	hTPCdEdxSigbefore_data->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hTPCdEdxSigbefore_data,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x}","",
						 kTRUE,0.0, 200.0,
						 kTRUE, 0.1,10.);
	 
	 
	TLatex *process03 = new TLatex(0.53,0.925,collisionSystem.Data());
        process03->SetNDC();
        process03->SetTextColor(1);
        process03->SetTextSize(0.040);
	process03->SetTextFont(42);
        process03->Draw();
	
	TLatex *DateTex03 = new TLatex(0.53,0.875,"Data");
	DateTex03->SetNDC();
	DateTex03->SetTextColor(1);
	DateTex03->SetTextSize(0.04);
	DateTex03->SetTextFont(42);
	DateTex03->Draw();
	
	
	
	TLatex *DateTexElectron = new TLatex(0.28,0.47,"e");
	DateTexElectron->SetNDC();
	DateTexElectron->SetTextColor(1);
	DateTexElectron->SetTextSize(0.04);
	DateTexElectron->SetTextFont(42);
	DateTexElectron->Draw();
	
	TLatex *DateTexKaon = new TLatex(0.33,0.755,"K");
	DateTexKaon->SetNDC();
	DateTexKaon->SetTextColor(1);
	DateTexKaon->SetTextSize(0.04);
	DateTexKaon->SetTextFont(42);
	DateTexKaon->Draw();
	
	TLatex *DateTexProton = new TLatex(0.43,0.755,"p");
	DateTexProton->SetNDC();
	DateTexProton->SetTextColor(1);
	DateTexProton->SetTextSize(0.04);
	DateTexProton->SetTextFont(42);
	DateTexProton->Draw();
	
	TLatex *DateTexPion = new TLatex(0.25,0.36,"#pi");
	DateTexPion->SetNDC();
	DateTexPion->SetTextColor(1);
	DateTexPion->SetTextSize(0.04);
	DateTexPion->SetTextFont(42);
	DateTexPion->Draw();
	
		
	
	
	canvasElectronPositronTPCdEdxSignalVsPBefore_data->Update();
	canvasElectronPositronTPCdEdxSignalVsPBefore_data->SaveAs(Form("%s/ElectronPositron_TPCdEdxSignalVsP_Before_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxSignalVsPBefore_data;
	
	cout<<"Pedrito lleno el material"<<endl;
	
	}
	
	
	if( hTPCEledEdxSigafter_data ) {
	
	
	cout<<"Llego a preparar el plot"<<endl;  
	
	
	
	cout<<"No entiendo que pasa aqui"<<endl;
	
	
	TCanvas * canvasElectronPositronTPCdEdxSignalVsPAfter_data = new TCanvas("canvasElectronPositronTPCdEdxSignalVsPAfter_data","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxSignalVsPAfter_data->SetLeftMargin(0.13);
	canvasElectronPositronTPCdEdxSignalVsPAfter_data->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxSignalVsPAfter_data->SetTopMargin(0.02);
	canvasElectronPositronTPCdEdxSignalVsPAfter_data->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxSignalVsPAfter_data->cd();
	canvasElectronPositronTPCdEdxSignalVsPAfter_data->SetLogz(1);
	canvasElectronPositronTPCdEdxSignalVsPAfter_data->SetLogx(1);
	hTPCEledEdxSigafter_data->SetStats(kFALSE);
	hTPCEledEdxSigafter_data->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hTPCEledEdxSigafter_data,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x}","",
						 kTRUE,0.0, 200.0,
						 kTRUE, 0.1,10.);
	 
	 
	TLatex *process03 = new TLatex(0.53,0.925,collisionSystem.Data());
        process03->SetNDC();
        process03->SetTextColor(1);
        process03->SetTextSize(0.040);
	process03->SetTextFont(42);
        process03->Draw();
	
	TLatex *DateTex03 = new TLatex(0.53,0.875,"Data");
	DateTex03->SetNDC();
	DateTex03->SetTextColor(1);
	DateTex03->SetTextSize(0.04);
	DateTex03->SetTextFont(42);
	DateTex03->Draw();
	
	
	
	TLatex *DateTexElectron = new TLatex(0.28,0.47,"e");
	DateTexElectron->SetNDC();
	DateTexElectron->SetTextColor(1);
	DateTexElectron->SetTextSize(0.04);
	DateTexElectron->SetTextFont(42);
	DateTexElectron->Draw();
	
	TLatex *DateTexKaon = new TLatex(0.37,0.555,"K");
	DateTexKaon->SetNDC();
	DateTexKaon->SetTextColor(1);
	DateTexKaon->SetTextSize(0.04);
	DateTexKaon->SetTextFont(42);
	DateTexKaon->Draw();
	
	TLatex *DateTexProton = new TLatex(0.47,0.555,"p");
	DateTexProton->SetNDC();
	DateTexProton->SetTextColor(1);
	DateTexProton->SetTextSize(0.04);
	DateTexProton->SetTextFont(42);
	DateTexProton->Draw();
	
	TLatex *DateTexPion = new TLatex(0.25,0.36,"#pi");
	DateTexPion->SetNDC();
	DateTexPion->SetTextColor(1);
	DateTexPion->SetTextSize(0.04);
	DateTexPion->SetTextFont(42);
	DateTexPion->Draw();
	
		
	
	
	canvasElectronPositronTPCdEdxSignalVsPAfter_data->Update();
	canvasElectronPositronTPCdEdxSignalVsPAfter_data->SaveAs(Form("%s/ElectronPositron_TPCdEdxSignalVsP_After_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxSignalVsPAfter_data;
	
	cout<<"Pedrito lleno el material"<<endl;
	
	}
	
	
	if( hTPCEledEdxafter_data ) {
	
	
	cout<<"Llego a preparar el plot"<<endl;  
	
	
	
	cout<<"No entiendo que pasa aqui"<<endl;
	
	
	TCanvas * canvasElectronPositronTPCdEdxVsPAfter_data = new TCanvas("canvasElectronPositronTPCdEdxVsPAfter_data","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxVsPAfter_data->SetLeftMargin(0.13);
	canvasElectronPositronTPCdEdxVsPAfter_data->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxVsPAfter_data->SetTopMargin(0.02);
	canvasElectronPositronTPCdEdxVsPAfter_data->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxVsPAfter_data->cd();
	canvasElectronPositronTPCdEdxVsPAfter_data->SetLogz(1);
	canvasElectronPositronTPCdEdxVsPAfter_data->SetLogx(1);
	hTPCEledEdxafter_data->SetStats(kFALSE);
	hTPCEledEdxafter_data->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hTPCEledEdxafter_data,
						 "",pLabel.Data(), "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e} (#sigma)","",
						 kTRUE,-6.5, 6.5,
						 kTRUE, 0.1,10.);
	 
		
	
	 TPaveText *paveTextElectronPositronTPCdEdxVsPAfter = new TPaveText(0.55,0.90,0.89,0.98,"b1NDC");
	 paveTextElectronPositronTPCdEdxVsPAfter->SetBorderSize(1);
	 paveTextElectronPositronTPCdEdxVsPAfter->SetFillColor(kWhite);
	 paveTextElectronPositronTPCdEdxVsPAfter->AddText(collisionSystem.Data());
	 paveTextElectronPositronTPCdEdxVsPAfter->AddText("This thesis");
	 
	 paveTextElectronPositronTPCdEdxVsPAfter->Draw();
	 
	
	
	
	canvasElectronPositronTPCdEdxVsPAfter_data->Update();
	canvasElectronPositronTPCdEdxVsPAfter_data->SaveAs(Form("%s/ElectronPositron_TPCdEdxVsP_After_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxVsPAfter_data;
	
	cout<<"Pedrito lleno el material"<<endl;
	
	}
	
	if( hTPCEledEdxafter_mc ) {
	
	
	cout<<"Llego a preparar el plot"<<endl;  
	
	
	
	cout<<"No entiendo que pasa aqui"<<endl;
	
	
	TCanvas * canvasElectronPositronTPCdEdxVsPAfter_mc = new TCanvas("canvasElectronPositronTPCdEdxVsPAfter_mc","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxVsPAfter_mc->SetLeftMargin(0.13);
	canvasElectronPositronTPCdEdxVsPAfter_mc->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxVsPAfter_mc->SetTopMargin(0.02);
	canvasElectronPositronTPCdEdxVsPAfter_mc->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxVsPAfter_mc->cd();
	canvasElectronPositronTPCdEdxVsPAfter_mc->SetLogz(1);
	canvasElectronPositronTPCdEdxVsPAfter_mc->SetLogx(1);
	hTPCEledEdxafter_mc->SetStats(kFALSE);
	hTPCEledEdxafter_mc->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hTPCEledEdxafter_mc,
						 "",pLabel.Data(), "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e} (#sigma)","",
						 kTRUE,-6.5, 6.5,
						 kTRUE, 0.1,10.);
	 
		
	
	 TPaveText *paveTextElectronPositronTPCdEdxVsPAfter = new TPaveText(0.55,0.90,0.89,0.98,"b1NDC");
	 paveTextElectronPositronTPCdEdxVsPAfter->SetBorderSize(1);
	 paveTextElectronPositronTPCdEdxVsPAfter->SetFillColor(kWhite);
	 paveTextElectronPositronTPCdEdxVsPAfter->AddText(collisionSystem.Data());
	 paveTextElectronPositronTPCdEdxVsPAfter->AddText("MC");
	 
	 paveTextElectronPositronTPCdEdxVsPAfter->Draw();
	 
	
	
	
	canvasElectronPositronTPCdEdxVsPAfter_mc->Update();
	canvasElectronPositronTPCdEdxVsPAfter_mc->SaveAs(Form("%s/ElectronPositron_TPCdEdxVsP_After_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxVsPAfter_mc;
	
	cout<<"Pedrito lleno el material"<<endl;
	
	}
	
	
	
	if( hESDDalitzElectronAfterTPCdEdxSignalVsP_mc ) {
	
	
	TH2F* hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_mc = (TH2F*)hESDDalitzElectronAfterTPCdEdxSignalVsP_mc->Clone("hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_mc");
	
	hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_mc->Add(hESDDalitzPositronAfterTPCdEdxSignalVsP_mc,1.0);
	
	TCanvas * canvasElectronPositronTPCdEdxSignalVsP_mc = new TCanvas("canvasElectronPositronTPCdEdxSignalVsP_mc","",10,10,500,500);  // gives the page size
	canvasElectronPositronTPCdEdxSignalVsP_mc->SetLeftMargin(0.13);
	canvasElectronPositronTPCdEdxSignalVsP_mc->SetRightMargin(0.11);
	canvasElectronPositronTPCdEdxSignalVsP_mc->SetTopMargin(0.02);
	canvasElectronPositronTPCdEdxSignalVsP_mc->SetBottomMargin(0.11);
	canvasElectronPositronTPCdEdxSignalVsP_mc->cd();
	canvasElectronPositronTPCdEdxSignalVsP_mc->SetLogz(1);
	canvasElectronPositronTPCdEdxSignalVsP_mc->SetLogx(1);
	hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_mc->SetStats(kFALSE);
	hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_mc->SetTitle("");
       
	
	DrawAutoGammaHisto2D(	hESDDalitzElectronPositronAfterTPCdEdxSignalVsP_mc,
						 "", "p (GeV/c)", "dE/dx signal in TPC (arb. units)","",
						 kTRUE,0.0, 180.0,
						 kTRUE, 0.,10.);
	 
	 
	TLatex *process03 = new TLatex(0.50,0.905,collisionSystem.Data());
        process03->SetNDC();
        process03->SetTextColor(1);
        process03->SetTextSize(0.040);
	process03->SetTextFont(42);
        process03->Draw();
	
	TLatex *DateTex03 = new TLatex(0.50,0.855,"MC");
	DateTex03->SetNDC();
	DateTex03->SetTextColor(1);
	DateTex03->SetTextSize(0.04);
	DateTex03->SetTextFont(42);
	DateTex03->Draw();
	
	
	 
	
	
	
	canvasElectronPositronTPCdEdxSignalVsP_mc->Update();
	canvasElectronPositronTPCdEdxSignalVsP_mc->SaveAs(Form("%s/ElectronPositron_After_MassCut_TPCdEdxSignalVsP_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronPositronTPCdEdxSignalVsP_mc;
	
	}
	
	
	
	
	
	
	if( hESDDalitzElectronAfterTPCdEdxVsEta_data ) {
	
	TCanvas * canvasElectronTPCdEdxVsEta_data = new TCanvas("canvasElectronTPCdEdxVsEta_data","",10,10,500,500);  // gives the page size
	canvasElectronTPCdEdxVsEta_data->SetLeftMargin(0.16);
	canvasElectronTPCdEdxVsEta_data->cd();
	hESDDalitzElectronAfterTPCdEdxVsEta_data->SetStats(kFALSE);
	canvasElectronTPCdEdxVsEta_data->SetLogz(1);
	
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterTPCdEdxVsEta_data,
						 "", "#eta", "dE/dx","Data",
						 kTRUE,-10., 10.,
						 kTRUE, -2.,2);
	
	
	canvasElectronTPCdEdxVsEta_data->Update();
	canvasElectronTPCdEdxVsEta_data->SaveAs(Form("%s/Electron_After_MassCut_TPCdEdxVsEta_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronTPCdEdxVsEta_data;
	
	}
	
	
	if( hESDDalitzElectronAfterTPCdEdxVsEta_mc ) {
	
	TCanvas * canvasElectronTPCdEdxVsEta_mc = new TCanvas("canvasElectronTPCdEdxVsEta_mc","",10,10,500,500);  // gives the page size
	canvasElectronTPCdEdxVsEta_mc->SetLeftMargin(0.16);
	canvasElectronTPCdEdxVsEta_mc->cd();
	hESDDalitzElectronAfterTPCdEdxVsEta_mc->SetStats(kFALSE);
	canvasElectronTPCdEdxVsEta_mc->SetLogz(1);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterTPCdEdxVsEta_mc,
						 "", "#eta", "dE/dx","MC",
						 kTRUE,-10., 10.,
						 kTRUE, -2.,2);
	
	
	canvasElectronTPCdEdxVsEta_mc->Update();
	canvasElectronTPCdEdxVsEta_mc->SaveAs(Form("%s/Electron_After_MassCut_TPCdEdxVsEta_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronTPCdEdxVsEta_mc;
	
	}
	
	
	if( hESDDalitzPositronAfterTPCdEdxVsEta_data ) {
	
	
	TCanvas * canvasPositronTPCdEdxVsEta_data = new TCanvas("canvasPositronTPCdEdxVsEta_data","",10,10,500,500);  // gives the page size
	canvasPositronTPCdEdxVsEta_data->SetLeftMargin(0.16);
	canvasPositronTPCdEdxVsEta_data->cd();
	hESDDalitzPositronAfterTPCdEdxVsEta_data->SetStats(kFALSE);
	canvasPositronTPCdEdxVsEta_data->SetLogz(1);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterTPCdEdxVsEta_data,
						 "", "#eta", "dE/dx","Data",
						 kTRUE,-10., 10.,
						 kTRUE, -2.,2);
	
	
	canvasPositronTPCdEdxVsEta_data->Update();
	canvasPositronTPCdEdxVsEta_data->SaveAs(Form("%s/Positron_After_MassCut_TPCdEdxVsEta_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronTPCdEdxVsEta_data;
	
	}
	
	
	if( hESDDalitzPositronAfterTPCdEdxVsEta_mc ) {
	
	TCanvas * canvasPositronTPCdEdxVsEta_mc = new TCanvas("canvasPositronTPCdEdxVsEta_mc","",10,10,500,500);  // gives the page size
	canvasPositronTPCdEdxVsEta_mc->SetLeftMargin(0.16);
	canvasPositronTPCdEdxVsEta_mc->cd();
	hESDDalitzPositronAfterTPCdEdxVsEta_mc->SetStats(kFALSE);
	canvasPositronTPCdEdxVsEta_mc->SetLogz(1);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterTPCdEdxVsEta_mc,
						 "", "#eta", "dE/dx","MC",
						 kTRUE,-10., 10.,
						 kTRUE, -2.,2);
	
	
	canvasPositronTPCdEdxVsEta_mc->Update();
	canvasPositronTPCdEdxVsEta_mc->SaveAs(Form("%s/Positron_After_MassCut_TPCdEdxVsEta_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronTPCdEdxVsEta_mc;
	}
	
	if( hESDDalitzElectronAfterTPCdEdxVsPhi_data ) {
	
	TCanvas * canvasElectronAfterTPCdEdxVsPhi_data = new TCanvas("canvasElectronAfterTPCdEdxVsPhi_data","",10,10,500,500);  // gives the page size
	canvasElectronAfterTPCdEdxVsPhi_data->SetLeftMargin(0.16);
	canvasElectronAfterTPCdEdxVsPhi_data->cd();
	canvasElectronAfterTPCdEdxVsPhi_data->SetLogz(1);
	hESDDalitzElectronAfterTPCdEdxVsPhi_data->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzElectronAfterTPCdEdxVsPhi_data,
						 "", "#phi", "dE/dx","Data",
						 kTRUE,-10., 10.,
						 kTRUE, -2.,2*TMath::Pi());
	
	
	canvasElectronAfterTPCdEdxVsPhi_data->Update();
	canvasElectronAfterTPCdEdxVsPhi_data->SaveAs(Form("%s/Electron_After_MassCut_TPCdEdxVsPhi_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronAfterTPCdEdxVsPhi_data;
	}
	
	
	if( hESDDalitzElectronAfterTPCdEdxVsPhi_mc ) {
	
	TCanvas * canvasElectronAfterTPCdEdxVsPhi_mc = new TCanvas("canvasElectronAfterTPCdEdxVsPhi_mc","",10,10,500,500);  // gives the page size
	canvasElectronAfterTPCdEdxVsPhi_mc->SetLeftMargin(0.16);
	canvasElectronAfterTPCdEdxVsPhi_mc->cd();
	canvasElectronAfterTPCdEdxVsPhi_mc->SetLogz(1);
	
	hESDDalitzElectronAfterTPCdEdxVsPhi_mc->SetStats(kFALSE);
       
	
	DrawAutoGammaHisto2D(	hESDDalitzElectronAfterTPCdEdxVsPhi_mc,
						 "", "#phi", "dE/dx","MC",
						 kTRUE,-10., 10.,
						 kTRUE, -2.,2*TMath::Pi());
	
	
	canvasElectronAfterTPCdEdxVsPhi_mc->Update();
	canvasElectronAfterTPCdEdxVsPhi_mc->SaveAs(Form("%s/Electron_TPCdEdxVsPhi_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasElectronAfterTPCdEdxVsPhi_mc;
	
	}
	
	
	if( hESDDalitzPositronAfterTPCdEdxVsPhi_data ) {
	
	TCanvas * canvasPositronAfterTPCdEdxVsPhi_data = new TCanvas("canvasPositronAfterTPCdEdxVsPhi_data","",10,10,500,500);  // gives the page size
	canvasPositronAfterTPCdEdxVsPhi_data->SetLeftMargin(0.16);
	canvasPositronAfterTPCdEdxVsPhi_data->cd();
	canvasPositronAfterTPCdEdxVsPhi_data->SetLogz(1);
	hESDDalitzPositronAfterTPCdEdxVsPhi_data->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterTPCdEdxVsPhi_data,
						 "", "#phi", "dE/dx","Data",
						 kTRUE,-10., 10.,
						 kTRUE, -2.,2*TMath::Pi());
	
	
	canvasPositronAfterTPCdEdxVsPhi_data->Update();
	canvasPositronAfterTPCdEdxVsPhi_data->SaveAs(Form("%s/Positron_After_MassCut_TPCdEdxVsPhi_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronAfterTPCdEdxVsPhi_data;
	
	}
	
	
	if( hESDDalitzPositronAfterTPCdEdxVsPhi_mc ) {
	
	
	TCanvas * canvasPositronAfterTPCdEdxVsPhi_mc = new TCanvas("canvasPositronAfterTPCdEdxVsPhi_mc","",10,10,500,500);  // gives the page size
	canvasPositronAfterTPCdEdxVsPhi_mc->SetLeftMargin(0.16);
	canvasPositronAfterTPCdEdxVsPhi_mc->cd();
	canvasPositronAfterTPCdEdxVsPhi_mc->SetLogz(1);
	hESDDalitzPositronAfterTPCdEdxVsPhi_mc->SetStats(kFALSE);
       
	
	 DrawAutoGammaHisto2D(	hESDDalitzPositronAfterTPCdEdxVsPhi_mc,
						 "", "#phi", "dE/dx","MC",
						 kTRUE,-10., 10.,
						 kTRUE, -2.,2*TMath::Pi());
	 
	
	
	canvasPositronAfterTPCdEdxVsPhi_mc->Update();
	canvasPositronAfterTPCdEdxVsPhi_mc->SaveAs(Form("%s/Positron_After_MassCut_TPCdEdxVsPhi_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasPositronAfterTPCdEdxVsPhi_mc;
	
	}
	

	
	if( hTPCdEdxSigafter_data ) {
	
	TCanvas * canvasGammaTPCdEdxSigVsP_data = new TCanvas("canvasGammaTPCdEdxSigVsP_data","",10,10,500,500);  // gives the page size
	canvasGammaTPCdEdxSigVsP_data->SetLeftMargin(0.13);
	canvasGammaTPCdEdxSigVsP_data->SetRightMargin(0.12);
	canvasGammaTPCdEdxSigVsP_data->SetTopMargin(0.02);
	canvasGammaTPCdEdxSigVsP_data->SetBottomMargin(0.1);
	
	canvasGammaTPCdEdxSigVsP_data->cd();
	canvasGammaTPCdEdxSigVsP_data->SetLogz(1);
	canvasGammaTPCdEdxSigVsP_data->SetLogx(1);
	hTPCdEdxSigafter_data->SetStats(kFALSE);
	hTPCdEdxSigafter_data->SetTitle("");
	
	DrawAutoGammaHisto2D(	hTPCdEdxSigafter_data,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e} (#sigma)","",
						 kTRUE,-6.1, 6.5,
						 kTRUE, 0.,200,1.0);
	
	TLatex *process03 = new TLatex(0.50,0.905,collisionSystem.Data());
        process03->SetNDC();
        process03->SetTextColor(1);
        process03->SetTextSize(0.040);
	process03->SetTextFont(42);
        process03->Draw();
	
	TLatex *DateTex03 = new TLatex(0.50,0.855,"Data");
	DateTex03->SetNDC();
	DateTex03->SetTextColor(1);
	DateTex03->SetTextSize(0.04);
	DateTex03->SetTextFont(42);
	DateTex03->Draw();
	
	
	
	
	
	
	canvasGammaTPCdEdxSigVsP_data->Update();
	canvasGammaTPCdEdxSigVsP_data->SaveAs(Form("%s/Gamma_TPCdEdxSigVsP_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasGammaTPCdEdxSigVsP_data;	
	
	}
	
	
	if( hTPCdEdxafter_data ) {
	
	TCanvas * canvasGammaTPCdEdxVsP_data = new TCanvas("canvasGammaTPCdEdxVsP_data","",10,10,500,500);  // gives the page size
	canvasGammaTPCdEdxVsP_data->SetLeftMargin(0.13);
	canvasGammaTPCdEdxVsP_data->SetRightMargin(0.12);
	canvasGammaTPCdEdxVsP_data->SetTopMargin(0.02);
	canvasGammaTPCdEdxVsP_data->SetBottomMargin(0.1);
	
	canvasGammaTPCdEdxVsP_data->cd();
	canvasGammaTPCdEdxVsP_data->SetLogz(1);
	canvasGammaTPCdEdxVsP_data->SetLogx(1);
	hTPCdEdxafter_data->SetStats(kFALSE);
	hTPCdEdxafter_data->SetTitle("");
	
	
	 DrawAutoGammaHisto2D(	hTPCdEdxafter_data,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x}","",
						 kTRUE,0.,200.,
						 kTRUE, 0.,20.,1.0);
	 
	 
	TLatex *process03 = new TLatex(0.50,0.905,collisionSystem.Data());
        process03->SetNDC();
        process03->SetTextColor(1);
        process03->SetTextSize(0.040);
	process03->SetTextFont(42);
        process03->Draw();
	
	TLatex *DateTex03 = new TLatex(0.50,0.855,"Data");
	DateTex03->SetNDC();
	DateTex03->SetTextColor(1);
	DateTex03->SetTextSize(0.04);
	DateTex03->SetTextFont(42);
	DateTex03->Draw();
	
	
	
	canvasGammaTPCdEdxVsP_data->Update();
	canvasGammaTPCdEdxVsP_data->SaveAs(Form("%s/Gamma_TPCdEdxVsP_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasGammaTPCdEdxVsP_data;	
	
	}
	
	
	if( hTPCdEdxSigafter_mc ) {
	
	TCanvas * canvasGammaTPCdEdxSigVsP_mc = new TCanvas("canvasGammaTPCdEdxSigVsP_mc","",10,10,500,500);  // gives the page size
	canvasGammaTPCdEdxSigVsP_mc->SetLeftMargin(0.13);
	canvasGammaTPCdEdxSigVsP_mc->SetRightMargin(0.12);
	canvasGammaTPCdEdxSigVsP_mc->SetTopMargin(0.05);
	canvasGammaTPCdEdxSigVsP_mc->SetBottomMargin(0.1);
	
	canvasGammaTPCdEdxSigVsP_mc->cd();
	canvasGammaTPCdEdxSigVsP_mc->SetLogz(1);
	canvasGammaTPCdEdxSigVsP_mc->SetLogx(1);
	hTPCdEdxSigafter_mc->SetStats(kFALSE);
	hTPCdEdxSigafter_mc->SetTitle("");
        //hTPCdEdxSigSigafter_mc->GetXaxis()->SetTitleOffset(1.5);
	
	 DrawAutoGammaHisto2D(	hTPCdEdxSigafter_mc,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x} - #LTd#it{E}/d#it{x}#GT_{e} (#sigma)","MC",
						 kTRUE,-4.1, 5.1,
						 kTRUE, 0.,200,1.0);
	
	
	canvasGammaTPCdEdxSigVsP_mc->Update();
	canvasGammaTPCdEdxSigVsP_mc->SaveAs(Form("%s/Gamma_TPCdEdxSigVsP_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasGammaTPCdEdxSigVsP_mc;		
	
	}
	
	
	if( hTPCdEdxafter_mc ) {
	
	TCanvas * canvasGammaTPCdEdxVsP_mc = new TCanvas("canvasGammaTPCdEdxVsP_mc","",10,10,500,500);  // gives the page size
	canvasGammaTPCdEdxVsP_mc->SetLeftMargin(0.13);
	canvasGammaTPCdEdxVsP_mc->SetRightMargin(0.12);
	canvasGammaTPCdEdxVsP_mc->SetTopMargin(0.05);
	canvasGammaTPCdEdxVsP_mc->SetBottomMargin(0.1);
	
	canvasGammaTPCdEdxVsP_mc->cd();
	canvasGammaTPCdEdxVsP_mc->SetLogz(1);
	canvasGammaTPCdEdxVsP_mc->SetLogx(1);
	hTPCdEdxafter_mc->SetStats(kFALSE);
	hTPCdEdxafter_mc->SetTitle("");
	
	
	 DrawAutoGammaHisto2D(	hTPCdEdxafter_mc,
						 "", "p (GeV/c)", "TPC d#it{E}/d#it{x}","",
						 kTRUE,0.,200.,
						 kTRUE, 0.,20.,1.0);
	 
	 
	 TLatex *process03 = new TLatex(0.50,0.905,collisionSystem.Data());
        process03->SetNDC();
        process03->SetTextColor(1);
        process03->SetTextSize(0.040);
	process03->SetTextFont(42);
        process03->Draw();
	
	TLatex *DateTex03 = new TLatex(0.50,0.855,"MC");
	DateTex03->SetNDC();
	DateTex03->SetTextColor(1);
	DateTex03->SetTextSize(0.04);
	DateTex03->SetTextFont(42);
	DateTex03->Draw();
	
	
	canvasGammaTPCdEdxVsP_mc->Update();
	canvasGammaTPCdEdxVsP_mc->SaveAs(Form("%s/Gamma_TPCdEdxVsP_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	delete canvasGammaTPCdEdxVsP_mc;	
	}
	
	
	if( hESDEposEnegInvMassPt_mc && hESDEposEnegTruePi0DalitzInvMassPt ) {
	
	
	TH1F* hESDEposEnegInvMass   = (TH1F*)hESDEposEnegInvMassPt_mc->ProjectionX("ESD_EposEneg_InvMass");
	TH1F* hESDEposEnegTruePi0DalitzInvMass = (TH1F*)hESDEposEnegTruePi0DalitzInvMassPt->ProjectionX("ESD_EposEneg_TruePi0Dalitz_InvMass");
	
	
	
	
	hMCPi0EposEnegInvMass            = (TH1F*)hMCPi0EposEnegInvMassPt_mc->ProjectionX("MC_Pi0_EposEneg_InvMass");

	
	
	
	
	TCanvas * canvasElectronDalitzMass = new TCanvas("canvasElectronDalitzMass","",10,10,500,500);  // gives the page size
	
	TH1F* histoElectronDalitzMass = new TH1F("histoElectronDalitzMass","histoElectronDalitzMass",1000,0.0,0.5);
     

	
       SetStyleHistoTH1ForGraphs(histoElectronDalitzMass,"#it{M}_{e^{+}e^{-}} (GeV/#it{c^{2}})","entries per 1 MeV/#it{c^{2}}", 0.032,0.04, 0.032,0.04, 1.0,1.1);
	
	
	Double_t fElectronMass = TDatabasePDG::Instance()->GetParticle(  ::kElectron   )->Mass(); 
	//Double_t fPi0Mass      = .1349766;
	
	Double_t fPi0Mass      = 0.135;
	
	
	
	
	
	TF1 * f1= new TF1("f1", Form("(abs( 1/(1-(x/0.43)**2) )**2 )*(1.-(x/%f)**2)**3*(1+2*%f**2/x**2)*(1-4*%f**2/x**2)**0.5*1/x",fPi0Mass,fElectronMass,fElectronMass),0.05,2.05);
	TF1 * f2= new TF1("f2", Form("(abs( 1/(1-(x/0.43)**2) )**2 )*(1.-(x/%f)**2)**3*(1+2*%f**2/x**2)*(1-4*%f**2/x**2)**0.5*1/x",fPi0Mass,fElectronMass,fElectronMass),0.0,0.5);
	
	
	f1->Print();
	
	
	
	TAxis *XAxisTemp   =  hMCPi0EposEnegInvMass->GetXaxis();
	
	XAxisTemp->Print();
	
	
	


	TH1F* KrollWada  =  new TH1F("KrollWada1","Kroll-Wada",100,0.0,0.5);
	//TH1F* KrollWada =   new TH1F("KrollWada1","Kroll-Wada",hMCPi0EposEnegInvMass->GetNbinsX(),hMCPi0EposEnegInvMass->GetXaxis()->GetXmin(),hMCPi0EposEnegInvMass->GetXaxis()->GetXmax());
	
	
	TH1F* KrollWada2 =  new TH1F("KrollWada2","Kroll-Wada2",hMCPi0EposEnegInvMass->GetNbinsX(),hMCPi0EposEnegInvMass->GetXaxis()->GetXmin(),hMCPi0EposEnegInvMass->GetXaxis()->GetXmax());
	
	
	TAxis *XAxis      = KrollWada->GetXaxis(); 
	Int_t       bins  = XAxis->GetNbins();
	Double_t    from  = XAxis->GetXmin();
	Double_t    to    = XAxis->GetXmax();
	Double_t *newBins = new Double_t[bins+1];
	 
	newBins[0] = from;
	Double_t factor = TMath::Power(to/from, 1./bins);
	
	cout<<"****************"<<endl;
	cout<<"to/from: "<<to/from<<" factor "<<factor<<"  1./bins "<<1./bins<<endl;
      
	for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];
        XAxis->Set(bins,newBins);
	
	
	
	
	
	
	
	
	for(Int_t i = 1; i < bins+1;  i++){
	  
	     Double_t x = KrollWada->GetXaxis()->GetBinCenter(i);
	     
	     if ( x > 0.135 ) continue;
	    
	     Double_t eval = f1->Eval(x);
	     cout<< "Eval "<< x << " "<<eval<<" binWitdth "<<KrollWada->GetBinWidth(i)<<" bin orig"<<hMCPi0EposEnegInvMass->GetBinWidth(i)<<endl;
	      
	     KrollWada->SetBinContent(i,eval);
	  
	}
	
	
	
	//TAxis *XAxisTemp2 =  hMCPi0EposEnegInvMass->GetXaxis();
	
	//KrollWada2->SetXaxiss(XAxi);
	
	for(Int_t i = 1; i <= KrollWada2->GetXaxis()->GetNbins(); i++){
	  
	    cout<<"Entro pedro"<<endl;
	  
	     
	    
	     Double_t xBin = KrollWada2->GetXaxis()->GetBinCenter(i);
	     
	     if ( xBin > 0.135  || i <=2 ) {
	        KrollWada2->SetBinContent(i,0);
	       
	       continue;
	     }
	    
	     Double_t eval = f2->Eval(xBin);
	    
	      
	     Int_t Xbin = KrollWada2->GetXaxis()->FindBin(xBin-0.001);
	     //cout<<" Pedrito xBin "<<xBin<<" Xbin "<<Xbin<<" eval "<<eval<<endl;
	     
	      cout<< "2Eval "<< xBin << " "<<eval<<" binWitdth "<<KrollWada2->GetBinWidth(i)<<" bin orig"<<hMCPi0EposEnegInvMass->GetBinWidth(i)<<endl;
	    
	     KrollWada2->SetBinContent(i,eval);
	   
	  
	}
	
	  for(Int_t i = 0; i<KrollWada2->GetXaxis()->GetNbins(); i++)
	  cout<<"X: "<<	KrollWada2->GetXaxis()->GetBinCenter(i)<<" "<<hMCPi0EposEnegInvMass->GetXaxis()->GetBinCenter(i)<<endl; 
	
	
	
	
	TH1F* KrollWadaScale0   = (TH1F*) KrollWada2->Clone("KrollWadaScale0");
	
	
	hMCPi0EposEnegInvMass->Rebin(2);
	KrollWadaScale0->Rebin(2);
	
	
	Int_t startBinIntegral = KrollWadaScale0->GetXaxis()->FindBin(0.005);
	Int_t endBinIntegral   = KrollWadaScale0->GetXaxis()->FindBin(0.010);
	
	
	Double_t integralKrollWada = KrollWadaScale0->Integral( startBinIntegral,endBinIntegral);
	
	startBinIntegral = hMCPi0EposEnegInvMass->GetXaxis()->FindBin(0.005);
	endBinIntegral   = hMCPi0EposEnegInvMass->GetXaxis()->FindBin(0.010);
	
	Double_t integralTruePi0Dalitz = hMCPi0EposEnegInvMass->Integral( startBinIntegral,endBinIntegral);
	
	 
	 
	Double_t scaling = integralTruePi0Dalitz/integralKrollWada;
	//KrollWadaScale0->Sumw2();
	KrollWadaScale0->Scale( scaling);
	
    
        canvasElectronDalitzMass->SetTickx();
        canvasElectronDalitzMass->SetTicky();
        

        canvasElectronDalitzMass->SetLogy(1);
	canvasElectronDalitzMass->SetTopMargin(0.02);
        canvasElectronDalitzMass->SetLeftMargin(0.1);
	canvasElectronDalitzMass->SetRightMargin(0.02);
	
	
        canvasElectronDalitzMass->cd();
        
        hMCPi0EposEnegInvMass->SetStats(kFALSE);
	hMCPi0EposEnegInvMass->GetXaxis()->SetTitleOffset(1.3);
	hMCPi0EposEnegInvMass->GetYaxis()->SetTitleOffset(1.4);
	hMCPi0EposEnegInvMass->GetYaxis()->SetLabelFont(42);
	hMCPi0EposEnegInvMass->GetXaxis()->SetLabelFont(42);
	hMCPi0EposEnegInvMass->GetYaxis()->SetTitleFont(62);
	hMCPi0EposEnegInvMass->GetXaxis()->SetTitleFont(62);
	hMCPi0EposEnegInvMass->GetYaxis()->SetLabelSize(0.030);
	hMCPi0EposEnegInvMass->GetYaxis()->SetTitleSize(0.035);	
	hMCPi0EposEnegInvMass->GetYaxis()->SetDecimals();
    	hMCPi0EposEnegInvMass->GetXaxis()->SetTitleSize(0.035);	
	hMCPi0EposEnegInvMass->GetXaxis()->SetLabelSize(0.030);

	
	//hMCPi0EposEnegInvMass->Rebin(5);
	//KrollWadaScale0->Rebin(5);
	
	
	
        hMCPi0EposEnegInvMass->Sumw2();
	histoElectronDalitzMass->GetXaxis()->SetRangeUser(0.001,0.15);
	//KrollWadaScale0->GetXaxis()->SetRangeUser(0.001,0.15);
	histoElectronDalitzMass->GetYaxis()->SetRangeUser(1,5e+07);
	//KrollWadaScale0->GetYaxis()->SetRangeUser(1,5e+06);
	
	//hMCPi0EposEnegInvMass->Rebin(5);
	//KrollWadaScale0->Rebin(5);
	histoElectronDalitzMass->Draw();
	hMCPi0EposEnegInvMass->Draw("same");
	KrollWadaScale0->Draw("c hist same");
	//KrollWadaScale0->Draw();
	
        
	
	
	//KrollWadaScale0->Draw();
	KrollWadaScale0->SetLineWidth(3);
	
	
	cout<<"hMCPi0EposEnegInvMass: "<<hMCPi0EposEnegInvMass->GetNbinsX()<<endl;
	

        DrawGammaSetMarker(hMCPi0EposEnegInvMass, 21, 0.8, kRed, kRed);
	//DrawGammaSetMarker(KrollWada, 3, 0.8, kBlue, kBlue);
	
	//DrawAliceLogoPerformancePlot(0.3,0.83,0.1,0.02,0.0,"06/05/2014","pPb","", "",1.0,1.0);
	
	//TLatex *ALICEsimulationText00 = new TLatex(0.70,0.91,"ALICE Simulation"); 
        TLatex *process00 = new TLatex(0.17,0.91,Form("HIJING, %s",collisionSystem.Data()));
        //ALICEsimulationText00->SetNDC();
	//ALICEsimulationText00->SetTextFont(42);
        //ALICEsimulationText00->SetTextColor(1);
        //ALICEsimulationText00->SetTextSize(0.035);
        //ALICEsimulationText00->Draw();
        process00->SetNDC();
        process00->SetTextColor(1);
        process00->SetTextSize(0.04);
	process00->SetTextFont(42);
        process00->Draw();
	

        TLegend* leg1 = new TLegend( 0.30,0.38,0.70,0.48);
        leg1->SetTextSize(0.04);                        
        leg1->SetFillColor(0);
	leg1->SetTextFont(42);
	
	leg1->SetFillStyle(0);
	leg1->SetLineColor(0);
        leg1->AddEntry(hMCPi0EposEnegInvMass,"#pi^{0} Dalitz");
	leg1->AddEntry(KrollWadaScale0,"Kroll-Wada formula");
        leg1->Draw();

      
	
        
        canvasElectronDalitzMass->Update();
        canvasElectronDalitzMass->SaveAs(Form("%s/Electron_Dalitz_InvMass_MC%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
        delete canvasElectronDalitzMass;          
        //delete hMCPi0EposEnegInvMass;
	

	
	
	TCanvas * canvasElectronDalitzMassMC= new TCanvas("canvasElectronDalitzMassMC","",10,10,500,480);  // gives the page size
	
	
	TH1F* KrollWadaScale1  = (TH1F*) KrollWada->Clone("KrollWadaScale1");
	

	
	
	 startBinIntegral = KrollWadaScale1->GetXaxis()->FindBin(0.05);
	 endBinIntegral   = KrollWadaScale1->GetXaxis()->FindBin(0.10);
	
	
	integralKrollWada = KrollWadaScale1->Integral( startBinIntegral,endBinIntegral);
	
	startBinIntegral = hESDEposEnegTruePi0DalitzInvMass->GetXaxis()->FindBin(0.05);
	endBinIntegral   = hESDEposEnegTruePi0DalitzInvMass->GetXaxis()->FindBin(0.10);
	
	integralTruePi0Dalitz = hESDEposEnegTruePi0DalitzInvMass->Integral( startBinIntegral,endBinIntegral);
	
	 
	 
	scaling = integralTruePi0Dalitz/integralKrollWada;
	//KrollWadaScale1->Sumw2();
	KrollWadaScale1->Scale( scaling);
	
	//hESDEposEnegTruePi0DalitzInvMass->Scale(1.0/scaling);
	
    
        canvasElectronDalitzMassMC->SetTickx();
        canvasElectronDalitzMassMC->SetTicky();
        

        canvasElectronDalitzMassMC->SetLogy(1);
        canvasElectronDalitzMassMC->SetLogx(1);
	canvasElectronDalitzMassMC->SetLogz(1);
        canvasElectronDalitzMassMC->SetLeftMargin(0.12);
	canvasElectronDalitzMassMC->SetRightMargin(0.03);
	canvasElectronDalitzMassMC->SetTopMargin(0.03);
       
        
        canvasElectronDalitzMassMC->cd();
        
        hESDEposEnegInvMass->SetStats(kFALSE);
	hESDEposEnegInvMass->GetXaxis()->SetTitleOffset(1.3);
	hESDEposEnegInvMass->GetYaxis()->SetTitleOffset(1.4);
	hESDEposEnegInvMass->GetYaxis()->SetLabelFont(42);
	hESDEposEnegInvMass->GetXaxis()->SetLabelFont(42);
	hESDEposEnegInvMass->GetYaxis()->SetTitleFont(62);
	hESDEposEnegInvMass->GetXaxis()->SetTitleFont(62);
	
	
        hESDEposEnegTruePi0DalitzInvMass->SetStats(kFALSE);
        hESDEposEnegInvMass->Sumw2();
        hESDEposEnegTruePi0DalitzInvMass->Sumw2();
	
	hESDEposEnegInvMass->GetXaxis()->SetRangeUser(0.001,0.21);
	hESDEposEnegInvMass->GetYaxis()->SetRangeUser(0.1,1e+08);
	hESDEposEnegTruePi0DalitzInvMass->GetXaxis()->SetRangeUser(0.055,0.21);
	hESDEposEnegTruePi0DalitzInvMass->GetYaxis()->SetRangeUser(0.1,1e+05);

	
        //hESDEposEnegInvMass->Draw();
	
	//hESDEposEnegTruePi0DalitzInvMass->Draw("same");
	hESDEposEnegTruePi0DalitzInvMass->Draw();
        hESDEposEnegInvMass->SetTitle("");
        hESDEposEnegTruePi0DalitzInvMass->SetTitle("");
        hESDEposEnegInvMass->SetYTitle("entries per 1 MeV/c^{2}");
        hESDEposEnegInvMass->SetXTitle("M_{e^{+}e^{-}} (GeV/c^{2})");
	KrollWadaScale1->SetLineWidth(3);
	
	KrollWadaScale1->Draw("c same");
	
	
	cout<<"hESDEposEnegInvMass: "<<hESDEposEnegInvMass->GetNbinsX()<<endl;
	cout<<"hMCPi0EposEnegInvMass: "<<hESDEposEnegTruePi0DalitzInvMass->GetNbinsX()<<endl;
	

        DrawGammaSetMarker(hESDEposEnegInvMass, 21, 0.3, kBlack, kBlack);
        DrawGammaSetMarker(hESDEposEnegTruePi0DalitzInvMass, 21, 0.3, kRed, kRed);
	 //DrawGammaSetMarker(KrollWadaScale1, 21, 0.3, kBlue, kBlue);
	

        TLegend* leg2 = new TLegend( 0.22,0.18,0.54,0.28);
        leg2->SetTextSize(0.027);                        
        leg2->SetFillColor(0);
        //leg2->AddEntry(hESDEposEnegInvMass,"All candidates e^{+}e^{-} pair");
        leg2->AddEntry(hESDEposEnegTruePi0DalitzInvMass,"True Dalitz e^{+}e^{-} pair");
	leg2->AddEntry(KrollWadaScale1,"Kroll-Wada formula");
        leg2->Draw();
	
	
	TLatex *ALICEsimulationText01 = new TLatex(0.18,0.44,"ALICE Simulation"); 
        TLatex *process01 = new TLatex(0.18,0.40,"HIJING p-Pb, #sqrt{s_{NN}} = 5.02 TeV");
	TLatex *DateTex01 = new TLatex(0.18,0.36,"10/05/2014");
	
	DateTex01->SetNDC();
	DateTex01->SetTextColor(1);
	DateTex01->SetTextSize(0.035);
	DateTex01->Draw();
	
	
	
        ALICEsimulationText01->SetNDC();
        ALICEsimulationText01->SetTextColor(1);
        ALICEsimulationText01->SetTextSize(0.035);
        ALICEsimulationText01->Draw();
        process01->SetNDC();
        process01->SetTextColor(1);
        process01->SetTextSize(0.035);
        process01->Draw();
	
	

        
        canvasElectronDalitzMassMC->Update();
        canvasElectronDalitzMassMC->SaveAs(Form("%s/Electron_Dalitz_InvMass%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
        delete canvasElectronDalitzMassMC;          
        delete hESDEposEnegInvMass;
        delete hESDEposEnegTruePi0DalitzInvMass;
	
	}
	
	//}
	
	if( hESDEposEnegInvMassPi0Pt_mc ){
	  
	   
	      
	     Int_t binMin, binMax;
	  
	     binMin = hESDEposEnegInvMassPi0Pt_mc->GetYaxis()->FindBin(0.0 + 0.02);
	     binMax = hESDEposEnegInvMassPi0Pt_mc->GetYaxis()->FindBin(1.0 - 0.02);
	     
	     TH1F* hESDEposEnegInvMassAll = (TH1F*) hESDEposEnegInvMassPi0Pt_mc->ProjectionX("ESD_EposEneg_InvMassAll");
	     
	     
	     hESDEposEnegInvMassPi0Below = (TH1F*) hESDEposEnegInvMassPi0Pt_mc->ProjectionX("ESD_EposEneg_InvMassPi0_Below",binMin,binMax);
	     
	     binMin = hESDEposEnegInvMassPi0Pt_mc->GetYaxis()->FindBin(1.0  + 0.02);
	     binMax = hESDEposEnegInvMassPi0Pt_mc->GetYaxis()->FindBin(10.0 - 0.02);
	    
	     
	     
	     hESDEposEnegInvMassPi0Above = (TH1F*) hESDEposEnegInvMassPi0Pt_mc->ProjectionX("ESD_EposEneg_InvMassPi0_Above",binMin,binMax);
	      
	     
	     
	    Double_t IntegralBelow = 0;
	    Double_t IntegralAbove = 0;
	    
	    binMin = hESDEposEnegInvMassPi0Below->GetXaxis()->FindBin(0.0   + 0.0001 );
	    binMax = hESDEposEnegInvMassPi0Below->GetXaxis()->FindBin(0.015 - 0.0001 );
	    
	    IntegralBelow = hESDEposEnegInvMassPi0Below->Integral( binMin,binMax );
	    
	    
	    binMin = hESDEposEnegInvMassPi0Above->GetXaxis()->FindBin(0.0   + 0.0001 );
	    binMax = hESDEposEnegInvMassPi0Above->GetXaxis()->FindBin(0.035 - 0.0001 );
	    
	    
	    IntegralAbove = hESDEposEnegInvMassPi0Above->Integral( binMin,binMax );
	    
	    Double_t IntegralBelowAbove = IntegralBelow + IntegralAbove;
	    
	    
	    Double_t IntegralTotal = hESDEposEnegInvMassAll->Integral();
	    
	    
	    cout<<"Total IntegralBelowAbove: "<<IntegralBelowAbove<<endl;
	    cout<<"Total IntegralTotal: "<<IntegralTotal<<endl;
	    
	    cout<<"Ratio: "<< IntegralBelowAbove / IntegralTotal<<endl;
	     
	     
	     
	  
	}
	
	if( sESDConvGammaZR_mc ){
	  
	TCanvas * canvasGammaZR_mc = new TCanvas("canvasGammaZR_mc","",10,10,500,500);  // gives the page size
	//rightMargin-0.02
	canvasGammaZR_mc->SetLeftMargin(0.11);
	canvasGammaZR_mc->SetRightMargin(0.10);
	canvasGammaZR_mc->SetTopMargin(0.02);
	canvasGammaZR_mc->SetBottomMargin(0.08);
	
	TH2D* f1 = sESDConvGammaZR_mc->Projection(1,0);
	
	f1->Draw("colz");
	
	canvasGammaZR_mc->Update();
          
	canvasGammaZR_mc->SaveAs(Form("%s/Gamma_Z_R_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
       
	cout <<"Entro al QA"<<endl;
	  
	}
	
	if( sESDConvGammaZR_data ){
	  
	TCanvas * canvasGammaZR_data = new TCanvas("canvasGammaZR_data","",500,440);  // gives the page size
	//canvasGammaZR_mc->SetLeftMargin(0.11);
	//canvasGammaZR_mc->SetRightMargin(0.10);
	//canvasGammaZR_mc->SetTopMargin(0.02);
	//canvasGammaZR_mc->SetBottomMargin(0.08);
	
	canvasGammaZR_data->SetLeftMargin(0.11);
	canvasGammaZR_data->SetRightMargin(0.10);
	canvasGammaZR_data->SetTopMargin(0.02);
	canvasGammaZR_data->SetBottomMargin(0.08);
	canvasGammaZR_data->SetLogz();
	
	TH2F* plot2DCanvas01 = new TH2F("plot2DCanvas01","plot2DCanvas01",1200,-150.0,150.0,800,0.0,200.0);
        SetStyleHistoTH2ForGraphs(plot2DCanvas01,"Z (cm)","R (cm)", 0.032,0.04, 0.032,0.04, 1.0,1.2);
	
	
	hESDConvGammaZR_data = (TH2D*)sESDConvGammaZR_data->Projection(1,0);
	hESDConvGammaZR_data->Scale(NormFactorData);
	
	
	
	plot2DCanvas01->GetYaxis()->SetRangeUser(0,100);
	
	
	
	plot2DCanvas01->Draw();
	hESDConvGammaZR_data->Draw("same,colz,p2");
	
	
	Double_t xMinTL = 0.49;
	Double_t xMaxTL = 0.85;
	
	if ( energy.CompareTo("7TeV") == 0 ){
	  
	  xMinTL = 0.60;
	  xMaxTL = 0.88;
	
	  
	} else if ( energy.CompareTo("2.76TeV") == 0 )  {
	  
	  xMinTL = 0.55;
	  xMaxTL = 0.85;
	  
	}
	
	
	TLegend*  legendGammaZR = new TLegend(xMinTL,0.15,xMaxTL,0.25);
	legendGammaZR->SetTextSize(0.04);
	legendGammaZR->SetTextFont(42);
	legendGammaZR->SetBorderSize(0);
	legendGammaZR->SetLineColor(-1);
	legendGammaZR->SetFillStyle(1);
	
	legendGammaZR->AddEntry((TObject*)0,"This thesis","");
	legendGammaZR->AddEntry((TObject*)0,collisionSystem.Data(),"");
	
	legendGammaZR->Draw("same");
	
	DrawStructureZRNew();
	
	canvasGammaZR_data->Update();
          
	canvasGammaZR_data->SaveAs(Form("%s/Gamma_Z_R_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
       
	cout <<"Entro al QA"<<endl;
	  
	}
	
	
	if( sESDConvGammaXY_data ){
	  
	TCanvas * canvasGammaXY_data = new TCanvas("canvasGammaXY_data","",500,440);  // gives the page size
	
	
	canvasGammaXY_data->SetLeftMargin(0.11);
	canvasGammaXY_data->SetRightMargin(0.12);
	canvasGammaXY_data->SetTopMargin(0.02);
	canvasGammaXY_data->SetBottomMargin(0.08);
	canvasGammaXY_data->SetLogz();
	
	//TH2F* plot2DCanvas01 = new TH2F("plot2DCanvas01","plot2DCanvas01",1200,-250.0,250.0,1200,-250.0,250.0);
        //SetStyleHistoTH2ForGraphs(plot2DCanvas01,"X (cm)","Y (cm)", 0.032,0.04, 0.032,0.04, 1.0,1.2);
	
	
	hESDConvGammaXY_data = (TH2D*)sESDConvGammaXY_data->Projection(1,0);
	hESDConvGammaXY_data->Scale(NormFactorData);
	
	
	
	DrawAutoGammaHisto2D(   hESDConvGammaXY_data,
									"", "X (cm)", "Y (cm)", "",
									kTRUE, -120., 120.,
									kTRUE, -120., 120.,0.9,1.2);
		
	Double_t xMinTL = 0.47;
	Double_t xMaxTL = 0.9;
	Color_t colorXY = kBlack;
	
	if ( energy.CompareTo("7TeV") == 0 ){
	  
	  xMinTL = 0.57;
	  xMaxTL = 0.88;
	  
	} else if ( energy.CompareTo("2.76TeV") == 0 )  {
	  
	  xMinTL = 0.45;
	  xMaxTL = 0.85;
	  //colorXY = 800;
	  
	}
	
	//Color_t colorXY = kYellow - 4;
	
	
	TLegend*  legendGammaXY = new TLegend(xMinTL,0.89,xMaxTL,0.98);
	legendGammaXY->SetTextSize(0.04);
	legendGammaXY->SetTextFont(42);
	legendGammaXY->SetTextColor(colorXY);
	//legendGammaXY->SetBorderSize(0);
	//legendGammaXY->SetLineColor(-1);
	//legendGammaXY->SetFillStyle(0);
	
	legendGammaXY->AddEntry((TObject*)0,"This thesis","");
	legendGammaXY->AddEntry((TObject*)0,collisionSystem.Data(),"");
	
	
	legendGammaXY->Draw("same");
	
	DrawStructureNew(colorXY);
	
	canvasGammaXY_data->Update();
          
	canvasGammaXY_data->SaveAs(Form("%s/Gamma_X_Y_data%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
       
	cout <<"Entro al QA"<<endl;
	  
	}
	
	
	
	
	if ( sESDConvGammaZR_data && sESDConvGammaZR_mc && hESDTrueConvGammaR ){
	  
	TCanvas * canvasGammaRadius = new TCanvas("canvasGammaRadius","",500,500);  // gives the page size
	canvasGammaRadius->SetLeftMargin(0.15);
	canvasGammaRadius->SetRightMargin(0.01);
	canvasGammaRadius->SetTopMargin(0.02);
	canvasGammaRadius->SetBottomMargin(0.08);
	canvasGammaRadius->SetLogy();
	
	
	  
	TH1D* hESDConvGammaR_data = (TH1D*)sESDConvGammaZR_data->Projection(1);
	
	
	
	hESDConvGammaR_data->SetLineColor(kBlack);
	hESDConvGammaR_data->SetLineWidth(4.0);
	
	TH1D* hESDConvGammaR_mc = (TH1D*)sESDConvGammaZR_mc->Projection(1);
	
	hESDConvGammaR_mc->SetLineColor(kRed);
	hESDConvGammaR_mc->SetLineWidth(4.0);
	
	hESDTrueConvGammaR->SetLineColor(kYellow-7);
	
	
	hESDTrueConvGammaR->SetLineColor(kYellow-7);
	hESDTrueConvGammaR->SetFillStyle(4050);
	hESDTrueConvGammaR->SetFillColor(kYellow-7);
	
	
	if ( energy.CompareTo("pPb_5.023TeV") == 0 ){
	  
	  
	  hESDConvGammaR_data->Scale(NormFactorData);
	  hESDConvGammaR_mc->Scale(NormFactorMC*1.1);
	  hESDTrueConvGammaR->Scale(NormFactorMC*1.1);
	  
	  
	} else {
	hESDConvGammaR_data->Scale(NormFactorData);
	hESDConvGammaR_mc->Scale(NormFactorMC);
	hESDTrueConvGammaR->Scale(NormFactorMC);
	
	
	
	}
	
	SetStyleHistoTH1ForGraphs(hESDTrueConvGammaR,"R (cm)","#frac{N_{#gamma}}{N_{ch}}", 0.032,0.04, 0.032,0.04, 1.0,1.3);
	
	
	//plot2DCanvas01->GetYaxis()->SetRangeUser(1e-6,7e-4);
	hESDTrueConvGammaR->GetXaxis()->SetRangeUser(0,185.0);
	
	Double_t xMinTL = 0.55;
	Double_t xMaxTL = 0.90;
	
	if ( energy.CompareTo("7TeV") == 0 ){
	  
	  xMinTL = 0.60;
	  xMaxTL = 0.90;
	  
	} else if ( energy.CompareTo("2.76TeV") == 0 )  {
	  
	  xMinTL = 0.60;
	  xMaxTL = 0.90;
	  
	}
	
	TLegend*  legendGammaRadius = new TLegend(xMinTL,0.87,xMaxTL,0.93);
	legendGammaRadius->SetTextSize(0.03);
	legendGammaRadius->SetTextFont(42);
	legendGammaRadius->SetBorderSize(0);
	legendGammaRadius->SetLineColor(-1);
	legendGammaRadius->SetFillStyle(0);
	
	legendGammaRadius->AddEntry((TObject*)0,collisionSystem.Data(),"");
	
	
	TLegend*  legendGammaRadiusSources = new TLegend(xMinTL,0.75,xMaxTL,0.85);
	legendGammaRadiusSources->SetTextSize(0.025);
	legendGammaRadiusSources->SetTextFont(42);
	legendGammaRadiusSources->SetBorderSize(0);
	legendGammaRadiusSources->SetLineColor(-1);
	legendGammaRadiusSources->SetFillStyle(0);
	
	legendGammaRadiusSources->AddEntry(hESDConvGammaR_data,"Data conversion candidates","l");
	legendGammaRadiusSources->AddEntry(hESDConvGammaR_mc,"MC conversion candidates","l");
	legendGammaRadiusSources->AddEntry(hESDTrueConvGammaR,"MC true conversions","f");
	
	
	
	hESDTrueConvGammaR->Draw();
	hESDConvGammaR_data->Draw("same");
	hESDConvGammaR_mc->Draw("same");
	legendGammaRadius->Draw();
	legendGammaRadiusSources->Draw();
	
	for( Int_t i = 0; i < 12; i++ ){
	 
	  TLatex *arrayLabels = new TLatex(0.18,0.5,arrayNamesRBins[i]);
	  
		 arrayLabels->SetNDC(kTRUE);
		 arrayLabels->SetTextSize(0.02);
		 arrayLabels->SetTextAngle(90);
		 //arrayLabels->Draw();
		
	}
	
	canvasGammaRadius->SaveAs(Form("%s/Gamma_Radius_data_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
	
	}
      
	
	 if( sESDMotherDalitzPlot_mc ){
	  
	TCanvas * canvasDalitzPlot_mc = new TCanvas("canvasDalitzPlot_mc","",10,10,500,500);  // gives the page size
	
	canvasDalitzPlot_mc->SetLeftMargin(0.11);
	canvasDalitzPlot_mc->SetRightMargin(0.10);
	canvasDalitzPlot_mc->SetTopMargin(0.02);
	canvasDalitzPlot_mc->SetBottomMargin(0.1);
	
	TH2D* f1 = sESDMotherDalitzPlot_mc->Projection(1,0);
	//f1->SetXTitle("Mass_{#gamma,e^{+}}^2");
	//f1->SetYTitle("Mass_{#gamma,e^{-}}^2");
	
	TLegend*  legendEposEnegEnergyInvMass = new TLegend(0.1,0.85,0.5,0.96);
	  legendEposEnegEnergyInvMass->SetTextSize(0.03);
	  legendEposEnegEnergyInvMass->SetTextFont(42);
	  legendEposEnegEnergyInvMass->SetBorderSize(0);
	  legendEposEnegEnergyInvMass->SetLineColor(-1);
	  legendEposEnegEnergyInvMass->SetFillStyle(0);
	 
	  legendEposEnegEnergyInvMass->AddEntry((TObject*)0,collisionSystem.Data(),"");
	  legendEposEnegEnergyInvMass->AddEntry((TObject*)0,"MC","");
	  legendEposEnegEnergyInvMass->AddEntry((TObject*)0,"DalitzPlot","");
	  
	SetStyleHistoTH2ForGraphs(f1,"Mass_{#gamma,e^{+}}^{2} GeV/#it{c}","Mass_{#gamma,e^{-}}^{2} GeV/#it{c}", 0.032,0.04, 0.032,0.04, 1.0,1.2);
	
	
	f1->Draw("colz");
	
	legendEposEnegEnergyInvMass->Draw();
	
	canvasDalitzPlot_mc->Update();
          
	canvasDalitzPlot_mc->SaveAs(Form("%s/DalitzPlot_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
       
	
	  
	}
	
	
	if( sESDTruePi0DalitzDalitzPlot_mc ){
	  
	TCanvas * canvasDalitzPlot_mc = new TCanvas("canvasDalitzPlot_mc","",10,10,500,500);  // gives the page size
	
	canvasDalitzPlot_mc->SetLeftMargin(0.11);
	canvasDalitzPlot_mc->SetRightMargin(0.10);
	canvasDalitzPlot_mc->SetTopMargin(0.02);
	canvasDalitzPlot_mc->SetBottomMargin(0.1);
	
	TH2D* f1 = sESDTruePi0DalitzDalitzPlot_mc->Projection(1,0);
	SetStyleHistoTH2ForGraphs(f1,"Mass_{#gamma,e^{+}}^{2}","Mass_{#gamma,e^{-}}^{2}", 0.032,0.04, 0.032,0.04, 1.0,1.2);
	
	
	TLegend*  legendEposEnegEnergyInvMass = new TLegend(0.1,0.85,0.5,0.96);
	  legendEposEnegEnergyInvMass->SetTextSize(0.03);
	  legendEposEnegEnergyInvMass->SetTextFont(42);
	  legendEposEnegEnergyInvMass->SetBorderSize(0);
	  legendEposEnegEnergyInvMass->SetLineColor(-1);
	  legendEposEnegEnergyInvMass->SetFillStyle(0);
	 
	  legendEposEnegEnergyInvMass->AddEntry((TObject*)0,collisionSystem.Data(),"");
	  legendEposEnegEnergyInvMass->AddEntry((TObject*)0,"MC: True #pi^{0}","");
	  legendEposEnegEnergyInvMass->AddEntry((TObject*)0,"DalitzPlot","");
	  
	
	f1->Draw("colz");
	legendEposEnegEnergyInvMass->Draw();
	
	canvasDalitzPlot_mc->Update();
          
	canvasDalitzPlot_mc->SaveAs(Form("%s/TruePi0DalitzPlot_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
       
	
	  
	}
	
	
	if( sMCPi0DalitzDalitzPlot_mc ) {
	  
	  
	 TCanvas * canvasDalitzPlot_mc = new TCanvas("canvasDalitzPlot_mc","",10,10,500,500);  // gives the page size
	
	canvasDalitzPlot_mc->SetLeftMargin(0.11);
	canvasDalitzPlot_mc->SetRightMargin(0.10);
	canvasDalitzPlot_mc->SetTopMargin(0.02);
	canvasDalitzPlot_mc->SetBottomMargin(0.08);
	
	TH2D* f1 = sMCPi0DalitzDalitzPlot_mc->Projection(1,0);
	
	f1->Draw("colz");
	
	canvasDalitzPlot_mc->Update();
          
	canvasDalitzPlot_mc->SaveAs(Form("%s/MC_Pi0DalitzPlot_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
       
	  
	  
	  
	}
	
	
	if(    hESDEposEnegTruePi0DalitzInvMassPt && hESDEposEnegTrueEtaDalitzInvMassPt && hESDEposEnegTruePhotonInvMassPt && hESDEposEnegInvMassPt 
	    && hESDTrueMotherEposEnegInvMassPt ){
	
	  cout<<"Entro a la masa invariante"<<endl;
	  
	  TCanvas * canvasEposEnegPt = new TCanvas("canvasEposEnegPt","",500,500);  // gives the page size
	  canvasEposEnegPt->SetLeftMargin(0.12);
	  canvasEposEnegPt->SetRightMargin(0.02);
	  canvasEposEnegPt->SetTopMargin(0.02);
	  canvasEposEnegPt->SetBottomMargin(0.1);
	  canvasEposEnegPt->SetLogy();
	  
	  
	   
	  
	  hESDEposEnegTruePi0DalitzPt = (TH1F*)hESDEposEnegTruePi0DalitzInvMassPt->ProjectionY("ESD_EposEneg_TruePi0Dalitz_Pt");
	  hESDEposEnegTrueEtaDalitzPt = (TH1F*)hESDEposEnegTrueEtaDalitzInvMassPt->ProjectionY("ESD_EposEneg_TrueEtaDalitz_Pt");
	  hESDEposEnegTruePhotonPt    = (TH1F*)hESDEposEnegTruePhotonInvMassPt->ProjectionY("ESD_EposEnegTruePhoton_Pt");
	  hESDEposEnegPt              = (TH1F*)hESDEposEnegInvMassPt->ProjectionY("ESD_EposEneg_Pt");
	  
	  hESDEposEnegTruePi0DalitzInvMass  = (TH1F*)hESDEposEnegTruePi0DalitzInvMassPt->ProjectionX("ESD_EposEneg_TruePi0Dalitz_Pt");
	  hESDEposEnegTrueEtaDalitzInvMass  = (TH1F*)hESDEposEnegTrueEtaDalitzInvMassPt->ProjectionX("ESD_EposEneg_TrueEtaDalitz_Pt");
	  hESDEposEnegTruePhotonInvMass     = (TH1F*)hESDEposEnegTruePhotonInvMassPt->ProjectionX("ESD_EposEnegTruePhoton_Pt");
	  hESDEposEnegInvMass               = (TH1F*)hESDEposEnegInvMassPt->ProjectionX("ESD_EposEneg_Pt");
	  hESDTrueEposEnegInvMass	    = (TH1F*)hESDEposEnegTrueInvMassPt->ProjectionX("ESD_TrueEPosEneg_InvMass");
	  hESDTrueMotherEposEnegInvMass     = (TH1F*)hESDTrueMotherEposEnegInvMassPt->ProjectionX("ESD_EposEneg_TrueMother_InvMass");
	  
	  SetStyleHistoTH1ForGraphs(hESDEposEnegPt,pTLabel.Data(),"counts per 100 MeV/#it{c}", 0.032,0.04, 0.032,0.04, 1.1,1.3);
	  hESDEposEnegPt->GetYaxis()->SetRangeUser(1e01,2e08);
	  
	  
	
	  
	  hESDEposEnegPt->SetFillStyle(4050);
	  hESDEposEnegPt->SetFillColor(18);
	  hESDEposEnegPt->SetLineColor(kBlack);
	  
	  hESDEposEnegTruePi0DalitzPt->SetLineColor(kRed+2);
	  hESDEposEnegTruePi0DalitzPt->SetFillStyle(4050);
	  hESDEposEnegTruePi0DalitzPt->SetFillColor(kRed-7);
	  
	  
	  hESDEposEnegTrueEtaDalitzPt->SetLineColor(kYellow+2);
	  hESDEposEnegTrueEtaDalitzPt->SetFillStyle(4050);
	  hESDEposEnegTrueEtaDalitzPt->SetFillColor(kYellow-7);
	  
	  
	  hESDEposEnegTruePhotonPt->SetLineColor(kGreen+2);
	  hESDEposEnegTruePhotonPt->SetFillStyle(4050);
	  hESDEposEnegTruePhotonPt->SetFillColor(kGreen-7);
	  
	  
	  hESDEposEnegPt->Draw();
	  hESDEposEnegTruePi0DalitzPt->Draw("same");
	  hESDEposEnegTrueEtaDalitzPt->Draw("same");
	  hESDEposEnegTruePhotonPt->Draw("same");
	  
	  TLegend*  legendEposEnegEnergy = new TLegend(0.6,0.9,0.9,0.96);
	  legendEposEnegEnergy->SetTextSize(0.03);
	  legendEposEnegEnergy->SetTextFont(42);
	  legendEposEnegEnergy->SetBorderSize(0);
	  legendEposEnegEnergy->SetLineColor(-1);
	  legendEposEnegEnergy->SetFillStyle(0);
	  legendEposEnegEnergy->AddEntry((TObject*)0,collisionSystem.Data(),"");
	  
	  legendEposEnegEnergy->Draw("same");
	
	  
	  
	  TLegend*  legendEposEnegSources = new TLegend(0.6,0.75,0.9,0.90);
	  legendEposEnegSources->SetTextSize(0.025);
	  legendEposEnegSources->SetTextFont(42);
	  legendEposEnegSources->SetBorderSize(0);
	  legendEposEnegSources->SetLineColor(-1);
	  legendEposEnegSources->SetFillStyle(0);
	
	  legendEposEnegSources->AddEntry(hESDEposEnegPt,"e^{+}e^{-} candidates","f");
	  legendEposEnegSources->AddEntry(hESDEposEnegTruePi0DalitzPt,"MC true #pi^{0} Dalitz","f");
	  legendEposEnegSources->AddEntry(hESDEposEnegTrueEtaDalitzPt,"MC true #eta Dalitz","f");
	  legendEposEnegSources->AddEntry(hESDEposEnegTruePhotonPt,"MC true photons","f");
	  
	  legendEposEnegSources->Draw("same");
		  
	  
	  canvasEposEnegPt->SaveAs(Form("%s/ElectronPositron_Pt_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
      
	  /////////////////////////////7
	  
	  TCanvas * canvasEposEnegInvMass = new TCanvas("canvasEposEnegInvMass","",500,500);  // gives the page size
	  canvasEposEnegInvMass->SetLeftMargin(0.12);
	  canvasEposEnegInvMass->SetRightMargin(0.02);
	  canvasEposEnegInvMass->SetTopMargin(0.02);
	  canvasEposEnegInvMass->SetBottomMargin(0.1);
	  canvasEposEnegInvMass->SetLogy();
	  canvasEposEnegInvMass->SetLogx();
	  
	  
	  SetStyleHistoTH1ForGraphs(hESDEposEnegInvMass,"#it{M}_{e^{+}e^{-}}","counts per 1 MeV/#it{c}", 0.032,0.04, 0.032,0.04, 1.1,1.3);
	  hESDEposEnegInvMass->GetYaxis()->SetRangeUser(5e01,6e06);
	  
	  
	
	  
	  hESDEposEnegInvMass->SetFillStyle(4050); //4050
	  hESDEposEnegInvMass->SetFillColor(18);
	  hESDEposEnegInvMass->SetLineColor(kBlack);
	  
	  hESDEposEnegTruePhotonInvMass->SetLineColor(kRed+2);
	  hESDEposEnegTruePhotonInvMass->SetFillStyle(3001);
	  hESDEposEnegTruePhotonInvMass->SetFillColor(kRed-7);
	  
	  hESDEposEnegTruePi0DalitzInvMass->SetLineColor(kGreen+2);
	  hESDEposEnegTruePi0DalitzInvMass->SetFillStyle(3002);
	  hESDEposEnegTruePi0DalitzInvMass->SetFillColor(kGreen-4);
	  
	  hESDEposEnegTrueEtaDalitzInvMass->SetLineColor(kYellow+2);
	  hESDEposEnegTrueEtaDalitzInvMass->SetFillStyle(3001);
	  hESDEposEnegTrueEtaDalitzInvMass->SetFillColor(kYellow-7);
	  
	  
	  hESDTrueEposEnegInvMass->SetLineColor(kYellow+2);
	  hESDTrueEposEnegInvMass->SetFillStyle(3001);
	  hESDTrueEposEnegInvMass->SetFillColor(kYellow+2);
	  
	  
	  
	  hESDEposEnegInvMass->Draw();
	  hESDTrueEposEnegInvMass->Draw("same");
	  hESDEposEnegTruePi0DalitzInvMass->Draw("same");
	  
	  
	  hESDEposEnegTrueEtaDalitzInvMass->Draw("same");
	  hESDEposEnegTruePhotonInvMass->Draw("same");
	  
	  //hESDTrueMotherEposEnegInvMass->Draw("same");
	  
	  TLegend*  legendEposEnegEnergyInvMass = new TLegend(0.6,0.90,0.9,0.96);
	  legendEposEnegEnergyInvMass->SetTextSize(0.03);
	  legendEposEnegEnergyInvMass->SetTextFont(42);
	  legendEposEnegEnergyInvMass->SetBorderSize(0);
	  legendEposEnegEnergyInvMass->SetLineColor(-1);
	  legendEposEnegEnergyInvMass->SetFillStyle(0);
	  legendEposEnegEnergyInvMass->AddEntry((TObject*)0,collisionSystem.Data(),"");
	  
	  legendEposEnegEnergyInvMass->Draw("same");
	
	  
	  
	  TLegend*  legendEposEnegSourcesInvMass = new TLegend(0.6,0.75,0.9,0.90);
	  legendEposEnegSourcesInvMass->SetTextSize(0.025);
	  legendEposEnegSourcesInvMass->SetTextFont(42);
	  legendEposEnegSourcesInvMass->SetBorderSize(0);
	  legendEposEnegSourcesInvMass->SetLineColor(-1);
	  legendEposEnegSourcesInvMass->SetFillStyle(0);
	
	  legendEposEnegSourcesInvMass->AddEntry(hESDEposEnegInvMass,"e^{+}e^{-} candidates","f");
	  legendEposEnegSourcesInvMass->AddEntry(hESDTrueEposEnegInvMass,"MC true e^{+}e^{-}","f");
	  //legendEposEnegSourcesInvMass->AddEntry(
	  legendEposEnegSourcesInvMass->AddEntry(hESDEposEnegTruePi0DalitzInvMass,"MC true #pi^{0} Dalitz","f");
	  legendEposEnegSourcesInvMass->AddEntry(hESDEposEnegTrueEtaDalitzInvMass,"MC true #eta Dalitz","f");
	  legendEposEnegSourcesInvMass->AddEntry(hESDEposEnegTruePhotonInvMass,"MC true photons","f");
	  
	  legendEposEnegSourcesInvMass->Draw("same");
	  
	  //compute contamination
	  //Integral 0 - 0.135 MeV
	  
	  Int_t startBinIntegral = hESDEposEnegInvMass->GetXaxis()->FindBin(0.0);
	  Int_t endBinIntegral   = hESDEposEnegInvMass->GetXaxis()->FindBin(0.135);
	    	
	  Double_t integralEposEnegInvMass = hESDEposEnegInvMass->Integral( startBinIntegral,endBinIntegral);
	 
	  
	  startBinIntegral = hESDEposEnegTruePi0DalitzInvMass->GetXaxis()->FindBin(0.0);
	  endBinIntegral   = hESDEposEnegTruePi0DalitzInvMass->GetXaxis()->FindBin(0.135);
	    	
	  Double_t integralTruePi0Dalitz = hESDEposEnegTruePi0DalitzInvMass->Integral( startBinIntegral,endBinIntegral);
	  
	  startBinIntegral = hESDEposEnegTruePhotonInvMass->GetXaxis()->FindBin(0.0);
	  endBinIntegral   = hESDEposEnegTruePhotonInvMass->GetXaxis()->FindBin(0.135);
	    	
	  Double_t integralTruePhoton = hESDEposEnegTruePhotonInvMass->Integral( startBinIntegral,endBinIntegral);
	
	
	  startBinIntegral = hESDEposEnegTrueEtaDalitzInvMass->GetXaxis()->FindBin(0.0);
	  endBinIntegral   = hESDEposEnegTrueEtaDalitzInvMass->GetXaxis()->FindBin(0.135);
	    	
	  Double_t integralTrueEtaDalitz = hESDEposEnegTrueEtaDalitzInvMass->Integral( startBinIntegral,endBinIntegral);
	
	  
	  startBinIntegral = hESDTrueEposEnegInvMass->GetXaxis()->FindBin(0.0);
	  endBinIntegral   = hESDTrueEposEnegInvMass->GetXaxis()->FindBin(0.135);
	    	
	  Double_t integralTrueEposEneg = hESDTrueEposEnegInvMass->Integral( startBinIntegral,endBinIntegral);
	
	  cout<<"Integral InvMass Pi0 Dalitz"<<endl;
	  cout<<"e+e- candidates: "<< integralEposEnegInvMass<<endl;
	  cout<<"Pi0 Dalitz: "<<integralTruePi0Dalitz<<endl;
	  cout<<"Eta Dalitz: "<<integralTrueEtaDalitz<<endl;
	  cout<<"Photon: "<<integralTruePhoton<<endl;
	  cout<<"True electrons: "<<integralTrueEposEneg<<endl;
	 
	  
	  
	  
	  
	  
	  canvasEposEnegInvMass->SaveAs(Form("%s/ElectronPositron_InvMass_mc%s_%s.%s",outputDir.Data(),period.Data(),fCutSelection.Data(),Suffix.Data()));
      
	  
	  
	  
	}  
	
	
	
	const char* nameOutput = Form("%s/ElectronQAHistos%s_%s.root",outputDir.Data(),period.Data(),fCutSelection.Data());
	TFile* fOutput = new TFile(nameOutput,"RECREATE");
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	cout<<nameOutput<<endl;
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	
	
	if( hESDTrueConvGammaRvsEffi    ) hESDTrueConvGammaRvsEffi->Write();
	if( hESDTrueConvGammaPtvsEffi   ) hESDTrueConvGammaPtvsEffi->Write();
	if( hMCGammaPtvsConvProb        ) hMCGammaPtvsConvProb->Write();
	if( hESDTruePositronPtvsEffi    ) hESDTruePositronPtvsEffi->Write();
	if( hESDTrueElectronPtvsEffi    ) hESDTrueElectronPtvsEffi->Write();
	if( hMCAllPositronsPtScaled     ) hMCAllPositronsPtScaled->Write();
        if( hMCAllElectronsPtScaled     ) hMCAllElectronsPtScaled->Write();
	if( hESDTruePositronPtScaled    ) hESDTruePositronPtScaled->Write();
        if( hESDTrueElectronPtScaled    ) hESDTrueElectronPtScaled->Write();
	if( hMCAllPositronsPtScaledBR   ) hMCAllPositronsPtScaledBR->Write();
        if( hMCAllElectronsPtScaledBR   ) hMCAllElectronsPtScaledBR->Write();
	if( hESDTruePositronPtScaledBR  ) hESDTruePositronPtScaledBR->Write();
        if( hESDTrueElectronPtScaledBR  ) hESDTrueElectronPtScaledBR->Write();
	if( hESDTruePi0DalitzElectronPtvsEffi  ) hESDTruePi0DalitzElectronPtvsEffi->Write();
	if( hESDTruePi0DalitzPositronPtvsEffi  ) hESDTruePi0DalitzPositronPtvsEffi->Write();
	if( hESDTruePi0DalitzConvGammaPtvsEffi ) hESDTruePi0DalitzConvGammaPtvsEffi->Write();
	if( hHistoTruePi0DalitzClusGammaPtvsEffi ) hHistoTruePi0DalitzClusGammaPtvsEffi->Write();
	if( hESDTruePi0DalitzGammaPtvsEffi ) hESDTruePi0DalitzGammaPtvsEffi->Write();
	if( hESDEposEnegInvMassPi0Pt_mc  ) hESDEposEnegInvMassPi0Pt_mc->Write();
	if( hESDEposEnegInvMassPi0Below ) hESDEposEnegInvMassPi0Below->Write();
	if( hESDEposEnegInvMassPi0Above ) hESDEposEnegInvMassPi0Above->Write();
	if( hMCPi0EposEnegInvMass ) hMCPi0EposEnegInvMass->Write();
	if ( hESDConvGammaZR_data ) hESDConvGammaZR_data->Write("ESD_ConvGamma_ZR_data");
	
       
	
	for(UInt_t iBin = 0; iBin < nBins; iBin++){
	  
	  if( hMCConvGammaEffiVsPtProb[iBin] ) hMCConvGammaEffiVsPtProb[iBin]->Write();
	  
	  
	}
	
	fOutput->Write();
	fOutput->Close();
	
	
		
    
}



void DrawAliceLogoPerformancePlot(Float_t startX, Float_t startY, Float_t widthLogo, Float_t textHeight, Float_t decrease, TString date, TString collisionSystem, TString textGenerator, TString textPeriod, Double_t xLengthCanvas, Double_t yLengthCanvas){
	
	
	Double_t widthLogoPix = xLengthCanvas*widthLogo;
	Double_t heightLogoPix = widthLogoPix/0.73447;
	Double_t totalXLogo = (startX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
	Double_t totalYLogo = (startY*yLengthCanvas - heightLogoPix)/yLengthCanvas;
	
	cout<<"starX       "<<startX<<endl;
	cout<<"totalYLogo  "<<totalYLogo<<endl;
	cout<<"totalXLogo  "<<totalXLogo<<endl;
	cout<<"startY      "<<startY<<endl;

	Float_t aliceStartY = totalYLogo - 0.00;  
	/*TLatex *alice = new TLatex((startX-2*decrease),aliceStartY,"ALICE Performance"); // Bo: this was modified
	alice->SetNDC();
	alice->SetTextColor(1);
	alice->SetTextFont(62);
	alice->SetTextSize(textHeight);
	alice->SetLineWidth(2);
	alice->Draw("same");
	TLatex *pp7 = new TLatex((startX-decrease),(aliceStartY-textHeight*1.1),collisionSystem.Data()); // Bo: this was modified
	pp7->SetNDC();
	pp7->SetTextColor(1);
	pp7->SetTextFont(62);	
	pp7->SetTextSize(textHeight);
	pp7->SetLineWidth(2);
	pp7->Draw("same");*/
	TLatex *today = new TLatex((startX-decrease),(aliceStartY-2*textHeight*1.1),date.Data()); // Bo: this was modified
	
	cout<<"(startX-decrease)      	       "<<(startX-decrease)<<endl;
	cout<<"(aliceStartY-2*textHeight*1.1)  "<<(aliceStartY-2*textHeight*1.1)<<endl;
	
	
	
	today->SetNDC();
	today->SetTextColor(1);
	today->SetTextFont(62);
	today->SetTextSize(textHeight);
	today->SetLineWidth(2);
	today->Draw("same");
	
	


	TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startX , totalYLogo, totalXLogo,startY);
	//  myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
	myPadLogo->SetBorderMode(0);
	myPadLogo->SetBorderSize(2);
	myPadLogo->SetFrameBorderMode(0);
	myPadLogo->SetLeftMargin(0.0);
	myPadLogo->SetTopMargin(0.0);
	myPadLogo->SetBottomMargin(0.0);
	myPadLogo->SetRightMargin(0.0);
	TASImage *myAliceLogo = new TASImage("2012-Jul-04-Performance_Logo.png");
	myPadLogo->Draw();  // to take out for not using a logo.
	myPadLogo->cd();
	myAliceLogo->Draw("same");
	
	
	
	
	
	
	
	
	
}

void Pal1()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;
   Double_t Red[3]    = { 1.00, 0.00, 0.00};
   Double_t Green[3]  = { 0.00, 1.00, 0.00};
   Double_t Blue[3]   = { 1.00, 0.00, 1.00};
   Double_t Length[3] = { 0.00, 0.50, 1.00 };
   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetPalette(50,colors);
}
void Pal2()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;
   Double_t Red[3]    = { 1.00, 0.50, 0.00};
   Double_t Green[3]  = { 0.50, 0.00, 1.00};
   Double_t Blue[3]   = { 1.00, 0.00, 0.50};
   Double_t Length[3] = { 0.00, 0.50, 1.00 };
   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetPalette(50,colors);
}






