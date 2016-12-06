//***********************************************************************************************
//**************************** CutStudiesOverview ***********************************************
//***********************************************************************************************
/************************************************************************************************
 ******     provided by Gamma Conversion Group, PWG4,                                         *****
 ******        Ana Marin, marin@physi.uni-heidelberg.de                                        *****
 ******        Friederike Bock, friederike.bock@cern.ch                                        *****
 ******        Hikari Murakami, Pedro Gonzalez                                        *****
 ***********************************************************************************************/

#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
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
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

TList* GetMCTList(TString path,TString CutSelection, TString topDirName);

Int_t  ExtractHistosbyPtHardBins(){
  
  
  TString dateForOutput                               = ReturnDateStringForOutput();
  
  TString outputDir = Form("JetJetStudies/InputFiles/%s",dateForOutput.Data());
  TString fCutSelectionRead = "00003113_00200009366300003800000000_0163103100900000";

  gSystem->Exec("mkdir -p "+outputDir);
  
  gROOT->SetStyle("Plain");
  TH1::AddDirectory(kFALSE);

  StyleSettingsThesis("pdf");

  Int_t colorsArray[4] = {kBlack, kBlue,kRed,kMagenta};
  Int_t nPtHardBins = 19;
  Int_t startPtHardBin = 5;
  TString pathDataSets[19];
  TString mainDir[19];
  TString pattern="LHC15g1a";
  TString nameMainDir[19];
  TH1D* MC_Pi0InAcc_Pt[19];
  TH1D* MC_Pi0_InAcc_Pt_merged;
  
  //  TString inputFileDir    = "/home/pgonzale/Dalitz/TRAIN_OUTPUTS/Legotrain-PCMEMCal-20161110-1213/LHC15g1aFineBins";
  //  TString inputMergedFile = "/home/pgonzale/Dalitz/TRAIN_OUTPUTS/Legotrain-PCMEMCal-20161110-1213/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_1.root";
  //**************************use all hard pT bins 2016/11/25,26******************************************
  //  TString inputFileDir    = "/afs/cern.ch/user/h/hmurakam/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161125-ConvTesting/LHC15g1aFineBins";//w weight
  //  TString inputMergedFile = "/afs/cern.ch/user/h/hmurakam/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161125-ConvTesting/GammaConv_MC_LHC15g1aFinerPtHardBins_139.root";//w weight
  //  TString inputFileDir    = "/afs/cern.ch/user/h/hmurakam/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161126-ConvTesting/LHC15g1aFineBins";//wo weight
  //  TString inputMergedFile = "/afs/cern.ch/user/h/hmurakam/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161126-ConvTesting/GammaConv_MC_LHC15g1aFinerPtHardBins_7.root";//wo weight
  //**************************removed hard bins :10,11,12,13 2016/11/28***********************************
  TString inputFileDir    = "/afs/cern.ch/user/h/hmurakam/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161128-ConvTesting_JJ/LHC15g1aFineBins";//wo weight
  TString inputMergedFile = "/afs/cern.ch/user/h/hmurakam/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161128-ConvTesting_JJ/GammaConv_MC_LHC15g1aFinerPtHardBins_7.root";//wo weight

 
  TList* fMCContainerMergedFile = (TList*)GetMCTList(inputMergedFile,fCutSelectionRead.Data(),"GammaConvV1_7");
  //  TList* fMCContainerMergedFile = (TList*)GetMCTList(inputMergedFile,fCutSelectionRead.Data(),"GammaConvV1_139");
  
  
  if( !fMCContainerMergedFile ) {
    cout<<"ERROR";
    return 0;
  }
  fMCContainerMergedFile->Print();
  TH1D* temp = (TH1D*) fMCContainerMergedFile->FindObject("MC_Pi0InAcc_Pt");
  MC_Pi0_InAcc_Pt_merged = (TH1D*)temp->Clone("MC_Pi0InAcc_Pt_merged");
  
  
  for(Int_t iPtHardBin=startPtHardBin; iPtHardBin < nPtHardBins; iPtHardBin++){
       pathDataSets[iPtHardBin] = Form("%s/GammaConv_MC_%s%d_7.root",inputFileDir.Data(),pattern.Data(),iPtHardBin);
       //       pathDataSets[iPtHardBin] = Form("%s/GammaConv_MC_%s%d_139.root",inputFileDir.Data(),pattern.Data(),iPtHardBin);
       cout<<iPtHardBin<<"\t"<<pathDataSets[iPtHardBin].Data()<<endl;
  }
  
  
  
  for(Int_t iPtHardBin=startPtHardBin; iPtHardBin<nPtHardBins; iPtHardBin++){
	TList *fMCContainer = NULL;
	  fMCContainer = (TList*)GetMCTList(pathDataSets[iPtHardBin].Data(),fCutSelectionRead.Data(),"GammaConvV1_7");
	//	fMCContainer = (TList*)GetMCTList(pathDataSets[iPtHardBin].Data(),fCutSelectionRead.Data(),"GammaConvV1_139");
 	if( !fMCContainer ){
	  cout<<"ERROR"<<endl;
 	  return 0;
 	}
 	          
         TH1D* temp                 = (TH1D*) fMCContainer->FindObject("MC_Pi0InAcc_Pt");
         MC_Pi0InAcc_Pt[iPtHardBin] = (TH1D*)temp->Clone(Form("MC_Pi0InAcc_Pt_%d",iPtHardBin));
	 delete fMCContainer;
       
       
			
  }
  
  
 TFile* fFile = new TFile(Form("%s/outputPtBins.root",outputDir.Data()),"Update");
 if(startPtHardBin == 5) MC_Pi0_InAcc_Pt_merged->Write();
 
 for(Int_t iPtHardBin=startPtHardBin;iPtHardBin<nPtHardBins;iPtHardBin++){
   if( MC_Pi0InAcc_Pt[iPtHardBin] )
   MC_Pi0InAcc_Pt[iPtHardBin]->Write();
 }
   
fFile->Clone();

  
}
  
TList* GetMCTList(TString path,TString CutSelection, TString topDirName){
  
  
  TFile* fFile = new TFile(path.Data(),"READ");
      
  if(fFile->IsZombie()){cout << "ERROR: File " << path.Data() << " could not be openend! Returning..." << endl; return 0;}
      
       	
  TList* topDir = (TList*)fFile->Get(topDirName.Data());
      
  //topDir->Print();
  TList* anaDir = (TList*)topDir->FindObject(Form("Cut Number %s",CutSelection.Data()) );
  if( !anaDir ){		
    cout<<Form("Cut Number %s",CutSelection.Data())<<" is not found in the file"<<endl;
    return 0;
  }
  TList* temp        = (TList*) anaDir->FindObject( Form("%s MC histograms",   CutSelection.Data()));
  if( !temp ){
    cout<<Form("%s MC histograms",CutSelection.Data())<<" is not found in the file"<<endl;
    return 0;
  }
      
     
      
  TList* fMCC = (TList*) temp->Clone("fMCC");
      
  fFile->Close();
  delete fFile;
      
    
  return fMCC;
  
}
  
  
  
  
  
