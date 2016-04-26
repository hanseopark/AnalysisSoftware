//***********************************************************************************************
//**************************** CutStudiesOverview ***********************************************
//***********************************************************************************************
/***************************************************************************************************
 ****** 	provided by Gamma Conversion Group, PWG4, 									*****
 ******		Ana Marin, marin@physi.uni-heidelberg.de								*****
 ******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 								*****
 ******		Friederike Bock, friederike.bock@cern.ch								*****
 ****************************************************************************************************/

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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

struct SysErrorConversion {
   Double_t value;
   Double_t error;
   //	TString name;
};

void CutStudiesDalitzOverview(const char* CombineCutsName = "CombineCuts.dat", const char *suffix = "gif", TString meson = "", TString isMC = "", const TString optionMult="", TString optionEnergy="", TString cutVariationName ="", const Int_t NumberOfCuts =1, Bool_t optGammaOn =1, TString fDalitz = "", Int_t mode=1){

   TString 	FileDirectory[50];
   TString 	cutNumber[50];
   TString 	cutNumberAdv[50];
   TString      cutName[50];
   char* 	Input[256];
   TString 	Something;
   TString 	name[10];
   TString 	prefix2;
   Double_t 	nColls[50];

   StyleSettingsThesis();

   SetPlotStyle();




   Bool_t pictDrawingOptions[4] = {kFALSE, kFALSE, kFALSE, kTRUE};

    if (cutVariationName.CompareTo("None")==0){
        cutVariationName = "";
    }



   TString outputDir = Form("CutStudies/%s",optionEnergy.Data());
   TString outputDir1 = "CutStudies";
   gSystem->Exec("mkdir "+outputDir1);
   gSystem->Exec("mkdir "+outputDir);


   TString textMeson;
   if ( meson.CompareTo("Pi0") == 0 || meson.CompareTo("Pi0EtaBinning" ) == 0 ){
        textMeson = "#pi^{0}";
        pictDrawingOptions[3] = kTRUE;
   } else {
        pictDrawingOptions[3] = kFALSE;
        textMeson = "#eta";
   }
   if (isMC.CompareTo("kTRUE") ==0){
        prefix2 = "MC";
        pictDrawingOptions[1] = kTRUE;
   }
   else {
      prefix2 = "data";
      pictDrawingOptions[1] = 	kFALSE;
   }

   Float_t 	pictDrawingCoordinates[9]     = 	{0.55, 0.8,  0.25, 0.04, 0.7,  0.5,  0.18, 0.035, 0.};
   Float_t 	pictDrawingCoordinatesSB[9]   = 	{0.6,  0.80, 0.35, 0.04, 0.15, 0.68, 0.18, 0.035, 0.};
   Float_t 	pictDrawingCoordinatesSign[9] = 	{0.6,  0.80, 0.35, 0.04, 0.7,  0.51, 0.18, 0.035, 0.};
   Float_t 	pictDrawingCoordinatesEff[9]  = 	{0.55, 0.8,  0.18, 0.04, 0.3,  0.2,  0.18, 0.035, 0.};
   Float_t 	pictDrawingCoordinatesAcc[9]  = 	{0.63, 0.2,  0.40, 0.04, 0.7,  0.33, 0.18, 0.035, 0.};
   Float_t 	pictDrawingCoordinatesInv[9]  = 	{0.63, 0.8,  0.40, 0.04, 0.7,  0.5,  0.18, 0.035, 0.};
   Float_t 	pictDrawingCoordinatesMass[9] = 	{0.4,  0.8,  0.75, 0.04, 0.2,  0.7,  0.18, 0.035, 0.};
   Float_t 	pictDrawingCoordinatesFWHM[9] = 	{0.4,  0.8,  0.75, 0.04, 0.2,  0.68, 0.18, 0.035, 0.};
   Float_t	pictDrawingCoordinatesNLO[9]  = 	{0.55, 0.8,  0.25, 0.04, 0.7,  0.5,  0.18, 0.035, 0.};


  TString collisionSystem = ReturnFullCollisionsSystem(optionEnergy);


   if (collisionSystem.CompareTo("") == 0){

      cout << "No correct collision system specification, has been given" << endl;
      return;

   }



   Int_t color[20] = {1,860,894,620,880,591,432,422,435,420,407,416,830,404,403,802,634,608,920,923};

   ifstream in(CombineCutsName);
   cout<<"Available Cuts:"<<endl;
   string TempCutNumber;
   Int_t Number = 0;

   while(getline(in, TempCutNumber)){

        TString tempCutNumberAdv = TempCutNumber;
        cutNumberAdv[Number] = tempCutNumberAdv;

        TObjArray *arr;
        arr = tempCutNumberAdv.Tokenize("_");

	    TObjString* objstrEvent    = 0x0;
            TObjString* objstrGamma    = 0x0;
            TObjString* objstrElectron = 0x0;
	    TObjString* objstrCluster  = 0x0;
            TObjString* objstrMeson    = 0x0;
	    TObjString* objstrName     = 0x0;
	    
	    TObjString* temp = (TObjString*)arr->At(0);
	    
	    
	    cout<<"Length: "<<temp->GetString().Length()<<endl;
	    
	    //if(  temp->GetString().Length() == 7 && ( temp->GetString().Contains("8000011") || temp->GetString().Contains("0000011") ) ){
	      
	    if( mode == 1 ) {  
	    
	    
		objstrEvent    = (TObjString*)arr->At(0);
		objstrGamma    = (TObjString*)arr->At(1);
		objstrElectron = (TObjString*)arr->At(2);
		objstrMeson    = (TObjString*)arr->At(3);
		objstrName     = (TObjString*)arr->At(4);
		cutNumber[Number] = objstrEvent->GetString()+"_"+objstrGamma->GetString()+"_"+ objstrElectron->GetString()+"_"+objstrMeson->GetString();
		
		if( objstrName ) {
		  
		  cutName[Number] = objstrName->GetString();
		} else {
		  cutName[Number] = "none";
		}
	    
	    } else if( mode == 6 || mode == 7 ){
	      
	        objstrEvent    	= (TObjString*)arr->At(0);
		objstrGamma    	= (TObjString*)arr->At(1);
		objstrCluster  	= (TObjString*)arr->At(2);
		objstrElectron 	= (TObjString*)arr->At(3);
		objstrMeson    	= (TObjString*)arr->At(4);
		objstrName      = (TObjString*)arr->At(5);
		cutNumber[Number] = objstrEvent->GetString()+"_"+objstrGamma->GetString()+"_"+objstrCluster->GetString()+"_"+objstrElectron->GetString()+"_"+objstrMeson->GetString();
		
		if( objstrName ) {
		  
		  cutName[Number] = objstrName->GetString();
		} else {
		  cutName[Number] = "none";
		}
		
	    
	      
	    } else {
	      
		objstrGamma    = (TObjString*)arr->At(0);
		objstrElectron = (TObjString*)arr->At(1);
		objstrMeson    = (TObjString*)arr->At(2);
		objstrName     = (TObjString*)arr->At(3);
		cutNumber[Number] = objstrGamma->GetString()+"_"+ objstrElectron->GetString()+"_"+objstrMeson->GetString();
		if( objstrName ) {
		  
		  cutName[Number] = objstrName->GetString();
		} else {
		  cutName[Number] = "none";
		}
		
		
		
	    }
	    
        
        cout<< cutNumber[Number]<<endl;
        Number++;
   }


   cout<<"=========================="<<endl;


   if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0) {
      for (Int_t i=0; i< NumberOfCuts; i++){
         nColls[i] = GetNCollFromCutNumber(cutNumber[i]);
      }
   } else {
      for (Int_t i=0; i< NumberOfCuts; i++) nColls[i] = 1.;
   }



  const Int_t ConstNumberOfCuts = NumberOfCuts;
  TString FileNameCorrected[ConstNumberOfCuts];
  TString FileNameUnCorrected[ConstNumberOfCuts];

  TFile *Cutcorrfile[ConstNumberOfCuts];
  TFile *Cutuncorrfile[ConstNumberOfCuts];


  TH1D *histoCorrectedYieldCut[ConstNumberOfCuts];
  TH1D *histoTrueEffiCut[ConstNumberOfCuts];
  TH1D *histoRawYieldCut[ConstNumberOfCuts];
  TH1D *histoAcceptanceMesonCut[ConstNumberOfCuts];
  TH1D *histoSignificanceMesonCut[ConstNumberOfCuts];
  TH1D *histoSBMesonCut[ConstNumberOfCuts];
  TH1D *histoGGFracCont[ConstNumberOfCuts];
  
    
  TH1D *histoRatioCorrectedYieldCut[ConstNumberOfCuts];
  TH1D *histoRatioTrueEffiCut[ConstNumberOfCuts];
  TH1D *histoRatioRawYieldCut[ConstNumberOfCuts];
  TH1D *histoRatioAcceptanceCut[ConstNumberOfCuts];
  TH1D *histoRatioSignificanceCut[ConstNumberOfCuts];
  TH1D *histoRatioSBCut[ConstNumberOfCuts];
  TH1D *histoRatioGGFracCont[ConstNumberOfCuts];
  
  TH1D *histoNumberOfGoodESDTracksVtxCut[ConstNumberOfCuts];
  TString cutStringsName[ConstNumberOfCuts];



  Float_t nEvt;
  TH1F *histoEventQuality;
  TGraphAsymmErrors* systErrGraphNegYieldExt;
  TGraphAsymmErrors* systErrGraphPosYieldExt;
  TString centralityString;
  Double_t maxPt	= 0;

  for (Int_t i=0; i< NumberOfCuts; i++){

      FileNameCorrected[i] = Form("%s/%s/%s_%s_GammaConvV1%sCorrection_%s.root",cutNumberAdv[i].Data(),optionEnergy.Data(), meson.Data(),prefix2.Data(),fDalitz.Data(), cutNumber[i].Data());
      cout<< FileNameCorrected[i] << endl;
      Cutcorrfile[i] = new TFile(FileNameCorrected[i]);
      FileNameUnCorrected[i] = Form("%s/%s/%s_%s_GammaConvV1%sWithoutCorrection_%s.root",cutNumberAdv[i].Data(),optionEnergy.Data(),meson.Data(),prefix2.Data(),fDalitz.Data(),cutNumber[i].Data());
      cout<< FileNameUnCorrected[i] << endl;
      Cutuncorrfile[i] = new TFile(FileNameUnCorrected[i]);

      TString fGammaCutSelection;
      TString fElectronCutSelection;
      TString fMesonCutSelection;    
      TString fEventCutSelection;
      TString fClusterCutSelection;
      
         
      if ( mode == 9 ){
		
		cout<<"Entro aqui antiguo"<<endl;
	      
		ReturnSeparatedCutNumber(cutNumberAdv[i].Data(), fGammaCutSelection, fElectronCutSelection,fMesonCutSelection,kTRUE);
		fEventCutSelection = fGammaCutSelection(0,7);
		fGammaCutSelection = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
		cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
	} else {
	  
		  
	  
		ReturnSeparatedCutNumberAdvanced(cutNumberAdv[i].Data(),fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
	}
      
      if(cutVariationName.Contains("Comparison")){
	
	cutStringsName[i] = cutName[i];
	
	   if( cutName[i].Contains("LHC13b2EMCAL") ){
	     cutStringsName[i] = "Dalitz-EMCAL";
	   } else if ( cutName[i].Contains("LHC13b2") ) {
	     cutStringsName[i] = "Dalitz-PHOS";
	     
	   }
	
      }

      else if (cutVariationName.Contains("V0Reader")){
         TString fV0Reader = fGammaCutSelection(0,1);
         cutStringsName[i] = AnalyseV0ReaderCut(fV0Reader.Atoi());      
      } else if (cutVariationName.Contains("Eta")){
	 TString fEtaCut = fGammaCutSelection(1,1);
         cutStringsName[i] = AnalyseEtaCut(fEtaCut.Atoi());
      } else if (cutVariationName.Contains("RCut")){
         TString fRCut = fGammaCutSelection(2,1);
         cutStringsName[i] = AnalyseRCut(fRCut.Atoi());
      } else if (cutVariationName.Contains("SinglePt")){
         TString fSinglePtCut = fGammaCutSelection(3,1);
         cutStringsName[i] = AnalyseSinglePtCut(fSinglePtCut.Atoi());
      } else if (cutVariationName.Contains("Cluster")){
         TString fClusterCut = fGammaCutSelection(4,1);
         cutStringsName[i] = AnalyseTPCClusterCut(fClusterCut.Atoi());
         cout << i << "\t" << fClusterCut.Data() << "\t" << cutStringsName[i].Data()<< endl;
      } else if (cutVariationName.Contains("dEdxE")){	 
         TString fdEdxCut = fGammaCutSelection(5,1);
         cutStringsName[i] = AnalyseTPCdEdxCutElectronLine(fdEdxCut.Atoi());
      } else if (cutVariationName.Contains("dEdxPi")){    
         TString fdEdxCut = fGammaCutSelection(6,3);
         cutStringsName[i] = AnalyseTPCdEdxCutPionLine(fdEdxCut.Data());      
      } else if (cutVariationName.Contains("Qt")){
         TString fQtCut = fGammaCutSelection(11,1);
         cutStringsName[i] = AnalyseQtMaxCut(fQtCut.Atoi());
      } else if (cutVariationName.Contains("Chi2")){
         TString fChi2Cut = fGammaCutSelection(12,1);
         cutStringsName[i] = AnalyseChi2GammaCut(fChi2Cut.Atoi());
      } else if (cutVariationName.Contains("Rapidity")){
         TString fRapidityCut = fMesonCutSelection(4,1);
         cutStringsName[i] = AnalyseRapidityMesonCut(fRapidityCut.Atoi());      
      } else if (cutVariationName.Contains("Alpha")){
         TString fAlphaCut = fMesonCutSelection(6,1);
         cutStringsName[i] = AnalyseAlphaMesonCut(fAlphaCut.Atoi());
	 
      } else if (cutVariationName.Contains("PsiPairDelthaPhi") ){
	
	 TString fPsiPairDelthaPhi = fElectronCutSelection(11,1);
	 
	 cout<<"String: "<< fElectronCutSelection.Data() <<" Code: "<<fPsiPairDelthaPhi.Data()<<endl;
	 
	 cutStringsName[i] = AnalysePsiPairDelthaPhiCut( fPsiPairDelthaPhi.Atoi() );
	 
      } else {
         cutStringsName[i] = cutNumberAdv[i].Data();
      }
      
      if(i == 0){

         histoEventQuality = (TH1F*)Cutcorrfile[i]->Get("NEvents");
         nEvt = GetNEvents(histoEventQuality);
         pictDrawingCoordinates[8] = nEvt;
         centralityString = GetCentralityString(cutNumberAdv[i].Data());
         cout << centralityString.Data()<< endl;
         systErrGraphNegYieldExt = (TGraphAsymmErrors*)Cutcorrfile[i]->Get(Form("%s_SystErrorRelNeg_YieldExtraction_%s",meson.Data(),centralityString.Data()));
         Double_t* negErrorYield = systErrGraphNegYieldExt->GetY();
         for (Int_t j = 0; j < systErrGraphNegYieldExt->GetN(); j++){
         negErrorYield[j] = -1*negErrorYield[j];
         }
         systErrGraphPosYieldExt = (TGraphAsymmErrors*)Cutcorrfile[i]->Get(Form("%s_SystErrorRelPos_YieldExtraction_%s",meson.Data(),centralityString.Data()));
      }

      TString nameCorrectedYield;
      TString nameEfficiency;


     

      nameCorrectedYield = "CorrectedYieldTrueEff";
      nameEfficiency     = "TrueMesonEffiPt";

      histoCorrectedYieldCut[i] =(TH1D*)Cutcorrfile[i]->Get(nameCorrectedYield.Data());
      histoCorrectedYieldCut[i]->SetName(Form("CorrectedYieldTrueEff_%s",cutNumber[i].Data()));
      histoTrueEffiCut[i] =(TH1D*)Cutcorrfile[i]->Get(nameEfficiency.Data());
      histoTrueEffiCut[i]->SetName(Form("TrueMesonEffiPt_%s",cutNumber[i].Data()));
      histoRawYieldCut[i] =(TH1D*)Cutuncorrfile[i]->Get("histoYieldMesonPerEvent");
      histoRawYieldCut[i]->SetName(Form("histoYieldMesonPerEvent_%s",cutNumber[i].Data()));


      histoAcceptanceMesonCut[i] = (TH1D*)Cutcorrfile[i]->Get("fMCMesonAccepPt");
      histoAcceptanceMesonCut[i]->SetName(Form("fMCMesonAccepPt_%s",cutNumber[i].Data()));

      histoSignificanceMesonCut[i] = (TH1D*)Cutuncorrfile[i]->Get("histoSignMeson");
      histoSignificanceMesonCut[i]->SetName(Form("histoSignMeson_%s",cutNumber[i].Data()));

      histoSBMesonCut[i] = (TH1D*)Cutuncorrfile[i]->Get("histoSBMeson");
      histoSBMesonCut[i]->SetName(Form("histoSBMesonCut_%s",cutNumber[i].Data()));
      
      histoGGFracCont[i] = (TH1D*)Cutcorrfile[i]->Get("TrueGGFracForData");
      histoGGFracCont[i]->SetName(Form("TrueGGFrac_%s",cutNumber[i].Data()));


      histoRawYieldCut[i]->Scale(1./nColls[i]);
      
      
      
      if (i > 0){
         Double_t maxPt= histoCorrectedYieldCut[i]->GetBinCenter(histoCorrectedYieldCut[i]->GetNbinsX()) + 0.5* histoCorrectedYieldCut[i]->GetBinWidth(histoCorrectedYieldCut[i]->GetNbinsX());

      }
      
      
      histoRatioCorrectedYieldCut[i] = (TH1D*) histoCorrectedYieldCut[i]->Clone(Form("histoRatioCorrectedYieldCut_%s",cutNumber[i].Data()));
      histoRatioCorrectedYieldCut[i]->Divide(histoRatioCorrectedYieldCut[i],histoCorrectedYieldCut[0],1.,1.,"B");

      histoRatioTrueEffiCut[i] = (TH1D*) histoTrueEffiCut[i]->Clone(Form("histoRatioTrueEffiCut_%s",cutNumber[i].Data()));
      histoRatioTrueEffiCut[i]->Divide(histoRatioTrueEffiCut[i],histoTrueEffiCut[0],1.,1.,"E");
		
      histoRatioRawYieldCut[i] = (TH1D*) histoRawYieldCut[i]->Clone(Form("histoRatioRawYieldCut_%s",cutNumber[i].Data()));
      histoRatioRawYieldCut[i]->Divide(histoRatioRawYieldCut[i],histoRawYieldCut[0],1.,1.,"E");
      
      histoRatioAcceptanceCut[i] = (TH1D*) histoAcceptanceMesonCut[i]->Clone(Form("histoRatioAcceptanceCut_%s",cutNumber[i].Data()));
      histoRatioAcceptanceCut[i]->Divide(  histoRatioAcceptanceCut[i], histoAcceptanceMesonCut[0],1.,1.,"E");
      
      histoRatioSignificanceCut[i] = (TH1D*) histoSignificanceMesonCut[i]->Clone(Form("histoRatioSignificanceCut_%s",cutNumber[i].Data()));
      histoRatioSignificanceCut[i]->Divide(  histoRatioSignificanceCut[i], histoSignificanceMesonCut[0],1.,1.,"E");
      
      histoRatioSBCut[i] = (TH1D*) histoSBMesonCut[i]->Clone(Form("histoSBMesonCut_%s",cutNumber[i].Data()));
      histoRatioSBCut[i]->Divide(  histoRatioSBCut[i], histoSBMesonCut[0],1.,1.,"E");
      
      histoRatioGGFracCont[i] = (TH1D*) histoGGFracCont[i]->Clone(Form("histoRatioTrueGGFrac_%s",cutNumber[i].Data()));
      histoRatioGGFracCont[i]->Divide(  histoRatioGGFracCont[i], histoGGFracCont[0],1.,1.,"E");
      

      histoNumberOfGoodESDTracksVtxCut[i] = (TH1D*)Cutcorrfile[i]->Get("GoodESDTracks");
      histoNumberOfGoodESDTracksVtxCut[i]->SetName(Form("ESD_NumberOfGoodESDTracksVtx_%s",cutNumber[i].Data()));
  }
  cout<<"=========================="<<endl;



   //**************************************************************************************
   //*************************** Plotting Gamma PurityEff *********************************
   //**************************************************************************************


   //**************************************************************************************
   //********************* Plotting RAW-Yield *********************************************
   //**************************************************************************************
  
  
  
  
  
  
  

   TCanvas* canvasRawYieldMeson = new TCanvas("canvasRawYieldMeson","",1350,1500);  // gives the page size
   canvasRawYieldMeson->SetTickx();
   canvasRawYieldMeson->SetTicky();
   canvasRawYieldMeson->SetGridx(0);
   canvasRawYieldMeson->SetGridy(0);
   canvasRawYieldMeson->SetLogy(0);
   canvasRawYieldMeson->SetLeftMargin(0.13);
   canvasRawYieldMeson->SetRightMargin(0.02);
   canvasRawYieldMeson->SetTopMargin(0.02);
   canvasRawYieldMeson->SetFillColor(0);

   TPad* padRawYield = new TPad("padRawYield", "", 0., 0.25, 1., 1.,-1, -1, -2);
   padRawYield->SetFillColor(0);
   padRawYield->GetFrame()->SetFillColor(0);
   padRawYield->SetBorderMode(0);
   padRawYield->SetLeftMargin(0.12);
   padRawYield->SetBottomMargin(0.);
   padRawYield->SetRightMargin(0.02);
   padRawYield->SetTopMargin(0.04);
   padRawYield->SetLogy();
   padRawYield->Draw();

   TPad* padRawYieldRatios = new TPad("padRawYieldRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
   padRawYieldRatios->SetFillColor(0);
   padRawYieldRatios->GetFrame()->SetFillColor(0);
   padRawYieldRatios->SetBorderMode(0);
   padRawYieldRatios->SetLeftMargin(0.12);
   padRawYieldRatios->SetBottomMargin(0.2);
   padRawYieldRatios->SetRightMargin(0.02);
   padRawYieldRatios->SetTopMargin(0.);
   padRawYieldRatios->Draw();

   padRawYield->cd();
   padRawYield->SetTickx();
   padRawYield->SetTicky();
   //padRawYield->SetLogy(1);


   TLegend* legendRawMeson = new TLegend(0.15,0.02,0.3,0.2);
   legendRawMeson->SetTextSize(0.02);
   legendRawMeson->SetFillColor(0);
   legendRawMeson->SetLineColor(0);
   
   
   
   TLegend* legendRAWYieldMesonFits = new TLegend(0.60,0.6,0.70,0.8);
   legendRAWYieldMesonFits->SetTextSize(0.03);
   legendRAWYieldMesonFits->SetFillColor(0);
   legendRAWYieldMesonFits->SetLineColor(0);
   
   TF1* fitRatioRAWYieldComparison[NumberOfCuts];
   
   
  
       
   
   Double_t minPtRatioRAWYield   = 1.40;
   Double_t maxPtRatioRAWYield   = 10.0;
   Double_t fitParameterRAWYield = 1.05;
   
   
   for(Int_t i = 0; i< NumberOfCuts; i++){
      if(i == 0){
         DrawGammaSetMarker(histoRawYieldCut[i], 20, 1., color[0], color[0]);
         DrawAutoGammaMesonHistos( histoRawYieldCut[i],
                                   "", "p_{t} (GeV/c)", Form("%s RAW Yield/(N_{ev} N_{coll})",textMeson.Data()),
                                   kTRUE, 5., 10e-10, kTRUE,
                                   kFALSE, 0.0, 0.030,
                                   kFALSE, 0., 10.);
         legendRawMeson->AddEntry(histoRawYieldCut[i],Form("standard: %s",cutStringsName[i].Data()));
      }
      else {
	
	 fitRatioRAWYieldComparison[i] = new TF1(Form("fitRatioRAWYieldComparison%03d",i),"pol0");
	 fitRatioRAWYieldComparison[i]->SetParameter(0,fitParameterRAWYield);
	 fitRatioRAWYieldComparison[i]->SetRange(minPtRatioRAWYield,maxPtRatioRAWYield);
	 histoRatioRawYieldCut[i]->Fit(fitRatioRAWYieldComparison[i],"SINRME+","",minPtRatioRAWYield,maxPtRatioRAWYield);
	 legendRAWYieldMesonFits->AddEntry(fitRatioRAWYieldComparison[i], Form("P_{0}: %f #pm %f",fitRatioRAWYieldComparison[i]->GetParameter(0),fitRatioRAWYieldComparison[i]->GetParError(0)));
   	
         if(i<20){
            DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i],color[i]);
	    fitRatioRAWYieldComparison[i]->SetLineColor(color[i]);
         } else {
            DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i-20],color[i]);
	    fitRatioRAWYieldComparison[i]->SetLineColor(color[i-20]);
         }
         histoRawYieldCut[i]->DrawCopy("same,e1,p");
         legendRawMeson->AddEntry(histoRawYieldCut[i],cutStringsName[i].Data());

	 if( NumberOfCuts < 3 ) {

                fitRatioRAWYieldComparison[i]->SetLineColor(kRed);
                fitRatioRAWYieldComparison[i]->SetLineWidth(3);       

        }
	

      }
      
   }
   
   legendRawMeson->Draw();
   
    if ( NumberOfCuts > 1 )legendRAWYieldMesonFits->Draw("sames");
   
   
   if (optionEnergy.CompareTo("pPb_5.023TeV")==0){
      TLatex *labelCentrality = new TLatex(0.65,0.93,"p-Pb #sqrt{s_{NN}} = 5.023 TeV");
      SetStyleTLatex( labelCentrality, 0.038,4);
      labelCentrality->Draw();
   } 

   
   padRawYieldRatios->cd();
   padRawYieldRatios->SetTickx();
   padRawYieldRatios->SetTicky();

 
   for(Int_t i = 0; i< NumberOfCuts; i++){

      if(i==0){

         histoRatioRawYieldCut[i]->SetYTitle("#frac{modified}{standard}");

         if( optionEnergy.Contains("Pb") ) histoRatioRawYieldCut[i]->GetYaxis()->SetRangeUser(0.9,1.1);
         else histoRatioRawYieldCut[i]->GetYaxis()->SetRangeUser(0.5,1.5);
         if (cutVariationName.Contains("Cent")) histoRatioRawYieldCut[i]->GetYaxis()->SetRangeUser(-0.2,3.2);

         histoRatioRawYieldCut[i]->GetYaxis()->SetLabelSize(0.07);
         histoRatioRawYieldCut[i]->GetYaxis()->SetNdivisions(505);
         histoRatioRawYieldCut[i]->GetYaxis()->SetTitleSize(0.1);
         histoRatioRawYieldCut[i]->GetYaxis()->SetDecimals();
         histoRatioRawYieldCut[i]->GetYaxis()->SetTitleOffset(0.5);
         histoRatioRawYieldCut[i]->GetXaxis()->SetTitleSize(0.11);
         histoRatioRawYieldCut[i]->GetXaxis()->SetLabelSize(0.08);
         histoRatioRawYieldCut[i]->SetMarkerStyle(22);
         histoRatioRawYieldCut[i]->SetMarkerSize(1.);
         histoRatioRawYieldCut[i]->SetMarkerColor(color[0]);
         histoRatioRawYieldCut[i]->SetLineColor(color[0]);
         DrawGammaSetMarker(histoRatioRawYieldCut[i], 20, 1., 1, 1);
	 histoRatioRawYieldCut[i]->Reset("ICES");
         histoRatioRawYieldCut[i]->DrawCopy("HIST");


      } else{
	
	 if(i<20){
            DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i],color[i]);	    
         } else {
            DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i-20],color[i]);	    
         }

         histoRatioRawYieldCut[i]->DrawCopy("same,e1,p");
	 fitRatioRAWYieldComparison[i]->Draw("same");
      }

      DrawGammaLines(0., maxPt,1., 1.,0.1);
   }

   canvasRawYieldMeson->Update();

   canvasRawYieldMeson->SaveAs(Form("%s/%s_%s_RAWYield.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix));
   delete canvasRawYieldMeson;




   //*****************************************************************************************
   //******************* Compare Corrected Yields ********************************************
   //*****************************************************************************************
   TCanvas* canvasCorrectedYieldMeson = new TCanvas("canvasCorrectedYieldMeson","",1350,1500);  // gives the page size
   canvasCorrectedYieldMeson->SetTickx();
   canvasCorrectedYieldMeson->SetTicky();
   canvasCorrectedYieldMeson->SetGridx(0);
   canvasCorrectedYieldMeson->SetGridy(0);
   canvasCorrectedYieldMeson->SetLogy(0);
   canvasCorrectedYieldMeson->SetLeftMargin(0.13);
   canvasCorrectedYieldMeson->SetRightMargin(0.02);
   canvasCorrectedYieldMeson->SetTopMargin(0.02);
   canvasCorrectedYieldMeson->SetFillColor(0);

   TPad* padCorrectedYield = new TPad("padCorrectedYield", "", 0., 0.33, 1., 1.,-1, -1, -2);
   padCorrectedYield->SetFillColor(0);
   padCorrectedYield->GetFrame()->SetFillColor(0);
   padCorrectedYield->SetBorderMode(0);
   padCorrectedYield->SetLeftMargin(0.12);
   padCorrectedYield->SetBottomMargin(0.);
   padCorrectedYield->SetRightMargin(0.02);
   padCorrectedYield->SetTopMargin(0.04);
   padCorrectedYield->SetLogy(1);
   padCorrectedYield->Draw();

   TPad* padCorrectedYieldRatios = new TPad("padCorrectedYieldRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
   padCorrectedYieldRatios->SetFillColor(0);
   padCorrectedYieldRatios->GetFrame()->SetFillColor(0);
   padCorrectedYieldRatios->SetBorderMode(0);
   padCorrectedYieldRatios->SetLeftMargin(0.12);
   padCorrectedYieldRatios->SetBottomMargin(0.2);
   padCorrectedYieldRatios->SetRightMargin(0.02);
   padCorrectedYieldRatios->SetTopMargin(0.);
   padCorrectedYieldRatios->SetLogy(1);
   padCorrectedYieldRatios->Draw();

   padCorrectedYield->cd();
   padCorrectedYield->SetTickx();
   padCorrectedYield->SetTicky();
   padCorrectedYield->SetLogy(1);

	
	
   TLegend* legendCorrectedYieldMeson = new TLegend(0.15,0.02,0.3,0.2);
   legendCorrectedYieldMeson->SetTextSize(0.02);
   legendCorrectedYieldMeson->SetFillColor(0);
   legendCorrectedYieldMeson->SetLineColor(0);
   
   TLegend* legendCorrectedYieldMesonFits = new TLegend(0.60,0.6,0.70,0.8);
   legendCorrectedYieldMesonFits->SetTextSize(0.03);
   legendCorrectedYieldMesonFits->SetFillColor(0);
   legendCorrectedYieldMesonFits->SetLineColor(0);
   
   
   
   
   TF1* fitRatioComparison[NumberOfCuts];
   
   TPaveText* fitParamRatioComparison[NumberOfCuts];
   
  
       
   
   Double_t minPtRatioCorrectYield   = 0.60;
   Double_t maxPtRatioCorrectYield   = 10.0;
   Double_t fitParameterCorrectYield = 1.00;
   


   for(Int_t i = 0; i< NumberOfCuts; i++){
     
      if(i == 0){
         DrawAutoGammaMesonHistos( histoCorrectedYieldCut[i],
                                   "", "p_{t} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{t}dp_{t}dy} (c/GeV)",
                                   kTRUE, 5., 5e-10,kTRUE,
                                   kFALSE, 0.0, 0.030,
                                   kFALSE, 0., 10.);
         DrawGammaSetMarker(histoCorrectedYieldCut[i], 20, 1., color[0], color[0]);
         histoCorrectedYieldCut[i]->DrawCopy("e1,p");
         legendCorrectedYieldMeson->AddEntry(histoCorrectedYieldCut[i], Form("standard: %s",cutStringsName[i].Data()));
      }

      else{
	
	 fitRatioComparison[i] = new TF1(Form("fitRatioComparison%03d",i),"pol0");
	 fitRatioComparison[i]->SetParameter(0,fitParameterCorrectYield);
	 fitRatioComparison[i]->SetRange(minPtRatioCorrectYield,maxPtRatioCorrectYield);
	 histoRatioCorrectedYieldCut[i]->Fit(fitRatioComparison[i],"SINRME+","",minPtRatioCorrectYield,maxPtRatioCorrectYield);
	 legendCorrectedYieldMesonFits->AddEntry(fitRatioComparison[i], Form("P_{0}: %f #pm %f",fitRatioComparison[i]->GetParameter(0),fitRatioComparison[i]->GetParError(0)));
   
  
	
         if(i<20){
            DrawGammaSetMarker(histoCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
	    fitRatioComparison[i]->SetLineColor(color[i]);
         } else {
            DrawGammaSetMarker(histoCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i]);
	    fitRatioComparison[i]->SetLineColor(color[i-20]);
         }

	 if( NumberOfCuts < 3 ) {

                fitRatioComparison[i]->SetLineColor(kRed);
                fitRatioComparison[i]->SetLineWidth(3);       

         }



         histoCorrectedYieldCut[i]->DrawCopy("same,e1,p");
         legendCorrectedYieldMeson->AddEntry(histoCorrectedYieldCut[i], cutStringsName[i].Data());
      }
      
      
   }
   
   legendCorrectedYieldMeson->Draw();
   
   if ( NumberOfCuts > 1 )legendCorrectedYieldMesonFits->Draw("sames");
   

   if (optionEnergy.CompareTo("pPb_5.023TeV")==0){
      TLatex *labelCentrality = new TLatex(0.65,0.93,"p-Pb #sqrt{s_{NN}} = 5.023 TeV");
      SetStyleTLatex( labelCentrality, 0.038,4);
      labelCentrality->Draw();
   }
   
   padCorrectedYieldRatios->cd();
   padCorrectedYieldRatios->SetTickx();
   padCorrectedYieldRatios->SetTicky();
   padCorrectedYieldRatios->SetLogy(0);

   for(Int_t i = 0; i< NumberOfCuts; i++){

      if(i==0){
         histoRatioCorrectedYieldCut[i]->SetYTitle("#frac{modified}{standard}");
         histoRatioCorrectedYieldCut[i]->SetXTitle("p_{t}");
         if(optionMult.CompareTo("Mult")==0) {
            if( optionEnergy.Contains("Pb")) histoRatioCorrectedYieldCut[i]->GetYaxis()->SetRangeUser(0.005,2.1);
               else histoRatioCorrectedYieldCut[i]->GetYaxis()->SetRangeUser(-1.05,12.5);
         } else {
            if( optionEnergy.Contains("PbPb")) histoRatioCorrectedYieldCut[i]->GetYaxis()->SetRangeUser(0.45,1.55);
               else if ( optionEnergy.Contains("pPb")) histoRatioCorrectedYieldCut[i]->GetYaxis()->SetRangeUser(0.94,1.06);
               else histoRatioCorrectedYieldCut[i]->GetYaxis()->SetRangeUser(0.75,1.15);
         }
         if (cutVariationName.Contains("Cent")) histoRatioCorrectedYieldCut[i]->GetYaxis()->SetRangeUser(0.005,3.2);
         
         histoRatioCorrectedYieldCut[i]->GetYaxis()->SetLabelSize(0.07);
         histoRatioCorrectedYieldCut[i]->GetYaxis()->SetNdivisions(505);
         histoRatioCorrectedYieldCut[i]->GetYaxis()->SetTitleSize(0.1);
         histoRatioCorrectedYieldCut[i]->GetYaxis()->SetDecimals();
         histoRatioCorrectedYieldCut[i]->GetYaxis()->SetTitleOffset(0.5);
         histoRatioCorrectedYieldCut[i]->GetXaxis()->SetTitleSize(0.11);
         histoRatioCorrectedYieldCut[i]->GetXaxis()->SetLabelSize(0.08);
         histoRatioCorrectedYieldCut[i]->SetMarkerStyle(22);
         histoRatioCorrectedYieldCut[i]->SetMarkerSize(1.);
         histoRatioCorrectedYieldCut[i]->SetMarkerColor(color[0]);
         histoRatioCorrectedYieldCut[i]->SetLineColor(color[0]);

         DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20, 1., 1, 1);
         histoRatioCorrectedYieldCut[i]->Reset("ICES");
         histoRatioCorrectedYieldCut[i]->DrawCopy("HIST");
		
      }
      else{
         if(i<20){
            DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
// 				fitRatioCorrectedYieldCut[i]->SetLineColor(color[i]);	
// 				fitRatioConstCorrectedYieldCut[i]->SetLineColor(color[i]);	
// 				fitRatioConstCorrectedYieldCut[i]->SetLineStyle(2);	
         } else {
            DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i]);
// 				fitRatioCorrectedYieldCut[i]->SetLineColor(color[i-20]);	
// 				fitRatioConstCorrectedYieldCut[i]->SetLineColor(color[i-20]);	
// 				fitRatioConstCorrectedYieldCut[i]->SetLineStyle(2);	
         }
         histoRatioCorrectedYieldCut[i]->DrawCopy("same,e1,p");
	 fitRatioComparison[i]->Draw("same");
				
// 			fitRatioCorrectedYieldCut[i]->Draw("same");
// 			fitRatioConstCorrectedYieldCut[i]->Draw("same");
      }
   }
   DrawGammaLines(0., maxPt,1., 1.,0.1);

   

   canvasCorrectedYieldMeson->Update();
   cout<<"save"<<endl;
   canvasCorrectedYieldMeson->SaveAs(Form("%s/%s_%s_CorrectedYield.%s",outputDir.Data(), meson.Data(),prefix2.Data(),suffix));
   delete canvasCorrectedYieldMeson;
   
   
   //********************************************Compare Acceptances**********************************************************//
   //															      //		
   //*************************************************************************************************************************//
   
   
   
   TCanvas* canvasAcceptanceMeson = new TCanvas("canvasAcceptanceMeson","",1350,1500);  // gives the page size
   canvasAcceptanceMeson->SetTickx();
   canvasAcceptanceMeson->SetTicky();
   canvasAcceptanceMeson->SetGridx(0);
   canvasAcceptanceMeson->SetGridy(0);
   canvasAcceptanceMeson->SetLogy(0);
   canvasAcceptanceMeson->SetLeftMargin(0.13);
   canvasAcceptanceMeson->SetRightMargin(0.02);
   canvasAcceptanceMeson->SetTopMargin(0.02);
   canvasAcceptanceMeson->SetFillColor(0);

   TPad* padAcceptance = new TPad("padAcceptance", "", 0., 0.33, 1., 1.,-1, -1, -2);
   padAcceptance->SetFillColor(0);
   padAcceptance->GetFrame()->SetFillColor(0);
   padAcceptance->SetBorderMode(0);
   padAcceptance->SetLeftMargin(0.12);
   padAcceptance->SetBottomMargin(0.);
   padAcceptance->SetRightMargin(0.02);
   padAcceptance->SetTopMargin(0.02);
   padAcceptance->SetLogy();
   padAcceptance->Draw();

   TPad* padAcceptanceRatios = new TPad("padAcceptanceRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
   padAcceptanceRatios->SetFillColor(0);
   padAcceptanceRatios->GetFrame()->SetFillColor(0);
   padAcceptanceRatios->SetBorderMode(0);
   padAcceptanceRatios->SetLeftMargin(0.12);
   padAcceptanceRatios->SetBottomMargin(0.2);
   padAcceptanceRatios->SetRightMargin(0.02);
   padAcceptanceRatios->SetTopMargin(0.);
   padAcceptanceRatios->Draw();

   padAcceptance->cd();
   padAcceptance->SetTickx();
   padAcceptance->SetTicky();


   TLegend* legendAcceptanceMeson = new TLegend(0.15,0.02,0.3,0.2);
   legendAcceptanceMeson->SetTextSize(0.02);
   legendAcceptanceMeson->SetFillColor(0);
   legendAcceptanceMeson->SetLineColor(0);
   
   for(Int_t i = 0; i< NumberOfCuts; i++){
     
     cout<<"Number Of Cuts : "<<NumberOfCuts<<" ///////////////////////////////////////////////////////////////////////////"<<endl;
      
      if(i == 0){
         DrawGammaSetMarker(histoAcceptanceMesonCut[i], 20, 1., color[0], color[0]);
         DrawAutoGammaMesonHistos( histoAcceptanceMesonCut[i],
                                   "", "p_{t} (GeV/c)", "Acceptance",
                                   kTRUE, 5., 10e-10, kTRUE,
                                   kFALSE, 0.0, 0.030,
                                   kFALSE, 0., 10.);
	  histoAcceptanceMesonCut[i]->DrawCopy("e1,p");
          legendAcceptanceMeson->AddEntry(histoAcceptanceMesonCut[i],Form("standard: %s",cutStringsName[i].Data()));
	 
      }
      else {
         if(i<20){
            DrawGammaSetMarker(histoAcceptanceMesonCut[i], 20+i, 1.,color[i],color[i]);
         } else {
            DrawGammaSetMarker(histoAcceptanceMesonCut[i], 20+i, 1.,color[i-20],color[i]);
         }
         histoAcceptanceMesonCut[i]->DrawCopy("same,e1,p");
         legendAcceptanceMeson->AddEntry(histoAcceptanceMesonCut[i],cutStringsName[i].Data());
      }
      
   }
   legendAcceptanceMeson->Draw();
   if (optionEnergy.CompareTo("pPb_5.023TeV")==0){
      TLatex *labelCentrality = new TLatex(0.65,0.93,"p-Pb #sqrt{s_{NN}} = 5.023 TeV");
      SetStyleTLatex( labelCentrality, 0.038,4);
      labelCentrality->Draw();
   } 

   
   padAcceptanceRatios->cd();
   padAcceptanceRatios->SetTickx();
   padAcceptanceRatios->SetTicky();
 
   for(Int_t i = 0; i< NumberOfCuts; i++){
      if(i==0){
         histoRatioAcceptanceCut[i]->SetYTitle("#frac{modified}{standard}");

         if( optionEnergy.Contains("Pb") ) histoRatioAcceptanceCut[i]->GetYaxis()->SetRangeUser(0.7,1.3);
         else histoRatioAcceptanceCut[i]->GetYaxis()->SetRangeUser(0.8,1.2);
         if (cutVariationName.Contains("Cent")) histoRatioAcceptanceCut[i]->GetYaxis()->SetRangeUser(0.4,1.6);


         histoRatioAcceptanceCut[i]->GetYaxis()->SetLabelSize(0.07);
         histoRatioAcceptanceCut[i]->GetYaxis()->SetNdivisions(505);
         histoRatioAcceptanceCut[i]->GetYaxis()->SetTitleSize(0.1);
         histoRatioAcceptanceCut[i]->GetYaxis()->SetDecimals();
         histoRatioAcceptanceCut[i]->GetYaxis()->SetTitleOffset(0.5);
         histoRatioAcceptanceCut[i]->GetXaxis()->SetTitleSize(0.11);
         histoRatioAcceptanceCut[i]->GetXaxis()->SetLabelSize(0.08);
         histoRatioAcceptanceCut[i]->SetMarkerStyle(22);
         histoRatioAcceptanceCut[i]->SetMarkerSize(1.);
         histoRatioAcceptanceCut[i]->SetMarkerColor(color[0]);
         histoRatioAcceptanceCut[i]->SetLineColor(color[0]);
         DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20, 1., 1, 1);
         histoRatioAcceptanceCut[i]->DrawCopy("HIST");

      } else{

         if(i<20){
            DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20+i, 1.,color[i],color[i]);
         } else {
            DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20+i, 1.,color[i-20],color[i]);
         }
         histoRatioAcceptanceCut[i]->DrawCopy("same,e1,p");
      }

      DrawGammaLines(0., maxPt,1., 1.,0.1);
   }

   canvasAcceptanceMeson->Update();

   canvasAcceptanceMeson->SaveAs(Form("%s/%s_%s_AcceptanceMeson.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix));
   delete canvasAcceptanceMeson;
   
  
   //**********************************Significance*******************************************************
   //
   //*****************************************************************************************************
   
   TCanvas* canvasSignificanceMeson = new TCanvas("canvasSignificanceMeson","",1350,1500);  // gives the page size
   canvasSignificanceMeson->SetTickx();
   canvasSignificanceMeson->SetTicky();
   canvasSignificanceMeson->SetGridx(0);
   canvasSignificanceMeson->SetGridy(0);
   canvasSignificanceMeson->SetLogy(0);
   canvasSignificanceMeson->SetLeftMargin(0.13);
   canvasSignificanceMeson->SetRightMargin(0.02);
   canvasSignificanceMeson->SetTopMargin(0.02);
   canvasSignificanceMeson->SetFillColor(0);

   TPad* padSignificance = new TPad("padSignificance", "", 0., 0.33, 1., 1.,-1, -1, -2);
   padSignificance->SetFillColor(0);
   padSignificance->GetFrame()->SetFillColor(0);
   padSignificance->SetBorderMode(0);
   padSignificance->SetLeftMargin(0.12);
   padSignificance->SetBottomMargin(0.);
   padSignificance->SetRightMargin(0.02);
   padSignificance->SetTopMargin(0.02);
   padSignificance->SetLogy();
   padSignificance->Draw();

   TPad* padSignificanceRatios = new TPad("padSignificanceRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
   padSignificanceRatios->SetFillColor(0);
   padSignificanceRatios->GetFrame()->SetFillColor(0);
   padSignificanceRatios->SetBorderMode(0);
   padSignificanceRatios->SetLeftMargin(0.12);
   padSignificanceRatios->SetBottomMargin(0.2);
   padSignificanceRatios->SetRightMargin(0.02);
   padSignificanceRatios->SetTopMargin(0.);
   padSignificanceRatios->Draw();

   padSignificance->cd();
   padSignificance->SetTickx();
   padSignificance->SetTicky();


   TLegend* legendSignificanceMeson = new TLegend(0.15,0.02,0.3,0.2);
   legendSignificanceMeson->SetTextSize(0.02);
   legendSignificanceMeson->SetFillColor(0);
   legendSignificanceMeson->SetLineColor(0);

   for(Int_t i = 0; i< NumberOfCuts; i++){

     cout<<"Number Of Cuts : "<<NumberOfCuts<<" ///////////////////////////////////////////////////////////////////////////"<<endl;

      if(i == 0){

         DrawGammaSetMarker(histoSignificanceMesonCut[i], 20, 1., color[0], color[0]);
         DrawAutoGammaMesonHistos( histoSignificanceMesonCut[i],
                                   "", "p_{t} (GeV/c)", "#frac{Signal}{#sqrt{Background}}",
                                   kTRUE, 5., 10e-10, kTRUE,
                                   kFALSE, 0.0, 0.030,
                                   kFALSE, 0., 10.);
         histoSignificanceMesonCut[i]->DrawCopy("e1,p");
         legendSignificanceMeson->AddEntry(histoSignificanceMesonCut[i],Form("standard: %s",cutStringsName[i].Data()));
      }
      else {

         if(i<20){
            DrawGammaSetMarker(histoSignificanceMesonCut[i], 20+i, 1.,color[i],color[i]);
         } else {
            DrawGammaSetMarker(histoSignificanceMesonCut[i], 20+i, 1.,color[i-20],color[i]);
         }
         histoSignificanceMesonCut[i]->DrawCopy("same,e1,p");
         legendSignificanceMeson->AddEntry(histoSignificanceMesonCut[i],cutStringsName[i].Data());
      }
      
   }
   legendSignificanceMeson->Draw();

   if (optionEnergy.CompareTo("pPb_5.023TeV")==0){
      TLatex *labelCentrality = new TLatex(0.65,0.93,"p-Pb #sqrt{s_{NN}} = 5.023 TeV");
      SetStyleTLatex( labelCentrality, 0.038,4);
      labelCentrality->Draw();
   } 

   
   padSignificanceRatios->cd();
   padSignificanceRatios->SetTickx();
   padSignificanceRatios->SetTicky();
   // if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0) padRawYieldRatios->SetLogy(0);
   // 	else padRawYieldRatios->SetLogy(0);

   for(Int_t i = 0; i< NumberOfCuts; i++){
      if(i==0){
         histoRatioSignificanceCut[i]->SetYTitle("#frac{modified}{standard}");
         if( optionEnergy.Contains("Pb") ) histoRatioSignificanceCut[i]->GetYaxis()->SetRangeUser(0.5,1.5);
            else histoRatioSignificanceCut[i]->GetYaxis()->SetRangeUser(0.8,1.2);
         if (cutVariationName.Contains("Cent")) histoRatioSignificanceCut[i]->GetYaxis()->SetRangeUser(0.4,1.6);
         histoRatioSignificanceCut[i]->GetYaxis()->SetLabelSize(0.07);
         histoRatioSignificanceCut[i]->GetYaxis()->SetNdivisions(505);
         histoRatioSignificanceCut[i]->GetYaxis()->SetTitleSize(0.1);
         histoRatioSignificanceCut[i]->GetYaxis()->SetDecimals();
         histoRatioSignificanceCut[i]->GetYaxis()->SetTitleOffset(0.5);
         histoRatioSignificanceCut[i]->GetXaxis()->SetTitleSize(0.11);
         histoRatioSignificanceCut[i]->GetXaxis()->SetLabelSize(0.08);
         histoRatioSignificanceCut[i]->SetMarkerStyle(22);
         histoRatioSignificanceCut[i]->SetMarkerSize(1.);
         histoRatioSignificanceCut[i]->SetMarkerColor(color[0]);
         histoRatioSignificanceCut[i]->SetLineColor(color[0]);
         DrawGammaSetMarker(histoRatioSignificanceCut[i], 20, 1., 1, 1);
         histoRatioSignificanceCut[i]->DrawCopy("HIST");
      } else{
         if(i<20){
            DrawGammaSetMarker(histoRatioSignificanceCut[i], 20+i, 1.,color[i],color[i]);
         } else {
            DrawGammaSetMarker(histoRatioSignificanceCut[i], 20+i, 1.,color[i-20],color[i]);
         }
         histoRatioSignificanceCut[i]->DrawCopy("same,e1,p");
      }

      DrawGammaLines(0., maxPt,1., 1.,0.1);
   }

   canvasSignificanceMeson->Update();

   canvasSignificanceMeson->SaveAs(Form("%s/%s_%s_SignificanceMeson.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix));
   delete canvasSignificanceMeson;
   
   
   //*************************************S/B*********************************************************
   //
   //*************************************************************************************************
   TCanvas* canvasSBMeson = new TCanvas("canvasSBMeson","",1350,1500);  // gives the page size
   canvasSBMeson->SetTickx();
   canvasSBMeson->SetTicky();
   canvasSBMeson->SetGridx(0);
   canvasSBMeson->SetGridy(0);
   canvasSBMeson->SetLogy(0);
   canvasSBMeson->SetLeftMargin(0.13);
   canvasSBMeson->SetRightMargin(0.02);
   canvasSBMeson->SetTopMargin(0.02);
   canvasSBMeson->SetFillColor(0);

   TPad* padSB = new TPad("padSB", "", 0., 0.33, 1., 1.,-1, -1, -2);
   padSB->SetFillColor(0);
   padSB->GetFrame()->SetFillColor(0);
   padSB->SetBorderMode(0);
   padSB->SetLeftMargin(0.12);
   padSB->SetBottomMargin(0.);
   padSB->SetRightMargin(0.02);
   padSB->SetTopMargin(0.02);
   padSB->SetLogy();
   padSB->Draw();

   TPad* padSBRatios = new TPad("padSBRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
   padSBRatios->SetFillColor(0);
   padSBRatios->GetFrame()->SetFillColor(0);
   padSBRatios->SetBorderMode(0);
   padSBRatios->SetLeftMargin(0.12);
   padSBRatios->SetBottomMargin(0.2);
   padSBRatios->SetRightMargin(0.02);
   padSBRatios->SetTopMargin(0.);
   padSBRatios->Draw();

   padSB->cd();
   padSB->SetTickx();
   padSB->SetTicky();


   TLegend* legendSBMeson = new TLegend(0.15,0.02,0.3,0.2);
   legendSBMeson->SetTextSize(0.02);
   legendSBMeson->SetFillColor(0);
   legendSBMeson->SetLineColor(0);
   
   for(Int_t i = 0; i< NumberOfCuts; i++){
     
     cout<<"Number Of Cuts : "<<NumberOfCuts<<" ///////////////////////////////////////////////////////////////////////////"<<endl;
      
      if(i == 0){
         DrawGammaSetMarker(histoSBMesonCut[i], 20, 1., color[0], color[0]);
         DrawAutoGammaMesonHistos( histoSBMesonCut[i],
                                   "", "p_{t} (GeV/c)", "Signal/Background",
                                   kTRUE, 5., 10e-10, kTRUE,
                                   kFALSE, 0.0, 0.030,
                                   kFALSE, 0., 10.);
	  histoSBMesonCut[i]->DrawCopy("e1,p");
          legendSBMeson->AddEntry(histoSBMesonCut[i],Form("standard: %s",cutStringsName[i].Data()));
	 
      }
      else {
         if(i<20){
            DrawGammaSetMarker(histoSBMesonCut[i], 20+i, 1.,color[i],color[i]);
         } else {
            DrawGammaSetMarker(histoSBMesonCut[i], 20+i, 1.,color[i-20],color[i]);
         }
         histoSBMesonCut[i]->DrawCopy("same,e1,p");
         legendSBMeson->AddEntry(histoSBMesonCut[i],cutStringsName[i].Data());
      }
      
   }
   legendSBMeson->Draw();
   if (optionEnergy.CompareTo("pPb_5.023TeV")==0){
      TLatex *labelCentrality = new TLatex(0.65,0.93,"p-Pb #sqrt{s_{NN}} = 5.023 TeV");
      SetStyleTLatex( labelCentrality, 0.038,4);
      labelCentrality->Draw();
   } 

   
   padSBRatios->cd();
   padSBRatios->SetTickx();
   padSBRatios->SetTicky();
   // if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0) padRawYieldRatios->SetLogy(0);
   // 	else padRawYieldRatios->SetLogy(0);

   for(Int_t i = 0; i< NumberOfCuts; i++){
      if(i==0){
         histoRatioSBCut[i]->SetYTitle("#frac{modified}{standard}");
         if( optionEnergy.Contains("Pb") ) histoRatioSBCut[i]->GetYaxis()->SetRangeUser(0.5,1.5);
            else histoRatioSBCut[i]->GetYaxis()->SetRangeUser(0.8,1.2);
         if (cutVariationName.Contains("Cent")) histoRatioSBCut[i]->GetYaxis()->SetRangeUser(0.4,1.6);
         histoRatioSBCut[i]->GetYaxis()->SetLabelSize(0.07);
         histoRatioSBCut[i]->GetYaxis()->SetNdivisions(505);
         histoRatioSBCut[i]->GetYaxis()->SetTitleSize(0.1);
         histoRatioSBCut[i]->GetYaxis()->SetDecimals();
         histoRatioSBCut[i]->GetYaxis()->SetTitleOffset(0.5);
         histoRatioSBCut[i]->GetXaxis()->SetTitleSize(0.11);
         histoRatioSBCut[i]->GetXaxis()->SetLabelSize(0.08);
         histoRatioSBCut[i]->SetMarkerStyle(22);
         histoRatioSBCut[i]->SetMarkerSize(1.);
         histoRatioSBCut[i]->SetMarkerColor(color[0]);
         histoRatioSBCut[i]->SetLineColor(color[0]);
         DrawGammaSetMarker(histoRatioSBCut[i], 20, 1., 1, 1);
         histoRatioSBCut[i]->DrawCopy("HIST");
      } else{
         if(i<20){
            DrawGammaSetMarker(histoRatioSBCut[i], 20+i, 1.,color[i],color[i]);
         } else {
            DrawGammaSetMarker(histoRatioSBCut[i], 20+i, 1.,color[i-20],color[i]);
         }
         histoRatioSBCut[i]->DrawCopy("same,e1,p");
      }

      DrawGammaLines(0., maxPt,1., 1.,0.1);
   }

   canvasSBMeson->Update();

   canvasSBMeson->SaveAs(Form("%s/%s_%s_SB.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix));
   delete canvasSBMeson;
   
   
   //*************************************************************************************************
   //******************** Output of the systematic Error due to Signal extraction for Pi0 ************
   //*************************************************************************************************
   Int_t NBinsPt = histoCorrectedYieldCut[0]->GetNbinsX();
   const Int_t NBinstPtConst = NBinsPt+1;
	
   Double_t  BinsXCenter[NBinstPtConst];
   Double_t  BinsXWidth[NBinstPtConst];
   Double_t BinValue[NBinstPtConst];
   BinsXCenter[0] = 0;
   BinsXWidth[0]=0.;
   BinValue[0]=0.;
   for (Int_t i = 1; i < NBinsPt +1; i++){
      BinsXCenter[i] = histoCorrectedYieldCut[0]->GetBinCenter(i);
      BinsXWidth[i]= histoCorrectedYieldCut[0]->GetBinWidth(i)/2.;
   }

   SysErrorConversion SysErrCut[ConstNumberOfCuts][NBinstPtConst];
   SysErrorConversion SysErrCutRaw[ConstNumberOfCuts][NBinstPtConst];

   for (Int_t j = 0; j < NumberOfCuts; j++){
      for (Int_t i = 1; i < NBinsPt +1; i++){
         BinValue[i]= histoCorrectedYieldCut[0]->GetBinContent(i);
         SysErrCut[j][i].value = histoCorrectedYieldCut[j]->GetBinContent(i);
         SysErrCut[j][i].error = histoCorrectedYieldCut[j]->GetBinError(i);
	 SysErrCutRaw[j][i].value = histoRawYieldCut[j]->GetBinContent(i);
         SysErrCutRaw[j][i].error = histoRawYieldCut[j]->GetBinError(i);
      }
   }

   Double_t DifferenceCut[ConstNumberOfCuts][NBinstPtConst];
   Double_t DifferenceErrorCut[ConstNumberOfCuts][NBinstPtConst];

   Double_t LargestDiffNeg[NBinstPtConst];
   Double_t LargestDiffPos[NBinstPtConst];
   Double_t LargestDiffErrorNeg[NBinstPtConst];
   Double_t LargestDiffErrorPos[NBinstPtConst];

	Double_t LargestDiffRelNeg[NBinstPtConst];
   Double_t LargestDiffRelPos[NBinstPtConst];
   Double_t LargestDiffRelErrorNeg[NBinstPtConst];
   Double_t LargestDiffRelErrorPos[NBinstPtConst];

   Double_t RelDifferenceCut[ConstNumberOfCuts][NBinstPtConst];
   Double_t RelDifferenceErrorCut[ConstNumberOfCuts][NBinstPtConst];
	Double_t RelDifferenceRawCut[ConstNumberOfCuts][NBinstPtConst];
	
   for (Int_t j = 1; j < NumberOfCuts; j++){
      for ( Int_t i = 1; i < NBinstPtConst; i++) {
         DifferenceCut[j][i]=0.;
         DifferenceErrorCut[j][i]=0.;
         LargestDiffNeg[i]=0.;
         LargestDiffPos[i]=0.;
         LargestDiffErrorNeg[i]=0.;
         LargestDiffErrorPos[i]=0.;
         RelDifferenceCut[j][i]=0.;
			RelDifferenceRawCut[j][i]=0.;
         RelDifferenceErrorCut[j][i]=0.;
      }
   }

   for(Int_t j = 1; j < NumberOfCuts; j++){
      for (Int_t i = 1; i < NBinsPt +1; i++){
         //Calculate differences
         DifferenceCut[j][i] = SysErrCut[j][i].value - SysErrCut[0][i].value;
         DifferenceErrorCut[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCut[j][i].error,2)-TMath::Power(SysErrCut[0][i].error,2)));
         if(SysErrCut[0][i].value != 0){
            RelDifferenceCut[j][i] = DifferenceCut[j][i]/SysErrCut[0][i].value*100. ;
            RelDifferenceErrorCut[j][i] = DifferenceErrorCut[j][i]/SysErrCut[0][i].value*100. ;
         } else {
            RelDifferenceCut[j][i] = -10000.;
            RelDifferenceErrorCut[j][i] = 100. ;
         }
 	if(SysErrCutRaw[0][i].value != 0){
		RelDifferenceRawCut[j][i] = (SysErrCutRaw[j][i].value - SysErrCutRaw[0][i].value)/SysErrCutRaw[0][i].value*100. ;
	} else {
		RelDifferenceRawCut[j][i] = -10000.;
	}
			
	if(DifferenceCut[j][i] < 0){
	  if (TMath::Abs(LargestDiffNeg[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
	      LargestDiffNeg[i] = DifferenceCut[j][i];
	      LargestDiffErrorNeg[i] = DifferenceErrorCut[j][i];
	  }
	}else{
	  if (TMath::Abs(LargestDiffPos[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
	      LargestDiffPos[i] = DifferenceCut[j][i];
	      LargestDiffErrorPos[i] = DifferenceErrorCut[j][i];
	  }
	}
      }
   }

   cout << "done filling" << endl;
   const char *SysErrDatname = Form("%s/%s_%s_SystematicErrorCutStudies.dat",outputDir.Data(),meson.Data(),prefix2.Data());
   fstream SysErrDat;
   SysErrDat.open(SysErrDatname, ios::out);
   SysErrDat << "Calculation of the systematic error due to the yield cuts" << endl;

   cout << "works" << endl;
   for (Int_t l=0; l< NumberOfCuts; l++){
      if (l == 0) {
         SysErrDat << endl <<"Bin" << "\t" << cutNumber[l] << "\t" <<endl;
         for(Int_t i = 1; i < (NBinsPt +1); i++){
	    SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << endl;	
         }
      } else{
         SysErrDat << endl <<"Bin" << "\t" << cutNumber[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
         for(Int_t i = 1; i < (NBinsPt +1); i++){
	    if (RelDifferenceRawCut[l][i] > -75.){
	      SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << "\t" <<  DifferenceCut[l][i] << "\t"<< DifferenceErrorCut[l][i] << "\t"<< RelDifferenceCut[l][i] <<  "\t" << RelDifferenceErrorCut[l][i] <<"\t" << RelDifferenceRawCut[l][i]<< endl;
	    } else {
	      SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << "\t" <<  DifferenceCut[l][i] << "\t"<< DifferenceErrorCut[l][i] << "\t"<< RelDifferenceCut[l][i] <<  "\t" << RelDifferenceErrorCut[l][i] <<"\t" << RelDifferenceRawCut[l][i]  <<"\t not considered in largest dev" <<endl;
	    }
         }
      }
//       if (l > 0){
// 	SysErrDat << WriteParameterToFile(fitRatioCorrectedYieldCut [l])<< endl;
// 	SysErrDat << WriteParameterToFile(fitRatioConstCorrectedYieldCut[l])<< endl;
//       }
   }


   SysErrDat << endl;
   SysErrDat << endl;
   SysErrDat << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
   for(Int_t i = 1; i < (NBinsPt +1); i++){
      SysErrDat << BinsXCenter[i]  << "\t" << LargestDiffNeg[i] << "\t" <<LargestDiffErrorNeg[i]<< "\t" << LargestDiffPos[i] << "\t" << LargestDiffErrorPos[i]<<endl;
   }
   SysErrDat << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
   for(Int_t i = 0; i < (NBinsPt +1); i++){
      if ( SysErrCut[0][i].value != 0.){
	LargestDiffRelNeg[i] = - LargestDiffNeg[i]/SysErrCut[0][i].value*100.;
	LargestDiffRelPos[i] = LargestDiffPos[i]/SysErrCut[0][i].value*100.;
	LargestDiffRelErrorNeg[i] = - LargestDiffErrorNeg[i]/SysErrCut[0][i].value*100.;
	LargestDiffRelErrorPos[i] = LargestDiffErrorPos[i]/SysErrCut[0][i].value*100.;
	if (i > 0) SysErrDat << BinsXCenter[i] << "\t" << LargestDiffNeg[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffErrorNeg[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffPos[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffErrorPos[i]/SysErrCut[0][i].value*100.<<endl;
      } else {
	LargestDiffRelNeg[i] = 0.;
	LargestDiffRelPos[i] = 0.;
	LargestDiffRelErrorNeg[i] = 0.;
	LargestDiffRelErrorPos[i] = 0.;
      }
   }
    TGraphAsymmErrors* SystErrGraphNeg = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffRelNeg, BinsXWidth, BinsXWidth, LargestDiffRelErrorNeg, LargestDiffRelErrorNeg);
    SystErrGraphNeg->SetName(Form("%s_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()));
    TGraphAsymmErrors* SystErrGraphPos = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffRelPos, BinsXWidth, BinsXWidth, LargestDiffRelErrorPos, LargestDiffRelErrorPos);
    SystErrGraphPos->SetName(Form("%s_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()));

    SysErrDat.close();

    const char* Outputname = Form("%s/%s_%s_SystematicErrorCuts.root",outputDir.Data(),meson.Data(),prefix2.Data());
    TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");

    SystErrGraphPos->Write(Form("%s_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
    SystErrGraphNeg->Write(Form("%s_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
    systErrGraphNegYieldExt->Write(Form("%s_SystErrorRelNeg_YieldExtraction_%s",meson.Data(),centralityString.Data()),TObject::kOverwrite);
    systErrGraphPosYieldExt->Write(Form("%s_SystErrorRelPos_YieldExtraction_%s",meson.Data(),centralityString.Data()),TObject::kOverwrite);
    SystematicErrorFile->Write();
    SystematicErrorFile->Close();
	
    if (cutVariationName.Contains("V0")){
      const char* Outputname2 = Form("%s/%s_RatioRawYields.root",outputDir.Data(),meson.Data());
      TFile* RatioRawYieldsFile = new TFile(Outputname2,"UPDATE");
// 	fitRatioCorrectedYieldCut[1]->Write(Form("%s_fitRatioCorrectedYieldCut%s_%s",meson.Data(),cutVariationName.Data(),prefix2.Data()),TObject::kOverwrite);
// 	fitRatioConstCorrectedYieldCut[1]->Write(Form("%s_fitRatioConstCorrectedYieldCut%s_%s",meson.Data(),cutVariationName.Data(),prefix2.Data()),TObject::kOverwrite);
	histoCorrectedYieldCut[1]->Write(Form("%s_histoCorrectedYield_%s_%s",meson.Data(),cutVariationName.Data(),prefix2.Data()),TObject::kOverwrite);
	histoRatioCorrectedYieldCut[1]->Write(Form("%s_histoRatioCorrectedYield_%s_%s",meson.Data(),cutVariationName.Data(),prefix2.Data()),TObject::kOverwrite);
	histoRatioRawYieldCut[1]->Write(Form("%s_histoRatioRawYield_%s_%s",meson.Data(),cutVariationName.Data(),prefix2.Data()),TObject::kOverwrite);
      RatioRawYieldsFile->Write();
      RatioRawYieldsFile->Close();
    }

	
   //**************************************************************************************
   //********************* Plotting RAW-Yield *********************************************
   //**************************************************************************************
	
   TCanvas* canvasTrueEffiMeson = new TCanvas("canvasTrueEffiMeson","",1350,1500);  // gives the page size
   canvasTrueEffiMeson->SetTickx();
   canvasTrueEffiMeson->SetTicky();
   canvasTrueEffiMeson->SetGridx(0);
   canvasTrueEffiMeson->SetGridy(0);
   canvasTrueEffiMeson->SetLogy(0);
   canvasTrueEffiMeson->SetLeftMargin(0.13);
   canvasTrueEffiMeson->SetRightMargin(0.02);
   canvasTrueEffiMeson->SetTopMargin(0.02);
   canvasTrueEffiMeson->SetFillColor(0);

   TPad* padTrueEffi = new TPad("padTrueEffi", "", 0., 0.25, 1., 1.,-1, -1, -2);
   padTrueEffi->SetFillColor(0);
   padTrueEffi->GetFrame()->SetFillColor(0);
   padTrueEffi->SetBorderMode(0);
   padTrueEffi->SetLeftMargin(0.12);
   padTrueEffi->SetBottomMargin(0.);
   padTrueEffi->SetRightMargin(0.02);
   padTrueEffi->SetTopMargin(0.04);
   padTrueEffi->Draw();

   TPad* padTrueEffiRatios = new TPad("padTrueEffiRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
   padTrueEffiRatios->SetFillColor(0);
   padTrueEffiRatios->GetFrame()->SetFillColor(0);
   padTrueEffiRatios->SetBorderMode(0);
   padTrueEffiRatios->SetLeftMargin(0.12);
   padTrueEffiRatios->SetBottomMargin(0.2);
   padTrueEffiRatios->SetRightMargin(0.02);
   padTrueEffiRatios->SetTopMargin(0.);
   padTrueEffiRatios->Draw();

   padTrueEffi->cd();
   padTrueEffi->SetTickx();
   padTrueEffi->SetTicky();
   padTrueEffi->SetLogy(0);


   TF1*       fitRatioComparisonTrueEffi[NumberOfCuts];
   TPaveText* fitParamRatioComparisonTrueEffi[NumberOfCuts];
   
   Double_t minPtRatioEfficiency   = 1.40;
   Double_t maxPtRatioEfficiency   = 4.0;
   Double_t fitParameterEfficiency = 1.05;
   
   TLegend* legendEffiMeson = new TLegend(0.15,0.65,0.3,0.85);
   legendEffiMeson->SetTextSize(0.02);
   legendEffiMeson->SetFillColor(0);
   legendEffiMeson->SetLineColor(0);
   
   TLegend* legendEffiMesonFits = new TLegend(0.15,0.40,0.3,0.60);
   legendEffiMesonFits->SetTextSize(0.03);
   legendEffiMesonFits->SetFillColor(0);
   legendEffiMesonFits->SetLineColor(0);
   
   
   for(Int_t i = 0; i< NumberOfCuts; i++){

      if(i == 0){
			
         DrawGammaSetMarker(histoTrueEffiCut[i], 20, 1., color[0], color[0]);
         DrawAutoGammaMesonHistos( histoTrueEffiCut[i],
                                   "", "p_{t} (GeV/c)", Form("%s Efficiency",textMeson.Data()),
                                   kTRUE, 5., 10e-10,kFALSE,
                                   kTRUE, -0.1, 0.0030,
                                   kFALSE, 0., 10.);
         histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(0.0,0.06);
         histoTrueEffiCut[i]->DrawCopy("e1,p");
         legendEffiMeson->AddEntry(histoTrueEffiCut[i],Form("standard: %s",cutStringsName[i].Data()));

      }
      else {
	
	fitRatioComparisonTrueEffi[i] = new TF1(Form("fitRatioComparisonTrueEffi%03d",i),"pol0");
	fitRatioComparisonTrueEffi[i]->SetParameter(0,fitParameterEfficiency);
	fitRatioComparisonTrueEffi[i]->SetRange(minPtRatioEfficiency,maxPtRatioEfficiency);
	histoRatioTrueEffiCut[i]->Fit(fitRatioComparisonTrueEffi[i],"SINRME+","",minPtRatioEfficiency,maxPtRatioEfficiency);
	legendEffiMesonFits->AddEntry(fitRatioComparisonTrueEffi[i], Form("P_{0}: %f #pm %f",fitRatioComparisonTrueEffi[i]->GetParameter(0),fitRatioComparisonTrueEffi[i]->GetParError(0)));
  	
	
         if(i<20){
            DrawGammaSetMarker(histoTrueEffiCut[i], 20+i, 1.,color[i],color[i]);
	    fitRatioComparisonTrueEffi[i]->SetLineColor(color[i]);
           
         } else {
            DrawGammaSetMarker(histoTrueEffiCut[i], 20+i, 1.,color[i-20],color[i]);
	    fitRatioComparisonTrueEffi[i]->SetLineColor(color[i-20]);
         }

           if( NumberOfCuts < 3 ) {

		fitRatioComparisonTrueEffi[i]->SetLineColor(kRed);
		fitRatioComparisonTrueEffi[i]->SetLineWidth(3);		
	
	   }
	
         histoTrueEffiCut[i]->DrawCopy("same,e1,p");
         legendEffiMeson->AddEntry(histoTrueEffiCut[i],cutStringsName[i].Data());
      }
   }
   legendEffiMeson->Draw();
   legendEffiMesonFits->Draw("same");


  if (optionEnergy.CompareTo("pPb_5.023TeV") == 0){
      TLatex *labelCentrality = new TLatex(0.65,0.91,"p-Pb #sqrt{s_{NN}} = 5.023 TeV");
      SetStyleTLatex( labelCentrality, 0.038,4);
      labelCentrality->Draw();
  }

  padTrueEffiRatios->cd();
  padTrueEffiRatios->SetTickx();
  padTrueEffiRatios->SetTicky();

  if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0) padTrueEffiRatios->SetLogy(0);
  else padTrueEffiRatios->SetLogy(0);

   for(Int_t i = 0; i< NumberOfCuts; i++){

      if(i==0){

         if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0){

               histoRatioTrueEffiCut[i]->SetYTitle(Form("R_{CP}(%s)",textMeson.Data()));
               histoRatioTrueEffiCut[i]->GetYaxis()->SetRangeUser(0.4,1.1);

          } else {
               histoRatioTrueEffiCut[i]->SetYTitle("#frac{modified}{standard}");
               histoRatioTrueEffiCut[i]->GetYaxis()->SetRangeUser(0.5,1.5);
          }
               histoRatioTrueEffiCut[i]->GetYaxis()->SetLabelSize(0.07);
               histoRatioTrueEffiCut[i]->GetYaxis()->SetNdivisions(505);
               histoRatioTrueEffiCut[i]->GetYaxis()->SetTitleSize(0.1);
               histoRatioTrueEffiCut[i]->GetYaxis()->SetDecimals();
               histoRatioTrueEffiCut[i]->GetYaxis()->SetTitleOffset(0.5);
               histoRatioTrueEffiCut[i]->GetXaxis()->SetTitleSize(0.11);
               histoRatioTrueEffiCut[i]->GetXaxis()->SetLabelSize(0.08);
               histoRatioTrueEffiCut[i]->SetMarkerStyle(22);
               histoRatioTrueEffiCut[i]->SetMarkerSize(1.);
               histoRatioTrueEffiCut[i]->SetMarkerColor(color[0]);
               histoRatioTrueEffiCut[i]->SetLineColor(color[0]);
               DrawGammaSetMarker(histoRatioTrueEffiCut[i], 20, 1., 1, 1);
  	       histoRatioTrueEffiCut[i]->Reset("ICES");
               histoRatioTrueEffiCut[i]->DrawCopy("HIST");
      } else{
	
	
	

         if(i<20){
            DrawGammaSetMarker(histoRatioTrueEffiCut[i], 20+i, 1.,color[i],color[i]);
         } else {
            DrawGammaSetMarker(histoRatioTrueEffiCut[i], 20+i, 1.,color[i-20],color[i]);
         }
         histoRatioTrueEffiCut[i]->DrawCopy("same,e1,p");
	 fitRatioComparisonTrueEffi[i]->Draw("same");
      }

      DrawGammaLines(0., maxPt,1., 1.,0.1);
   }

   
   canvasTrueEffiMeson->Update();

   canvasTrueEffiMeson->SaveAs(Form("%s/%s_%s_Efficiencies.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix));
   delete canvasTrueEffiMeson;
   


   TCanvas* canvasGGFracCont = new TCanvas("canvasGGFracCont","",1350,1500);  // gives the page size
   canvasGGFracCont->SetTickx();
   canvasGGFracCont->SetTicky();
   canvasGGFracCont->SetGridx(0);
   canvasGGFracCont->SetGridy(0);
   canvasGGFracCont->SetLogy(0);
   canvasGGFracCont->SetLeftMargin(0.13);
   canvasGGFracCont->SetRightMargin(0.02);
   canvasGGFracCont->SetTopMargin(0.02);
   canvasGGFracCont->SetFillColor(0);

   TPad* padGGFracCont = new TPad("padGGFracCont", "", 0., 0.18, 1., 1.,-1, -1, -2);
   padGGFracCont->SetFillColor(0);
   padGGFracCont->GetFrame()->SetFillColor(0);
   padGGFracCont->SetBorderMode(0);
   padGGFracCont->SetLeftMargin(0.12);
   padGGFracCont->SetBottomMargin(0.);
   padGGFracCont->SetRightMargin(0.02);
   padGGFracCont->SetTopMargin(0.04);
   padGGFracCont->Draw();
   
       
   TPad* padGGFracContRatios = new TPad("padGGFracContRatios", "", 0., 0., 1., 0.18,-1, -1, -2);
   padGGFracContRatios->SetFillColor(0);
   padGGFracContRatios->GetFrame()->SetFillColor(0);
   padGGFracContRatios->SetBorderMode(0);
   padGGFracContRatios->SetLeftMargin(0.12);
   padGGFracContRatios->SetBottomMargin(0.3);
   padGGFracContRatios->SetRightMargin(0.02);
   padGGFracContRatios->SetTopMargin(0.01);
   padGGFracContRatios->Draw();

   padGGFracCont->cd();
   padGGFracCont->SetTickx();
   padGGFracCont->SetTicky();
   padGGFracCont->SetLogy(0);


   TLegend* legendGGFracCont = new TLegend(0.15,0.75,0.3,0.95);
   legendGGFracCont->SetTextSize(0.02);
   legendGGFracCont->SetFillColor(0);
   legendGGFracCont->SetLineColor(0);
   
   
   TLegend* legendGGFracContFits = new TLegend(0.15,0.50,0.3,0.70);
   legendGGFracContFits->SetTextSize(0.03);
   legendGGFracContFits->SetFillColor(0);
   legendGGFracContFits->SetLineColor(0);
   
   
   TF1*       fitRatioComparisonGGFracCont[NumberOfCuts];
   TPaveText* fitParamRatioComparisonGGFracCont[NumberOfCuts];

   Float_t minPtRatioGGFrac = 0.6;
   Float_t maxPtRatioGGFrac = 5.0;  
   Float_t fitParameterGGFrac=1.0; 
   
   
  
   
   
   for(Int_t i = 0; i< NumberOfCuts; i++){

      if(i == 0){
			
         DrawGammaSetMarker(histoGGFracCont[i], 20, 1., color[0], color[0]);
         DrawAutoGammaMesonHistos( histoGGFracCont[i],
                                   "", "p_{t} (GeV/c)","Contamination fraction",
				   
                                   kTRUE, 1., 1e-1,kFALSE,
                                   kTRUE, -0.1, 1.0,
                                   kFALSE, 0., 10.);
         histoGGFracCont[i]->GetYaxis()->SetRangeUser(-0.05,0.3);
         histoGGFracCont[i]->DrawCopy("e1,p");
         legendGGFracCont->AddEntry(histoGGFracCont[i],Form("standard: %s",cutStringsName[i].Data()));

      }
      else {
	
	fitRatioComparisonGGFracCont[i] = new TF1(Form("fitRatioComparisonGGFracCont%03d",i),"pol0");
	fitRatioComparisonGGFracCont[i]->SetParameter(0,fitParameterGGFrac);
	fitRatioComparisonGGFracCont[i]->SetRange(minPtRatioGGFrac,maxPtRatioGGFrac);
	histoRatioGGFracCont[i]->Fit(fitRatioComparisonGGFracCont[i],"SINRME+","",minPtRatioGGFrac,maxPtRatioGGFrac);
	legendGGFracContFits->AddEntry(fitRatioComparisonGGFracCont[i], Form("P_{0}: %f #pm %f",fitRatioComparisonGGFracCont[i]->GetParameter(0),fitRatioComparisonGGFracCont[i]->GetParError(0)));
  	
         if(i<20){
            DrawGammaSetMarker(histoGGFracCont[i], 20+i, 1.,color[i],color[i]);
	    fitRatioComparisonGGFracCont[i]->SetLineColor(color[i]);
            //fitRatioComparisonGGFracCont[i]->SetLineWidth(2);
         } else {
            DrawGammaSetMarker(histoGGFracCont[i], 20+i, 1.,color[i-20],color[i]);
	    fitRatioComparisonGGFracCont[i]->SetLineColor(color[i-20]);
	    //fitRatioComparisonGGFracCont[i]->SetLineWidth(2);	
         }

	   if( NumberOfCuts < 3 )fitRatioComparisonGGFracCont[i]->SetLineColor(kRed);		

          histoGGFracCont[i]->DrawCopy("same,e1,p");
         legendGGFracCont->AddEntry(histoGGFracCont[i],cutStringsName[i].Data());
      }
   }
     legendGGFracCont->Draw();
   //legendGGFracContFits->Draw("same");
   
  

  if (optionEnergy.CompareTo("pPb_5.02TeV") == 0){
      TLatex *labelCentrality = new TLatex(0.65,0.91,"p-Pb #sqrt{s_{NN}} = 5.02 TeV");
      SetStyleTLatex( labelCentrality, 0.038,4);
      labelCentrality->Draw();
  }

  padGGFracContRatios->cd();
  padGGFracContRatios->SetTickx();
  padGGFracContRatios->SetTicky();

  if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0) padGGFracContRatios->SetLogy(0);
  else padGGFracContRatios->SetLogy(0);

   for(Int_t i = 0; i< NumberOfCuts; i++){

      if(i==0){

         if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0){

               histoRatioGGFracCont[i]->SetYTitle(Form("R_{CP}(%s)",textMeson.Data()));
               histoRatioGGFracCont[i]->GetYaxis()->SetRangeUser(-0.1,2.0);

          } else {
               histoRatioGGFracCont[i]->SetYTitle("#frac{modified}{standard}");
               histoRatioGGFracCont[i]->GetYaxis()->SetRangeUser(-0.1,3.0);
          }
               histoRatioGGFracCont[i]->GetYaxis()->SetLabelSize(0.07);
               histoRatioGGFracCont[i]->GetYaxis()->SetNdivisions(505);
               histoRatioGGFracCont[i]->GetYaxis()->SetTitleSize(0.1);
               histoRatioGGFracCont[i]->GetYaxis()->SetDecimals();
               histoRatioGGFracCont[i]->GetYaxis()->SetTitleOffset(0.5);
               histoRatioGGFracCont[i]->GetXaxis()->SetTitleSize(0.11);
               histoRatioGGFracCont[i]->GetXaxis()->SetLabelSize(0.08);
               histoRatioGGFracCont[i]->SetMarkerStyle(22);
               histoRatioGGFracCont[i]->SetMarkerSize(1.);
               histoRatioGGFracCont[i]->SetMarkerColor(color[0]);
               histoRatioGGFracCont[i]->SetLineColor(color[0]);
               DrawGammaSetMarker(histoRatioGGFracCont[i], 20, 1., 1, 1);
               histoRatioGGFracCont[i]->Reset("ICES");
               histoRatioGGFracCont[i]->DrawCopy("HIST");
      } else{

         if(i<20){
            DrawGammaSetMarker(histoRatioGGFracCont[i], 20+i, 1.,color[i],color[i]);
         } else {
            DrawGammaSetMarker(histoRatioGGFracCont[i], 20+i, 1.,color[i-20],color[i]);
         }
           histoRatioGGFracCont[i]->DrawCopy("same,e1,p");
	 //fitRatioComparisonGGFracCont[i]->Draw("same");
	 
      }

      DrawGammaLines(0., maxPt,1., 1.,0.1);
   }

    canvasGGFracCont->Update();
   canvasGGFracCont->SaveAs(Form("%s/%s_%s_GGFracCont.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix));
   
   
   delete canvasGGFracCont;

}




