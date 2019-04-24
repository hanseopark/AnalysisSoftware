/****************************************************************************************************************************
******      Friederike Bock, friederike.bock@cern.ch                                      *****
******      Hikari Murakami, Pedro Gonzalez
******      Efficiency study 
******
*****************************************************************************************************************************/

#include <Riostream.h>
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
#include "TMinuit.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/CombinationFunctions.h"

extern TSystem* gSystem;


// run afterburner with data and from MC only added particles (specify path in baseDir)
// read CorrectedYieldTrueEff from Data afterburner output
// read MCYield_Meson_oldBin from MC afterurner output


//**********************************************************************************************************************
//***************************** Main function **************************************************************************
//**********************************************************************************************************************
void ExtractInputForWeightsPbPb(    TString suffix                  = "eps", 
				    TString fOptEnergy              = "PbPb_5.02TeV", 
				    Int_t iteration                 = 5
                            ){

  //**********************************************************************************************************************
  //****************************************** settings and declarations *************************************************
  //**********************************************************************************************************************

  Bool_t WriteFile    = kFALSE;
  Bool_t producePlots = kTRUE;
  Bool_t crossCheck   = kFALSE;

  gROOT->Reset();   
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();  
  SetPlotStyle();

  TString dateForOutput = ReturnDateStringForOutput();
  
  //TString outputFileName = Form("MCSpectraInputPbPb_%s.root",dateForOutput.Data());
  TString outputFileName = "MCSpectraInputPbPb_it5.root";

  TString fileNameTrainInput;
  if (iteration == 1)      fileNameTrainInput = "/home/meike/alice/AliPhysics/PWGGA/GammaConv/macros/data/MCSpectraInputPbPb_it1.root";
  else if (iteration == 2) fileNameTrainInput = "/home/meike/alice/AliPhysics/PWGGA/GammaConv/macros/data/MCSpectraInputPbPb_it2.root";
  else if (iteration == 3) fileNameTrainInput = "/home/meike/alice/AliPhysics/PWGGA/GammaConv/macros/data/MCSpectraInputPbPb_it3.root";
  else if (iteration == 5) fileNameTrainInput = "/home/meike/alice/AliPhysics/PWGGA/GammaConv/macros/data/MCSpectraInputPbPb_it4.root"; 
  
  // For calculation the ratios to the data-fit
  // bin content / fit function evaluated at bin center
  Bool_t integrateFunction = kFALSE;
  // if kTRUE: integrate function to calculate ratio

  Bool_t plotWeightedMCWithOldBin = kTRUE;
  
  const Int_t nCentClasses     = 9;
  const Int_t nCentClassesUsed = 9;
  TString baseDir;
    
  TString cent[nCentClasses]  = {"0-10%","0-20%","60-80%","0-5%","5-10%","10-20%","40-60%","20-40%","20-50%"};
  
  TString cutStringsAdded[nCentClasses] = {      
    "10110a23_00200009247602008250404000_0152501500000000",
    "10210a23_00200009247602008250404000_0152501500000000",
    "16810a23_00200009247602008250404000_0152501500000000",
    "30110a23_00200009247602008250404000_0152501500000000",
    "31210a23_00200009247602008250404000_0152501500000000",
    "11210a23_00200009247602008250404000_0152501500000000",
    "14610a23_00200009247602008250404000_0152501500000000",
    "12410a23_00200009247602008250404000_0152501500000000",
    "12510a23_00200009247602008250404000_0152501500000000" 
  };

  TString cutStrings[nCentClasses] = {
    "10110a13_00200009247602008250404000_0152501500000000",
    "10210a13_00200009247602008250404000_0152501500000000",
    "16810a13_00200009247602008250404000_0152501500000000",
    "30110a13_00200009247602008250404000_0152501500000000",
    "31210a13_00200009247602008250404000_0152501500000000",
    "11210a13_00200009247602008250404000_0152501500000000",
    "14610a13_00200009247602008250404000_0152501500000000",
    "12410a13_00200009247602008250404000_0152501500000000",
    "12510a13_00200009247602008250404000_0152501500000000"
  }; 

  // qcd, doubqcd, qmpt, rad, modkfunc, h, oHag, doHag
  TString fitFunctionsPi0[nCentClasses] = {"rad",      "rad",  "doHag", "doHag", "doHag", "rad", "oHag", "rad",   "doHag"};
  TString fitFunctionsEta[nCentClasses] = {"modkfunc", "oHag", "oHag",  "doHag", "doHag", "rad", "oHag", "doHag", "oHag"};
 
  const char* fitOptionPi0 = "NRME+";
  const char* fitOptionEta = "NRME+";
  // char* fitOptionEta = "WNRME+";  // option W to ignore error bars
  
  TString fileNameDataPi0[nCentClasses];
  TString fileNameDataEta[nCentClasses];
  TString fileNameAddedPi0[nCentClasses];
  TString fileNameAddedEta[nCentClasses];
  TString fileNameMBMCPi0[nCentClasses];
  TString fileNameMBMCEta[nCentClasses];
  //TString fileNameMergedMCPi0[nCentClasses];
  //TString fileNameMergedMCEta[nCentClasses];
   
  Color_t colorData           = kBlack;     //GetColorDefaultColor( fOptEnergy, "", cent[i], kFALSE);
  Color_t colorMC             = kGreen+2;    //GetColorDefaultColor( fOptEnergy, "", cent[i], kFALSE);
  Color_t colorMBMC           = kRed-4;
  Color_t colorFit            = kOrange+7;
  Color_t colorOldFit         = kBlue+1;
  
  const Int_t nMCs = 15;
  TString MCs[nMCs] = {"merged", "LHC18e1", "LHC18e1a", "LHC18e1b", "LHC18e1c", "LHC16h4", "LHC16i1a", "LHC16i1b", "LHC16i1c", "LHC16i2a", "LHC16i2b", "LHC16i2c", "LHC16i3a", "LHC16i3b", "LHC16i3c"};
 
  Double_t minPtPlot    = 0.3;
  Double_t maxPtPlot    = 20.0;
  Double_t minPtFitPi0  = 0.4;
  Double_t maxPtFitPi0  = 20.0;
  Double_t minPtFitEta  = 1.0;
  Double_t maxPtFitEta  = 10.0;

  TString collisionSystem              = ReturnFullCollisionsSystem(fOptEnergy);
  TString collisionSystemForWriting    = ReturnCollisionEnergyStringForTheory(fOptEnergy);
  
  Style_t markerStyleData    = 20;
  Size_t  markerSizeSpectrum = 1.0;
  Style_t lineStyleDataFit   = 1;
  Double_t lineWidthFit      = 1.0;
  
  Double_t markerSizeMC      = 1.0;
  Double_t markerSizeMBMC    = 1.0;
  Int_t markerStyleMC        = 20;
  Int_t markerStyleMBMC      = 24;
    
  TFile* filePi0DataInput[nCentClasses];
  TFile* fileEtaDataInput[nCentClasses];
  TFile* filePi0AddedInput[nCentClasses];
  TFile* fileEtaAddedInput[nCentClasses];
  TFile* filePi0MBMCInput[nCentClasses];
  TFile* fileEtaMBMCInput[nCentClasses];
  //TFile* filePi0MergedMC[nCentClasses];
  //TFile* fileEtaMergedMC[nCentClasses];
  
  TH1D* histoYieldDataPi0[nCentClasses];                // Data histo
  TH1D* histoYieldDataEta[nCentClasses];
  TGraphAsymmErrors* graphYieldDataPi0[nCentClasses];   // Data graph
  TGraphAsymmErrors* graphYieldDataEta[nCentClasses];
  TF1* fitPi0DataYield[nCentClasses];                   // Data fit
  TF1* fitEtaDataYield[nCentClasses]; 
  
  TH1D* histoPi0InputMCWOWeights[nCentClasses];       // MC histo without weights (added particles)
  TH1D* histoPi0InputMCWWeights[nCentClasses];        // MC histo with    weights
  //TH1D* histoPi0Efficiency[nCentClasses];
  
  TH1D* histoPi0InputMBMCWOWeights[nCentClasses];     // MC histo without weights (normal particles)
  TH1D* histoPi0InputMBMCWWeights[nCentClasses];      // MC histo with    weights
  //TH1D* histoPi0MBEfficiency[nCentClasses];

  TH1D* histoPi0MergedEfficiency[nCentClasses];       // effi merged from weighted 
  
  TH1D* histoEtaInputMCWOWeights[nCentClasses];
  TH1D* histoEtaInputMCWWeights[nCentClasses];
  //TH1D* histoEtaEfficiency[nCentClasses];
  TH1D* histoEtaInputMBMCWOWeights[nCentClasses];  
  TH1D* histoEtaInputMBMCWWeights[nCentClasses];     
  //TH1D* histoEtaMBEfficiency[nCentClasses];  
  TH1D* histoEtaMergedEfficiency[nCentClasses];       

  
  TH1D* histoPi0RatioDataToFit[nCentClasses];   // Data histo to new Data fit
  TH1D* histoEtaRatioDataToFit[nCentClasses];

  TH1D*  histoPi0RatioDataToOldFit[nCentClasses];   // Data histo to old Data fit
  TH1D*  histoEtaRatioDataToOldFit[nCentClasses];   
    
  TH1D* histoPi0Ratio[nCentClasses];            // ratio MC(Added) without weight histo to Data fit
  TH1D* histoEtaRatio[nCentClasses];            
  
  TH1D* histoPi0MBRatio[nCentClasses];            // ratio MC(Added) without weight histo to Data fit
  TH1D* histoEtaMBRatio[nCentClasses];

  TFile *trainInputFile[nCentClasses];           // file which was used for weighting on lego train
  TF1 *fitPi0FromFile[nCentClasses];             // fits used for weighting
  TF1 *fitEtaFromFile[nCentClasses];
  
  TString folderName = "";
  

  if(crossCheck){

    for(Int_t i=0; i<nCentClassesUsed; i++){
    
      TString cutNumberAdded  = cutStringsAdded[i];
      TString eventCutAdded   = cutNumberAdded(0,8);
      TString cutNumber       = cutStrings[i];
      TString eventCut        = cutNumber(0,8);   
      TString eventCutShort   = cutNumber(0,6);
      
      Double_t functionResultMC = 1.;
      Double_t mesonPt = 0.;  // would have to go through all values in small steps smaller than bin width
      
      trainInputFile[i] = new TFile(fileNameTrainInput);

      // read data
      fitPi0FromFile[i] = (TF1*)trainInputFile[i]->Get(Form("Pi0_Data_%s_%s", collisionSystemForWriting.Data(), eventCutShort.Data()));

      for(Int_t j=0; j<nMCs; j++){

	// read MB MC
	histoPi0InputMBMCWOWeights[i] = (TH1D*)trainInputFile[i]->Get(Form("Pi0_%s_%s_%s",MCs[j].Data(),collisionSystemForWriting.Data(),eventCut.Data()));
	// read added particles MC
	//histoPi0InputMCWOWeights[i]   = (TH1D*)trainInputFile[i]->Get(Form("Pi0_LHC16h4_%s_%s",collisionSystemForWriting.Data(),eventCutAdded.Data()));

	// get MC value
	functionResultMC = histoPi0InputMBMCWOWeights[i]->Interpolate(mesonPt);

	
      } // end of MCs loop
    } // end of cent classes loop
  } // end of cross check
  
//**********************************************************************************************************************
//************************************************* writeFile **********************************************************
//**********************************************************************************************************************

  if(WriteFile){

    Bool_t writeData = kTRUE;    // input from previous iteration
    Bool_t writeMC   = kFALSE;    // needed only before the first iteration

    cout << "Write file..." << endl;
      
    TFile fMCSpectraInput(outputFileName, "UPDATE");
  
    for(Int_t i=0; i<nCentClassesUsed; i++){
    
      TString cutNumberAdded  = cutStringsAdded[i];
      TString eventCutAdded   = cutNumberAdded(0,8);
      TString cutNumber       = cutStrings[i];
      TString eventCut        = cutNumber(0,8);   
      TString eventCutShort   = cutNumber(0,6);
      TString periodNameMC;
    
      cout << "writing everything for " << cent[i].Data() << endl;
    
      //**********************************************************************************************************************
      //************************************** read in data and MC files Pi0 and Eta *****************************************
      //**********************************************************************************************************************

      // +++ write data +++
    
      if(writeData){

	if(iteration == 0){        // for first iteration
	  baseDir               = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeights"; 
	  folderName            = "afterburner_AOD_MCmerged";
	} else if (iteration == 1){// have first iteration, for second iteration
	  baseDir               = "/home/meike/analysis/results/photonconvResults/PbPb/final";       
	  folderName            = "afterburner_AOD_MCmergedWithAddSig";
	} else if (iteration == 3){  // add cent classes 20-40% and 20-50%
	  baseDir               = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeights3rd";       
	  folderName            = "afterburner_AOD_MCLHC16g1";  // could not merge with other MCs because of light output complication
	} else if (iteration == 4){  // added cent class 0-20% (now still without mult weights), new 18e1* MCs, new pT binning
	  baseDir               = "/home/meike/analysis/results/photonconvResults/PbPb/newMCs";       
	  folderName            = "afterburner_AOD_MCmerged";
	} else if (iteration == 5){
	  baseDir               = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeightsNewMCs";
	  folderName            = "";
	} else
	  cout << "Iteration not correctly defined" << endl;
	    
	fileNameDataPi0[i]    = Form("%s/%s/%s/%s/Pi0_data_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
	fileNameDataEta[i]    = Form("%s/%s/%s/%s/Eta_data_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
	filePi0DataInput[i]            = new TFile(fileNameDataPi0[i]);
	fileEtaDataInput[i]            = new TFile(fileNameDataEta[i]);
	histoYieldDataPi0[i]           = (TH1D*)filePi0DataInput[i]->Get("CorrectedYieldTrueEff");
	histoYieldDataEta[i]           = (TH1D*)fileEtaDataInput[i]->Get("CorrectedYieldTrueEff");
	// convert histo to graph
	graphYieldDataPi0[i]           = new TGraphAsymmErrors(histoYieldDataPi0[i]);
	graphYieldDataEta[i]           = new TGraphAsymmErrors(histoYieldDataEta[i]);
	// fit data
	fitPi0DataYield[i] = FitObject(fitFunctionsPi0[i].Data(),"fitPi0DataYield","Pi0", NULL, minPtFitPi0, maxPtFitPi0);

	//fitPi0DataYield[i]->SetParameter(i,x);
	
	graphYieldDataPi0[i]->Fit(fitPi0DataYield[i],fitOptionPi0, fitOptionPi0, minPtFitPi0, maxPtFitPi0);
	//fitPi0DataYield[i]->SetRange(minPtFitPi0, maxPtFitPi0);
	fitPi0DataYield[i]->SetRange(minPtPlot, maxPtPlot);
	fitPi0DataYield[i]->SetLineColor(colorFit);
	fitPi0DataYield[i]->SetLineStyle(lineStyleDataFit);
	fitPi0DataYield[i]->SetLineWidth(lineWidthFit);
	fitEtaDataYield[i]     = FitObject(fitFunctionsEta[i].Data(),"fitEtaDataYield","Eta", NULL, minPtFitEta, maxPtFitEta);
	graphYieldDataEta[i]->Fit(fitEtaDataYield[i],fitOptionEta, fitOptionEta, minPtFitEta, maxPtFitEta);
	//fitEtaDataYield[i]->SetRange(minPtFitEta, maxPtFitEta);
	fitEtaDataYield[i]->SetRange(minPtPlot, maxPtPlot);
	fitEtaDataYield[i]->SetLineColor(colorFit);
	fitEtaDataYield[i]->SetLineStyle(lineStyleDataFit);
	fitEtaDataYield[i]->SetLineWidth(lineWidthFit);
    
	// write to file
	fMCSpectraInput.cd();
	if(fitPi0DataYield[i]){
	  fitPi0DataYield[i]->SetTitle(Form("Pi0_Data_%s_%s", collisionSystemForWriting.Data(), eventCutShort.Data()));
	  fitPi0DataYield[i]->Write(Form("Pi0_Data_%s_%s", collisionSystemForWriting.Data(), eventCutShort.Data()),TObject::kOverwrite);
	} else cout << "ERROR "  << endl;
	if(fitEtaDataYield[i]){
	  cout << "Eta fit " << i << " " << fitEtaDataYield[i]->GetName() << endl;
	  fitEtaDataYield[i]->SetTitle(Form("Eta_Data_%s_%s", collisionSystemForWriting.Data(), eventCutShort.Data()));
	  fitEtaDataYield[i]->Write(Form("Eta_Data_%s_%s", collisionSystemForWriting.Data(), eventCutShort.Data()),TObject::kOverwrite);
	} else cout << "ERROR "  << endl;

	//cleanup
	delete filePi0DataInput[i];
	delete graphYieldDataPi0[i];
	delete graphYieldDataEta[i];
	delete fileEtaDataInput[i];
      }  // end of writeData

      if(writeMC){
		
	// +++ added MC particles +++++++++++++++++++++++++++++++++++++++++++++++++
	
	if(iteration==0){
	  
	  baseDir              = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeights";
	  folderName             = "afterburner_AOD_MCLHC16h4";
	  
	  fileNameAddedPi0[i]    = Form("%s/%s/%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStringsAdded[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStringsAdded[i].Data());
	  fileNameAddedEta[i]    = Form("%s/%s/%s/%s/Eta_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStringsAdded[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStringsAdded[i].Data());
	  
	  filePi0AddedInput[i]         = new TFile(fileNameAddedPi0[i]);
	  fileEtaAddedInput[i]         = new TFile(fileNameAddedEta[i]);
	  
	  histoPi0InputMCWOWeights[i]    = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson_oldBinWOWeights");
	  histoEtaInputMCWOWeights[i]    = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson_oldBinWOWeights");

	  // cleanup
	  //delete filePi0AddedInput[i];
	  //delete fileEtaAddedInput[i];
	  // segmentation violation ?!
	  
	} else if (iteration==3){

	  baseDir    = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeights3rd";
	  folderName = "afterburner_AOD_MCLHC16h4WithAddSig";  // did run MB and added particles within the same afterburner run

	  fileNameAddedPi0[i]    = Form("%s/%s/%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
	  fileNameAddedEta[i]    = Form("%s/%s/%s/%s/Eta_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
	  
	  filePi0AddedInput[i]         = new TFile(fileNameAddedPi0[i]);  
	  fileEtaAddedInput[i]         = new TFile(fileNameAddedEta[i]);
	  
	  histoPi0InputMCWOWeights[i]    = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson_oldBinWOWeights_AddedSig");  
	  histoEtaInputMCWOWeights[i]    = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson_oldBinWOWeights_AddedSig");

	} else if (iteration==4){

	  baseDir    = "/home/meike/analysis/results/photonconvResults/PbPb/newMCs";
	  folderName = "afterburner_AOD_MCLHC16h4WithAddSig";  // code below is the same as for it 3

	  fileNameAddedPi0[i]    = Form("%s/%s/%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
	  fileNameAddedEta[i]    = Form("%s/%s/%s/%s/Eta_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
	  
	  filePi0AddedInput[i]         = new TFile(fileNameAddedPi0[i]);  
	  fileEtaAddedInput[i]         = new TFile(fileNameAddedEta[i]);
	  
	  histoPi0InputMCWOWeights[i]    = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson_oldBinWOWeights_AddedSig");  
	  histoEtaInputMCWOWeights[i]    = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson_oldBinWOWeights_AddedSig");
	  
	} 
	
	// write to file
	fMCSpectraInput.cd();
	histoPi0InputMCWOWeights[i]->SetTitle(Form("Pi0_LHC16h4_%s_%s", collisionSystemForWriting.Data(),eventCutAdded.Data()));
	histoPi0InputMCWOWeights[i]->Write(Form("Pi0_LHC16h4_%s_%s", collisionSystemForWriting.Data(),eventCutAdded.Data()),TObject::kOverwrite);
	histoEtaInputMCWOWeights[i]->SetTitle(Form("Eta_LHC16h4_%s_%s", collisionSystemForWriting.Data(),eventCutAdded.Data()));
	histoEtaInputMCWOWeights[i]->Write(Form("Eta_LHC16h4_%s_%s", collisionSystemForWriting.Data(),eventCutAdded.Data()),TObject::kOverwrite);

	cout << "done writing added particles" << endl;
	
	// +++  normal MC particles ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
	for(Int_t j=1; j<nMCs; j++){    // loop over all MCs excluding 'merged' and including 'LHC16h4'

	  if(iteration==0){
	    baseDir              = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeights";
	    folderName            = Form("afterburner_AOD_MC%s", MCs[j].Data());
	  } else if(iteration==3){
	    baseDir              = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeights3rd";
	    if(MCs[j].CompareTo("LHC16h4")==0) folderName = "afterburner_AOD_MCLHC16h4WithAddSig";    // did run MB and added particles within the same afterburner run
	    else                               folderName = Form("afterburner_AOD_MC%s", MCs[j].Data());
	  } else if(iteration==4){
	    baseDir              = "/home/meike/analysis/results/photonconvResults/PbPb/newMCs";
	    if(MCs[j].CompareTo("LHC16h4")==0) folderName = "afterburner_AOD_MCLHC16h4WithAddSig";    
	    else                               folderName = Form("afterburner_AOD_MC%s", MCs[j].Data());
	  } 
	  
	  fileNameMBMCPi0[i]    = Form("%s/%s/%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
	  fileNameMBMCEta[i]    = Form("%s/%s/%s/%s/Eta_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
	  filePi0MBMCInput[i]         = new TFile(fileNameMBMCPi0[i]);
	  fileEtaMBMCInput[i]         = new TFile(fileNameMBMCEta[i]);
	  histoPi0InputMBMCWOWeights[i]  = (TH1D*)filePi0MBMCInput[i]->Get("MCYield_Meson_oldBinWOWeights");
	  histoEtaInputMBMCWOWeights[i]  = (TH1D*)fileEtaMBMCInput[i]->Get("MCYield_Meson_oldBinWOWeights");
	  
	  // write to file
	  fMCSpectraInput.cd();
	  if(histoPi0InputMBMCWOWeights[i]){
	    histoPi0InputMBMCWOWeights[i]->SetTitle(Form("Pi0_%s_%s_%s",MCs[j].Data(),collisionSystemForWriting.Data(),eventCut.Data()));
	    histoPi0InputMBMCWOWeights[i]->Write(Form("Pi0_%s_%s_%s",MCs[j].Data(),collisionSystemForWriting.Data(),eventCut.Data()),TObject::kOverwrite);
	  } else cout << "ERROR: histo does not exist" << endl;
	  if(histoEtaInputMBMCWOWeights[i]){
	    histoEtaInputMBMCWOWeights[i]->SetTitle(Form("Eta_%s_%s_%s",MCs[j].Data(),collisionSystemForWriting.Data(),eventCut.Data()));
	    histoEtaInputMBMCWOWeights[i]->Write(Form("Eta_%s_%s_%s",MCs[j].Data(),collisionSystemForWriting.Data(),eventCut.Data()),TObject::kOverwrite);
	  } else cout << "ERROR: histo does not exist" << endl;

	  // cleanup
	  delete filePi0MBMCInput[i];
	  delete fileEtaMBMCInput[i];
      
	} // end of loop over MCs
      } // end of writeMC
    } // end of loop over cent classes
    
    fMCSpectraInput.Close();
    cout << "Done writing to file" << endl;

  } // end of writeFile
    
//**********************************************************************************************************************
//******************************************************* producePlots *************************************************
//**********************************************************************************************************************

  if(producePlots){
  
    for(Int_t i=0; i<nCentClassesUsed; i++){

      cout << endl;    
      cout << "========= Producing plots for " << cent[i].Data() << endl;
      cout << endl;
      
      TString cutNumberAdded  = cutStringsAdded[i];
      TString eventCutAdded   = cutNumberAdded(0,8);
      TString cutNumber       = cutStrings[i];
      TString eventCut        = cutNumber(0,8);   
      TString eventCutShort   = cutNumber(0,6);

      cout << eventCutShort << endl;
            
    //**********************************************************************************************************************
    //************************************** read in data and MC files Pi0 and Eta *****************************************
    //**********************************************************************************************************************
    //if ( iteration == 0 || iteration == 1 ) 
      // baseDir                = "/home/meike/analysis/results/photonconvResults/PbPb/final";
      //else if ( iteration == 2 )
      if(iteration < 3)        baseDir = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeights2nd";
      else if(iteration == 3)  baseDir = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeights3rd";
      else if(iteration == 4)  baseDir = "/home/meike/analysis/results/photonconvResults/PbPb/newMCs";
      else if(iteration == 5)  baseDir = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeightsNewMCs";
      
    // +++ data +++
    if (iteration == 0)      folderName = "afterburner_AOD_MCmerged_WoPTW";              // corrected with effi from normal un-weighted MC
    //else if (iteration == 1) folderName = "afterburner_AOD_MCmergedWithAddSig";        // corrected with effi from normal and added particles merged, both with pT-weighting
    else if (iteration == 1) folderName = "afterburner_AOD_MCmergedWithAddSig_WPTW1";    // corrected with effi from normal and added particles merged, both with pT-weighting
    else if (iteration == 2) folderName = "afterburner_AOD_MCmergedWithAddSig_WPTW2";    // corrected with effi from normal and added particles merged, both with pT-weighting
    else if (iteration == 3) folderName = "afterburner_AOD_MCLHC16g1";
    else if (iteration == 4) folderName = "afterburner_AOD_MCmerged";                    // corrected with effi from normal un-weighted MC
    else if (iteration == 5) folderName = "";                                            // corrected with effi from normal and added particles merged, both with pT-weighting
    
    fileNameDataPi0[i]    = Form("%s/%s/%s/%s/Pi0_data_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
    fileNameDataEta[i]    = Form("%s/%s/%s/%s/Eta_data_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
    // efficiencies from the same file; normal and added MC particles merged
    //fileNameMergedMCPi0[i] = fileNameDataPi0[i];
    //fileNameMergedMCEta[i] = fileNameDataEta[i];

    
    // +++ normal MC particles +++ 
    if ( iteration == 0 )     folderName = "afterburner_AOD_MCmerged_WoPTW";
    else if (iteration == 1 ) folderName = "afterburner_AOD_MCmerged_WPTW1";  //folderName  = "afterburner_AOD_MCmerged";
    else if (iteration == 2)  folderName = "afterburner_AOD_MCmerged_WPTW2";
    else if (iteration == 3)  folderName = "afterburner_AOD_MCLHC16g1";
    else if (iteration == 4)  folderName = "afterburner_AOD_MCmerged";
    else if (iteration == 5)  folderName = "";                                 // read normal and added particles from the same file
    
    fileNameMBMCPi0[i]    = Form("%s/%s/%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
    fileNameMBMCEta[i]    = Form("%s/%s/%s/%s/Eta_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());

    // +++ added particles +++ 
    if ( iteration == 0 )    folderName = "afterburner_AOD_MCLHC16h4_WoPTW";
    else if (iteration == 1) folderName = "afterburner_AOD_MCLHC16h4_WPTW1";      // folderName = "afterburner_AOD_MCLHC16h4";
    else if (iteration == 2) folderName = "afterburner_AOD_MCLHC16h4_WPTW2";    
    else if (iteration == 3) folderName = "afterburner_AOD_MCLHC16h4WithAddSig";  
    else if (iteration == 4) folderName = "afterburner_AOD_MCLHC16h4WithAddSig";
    else if (iteration == 5)  folderName = "";                                   // read normal and added particles from the same file
    
    if( iteration < 3){
      fileNameAddedPi0[i]    = Form("%s/%s/%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStringsAdded[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStringsAdded[i].Data());
      fileNameAddedEta[i]    = Form("%s/%s/%s/%s/Eta_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStringsAdded[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStringsAdded[i].Data());
    } else {
      fileNameAddedPi0[i]    = Form("%s/%s/%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
      fileNameAddedEta[i]    = Form("%s/%s/%s/%s/Eta_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStrings[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStrings[i].Data());
    }
    // read fit from lego train input file

    if(iteration == 1 || iteration == 2 || iteration == 3 || iteration == 5){
      trainInputFile[i] = new TFile(fileNameTrainInput);
      fitPi0FromFile[i] = (TF1*)trainInputFile[i]->Get(Form("Pi0_Data_%s_%s", collisionSystemForWriting.Data(), eventCutShort.Data()));
      fitEtaFromFile[i] = (TF1*)trainInputFile[i]->Get(Form("Eta_Data_%s_%s", collisionSystemForWriting.Data(), eventCutShort.Data()));
    }
    
    filePi0DataInput[i]          = new TFile(fileNameDataPi0[i]);
    fileEtaDataInput[i]          = new TFile(fileNameDataEta[i]);
    filePi0AddedInput[i]         = new TFile(fileNameAddedPi0[i]);
    fileEtaAddedInput[i]         = new TFile(fileNameAddedEta[i]);
    filePi0MBMCInput[i]          = new TFile(fileNameMBMCPi0[i]);
    fileEtaMBMCInput[i]          = new TFile(fileNameMBMCEta[i]);
    //filePi0MergedMC[i]           = new TFile(fileNameMergedMCPi0[i]);
    //fileEtaMergedMC[i]           = new TFile(fileNameMergedMCEta[i]);
    
    histoYieldDataPi0[i]            = (TH1D*)filePi0DataInput[i]->Get("CorrectedYieldTrueEff");

    if(iteration == 0){
      histoPi0InputMCWOWeights[i]    = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson_oldBinWOWeights");   // there is no MCYield_Meson_WOWeights
      histoPi0InputMBMCWOWeights[i]   = (TH1D*)filePi0MBMCInput[i]->Get("MCYield_Meson_oldBinWOWeights");  
    } else if(iteration==1 || iteration==2){
      if(plotWeightedMCWithOldBin){
	histoPi0InputMCWWeights[i]     = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson_oldBin");    
	histoPi0InputMBMCWWeights[i]    = (TH1D*)filePi0MBMCInput[i]->Get("MCYield_Meson_oldBin");
      } else {
	histoPi0InputMCWWeights[i]     = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson");
	histoPi0InputMBMCWWeights[i]    = (TH1D*)filePi0MBMCInput[i]->Get("MCYield_Meson");
      }
    } else if(iteration == 3  || iteration == 4){
      histoPi0InputMCWOWeights[i]    = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson_oldBinWOWeights_AddedSig"); 
      histoPi0InputMBMCWOWeights[i]   = (TH1D*)filePi0MBMCInput[i]->Get("MCYield_Meson_oldBinWOWeights");  
    } else if (iteration==5){
      if(plotWeightedMCWithOldBin){
	histoPi0InputMCWWeights[i]     = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson_oldBin_AddedSig");    
	histoPi0InputMBMCWWeights[i]    = (TH1D*)filePi0MBMCInput[i]->Get("MCYield_Meson_oldBin");
      } else {
	cout << "Not working: there is no MCYield_Meson_AddedSig in file" << endl;
	return;
      }
    }

    histoYieldDataEta[i]            = (TH1D*)fileEtaDataInput[i]->Get("CorrectedYieldTrueEff");

    if(iteration == 0){
      histoEtaInputMCWOWeights[i]    = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson_oldBinWOWeights");   
      histoEtaInputMBMCWOWeights[i]   = (TH1D*)fileEtaMBMCInput[i]->Get("MCYield_Meson_oldBinWOWeights");  
    } else if(iteration==1 || iteration==2 || iteration == 5){
      if(plotWeightedMCWithOldBin){
	histoEtaInputMCWWeights[i]     = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson_oldBin");    
	histoEtaInputMBMCWWeights[i]    = (TH1D*)fileEtaMBMCInput[i]->Get("MCYield_Meson_oldBin"); 
      } else{
	histoEtaInputMCWWeights[i]     = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson");    
	histoEtaInputMBMCWWeights[i]    = (TH1D*)fileEtaMBMCInput[i]->Get("MCYield_Meson");
      }
    } else if(iteration == 3  || iteration == 4){
      histoEtaInputMCWOWeights[i]    = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson_oldBinWOWeights_AddedSig"); 
      histoEtaInputMBMCWOWeights[i]   = (TH1D*)fileEtaMBMCInput[i]->Get("MCYield_Meson_oldBinWOWeights");  
    } else if (iteration==5){
      if(plotWeightedMCWithOldBin){
	histoEtaInputMCWWeights[i]     = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson_oldBin_AddedSig");    
	histoEtaInputMBMCWWeights[i]    = (TH1D*)fileEtaMBMCInput[i]->Get("MCYield_Meson_oldBin");
      } else {
	cout << "Not working: there is no MCYield_Meson_AddedSig in file" << endl;
	return;
      }
    }

    
    // convert histo to graph
    graphYieldDataPi0[i]           = new TGraphAsymmErrors(histoYieldDataPi0[i]);
    graphYieldDataEta[i]           = new TGraphAsymmErrors(histoYieldDataEta[i]);
    
    // fit data
    cout << endl;
    cout << "=========== Fit Pi0 =============" << endl;
    cout << endl;
    fitPi0DataYield[i] = FitObject(fitFunctionsPi0[i].Data(),"fitPi0DataYield","Pi0", NULL, minPtFitPi0, maxPtFitPi0);
    graphYieldDataPi0[i]->Fit(fitPi0DataYield[i], fitOptionPi0, fitOptionPi0, minPtFitPi0, maxPtFitPi0);
    //fitPi0DataYield[i]->SetRange(minPtFitPi0, maxPtFitPi0);
    fitPi0DataYield[i]->SetRange(minPtPlot, maxPtPlot);
    fitPi0DataYield[i]->SetLineColor(colorFit);
    fitPi0DataYield[i]->SetLineStyle(lineStyleDataFit);
    fitPi0DataYield[i]->SetLineWidth(lineWidthFit);
    cout << endl;
    cout << "=========== Fit Eta =============" << endl;
    cout << endl;
    fitEtaDataYield[i]     = FitObject(fitFunctionsEta[i].Data(),"fitEtaDataYield","Eta", NULL, minPtFitEta, maxPtFitEta);
    graphYieldDataEta[i]->Fit(fitEtaDataYield[i],fitOptionEta, fitOptionEta, minPtFitEta, maxPtFitEta);
    //fitEtaDataYield[i]->SetRange(minPtFitEta, maxPtFitEta);
    fitEtaDataYield[i]->SetRange(minPtPlot, maxPtPlot);
    fitEtaDataYield[i]->SetLineColor(colorFit);
    fitEtaDataYield[i]->SetLineStyle(lineStyleDataFit);
    fitEtaDataYield[i]->SetLineWidth(lineWidthFit);

    cout << endl;
    cout << "======== Calculating ratios ==========" << endl;
    cout << endl;
    
    //if(histoPi0InputMCWOWeights[i]) cout << "Pi0 histo ok" << endl;
    //if(fitPi0DataYield[i])          cout << "Pi0 fit   ok" << endl;
    
    // calculate ratio added signals MC histo to Data fit
    if(iteration == 0 || iteration == 3 || iteration == 4){
      histoPi0Ratio[i]  = CalculateHistoRatioToFit(histoPi0InputMCWOWeights[i], fitPi0DataYield[i],integrateFunction);
      histoEtaRatio[i]  = CalculateHistoRatioToFit(histoEtaInputMCWOWeights[i], fitEtaDataYield[i],integrateFunction);
    } else{
      histoPi0Ratio[i]  = CalculateHistoRatioToFit(histoPi0InputMCWWeights[i], fitPi0DataYield[i],integrateFunction);
      histoEtaRatio[i]  = CalculateHistoRatioToFit(histoEtaInputMCWWeights[i], fitEtaDataYield[i],integrateFunction);
    }

    // calculate ratio normal MC particles histo to Data fit
    if(iteration == 0 || iteration == 3 || iteration == 4){
      histoPi0MBRatio[i]  = CalculateHistoRatioToFit(histoPi0InputMBMCWOWeights[i], fitPi0DataYield[i],integrateFunction);
      histoEtaMBRatio[i]  = CalculateHistoRatioToFit(histoEtaInputMBMCWOWeights[i], fitEtaDataYield[i],integrateFunction);
    } else {
      histoPi0MBRatio[i]  = CalculateHistoRatioToFit(histoPi0InputMBMCWWeights[i], fitPi0DataYield[i],integrateFunction);
      histoEtaMBRatio[i]  = CalculateHistoRatioToFit(histoEtaInputMBMCWWeights[i], fitEtaDataYield[i],integrateFunction);
    }

    // calculate ratio Data histo to new fit
    histoPi0RatioDataToFit[i] = CalculateHistoRatioToFit(histoYieldDataPi0[i], fitPi0DataYield[i], integrateFunction);
    histoEtaRatioDataToFit[i] = CalculateHistoRatioToFit(histoYieldDataEta[i], fitEtaDataYield[i], integrateFunction);

    // calculate ratio Data histo to old fit
    if(iteration == 5){
      histoPi0RatioDataToOldFit[i] = CalculateHistoRatioToFit(histoYieldDataPi0[i], fitPi0FromFile[i], integrateFunction);
      histoEtaRatioDataToOldFit[i] = CalculateHistoRatioToFit(histoYieldDataEta[i], fitEtaFromFile[i], integrateFunction);
    }
    
    } // end of loop over cent classes
  


  //**********************************************************************************************************************
  //************************************************** PLOTTING **********************************************************
  //**********************************************************************************************************************
    cout << endl;
    cout << "============ Plotting ================" << endl;
    cout << endl;
    
    TString outputDir               = Form("../%s/ExtractInputForWeights/%s",suffix.Data(),dateForOutput.Data());
    cout << "output dir: " << outputDir.Data() << endl;
    gSystem->Exec("mkdir -p "+outputDir);
    
    // canvas
    TCanvas* canvas  = new TCanvas("canvas", "", 200, 10, 1500, 500);  // gives the page size
    
    // legends
    TLegend* legendSpectra            = GetAndSetLegend2(0.2, 0.15, 0.8, 0.15+0.035*6, 0.035, 1, "", 42, 0.15);
    TLegend* legendRatio1             = GetAndSetLegend2(0.45, 0.15, 0.9, 0.15+0.035*3, 0.035, 1, "", 42, 0.15);
    TLegend* legendRatio              = GetAndSetLegend2(0.45, 0.15, 0.9, 0.15+0.035*3, 0.035, 1, "", 42, 0.15);
    
    // labels
    TLatex *labelSpectraPi0Label        = new TLatex(0.925,0.87,"#pi^{0} #rightarrow #gamma #gamma");
    SetStyleTLatex( labelSpectraPi0Label, 0.035, 4, 1, 42, kTRUE, 31);
    TLatex *labelSpectraEtaLabel        = new TLatex(0.925,0.87,"#eta #rightarrow #gamma #gamma");
    SetStyleTLatex( labelSpectraEtaLabel, 0.035, 4, 1, 42, kTRUE, 31);
    TLatex *labelSpectraPi0LabelRatio   = new TLatex(0.17,0.87,"#pi^{0} #rightarrow #gamma #gamma");
    SetStyleTLatex( labelSpectraPi0LabelRatio, 0.035, 4, 1, 42, kTRUE, 11);
    TLatex *labelSpectraEtaLabelRatio   = new TLatex(0.17,0.87,"#eta #rightarrow #gamma #gamma");
    SetStyleTLatex( labelSpectraEtaLabelRatio, 0.035, 4, 1, 42, kTRUE, 11);

    //TLatex *labelIteration = new TLatex(0.2, 0.15+0.035*6+0.02, Form("After iteration %d",iteration));
    //SetStyleTLatex(labelIteration, 0.035, 4, 1, 42, kTRUE, 11);
    
    // default histos
    TH1F * histo1DSpectra;
    histo1DSpectra          = new TH1F("histo1DSpectra", "histo1DSpectra",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( histo1DSpectra, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
                            0.03, 0.035, 0.03, 0.035, 0.83, 2.0);
    histo1DSpectra->GetYaxis()->SetRangeUser(1E-8, 5E3);
    histo1DSpectra->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DSpectra->GetXaxis()->SetLabelOffset(-0.01);
    histo1DSpectra->GetYaxis()->SetLabelOffset(0.01);

    // for ratio data/fit
    TH1F * histo1DRatio;
    histo1DRatio          = new TH1F("histo1DRatio", "histo1DRatio",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( histo1DRatio, "#it{p}_{T} (GeV/#it{c})", "Ratio data histo to data fit",
                            0.03, 0.035, 0.03, 0.035, 0.83, 1.6);
    histo1DRatio->GetYaxis()->SetRangeUser(0.61, 1.5);
    histo1DRatio->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DRatio->GetXaxis()->SetLabelOffset(-0.01);
    histo1DRatio->GetYaxis()->SetLabelOffset(0.01);

    // for ratio MC/data
    TH1F * histo1DRatio2Pi0;
    histo1DRatio2Pi0          = new TH1F("histo1DRatio2Pi0", "histo1DRatio2Pi0",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( histo1DRatio2Pi0, "#it{p}_{T} (GeV/#it{c})", "Ratio MC to new data fit",
			       0.03, 0.035, 0.03, 0.035, 0.83, 1.5);
    if(iteration == 0 || iteration == 3 || iteration == 4) histo1DRatio2Pi0->GetYaxis()->SetRangeUser(1E-4, 5E5);  
    else               histo1DRatio2Pi0->GetYaxis()->SetRangeUser(0.85,1.2);  
    histo1DRatio2Pi0->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DRatio2Pi0->GetXaxis()->SetLabelOffset(-0.01);
    histo1DRatio2Pi0->GetYaxis()->SetLabelOffset(0.01);

    TH1F * histo1DRatio2Eta;
    histo1DRatio2Eta          = new TH1F("histo1DRatio2Eta", "histo1DRatio2Eta",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( histo1DRatio2Eta, "#it{p}_{T} (GeV/#it{c})", "Ratio MC to new data fit",
			       0.03, 0.035, 0.03, 0.035, 0.83, 1.5);
    if(iteration == 0 || iteration == 3 || iteration == 4) histo1DRatio2Eta->GetYaxis()->SetRangeUser(0.00000005, 900000000.);   
    else               histo1DRatio2Eta->GetYaxis()->SetRangeUser(0.1,2.9);  
    histo1DRatio2Eta->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DRatio2Eta->GetXaxis()->SetLabelOffset(-0.01);
    histo1DRatio2Eta->GetYaxis()->SetLabelOffset(0.01);
    
    for(Int_t i=0; i<nCentClassesUsed; i++){

      // labels
      TLatex *labelCollisionSystem        = new TLatex(0.93,0.91,Form("%s  %s", collisionSystem.Data(), cent[i].Data()));
      SetStyleTLatex( labelCollisionSystem, 0.035, 4, 1, 42, kTRUE, 31);
      TLatex *labelCollisionSystemRatio     = new TLatex(0.17,0.91,Form("%s  %s", collisionSystem.Data(), cent[i].Data()));
      SetStyleTLatex( labelCollisionSystemRatio, 0.035 , 4, 1, 42, kTRUE, 11); 

      
      //**********************************************************************************************************************
      //********************************************* PLOTTING SPECTRA Pi0 ***************************************************
      //**********************************************************************************************************************
      cout << "plotting spectra Pi0..." << endl;
      
      // start drawing plots for cent classes separately
      canvas->Clear();   // does not work without this line
      DrawGammaCanvasSettings( canvas,  0.13, 0.01, 0.015, 0.08);
      canvas->Divide(3,1);
	
      // +++ Pi0 +++ Data  histo & MC Added histo UN-weighted +++
      canvas->cd(1);
      gPad->SetLeftMargin(0.16);
      gPad->SetRightMargin(0.0);
      gPad->SetLogy(1);
      gPad->SetLogx(1);
      histo1DSpectra->DrawCopy();

      histoYieldDataPi0[i]->Draw("p,same,e1");   
      DrawGammaSetMarker(histoYieldDataPi0[i], markerStyleData, markerSizeSpectrum, colorData , colorData);
      legendSpectra->AddEntry(histoYieldDataPi0[i], " Data, corr. with new weighted effi", "pe");
      fitPi0DataYield[i]->Draw("same");
      legendSpectra->AddEntry(fitPi0DataYield[i], Form(" New fit %s to data", fitFunctionsPi0[i].Data()),"l");
      if(iteration == 0 || iteration == 3  || iteration == 4){
	DrawGammaSetMarker(histoPi0InputMCWOWeights[i], markerStyleMC, 1.0, colorMC, colorMC);        // marker style, size, color, line color
	histoPi0InputMCWOWeights[i]->Draw("same");
	legendSpectra->AddEntry(histoPi0InputMCWOWeights[i], " Added MC particles WoPTW ","p");
	DrawGammaSetMarker(histoPi0InputMBMCWOWeights[i], markerStyleMBMC, 1.0, colorMBMC, colorMBMC);      
	histoPi0InputMBMCWOWeights[i]->Draw("same");
	legendSpectra->AddEntry(histoPi0InputMBMCWOWeights[i], " Normal MC particles WoPTW","p");
      } else {
	DrawGammaSetMarker(histoPi0InputMCWWeights[i], markerStyleMC, 1.0, colorMC, colorMC);       
	histoPi0InputMCWWeights[i]->Draw("same");
	legendSpectra->AddEntry(histoPi0InputMCWWeights[i], " Added MC particles WPTW","p");
	DrawGammaSetMarker(histoPi0InputMBMCWWeights[i], markerStyleMBMC, 1.0, colorMBMC, colorMBMC);        	
	histoPi0InputMBMCWWeights[i]->Draw("same");
	legendSpectra->AddEntry(histoPi0InputMBMCWWeights[i], " Normal MC particles WPTW","p");

	fitPi0FromFile[i]->SetLineColor(colorOldFit); 
	fitPi0FromFile[i]->Draw("same");
	legendSpectra->AddEntry(fitPi0FromFile[i], "fit used for weighting","l");
	
      }          
      histo1DSpectra->Draw("p,same,e1");
      labelCollisionSystem->Draw();
      labelSpectraPi0Label->Draw();
      //labelIteration->Draw();
      legendSpectra->Draw();   

    // **********************************************************************************************************************
    // ***************************************** PLOTTING RATIOS PI0 ********************************************************
    // **********************************************************************************************************************
      cout << "plotting ratios Pi0..." << endl;

      // +++ Pi0 +++ Data histo / Data Fit +++
      canvas->cd(2);
      gPad->SetLeftMargin(0.1);
      gPad->SetRightMargin(0.0);
      gPad->SetLogx(1);
      gPad->SetLogy(0);
      histo1DRatio->DrawCopy();
      DrawGammaSetMarker(histoPi0RatioDataToFit[i], 20, 1.0, colorFit, colorFit);        // marker style, size, color, line color
      histoPi0RatioDataToFit[i]->Draw("same");
      if(iteration == 5){
	DrawGammaSetMarker(histoPi0RatioDataToOldFit[i], 24, 1.5, colorOldFit, colorOldFit);   
	histoPi0RatioDataToOldFit[i]->Draw("same");
	legendRatio1->AddEntry(histoPi0RatioDataToFit[i], "to new fit", "p");
	legendRatio1->AddEntry(histoPi0RatioDataToOldFit[i], "to old fit", "p");
	legendRatio1->Draw("same");
      }
      histo1DRatio->DrawCopy("SAME");       
      labelCollisionSystemRatio->Draw();  
      labelSpectraPi0LabelRatio->Draw();
      DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
      
      // +++ Pi0 +++ MC / Data fit +++
      canvas->cd(3);
      gPad->SetLeftMargin(0.1);
      gPad->SetRightMargin(0.0);
      if(iteration == 0 || iteration == 3  || iteration == 4){
	gPad->SetLogx(1);
	gPad->SetLogy(1);
      } else {
	gPad->SetLogx(1);
	gPad->SetLogy(0);
      }

      histo1DRatio2Pi0->DrawCopy();

      DrawGammaSetMarker(histoPi0Ratio[i], markerStyleMC, 1.0, colorMC, colorMC);        // marker style, size, color, line color
      histoPi0Ratio[i]->Draw("same");
      if(iteration == 0 || iteration == 3  || iteration == 4) legendRatio->AddEntry(histoPi0Ratio[i], " Added MC particles WoPTW","p");
      else               legendRatio->AddEntry(histoPi0Ratio[i], " Added MC particles WPTW","p");

      DrawGammaSetMarker(histoPi0MBRatio[i], markerStyleMBMC, 1.0, colorMBMC, colorMBMC); 
      histoPi0MBRatio[i]->Draw("same");
      if (iteration == 0  || iteration == 4) legendRatio->AddEntry(histoPi0MBRatio[i], " Normal MC particles WoPTW","p");
      else                                   legendRatio->AddEntry(histoPi0MBRatio[i], " Normal MC particles WPTW","p");
      
      histo1DRatio2Pi0->DrawCopy("SAME");
      legendRatio->Draw("same");
      labelCollisionSystemRatio->Draw();  
      labelSpectraPi0LabelRatio->Draw();
      DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);


      // +++ save +++
      canvas->Update();
      canvas->SaveAs(Form("%s/Pi0_%s_%s_it%d.%s",outputDir.Data(),collisionSystemForWriting.Data(),cent[i].Data(),iteration,suffix.Data()));

      // cleanup
      canvas->Clear();
      DrawGammaCanvasSettings( canvas,  0.13, 0.01, 0.015, 0.08);
      canvas->Divide(3,1);
      legendSpectra->Clear();
      legendRatio->Clear();
      
      //**********************************************************************************************************************
      //********************************************* PLOTTING SPECTRA ETA ***************************************************
      //**********************************************************************************************************************
      cout << "plotting spectra Eta..." << endl;
      
      // +++ Eta +++ Data  histo & MC Added histo UN-weighted +++
      canvas->cd(1);
      gPad->SetLeftMargin(0.16);
      gPad->SetRightMargin(0.0);
      gPad->SetLogy(1);
      gPad->SetLogx(1);
      histo1DSpectra->DrawCopy(); 
      histoYieldDataEta[i]->Draw("p,same,e1");
      DrawGammaSetMarker(histoYieldDataEta[i], markerStyleData, markerSizeSpectrum, colorData , colorData);
      legendSpectra->AddEntry(histoYieldDataEta[i], " Data, corr. with new weighted effi", "pe");
      fitEtaDataYield[i]->Draw("same");
      legendSpectra->AddEntry(fitEtaDataYield[i], Form(" New fit %s to data", fitFunctionsEta[i].Data()),"l");
      if(iteration == 0 || iteration == 3 || iteration == 4){
	DrawGammaSetMarker(histoEtaInputMCWOWeights[i], markerStyleMC, 1.0, colorMC, colorMC);        // marker style, size, color, line color
	histoEtaInputMCWOWeights[i]->Draw("same");
	legendSpectra->AddEntry(histoEtaInputMCWOWeights[i], " Added MC particles WoPTW ","p");
	DrawGammaSetMarker(histoEtaInputMBMCWOWeights[i], markerStyleMBMC, 1.0, colorMBMC, colorMBMC);        // marker style, size, color, line color
	histoEtaInputMBMCWOWeights[i]->Draw("same");
	legendSpectra->AddEntry(histoEtaInputMBMCWOWeights[i], " Normal MC particles WoPTW","p");
      } else {
	DrawGammaSetMarker(histoEtaInputMCWWeights[i], markerStyleMC, 1.0, colorMC, colorMC);        // marker style, size, color, line color
	histoEtaInputMCWWeights[i]->Draw("same");
	legendSpectra->AddEntry(histoEtaInputMCWWeights[i], " Added MC particles WPTW","p");
	DrawGammaSetMarker(histoEtaInputMBMCWWeights[i], markerStyleMBMC, 1.0, colorMBMC, colorMBMC);        // marker style, size, color, line color
	histoEtaInputMBMCWWeights[i]->Draw("same");
	legendSpectra->AddEntry(histoEtaInputMBMCWWeights[i], " Normal MC particles WPTW","p");

	fitEtaFromFile[i]->SetLineColor(colorOldFit);  
	fitEtaFromFile[i]->Draw("same");
	legendSpectra->AddEntry(fitEtaFromFile[i], "fit used for weighting","l");
	
      }          
      histo1DSpectra->Draw("p,same,e1");
      labelCollisionSystem->Draw();
      labelSpectraEtaLabel->Draw();
      //labelIteration->Draw();
      legendSpectra->Draw();  
  
           
      // **********************************************************************************************************************
      // ***************************************** PLOTTING RATIOS ETA ********************************************************
      // **********************************************************************************************************************
      cout << "plotting ratios Eta..." << endl;

      // +++ Eta +++ Data histo / Data Fit +++
      canvas->cd(2);
      gPad->SetLeftMargin(0.1);
      gPad->SetRightMargin(0.0);
      gPad->SetLogx(1);
      gPad->SetLogy(0);
      histo1DRatio->DrawCopy(); 
      DrawGammaSetMarker(histoEtaRatioDataToFit[i], 20, 1.0, colorFit, colorFit);        // marker style, size, color, line color
      histoEtaRatioDataToFit[i]->Draw("same");

      if(iteration == 5){
	DrawGammaSetMarker(histoEtaRatioDataToOldFit[i], 24, 1.5, colorOldFit, colorOldFit);        // marker style, size, color, line color
	histoEtaRatioDataToOldFit[i]->Draw("same");
	legendRatio1->AddEntry(histoEtaRatioDataToFit[i], "to new fit", "p");
	legendRatio1->AddEntry(histoEtaRatioDataToOldFit[i], "to old fit", "p");
	legendRatio1->Draw("same");
      }

      histo1DRatio->DrawCopy("SAME");
      labelCollisionSystemRatio->Draw();  
      labelSpectraEtaLabelRatio->Draw();
      DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);

      // +++ Eta +++ MC / Data fit +++
      canvas->cd(3);
      gPad->SetLeftMargin(0.1);
      gPad->SetRightMargin(0.0);
      if(iteration == 0 || iteration == 3 || iteration == 4){
	gPad->SetLogx(1);
	gPad->SetLogy(1);
      } else {
	gPad->SetLogx(1);
	gPad->SetLogy(0);
      }
      
      histo1DRatio2Eta->DrawCopy();
      

      DrawGammaSetMarker(histoEtaRatio[i], markerStyleMC, 1.0, colorMC, colorMC);        // marker style, size, color, line color
      histoEtaRatio[i]->Draw("same");
      if(iteration == 0 || iteration == 3 || iteration == 4) legendRatio->AddEntry(histoEtaRatio[i], " Added MC particles WoPTW","p");
      else                                                   legendRatio->AddEntry(histoEtaRatio[i], " Added MC particles WPTW","p");
      
      DrawGammaSetMarker(histoEtaMBRatio[i], markerStyleMBMC, 1.0, colorMBMC, colorMBMC);  
      histoEtaMBRatio[i]->Draw("same");
      if (iteration == 0  || iteration == 4) legendRatio->AddEntry(histoEtaMBRatio[i], " Normal MC particles WoPTW","p");
      else                                   legendRatio->AddEntry(histoEtaMBRatio[i], " Normal MC particles WPTW","p");
      
      histo1DRatio2Eta->DrawCopy("SAME");
      legendRatio->Draw("same");
      labelCollisionSystemRatio->Draw();  
      labelSpectraEtaLabelRatio->Draw();
      DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
      
      // +++ save +++
      canvas->Update();
      canvas->SaveAs(Form("%s/Eta_%s_%s_it%d.%s",outputDir.Data(),collisionSystemForWriting.Data(),cent[i].Data(),iteration,suffix.Data()));

      // cleanup
      canvas->Clear();
      DrawGammaCanvasSettings( canvas,  0.13, 0.01, 0.015, 0.08);
      canvas->Divide(3,1);
      legendSpectra->Clear();
      legendRatio->Clear();

      delete labelCollisionSystem;
      delete labelCollisionSystemRatio;
	  
    } // end of loop over cent classes

    // cleanup

    delete histo1DSpectra;
    delete labelSpectraPi0Label;
    delete labelSpectraEtaLabel;
    delete canvas;
    delete histo1DRatio;
    delete histo1DRatio2Pi0;
    delete histo1DRatio2Eta;
    delete labelSpectraPi0LabelRatio;
    delete labelSpectraEtaLabelRatio;
    cout << "done plotting" << endl;
    
  } // end of producePlots         
} // end of main


/*
// set point 0,0
    if(kFALSE){
      graphYieldDataEta[i]->Print();
      TGraphAsymmErrors *graphYieldDataEtaNew = new TGraphAsymmErrors;
      graphYieldDataEtaNew->SetPoint(1, 0.7, 1E-4);
      graphYieldDataEtaNew->SetPoint(1, 0., 0.);
      graphYieldDataEtaNew->SetPointError(1, 0.5, 0.5, 1E-4, 1E-4);
      Double_t x = -1;
      Double_t y = -1;
      for(int p=1; p<=graphYieldDataEta[i]->GetN(); p++){
	graphYieldDataEta[i]->GetPoint(p, x, y);
	graphYieldDataEtaNew->SetPoint(p+1, x, y);
	cout << "Point " << p << " of graph "<< i << ": " << x << " "<< y << endl;
      }
      graphYieldDataEta[i]=graphYieldDataEtaNew;
      graphYieldDataEta[i]->GetPoint(1, x, y);
      cout << "Point 1 of graph "<< i << ": " << x << " "<< y << endl;
    }*/
