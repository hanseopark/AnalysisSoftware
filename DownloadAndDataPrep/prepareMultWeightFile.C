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
#include <vector>
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

void prepareMultWeightFile(){

  // this macro reads GoodESDTracks histograms from GammaConv files
  // for different centrality classes, for data and different MCs
  // and saves them in a root file which has to be handed over to the addtask to apply multiplicity weighting
  // has to be compiled!
  // current settings are for LHC15o

  //######### settings#########################
  const Int_t nProductions                 = 17;
  //const Int_t nCentralityClasses         = 9;  // nCentClassesInFirstConfig + nCentClassesInSecondConfig + nCentClassesAdditional
  const Int_t nCentClassesInFirstConfig    = 4;
  const Int_t nCentClassesInSecondConfig   = 4;
  const Int_t nCentClassesAdditional       = 1;

  const Int_t nProductionsUsed               = 17;
  const Int_t nCentralityClassesUsed         = 1;
  const Int_t nCentClassesInFirstConfigUsed  = 0;
  const Int_t nCentClassesInSecondConfigUsed = 0;
  const Int_t nCentClassesAdditionalUsed     = 1;

  const TString productionNames[nProductions] = {"LHC18e1", "LHC18e1a", "LHC18e1b", "LHC18e1b_extra", "LHC18e1c", "LHC18e1c_extra",
						 "LHC15o",
						 "LHC16h4",
						 "LHC16i1a", "LHC16i1b", "LHC16i1c",
						 "LHC16i2a", "LHC16i2b", "LHC16i2c", 
						 "LHC16i3a", "LHC16i3b", "LHC16i3c"};
  // first three digits of event cut number:
  const TString centralityClasses1[nCentClassesInFirstConfig]  = {"101", "112", "124", "301"};    // config 294
  const TString centralityClasses2[nCentClassesInSecondConfig] = {"146", "125", "168", "312"};    // config 295
  const TString centralityClasses3[nCentClassesAdditional]     = {"102"};                         // config 295 starting Feb 2019
  
  const TString eventCutNo                             = "10a13";   // rest after first three digits
  // TString eventCutNo                                = "10a23";   // for added particles
  const TString inputPath               = "/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD";
  const TString photonAndMesonCutNo     = "00200009247602008250404000_0152501500000000";
  //TString dateForOutput 	      	= ReturnDateStringForOutput();
  TString dateForOutput 	      	= "2019_03_14";
  const TString targetFileName          = Form("histosForMultWeighting_%s.root",dateForOutput.Data());
  const TString optionFileHandling      = "UPDATE";

  // file names for config 294
  const TString fileNames1[nProductions] = {"GammaConvV1_294_AOD198_highIR_train1197.root",  // LHC18e1
					    "GammaConvV1_294_AOD198_highIR_train1198.root",  // LHC18e1a
					    "GammaConvV1_294_AOD198_highIR_train1199.root",  // LHC18e1b
					    "GammaConvV1_294_AOD198_highIR_train1199.root",  // LHC18e1b_extra
					    "GammaConvV1_294_AOD198_highIR_train1200.root",  // LHC18e1c
					    "GammaConvV1_294_AOD198_highIR_train1200.root",  // LHC18e1c_extra
					    "GammaConvV1_294_AOD194_highIR_train479.root",   // LHC15o 
					    "GammaConvV1_294_AOD198_highIR_train1050.root",  // LHC16h4
					    "GammaConvV1_294_AOD198_highIR_train1051.root",  // LHC16i1a
					    "GammaConvV1_294_AOD198_highIR_train1052.root",  // LHC16i1b
					    "GammaConvV1_294_AOD198_highIR_train1053.root",  // LHC16i1c
					    "GammaConvV1_294_AOD198_highIR_train1054.root",  // LHC16i2a
					    "GammaConvV1_294_AOD198_highIR_train1055.root",  // LHC16i2b
					    "GammaConvV1_294_AOD198_highIR_train1056.root",  // LHC16i2c
					    "GammaConvV1_294_AOD198_highIR_train1057.root",  // LHC16i3a
					    "GammaConvV1_294_AOD198_highIR_train1058.root",  // LHC16i3b
					    "GammaConvV1_294_AOD198_highIR_train1059.root"}; // LHC16i3c 

  // file names for config 295
  const TString fileNames2[nProductions] = {"GammaConvV1_295_AOD198_highIR_train1197.root",  // LHC18e1
					    "GammaConvV1_295_AOD198_highIR_train1198.root",  // LHC18e1a
					    "GammaConvV1_295_AOD198_highIR_train1199.root",  // LHC18e1b
					    "GammaConvV1_295_AOD198_highIR_train1199.root",  // LHC18e1b_extra
					    "GammaConvV1_295_AOD198_highIR_train1200.root",  // LHC18e1c
					    "GammaConvV1_295_AOD198_highIR_train1200.root",  // LHC18e1c_extra
					    "GammaConvV1_295_AOD194_highIR_train479.root",   // LHC15o 
					    "GammaConvV1_295_AOD198_highIR_train1050.root",  // LHC16h4
					    "GammaConvV1_295_AOD198_highIR_train1051.root",  // LHC16i1a
					    "GammaConvV1_295_AOD198_highIR_train1052.root",  // LHC16i1b
					    "GammaConvV1_295_AOD198_highIR_train1053.root",  // LHC16i1c
					    "GammaConvV1_295_AOD198_highIR_train1054.root",  // LHC16i2a
					    "GammaConvV1_295_AOD198_highIR_train1055.root",  // LHC16i2b
					    "GammaConvV1_295_AOD198_highIR_train1056.root",  // LHC16i2c
					    "GammaConvV1_295_AOD198_highIR_train1057.root",  // LHC16i3a
					    "GammaConvV1_295_AOD198_highIR_train1058.root",  // LHC16i3b
					    "GammaConvV1_295_AOD198_highIR_train1059.root"}; // LHC16i3c 

  // file names for 0-20% without mult weights
  const TString fileNames3[nProductions] = {"GammaConvV1_295_AOD198_highIR_train1218.root",  // LHC18e1
					    "GammaConvV1_295_AOD198_highIR_train1214.root",  // LHC18e1a
					    "GammaConvV1_295_AOD198_highIR_train1215.root",  // LHC18e1b
					    "GammaConvV1_295_AOD198_highIR_train1215.root",  // LHC18e1b_extra
					    "GammaConvV1_295_AOD198_highIR_train1216.root",  // LHC18e1c
					    "GammaConvV1_295_AOD198_highIR_train1216.root",  // LHC18e1c_extra
					    "GammaConvV1_295_AOD194_highIR_train501.root",   // LHC15o 
					    "GammaConvV1_295_AOD198_highIR_train1217.root",  // LHC16h4
					    "GammaConvV1_295_AOD198_highIR_train1221.root",  // LHC16i1a
					    "GammaConvV1_295_AOD198_highIR_train1222.root",  // LHC16i1b
					    "GammaConvV1_295_AOD198_highIR_train1223.root",  // LHC16i1c
					    "GammaConvV1_295_AOD198_highIR_train1224.root",  // LHC16i2a
					    "GammaConvV1_295_AOD198_highIR_train1225.root",  // LHC16i2b
					    "GammaConvV1_295_AOD198_highIR_train1226.root",  // LHC16i2c
					    "GammaConvV1_295_AOD198_highIR_train1227.root",  // LHC16i3a
					    "GammaConvV1_295_AOD198_highIR_train1228.root",  // LHC16i3b
					    "GammaConvV1_295_AOD198_highIR_train1229.root"}; // LHC16i3c 

  
    //########## start ########################
    TFile targetFile(targetFileName,optionFileHandling);

    // go through productions j
    for(Int_t j=0; j<nProductionsUsed; j++){
      cout << productionNames[j].Data() << endl;
      // go through centrality classes n
      for(Int_t n=0; n<nCentralityClassesUsed; n++){ 
	Int_t i;  // centrality class array position counter
	TString fileNames;
	TString centralityClass;
	if(n<nCentClassesInFirstConfigUsed) {
	  i = n;
	  //cout << i << endl;
	  fileNames       = fileNames1[j];
	  centralityClass = centralityClasses1[i];
	} else if(n<(nCentClassesInFirstConfigUsed+nCentClassesInSecondConfigUsed)) {
	  i = n-nCentClassesInFirstConfigUsed;
	  //cout << i << endl;
	  fileNames       = fileNames2[j];
	  centralityClass = centralityClasses2[i];
	} else {
	  i = n-nCentClassesInFirstConfigUsed-nCentClassesInSecondConfigUsed;
	  //cout << i << endl;
	  fileNames       = fileNames3[j];
	  centralityClass = centralityClasses3[i];
	}

	TString fullEventCutNo = Form("%s%s", centralityClass.Data(), eventCutNo.Data());
	TString fullCutString  = Form("%s_%s", fullEventCutNo.Data(), photonAndMesonCutNo.Data());
	TString centForTitle   = GetCentralityString(fullEventCutNo);
	cout << "\t" << centForTitle << endl;
	TString inputFileName = Form("%s/%s/%s", inputPath.Data(), productionNames[j].Data(), fileNames.Data());
	//gSystem->Exec(Form("ls %s", inputFileName.Data()));
	TFile inputFile(inputFileName, "READ");
	TList *topList        = (TList*)inputFile.Get("GammaConvV1");
	TList *cutList        = (TList*)topList->FindObject(Form("Cut Number %s",fullCutString.Data()));
	TList *subList        = (TList*)cutList->FindObject(Form("%s ESD histograms", fullCutString.Data()));
	TH1D *histoTracks     = (TH1D*)subList->FindObject("GoodESDTracks");
	topList->SetOwner(kTRUE);  // enables recursive free of memory
	cutList->SetOwner(kTRUE);
	subList->SetOwner(kTRUE);
	TH1D* histoTracksSave = (TH1D*)histoTracks->Clone();
	histoTracksSave->SetTitle(Form("GoodESDTracks %s %s", productionNames[j].Data(), centForTitle.Data()));
	histoTracksSave->SetName(Form("%s_%s", productionNames[j].Data(), centralityClass.Data()));
	targetFile.cd();
	histoTracksSave->Write();
	cout << "\t done"<< endl;
	// clean up
	inputFile.Close();
	delete topList;   // needs to be deleted because Close() does not free list memory correctly
      } // end of loop over cent classes
    } // end of loop over productions
    cout << "to " <<  targetFileName.Data() << endl;
    targetFile.Close(); 
}
