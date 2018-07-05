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

  //######### settings#################
  const Int_t nProductions       = 17;
  const Int_t nProductionsUsed   = 17;
  const Int_t nCentralityClasses = 10;

  const TString productionNames[nProductions] = {"LHC15o", "LHC16g1", 
						 "LHC16g1a", "LHC16g1b", "LHC16g1b_extra", "LHC16g1c", "LHC16g1c_extra",
						 "LHC16h4",
						 "LHC16i1a", "LHC16i1b", "LHC16i1c",
						 "LHC16i2a", "LHC16i2b", "LHC16i2c", 
						 "LHC16i3a", "LHC16i3b", "LHC16i3c"};
  // first three digits of event cut number:
  const TString centralityClasses1[nCentralityClasses] = {"301","312","112","123","134"};      // config 283
  const TString centralityClasses2[nCentralityClasses] = {"146","168","189", "101", "109"};    // config 285
  const TString eventCutNo                             = "10a13";   // rest after first three digits
  // TString eventCutNo                                = "10a23"; // for added particles
  const TString inputPath               = "/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD";
  const TString photonAndMesonCutNo     = "00200009247602008250404000_0652501500000000";
  const TString targetFileName          = "histosForMultWeighting.root";

  // file names for centrality classes 0-5-10-20-30-40:
  const TString fileNames1[nProductions] = {"GammaConvV1_283_AOD194_highIR_train461_noMissingTracks.root",  // LHC15o 
					    "GammaConvV1_283_AOD198_highIR_train942.root",  // LHC16g1
					    "GammaConvV1_283_AOD198_highIR_train943.root",  // LHC16g1a
					    "GammaConvV1_283_AOD198_highIR_train944.root",  // LHC16g1b
					    "GammaConvV1_283_AOD_highIR_train945.root",     // LHC16g1b_extra
					    "GammaConvV1_283_AOD198_highIR_train946.root",  // LHC16g1c
					    "GammaConvV1_283_AOD_highIR_train947.root",     // LHC16g1c_extra
					    "GammaConvV1_283_AOD198_highIR_train948.root",  // LHC16h4
					    "GammaConvV1_283_AOD198_highIR_train949.root",  // LHC16i1a
					    "GammaConvV1_283_AOD198_highIR_train950.root",  // LHC16i1b
					    "GammaConvV1_283_AOD198_highIR_train951.root",  // LHC16i1c
					    "GammaConvV1_283_AOD198_highIR_train952.root",  // LHC16i2a
					    "GammaConvV1_283_AOD198_highIR_train953.root",  // LHC16i2b
					    "GammaConvV1_283_AOD198_highIR_train954.root",  // LHC16i2c
					    "GammaConvV1_283_AOD198_highIR_train955.root",  // LHC16i3a
					    "GammaConvV1_283_AOD198_highIR_train956.root",  // LHC16i3b
					    "GammaConvV1_283_AOD198_highIR_train957.root"}; // LHC16i3c 

  // file names for centrality classes 40-60-80-90:
  const TString fileNames2[nProductions] = {"GammaConvV1_285_AOD194_highIR_train461_noMissingTracks.root",  // LHC15o
					    "GammaConvV1_285_AOD198_highIR_train942.root",  // LHC16g1
					    "GammaConvV1_285_AOD198_highIR_train943.root",  // LHC16g1a
					    "GammaConvV1_285_AOD198_highIR_train944.root",  // LHC16g1b
					    "GammaConvV1_285_AOD_highIR_train945.root",     // LHC16g1b_extra
					    "GammaConvV1_285_AOD198_highIR_train946.root",  // LHC16g1c
					    "GammaConvV1_285_AOD_highIR_train947.root",     // LHC16g1c_extra
					    "GammaConvV1_285_AOD198_highIR_train948.root",  // LHC16h4
					    "GammaConvV1_285_AOD198_highIR_train949.root",  // LHC16i1a
					    "GammaConvV1_285_AOD198_highIR_train950.root",  // LHC16i1b
					    "GammaConvV1_285_AOD198_highIR_train951.root",  // LHC16i1c
					    "GammaConvV1_285_AOD198_highIR_train952.root",  // LHC16i2a
					    "GammaConvV1_285_AOD198_highIR_train953.root",  // LHC16i2b
					    "GammaConvV1_285_AOD198_highIR_train954.root",  // LHC16i2c
					    "GammaConvV1_285_AOD198_highIR_train955.root",  // LHC16i3a
					    "GammaConvV1_285_AOD198_highIR_train956.root",  // LHC16i3b
					    "GammaConvV1_285_AOD198_highIR_train957.root"}; // LHC16i3c 

    //########## start ########################
    TFile targetFile(targetFileName,"RECREATE");

    // go through productions j
    for(Int_t j=0; j<nProductionsUsed; j++){
      cout << productionNames[j].Data() << endl;
      // go through centrality classes n
      for(Int_t n=0; n<nCentralityClasses; n++){ 
	Int_t i;  // centrality class array position counter
	TString fileNames;
	TString centralityClass;
	if(n<5) {
	  i = n;
	  fileNames = fileNames1[j];
	  centralityClass   = centralityClasses1[i]; 
	} else {
	  i = n-5;
	  fileNames = fileNames2[j];
	  centralityClass   = centralityClasses2[i];
	}
	TString fullEventCutNo    = Form("%s%s", centralityClass.Data(), eventCutNo.Data());
	TString fullCutString = Form("%s_%s", fullEventCutNo.Data(), photonAndMesonCutNo.Data());
	TString centForTitle = GetCentralityString(fullEventCutNo);
	cout << "\t" << centForTitle;
	TString inputFileName = Form("%s/%s/%s", inputPath.Data(), productionNames[j].Data(), fileNames.Data());
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
      }
    }
    cout << "to " <<  targetFileName.Data() << endl;
    targetFile.Close(); 
}
