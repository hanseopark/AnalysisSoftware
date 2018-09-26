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
void prepareMultWeightFile_pp(){

  // Macro copied  from prepareMultWeightFile() that is for PbPb. Meike Danisch
  // Adapted for pp. Ana Marin
  // this macro reads GoodESDTracks histograms from GammaConv files
  // for different centrality classes, for data and different MCs
  // and saves them in a root file which has to be handed over to the addtask to apply multiplicity weighting
  // has to be compiled!
  // current settings are for LHC15o
  // cureent settings are for LHC16d


  //######### settings#################
  const Int_t nProductions       = 3;
  const Int_t nProductionsUsed   = 3;
  const Int_t nCentralityClasses = 1;

  const TString productionNames[nProductions] = {"LHC16d", "LHC16P1Pyt8","LHC18d6a2" };

  // first three digits of event cut number:
  const TString centralityClasses1[nCentralityClasses] = {"000"};   // config 23
  //  const TString centralityClasses2[nCentralityClasses] = {"000","000","000","000","000"};   // config 23


  const TString eventCutNo                             = "10a13";   // rest after first three digits

  const TString inputPath               = "/Users/marin/analysis/multWeightFile";
  const TString eventAndphotonCutNo             = "00010103_0d000009266300008850404000";
  const TString targetFileName          = "histosForMultWeighting_pp.root";


  // file names for Pythia - Phojet- all  MB
  const TString fileNames1[nProductions] = {"GammaConv_Material_LHC16d_23.root",  // LHC16d
					    "GammaConv_Material_MC_LHC17f6total_23.root", // Pythia 13 TeV
					    "GammaConv_Material_MC_LHC18d6a2_23.root" };  // Phojet 13 TeV 

  TString fileNames;
  TString centralityClass;

  //########## start ########################
  TFile targetFile( Form("%s/%s", inputPath.Data(),targetFileName.Data()),"RECREATE");
  
    // go through productions j
    for(Int_t j=0; j<nProductionsUsed; j++){
      cout << productionNames[j].Data() << endl;
      fileNames = fileNames1[j];
      centralityClass   = centralityClasses1[0]; 
      TString fullCutString = Form("%s", eventAndphotonCutNo.Data());
      TString inputFileName = Form("%s/%s/%s", inputPath.Data(), productionNames[j].Data(), fileNames.Data());
      cout<< inputFileName.Data()<< endl;
      cout<< fullCutString.Data() << endl;
      TFile inputFile(inputFileName, "READ");
      //	TList *topList        = (TList*)inputFile.Get("GammaConvV1");
      TList *topList        = (TList*)inputFile.Get("GammaConvMaterial");
      TList *cutList        = (TList*)topList->FindObject(Form("Cut Number %s",fullCutString.Data()));
      TList *subList        = (TList*)cutList->FindObject(Form("%s ESD histograms", fullCutString.Data()));
      TH1D *histoTracks     = (TH1D*)subList->FindObject("GoodESDTracksEta08");
      topList->SetOwner(kTRUE);  // enables recursive free of memory
      cutList->SetOwner(kTRUE);
      subList->SetOwner(kTRUE);
      TH1D* histoTracksSave = (TH1D*)histoTracks->Clone();
      histoTracksSave->Sumw2();
      cout<< histoTracksSave->GetEntries()<< endl;
      if(histoTracksSave->GetEntries()>0) histoTracksSave->Scale(1./histoTracksSave->GetEntries());
      //	histoTracksSave->SetTitle(Form("GoodESDTracks %s %s", productionNames[j].Data(), centForTitle.Data()));
      histoTracksSave->SetTitle(Form("GoodESDTracks08 %s", productionNames[j].Data() ));
      histoTracksSave->SetName(Form("%s_%s", productionNames[j].Data(), centralityClass.Data()));
      targetFile.cd();
      histoTracksSave->Write();
      cout << "\t done"<< endl;
      // clean up
      inputFile.Close();
      delete topList;   // needs to be deleted because Close() does not free list memory correctly

    }
    cout << "to " <<  targetFileName.Data() << endl;
    targetFile.Close(); 
}
