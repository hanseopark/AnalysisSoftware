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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"

void prepareGammaPtWeightFile_ppV2(TString inputPath = "/Users/marin/analysis/gammaPtWeightFile",
                              TString eventAndphotonCutNo = "00010103_0d000009266300008850404000",
                              TString targetFileName = "histosForGammaPtWeighting_ppV2.root",
                              TString energy = ""){

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

    TString centralityClasses1[nCentralityClasses] = {"000"};   // config 23

    // file names for Pythia - Phojet - all  MB
    TString productionNames13TeV[nProductions] = {"LHC16d", "LHC16P1Pyt8","LHC18d6a2" };
    TString fileNames13TeV[nProductions]       = {"AdditionalMBHistosLHC16_Pythia_00010103_0d000009266300008850404000.root",  // LHC16d
                                                  "AdditionalMBHistosLHC16_Pythia_00010103_0d000009266300008850404000.root", // Pythia 13 TeV
                                                  "AdditionalMBHistosLHC16_Phojet_00010103_0d000009266300008850404000.root" };  // Phojet 13 TeV

    TString productionNames5TeV[nProductions] = {"LHC17p", "LHC17l3b","LHC18d6b" };
    TString fileNames5TeV[nProductions]       = {"MaterialBudget_LHC17p_fastwoSDD_LowInt_23.root",  // LHC16d
                                                 "MaterialBudget_PythiaLowInt_fastwoSDD_23.root", // Pythia 13 TeV
                                                 "MaterialBudget_PhojetLowInt_fastwoSDD_23.root" };  // Phojet 13 TeV

    TString productions;
    TString fileNames;
    TString centralityClass;

    const int nBinsPt = 69;
    Double_t fMaxPt = 20.;
    Double_t arrayPtBins[nBinsPt]; 

    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 40) arrayPtBins[i]              = 0.05*i;
      else if(i<50) arrayPtBins[i]           = 2.+ 0.1*(i-40);
      else if(i<64) arrayPtBins[i]          = 3. + 0.5*(i-50);
      else if(i<69) arrayPtBins[i]          = 10.+ 2.0*(i-64);
      else arrayPtBins[i]                    = fMaxPt;
    }
    

    //########## start ########################
    TFile targetFile( Form("%s/%s", inputPath.Data(),targetFileName.Data()),"RECREATE");
    TH1F* histoGammaPtRebin=0x0;
    // go through productions j
    for(Int_t j=0; j<nProductionsUsed; j++){
        if(energy.Contains("5TeV")){
            productions = productionNames5TeV[j];
            fileNames   = fileNames5TeV[j];
        } else if (energy.Contains("13TeV")){
            productions = productionNames13TeV[j];
            fileNames   = fileNames13TeV[j];
        }
        cout << productions.Data() << endl;
        centralityClass   = centralityClasses1[0];

        TString fullCutString = Form("%s", eventAndphotonCutNo.Data());
        TString inputFileName = Form("%s/%s/%s", inputPath.Data(), productions.Data(), fileNames.Data());
        if(energy.Contains("5TeV")) inputFileName = Form("%s/%s", inputPath.Data(), fileNames.Data());
        cout<< inputFileName.Data()<< endl;
        cout<< fullCutString.Data() << endl;

        TFile inputFile(inputFileName, "READ");
 	if (j==0){
	  histoGammaPtRebin= (TH1F*)inputFile.Get("histoPtSumRBinDataSecSubtractedUsingCocktailRebin");
	} else{
	  histoGammaPtRebin= (TH1F*)inputFile.Get("histoPtSumRMCSecSubtractedRebin");
	}
	histoGammaPtRebin->Draw();
	//        topList->SetOwner(kTRUE);  // enables recursive free of memory
	//        cutList->SetOwner(kTRUE);
	//        subList->SetOwner(kTRUE);

	//        TH1D* histoTracksSave = (TH1D*)histoTracks->Clone();
	TH1D* histoGammaPtRebinSave = (TH1D*)histoGammaPtRebin->Clone();

	histoGammaPtRebinSave->Sumw2();
	histoGammaPtRebinSave->SetName(Form("pT_%s_%s", productions.Data(), centralityClass.Data()));
        targetFile.cd();
	histoGammaPtRebinSave->Write();
        cout << "\t done"<< endl;
        // clean up
        inputFile.Close();

    }
    cout << "to " <<  targetFileName.Data() << endl;
    targetFile.Close();
}
