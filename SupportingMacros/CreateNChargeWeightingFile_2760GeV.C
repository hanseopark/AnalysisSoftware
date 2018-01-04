// Hikari Murakmi, Pedro Gonzalez
// for Efficiency study for pp 2.76TeV 
// To calculate Charged track multiplicity weighing for JetJet

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "iostream"

void CreateNChargeWeightingFile_2760GeV(){

  // 2.76 TeV
  const Int_t nFiles=10;

  TString mainDir[nFiles];
  TString path[nFiles]={
    //    "~/cernbox/Train/GA_pp_MC/2569_20161024_JJ/GammaConvV1_MC_LHC15g1a_JJ_2760GeV_ancLHC11a.root"

    //********This is for 0th iteration, created input file ESDTracks_ChargedWeight_PCM_temp.root **********************
    // "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161106-ConvTesting/GammaConvV1_LHC11a-pass4-WSDD_7.root",
    // "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161106-ConvTesting/GammaConvV1_LHC11a-pass4-WSDD_7.root",
    // "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161106-ConvTesting/GammaConvV1_MC_LHC12f1a_7.root",
    // "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161106-ConvTesting/GammaConvV1_MC_LHC12f1a_7.root",
    // "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161106-ConvTesting/GammaConvV1_MC_LHC12f1b_7.root",
    // "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161106-ConvTesting/GammaConvV1_MC_LHC12f1b_7.root",
    // "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161106-ConvTesting/GammaConvV1_MC_LHC12i3-WSDD_7.root",
    // "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161106-ConvTesting/GammaConvV1_MC_LHC12i3-WSDD_7.root",
    // "~/cernbox/Train/GA_pp_MC/2611_20161103_JJ/GammaConvV1_LHC15g1a_JJ_WSDD_7.root",
    // "~/cernbox/Train/GA_pp_MC/2611_20161103_JJ/GammaConvV1_LHC15g1a_JJ_WSDD_7.root"
    //********This is for 0th iteration, created input file ESDTracks_ChargedWeight_PCM_temp.root **********************
    "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161115-ConvTesting/GammaConvV1_LHC11a-pass4-WSDD_7.root",
    "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161115-ConvTesting/GammaConvV1_LHC11a-pass4-WSDD_7.root",
    "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161115-ConvTesting/GammaConvV1_MC_LHC12f1a_139.root",
    "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161115-ConvTesting/GammaConvV1_MC_LHC12f1a_139.root",
    "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161115-ConvTesting/GammaConvV1_MC_LHC12f1b_139.root",
    "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161115-ConvTesting/GammaConvV1_MC_LHC12f1b_139.root",
    "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161115-ConvTesting/GammaConvV1_MC_LHC12i3-WSDD_139.root",
    "~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161115-ConvTesting/GammaConvV1_MC_LHC12i3-WSDD_139.root",
    "~/cernbox/Train/GA_pp_MC/2640_20161110-1213_JJ/GammaConvV1_LHC15g1a_JJ_WSDD_139.root",
    "~/cernbox/Train/GA_pp_MC/2640_20161110-1213_JJ/GammaConvV1_LHC15g1a_JJ_WSDD_139.root"

  };
  TString name[nFiles]={
    "LHC11a_00",
    "LHC11a_13",// "_13"...temporary copy 2016/11/09 
    "LHC12f1a_00",
    "LHC12f1a_13",
    "LHC12f1b_00",
    "LHC12f1b_13",
    "LHC12i3_00",
    "LHC12i3_13",
    "LHC15g1a_00",
    "LHC15g1a_13"
  };
  Int_t cutNr[nFiles]={0,0,0,0,0,0,0,0,0,0};

  // // 5.02 TeV
  // const Int_t nFiles=2;

  // TString mainDir[nFiles];
  // TString path[nFiles]={
  //   //    "~/cernbox/Train/GA_pp_MC/2569_20161024_JJ/GammaConvV1_MC_LHC15g1a_JJ_2760GeV_ancLHC11a.root"
  //   "~/cernbox/Train/GA_pp/1871_20161006-1806_15n_p2/GammaConvV1_LHC15n-pass2_28.root",
  //   "~/cernbox/Train/GA_pp_MC/2519_20160913-0927_JJ_ancLHC15n/244628/GammaConvV1_MC_LHC16h3b_JJ_5020G_ancLHC15n_28.root",
  // };
  // TString name[nFiles]={
  //   "LHC15n_00",
  //   "LHC16h3b_00"
  // };
  // Int_t cutNr[nFiles]={0,0};

  TDirectory::AddDirectory(0);
  TFile::AddDirectory(0);

  for(Int_t i=0; i<nFiles; i++){
    TFile* fFile = new TFile(path[i].Data(),"READ");
    if(fFile->IsZombie()){cout << "ERROR: File " << path[i].Data() << " could not be openend! Returning..." << endl; return;}
    else{
      cout << "Reading file: " << path[i].Data();
      TKey *key;
      TIter next(fFile->GetListOfKeys());
      while ((key=(TKey*)next())) {cout << Form(" - found TopDir: %s",key->GetName()); mainDir[i] = key->GetName();}
      cout << endl;
    }

    TList *listInput =(TList*)fFile->Get(mainDir[i].Data());
    listInput->SetOwner(kTRUE);
    TList *listCuts = (TList*)listInput->At(cutNr[i]);
    listCuts->SetOwner(kTRUE);
    TString nameCuts = listCuts->GetName();
    nameCuts.Replace(0,11,"");
    delete listInput;

    TList* TopDir = (TList*) fFile->Get(mainDir[i].Data());
    TopDir->SetOwner(kTRUE);
    if(TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}
    TList* TopContainer= (TList*) TopDir->FindObject(Form("Cut Number %s",nameCuts.Data()));
    TopContainer->SetOwner(kTRUE);
    if(TopContainer == NULL) {cout << "ERROR: " << Form("Cut Number %s",nameCuts.Data()) << " not found in File" << endl; return;}
    TList* ESDContainer = (TList*) TopContainer->FindObject(Form("%s ESD histograms",nameCuts.Data()));
    ESDContainer->SetOwner(kTRUE);
    if(ESDContainer == NULL) {cout << "ERROR: " << Form("%s ESD histograms",nameCuts.Data()) << " not found in File" << endl; return;}

    TH1D* fHistNGoodTracks = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
    fHistNGoodTracks->Sumw2();
    Double_t integral = fHistNGoodTracks->Integral();
    fHistNGoodTracks->Scale(1./integral);

    TFile* fOutput;
    if(i==0) fOutput = new TFile("ESDTracks.root","RECREATE");
    else fOutput = new TFile("ESDTracks.root","UPDATE");

    fHistNGoodTracks->Write(name[i].Data());

    fOutput->Close();
    delete fOutput;

    delete TopDir;
    fFile->Close();
    delete fFile;
  }

  return;
}
