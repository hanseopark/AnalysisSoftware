#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "iostream"

void CreateNChargeWeightingFile(){

// 2.76 TeV
//	const Int_t nFiles=15;

//	TString mainDir[nFiles];
//	TString path[nFiles]={
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_LHC11a-pass4_1.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_LHC11a-pass4_1.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_LHC13g-pass3_96.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_LHC13g-pass3_96.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_LHC13g-pass3_95.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_LHC13g-pass3_95.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_MC_LHC12f1a_1.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_MC_LHC12f1a_1.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_MC_LHC12f1b_1.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_MC_LHC12f1b_1.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20150825-ConvCalo/GammaConvCalo_MC_LHC12i3_1.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_MC_LHC15g2_96.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_MC_LHC15g2_96.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_MC_LHC15g2_95.root",
//		"/home/daniel/data/work/Grid/Legotrain-vAN-20151019-2.76TeV-QA/GammaConvCalo_MC_LHC15g2_95.root"
//	};
//	TString name[nFiles]={
//		"LHC11a_00",
//		"LHC11a_51",
//		"LHC13g_00",
//		"LHC13g_52",
//		"LHC13g_83",
//		"LHC13g_85",
//		"LHC12f1a_00",
//		"LHC12f1a_51",
//		"LHC12f1b_00",
//		"LHC12f1b_51",
//		"LHC12i3_00",
//		"LHC15g2_00",
//		"LHC15g2_52",
//		"LHC15g2_83",
//		"LHC15g2_85"
//	};
//	Int_t cutNr[nFiles]={0,1,0,1,0,1,0,1,0,1,0,0,1,0,1};

  // 8TeV
  const Int_t nFiles=63;

  TString mainDir[nFiles];
  TString path[nFiles]={
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12a_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12b_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12c_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12d_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12f_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12h_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12i_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12a_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12b-kEMC7_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12c-kEMC7_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12d-kEMC7_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12f-kEMC7_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12h-kEMC7_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12i-kEMC7_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12a_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC12b_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12c-kEMCEGA_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12d-kEMCEGA_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12f-kEMCEGA_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12h-kEMCEGA_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA-Trigger/LHC12i-kEMCEGA_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1a1_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1b_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1c_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1d_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1f_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1h_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1i_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2a_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2b_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2c_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2d_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2f_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2h_GammaConvCalo_101.root",
        "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2i_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1a1_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1b_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1c_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1d_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1f_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1h_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1i_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2a_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2b_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2c_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2d_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2f_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2h_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2i_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1a1_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1b_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1c_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1d_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1f_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1h_GammaConvCalo_101.root",
    "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h1i_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2a_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2b_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2c_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2d_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2f_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2h_GammaConvCalo_101.root",
      "/home/daniel/data/work/Grid/Legotrain-vAN-20160518-8TeV-QA/LHC15h2i_GammaConvCalo_101.root",
  };
  TString name[nFiles]={
      "LHC12a_10",
      "LHC12b_10",
      "LHC12c_10",
      "LHC12d_10",
      "LHC12f_10",
      "LHC12h_10",
      "LHC12i_10",
    "LHC12a_52",
    "LHC12b_52",
    "LHC12c_52",
    "LHC12d_52",
    "LHC12f_52",
    "LHC12h_52",
    "LHC12i_52",
    "LHC12a_81",
    "LHC12b_81",
    "LHC12c_81",
    "LHC12d_81",
    "LHC12f_81",
    "LHC12h_81",
    "LHC12i_81",
      "LHC15h1a1_10",
      "LHC15h1b_10",
      "LHC15h1c_10",
      "LHC15h1d_10",
      "LHC15h1f_10",
      "LHC15h1h_10",
      "LHC15h1i_10",
      "LHC15h2a_10",
      "LHC15h2b_10",
      "LHC15h2c_10",
      "LHC15h2d_10",
      "LHC15h2f_10",
      "LHC15h2h_10",
      "LHC15h2i_10",
    "LHC15h1a1_52",
    "LHC15h1b_52",
    "LHC15h1c_52",
    "LHC15h1d_52",
    "LHC15h1f_52",
    "LHC15h1h_52",
    "LHC15h1i_52",
    "LHC15h2a_52",
    "LHC15h2b_52",
    "LHC15h2c_52",
    "LHC15h2d_52",
    "LHC15h2f_52",
    "LHC15h2h_52",
    "LHC15h2i_52",
    "LHC15h1a1_81",
    "LHC15h1b_81",
    "LHC15h1c_81",
    "LHC15h1d_81",
    "LHC15h1f_81",
    "LHC15h1h_81",
    "LHC15h1i_81",
    "LHC15h2a_81",
    "LHC15h2b_81",
    "LHC15h2c_81",
    "LHC15h2d_81",
    "LHC15h2f_81",
    "LHC15h2h_81",
    "LHC15h2i_81"
  };
  Int_t cutNr[nFiles]={0,0,0,0,0,0,0,
                       0,1,1,1,1,1,1,
                       0,0,2,2,2,2,2,
                       0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,
                       1,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,
                       2,2,2,2,2,2,2,
                       2,2,2,2,2,2,2
                      };

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
