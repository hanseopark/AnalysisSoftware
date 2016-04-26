#include "TString.h"

void create_epos_root_file() {

  TFile* fEpos = new TFile("Epos_pi0_PbPb_2760GeV.root","recreate");

  const Int_t nFiles = 6;

  TString label[nFiles] = {"Epos_pi0_PbPb_2760GeV_00_to_05", "Epos_pi0_PbPb_2760GeV_05_to_10", "Epos_pi0_PbPb_2760GeV_10_to_20", 
			       "Epos_pi0_PbPb_2760GeV_20_to_40", "Epos_pi0_PbPb_2760GeV_40_to_60", "Epos_pi0_PbPb_2760GeV_60_to_80"};

  TGraphErrors* gEpos[nFiles];
  
  for (int iFile=0; iFile<nFiles; iFile++) {

    const Int_t nPoints = 70;
    Double_t pt[nPoints], dndptpy[nPoints], dndptpy_staterr[nPoints];

    TString filename = label[iFile] + ".txt";
    cout << filename << endl;

    ifstream fEposTxt(filename);
       
    Int_t iPoint = 0;
       
    while(!fEposTxt.eof() && iPoint < nPoints){ 
      fEposTxt >> pt[iPoint] >> dndptpy[iPoint] >> dndptpy_staterr[iPoint];
      // cout << pt[iPoint] << "   " << dndptpy[iPoint] << "   " << dndptpy_staterr[iPoint] << endl;

      if (fEposTxt.fail()) {
	fEposTxt.clear();
	fEposTxt.ignore(0xFFFFFF,'\n');
	continue;
      }

      iPoint++;
    }

    fEposTxt.close();
  
    gEpos[iFile] = new TGraphErrors(nPoints, pt, dndptpy, 0, dndptpy_staterr);
    gEpos[iFile]->Print();
    gEpos[iFile]->Write("g" + label[iFile]);
  }

  fEpos->Close();

}
