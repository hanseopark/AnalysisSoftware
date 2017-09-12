void purity_ratio() {

  TString fn020 = "Purity_InclusivePhotonSample_020_50200013_00200009007000008250400000_Data LHC10h.root";
  TString fn2040 = "Purity_InclusivePhotonSample_2040_52400013_00200009007000008250400000_Data LHC10h.root";
  TFile* f020 = new TFile(fn020);
  TFile* f2040 = new TFile(fn2040);

  TGraphErrors* purity_data_driven_020 = (TGraphErrors*) f020->Get("Gamma_Purity_020");
  TGraphErrors* purity_MC_driven_020 = (TGraphErrors*) f020->Get("Gamma_Purity_MCdriven_020");

  TGraphErrors* purity_data_driven_2040 = (TGraphErrors*) f2040->Get("Gamma_Purity_2040");
  TGraphErrors* purity_MC_driven_2040 = (TGraphErrors*) f2040->Get("Gamma_Purity_MCdriven_2040");

  cout << "0-20%: ";
  for (Int_t i=0; i<16; i++) {
    Double_t purity_data = purity_data_driven_020->GetY()[i];
    Double_t purity_MC = purity_MC_driven_020->GetY()[i];
    Double_t ratio = purity_data / purity_MC;
    cout << ratio << ", ";
  }
  cout << endl << endl;

  cout << "20-40%: ";
  for (Int_t i=0; i<16; i++) {
    Double_t purity_data = purity_data_driven_2040->GetY()[i];
    Double_t purity_MC = purity_MC_driven_2040->GetY()[i];
    Double_t ratio = purity_data / purity_MC;
    cout << ratio << ", ";
  }

}
