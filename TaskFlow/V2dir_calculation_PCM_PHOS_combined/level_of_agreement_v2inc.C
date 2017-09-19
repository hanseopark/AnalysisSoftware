//
// quantify the level of agreement between v2,inc measurements from
// PCM and PHOS talking into account pT-correlated uncertainties
//
// Approach: check whether difference between v2,inc,PCM and v2,inc PHOS
// is consistent with zero
//
// Klaus Reygers, 18 September 2017
//

Double_t DotProduct(const TVectorD& a, const TVectorD& b);
Double_t chi2(const TVectorD& x, const TVectorD& mu, const TMatrixDSym& V);

void level_of_agreement_v2inc() {

  // TString centr = "00-20";
  TString centr = "20-40";

  TString fn_pcm = "data_pcm/data_pcm.root";
  TFile f_pcm(fn_pcm.Data());
  TDirectory* dir_pcm = (TDirectory*) f_pcm.Get(centr.Data());

  TString fn_phos = "data_phos/data_phos.root";
  TFile f_phos(fn_phos.Data());
  TDirectory* dir_phos = (TDirectory*) f_phos.Get(centr.Data());

  TVectorD v2_inc_meas_values_pcm = *(TVectorD*) dir_pcm->Get("v2_inc_meas_values");
  TMatrixDSym cov_v2_inc_toterr_pcm = *(TMatrixDSym*) dir_pcm->Get("cov_v2_inc_toterr");
  TMatrixDSym cov_v2_inc_staterr_pcm = *(TMatrixDSym*) dir_pcm->Get("cov_v2_inc_staterr");
  TMatrixDSym cov_v2_inc_syserr_pcm = *(TMatrixDSym*) dir_pcm->Get("cov_v2_inc_syserr");

  TVectorD v2_inc_meas_values_phos = *(TVectorD*) dir_phos->Get("v2_inc_meas_values");
  TMatrixDSym cov_v2_inc_toterr_phos = *(TMatrixDSym*) dir_phos->Get("cov_v2_inc_toterr");
  TMatrixDSym cov_v2_inc_syserr_phos = *(TMatrixDSym*) dir_phos->Get("cov_v2_inc_syserr");
  TMatrixDSym cov_v2_inc_staterr_phos = *(TMatrixDSym*) dir_phos->Get("cov_v2_inc_staterr");

  // check: scale factor for uncertainties
  const Double_t sf = 1.2;
  cov_v2_inc_toterr_pcm *= sf*sf;
  cov_v2_inc_toterr_phos *= sf*sf;

  TVectorD v2_inc_diff = v2_inc_meas_values_phos - v2_inc_meas_values_pcm;
  TMatrixDSym cov_v2_inc_diff = cov_v2_inc_toterr_phos + cov_v2_inc_toterr_pcm;

  TVectorD x(16); // contains only zeros
  Double_t chi2val = chi2(x, v2_inc_diff, cov_v2_inc_diff);
  cout << centr << "%: ";
  cout << "chi2 = " << chi2val;
  cout << ", p-value = " << TMath::Prob(chi2val, 16) << endl;

}

Double_t DotProduct(const TVectorD& a, const TVectorD& b) {
  Int_t n = a.GetNrows();
  Double_t dotprod = 0;
  for (Int_t i=0; i<n; i++) dotprod += a(i)*b(i);
  return dotprod;
}

Double_t chi2(const TVectorD& x, const TVectorD& mu, const TMatrixDSym& V) {
  TVectorD d = x - mu;
  TMatrixDSym V_copy = V;
  TMatrixDSym Vinv = V_copy.Invert();
  return DotProduct(d, Vinv*d);
}
