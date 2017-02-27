void read_phos() {

  // file with PHOS results
  TFile file_phos("data_phos/PHOSFlow_v2_final_161116.root");

  // open output file
  TFile f_out("data_phos/data_phos.root","recreate");

  const Int_t n_pt_bins = 16;

  // loop over centralities
  // centrality, 0-20% = cen6, 20-40% = cen3, 40-80% = cen8
  const Int_t n_centr = 3;
  TString centr[n_centr] = {"00-20", "20-40", "40-80"};
  TString centr_cocktail[n_centr] = {"cen6", "cen3", "cen8"}; 
  TString centr_Rgamma[n_centr] = {"0-20%", "20-40%", "40-80%"};

  // inclusive photon vn
  TVectorD v2_inc_meas_values(n_pt_bins); // data
  TMatrixDSym cov_v2_inc_staterr(n_pt_bins); // covariance matrix
  TMatrixDSym cov_v2_inc_syserr(n_pt_bins); // covariance matrix
  TMatrixDSym cov_v2_inc_toterr(n_pt_bins); // covariance matrix
  
  // decay photon vn
  TVectorD v2_dec_meas_values_staterr(n_pt_bins); // data
  TMatrixDSym cov_v2_dec_staterr(n_pt_bins); // covariance matrix
  TVectorD v2_dec_meas_values_syserr(n_pt_bins); // data
  TMatrixDSym cov_v2_dec_syserr(n_pt_bins); // covariance matrix
  
  // Rgamma
  TVectorD Rgamma_meas_vec(n_pt_bins);
  TMatrixDSym cov_Rgamma(n_pt_bins);

  for (Int_t i_centr=0; i_centr<n_centr; i_centr++) {

    TString dir = "PHOS_v2SP_PbPb_2760_Centrality_" + centr[i_centr] + "/";
    
    // get histograms
    TString name_histo_v2inc_staterr = dir + "hPHOS_v2SP_gammaIncl_noK0s_PbPb_cen" + centr[i_centr] + "_Stat";
    TString name_histo_v2inc_syserr = dir + "hPHOS_v2SP_gammaIncl_noK0s_PbPb_cen" + centr[i_centr] + "_Syst";
    TString name_histo_v2dec = dir + "hPHOS_v2SP_Cocktail_noK0s_PbPb_cen" + centr[i_centr] + "_Stat";
    TString name_graph_v2dec_syserr = dir + "hPHOS_v2SP_Cocktail_noK0s_PbPb_cen" + centr[i_centr] + "_Syst";
    TString name_histo_Rgamma = dir + "hPHOS_DoubleRatio_PbPb_cen" + centr[i_centr] + "_Stat";
	    
    TH1D* h_v2inc_staterr = (TH1D*) file_phos.Get(name_histo_v2inc_staterr);
    TH1D* h_v2inc_syserr = (TH1D*) file_phos.Get(name_histo_v2inc_syserr);

    //    for (Int_t ibin=1; ibin<=h_v2inc_toterr->GetNbinsX(); ibin++) {
    //   Double_t sys_err =  h_v2inc_syserr->GetBinError(ibin);
    //   Double_t stat_err =  h_v2inc_staterr->GetBinError(ibin);
    //   Double_t 
    // };

    
    TH1D* h_v2dec = (TH1D*) file_phos.Get(name_histo_v2dec);
    TGraphAsymmErrors* g_v2dec_syserr = (TGraphAsymmErrors*) file_phos.Get(name_graph_v2dec_syserr);
    TH1D* h_Rgamma = (TH1D*) file_phos.Get(name_histo_Rgamma);

    if (!h_v2dec) cout << "Error" << endl;
    
    Double_t pt_bin_center_first_bin = 1.; // GeV
    Int_t ibin_first_v2inc = h_v2inc_staterr->GetXaxis()->FindBin(pt_bin_center_first_bin);
    Int_t ibin_first_v2dec = h_v2dec->GetXaxis()->FindBin(pt_bin_center_first_bin);
    Int_t ibin_first_Rgamma = h_Rgamma->GetXaxis()->FindBin(pt_bin_center_first_bin);
    
    for (Int_t i_pt_bin=0; i_pt_bin<n_pt_bins; ++i_pt_bin) {

      // Rgamma
      Rgamma_meas_vec(i_pt_bin) = h_Rgamma->GetBinContent(ibin_first_Rgamma + i_pt_bin);

      // Rgamma covariance matrix
      Double_t sigma_Rgamma = h_Rgamma->GetBinError(ibin_first_Rgamma + i_pt_bin);
      cov_Rgamma(i_pt_bin,i_pt_bin) = sigma_Rgamma * sigma_Rgamma;
      
      // inclusive photon v2
      v2_inc_meas_values(i_pt_bin) = h_v2inc_staterr->GetBinContent(ibin_first_v2inc + i_pt_bin);
      
      // inclusive photon covariance matrix
      Double_t sigma_v2inc_stat = h_v2inc_staterr->GetBinError(ibin_first_v2inc + i_pt_bin);
      Double_t sigma_v2inc_sys = h_v2inc_syserr->GetBinError(ibin_first_v2inc + i_pt_bin);
      Double_t sigma_v2inc_tot_squared = sigma_v2inc_stat*sigma_v2inc_stat + sigma_v2inc_sys*sigma_v2inc_sys;
      cov_v2_inc_staterr(i_pt_bin,i_pt_bin) = sigma_v2inc_stat * sigma_v2inc_stat;
      cov_v2_inc_syserr(i_pt_bin,i_pt_bin) = sigma_v2inc_sys * sigma_v2inc_sys;
      cov_v2_inc_toterr(i_pt_bin,i_pt_bin) = sigma_v2inc_tot_squared;

      // decay photon v2 (stat)
      v2_dec_meas_values_staterr(i_pt_bin) = h_v2dec->GetBinContent(ibin_first_v2dec + i_pt_bin);
      
      // decay photon covariance matrix (stat.)
      Double_t sigma_v2dec_staterr = h_v2dec->GetBinError(ibin_first_v2dec + i_pt_bin);
      if (sigma_v2dec_staterr == 0) {
	cout << "Warming: sigma_v2dec equal to 0" << endl;
	sigma_v2dec_staterr = 0.0001;
      }
      cov_v2_dec_staterr(i_pt_bin,i_pt_bin) = sigma_v2dec_staterr * sigma_v2dec_staterr;

      Int_t graph_offset_v2dec = 1;
      v2_dec_meas_values_syserr(i_pt_bin) = g_v2dec_syserr->GetY()[graph_offset_v2dec + i_pt_bin];

      Double_t sigma_v2dec_syserr =  g_v2dec_syserr->GetEYhigh()[graph_offset_v2dec + i_pt_bin];
      cov_v2_dec_syserr(i_pt_bin,i_pt_bin) = sigma_v2dec_syserr * sigma_v2dec_syserr;

      // sanity check: check whether central values of v2dec with stat and sys err are the same
      if (v2_dec_meas_values_staterr(i_pt_bin) != v2_dec_meas_values_syserr(i_pt_bin)) {
	cout << "WARNING: problem with v2dec values" << endl; 
      }
      
      
    }

    // cov_Rgamma.Print();
    v2_dec_meas_values_staterr.Print();
    v2_dec_meas_values_syserr.Print();
    
    TDirectory *cen_dir = f_out.mkdir(centr[i_centr].Data());
    cen_dir->cd();
    
    Rgamma_meas_vec.Write("Rgamma_meas_vec");
    cov_Rgamma.Write("cov_Rgamma");

    v2_inc_meas_values.Write("v2_inc_meas_values");
    cov_v2_inc_staterr.Write("cov_v2_inc_staterr");
    cov_v2_inc_syserr.Write("cov_v2_inc_syserr");
    cov_v2_inc_toterr.Write("cov_v2_inc_toterr");

    v2_dec_meas_values_staterr.Write("v2_dec_meas_values_staterr");
    cov_v2_dec_staterr.Write("cov_v2_dec_staterr");
    v2_dec_meas_values_syserr.Write("v2_dec_meas_values_syserr");
    cov_v2_dec_syserr.Write("cov_v2_dec_syserr");
       
  } // end loop over centrality

  f_out.Close();
  
}
