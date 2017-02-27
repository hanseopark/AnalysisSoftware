void read_pcm() {

  // read PCM inclusive photon v2, decay photon v2, and Rgamm
  // and prepare input for merging macro

  const Int_t n_pt_bins = 16;

  // define position in Tgraphs where data start
  const Int_t i_bin_offset = 10;

  // open output file
  TFile f_out("data_pcm/data_pcm.root","recreate");
  
  // open Rgamma file
  TString filename_Rgamma = "data_pcm/Gamma_CombResults_PbPb_2.76TeV.root";
  TFile file_Rgamma(filename_Rgamma);
  
  // open cocktail file
  TString filename_v2dec = "data_pcm/CocktailV2.root";
  TFile file_v2dec(filename_v2dec);

  // loop over centralities
  // centrality, 0-20% = cen6, 20-40% = cen3, 40-80% = cen8
  const Int_t n_centr = 3;
  TString centr[n_centr] = {"020", "2040", "4080"};
  TString centr_out[n_centr] = {"00-20", "20-40", "40-80"};
  TString centr_cocktail[n_centr] = {"cen6", "cen3", "cen8"}; 
  TString centr_Rgamma[n_centr] = {"0-20%", "20-40%", "40-80%"};

  for (Int_t i_centr=0; i_centr<n_centr; i_centr++) {
    
    // open v2 inc file
    TString filename_v2inc = "data_pcm/PCM_InclusivePhotonFlow_Syst_" + centr[i_centr] + ".root";
    TFile file_v2inc(filename_v2inc);
  
    // get graphs / histograms
    TString name_graph_v2inc_syserr = "graphTotalSystErrors_" + centr[i_centr];
    TString name_graph_v2inc_toterr = "graphTotalErrors_" + centr[i_centr];
    TString name_histo_v2dec = "v2gammaCocktail_" + centr_cocktail[i_centr];
    TString name_graph_Rgamma_comb_toterr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_comb_totErr";
    TString name_graph_Rgamma_comb_staterr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_comb_StatErr";
    TString name_histo_Rgamma_staterr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/hDR_PCM_StatErr";
    TString name_graph_Rgamma_sysAerr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_PCM_SysAErr";
    TString name_graph_Rgamma_sysBerr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_PCM_SysBErr";
    TString name_graph_Rgamma_sysCerr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_PCM_SysCErr";
  
    TGraphAsymmErrors* g_v2inc_syserr = (TGraphAsymmErrors*) file_v2inc.Get(name_graph_v2inc_syserr);
    TGraphAsymmErrors* g_v2inc_toterr = (TGraphAsymmErrors*) file_v2inc.Get(name_graph_v2inc_toterr);
    TH1D* h_v2dec = (TH1D*) file_v2dec.Get(name_histo_v2dec);
    h_v2dec->DrawCopy();
    g_v2inc_toterr->DrawClone("p");

    // find first bin in dec photon v2 histogram
    Double_t pt_first_bin = g_v2inc_toterr->GetX()[i_bin_offset];
    Int_t first_bin_h_v2dec = h_v2dec->GetXaxis()->FindBin(pt_first_bin);
    if (h_v2dec->GetBinCenter(first_bin_h_v2dec) != pt_first_bin) {
      cout << "ERROR: pt bin mismatch" << endl;
    }

    TGraphAsymmErrors* g_Rgamma_comb_toterr = (TGraphAsymmErrors*) file_Rgamma.Get(name_graph_Rgamma_comb_toterr);
    TGraphAsymmErrors* g_Rgamma_comb_staterr = (TGraphAsymmErrors*) file_Rgamma.Get(name_graph_Rgamma_comb_staterr);
    TGraphAsymmErrors* g_Rgamma_sysAerr = (TGraphAsymmErrors*) file_Rgamma.Get(name_graph_Rgamma_sysAerr);
    TGraphAsymmErrors* g_Rgamma_sysBerr = (TGraphAsymmErrors*) file_Rgamma.Get(name_graph_Rgamma_sysBerr);
    TGraphAsymmErrors* g_Rgamma_sysCerr = (TGraphAsymmErrors*) file_Rgamma.Get(name_graph_Rgamma_sysCerr);

    
    // pT values
    TVectorD pt(n_pt_bins);

    // inclusive photon vn
    TVectorD v2_inc_meas_values(n_pt_bins); // data
    TMatrixDSym cov_v2_inc_syserr(n_pt_bins); // covariance matrix
    TMatrixDSym cov_v2_inc_toterr(n_pt_bins); // covariance matrix
    TMatrixDSym cov_v2_inc_staterr(n_pt_bins); // covariance matrix
  
    // decay photon vn
    TVectorD v2_dec_meas_values(n_pt_bins); // data
    TMatrixDSym cov_v2_dec(n_pt_bins); // covariance matrix

    // Rgamma
    TVectorD Rgamma_meas_vec_comb(n_pt_bins); // PCM+PHOS combined
    TVectorD Rgamma_meas_vec(n_pt_bins);  // PCM
    TMatrixDSym cov_Rgamma_comb_toterr(n_pt_bins);
    TMatrixDSym cov_Rgamma_comb_staterr(n_pt_bins);
    TMatrixDSym cov_Rgamma_syserr(n_pt_bins);

    for (Int_t i_pt_bin=0; i_pt_bin<n_pt_bins; ++i_pt_bin) {

      // Rgamma
      Rgamma_meas_vec_comb(i_pt_bin) = g_Rgamma_comb_toterr->GetY()[i_pt_bin];
      Rgamma_meas_vec(i_pt_bin) = g_Rgamma_sysAerr->GetY()[i_pt_bin];
      
      // Rgamma uncertainty
      Double_t sigma_Rgamma_comb_toterr = g_Rgamma_comb_toterr->GetEYhigh()[i_pt_bin];
      Double_t sigma_Rgamma_comb_staterr = g_Rgamma_comb_staterr->GetEYhigh()[i_pt_bin];
      cov_Rgamma_comb_toterr(i_pt_bin, i_pt_bin) = sigma_Rgamma_comb_toterr * sigma_Rgamma_comb_toterr;
      cov_Rgamma_comb_staterr(i_pt_bin, i_pt_bin) = sigma_Rgamma_comb_staterr * sigma_Rgamma_comb_staterr;

      Double_t sys_err_A = g_Rgamma_sysAerr->GetEYhigh()[i_pt_bin];
      Double_t sys_err_B = g_Rgamma_sysBerr->GetEYhigh()[i_pt_bin];
      Double_t sys_err_C = g_Rgamma_sysCerr->GetEYhigh()[i_pt_bin];
      Double_t sys_err_squared = sys_err_A*sys_err_A + sys_err_B*sys_err_B + sys_err_C*sys_err_C;
      cov_Rgamma_syserr(i_pt_bin, i_pt_bin) = sys_err_squared;
      
      // v2inc
      v2_inc_meas_values(i_pt_bin) = g_v2inc_toterr->GetY()[i_bin_offset + i_pt_bin];
      if (v2_inc_meas_values(i_pt_bin) < 0) cout << "WARNING: v2inc < 0";

      // v2inc uncertainty
      Double_t sigma_v2inc_sys = g_v2inc_syserr->GetEYhigh()[i_bin_offset + i_pt_bin];
      Double_t sigma_v2inc_tot = g_v2inc_toterr->GetEYhigh()[i_bin_offset + i_pt_bin];
      Double_t sigma_v2inc_stat = TMath::Sqrt(sigma_v2inc_tot*sigma_v2inc_tot - sigma_v2inc_sys*sigma_v2inc_sys);
      cov_v2_inc_syserr(i_pt_bin, i_pt_bin) = sigma_v2inc_sys * sigma_v2inc_sys;
      cov_v2_inc_toterr(i_pt_bin, i_pt_bin) = sigma_v2inc_tot * sigma_v2inc_tot;
      cov_v2_inc_staterr(i_pt_bin, i_pt_bin) = sigma_v2inc_stat * sigma_v2inc_stat;

      // v2dec
      v2_dec_meas_values(i_pt_bin) = h_v2dec->GetBinContent(first_bin_h_v2dec + i_pt_bin);

      // v2dec uncertainty
      Double_t sigma_v2dec = h_v2dec->GetBinError(first_bin_h_v2dec + i_pt_bin);
      cov_v2_dec(i_pt_bin, i_pt_bin) = sigma_v2dec * sigma_v2dec;

      // pt
      pt(i_pt_bin) = g_v2inc_toterr->GetX()[i_bin_offset + i_pt_bin];
      
    }
  
    // v2_dec_meas_values.Print();
    // cov_v2_dec.Print();
    // Rgamma_meas_vec.Print();
    // cov_Rgamma.Print();

    TDirectory *cen_dir = f_out.mkdir(centr_out[i_centr].Data());
    cen_dir->cd();
    
    Rgamma_meas_vec_comb.Write("Rgamma_meas_vec_comb");  // PCM+PHOS combined
    Rgamma_meas_vec.Write("Rgamma_meas_vec"); // only PCM
    cov_Rgamma_comb_toterr.Write("cov_Rgamma_comb_toterr");
    cov_Rgamma_comb_staterr.Write("cov_Rgamma_comb_staterr");
    cov_Rgamma_syserr.Write("cov_Rgamma_syserr");
    
    v2_inc_meas_values.Write("v2_inc_meas_values");
    cov_v2_inc_syserr.Write("cov_v2_inc_syserr");
    cov_v2_inc_toterr.Write("cov_v2_inc_toterr");
    cov_v2_inc_staterr.Write("cov_v2_inc_staterr");
    
    v2_dec_meas_values.Write("v2_dec_meas_values");
    cov_v2_dec.Write("cov_v2_dec");

    pt.Write("pt");

  }  // end loop over centrality
  
  f_out.Close();

}
