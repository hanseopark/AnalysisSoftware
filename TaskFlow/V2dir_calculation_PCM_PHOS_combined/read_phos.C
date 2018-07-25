//
// read PHOS v2 data and prepare input for v2dir_pcm_phos_comb.C
// Klaus Reygers
//

void read_phos() {

    // file with PHOS results
    // TFile file_phos("data_phos/PHOSFlow_v2_final_161116.root");
    TFile file_phos("data_phos/PHOS_SP_final_16102017.root");
    TFile file_phos_rgamma("data_phos/PHOSFlow_v2_final_161116.root"); // only for Rgamma

    // open output file
    TString fn_out = "data_phos/data_phos.root";
    TFile f_out(fn_out, "recreate");

    // reduce uncertinaties to make projections for the yellow report
    bool projection_for_yellow_report = true;
    const Double_t reduction_factor_staterr = 1./3.;
    
    // define correlation coefficient for systematic uncertainties (0 = no correlation, 1 = fully positively correlated)
    const Double_t corr_coeff_v2inc = 1.;
    const Double_t corr_coeff_Rgam = 0.25;  // not relevant for combined Rgamma
    
    const Int_t n_pt_bins = 16;

    // loop over centralities
    // centrality, 0-20% = cen6, 20-40% = cen3, 40-80% = cen8
    // const Int_t n_centr = 3;
    // TString centr[n_centr] = {"00-20", "20-40", "40-80"};
    // TString centr_cocktail[n_centr] = {"cen6", "cen3", "cen8"};
    // TString centr_Rgamma[n_centr] = {"0-20%", "20-40%", "40-80%"};
    const Int_t n_centr = 2;
    TString centr[n_centr] = {"00-20", "20-40"};
    TString centr_cocktail[n_centr] = {"cen6", "cen3"};
    TString centr_Rgamma[n_centr] = {"0-20%", "20-40%"};

    // inclusive photon vn
    TVectorD v2_inc_meas_values(n_pt_bins);    // data
    TMatrixDSym cov_v2_inc_staterr(n_pt_bins); // covariance matrix
    TMatrixDSym cov_v2_inc_syserr(n_pt_bins);  // covariance matrix
    TMatrixDSym cov_v2_inc_toterr(n_pt_bins);  // covariance matrix

    // decay photon vn
    TVectorD v2_dec_meas_values(n_pt_bins);    // data
    TMatrixDSym cov_v2_dec_staterr(n_pt_bins); // covariance matrix
    TMatrixDSym cov_v2_dec_syserr(n_pt_bins);  // covariance matrix
    TMatrixDSym cov_v2_dec_toterr(n_pt_bins);  // covariance matrix

    // Rgamma
    TVectorD Rgamma_meas_vec(n_pt_bins);
    TMatrixDSym cov_Rgamma_staterr(n_pt_bins);
    TMatrixDSym cov_Rgamma_syserr(n_pt_bins);
    TMatrixDSym cov_Rgamma_toterr(n_pt_bins);

    for (Int_t i_centr = 0; i_centr < n_centr; i_centr++) {

        TString dir = "PHOS_v2SP_PbPb_2760_Centrality_" + centr[i_centr] + "/";

        // // get histograms (works with PHOSFlow_v2_final_161116.root)
        // TString name_histo_v2inc_staterr = dir + "hPHOS_v2SP_gammaIncl_noK0s_PbPb_cen" + centr[i_centr] + "_Stat";
        // TString name_histo_v2inc_syserr = dir + "hPHOS_v2SP_gammaIncl_noK0s_PbPb_cen" + centr[i_centr] + "_Syst";
        // TString name_histo_v2dec = dir + "hPHOS_v2SP_Cocktail_noK0s_PbPb_cen" + centr[i_centr] + "_Stat";
        // TString name_graph_v2dec_syserr = dir + "hPHOS_v2SP_Cocktail_noK0s_PbPb_cen" + centr[i_centr] + "_Syst";
        // TString name_histo_Rgamma_staterr = dir + "hPHOS_DoubleRatio_PbPb_cen" + centr[i_centr] + "_Stat";
        // TString name_histo_Rgamma_syserr = dir + "hPHOS_DoubleRatio_PbPb_cen" + centr[i_centr] + "_Syst";

        // get histograms (works with PHOS_SP_final_16102017.root)
        TString name_histo_v2inc_staterr = "v2SP_PHOS_gammaIncl_noK0s_PbPb_cen" + centr[i_centr] + "_Stat";
        TString name_histo_v2inc_syserr = "v2SP_PHOS_gammaIncl_noK0s_PbPb_cen" + centr[i_centr] + "_Syst";
        TString name_histo_Rgamma_staterr = dir + "hPHOS_DoubleRatio_PbPb_cen" + centr[i_centr] + "_Stat";
        TString name_histo_Rgamma_syserr = dir + "hPHOS_DoubleRatio_PbPb_cen" + centr[i_centr] + "_Syst";

        TH1D *h_v2inc_staterr = (TH1D *)file_phos.Get(name_histo_v2inc_staterr);
        TH1D *h_v2inc_syserr = (TH1D *)file_phos.Get(name_histo_v2inc_syserr);

        // TH1D* h_v2dec = (TH1D*) file_phos.Get(name_histo_v2dec);
        // TGraphAsymmErrors* g_v2dec_syserr = (TGraphAsymmErrors*) file_phos.Get(name_graph_v2dec_syserr);
        TH1D *h_Rgamma_staterr = (TH1D *)file_phos_rgamma.Get(name_histo_Rgamma_staterr);
        TH1D *h_Rgamma_syserr = (TH1D *)file_phos_rgamma.Get(name_histo_Rgamma_syserr);

        // if (!h_v2dec) cout << "Error" << endl;

        Double_t pt_bin_center_first_bin = 1.; // GeV
        Int_t ibin_first_v2inc = h_v2inc_staterr->GetXaxis()->FindBin(pt_bin_center_first_bin);
        // Int_t ibin_first_v2dec = h_v2dec->GetXaxis()->FindBin(pt_bin_center_first_bin);
        Int_t ibin_first_Rgamma = h_Rgamma_staterr->GetXaxis()->FindBin(pt_bin_center_first_bin);

        for (Int_t i_pt_bin = 0; i_pt_bin < n_pt_bins; ++i_pt_bin) {

            //
            // Rgamma
            //
            Rgamma_meas_vec(i_pt_bin) = h_Rgamma_staterr->GetBinContent(ibin_first_Rgamma + i_pt_bin);

            // Rgamma covariance matrix
            Double_t sigma_Rgamma_staterr = h_Rgamma_staterr->GetBinError(ibin_first_Rgamma + i_pt_bin);
            Double_t sigma_Rgamma_syserr = h_Rgamma_syserr->GetBinError(ibin_first_Rgamma + i_pt_bin);
            Double_t sigma_Rgamma_toterr =
                TMath::Sqrt(sigma_Rgamma_staterr * sigma_Rgamma_staterr + sigma_Rgamma_syserr * sigma_Rgamma_syserr);
            cov_Rgamma_staterr(i_pt_bin, i_pt_bin) = sigma_Rgamma_staterr * sigma_Rgamma_staterr;
            cov_Rgamma_syserr(i_pt_bin, i_pt_bin) = sigma_Rgamma_syserr * sigma_Rgamma_syserr;
            cov_Rgamma_toterr(i_pt_bin, i_pt_bin) = sigma_Rgamma_toterr * sigma_Rgamma_toterr;

            //
            // inclusive photon v2
            //
            v2_inc_meas_values(i_pt_bin) = h_v2inc_staterr->GetBinContent(ibin_first_v2inc + i_pt_bin);

            // inclusive photon covariance matrix
            Double_t sigma_v2inc_stat = h_v2inc_staterr->GetBinError(ibin_first_v2inc + i_pt_bin);
	    // Double_t sigma_v2inc_stat = h_v2inc_staterr->GetBinError(ibin_first_v2inc + i_pt_bin);  // just a check
            Double_t sigma_v2inc_sys = h_v2inc_syserr->GetBinError(ibin_first_v2inc + i_pt_bin);

	    if (projection_for_yellow_report) {
		cout << "ATTENTION: reduced uncertainties (projections for yellow report)" << endl;
		sigma_v2inc_stat *= reduction_factor_staterr;
	    }

            Double_t sigma_v2inc_tot_squared = sigma_v2inc_stat * sigma_v2inc_stat + sigma_v2inc_sys * sigma_v2inc_sys;
	    
            cov_v2_inc_staterr(i_pt_bin, i_pt_bin) = sigma_v2inc_stat * sigma_v2inc_stat;
            cov_v2_inc_syserr(i_pt_bin, i_pt_bin) = sigma_v2inc_sys * sigma_v2inc_sys;
            cov_v2_inc_toterr(i_pt_bin, i_pt_bin) = sigma_v2inc_tot_squared;
        }

        //
        // correlated uncertainties
        //
        bool correlated_uncertainties = true;
        if (correlated_uncertainties) {
            for (Int_t i_pt_bin = 0; i_pt_bin < n_pt_bins; ++i_pt_bin) {
                for (Int_t j_pt_bin = 0; j_pt_bin < n_pt_bins; ++j_pt_bin) {
                    if (i_pt_bin != j_pt_bin) {

                        // Rgamma
                        Double_t sig_i_Rgamma_syserr = TMath::Sqrt(cov_Rgamma_syserr(i_pt_bin, i_pt_bin));
                        Double_t sig_j_Rgamma_syserr = TMath::Sqrt(cov_Rgamma_syserr(j_pt_bin, j_pt_bin));
                        cov_Rgamma_syserr(i_pt_bin, j_pt_bin) = corr_coeff_Rgam * sig_i_Rgamma_syserr * sig_j_Rgamma_syserr;
                        cov_Rgamma_toterr(i_pt_bin, j_pt_bin) = corr_coeff_Rgam * sig_i_Rgamma_syserr * sig_j_Rgamma_syserr;

                        // v2inc (assumption: systematic uncertainties is fully correlated among pT bins)
                        Double_t sig_i_v2_inc_syserr = TMath::Sqrt(cov_v2_inc_syserr(i_pt_bin, i_pt_bin));
                        Double_t sig_j_v2_inc_syserr = TMath::Sqrt(cov_v2_inc_syserr(j_pt_bin, j_pt_bin));
                        cov_v2_inc_syserr(i_pt_bin, j_pt_bin) = corr_coeff_v2inc * sig_i_v2_inc_syserr * sig_j_v2_inc_syserr;
                        cov_v2_inc_toterr(i_pt_bin, j_pt_bin) = corr_coeff_v2inc * sig_i_v2_inc_syserr * sig_j_v2_inc_syserr;
                    }
                }
            }
        }

        // cov_Rgamma.Print();
        // v2_dec_meas_values.Print();
        // v2_dec_meas_values.Print();

        TDirectory *cen_dir = f_out.mkdir(centr[i_centr].Data());
        cen_dir->cd();

        Rgamma_meas_vec.Write("Rgamma_meas_vec");
        cov_Rgamma_staterr.Write("cov_Rgamma_staterr");
        cov_Rgamma_syserr.Write("cov_Rgamma_syserr");
        cov_Rgamma_toterr.Write("cov_Rgamma_toterr");

        v2_inc_meas_values.Write("v2_inc_meas_values");
        cov_v2_inc_staterr.Write("cov_v2_inc_staterr");
        cov_v2_inc_syserr.Write("cov_v2_inc_syserr");
        cov_v2_inc_toterr.Write("cov_v2_inc_toterr");

	TObjString correlation_coeff_v2inc(Form("correlation coefficient for v2inc: %4.2f", corr_coeff_v2inc));
	correlation_coeff_v2inc.Write("correlation_coeff_v2inc");

	TObjString correlation_coeff_Rgam(Form("correlation coefficient for Rgam: %4.2f", corr_coeff_Rgam));
	correlation_coeff_Rgam.Write("correlation_coeff_Rgam");

    } // end loop over centrality

    f_out.Close();

    cout << "read_phos.C: created " << fn_out << endl;
}
