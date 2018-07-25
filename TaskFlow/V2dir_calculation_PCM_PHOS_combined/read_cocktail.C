void read_cocktail() {

    // 16 pT bins: bin i. central pT value
    // 0  |1
    // 1  |1.2
    // 2  |1.4
    // 3  |1.6
    // 4  |1.8
    // 5  |2
    // 6  |2.2
    // 7  |2.4
    // 8  |2.6
    // 9  |2.85
    // 10 |3.15
    // 11 |3.5
    // 12 |3.9
    // 13 |4.35
    // 14 |5
    // 15 |5.8

    // reduce uncertinaties to make projections for the yellow report
    bool projection_for_yellow_report = true;
    const reduction_factor_v2dec_syserr = 1.; // could add reduction for v2dec cocktail here

    // open cocktail file
    // TString filename_v2dec = "cocktail/CocktailV2_FullSys_19092017.root";
    // // TString filename_v2dec = "cocktail/CocktailV2_19092017.root";
    TString filename_v2dec = "cocktail/CocktailV2_16102017.root";
    TFile file_v2dec(filename_v2dec);

    const Double_t corr_coeff_v2dec = 1; // not 1 in order to avoid numerical instabilities(!?!)

    // loop over centralities
    // centrality, 0-20% = cen6, 20-40% = cen3, 40-80% = cen8
    const Int_t n_centr = 2;
    TString centr[n_centr] = {"020", "2040"};
    TString centr_out[n_centr] = {"00-20", "20-40"};
    TString centr_cocktail[n_centr] = {"cen6", "cen3"};
    TString centr_Rgamma[n_centr] = {"0-20%", "20-40%"};

    // number of pt bins
    const Int_t n_pt_bins = 16;

    // first bin in graphs
    const Int_t index_first_point = 1; // for TGraph
    const Int_t first_bin = 2;         // for histogram with statistical uncertainties

    TVectorD v2_dec_meas_values(n_pt_bins);    // data
    TMatrixDSym cov_v2_dec_staterr(n_pt_bins); // covariance matrix, stat. error
    TMatrixDSym cov_v2_dec_syserr(n_pt_bins);  // covariance matrix, sys. error
    TMatrixDSym cov_v2_dec_toterr(n_pt_bins);  // covariance matrix, sys. error

    // open output file
    TString fn_out = "cocktail/cocktail.root";
    TFile f_out(fn_out, "recreate");

    for (Int_t i_centr = 0; i_centr < n_centr; i_centr++) {

        TString name_histo_v2dec_staterr = "v2gammaCocktailNoK0s_" + centr_cocktail[i_centr];
        TString name_graph_v2dec_syserr = "v2gammaCocktailNoK0s_sys_" + centr_cocktail[i_centr];

        // decay photon cocktail
        TH1D *h_v2dec_staterr = (TH1D *)file_v2dec.Get(name_histo_v2dec_staterr);
        TGraphAsymmErrors *g_v2dec_syserr = (TGraphAsymmErrors *)file_v2dec.Get(name_graph_v2dec_syserr);

        for (Int_t i_pt_bin = 0; i_pt_bin < n_pt_bins; ++i_pt_bin) {

            // v2dec
            v2_dec_meas_values(i_pt_bin) = g_v2dec_syserr->GetY()[index_first_point + i_pt_bin];

            // v2dec uncertainty
            Double_t sigma_v2dec_stat = h_v2dec_staterr->GetBinError(first_bin + i_pt_bin);

            // cout << "WARNING: just a test: increased stat. error by a factor 2 !!!" << endl;
            // sigma_v2dec_stat *= 2;

            Double_t sigma_v2dec_sys = g_v2dec_syserr->GetEYhigh()[index_first_point + i_pt_bin];

	    if (projection_for_yellow_report) {
		cout << "ATTENTION: reduced uncertainties (projections for yellow report)" << endl;
		sigma_v2dec_sys *= reduction_factor_v2dec_syserr;
	    }

	    Double_t sigma_v2dec_tot_squared = sigma_v2dec_stat * sigma_v2dec_stat + sigma_v2dec_sys * sigma_v2dec_sys;

            cov_v2_dec_staterr(i_pt_bin, i_pt_bin) = sigma_v2dec_stat * sigma_v2dec_stat;
            cov_v2_dec_syserr(i_pt_bin, i_pt_bin) = sigma_v2dec_sys * sigma_v2dec_sys;
            cov_v2_dec_toterr(i_pt_bin, i_pt_bin) = sigma_v2dec_tot_squared;
        }

        // off diagonal elements: correlated errors
        for (Int_t i_pt_bin = 0; i_pt_bin < n_pt_bins; ++i_pt_bin) {
            for (Int_t j_pt_bin = 0; j_pt_bin < n_pt_bins; ++j_pt_bin) {

                Double_t sigma_v2dec_sys_i = TMath::Sqrt(cov_v2_dec_syserr(i_pt_bin, i_pt_bin));
                Double_t sigma_v2dec_sys_j = TMath::Sqrt(cov_v2_dec_syserr(j_pt_bin, j_pt_bin));

                if (i_pt_bin != j_pt_bin) {
                    cov_v2_dec_syserr(i_pt_bin, j_pt_bin) = corr_coeff_v2dec * sigma_v2dec_sys_i * sigma_v2dec_sys_j;
                    cov_v2_dec_toterr(i_pt_bin, j_pt_bin) = cov_v2_dec_syserr(i_pt_bin, j_pt_bin);
                }
            }
        }

        // write output
        TDirectory *cen_dir = f_out.mkdir(centr_out[i_centr].Data());
        cen_dir->cd();

        v2_dec_meas_values.Write("v2_dec_meas_values");
        cov_v2_dec_staterr.Write("cov_v2_dec_staterr");
        cov_v2_dec_syserr.Write("cov_v2_dec_syserr");
        cov_v2_dec_toterr.Write("cov_v2_dec_toterr");

        // cov_v2_dec_toterr.Print();
    }

    f_out.Close();

    cout << "read_cocktail.C: created " << fn_out << endl;
}
