//
// read PCM v2 data and prepare input for v2dir_pcm_phos_comb.C
// Klaus Reygers
//

void read_pcm() {

    // read PCM inclusive photon v2 and Rgamma
    // and prepare input for merging macro

    // reduce uncertinaties to make projections for the yellow report
    bool projection_for_yellow_report = true;
    const Double_t reduction_factor_staterr = 1./3.;
    const Double_t reduction_factor_Rgamma_comb_syserr = 1./2.;
    
    // define correlation coefficient for systematic uncertainties (0 = no correlation, 1 = fully positively correlated)
    const Double_t corr_coeff_v2inc = 1.;
    // const Double_t corr_coeff_Rgam = 0.25; // estimated in a meeting with Dmitri and Mike on April 22, 2018
    // const Double_t corr_coeff_v2inc = 0.5; // estimated in a meeting with Dmitri and Mike on April 22, 2018

    const Double_t corr_coeff_Rgam = 1.; // a check

    Double_t pt_bin_min_Rgam_corr = 0; // default
    // Double_t pt_bin_min_Rgam_corr = 3; // no correlation in Rgamma for bins below this bin (just a test)
    
    const Int_t n_pt_bins = 16;

    // define position in Tgraphs where data start
    const Int_t i_bin_offset = 10;

    // open output file
    TString fn_out = "data_pcm/data_pcm.root";
    TFile f_out(fn_out, "recreate");

    // open Rgamma file
    TString filename_Rgamma = "data_pcm/Gamma_CombResults_PbPb_2.76TeV.root";
    TFile file_Rgamma(filename_Rgamma);

    // loop over centralities
    // centrality, 0-20% = cen6, 20-40% = cen3, 40-80% = cen8
    const Int_t n_centr = 2;
    TString centr[n_centr] = {"020", "2040"};
    TString centr_out[n_centr] = {"00-20", "20-40"};
    TString centr_cocktail[n_centr] = {"cen6", "cen3"};
    TString centr_Rgamma[n_centr] = {"0-20%", "20-40%"};

    for (Int_t i_centr = 0; i_centr < n_centr; i_centr++) {

        // open v2 inc file
        // TString filename_v2inc = "data_pcm/PCM_InclusivePhotonFlow_Syst_" + centr[i_centr] + "_170508.root";
        // TString filename_v2inc = "data_pcm/PCM_InclusivePhotonFlow_Syst_" + centr[i_centr] + "_171012.root";
        // TString filename_v2inc = "data_pcm/PCM_InclusivePhotonFlow_Syst_" + centr[i_centr] + "_171031.root";

        // switch to new convention with date first
        TString filename_v2inc = "data_pcm/171106_PCM_InclusivePhotonFlow_Syst_" + centr[i_centr] + ".root";    

        TFile file_v2inc(filename_v2inc);

        // get graphs / histograms
        TString name_graph_v2inc_syserr = "graphTotalSystErrors_" + centr[i_centr];
        TString name_graph_v2inc_toterr = "graphTotalErrors_" + centr[i_centr];
	
        TString name_graph_Rgamma_comb_toterr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_comb_totErr";
        TString name_graph_Rgamma_comb_staterr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_comb_StatErr";
        TString name_histo_Rgamma_staterr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/hDR_PCM_StatErr";
        TString name_graph_Rgamma_sysAerr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_PCM_SysAErr";
        TString name_graph_Rgamma_sysBerr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_PCM_SysBErr";
        TString name_graph_Rgamma_sysCerr = "Gamma_PbPb_2.76TeV_" + centr_Rgamma[i_centr] + "/DR_PCM_SysCErr";

        TGraphAsymmErrors *g_v2inc_syserr = (TGraphAsymmErrors *)file_v2inc.Get(name_graph_v2inc_syserr);
        TGraphAsymmErrors *g_v2inc_toterr = (TGraphAsymmErrors *)file_v2inc.Get(name_graph_v2inc_toterr);

        TGraphAsymmErrors *g_Rgamma_comb_toterr = (TGraphAsymmErrors *)file_Rgamma.Get(name_graph_Rgamma_comb_toterr);
        TGraphAsymmErrors *g_Rgamma_comb_staterr = (TGraphAsymmErrors *)file_Rgamma.Get(name_graph_Rgamma_comb_staterr);
        TGraphAsymmErrors *g_Rgamma_sysAerr = (TGraphAsymmErrors *)file_Rgamma.Get(name_graph_Rgamma_sysAerr);
        TGraphAsymmErrors *g_Rgamma_sysBerr = (TGraphAsymmErrors *)file_Rgamma.Get(name_graph_Rgamma_sysBerr);
        TGraphAsymmErrors *g_Rgamma_sysCerr = (TGraphAsymmErrors *)file_Rgamma.Get(name_graph_Rgamma_sysCerr);
        TH2F *h_Rgamma_staterr = (TH2F *)file_Rgamma.Get(name_histo_Rgamma_staterr);

        // pT values
        TVectorD pt(n_pt_bins);

        // inclusive photon vn
        TVectorD v2_inc_meas_values(n_pt_bins);    // data
        TMatrixDSym cov_v2_inc_syserr(n_pt_bins);  // covariance matrix
        TMatrixDSym cov_v2_inc_toterr(n_pt_bins);  // covariance matrix
        TMatrixDSym cov_v2_inc_staterr(n_pt_bins); // covariance matrix

        // Rgamma
        TVectorD Rgamma_meas_vec_comb(n_pt_bins); // PCM+PHOS combined
        TVectorD Rgamma_meas_vec_pcm(n_pt_bins);  // PCM
        TMatrixDSym cov_Rgamma_comb_toterr(n_pt_bins);
        TMatrixDSym cov_Rgamma_comb_staterr(n_pt_bins);
        TMatrixDSym cov_Rgamma_comb_syserr(n_pt_bins);
        TMatrixDSym cov_Rgamma_pcm_syserr(n_pt_bins);
        TMatrixDSym cov_Rgamma_pcm_staterr(n_pt_bins);
        TMatrixDSym cov_Rgamma_pcm_toterr(n_pt_bins);

        //
        // define central values and uncertainties of Rgamma, v2,inc, and v2,dec
        // (only diagonals of the covariance matrices will be filled here)
        //
        for (Int_t i_pt_bin = 0; i_pt_bin < n_pt_bins; ++i_pt_bin) {

            // Rgamma
            Rgamma_meas_vec_comb(i_pt_bin) = g_Rgamma_comb_toterr->GetY()[i_pt_bin];
            Rgamma_meas_vec_pcm(i_pt_bin) = g_Rgamma_sysAerr->GetY()[i_pt_bin];

            // Rgamma uncertainty
            Double_t sigma_Rgamma_comb_toterr = g_Rgamma_comb_toterr->GetEYhigh()[i_pt_bin];
            Double_t sigma_Rgamma_comb_staterr = g_Rgamma_comb_staterr->GetEYhigh()[i_pt_bin];
	    Double_t sigma_Rgamma_comb_syserr = sigma_Rgamma_comb_toterr * sigma_Rgamma_comb_toterr -
		sigma_Rgamma_comb_staterr * sigma_Rgamma_comb_staterr;
	    	    
	    if (projection_for_yellow_report) {
		cout << "ATTENTION: reduced uncertainties (projections for yellow report)" << endl;
		sigma_Rgamma_comb_staterr *= reduction_factor_staterr;
		sigma_Rgamma_comb_syserr *= reduction_factor_Rgamma_comb_syserr;
		sigma_Rgamma_comb_toterr =
		    TMath::Sqrt(sigma_Rgamma_comb_staterr * sigma_Rgamma_comb_staterr + sigma_Rgamma_comb_syserr * sigma_Rgamma_comb_syserr);
	    }

	    cov_Rgamma_comb_toterr(i_pt_bin, i_pt_bin) = sigma_Rgamma_comb_toterr * sigma_Rgamma_comb_toterr;
            cov_Rgamma_comb_staterr(i_pt_bin, i_pt_bin) = sigma_Rgamma_comb_staterr * sigma_Rgamma_comb_staterr;
            cov_Rgamma_comb_syserr(i_pt_bin, i_pt_bin) = sigma_Rgamma_comb_syserr * sigma_Rgamma_comb_syserr;

            Double_t sys_err_A = g_Rgamma_sysAerr->GetEYhigh()[i_pt_bin];
            Double_t sys_err_B = g_Rgamma_sysBerr->GetEYhigh()[i_pt_bin];
            Double_t sys_err_C = g_Rgamma_sysCerr->GetEYhigh()[i_pt_bin];
            Double_t sys_err_squared = sys_err_A * sys_err_A + sys_err_B * sys_err_B + sys_err_C * sys_err_C;
            cov_Rgamma_pcm_syserr(i_pt_bin, i_pt_bin) = sys_err_squared;

            Double_t ptRg = g_Rgamma_sysAerr->GetX()[i_pt_bin];
            Int_t ib = h_Rgamma_staterr->FindBin(ptRg);
            Double_t stat_err = h_Rgamma_staterr->GetBinError(ib);
            cov_Rgamma_pcm_staterr(i_pt_bin, i_pt_bin) = stat_err * stat_err;

            Double_t tot_err_squared = sys_err_squared + stat_err * stat_err;
            cov_Rgamma_pcm_toterr(i_pt_bin, i_pt_bin) = tot_err_squared;

            // v2inc
            v2_inc_meas_values(i_pt_bin) = g_v2inc_toterr->GetY()[i_bin_offset + i_pt_bin];
            if (v2_inc_meas_values(i_pt_bin) < 0)
                cout << "WARNING: v2inc < 0";

            // v2inc uncertainty
            Double_t sigma_v2inc_sys = g_v2inc_syserr->GetEYhigh()[i_bin_offset + i_pt_bin];
            Double_t sigma_v2inc_tot = g_v2inc_toterr->GetEYhigh()[i_bin_offset + i_pt_bin];
            Double_t sigma_v2inc_stat =
                TMath::Sqrt(sigma_v2inc_tot * sigma_v2inc_tot - sigma_v2inc_sys * sigma_v2inc_sys);

	    if (projection_for_yellow_report) {
		cout << "ATTENTION: reduced uncertainties (projections for yellow report)" << endl;
		sigma_v2inc_stat *= reduction_factor_staterr;
		sigma_v2inc_tot = TMath::Sqrt(sigma_v2inc_stat * sigma_v2inc_stat + sigma_v2inc_sys * sigma_v2inc_sys);
	    }
	    
            cov_v2_inc_syserr(i_pt_bin, i_pt_bin) = sigma_v2inc_sys * sigma_v2inc_sys;
            cov_v2_inc_toterr(i_pt_bin, i_pt_bin) = sigma_v2inc_tot * sigma_v2inc_tot;
            cov_v2_inc_staterr(i_pt_bin, i_pt_bin) = sigma_v2inc_stat * sigma_v2inc_stat;

            // pt
            pt(i_pt_bin) = g_v2inc_toterr->GetX()[i_bin_offset + i_pt_bin];
        }

        //
        // Correlated uncertainties: off-diagonal elements of the covariance matrices
        // Assumption: systematic are fully (positively) correlated among pT bins
        //
        bool correlated_uncertainties = true;
        if (correlated_uncertainties) {
            for (Int_t i_pt_bin =  0; i_pt_bin < n_pt_bins; ++i_pt_bin) {
                for (Int_t j_pt_bin = 0; j_pt_bin < n_pt_bins; ++j_pt_bin) {
                    if (i_pt_bin != j_pt_bin) {

                        // Rgamma
                        Double_t sig_i_Rgamma_syserr = TMath::Sqrt(cov_Rgamma_comb_syserr(i_pt_bin, i_pt_bin));
                        Double_t sig_j_Rgamma_syserr = TMath::Sqrt(cov_Rgamma_comb_syserr(j_pt_bin, j_pt_bin));
			if (i_pt_bin >= pt_bin_min_Rgam_corr && j_pt_bin >= pt_bin_min_Rgam_corr)
			    cov_Rgamma_comb_toterr(i_pt_bin, j_pt_bin) = corr_coeff_Rgam * sig_i_Rgamma_syserr * sig_j_Rgamma_syserr;
						
                        Double_t sig_i_Rgamma_pcm_syserr = TMath::Sqrt(cov_Rgamma_pcm_syserr(i_pt_bin, i_pt_bin));
                        Double_t sig_j_Rgamma_pcm_syserr = TMath::Sqrt(cov_Rgamma_pcm_syserr(j_pt_bin, j_pt_bin));
                        cov_Rgamma_pcm_toterr(i_pt_bin, j_pt_bin) = corr_coeff_Rgam * sig_i_Rgamma_pcm_syserr * sig_j_Rgamma_pcm_syserr;

                        // v2inc
                        Double_t sig_i_v2inc_syserr = TMath::Sqrt(cov_v2_inc_syserr(i_pt_bin, i_pt_bin));
                        Double_t sig_j_v2inc_syserr = TMath::Sqrt(cov_v2_inc_syserr(j_pt_bin, j_pt_bin));
                        cov_v2_inc_syserr(i_pt_bin, j_pt_bin) = corr_coeff_v2inc * sig_i_v2inc_syserr * sig_j_v2inc_syserr;
                        cov_v2_inc_toterr(i_pt_bin, j_pt_bin) = corr_coeff_v2inc * sig_i_v2inc_syserr * sig_j_v2inc_syserr;
                    }
                }
            }
        }

        // write output
        TDirectory *cen_dir = f_out.mkdir(centr_out[i_centr].Data());
        cen_dir->cd();

        Rgamma_meas_vec_comb.Write("Rgamma_meas_vec_comb"); // PCM+PHOS combined
        Rgamma_meas_vec_pcm.Write("Rgamma_meas_vec");       // only PCM
        cov_Rgamma_comb_toterr.Write("cov_Rgamma_comb_toterr");
        cov_Rgamma_comb_staterr.Write("cov_Rgamma_comb_staterr");
        cov_Rgamma_pcm_syserr.Write("cov_Rgamma_syserr");
        cov_Rgamma_pcm_staterr.Write("cov_Rgamma_staterr");
        cov_Rgamma_pcm_toterr.Write("cov_Rgamma_toterr");

        v2_inc_meas_values.Write("v2_inc_meas_values");
        cov_v2_inc_syserr.Write("cov_v2_inc_syserr");
        cov_v2_inc_toterr.Write("cov_v2_inc_toterr");
        cov_v2_inc_staterr.Write("cov_v2_inc_staterr");

	TObjString correlation_coeff_v2inc(Form("correlation coefficient for v2inc: %4.2f", corr_coeff_v2inc));
	correlation_coeff_v2inc.Write("correlation_coeff_v2inc");

	TObjString correlation_coeff_Rgam(Form("correlation coefficient for Rgam: %4.2f", corr_coeff_Rgam));
	correlation_coeff_Rgam.Write("correlation_coeff_Rgam");
	
        pt.Write("pt");

    } // end loop over centrality

    f_out.Close();

    cout << "read_pcm.C: created " << fn_out << endl;
}
