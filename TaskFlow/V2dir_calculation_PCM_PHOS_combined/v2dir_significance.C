//
// Significance of certain v2dir hypothesis based on marginalized likelihoods
//
// Method:
// Generate v2inc pseudo data for the hypothesis v2dir = 0 and compare with the actually
// measured v2inc data. 
// 
// Klaus Reygers, March 2018
//

// main setting: significance calculation with or without remapping of sampled Rgamma value
// bool old_method = true;  // for v2dir paper on arXiv (before QM 2018)
bool old_method = false;

// decide whether to write v2inc pseudo data into a tree
bool v2_inc_pd_tree = false;

//
// used functions
//
Double_t dot_product(const TVectorD &a, const TVectorD &b);

RooDataSet *generate_multivariate_pseudo_data(const TVectorD &mean_vect, const TMatrixDSym &cov_matrix,
                                              const Double_t &range_min, const Double_t &range_max,
                                              const Int_t n_samples);

void update_sums(RooDataSet *pd, Int_t &n_evt, TVectorD &v2inc_mean_sum, TMatrixDSym &v2inc_cov_sum);

void get_v2inc_pseudo_data(RooDataSet *pd, Double_t* v2inc_pd, const Int_t n_pt_bins);

Double_t p_value_to_n_sigma(const Double_t &pvalue);

Double_t cdf_standard_normal_distr(const Double_t& x);

Double_t inverse_cdf_standard_normal_distr(const Double_t& x);

// for visualization of the correlation of the v2inc pseudo-data for the first two pt bins
// TH2F hv2inc_1_2("hv2inc_1_2", "v2inc_pt_bins_2 vs. v2inc_pt_bins_1 (pseudo data)", 100, 0.05, 0.09, 100, 0.05, 0.11);
TH2F hv2inc_1_2("hv2inc_1_2", "v2inc_pt_bins_2 vs. v2inc_pt_bins_1 (pseudo data)", 100, 0.08, 0.22, 100, 0.08, 0.22);

//
// the main function - for single pt bin choose pt_bin_low = pt_bin_up
//
void v2dir_significance(Int_t hypothesis_id = 0, TString centr = "20-40", const Int_t pt_bin_low = 0,
                        const Int_t pt_bin_up = 15, Int_t n_samples = 1000, TString output_dir = ".") {

    if (pt_bin_low < 0 || pt_bin_up > 15 || pt_bin_low > pt_bin_up) {
        cout << "ERROR: bin numbers incorrect" << endl;
        exit(-1);
    }
    
    const Int_t n_pt_bins = pt_bin_up - pt_bin_low + 1;

    TString fn = Form("%s/v2dir_pcm_phos_comb_%s.root", output_dir.Data(), centr.Data());
    TFile f(fn.Data());

    // v2inc
    TGraphAsymmErrors g_v2_inc_comb_toterr = *(TGraphAsymmErrors *)f.Get("g_v2_inc_comb_toterr");
    TGraphAsymmErrors g_v2_inc_comb_staterr = *(TGraphAsymmErrors *)f.Get("g_v2_inc_comb_staterr");    
    TMatrixDSym cov_v2inc_all_bins = *(TMatrixDSym *)f.Get("cov_v2_inc_toterr_comb");

    //
    // just a test (playing with different correlations coefficients of the systematic uncertainties)
    //
    // for (Int_t i=0; i<16; i++) {
    //  	for (Int_t j=i+1; j<16; j++) {
    // 	    cov_v2inc_all_bins(i,j) = 0.99 * TMath::Sqrt(cov_v2inc_all_bins(i,i)) * TMath::Sqrt(cov_v2inc_all_bins(j,j)); 
    // 	    cov_v2inc_all_bins(i,j) = 0.;
    //  	    cov_v2inc_all_bins(j,i) = cov_v2inc_all_bins(i,j);
    // 	}
    // }

    TMatrixDSym cov_v2inc(n_pt_bins);
    cov_v2inc_all_bins.GetSub(pt_bin_low, pt_bin_up, pt_bin_low, pt_bin_up, cov_v2inc);
    TMatrixDSym cov_v2inc_inv = cov_v2inc;
    cov_v2inc_inv.Invert();
    // cov_v2inc.Print();
    
    // v2dec
    TGraphAsymmErrors g_v2_dec_comb = *(TGraphAsymmErrors *)f.Get("g_v2_dec_comb");
    TMatrixDSym cov_v2dec_all_bins = *(TMatrixDSym *)f.Get("cov_v2_dec_toterr");

    //
    // just a test (playing with different correlations coefficients of the systematic uncertainties)
    //
    // for (Int_t i=0; i<16; i++) {
    // 	for (Int_t j=i+1; j<16; j++) {
    //    	    // cov_v2dec_all_bins(i,j) = 0.99 * TMath::Sqrt(cov_v2dec_all_bins(i,i)) * TMath::Sqrt(cov_v2dec_all_bins(j,j));
    // 	    // cov_v2dec_all_bins(i,j) = 0.;
    //  	    // cov_v2dec_all_bins(j,i) = cov_v2dec_all_bins(i,j);
	    
    //  	}
    //  }

    TMatrixDSym cov_v2dec(n_pt_bins);
    cov_v2dec_all_bins.GetSub(pt_bin_low, pt_bin_up, pt_bin_low, pt_bin_up, cov_v2dec);

    // Rgamma
    TGraphErrors g_Rgamma_toterr = *(TGraphErrors *)f.Get("g_Rgamma_toterr");
    TMatrixDSym cov_Rgam_all_bins = *(TMatrixDSym *)f.Get("cov_Rgamma_comb_toterr");

    //
    // just a test (playing with different correlations coefficients of the systematic uncertainties)
    //
    // for (Int_t i=0; i<16; i++) {
    // 	for (Int_t j=i+1; j<16; j++) {
    // 	    cov_Rgam_all_bins(i,j) = 0.99 * TMath::Sqrt(cov_Rgam_all_bins(i,i)) * TMath::Sqrt(cov_Rgam_all_bins(j,j)); 
    // 	    cov_Rgam_all_bins(i,j) = 0.;
    // 	    cov_Rgam_all_bins(j,i) = cov_Rgam_all_bins(i,j);
    // 	}
    // }

    TMatrixDSym cov_Rgam(n_pt_bins);
    cov_Rgam_all_bins.GetSub(pt_bin_low, pt_bin_up, pt_bin_low, pt_bin_up, cov_Rgam);

    // fill measured v2inc, v2dec, and Rgam into TVectorD
    TVectorD v2inc_meas(n_pt_bins);
    TVectorD v2inc_meas_stat_err(n_pt_bins);
    TVectorD v2dec_calc(n_pt_bins);
    TVectorD Rgam_meas(n_pt_bins);
    TVectorD pt(n_pt_bins);
    for (Int_t i = 0; i < n_pt_bins; i++) {
        pt(i) = g_v2_inc_comb_toterr.GetX()[i + pt_bin_low];
        v2inc_meas(i) = g_v2_inc_comb_toterr.GetY()[i + pt_bin_low];
	v2inc_meas_stat_err(i) = g_v2_inc_comb_staterr.GetEYhigh()[i + pt_bin_low];
        v2dec_calc(i) = g_v2_dec_comb.GetY()[i + pt_bin_low];
        Rgam_meas(i) = g_Rgamma_toterr.GetY()[i + pt_bin_low];
        if (Rgam_meas(i) < 1.)
            Rgam_meas(i) = 1.;
    }

    // Paquet hydro prediction
    TVectorD v2dir_Paquet_00_20(16);
    v2dir_Paquet_00_20(0) = 0.048;
    v2dir_Paquet_00_20(1) = 0.052;
    v2dir_Paquet_00_20(2) = 0.056;
    v2dir_Paquet_00_20(3) = 0.057;
    v2dir_Paquet_00_20(4) = 0.056;
    v2dir_Paquet_00_20(5) = 0.055;
    v2dir_Paquet_00_20(6) = 0.051;
    v2dir_Paquet_00_20(7) = 0.047;
    v2dir_Paquet_00_20(8) = 0.042;
    v2dir_Paquet_00_20(9) = 0.035;

    TVectorD v2dir_Paquet_20_40(16);
    v2dir_Paquet_20_40(0) = 0.083;
    v2dir_Paquet_20_40(1) = 0.092;
    v2dir_Paquet_20_40(2) = 0.098;
    v2dir_Paquet_20_40(3) = 0.101;
    v2dir_Paquet_20_40(4) = 0.099;
    v2dir_Paquet_20_40(5) = 0.094;
    v2dir_Paquet_20_40(6) = 0.087;
    v2dir_Paquet_20_40(7) = 0.079;
    v2dir_Paquet_20_40(8) = 0.070;
    v2dir_Paquet_20_40(9) = 0.057;
    
    // print centrality class and selected pt range
    cout << "Centrality class " << centr.Data() << "%, selected pt range: " << pt(0) << " - " << pt(n_pt_bins - 1)
         << " GeV/c" << endl;

    // vector for v2dir hypothesis
    TVectorD v2dir_hyp(n_pt_bins);
    TString v2dir_hyp_name = "";
    TString v2dir_hyp_name_latex = "";

    // select hypothesis
    switch (hypothesis_id) {
    case 0:
	// hypothesis: v2dir = 0 for all pt bins
        v2dir_hyp_name = "v2dir_equals_0";
	v2dir_hyp_name_latex = "v_{2,dir} = 0";
        for (Int_t i = 0; i < n_pt_bins; i++)
            v2dir_hyp(i) = 0;
        break;
    case 1:
	// hypothesis: v2dir = v2dec
        v2dir_hyp_name = "v2dir_equals_v2dec";
        for (Int_t i = 0; i < n_pt_bins; i++)
            v2dir_hyp(i) = v2dec_calc(i);
        break;
    case 2:
	// hypothesis: v2dir from (Paquet et al)
	v2dir_hyp_name = "v2dir_equals_paquet";
	if (centr == "00-20") {
	    for (Int_t i = 0; i < n_pt_bins; i++) {
		v2dir_hyp(i) = v2dir_Paquet_00_20(i + pt_bin_low);
	    }
	}
	else if (centr == "20-40") {
	    for (Int_t i = 0; i < n_pt_bins; i++) {
		v2dir_hyp(i) = v2dir_Paquet_20_40(i + pt_bin_low);
	    }
	}
	else {
	    cout << "ERROR: undefined centrality class" << endl;
	    exit(-1);
	}
	break;
	
    }

    // histogram for test statistic
    TH1F h_test_statistic("h_test_statistic", "h_test_statistic", 1000, 0., 500.);

    Int_t n_events = 0;
    TVectorD v2inc_mean_sum(n_pt_bins);
    TMatrixDSym v2inc_cov_sum(n_pt_bins);

    // v2inc pseudo data
    RooDataSet *v2inc_pseudo_data_all = 0;

    // sample Rgam posterior distribution (new method: Rgam > 1 limit not realized via parameter limits)
    Double_t Rgam_lower_limit = 0;
    if (old_method) Rgam_lower_limit = 1.;
    RooDataSet *Rgam_data_set = generate_multivariate_pseudo_data(Rgam_meas, cov_Rgam, Rgam_lower_limit, 5., n_samples);
 
    // sample v2dec posterior distribution
    RooDataSet *v2dec_data_set = generate_multivariate_pseudo_data(v2dec_calc, cov_v2dec, -0.5, 0.5, n_samples);

    // plot some v2inc pseudo data points (for debugging purposes)
    Double_t v2inc_pd[n_pt_bins];
    TTree *tr = NULL;
    if (v2_inc_pd_tree) {
	tr = new TTree("tr","v2inc pseudo data");    
	for (Int_t i=0; i<n_pt_bins; ++i) {
	    TString sv2inc = Form("v2inc%d", i);
	    TString sv2inc2 = sv2inc + "/D";
	    tr->Branch(sv2inc.Data(), &v2inc_pd[i], sv2inc2.Data());
	}
    }

    // loop over Rgam and v2dec values
    for (Int_t k = 0; k < n_samples; k++) {

        TVectorD v2inc_sampled(n_pt_bins);

        const RooArgSet *Rgam_row = Rgam_data_set->get(k);
        const RooArgSet *v2dec_row = v2dec_data_set->get(k);

        TVectorD Rgam_sampled(n_pt_bins);
        TVectorD v2dec_sampled(n_pt_bins);

        for (Int_t i = 0; i < n_pt_bins; i++) {

            // is there a simpler way to retrieve the pseudo data?
            RooRealVar *v2dec_var = (RooRealVar *)v2dec_row->find(((RooArgList *)v2dec_row)->at(i)->GetName());
            v2dec_sampled(i) = v2dec_var->getVal();

            RooRealVar *Rgam_var = (RooRealVar *)Rgam_row->find(((RooArgList *)Rgam_row)->at(i)->GetName());

	    if (old_method) {
		Rgam_sampled(i) = Rgam_var->getVal();
	    }
	    else {
		// new method: take into account physical limit for Rgamma here (Rgamma > 1)
		Double_t RgamValNoLimit = Rgam_var->getVal();	    
		Double_t n_sigma_RgamValNoLimt =  (RgamValNoLimit - Rgam_meas(i)) / TMath::Sqrt(cov_Rgam(i,i));
		Double_t n_sigma_1 = (1. - Rgam_meas(i)) / TMath::Sqrt(cov_Rgam(i,i));
		Double_t RgamValNoLimitQuantile = (1. - cdf_standard_normal_distr(n_sigma_1)) * 
		    cdf_standard_normal_distr(n_sigma_RgamValNoLimt) + cdf_standard_normal_distr(n_sigma_1);
		Rgam_sampled(i) = TMath::Sqrt(cov_Rgam(i,i)) * inverse_cdf_standard_normal_distr(RgamValNoLimitQuantile) + Rgam_meas(i);
		if ( Rgam_sampled(i) < 1.) cout << "What?" << endl;
	    }
	    
	    v2inc_sampled(i) = ((Rgam_sampled(i) - 1.) * v2dir_hyp(i) + v2dec_sampled(i)) / Rgam_sampled(i);
	    
        }

        // generate v2inc pseudo data for the sampled v2inc prediction
        RooDataSet *v2inc_pseudo_data = generate_multivariate_pseudo_data(v2inc_sampled, cov_v2inc, -0.5, 0.5, 1000);

	// update sums needed for the calculation of the covariance matrix if the v2inc pseudo-data geenrated
	// under the selected hypothesis
        update_sums(v2inc_pseudo_data, n_events, v2inc_mean_sum, v2inc_cov_sum);

	if (v2_inc_pd_tree) {
	    get_v2inc_pseudo_data(v2inc_pseudo_data, v2inc_pd, n_pt_bins);
	    tr->Fill();
	}

        delete v2inc_pseudo_data;
    }

    delete Rgam_data_set;
    delete v2dec_data_set;

    //
    // calculate significance
    //

    // matrices needed to calculate covariance of v2inc pseudo data
    TVectorD v2inc_mean_norm = v2inc_mean_sum * (1. / n_events);
    TMatrixDSym v2inc_cov_norm = v2inc_cov_sum * (1. / n_events);

    TMatrixDSym v2inc_pseudo_data_cov(n_pt_bins);

    // prepare covariance matrix of v2inc pseudo data
    for (Int_t i = 0; i < n_pt_bins; ++i) {
        for (Int_t j = 0; j < n_pt_bins; ++j) {
            v2inc_pseudo_data_cov(i, j) = v2inc_cov_norm(i, j) - v2inc_mean_norm(i) * v2inc_mean_norm(j);
        }
    }

    TMatrixDSym v2inc_pseudo_data_cov_inv = v2inc_pseudo_data_cov;
    v2inc_pseudo_data_cov_inv.Invert();

    // calculate significance
    TVectorD d = v2inc_meas - v2inc_mean_norm;
    Double_t chi2 = dot_product(d, v2inc_pseudo_data_cov_inv * d);
    Double_t p_value = TMath::Prob(chi2, n_pt_bins);
    Double_t n_sigma = p_value_to_n_sigma(p_value);

    // print results
    cout << "hypothesis: " << v2dir_hyp_name << endl;
    cout << "chi2 = " << chi2 << endl;
    cout << "p-value: " << p_value << endl;
    cout << "n_sigma = " << n_sigma << endl;

    // cross check: Rgamma  significance
    TMatrixDSym cov_Rgam_inv = cov_Rgam;
    cov_Rgam_inv.Invert();
    TVectorD d_Rgam(n_pt_bins);
    for (Int_t i=0; i<n_pt_bins; ++i) d_Rgam(i) = Rgam_meas(i) - 1.;
    Double_t chi2_Rgam = dot_product(d_Rgam, cov_Rgam_inv * d_Rgam);
    Double_t p_value_Rgam = TMath::Prob(chi2_Rgam, n_pt_bins);
    Double_t n_sigma_Rgam = p_value_to_n_sigma(p_value_Rgam);
    cout << endl << "n_sigma_Rgam = " << n_sigma_Rgam << endl;
    
    // plot first two bins
       
    // style settings
    TStyle *myStyle = new TStyle("myStyle", "My root style");
    myStyle->SetOptStat(kFALSE);
    myStyle->SetLabelOffset(0.005, "x"); // 0.005 = root default
    myStyle->SetLabelOffset(0.005, "y"); // 0.005 = root default
    myStyle->SetTitleXSize(0.04);        // 0.04  = root default
    myStyle->SetTitleYSize(0.04);        // 0.04  = root default
    myStyle->SetTitleXOffset(1.2);
    myStyle->SetTitleYOffset(1.3);
    myStyle->SetPadLeftMargin(0.15);
    myStyle->SetPadRightMargin(0.12); // 0.1 = root default
    myStyle->SetPadTopMargin(0.1);
    myStyle->SetPadBottomMargin(0.12);
    myStyle->SetCanvasColor(0);
    myStyle->SetPadColor(0);
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetPadBorderMode(0);
    myStyle->SetOptTitle(0);

    // another style
    TStyle *myStyle2 = new TStyle("myStyle2", "Another root style");
    myStyle->SetOptStat(kFALSE);
    myStyle->SetLabelOffset(0.005, "x"); // 0.005 = root default
    myStyle->SetLabelOffset(0.005, "y"); // 0.005 = root default
    myStyle->SetTitleXSize(0.04);        // 0.04  = root default
    myStyle->SetTitleYSize(0.04);        // 0.04  = root default
    myStyle->SetTitleXOffset(1.2);
    myStyle->SetTitleYOffset(1.3);
    myStyle->SetPadLeftMargin(0.15);
    myStyle->SetPadRightMargin(0.12); // 0.1 = root default
    myStyle->SetPadTopMargin(0.1);
    myStyle->SetPadBottomMargin(0.12);
    myStyle->SetCanvasColor(0);
    myStyle->SetPadColor(0);
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetPadBorderMode(0);
    myStyle->SetOptTitle(0);
    myStyle->SetPadTickX(1);
    myStyle->SetPadTickY(1);
	
    gROOT->SetStyle("myStyle");

    hv2inc_1_2.SetXTitle("v_{2,inc,m}(p_{T} bin 0)");
    hv2inc_1_2.SetYTitle("v_{2,inc,m}(p_{T} bin 1)");
    TCanvas* c1 = new TCanvas("c1");
    if (n_pt_bins > 1) {
	hv2inc_1_2.DrawCopy("colz");
	TMarker m(v2inc_meas(0), v2inc_meas(1), 20);
	m.SetMarkerColor(kRed);
	m.DrawClone();
    }

    TString fn1 =  output_dir + "/v2inc_pseudo_data_correlation_plot_pt_bins_1_2_" + centr + "_" + output_dir + ".png";
    c1->SaveAs(fn1.Data());
    
    // TVectorD v2inc_mean_roofit = *(v2inc_pseudo_data_all->mean());
    // TMatrixDSym v2inc_pseuo_data_cov_roofit = *(v2inc_pseudo_data_all->covarianceMatrix());
    // v2inc_pseuo_data_cov_roofit.Print();

    //
    // draw mean and sigma of the v2inc pseudo data along with the measured v2inc points
    //
    
    TGraphErrors g_v2inc_pseudo_data(n_pt_bins);
    TGraphErrors g_v2inc_meas(n_pt_bins);
    for (Int_t i=0; i<n_pt_bins; i++) {
	g_v2inc_meas.SetPoint(i, pt(i), v2inc_meas(i));
	//  g_v2inc_meas.SetPointError(i, 0, 0);
	g_v2inc_meas.SetPointError(i, 0, v2inc_meas_stat_err(i));
	
	g_v2inc_pseudo_data.SetPoint(i, pt(i), v2inc_mean_norm(i));
	g_v2inc_pseudo_data.SetPointError(i, 0, TMath::Sqrt(v2inc_pseudo_data_cov(i,i)));	
    }
    
    TCanvas* c2 = new TCanvas("c2");
    // c2->cd();
    
    TH2F fr2("fr2", "fr2", 1, pt(0)-0.25, pt(n_pt_bins-1)+0.25, 1, 0., 0.25);
    fr2.SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fr2.SetYTitle("#it{v}_{2,inc}");
    
    fr2.DrawCopy();
    g_v2inc_pseudo_data.SetMarkerColor(kRed);
    g_v2inc_pseudo_data.SetMarkerStyle(kOpenCircle);

    g_v2inc_meas.SetMarkerColor(kBlue);
    g_v2inc_meas.SetMarkerStyle(kFullCircle);
    
    g_v2inc_pseudo_data.DrawClone("p");
    g_v2inc_meas.DrawClone("p");

    TLegend leg2(0.3, 0.15, 0.82, 0.25, NULL, "brNDC");
    TString legend_text_pseudo_data = Form("pseudo data for hypothesis %s", v2dir_hyp_name_latex.Data());
    leg2.AddEntry(&g_v2inc_meas, "real data", "p");
    leg2.AddEntry(&g_v2inc_pseudo_data, legend_text_pseudo_data.Data(), "p");
    leg2.SetBorderSize(0);
    leg2.SetTextFont(42);
    leg2.DrawClone();

    Double_t rho = -99; // take one asdefault
    if (output_dir == "output_corr_coeff_0") rho = 0;
    else if (output_dir == "output_corr_coeff_05") rho = 0.5;
    else if (output_dir == "output_corr_coeff_1") rho = 1;
    else cout << "WARNING: could not identify correlation coefficient" << endl;
    
    TLatex la2a;
    la2a.SetTextFont(42); 
    // la2a.SetTextSize(0.042);
    if (output_dir == "output_rgam_corr_coeff_025") {
	// the setting used for the final version of the v2dir paper
	la2a.DrawLatexNDC(0.2, 0.8, Form("%s%%, #rho(v_{2,inc}) = #rho(v_{2,dec}) = 1, #rho(R_{#gamma}) = 0.25", centr.Data()));
    }
    else {
	la2a.DrawLatexNDC(0.2, 0.8, Form("%s%%, #rho = %3.1f", centr.Data(), rho));
    }
    
    TLatex la2b;
    la2b.SetTextSize(0.042);
    la2b.SetTextFont(42); 
    la2b.DrawLatexNDC(0.2, 0.75, Form("#it{p}-value = %5.3f (%3.1f#sigma)", p_value, n_sigma));

    TString fn2 = output_dir + "/v2inc_pseudo_and_measured_data_" + centr + "_" + output_dir + ".pdf";
    c2->SaveAs(fn2.Data());

    // write output
    TString fn_out(output_dir + "/v2dir_significance_" + centr + ".root");
    TFile f_out(fn_out, "recreate");
    g_v2inc_pseudo_data.SetName("g_v2inc_pseudo_data");
    g_v2inc_pseudo_data.SetTitle("v2inc pseudo data");
    g_v2inc_pseudo_data.Write();
    g_v2inc_meas.SetName("g_v2inc_meas");
    g_v2inc_meas.SetTitle("v2inc measured data points");
    g_v2inc_meas.Write();
    v2inc_pseudo_data_cov.Write("v2inc_pseudo_data_cov");
    v2inc_pseudo_data_cov_inv.Write("v2inc_pseudo_data_cov_inv");
    if (v2_inc_pd_tree) tr->Write();
    f_out.Close();
    
}

void update_sums(RooDataSet *pd, Int_t &n_evt, TVectorD &v2inc_mean_sum, TMatrixDSym &v2inc_cov_sum) {

    const Int_t n_pt_bins = v2inc_mean_sum.GetNrows();

    // loop over generated v2inc pseudo data
    for (Int_t k = 0; k < pd->numEntries(); k++) {

        TVectorD v2inc_pseudo_data(n_pt_bins);

        const RooArgSet *v2inc_row = pd->get(k);

        // cout << mvgRgam.getVal() << endl;
        for (Int_t i = 0; i < n_pt_bins; i++) {

            // It works, but is there a simpler way to retrieve the pseudo data?
            RooRealVar *v2inc_var = (RooRealVar *)v2inc_row->find(((RooArgList *)v2inc_row)->at(i)->GetName());
            v2inc_pseudo_data(i) = v2inc_var->getVal();
        }

	// fill 2d histogram for the first two pt bins for visulization 
	if (n_pt_bins > 1) hv2inc_1_2.Fill(v2inc_pseudo_data(0), v2inc_pseudo_data(1));

	// update sum needed for the calculation of the covariance matrix
        for (Int_t i = 0; i < n_pt_bins; ++i) {
            v2inc_mean_sum(i) += v2inc_pseudo_data(i);
            for (Int_t j = i; j < n_pt_bins; ++j) {
                Double_t p_ij = v2inc_pseudo_data(i) * v2inc_pseudo_data(j);
                v2inc_cov_sum(i, j) += p_ij;
                if (i != j)
                    v2inc_cov_sum(j, i) += p_ij;
            }
        }

        n_evt++;
    }
}

void get_v2inc_pseudo_data(RooDataSet *pd, Double_t* v2inc_pd, const Int_t n_pt_bins) {

    // plot only the first entry
    TGraph g_v2inc_pd(n_pt_bins);

    const RooArgSet *v2inc_row = pd->get(0);  // first row

    // cout << mvgRgam.getVal() << endl;
    for (Int_t i = 0; i < n_pt_bins; i++) {
	// It works, but is there a simpler way to retrieve the pseudo data?
	RooRealVar *v2inc_var = (RooRealVar *)v2inc_row->find(((RooArgList *)v2inc_row)->at(i)->GetName());
	v2inc_pd[i] = v2inc_var->getVal();
    }

};

RooDataSet *generate_multivariate_pseudo_data(const TVectorD &mean_vect, const TMatrixDSym &cov_matrix,
                                              const Double_t &range_min, const Double_t &range_max,
                                              const Int_t n_samples = 100) {

    const Int_t n_pt_bins = mean_vect.GetNrows();

    RooRealVar *pseudo_data_var[n_pt_bins];
    RooRealVar *mean_var[n_pt_bins];

    RooArgList pseudo_data_list;
    RooArgList mean_list;

    // prepare RooFit variables for multivariate Gaussian
    for (Int_t i = 0; i < n_pt_bins; i++) {
        char *name_pseudo_data = Form("pseudo_data_pt_bin_%d", i);
        pseudo_data_var[i] = new RooRealVar(name_pseudo_data, name_pseudo_data, 0, range_min, range_max);
        pseudo_data_list.add(*pseudo_data_var[i]);

        char *name_mean = Form("mean_vect_pt_bin_%d", i);
        mean_var[i] = new RooRealVar(name_mean, name_mean, mean_vect(i), range_min, range_max);
        mean_list.add(*mean_var[i]);
    }

    RooMultiVarGaussian mvg("mvg", "mvg", pseudo_data_list, mean_list, cov_matrix);

    RooDataSet *ds = mvg.generate(pseudo_data_list, n_samples);

    // free memory
    for (Int_t i = 0; i < n_pt_bins; i++) {
        delete pseudo_data_var[i];
        delete mean_var[i];
    }

    return ds;
};

Double_t dot_product(const TVectorD &a, const TVectorD &b) {
    Int_t n = a.GetNrows();
    Double_t dotprod = 0;
    for (Int_t i = 0; i < n; i++)
        dotprod += a(i) * b(i);
    return dotprod;
}

Double_t p_value_to_n_sigma(const Double_t &pvalue) {

    // two-sided p-value calculation

    TF1 f("f", "1-TMath::Erf(x/sqrt(2.)) - [0]", -100, 100);
    f.SetParameter(0, pvalue);

    // create wrapper function
    ROOT::Math::WrappedTF1 wfrf(f);

    // create root finder
    ROOT::Math::BrentRootFinder brf;

    // set parameters of the method
    brf.SetFunction(wfrf, 0., 50.);
    brf.Solve();
    Double_t nsigma = brf.Root();

    return nsigma;
}


Double_t cdf_standard_normal_distr(const Double_t& x) {
    return 1. - 0.5 * TMath::Erfc(x / TMath::Sqrt(2));
}


Double_t inverse_cdf_standard_normal_distr(const Double_t& x) {
    if (x < 0 || x > 1) {
	cout << "ERROR in inverse_cdf_standard_normal_distr: argument out of range" << endl;
	return -999.;
    }
    else {
	return TMath::Sqrt(2) * TMath::ErfcInverse(2. * (1. - x));
    }
}
