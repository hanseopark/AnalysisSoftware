//
// Significance of certain v2dir hypothesis based on marginalized likelihoods
// Klaus Reygers, March 2018
//

//
// used functions
//
Double_t dot_product(const TVectorD &a, const TVectorD &b);

RooDataSet *generate_multivariate_pseudo_data(const TVectorD &mean_vect, const TMatrixDSym &cov_matrix,
                                              const Double_t &range_min, const Double_t &range_max,
                                              const Int_t n_samples);

void update_sums(RooDataSet *pd, Int_t &n_evt, TVectorD &v2inc_mean_sum, TMatrixDSym &v2inc_cov_sum);

Double_t p_value_to_n_sigma(const Double_t &pvalue);

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

    output_dir += "/";

    TString fn = Form("output/v2dir_pcm_phos_comb_%s.root", centr.Data());
    TFile f(fn.Data());

    // v2inc
    TGraphAsymmErrors g_v2_inc_comb_toterr = *(TGraphAsymmErrors *)f.Get("g_v2_inc_comb_toterr");
    TMatrixDSym cov_v2inc_all_bins = *(TMatrixDSym *)f.Get("cov_v2_inc_toterr_comb");
    TMatrixDSym cov_v2inc(n_pt_bins);
    cov_v2inc_all_bins.GetSub(pt_bin_low, pt_bin_up, pt_bin_low, pt_bin_up, cov_v2inc);
    TMatrixDSym cov_v2inc_inv = cov_v2inc;
    cov_v2inc_inv.Invert();

    // v2dec
    TGraphAsymmErrors g_v2_dec_comb = *(TGraphAsymmErrors *)f.Get("g_v2_dec_comb");
    TMatrixDSym cov_v2dec_all_bins = *(TMatrixDSym *)f.Get("cov_v2_dec_toterr");
    TMatrixDSym cov_v2dec(n_pt_bins);
    cov_v2dec_all_bins.GetSub(pt_bin_low, pt_bin_up, pt_bin_low, pt_bin_up, cov_v2dec);

    // Rgamma
    TGraphErrors g_Rgamma_toterr = *(TGraphErrors *)f.Get("g_Rgamma_toterr");
    TMatrixDSym cov_Rgam_all_bins = *(TMatrixDSym *)f.Get("cov_Rgamma_comb_toterr");
    TMatrixDSym cov_Rgam(n_pt_bins);
    cov_Rgam_all_bins.GetSub(pt_bin_low, pt_bin_up, pt_bin_low, pt_bin_up, cov_Rgam);

    // fill measured v2inc, v2dec, and Rgam into TVectorD
    TVectorD v2inc_meas(n_pt_bins);
    TVectorD v2dec_calc(n_pt_bins);
    TVectorD Rgam_meas(n_pt_bins);
    TVectorD pt(n_pt_bins);
    for (Int_t i = 0; i < n_pt_bins; i++) {
        pt(i) = g_v2_inc_comb_toterr.GetX()[i + pt_bin_low];
        v2inc_meas(i) = g_v2_inc_comb_toterr.GetY()[i + pt_bin_low];
        v2dec_calc(i) = g_v2_dec_comb.GetY()[i + pt_bin_low];
        Rgam_meas(i) = g_Rgamma_toterr.GetY()[i + pt_bin_low];
        if (Rgam_meas(i) < 1.)
            Rgam_meas(i) = 1.;
    }

    // print centrality class and selected pt range
    cout << "Centrality class " << centr.Data() << "%, selected pt range: " << pt(0) << " - " << pt(n_pt_bins - 1)
         << " GeV/c" << endl;

    // vector for v2dir hypothesis
    TVectorD v2dir_hyp(n_pt_bins);
    TString v2dir_hyp_name = "";

    // select hypothesis
    switch (hypothesis_id) {
    case 0:
        v2dir_hyp_name = "v2dir_equals_0";
        for (Int_t i = 0; i < n_pt_bins; i++)
            v2dir_hyp(i) = 0;
        break;
    case 1:
        v2dir_hyp_name = "v2dir_equals_v2dec";
        for (Int_t i = 0; i < n_pt_bins; i++)
            v2dir_hyp(i) = v2dec_calc(i);
        break;
    }

    bool ignore_Rgamma_uncertainty = false; // should be false, only true for tests
    if (ignore_Rgamma_uncertainty)
        cout << "TEST: ignore Rgamma uncertainty (NOT for final result)" << endl;

    // histogram for test statistic
    TH1F h_test_statistic("h_test_statistic", "h_test_statistic", 1000, 0., 500.);

    Int_t n_events = 0;
    TVectorD v2inc_mean_sum(n_pt_bins);
    TMatrixDSym v2inc_cov_sum(n_pt_bins);

    // calculate central value of the v2inc prediction
    TVectorD v2inc_pred_central_value(n_pt_bins);
    for (Int_t i = 0; i < n_pt_bins; i++) {
        v2inc_pred_central_value(i) = ((Rgam_meas(i) - 1.) * v2dir_hyp(i) + v2dec_calc(i)) / Rgam_meas(i);
    }

    // v2inc pseudo data
    RooDataSet *v2inc_pseudo_data_all = 0;

    // sample Rgam posterior distribution
    RooDataSet *Rgam_data_set = generate_multivariate_pseudo_data(Rgam_meas, cov_Rgam, 1., 5., n_samples);

    // cout << "TEST: Rgamma allowd to fluctuate below unity!!!" << endl;
    // RooDataSet* Rgam_data_set = generate_multivariate_pseudo_data(Rgam_meas, cov_Rgam, 0., 2., n_samples);

    // sample v2dec posterior distribution
    RooDataSet *v2dec_data_set = generate_multivariate_pseudo_data(v2dec_calc, cov_v2dec, -0.5, 0.5, n_samples);

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
            Rgam_sampled(i) = Rgam_var->getVal();

            if (!ignore_Rgamma_uncertainty) {
                // default
                v2inc_sampled(i) = ((Rgam_sampled(i) - 1.) * v2dir_hyp(i) + v2dec_sampled(i)) / Rgam_sampled(i);
            } else {
                // just a test
                // v2inc_sampled(i) = ((Rgam_meas(i) - 1.) * v2dir_hyp(i) + v2dec_sampled(i)) / Rgam_meas(i);
            }
        }

        // generate v2inc pseudo data for the sampled v2inc prediction
        RooDataSet *v2inc_pseudo_data = generate_multivariate_pseudo_data(v2inc_sampled, cov_v2inc, -0.5, 0.5, 1000);

        // fill_test_statistic(v2inc_pseudo_data, v2inc_pred_central_value, cov_v2inc_inv, h_test_statistic);
	// if (k==0) {
        //     v2inc_pseudo_data_all = (RooDataSet*) v2inc_pseudo_data->Clone();
        // }
        // else {
        //     v2inc_pseudo_data_all->append(*v2inc_pseudo_data);
        // }

        update_sums(v2inc_pseudo_data, n_events, v2inc_mean_sum, v2inc_cov_sum);

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

    TVectorD d = v2inc_meas - v2inc_mean_norm;
    Double_t chi2 = dot_product(d, v2inc_pseudo_data_cov_inv * d);
    Double_t p_value = TMath::Prob(chi2, n_pt_bins);
    Double_t n_sigma = p_value_to_n_sigma(p_value);
    cout << "chi2 = " << chi2 << endl;
    cout << "p-value: " << p_value << endl;
    cout << "n_sigma = " << n_sigma << endl;

    // TVectorD v2inc_mean_roofit = *(v2inc_pseudo_data_all->mean());
    // TMatrixDSym v2inc_pseuo_data_cov_roofit = *(v2inc_pseudo_data_all->covarianceMatrix());
    // v2inc_pseuo_data_cov_roofit.Print();
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
