// used function
Double_t p_value_to_n_sigma(const Double_t &pvalue);

void bayes(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
	   const Double_t &v2inc_meas, const Double_t &v2inc_err,
	   const Double_t &v2dec_meas, const Double_t &v2dec_err);

void v2_dir_central_value_and_error_from_v2dir_distribution(TH1 *h_v2dir_distr, Double_t &v2dir_central_value,
								 Double_t &v2dir_err_up, Double_t &v2dir_err_low);

void v2dir_significance_marginalized_likelihood(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
						    const Double_t &v2inc_meas, const Double_t &v2inc_err,
						    const Double_t &v2dec_meas, const Double_t &v2dec_err);

void v2dir_significance_3d_likelihood(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
				      const Double_t &v2inc_meas, const Double_t &v2inc_err,
				      const Double_t &v2dec_meas, const Double_t &v2dec_err);

void gaussian_error_propagation(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
				const Double_t &v2inc_meas, const Double_t &v2inc_err,
				const Double_t &v2dec_meas, const Double_t &v2dec_err);

void neyman_construction(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
 			 const Double_t &v2inc_meas, const Double_t &v2inc_err,
 			 const Double_t &v2dec_meas, const Double_t &v2dec_err);


void test_v2dir_significance_different_methods() {

    gRandom->SetSeed(0);

    // 20-40%, pt bin 2 (pT = 1.4 GeV)
    // const Double_t Rgam_meas = 1.06129;
    // const Double_t Rgam_meas = 1.05;
    // const Double_t Rgam_err = 0.0484595;
    // const Double_t Rgam_err = 0.;

    // const Double_t v2inc_meas = 0.181394;
    // const Double_t v2inc_err = 0.00279544;
    // const Double_t v2inc_err = 0;
    
    // const Double_t v2dec_meas = 0.182441;
    // const Double_t v2dec_err = 0.00463759;
    // const Double_t v2dec_err = 0;

    // toy values
    const Double_t Rgam_meas = 1.1;
    const Double_t Rgam_err = 0.04; // "2.5 sigma"
    const Double_t v2inc_meas = 0.1;
    const Double_t v2inc_err = v2inc_meas * 0.001;
    const Double_t v2dec_meas = 0.1;
    const Double_t v2dec_err = v2dec_meas * 0.001;
     
    // // Observed values
    // const Double_t v2inc_meas = 0.10 ;  //Inclusive photons
    // const Double_t v2dec_meas = 0.105 ; //Decay photons
    // const Double_t Rgam_meas  = 1.12 ;  //Rg

    // // Resolutions
    // const Double_t v2inc_err = 0.03*v2inc_meas ; //absolute width of Gaussian unc.
    // const Double_t v2dec_err = 0.02*v2dec_meas ; //absolute width of Gaussian unc.
    // const Double_t Rgam_err  = 0.03*Rgam_meas ; //absolute width of Gaussian unc.

    // pick a large Rgamm far away from the pole
    // const Double_t Rgam_meas = 1.3;
    // const Double_t Rgam_err = 0.001;
    // const Double_t Rgam_err = 0.01;
    // const Double_t v2inc_meas = 0.17;
    // const Double_t v2inc_err = 0.005;
    // const Double_t v2dec_meas = 0.2;
    // const Double_t v2dec_err = 0.005;

    const Double_t v2dir_hyp = 0.;

    gaussian_error_propagation(v2dir_hyp, Rgam_meas, Rgam_err,
			       v2inc_meas, v2inc_err,
    			       v2dec_meas, v2dec_err);
	
    bayes(v2dir_hyp, Rgam_meas, Rgam_err,
   	  v2inc_meas, v2inc_err,
    	  v2dec_meas, v2dec_err);

    v2dir_significance_marginalized_likelihood(v2dir_hyp, Rgam_meas, Rgam_err,
					       v2inc_meas, v2inc_err,
					       v2dec_meas, v2dec_err);

    // v2dir_significance_3d_likelihood(v2dir_hyp, Rgam_meas, Rgam_err,
    //  				     v2inc_meas, v2inc_err,
    // 				     v2dec_meas, v2dec_err);


    // neyman_construction(v2dir_hyp, Rgam_meas, Rgam_err,
    //			v2inc_meas, v2inc_err,
    //			v2dec_meas, v2dec_err);

}

void bayes(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
	   const Double_t &v2inc_meas, const Double_t &v2inc_err,
	   const Double_t &v2dec_meas, const Double_t &v2dec_err) {

    // TH1F h_v2dir_check("h_v2dir_check", "h_v2dir_check", 100, -0.5, 0.5);
    TH1F h_v2dir_check("h_v2dir_check", "h_v2dir_check", 1000, -2, 2);

    TF1 f_Rgam("f_Rgam", "gaus", 1., 3.);
    f_Rgam.SetParameters(1., Rgam_meas, Rgam_err);

    for (Int_t i=0; i<10000000; ++i) {

	Double_t Rgam_s = 0;
	if (Rgam_err != 0) Rgam_s = f_Rgam.GetRandom();
	else Rgam_s = Rgam_meas;

	Double_t v2inc_s = 0;
	if (v2inc_err != 0) v2inc_s = gRandom->Gaus(v2inc_meas, v2inc_err);
	else v2inc_s = v2inc_meas;

	Double_t v2dec_s = 0;
	if (v2dec_err != 0) v2dec_s = gRandom->Gaus(v2dec_meas, v2dec_err);
	else v2dec_s = v2dec_meas;

	Double_t v2dir_s = (Rgam_s * v2inc_s - v2dec_s) / (Rgam_s - 1.);

	//	if (TMath::Abs(v2dir_s) < 0.5) h_v2dir_check.Fill(v2dir_s);
	h_v2dir_check.Fill(v2dir_s);
	
    }

    Double_t v2dir_central_value, v2dir_err_low, v2dir_err_up;
    v2_dir_central_value_and_error_from_v2dir_distribution(&h_v2dir_check, v2dir_central_value,
							   v2dir_err_up, v2dir_err_low);

    cout << endl << "Bayesian method:" << endl;
    cout << "v2dir = " << v2dir_central_value << " + " << v2dir_err_up << " - " << v2dir_err_low << endl;
    cout << "v2dir_central_value / v2dir_err_low  = " <<  v2dir_central_value / v2dir_err_low << endl;

    Int_t bin_v2dir0 = h_v2dir_check.GetXaxis()->FindBin(0.);
    Double_t p_less_then_0 = h_v2dir_check.Integral(1, bin_v2dir0) / h_v2dir_check.Integral();
    cout << "p_less_then_0 = " << p_less_then_0 << ", 2 * p_less_then_0 = " << 2. * p_less_then_0 << endl;
    
    // draw distribution of the test statistic
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
    gROOT->SetStyle("myStyle");

    TCanvas* c2 = new TCanvas("c2");
    h_v2dir_check.DrawCopy();
    
};


void v2_dir_central_value_and_error_from_v2dir_distribution(TH1 *h_v2dir_distr, Double_t &v2dir_central_value,
                                                            Double_t &v2dir_err_up, Double_t &v2dir_err_low) {

    // determine median and errors
    const Int_t n_quant = 3;
    Double_t prob_integrals[n_quant] = {0.15865, 0.5, 0.84135};
    Double_t quantiles[n_quant];

    h_v2dir_distr->GetQuantiles(n_quant, quantiles, prob_integrals);
    v2dir_central_value = quantiles[1];
    v2dir_err_up = quantiles[2] - v2dir_central_value;
    v2dir_err_low = v2dir_central_value - quantiles[0];
}

// likelihood test statistic
void v2dir_significance_marginalized_likelihood(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
						const Double_t &v2inc_meas, const Double_t &v2inc_err,
						const Double_t &v2dec_meas, const Double_t &v2dec_err) {

    // histogram for test statistic
    TH1F h_v2inc_pseudo_data("h_v2inc_pseudo_data", "h_v2inc_pseudo_data", 10000, -0.5, 0.5);

    // central v2inc prediction
    // const Double_t v2inc_pred_central_value = ((Rgam_meas - 1.) * v2dir_hyp + v2dec_meas) / Rgam_meas;
    // cout << "v2inc_pred_central_value = " << v2inc_pred_central_value << endl;

    TF1 f_Rgam("f_Rgam", "gaus", 1., Rgam_meas + 5. * Rgam_err);
    f_Rgam.SetParameters(1., Rgam_meas, Rgam_err);

    TF1 f_v2dec("f_v2dec", "gaus", v2dec_meas - 5. * v2dec_err, v2dec_meas + 5. * v2dec_err);
    f_v2dec.SetParameters(1., v2dec_meas, v2dec_err);

    TF1 f_v2inc("f_v2inc", "gaus", 0., 0.3);
 
    for (Int_t i=0; i<5000000; ++i) {

	Double_t Rgam_sampled = 0;
	if (Rgam_err != 0) Rgam_sampled = f_Rgam.GetRandom();
	else Rgam_sampled = Rgam_meas;

	Double_t v2dec_sampled = 0;
	if (v2dec_err != 0) v2dec_sampled = f_v2dec.GetRandom();
	else v2dec_sampled = v2dec_meas;
		 
	Double_t v2inc_sampled = ((Rgam_sampled - 1.) * v2dir_hyp + v2dec_sampled) / Rgam_sampled;

	// generate toy data from v2inc_sampled
	f_v2inc.SetParameters(1., v2inc_sampled, v2inc_err);
	for (Int_t j=0; j<1; ++j) {
	    Double_t v2inc_pseudo_data = 0;
	    if (v2inc_err != 0) v2inc_pseudo_data = f_v2inc.GetRandom();
	    else v2inc_pseudo_data = v2inc_sampled;
	    h_v2inc_pseudo_data.Fill(v2inc_pseudo_data);
	}

    }

    // normalize test statistics histogram
    Double_t integral = h_v2inc_pseudo_data.Integral("width");
    h_v2inc_pseudo_data.Scale(1./integral);

    // calculate p-value and deviation in units of the standard deviation
    Int_t bin_v2inc_meas = h_v2inc_pseudo_data.GetXaxis()->FindBin(v2inc_meas);
    Double_t likelihood_v2inc_meas = h_v2inc_pseudo_data.GetBinContent(bin_v2inc_meas);

    cout << "likelihood_v2inc_meas = " << likelihood_v2inc_meas << endl;

    Double_t integral_less_than_observed_likelihood = 0;
    Double_t integral_total = 0;
        for (Int_t ib=1; ib<=h_v2inc_pseudo_data.GetNbinsX(); ++ib) {
	Double_t l = h_v2inc_pseudo_data.GetBinContent(ib);
	integral_total += l;
	if (l < likelihood_v2inc_meas) integral_less_than_observed_likelihood += l;
    }

    Double_t p_value = integral_less_than_observed_likelihood / integral_total;
    Double_t n_sigma = p_value_to_n_sigma(p_value);

    cout << endl << "Marginalized likelihood (with likelihood test statistic)" << endl;
    cout << "p-value = " << p_value << endl;
    cout << "number of standard deviations: " << n_sigma << endl;

    Double_t mean = h_v2inc_pseudo_data.GetMean();
    Double_t sigma = h_v2inc_pseudo_data.GetStdDev();
    Double_t dev = TMath::Abs(v2inc_meas - mean) / sigma;
    cout << "approx. number of standard deviations:" << dev << endl;

    const Int_t n_bins = h_v2inc_pseudo_data.GetNbinsX();
    Double_t p_v2inc_gt_v2inc_meas = h_v2inc_pseudo_data.Integral(bin_v2inc_meas, n_bins) / h_v2inc_pseudo_data.Integral(1, n_bins);
    cout << "p_v2inc_gt_v2inc_meas = " << p_v2inc_gt_v2inc_meas << endl;
    
    // draw distribution of the test statistic
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
    gROOT->SetStyle("myStyle");

    TCanvas* c1 = new TCanvas("c1");
    // c1->SetLogy();
    h_v2inc_pseudo_data.SetLineColor(kBlue);
		// h_v2inc_pseudo_data.SetXTitle("v_{2,inc}^{m}");
    // h_v2inc_pseudo_data.SetYTitle("L(v_{2,inc}^{m})");
		h_v2inc_pseudo_data.SetXTitle("v_{2,inc,m}");
    h_v2inc_pseudo_data.SetYTitle("L(v_{2,inc,m})");
    h_v2inc_pseudo_data.SetAxisRange(0.1, 0.25);
    // h_v2inc_pseudo_data.SetAxisRange(-0.1, 0.3);
    h_v2inc_pseudo_data.DrawCopy("hist");

    TLine* line_v2inc_meas = new TLine(v2inc_meas, 0, v2inc_meas, 30);
    line_v2inc_meas->SetLineWidth(2);
    line_v2inc_meas->SetLineColor(kRed);
    line_v2inc_meas->Draw();

    // highlight integration range to obtain the p-value
    Double_t v2inc_meas_lim1 = 0;
    Double_t v2inc_meas_lim2 = 0;
    Double_t L_obs = h_v2inc_pseudo_data.GetBinContent(h_v2inc_pseudo_data.FindBin(v2inc_meas));
    for (Int_t ib=2; ib<=h_v2inc_pseudo_data.GetNbinsX(); ib++) {
	Double_t L_current_bin = h_v2inc_pseudo_data.GetBinContent(ib);
	Double_t L_previous_bin = h_v2inc_pseudo_data.GetBinContent(ib-1);
	if (L_previous_bin < L_obs && L_current_bin >= L_obs) v2inc_meas_lim1 = h_v2inc_pseudo_data.GetBinCenter(ib-1);
	if (L_previous_bin >= L_obs && L_current_bin < L_obs) v2inc_meas_lim2 = h_v2inc_pseudo_data.GetBinCenter(ib);
    }
    cout << "v2inc_meas_lim2 = " << v2inc_meas_lim2 << endl;
    h_v2inc_pseudo_data.SetFillStyle(1001);
    h_v2inc_pseudo_data.SetFillColorAlpha(kRed, 0.3);
    h_v2inc_pseudo_data.SetAxisRange(0.14, v2inc_meas_lim1);
    h_v2inc_pseudo_data.DrawCopy("same hist");
    h_v2inc_pseudo_data.SetAxisRange(v2inc_meas_lim2, 0.2);
    h_v2inc_pseudo_data.DrawCopy("same hist");

}


void v2dir_significance_3d_likelihood(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
						const Double_t &v2inc_meas, const Double_t &v2inc_err,
						const Double_t &v2dec_meas, const Double_t &v2dec_err) {

    // random seed
    gRandom->SetSeed(0);

    // number of bins
    Int_t nb_v2inc_meas = 500;
    Int_t nb_v2dec_meas = 500;
    Int_t nb_Rgam_meas = 500;

    // ranges of measured values
    Double_t v2_meas_min = -0.1;
    Double_t v2_meas_max = 0.5;
    Double_t Rgam_meas_min = 0.9;
    Double_t Rgam_meas_max = 1.4;

    // ranges of true values
    Double_t v2_true_min = 0.;
    Double_t v2_true_max = 0.4;
    Double_t Rgam_true_min = 1.;
    Double_t Rgam_true_max = 1.5;

    // histogram for test statistic
    TH3F h_pseudo_data("h_pseudo_data", "pseudo data: v2inc, v2dec, Rgam",
		       nb_v2inc_meas, v2_meas_min, v2_meas_max, nb_v2dec_meas, v2_meas_min, v2_meas_max, nb_Rgam_meas, Rgam_meas_min, Rgam_meas_max);

    // for simulating the nuisance parameters
    TF1 f_Rgam_true("f_Rgam_true", "gaus", Rgam_true_min, Rgam_true_max);
    TF1 f_v2dec_true("f_v2dec_true", "gaus", v2_true_min, v2_true_max);

    f_Rgam_true.SetParameters(1., Rgam_meas, Rgam_err);
    f_v2dec_true.SetParameters(1., v2dec_meas, v2dec_err);

    // for simulating the measured values (pseudo-data)
    TF1 f_v2inc("f_v2inc", "gaus", v2_meas_min, v2_meas_max);
    TF1 f_v2dec("f_v2dec", "gaus", v2_meas_min, v2_meas_max);
    TF1 f_Rgam("f_Rgam", "gaus", Rgam_meas_min, Rgam_meas_max);

    for (Int_t i=0; i<10000000; ++i) {

	// uniform distribution of niusance parameters
	// Double_t Rgam_sampled = gRandom->Uniform(Rgam_true_min, Rgam_true_max);
	// Double_t v2dec_sampled = gRandom->Uniform(v2_true_min, v2_true_max);

	// Bayesian posterioer distribution for niusance parameters
	Double_t Rgam_sampled = f_Rgam_true.GetRandom();
	Double_t v2dec_sampled =  f_v2dec_true.GetRandom();

	Double_t v2inc_sampled = ((Rgam_sampled - 1.) * v2dir_hyp + v2dec_sampled) / Rgam_sampled;

	// generate toy data for v2inc, v2dec, and Rgam
	f_v2inc.SetParameters(1., v2inc_sampled, v2inc_err);
	f_v2dec.SetParameters(1., v2dec_sampled, v2dec_err);
	f_Rgam.SetParameters(1., Rgam_sampled, Rgam_err);

	for (Int_t j=0; j<500; ++j) {
	    Double_t v2inc_pseudo_data = f_v2inc.GetRandom();
	    Double_t v2dec_pseudo_data = f_v2dec.GetRandom();
	    Double_t Rgam_pseudo_data = f_Rgam.GetRandom();

	    h_pseudo_data.Fill(v2inc_pseudo_data, v2dec_pseudo_data, Rgam_pseudo_data);

	}

    }

    TCanvas* c2 = new TCanvas("c2");
    h_pseudo_data.SetXTitle("v2inc");
    h_pseudo_data.SetYTitle("v2dec");
    h_pseudo_data.SetZTitle("Rgam");
    h_pseudo_data.DrawCopy("box");

    // likelihood (not normalized) from observed values
    Double_t L_obs =  h_pseudo_data.GetBinContent(h_pseudo_data.FindBin(v2inc_meas, v2dec_meas, Rgam_meas));
    Double_t L_obs_up =  L_obs + TMath::Sqrt(L_obs);
    Double_t L_obs_low =  L_obs - TMath::Sqrt(L_obs);

    Double_t sum_total = 0;
    Double_t sum_smaller_than_L_obs = 0;
    Double_t sum_smaller_than_L_obs_up = 0;
    Double_t sum_smaller_than_L_obs_low = 0;

    for (Int_t ib_v2inc=1; ib_v2inc <= h_pseudo_data.GetNbinsX(); ++ib_v2inc) {
	for (Int_t ib_v2dec=1; ib_v2dec <= h_pseudo_data.GetNbinsY(); ++ib_v2dec) {
	    for (Int_t ib_Rgam=1; ib_Rgam <= h_pseudo_data.GetNbinsZ(); ++ib_Rgam) {
		Double_t L_bin = h_pseudo_data.GetBinContent(ib_v2inc, ib_v2dec, ib_Rgam);
		sum_total += L_bin;
		if (L_bin <= L_obs) sum_smaller_than_L_obs += L_bin;
		if (L_bin <= L_obs_up) sum_smaller_than_L_obs_up += L_bin;
		if (L_bin <= L_obs_low) sum_smaller_than_L_obs_low += L_bin;
	    }
	}
    }

    Double_t p_value = sum_smaller_than_L_obs / sum_total;
    Double_t p_value_up = sum_smaller_than_L_obs_up / sum_total;
    Double_t p_value_low = sum_smaller_than_L_obs_low / sum_total;

    Double_t n_sigma = p_value_to_n_sigma(p_value);

    cout << endl << "3d likelihood:" << endl;
    cout << "L_obs = " <<  L_obs << endl;
    cout << "p-value = " << p_value << " + " << p_value_up - p_value << " - " << p_value - p_value_low << endl;
    cout << "n_sigma = " << n_sigma << endl;
}


void gaussian_error_propagation(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
						const Double_t &v2inc_meas, const Double_t &v2inc_err,
						const Double_t &v2dec_meas, const Double_t &v2dec_err) {

    // central value
    Double_t v2dir = (Rgam_meas * v2inc_meas -  v2dec_meas) / (Rgam_meas - 1.);

    // uncertainty from Gaussian error propagation
    Double_t v2_dir_err_squared =
                (TMath::Power(Rgam_meas - 1., 2) *
                     (TMath::Power(v2dec_err, 2) + TMath::Power(Rgam_meas, 2) * TMath::Power(v2inc_err, 2)) +
                 TMath::Power(Rgam_err, 2) * TMath::Power(v2dec_meas - v2inc_meas, 2)) /
                TMath::Power(Rgam_meas - 1., 4);
    Double_t v2_dir_err = TMath::Sqrt(v2_dir_err_squared);
    Double_t sign_gaus = TMath::Abs(v2dir - v2dir_hyp) / v2_dir_err;

    cout << endl << "Gaussian error propagation:" << endl;
    cout << "v2dir = " << v2dir << " +/- " << v2_dir_err << endl;
    cout << "significance from Gaussian error prop. (number of standard deviations): " << sign_gaus << endl;

}


void neyman_construction(const Double_t &v2dir_hyp, const Double_t &Rgam_meas, const Double_t &Rgam_err,
			 const Double_t &v2inc_meas, const Double_t &v2inc_err,
			 const Double_t &v2dec_meas, const Double_t &v2dec_err) {

    const Int_t nRg=1000;    // Step size in Rg scan: (maxRg-1)/nRg
    // const Double_t maxRg=5.;
    // const Double_t maxRg=1.5;
    const Double_t MaxRg=2.;
    const Int_t nInc=1000;   //Step size in v2Incl scan: (0.5- (-0.5))/nInc
    const Int_t nDec=1000;   //Step size in v2Decay scan: (0.5- (-0.5))/nDec
    const Int_t nDir=100;
    const Double_t chi2cut = TMath::ChisquareQuantile(0.6827,3);
    const Double_t v2Max=0.5 ;

    Double_t v2DirMin=0.5;
    Double_t v2DirMax=-0.5;

   for(Int_t iDir=0; iDir<nDir; iDir++){

       Double_t v2DirTry = -v2Max + 2.*iDir*v2Max/nDir;
       bool v2DirAccepted = false;

       for(Int_t iRg=0; iRg<nRg; iRg++){ //Test all true values of Rg, consistent with measurement within unc.
	   Double_t RgTry=1.00001+iRg*(MaxRg-1.)/nRg;

	   for(Int_t iDec=0; iDec<nDec; iDec++){ // Test all V2Dec values, consistent with measurement

	       Double_t v2DecTry = -v2Max + 2.*iDec*v2Max/nDec;

	       Double_t v2IncTry = ((RgTry - 1.) * v2DirTry + v2DecTry) / RgTry;

	       Double_t chi2 = TMath::Power((v2inc_meas - v2IncTry) / v2inc_err, 2) +
		   TMath::Power((v2dec_meas - v2DecTry) / v2dec_err, 2) +
		   TMath::Power((Rgam_meas - RgTry) / Rgam_err, 2);

	       if(chi2 < chi2cut){
		   if(v2DirTry<v2DirMin) v2DirMin=v2DirTry;
		   if(v2DirTry>v2DirMax) v2DirMax=v2DirTry;
		   v2DirAccepted = true;
	       }

	   }
       }
       // cout << "v2Dir = " << v2DirTry;
       // if (v2DirAccepted) cout << ": accepted." << endl;
       // else cout << ": rejected." << endl;
   }

   cout << endl << "Neyman construction:" << endl;
   cout << "v2dir interval: [" << v2DirMin << "," << v2DirMax << "]" << endl;

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
