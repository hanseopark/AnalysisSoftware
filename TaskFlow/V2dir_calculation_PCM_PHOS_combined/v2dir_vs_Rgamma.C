//
// Plot v2dir vs. Rgamma to illustrate how one loses sensitivty
//
// Klaus Reygers, 22 January 2018
//

void v2_dir_central_value_and_error_from_v2dir_distribution(TH1 *h_v2dir_distr, Double_t &v2dir_central_value,
                                                            Double_t &v2dir_err_up, Double_t &v2dir_err_low);

const Int_t n_bins = 16;

// select pt bin (starting from 0)
// Int_t ibin = 2;  // 1.4 GeV
// Int_t ibin = 5;  // 2.0 GeV
Int_t ibin = 3;

void v2dir_vs_Rgamma() {

    TString cen_string = "20-40%";
    TFile *f = new TFile("./output_rgam_corr_coeff_025/v2dir_pcm_phos_comb_20-40.root");

    TGraphAsymmErrors *g_v2_inc_comb_toterr = (TGraphAsymmErrors *)f->Get("g_v2_inc_comb_toterr");
    TGraphAsymmErrors *g_v2_dec_comb = (TGraphAsymmErrors *)f->Get("g_v2_dec_comb");
    TGraphErrors *g_Rgamma_toterr = (TGraphErrors *) f->Get("g_Rgamma_toterr");

    g_v2_inc_comb_toterr->Print();

    Double_t v2inc = g_v2_inc_comb_toterr->GetY()[ibin];
    Double_t v2inc_err = g_v2_inc_comb_toterr->GetEYhigh()[ibin];
    Double_t v2dec = g_v2_dec_comb->GetY()[ibin];
    Double_t v2dec_err = g_v2_dec_comb->GetEYhigh()[ibin];
    Double_t pt = g_v2_inc_comb_toterr->GetX()[ibin];

    // get also the actual Rgamma for the given pT bin
    Double_t Rgamma_meas = g_Rgamma_toterr->GetY()[ibin];
    Double_t Rgamma_meas_err = g_Rgamma_toterr->GetEY()[ibin];
    TF1* fGaus = new TF1("fGaus", "gaus", 1., 1.3);
    fGaus->SetParameters(0.1, Rgamma_meas, Rgamma_meas_err);

    // print values
    cout << "Rgamma = " << Rgamma_meas << " +/- " << Rgamma_meas_err << endl;
    cout << "v2inc = " << v2inc << " +/- " << v2inc_err << endl;
    cout << "v2dec = " << v2dec << " +/- " << v2dec_err << endl;

    const Int_t n_Rgamma = 200;
    const Double_t Rgamma_min = 1.001;
    const Double_t delta_Rgamma = 0.001;

    TGraphErrors* g_v2dir = new TGraphErrors(n_Rgamma);

    // loop over Rgamma values
    for (Int_t i_Rgamma = 0; i_Rgamma < n_Rgamma; i_Rgamma++) {
        Double_t Rgamma = Rgamma_min + i_Rgamma * delta_Rgamma;
        Double_t v2dir = (Rgamma * v2inc - v2dec) / (Rgamma - 1);

        Double_t v2dir_err =
            1. / (Rgamma - 1.) * TMath::Sqrt(TMath::Power(Rgamma * v2inc_err, 2) + TMath::Power(v2dec_err, 2));
        g_v2dir->SetPoint(i_Rgamma, Rgamma, v2dir);
        g_v2dir->SetPointError(i_Rgamma, 0., v2dir_err);

    }

    // make the plot
    // style settings
    TStyle *myStyle = new TStyle("myStyle", "My root style");
    myStyle->SetOptStat(kFALSE);
    myStyle->SetLabelOffset(0.005, "x"); // 0.005 = root default
    myStyle->SetLabelOffset(0.005, "y"); // 0.005 = root default
    myStyle->SetTitleXSize(0.06);        // 0.04  = root default
    myStyle->SetTitleYSize(0.06);        // 0.04  = root default
    myStyle->SetTitleXOffset(0.8);
    myStyle->SetTitleYOffset(0.8);
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
    TH2F* fr = new TH2F("fr", "fr", 1, 1., 1.2, 1, 0., 0.4);
    fr->SetXTitle("#it{R}_{#gamma}");
    fr->SetYTitle("#it{v}_{2,dir}");
    fr->Draw();

    // Gaussian for measured value of Rgamma + uncertainty
    // +/- one sigma range in darker color
    fGaus->SetLineWidth(3);
    fGaus->SetLineColor(kBlue);
    fGaus->SetFillStyle(1001);
    fGaus->SetFillColorAlpha(kBlue, 0.3);
    fGaus->SetRange(Rgamma_meas - Rgamma_meas_err, Rgamma_meas + Rgamma_meas_err);
    fGaus->DrawCopy("same");
    fGaus->SetFillColorAlpha(kBlue, 0.10);
    fGaus->SetRange(1., Rgamma_meas - Rgamma_meas_err);
    fGaus->DrawCopy("same");
    fGaus->SetRange(Rgamma_meas + Rgamma_meas_err, 1.2);
    fGaus->DrawCopy("same");

    // vertical line at central value of Rgamma
    TLine *l_Rgamma = new TLine(Rgamma_meas, 0., Rgamma_meas, 0.1);
    l_Rgamma->SetLineStyle(2);
    l_Rgamma->SetLineColor(kBlue);
    l_Rgamma->Draw();

    // g_v2dir->SetMarkerStyle(20);
    g_v2dir->SetFillStyle(1001);
    g_v2dir->SetFillColorAlpha(kRed, 0.3);
    g_v2dir->Draw("3");
    g_v2dir->SetLineColor(kRed);
    g_v2dir->SetLineWidth(3);
    g_v2dir->Draw("L X");

    TString la1_string = cen_string + Form(", #it{p}_{T} = %4.1f GeV/#it{c}", pt);

    TLatex* la1 = new TLatex(1.1, 0.36, la1_string);
    la1->Draw("same");

    c1->SaveAs("v2dir_vs_Rgamma.pdf");

    //
    // plot posterior distribution for the same pt bin
    //
    TString c_v2dir_posterior_distr = Form("h_v2_dir_comb_toterr_ptbin_%d", ibin);
    // TString c_v2dir_posterior_distr = Form("h_v2_dir_comb_toterr_uncorr_ptbin_%d", ibin);
    TH1F *h_v2dir_posterior_distr = (TH1F *) f->Get(c_v2dir_posterior_distr.Data());

    // normalize posterior distribution so tht integral gives unity
    Double_t v2dir_distr_integral = h_v2dir_posterior_distr->Integral("width");
    h_v2dir_posterior_distr->Scale(1./v2dir_distr_integral);

    TCanvas* c2 = new TCanvas("c2");
    TH2F* fr2 = new TH2F("fr2", "fr2", 1, -0.1, 0.4, 1, 0., 8.5);
    fr2->SetXTitle("#it{v}_{2,dir}");
    fr2->SetYTitle("d#it{p}/d#it{v}_{2,dir}");
    fr2->Draw();

    Double_t v2dir_central_value, v2dir_err_low, v2dir_err_up;
    v2_dir_central_value_and_error_from_v2dir_distribution(h_v2dir_posterior_distr, v2dir_central_value,
                                                                v2dir_err_up, v2dir_err_low);

    h_v2dir_posterior_distr->SetLineColor(kRed);
    h_v2dir_posterior_distr->DrawCopy("same hist");
    h_v2dir_posterior_distr->SetAxisRange(v2dir_central_value - v2dir_err_low, v2dir_central_value + v2dir_err_up);
    h_v2dir_posterior_distr->SetFillStyle(1001);
    h_v2dir_posterior_distr->SetFillColorAlpha(kRed, 0.3);
    h_v2dir_posterior_distr->DrawCopy("same hist");

    TLine* l_v2dir_central_value = new TLine(v2dir_central_value, 0., v2dir_central_value, 4.);
    l_v2dir_central_value->SetLineColor(kRed);
    l_v2dir_central_value->SetLineWidth(3);
    l_v2dir_central_value->SetLineStyle(2);
    l_v2dir_central_value->Draw();

    TLatex* la2 = new TLatex(0.18, 0.83, la1_string);
    la2->SetNDC();
    la2->Draw("same");

    // cross check (probability corresponding to the +/- 2 sigma interval, 95.45% for a Gaussian)
    Double_t v2dir_bin_2sigma_min = h_v2dir_posterior_distr->FindBin(v2dir_central_value - 2. * v2dir_err_low);
    Double_t v2dir_bin_2sigma_max = h_v2dir_posterior_distr->FindBin(v2dir_central_value + 2. * v2dir_err_up);
    cout << "2 sigma interval: " << h_v2dir_posterior_distr->Integral(v2dir_bin_2sigma_min, v2dir_bin_2sigma_max, "width") << endl;

    c2->SaveAs("v2dir_posterior_distr.pdf");

    //
    // cross check: simple Monte Carlo error propagation
    //
    TH1F* h_v2dir_check = new TH1F("h_v2dir_check", "h_v2dir_check", 100, -0.5, 0.5);

    TF1 f_Rgam("f_Rgam", "gaus / x^10", 1., 5.);
    // TF1 f_Rgam("f_Rgam", "gaus", 1., 5.);
    f_Rgam.SetParameters(1., Rgamma_meas, Rgamma_meas_err);
    
    for (Int_t i=0; i<10000; ++i) {

	Double_t Rgam_s = f_Rgam.GetRandom();
	Double_t v2inc_s = gRandom->Gaus(v2inc, v2inc_err);
	Double_t v2dec_s = gRandom->Gaus(v2dec, v2dec_err);

	Double_t v2dir_s = (Rgam_s * v2inc_s - v2dec_s) / (Rgam_s - 1.);

	if (TMath::Abs(v2dir_s) < 0.5) h_v2dir_check->Fill(v2dir_s);

    }

    TCanvas* c3 = new TCanvas("c3");
    TH2F* fr3 = new TH2F("fr3", "fr3", 1, -0.1, 0.4, 1, 0., 8.5);
    fr3->SetXTitle("#it{v}_{2,dir}");
    fr3->SetYTitle("d#it{p}/d#it{v}_{2,dir}");
    fr3->Draw();

    Double_t v2dir_int = h_v2dir_check->Integral("width");
    h_v2dir_check->Scale(1./v2dir_int);
    h_v2dir_check->Draw("same");

    Double_t v2dir_check_central_value, v2dir_check_err_low, v2dir_check_err_up;
    v2_dir_central_value_and_error_from_v2dir_distribution(h_v2dir_check, v2dir_check_central_value,
                                                                v2dir_check_err_up, v2dir_check_err_low);

    h_v2dir_check->SetLineColor(kRed);
    h_v2dir_check->DrawCopy("same hist");
    h_v2dir_check->SetAxisRange(v2dir_check_central_value - v2dir_check_err_low, v2dir_check_central_value + v2dir_check_err_up);
    h_v2dir_check->SetFillStyle(1001);
    h_v2dir_check->SetFillColorAlpha(kRed, 0.3);
    h_v2dir_check->DrawCopy("same hist");

    l_v2dir_central_value->Draw();

}

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
