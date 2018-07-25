// global variables (used in chi2 function)
const Int_t pt_bin_low = 0;
const Int_t pt_bin_up = 5;
const Int_t n_pt_bins = pt_bin_up - pt_bin_low + 1;
TVectorD pt(n_pt_bins);
TVectorD v2inc(n_pt_bins);
TMatrixDSym cov_v2inc(n_pt_bins);
TMatrixDSym cov_v2inc_inv(n_pt_bins);

Double_t dot_product(const TVectorD &a, const TVectorD &b) {
    Int_t n = a.GetNrows();
    Double_t dotprod = 0;
    for (Int_t i = 0; i < n; i++)
        dotprod += a(i) * b(i);
    return dotprod;
}


Double_t f(Double_t *x, Double_t *par) {
    Double_t p0 = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t xv = x[0];
    return p2 * xv * xv + p1 * xv + p0;
}


// chi2 function
void chi2(Int_t &npar, Double_t *gin, Double_t &chi2, Double_t *par, Int_t iflag) {

    TVectorD v2inc_model(n_pt_bins);
    for (Int_t i=0; i<n_pt_bins; ++i) {
	Double_t ptv = pt(i);
	v2inc_model(i) = f(&ptv, par);
    }

    TVectorD diff = v2inc - v2inc_model;
    TVectorD tmp = cov_v2inc_inv * diff;
    chi2 = dot_product(diff, tmp);
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


// the main function
void test_fit_v2inc() {

    // TFile* file = new TFile("output_rgam_corr_coeff_025_2018_07_04/v2dir_pcm_phos_comb_20-40.root");
    TFile* file = new TFile("output_test/v2dir_pcm_phos_comb_20-40.root");

    TMatrixDSym cov_v2inc_all_bins = *(TMatrixDSym*) file->Get("cov_v2_inc_toterr_comb");

    TGraphAsymmErrors g_v2inc_all_bins = *(TGraphAsymmErrors*) file->Get("g_v2_inc_comb_staterr");
    
    cov_v2inc_all_bins.GetSub(pt_bin_low, pt_bin_up, pt_bin_low, pt_bin_up, cov_v2inc);

    // invert covariance matrix
    cov_v2inc_inv = cov_v2inc;
    cov_v2inc_inv.Invert();
    
    for (Int_t i=0; i<n_pt_bins; ++i) {
	pt(i) = g_v2inc_all_bins.GetX()[i + pt_bin_low];
	v2inc(i) = g_v2inc_all_bins.GetY()[i + pt_bin_low];
    }

    //
    // now fit the v2inc points
    //
    
    // prepare minuit
    Int_t nPar = 3; // number of fit parameters
    TMinuit m(nPar);
    m.SetFCN(chi2);
    m.SetPrintLevel(0); // -1 quiet, 0 normal, 1 verbose

    // 1 for chi2 fit, 0.5 for negative log-likelihood fir
    // see section 1.4.1 in MINUIT manual, e.g., http://hep.fi.infn.it/minuit.pdf
    m.SetErrorDef(1);

    // parameters:
    // parameter no., name, start value, step size, range min., range max.
    // range min = range max = 0 -> no limits
    m.DefineParameter(0, "p0", 0.05, 0.01, 0, 0);
    m.DefineParameter(1, "p1", 0.1, 0.01, 0, 0);
    m.DefineParameter(2, "p2", -0.05, 0.01, 0, 0);

    // now ready for minimization step
    m.Migrad();
    m.Command("SHOW COV"); // show covariance matrix

    // get the fit results
    Double_t fmin; // best function value so far
    Double_t fedm; // estimated vertical distance to the minimum
    Double_t errdef; // the value used in SetErrorDef() before
    Int_t npari; // the number of available parameters
    Int_t nparx; // the highest parameter number
    Int_t istat; // integer status about covriance matrix
    m.mnstat(fmin, fedm, errdef, npari, nparx, istat);
    
    // draw fit
    Double_t p0, p0_err, p1, p1_err,  p2, p2_err;
    m.GetParameter(0, p0, p0_err);
    m.GetParameter(1, p1, p1_err);
    m.GetParameter(2, p2, p2_err);

    // fit statistics
    Int_t ndf = n_pt_bins - nPar;
    Double_t chi2_per_ndf = fmin / ndf;
    Double_t p_value = TMath::Prob(fmin, ndf);
    cout << "chi2_per_ndf = " << chi2_per_ndf << ", p-value = " << p_value << endl;

    
    //
    // plot data and fit
    //
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
    
    TH2F fr1("fr1", "fr1", 1, 0.8, 2.2, 1, 0.13, 0.22);
    fr1.SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fr1.SetYTitle("#it{v}_{2,inc}");
    fr1.DrawCopy();

    g_v2inc_all_bins.SetMarkerStyle(24);
    g_v2inc_all_bins.DrawClone("p");

    TF1* ff = new TF1("ff", f, 0., 5., 3);
    ff->SetParameters(p0, p1, p2);
    ff->Draw("same");

    TLatex* la = new TLatex();
    TString text = Form("chi2_per_ndf = %4.2f, p-value = %4.2f", chi2_per_ndf, p_value);
    la->SetTextSize(0.035);
    la->DrawLatexNDC(0.3, 0.2, text);

    TLatex* la2 = new TLatex();
    la2->SetTextSize(0.025);
    la2->DrawLatexNDC(0.3, 0.3, "20-40% class");

    // do simple root fit in addition
    TF1* fs = new TF1("fs", "pol2", 0.9, 2.1);
    g_v2inc_all_bins.Fit("fs", "0", "", 0.9, 2.1);
    fs->SetLineStyle(2);
    fs->SetLineColor(kRed);
    fs->Draw("same");
    Double_t chi2_simple = fs->GetChisquare();
    Double_t chi2_simple_per_ndf = chi2_simple / ndf;
    Double_t p_value_simple = TMath::Prob(chi2_simple, ndf);
    cout << "chi2_simple_per_ndf = " << chi2_simple_per_ndf;
    cout << ", p-value_simple = " << p_value_simple << endl;

    TString text_simple = Form("chi2_simple_per_ndf = %4.2f, p-value_simple = %4.2f", chi2_simple_per_ndf, p_value_simple);
    la->SetTextColor(kRed);
    la->DrawLatexNDC(0.3, 0.25, text_simple);
    
    c1->SaveAs("test_fit_v2inc.pdf");
    
}


