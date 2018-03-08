//
// quantify the level of agreement between v2,inc measurements from
// PCM and PHOS talking into account pT-correlated uncertainties
//
// Approach: check whether difference between v2,inc,PCM and v2,inc PHOS
// is consistent with zero
//
// Klaus Reygers, 18 September 2017
//

const Int_t n_bins = 16;

Double_t DotProduct(const TVectorD &a, const TVectorD &b);
Double_t chi2(const TVectorD &x, const TVectorD &mu, const TMatrixDSym &V);

void level_of_agreement_v2inc() {

    // TString centr = "00-20";
    TString centr = "20-40";

    TString fn_pcm = "data_pcm/data_pcm.root";
    TFile f_pcm(fn_pcm.Data());
    TDirectory *dir_pcm = (TDirectory *)f_pcm.Get(centr.Data());

    TString fn_phos = "data_phos/data_phos.root";
    TFile f_phos(fn_phos.Data());
    TDirectory *dir_phos = (TDirectory *)f_phos.Get(centr.Data());

    TVectorD v2_inc_meas_values_pcm = *(TVectorD *)dir_pcm->Get("v2_inc_meas_values");
    TMatrixDSym cov_v2_inc_toterr_pcm = *(TMatrixDSym *)dir_pcm->Get("cov_v2_inc_toterr");
    TMatrixDSym cov_v2_inc_staterr_pcm = *(TMatrixDSym *)dir_pcm->Get("cov_v2_inc_staterr");
    TMatrixDSym cov_v2_inc_syserr_pcm = *(TMatrixDSym *)dir_pcm->Get("cov_v2_inc_syserr");

    TVectorD v2_inc_meas_values_phos = *(TVectorD *)dir_phos->Get("v2_inc_meas_values");
    TMatrixDSym cov_v2_inc_toterr_phos = *(TMatrixDSym *)dir_phos->Get("cov_v2_inc_toterr");
    TMatrixDSym cov_v2_inc_syserr_phos = *(TMatrixDSym *)dir_phos->Get("cov_v2_inc_syserr");
    TMatrixDSym cov_v2_inc_staterr_phos = *(TMatrixDSym *)dir_phos->Get("cov_v2_inc_staterr");

    TVectorD pt = *(TVectorD *)dir_pcm->Get("pt");

    // possible cross check: no correlations, i.e., all off-diagonal elements are zero
    bool no_correlations = false; // default
    // bool no_correlations = true; // cross check
    if (no_correlations) {
        for (Int_t i = 0; i < n_bins; ++i) {
            for (Int_t j = 0; j < n_bins; ++j) {
                if (i != j) {
                    cov_v2_inc_toterr_pcm(i, j) = 0;
                    cov_v2_inc_toterr_phos(i, j) = 0;
                }
            }
        }
    }

    TVectorD v2_inc_diff = v2_inc_meas_values_phos - v2_inc_meas_values_pcm;
    TMatrixDSym cov_v2_inc_diff = cov_v2_inc_toterr_phos + cov_v2_inc_toterr_pcm;

    TVectorD x(16); // contains only zeros
    Double_t chi2val = chi2(x, v2_inc_diff, cov_v2_inc_diff);
    cout << centr << "%: ";
    cout << "chi2 = " << chi2val;
    cout << ", p-value = " << TMath::Prob(chi2val, 16) << endl << endl;

    //
    // Plot comparison as a visual cross check
    //
    TGraphErrors g_v2_inc_PCM(n_bins);
    TGraphErrors g_v2_inc_PHOS(n_bins);
    for (Int_t i = 0; i < n_bins; ++i) {
        Double_t sigma_pcm = TMath::Sqrt(cov_v2_inc_toterr_pcm(i, i));
        g_v2_inc_PCM.SetPoint(i, pt(i), v2_inc_meas_values_pcm(i));
        g_v2_inc_PCM.SetPointError(i, 0.05, sigma_pcm);

        Double_t sigma_phos = TMath::Sqrt(cov_v2_inc_toterr_phos(i, i));
        g_v2_inc_PHOS.SetPoint(i, pt(i), v2_inc_meas_values_phos(i));
        g_v2_inc_PHOS.SetPointError(i, 0.05, sigma_phos);

        // print PHOS / PCM ratio
        cout << pt(i) << " GeV/c, v2,inc,PHOS / v2,inc,PCM = " << v2_inc_meas_values_phos(i) / v2_inc_meas_values_pcm(i)
             << endl;
    }

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

    TCanvas *c2b = new TCanvas("c2b", "c2b", 700, 10, 700, 500);
    TH2F frame("frame", "inclusive photon v_{2}", 1, 0., 6, 1, -0.05, 0.30);
    frame.SetXTitle("p_{T} (GeV/#it{c})");
    frame.SetYTitle("inclusive photon v_{2}");
    frame.DrawCopy();

    TLine ly0(0., 0., 6., 0.);
    ly0.DrawClone();

    g_v2_inc_PCM.SetFillColorAlpha(kRed, 0.2);
    g_v2_inc_PHOS.SetFillColorAlpha(kBlue, 0.2);

    g_v2_inc_PCM.SetMarkerStyle(kFullCircle);
    g_v2_inc_PHOS.SetMarkerStyle(kFullCircle);

    g_v2_inc_PCM.SetLineColor(kRed);
    g_v2_inc_PHOS.SetLineColor(kBlue);

    g_v2_inc_PCM.SetMarkerColor(kRed);
    g_v2_inc_PHOS.SetMarkerColor(kBlue);

    g_v2_inc_PCM.DrawClone("pE5");
    g_v2_inc_PHOS.DrawClone("pE5");

    TLatex l2b;
    l2b.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));

    TLegend leg2b(0.55, 0.72, 0.87, 0.88, NULL, "brNDC");
    leg2b.AddEntry(&g_v2_inc_PCM, "#it{v}_{2,inc}^{PCM}", "p");
    leg2b.AddEntry(&g_v2_inc_PHOS, "#it{v}_{2,inc}^{PHOS}", "p");
    leg2b.SetBorderSize(1);
    leg2b.SetTextSize(0.032);
    leg2b.DrawClone();
}

Double_t DotProduct(const TVectorD &a, const TVectorD &b) {
    Int_t n = a.GetNrows();
    Double_t dotprod = 0;
    for (Int_t i = 0; i < n; i++)
        dotprod += a(i) * b(i);
    return dotprod;
}

Double_t chi2(const TVectorD &x, const TVectorD &mu, const TMatrixDSym &V) {
    TVectorD d = x - mu;
    TMatrixDSym V_copy = V;
    TMatrixDSym Vinv = V_copy.Invert();
    eturn Product(d, Vinv * d);
}
