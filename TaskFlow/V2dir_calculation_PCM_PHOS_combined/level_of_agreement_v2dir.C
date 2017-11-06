//
// quantify the level of agreement between v2,inc measurements from
// PCM and PHOS talking into account pT-correlated uncertainties
//
// Approach: check whether difference between v2,inc,PCM and v2,inc PHOS
// is consistent with zero
//
// Klaus Reygers, 10 October 2017
//

Double_t DotProduct(const TVectorD &a, const TVectorD &b);
Double_t chi2(const TVectorD &x, const TVectorD &mu, const TMatrixDSym &V);

void level_of_agreement_v2dir() {

    // TString centr = "00-20";
    TString centr = "20-40";

    TString fn = Form("output/v2dir_pcm_phos_comb_%s.root", centr.Data());
    TFile f(fn.Data());

    TGraphAsymmErrors g_v2dir_PCM = *(TGraphAsymmErrors *)f.Get("g_v2_dir_pcm_Rpcm_toterr");
    TMatrixDSym cov_v2dir_PCM = *(TMatrixDSym *)f.Get("cov_v2_dir_pcm_Rpcm_toterr");

    TGraphAsymmErrors g_v2dir_PHOS = *(TGraphAsymmErrors *)f.Get("g_v2_dir_phos_Rphos_toterr");
    TMatrixDSym cov_v2dir_PHOS = *(TMatrixDSym *)f.Get("cov_v2_dir_phos_Rphos_toterr");

    // fill TVectorD from graph
    const Int_t n_pt_bins = 16;
    TVectorD v2dir_PCM(n_pt_bins);
    TVectorD v2dir_PHOS(n_pt_bins);

    for (Int_t i = 0; i < n_pt_bins; i++) {
        v2dir_PCM(i) = g_v2dir_PCM.GetY()[i];
        v2dir_PHOS(i) = g_v2dir_PHOS.GetY()[i];
    }

    // check: scale factor for uncertainties
    // const Double_t sf = 1.2;
    // cov_v2dir_PCM *= sf*sf;
    // cov_v2dir_PHOS *= sf*sf;

    TVectorD v2dir_diff = v2dir_PHOS - v2dir_PCM;
    TMatrixDSym cov_v2dir_diff = cov_v2dir_PHOS + cov_v2dir_PCM;

    TVectorD x(16); // contains only zeros
    Double_t chi2val = chi2(x, v2dir_diff, cov_v2dir_diff);
    cout << centr << "%: ";
    cout << "chi2 = " << chi2val;
    cout << ", p-value = " << TMath::Prob(chi2val, n_pt_bins) << endl;

    // plot v2dir as a visual check

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
    TH2F frame("frame", "direct photon v_{2}", 1, 0., 6, 1, -0.05, 0.30);
    frame.SetXTitle("p_{T} (GeV/#it{c})");
    frame.SetYTitle("direct photon v_{2}");
    frame.DrawCopy();

    TLine ly0(0., 0., 6., 0.);
    ly0.DrawClone();

    g_v2dir_PCM.SetFillColorAlpha(kRed, 0.2);
    g_v2dir_PHOS.SetFillColorAlpha(kBlue, 0.2);

    g_v2dir_PCM.SetMarkerStyle(kFullCircle);
    g_v2dir_PHOS.SetMarkerStyle(kFullCircle);

    g_v2dir_PCM.SetLineColor(kRed);
    g_v2dir_PHOS.SetLineColor(kBlue);

    g_v2dir_PCM.SetMarkerColor(kRed);
    g_v2dir_PHOS.SetMarkerColor(kBlue);

    g_v2dir_PCM.DrawClone("pE5");
    g_v2dir_PHOS.DrawClone("pE5");

    TLatex l2b;
    l2b.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));

    TLegend leg2b(0.55, 0.72, 0.87, 0.88, NULL, "brNDC");
    leg2b.AddEntry(&g_v2dir_PCM, "#it{v}_{2,dir}^{PCM}(#it{v}_{2,inc}^{PCM},#it{v}_{2,dec},#it{R}_{#gamma}^{comb})",
                   "p");
    leg2b.AddEntry(&g_v2dir_PHOS, "#it{v}_{2,dir}^{PHOS}(#it{v}_{2,inc}^{PHOS},#it{v}_{2,dec},#it{R}_{#gamma}^{comb})",
                   "p");
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
    return DotProduct(d, Vinv * d);
}
