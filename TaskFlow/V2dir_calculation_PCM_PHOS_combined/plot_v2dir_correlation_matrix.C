

const Int_t n_bins = 16;

Double_t DotProduct(const TVectorD &a, const TVectorD &b);
Double_t chi2(const TVectorD &x, const TVectorD &mu, const TMatrixDSym &V);

void plot_v2dir_correlation_matrix() {

    TString centr = "00-20";
    // TString centr = "20-40";

    TString fn = "./output/v2dir_pcm_phos_comb_" + centr + ".root";
	
    TFile f(fn.Data());
    TMatrixDSym corr_v2_dir_comb_toterr = *(TMatrixDSym *)f.Get("corr_v2_dir_comb_toterr");
    corr_v2_dir_comb_toterr.Print();

    TH2F h_corr_v2_dir_comb_toterr("h_corr_v2_dir_comb_toterr", "v2dir cov matrix",
				   n_bins, 0.5, n_bins + 0.5, n_bins, 0.5, n_bins + 0.5);
    for (Int_t i = 0; i < n_bins; ++i) {
	for (Int_t j = 0; j < n_bins; ++j) {
	    h_corr_v2_dir_comb_toterr.SetBinContent(i+1, j+1, corr_v2_dir_comb_toterr(i, j));
	}
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
    myStyle->SetPaintTextFormat("4.2f");
    gROOT->SetStyle("myStyle");

    TCanvas* c1 = new TCanvas("c1", "c1", 500, 500);
    h_corr_v2_dir_comb_toterr.SetMarkerColor(kOrange+7);
    h_corr_v2_dir_comb_toterr.DrawCopy("col");
    h_corr_v2_dir_comb_toterr.DrawCopy("text same");

    TString fout = "corr_matrix_v2_dir_comb_" + centr + ".pdf";
    c1->SaveAs(fout);

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
