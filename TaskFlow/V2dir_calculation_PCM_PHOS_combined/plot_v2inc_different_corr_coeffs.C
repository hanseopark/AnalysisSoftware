void plot_v2inc_different_corr_coeffs() {

    TString centr = "00-20";
    // TString centr = "20-40";

    TString fn_corr_coeff_1 = "./output_corr_coeff_1/v2dir_pcm_phos_comb_" + centr + ".root";
    TString fn_corr_coeff_05 = "./output_corr_coeff_05/v2dir_pcm_phos_comb_" + centr + ".root";
    TString fn_corr_coeff_0 = "./output_corr_coeff_0/v2dir_pcm_phos_comb_" + centr + ".root";

    TFile f_corr_coeff_1(fn_corr_coeff_1.Data());
    TFile f_corr_coeff_05(fn_corr_coeff_05.Data());
    TFile f_corr_coeff_0(fn_corr_coeff_0.Data());

    TGraphAsymmErrors g_v2_inc_comb_toterr_corr_coeff_1 = *(TGraphAsymmErrors*) f_corr_coeff_1.Get("g_v2_inc_comb_toterr");
    TGraphAsymmErrors g_v2_inc_comb_toterr_corr_coeff_05 = *(TGraphAsymmErrors*) f_corr_coeff_05.Get("g_v2_inc_comb_toterr");
    TGraphAsymmErrors g_v2_inc_comb_toterr_corr_coeff_0 = *(TGraphAsymmErrors*) f_corr_coeff_0.Get("g_v2_inc_comb_toterr");

    g_v2_inc_comb_toterr_corr_coeff_1.SetLineColor(kBlack);
    g_v2_inc_comb_toterr_corr_coeff_05.SetLineColor(kGreen+4);
    g_v2_inc_comb_toterr_corr_coeff_0.SetLineColor(kOrange+7);

    g_v2_inc_comb_toterr_corr_coeff_1.SetMarkerColor(kBlack);
    g_v2_inc_comb_toterr_corr_coeff_05.SetMarkerColor(kGreen+4);
    g_v2_inc_comb_toterr_corr_coeff_0.SetMarkerColor(kOrange+7);

    g_v2_inc_comb_toterr_corr_coeff_1.SetMarkerStyle(kFullCircle);
    g_v2_inc_comb_toterr_corr_coeff_05.SetMarkerStyle(29);
    g_v2_inc_comb_toterr_corr_coeff_0.SetMarkerStyle(33);

    // displace points on pt axis
    for (Int_t i=0; i<16; i++) {
	Double_t pt = g_v2_inc_comb_toterr_corr_coeff_1.GetX()[i];
	Double_t v2inc_corr_coeff_05 = g_v2_inc_comb_toterr_corr_coeff_05.GetY()[i];
	Double_t v2inc_corr_coeff_0 = g_v2_inc_comb_toterr_corr_coeff_0.GetY()[i];

	g_v2_inc_comb_toterr_corr_coeff_05.SetPoint(i, pt-0.08, v2inc_corr_coeff_05);
	g_v2_inc_comb_toterr_corr_coeff_0.SetPoint(i, pt+0.08, v2inc_corr_coeff_0);
       	
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

    TCanvas* c1 = new TCanvas("c1", "c1", 700, 500);
    TH2F fr("fr", "fr", 1, 0., 6.0, 1, 0., 0.3);
    fr.SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fr.SetYTitle("#it{v}_{2,inc}");
    fr.DrawCopy();

    g_v2_inc_comb_toterr_corr_coeff_1.DrawClone("p");
    g_v2_inc_comb_toterr_corr_coeff_05.DrawClone("p");
    g_v2_inc_comb_toterr_corr_coeff_0.DrawClone("p");

    TString cen = centr + "%";
    TLegend leg(0.7, 0.65, 0.82, 0.87, cen.Data(), "brNDC");
    leg.AddEntry(&g_v2_inc_comb_toterr_corr_coeff_1, "#rho = 1", "p");
    leg.AddEntry(&g_v2_inc_comb_toterr_corr_coeff_05, "#rho = 0.5", "p");
    leg.AddEntry(&g_v2_inc_comb_toterr_corr_coeff_0, "#rho = 0", "p");
    leg.SetBorderSize(0);
    leg.DrawClone();
    
    TString fout = "v2inc_different_corr_coeffs_" + centr + ".pdf";
    c1->SaveAs(fout);
	
    
}
