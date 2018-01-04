Bool_t kHeavyIon = kTRUE;

PLBDrawInvMass(TString fileInputPi0 = "",TString fileInputEta = "", Int_t centralityStart, Int_t centralityEnd, Int_t binPi0, Int_t binEta,Double_t scaleFactorPi0, Double_t scaleFactorEta)
{
  PPRstyle();
//   PLBDrawInvMassPHOS();
  if (centralityStart == 100 && centralityEnd == 100) kHeavyIon =kFALSE;
  PLBDrawInvMassPCM(fileInputPi0,fileInputEta,centralityStart, centralityEnd,binPi0,binEta,scaleFactorPi0,scaleFactorEta);
  
}

void SetStyleTLatex( TLatex* text, Size_t textSize, Width_t lineWidth, Color_t textColor = 1, Bool_t kNDC = kTRUE){
	if (kNDC) {text->SetNDC();}
	text->SetTextColor(textColor);
	text->SetTextSize(textSize);
	text->SetLineWidth(lineWidth);
}

//-----------------------------------------------------------------------------
// PLBDrawInvMassPHOS()
// {
//   fPi0 = TFile::Open("PLB_LHC10de_pi0_FitResult.root");
//   fEta = TFile::Open("PLB_LHC10de_Eta_FitResult.root");
// 
//   TH1F* hPi0Real = (TH1F*)fPi0->Get("Real");
//   TH1F* hEtaReal = (TH1F*)fEta->Get("Real");
//   TH1F* hPi0Peak = (TH1F*)fPi0->Get("hp2");
//   TH1F* hEtaPeak = (TH1F*)fEta->Get("hp2");
//   TF1 * fPi0Fit  = hPi0Peak->GetFunction("gs");
//   TF1 * fEtaFit  = hEtaPeak->GetFunction("gs");
// 
//   // hEtaReal->Rebin(2);
//   // hEtaPeak->Rebin(2);
// 
//   PlotInvMass("PHOS",hPi0Real, hPi0Peak, hEtaReal, hEtaPeak, fPi0Fit, fEtaFit);
// 
// }

//-----------------------------------------------------------------------------
PLBDrawInvMassPCM(TString fileInputPi02 = "",TString fileInputEta2 = "",Int_t centralityStart, Int_t centralityEnd, Int_t binPi02=0, Int_t binEta2=0, Double_t scaleFactorPi0, Double_t scaleFactorEta)
{
  fPi0 = TFile::Open(fileInputPi02.Data());
  fEta = TFile::Open(fileInputEta2.Data());

  TH1F* hPi0Real = (TH1F*)fPi0->Get(Form("Mapping_GG_InvMass_in_Pt_Bin0%i",binPi02));
  TH1F* hEtaReal = (TH1F*)fEta->Get(Form("Mapping_GG_InvMass_in_Pt_Bin0%i",binEta2));
  TH1F* hPi0Peak = (TH1F*)fPi0->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin0%i",binPi02));
  TH1F* hEtaPeak = (TH1F*)fEta->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin0%i",binEta2));
  TF1 * fPi0Fit  = (TF1* )fPi0->Get(Form("Signal_InvMassFit_in_Pt_Bin0%i",binPi02));
  TF1 * fEtaFit  = (TF1* )fEta->Get(Form("Signal_InvMassFit_in_Pt_Bin0%i",binEta2));

  PlotInvMass("PCM",hPi0Real, hPi0Peak, hEtaReal, hEtaPeak, fPi0Fit, fEtaFit,centralityStart,centralityEnd,scaleFactorPi0,scaleFactorEta);

}

//-----------------------------------------------------------------------------

PlotInvMass(const TString det,
	    const TH1F *hPi0Real, const TH1F *hPi0Peak, 
	    const TH1F *hEtaReal, const TH1F *hEtaPeak,
	    const TF1  *fPi0Fit,  const TF1  *fEtaFit,
		 Int_t centralityStart2, Int_t centralityEnd2, Double_t scaleFactorPi0 , Double_t scaleFactorEta
  			)
{
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(3);

  hPi0Real->SetLineWidth(2);
  hEtaReal->SetLineWidth(2);
  hPi0Peak->SetLineWidth(2);
  hEtaPeak->SetLineWidth(2);
  hPi0Peak->SetMarkerColor(kRed+1);
  hEtaPeak->SetMarkerColor(kRed+1);
  hPi0Peak->SetLineColor(kRed+1);
  hEtaPeak->SetLineColor(kRed+1);
  hPi0Peak->SetLineStyle(9);
  hEtaPeak->SetLineStyle(9);
  hPi0Peak->SetMarkerStyle(20);
  hEtaPeak->SetMarkerStyle(20);
  hPi0Peak->SetMarkerSize(1.5);
  hEtaPeak->SetMarkerSize(1.5);
//   hPi0Real->SetMarkerSize(1.5);
	
  fPi0Fit->SetLineColor(kBlue);
  fEtaFit->SetLineColor(kBlue);
  fPi0Fit->SetLineWidth(2);
  fEtaFit->SetLineWidth(2);

  hPi0Real->SetLabelFont(42,"XY");
  hEtaReal->SetLabelFont(42,"XY");
  hPi0Real->SetTitleFont(42,"XY");
  hEtaReal->SetTitleFont(42,"XY");
  hPi0Real->SetLabelSize(0.04,"XY");
  hEtaReal->SetLabelSize(0.04,"XY");
  hPi0Real->SetTitleSize(0.04,"XY");
  hEtaReal->SetTitleSize(0.04,"XY");

  hEtaReal->SetNdivisions(505,"XY");

  hPi0Real->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
  hEtaReal->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");

  Double_t mBinPi0 = hPi0Real->GetBinWidth(1);
  Double_t mBinEta = hEtaReal->GetBinWidth(1);

  hPi0Real->SetYTitle(Form("Events/%.0f MeV/c^{2}",mBinPi0*1000));
  hEtaReal->SetYTitle(Form("Events/%.0f MeV/c^{2}",mBinEta*1000));

  hPi0Real->SetTitleOffset(1.2,"Y");
  hEtaReal->SetTitleOffset(1.2,"Y");
  hPi0Real->SetLabelOffset(0.01,"Y");
  hEtaReal->SetLabelOffset(0.01,"Y");
  hPi0Real->SetTitleOffset(1.0,"X");
  hEtaReal->SetTitleOffset(1.0,"X");

  hPi0Real->SetMaximum(+hPi0Real->GetMaximum()*1.1);
  if (det.Contains("PCM"))
    hEtaReal->SetMaximum(+hEtaReal->GetMaximum()*1.5);
  hPi0Real->SetMinimum(-hPi0Real->GetMaximum()*0.02);
  hEtaReal->SetMinimum(-hEtaReal->GetMaximum()*0.01);

  hPi0Real->SetNdivisions(505,"X");
  hPi0Real->SetAxisRange(0.08,0.24,"X");
  hEtaReal->SetAxisRange(0.4,0.7,"X");

  TCanvas *c1 = new TCanvas(det,"InvMass",0,0,600*1.2,600);
  c1->SetRightMargin(0.03);
  c1->SetTopMargin(0.05);
  c1->SetLeftMargin(0.09);
  c1->SetBottomMargin(0.09);
  hPi0Real->GetXaxis()->SetRangeUser(0.1,0.198);
  hPi0Real->Draw("e,hist");
  hPi0Peak->Draw("ephist same");
  if (scaleFactorPi0 !=1){
    hPi0Peak->Scale(scaleFactorPi0);
//   fEtaFit ->Draw("same");
	TH1D* histoFitPi0 = (TH1D*)fPi0Fit->GetHistogram();
	histoFitPi0->Scale(scaleFactorPi0);
	histoFitPi0->Draw("same,lc");
	TLatex *labelScalingPi0 = new TLatex(0.47,0.3,Form("Signal #times %.0f",scaleFactorPi0));
	SetStyleTLatex( labelScalingPi0, 0.05,4,hEtaPeak->GetLineColor(),kTRUE);
	labelScalingPi0->Draw();
  } else {
	fPi0Fit ->Draw("same");
  }

	TLatex *labelMethod = new TLatex(0.65,0.89,Form("PCM"));
	SetStyleTLatex( labelMethod, 0.04,4,kBlack,kTRUE);
	labelMethod->Draw();

	if (kHeavyIon){	
		TLatex *labelRange = new TLatex(0.65,0.835,Form("#pi^{0}: 0.8<p_{T}<1.0 GeV/c"));
		SetStyleTLatex( labelRange, 0.04,4,kBlack,kTRUE);
		labelRange->Draw();
	} else {
		TLatex *labelRange = new TLatex(0.65,0.835,Form("#pi^{0}: 0.8<p_{T}<1.0 GeV/c"));
		SetStyleTLatex( labelRange, 0.04,4,kBlack,kTRUE);
		labelRange->Draw();
	}
  
  
  if (kHeavyIon){
	cout << centralityStart2 << "\t" <<centralityEnd2<< endl;
	cout << Form("PbPb #sqrt{#it{s}_{_{NN}}} 2.76 TeV, %s %% " ,Form(" %i - %i", centralityStart2,centralityEnd2)) << endl;
	TLatex *labelCentrality = new TLatex(0.65,0.775,Form("PbPb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV"));
	SetStyleTLatex( labelCentrality, 0.04,4,kBlack,kTRUE);
	labelCentrality->Draw();
	TLatex *labelCentrality2 = new TLatex(0.65,0.72,Form("%s %% " ,Form("%i - %i", centralityStart2,centralityEnd2)));
	SetStyleTLatex( labelCentrality2, 0.04,4,kBlack,kTRUE);
	labelCentrality2->Draw();
  } else {
	cout << Form("pp #sqrt{#it{s}} 2.76 TeV") << endl;
	
	TLatex *labelPP = new TLatex(0.65,0.77,Form("pp #sqrt{#it{s}} = 2.76 TeV"));
	SetStyleTLatex( labelPP, 0.04,4,kBlack,kTRUE);
	labelPP->Draw();
  }
  
  if (kHeavyIon){
		c1->Print(Form("%s_InvMassPi0Only_%i-%i.pdf",det.Data(),centralityStart2,centralityEnd2));
	} else {
		c1->Print(Form("%s_InvMassPi0Only_pp_2760GeV.pdf",det.Data()));
	}
// 	delete c1;
  cout << "here" << endl;
  
  TCanvas *c2 = new TCanvas("InvMassEta","InvMassEta",0,0,600*1.2,600);
  c2->SetRightMargin(0.03);
  c2->SetTopMargin(0.05);
  c2->SetLeftMargin(0.09);
  c2->SetBottomMargin(0.09);
//   c2->Draw();
  c2->cd();
//   hEtaReal->GetXaxis()->SetRangeUser(0.45,0.65);  
  hEtaReal->Draw("ehist");
  hEtaPeak->Draw("ephist same");
  hEtaPeak->Scale(scaleFactorEta);
// //   fEtaFit ->Draw("same");
  TH1D* histoFitEta = (TH1D*)fEtaFit->GetHistogram();
  histoFitEta->Scale(scaleFactorEta);
  histoFitEta->Draw("same,lc");

  if (scaleFactorEta!=1){
  	TLatex *labelScalingEta = new TLatex(0.57,0.35,Form("Signal #times %.0f",scaleFactorEta));
	SetStyleTLatex( labelScalingEta, 0.05,4,hEtaPeak->GetLineColor(),kTRUE);
	labelScalingEta->Draw();
  } 

  
  TPaveText *txtEta = new TPaveText(0.30,0.80,0.95,0.90,"NDC");
  txtEta->SetFillColor(kWhite);
  txtEta->SetBorderSize(0);
  txtEta->SetTextAlign(14);
  txtEta->SetTextSize(0.07);

  	TLatex *labelMethod = new TLatex(0.65,0.89,Form("PCM"));
	SetStyleTLatex( labelMethod, 0.04,4,kBlack,kTRUE);
	labelMethod->Draw();

	if (kHeavyIon && centralityStart2 < 40){	
		TLatex *labelRange = new TLatex(0.65,0.835,Form("#eta: 4<p_{T}<7 GeV/c"));
		SetStyleTLatex( labelRange, 0.04,4,kBlack,kTRUE);
		labelRange->Draw();
	} else if (kHeavyIon ){	
		TLatex *labelRange = new TLatex(0.65,0.835,Form("#eta: 2<p_{T}<4 GeV/c"));
		SetStyleTLatex( labelRange, 0.04,4,kBlack,kTRUE);
		labelRange->Draw();
	} else {
		TLatex *labelRange = new TLatex(0.65,0.835,Form("#eta: 1.4<p_{T}<1.8 GeV/c"));
		SetStyleTLatex( labelRange, 0.04,4,kBlack,kTRUE);
		labelRange->Draw();
	}
  
  
  if (kHeavyIon){
	cout << centralityStart2 << "\t" <<centralityEnd2<< endl;
	cout << Form("PbPb #sqrt{#it{s}_{_{NN}}} 2.76 TeV, %s %% " ,Form(" %i - %i", centralityStart2,centralityEnd2)) << endl;
	TLatex *labelCentrality = new TLatex(0.65,0.775,Form("PbPb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV"));
	SetStyleTLatex( labelCentrality, 0.04,4,kBlack,kTRUE);
	labelCentrality->Draw();
	TLatex *labelCentrality2 = new TLatex(0.65,0.72,Form("%s %% " ,Form("%i - %i", centralityStart2,centralityEnd2)));
	SetStyleTLatex( labelCentrality2, 0.04,4,kBlack,kTRUE);
	labelCentrality2->Draw();
  } else {
	cout << Form("pp #sqrt{#it{s}} 2.76 TeV") << endl;
	
	TLatex *labelPP = new TLatex(0.65,0.77,Form("pp #sqrt{#it{s}} = 2.76 TeV"));
	SetStyleTLatex( labelPP, 0.04,4,kBlack,kTRUE);
	labelPP->Draw();
  }


//   if (scaleFactorEta != 1){
// 	TLatex *labelScalingEta = new TLatex(0.77,0.3,Form("x %.0f",scaleFactorEta));
// 	SetStyleTLatex( labelScalingEta, 0.05,4,hEtaPeak->GetLineColor(),kTRUE);
// 	labelScalingEta->Draw();
//   }

    if (kHeavyIon){
		c2->Print(Form("%s_InvMassEtaOnly_%i-%i.pdf",det.Data(),centralityStart2,centralityEnd2));
	} else {
		c2->Print(Form("%s_InvMassEtaOnly_pp_2760GeV.pdf",det.Data()));
	}

  
}

//=============================================================================
PPRstyle()
{

  //////////////////////////////////////////////////////////////////////
  //
  // ROOT style macro for the TRD TDR
  //
  //////////////////////////////////////////////////////////////////////

  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(-1);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetCanvasColor(10);

  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameBorderMode(-1);
  gStyle->SetFrameLineWidth(1.2);
  gStyle->SetFrameLineColor(1);

  gStyle->SetHistFillColor(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetHistLineColor(1);

  gStyle->SetPadColor(10);
  gStyle->SetPadBorderSize(1);
  // gStyle->SetPadBorderMode(-1);

  gStyle->SetStatColor(10);
  gStyle->SetTitleColor(kBlack,"X");
  gStyle->SetTitleColor(kBlack,"Y");

  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetLabelSize(0.04,"Z");
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelFont(42,"Z");
  gStyle->SetStatFont(42);

  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.4,"Y");

  gStyle->SetTitleFillColor(kWhite);

  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

}
