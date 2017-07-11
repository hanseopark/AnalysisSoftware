// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

//This file is not supposed to be run on outputfiles of the GammaConv-Software before the 30th Sept 2010.

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "../CommonHeaders/PlottingMeson.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "ExtractSignalPiPlPiMiPiZeroTemplate.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "THnSparse.h"

Double_t TemplateBGCorrectionFitFunctionWOPeakRegion(Double_t *x, Double_t *par)
{
	// mode:	0 // new output PCM-PCM
	//			2 // new output PCM-EMCAL
	//			3 // new output PCM-PHOS
	//			4 // new output EMCAL-EMCAL
	//			5 // new output PHOS-PHOS
	//			9 // old output PCM-PCM
	Double_t fEtaPeakInterval[2], fOmegaPeakInterval[2];
	switch(fMode)
	{
		case 0: fEtaPeakInterval[0] = 0.53; fEtaPeakInterval[1] = 0.56; fOmegaPeakInterval[0] = 0.73; fOmegaPeakInterval[1] = 0.82; break;
		
		case 2: fEtaPeakInterval[0] = 0.5; fEtaPeakInterval[1] = 0.57; fOmegaPeakInterval[0] = 0.68; fOmegaPeakInterval[1] = 0.84; break;
		case 3: fEtaPeakInterval[0] = 0.5; fEtaPeakInterval[1] = 0.57; fOmegaPeakInterval[0] = 0.68; fOmegaPeakInterval[1] = 0.84; break;
		
		case 4: fEtaPeakInterval[0] = 0.47; fEtaPeakInterval[1] = 0.58; fOmegaPeakInterval[0] = 0.64; fOmegaPeakInterval[1] = 0.86; break;
		case 5: fEtaPeakInterval[0] = 0.47; fEtaPeakInterval[1] = 0.58; fOmegaPeakInterval[0] = 0.64; fOmegaPeakInterval[1] = 0.86; break;
		
		default: fEtaPeakInterval[0] = 0.48; fEtaPeakInterval[1] = 0.57; fOmegaPeakInterval[0] = 0.68; fOmegaPeakInterval[1] = 0.83; break;
	}

    if (((x[0]>fEtaPeakInterval[0]) && (x[0]<fEtaPeakInterval[1])) || ((x[0]>fOmegaPeakInterval[0]) && (x[0]<fOmegaPeakInterval[1]))) return 0.;
    else
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t TemplateBGCorrectionFitFunction(Double_t *x, Double_t *par)
{
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t FunctionBGExclusion(Double_t *x, Double_t *par)
{
    if (x[0] > fBGFitRangeLeft[1] && x[0] < fBGFitRange[0]) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}

void GetPurityHistoFrom2DHisto(TH1D* purityHisto, TH2D* histoReconstructed, TH2D* histoTrue, Double_t* RangeInMass, Double_t* RangeInPt)
{
	Int_t LeftMassBin = histoReconstructed -> GetXaxis() -> FindBin(RangeInMass[0]);
	Int_t RightMassBin = histoReconstructed -> GetXaxis() -> FindBin(RangeInMass[1]);
	
	TH1D* True1DHisto, *Reco1DHisto;
	True1DHisto = new TH1D("True1DHisto", "True1DHisto", histoReconstructed->GetNbinsY(), RangeInPt[0], RangeInPt[1]);
	Reco1DHisto = new TH1D("Reco1DHisto", "Reco1DHisto", histoReconstructed->GetNbinsY(), RangeInPt[0], RangeInPt[1]);
	
	Double_t errorTrue, errorReco;
	
	for (Int_t iPt = 1; iPt <= histoReconstructed->GetNbinsY(); iPt++)
	{
		Reco1DHisto->SetBinContent(iPt, histoReconstructed->IntegralAndError(LeftMassBin, RightMassBin, iPt, iPt, errorReco, ""));
		Reco1DHisto->SetBinError(iPt, errorReco);
		True1DHisto->SetBinContent(iPt, histoTrue->IntegralAndError(LeftMassBin, RightMassBin, iPt, iPt, errorTrue, ""));
		True1DHisto->SetBinError(iPt, errorTrue);
	}
	
	purityHisto->Sumw2();
	purityHisto->Divide(True1DHisto, Reco1DHisto, 1., 1., "");
	
	delete Reco1DHisto;
	delete True1DHisto;
}


//FUNCTION ADDED BY ME!!!
void CreateQAHistosMC(TString cutSelection, TString suffix, TString energy, TList* ESDContainer, TList* MCContainer, TList* TrueConversionContainer, Int_t fMode)
{
	TString fdate = ReturnDateString();
	cout<<"CREATING QA HISTOS MC"<<endl;
	
	//1D histos
	TH1D* histo_NEvents = (TH1D*)ESDContainer->FindObject("NEvents");
	TH1D* histo_GoodESDTracks = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
	TH1D* histo_ESD_ConvGamma_Pt = (TH1D*)ESDContainer->FindObject("ESD_ConvGamma_Pt");
	TH1D* histo_ESD_ConvGamma_Eta = (TH1D*)ESDContainer->FindObject("ESD_ConvGamma_Eta");
	TH1D* histo_PrimaryNegPions_Pt = (TH1D*)ESDContainer->FindObject("ESD_PrimaryNegPions_Pt");
	TH1D* histo_PrimaryPosPions_Pt = (TH1D*)ESDContainer->FindObject("ESD_PrimaryPosPions_Pt");
	TH1D* histo_PrimaryNegPions_Phi = (TH1D*)ESDContainer->FindObject("ESD_PrimaryNegPions_Phi");
	TH1D* histo_PrimaryPosPions_Phi = (TH1D*)ESDContainer->FindObject("ESD_PrimaryPosPions_Phi");
	TH1D* histo_PrimaryNegPions_Eta = (TH1D*)ESDContainer->FindObject("ESD_PrimaryNegPions_Eta");
	TH1D* histo_PrimaryPosPions_Eta = (TH1D*)ESDContainer->FindObject("ESD_PrimaryPosPions_Eta");
	//2D histos
	TH2D* histo_ESD_PrimaryNegPions_ClsTPC = (TH2D*)ESDContainer->FindObject("ESD_PrimaryNegPions_ClsTPC");
	TH2D* histo_ESD_PrimaryPosPions_ClsTPC = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPosPions_ClsTPC");
	TH2D* histo_ESD_PrimaryPions_DCAxy = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPions_DCAxy");
	TH2D* histo_ESD_PrimaryPions_DCAz = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPions_DCAz");
	TH2D* histo_ESD_PrimaryPions_TPCdEdx = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPions_TPCdEdx");
	TH2D* histo_ESD_PrimaryPions_TPCdEdxSignal = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPions_TPCdEdxSignal");
	TH2D* histo_ESD_PiPlusPiNeg_InvMassPt = (TH2D*)ESDContainer->FindObject("ESD_PiPlusPiNeg_InvMassPt");
	TH2D* histo_ESD_GammaGamma_InvMass_Pt = (TH2D*)ESDContainer->FindObject("ESD_GammaGamma_InvMass_Pt");
	TH2D* histo_ESD_Mother_InvMass_Pt = (TH2D*)ESDContainer->FindObject("ESD_Mother_InvMass_Pt");
	TH2D* histo_ESD_Background_InvMass_Pt = (TH2D*)ESDContainer->FindObject("ESD_Background_InvMass_Pt");
	
	//cout<<"histo_ESDConvGamma_Pt NbinsX = "<<histo_ESDConvGamma_Pt->GetNbinsX()<<endl;
	//cout<<"histo_ESDConvGamma_Pt Content in bin 6 = "<<histo_ESDConvGamma_Pt->GetBinContent(6)<<endl;
	
	Double_t RangeX[2],RangeY[2];
	
	//plot export
	{
	//1D histos export
	RangeX[0] = 0.; RangeX[1] = 8.;
	RangeY[0] = 0.1; RangeY[1] = 1.2e8;
	PlotSingle1DHistogram(histo_NEvents,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/NEvents.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"NEvents", " ", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 110.;
	RangeY[0] = 0.1; RangeY[1] = 4e6;
	PlotSingle1DHistogram(histo_GoodESDTracks,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/GoodESDTracks.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"GoodESDTracks", "Number of tracks", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	if (fMode!=4 && fMode!=5)
	PlotSingle1DHistogram(histo_ESD_ConvGamma_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_ConvGamma_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_ConvGamma_Pt", "p_{T,#gamma}", "Events",suffix,0,1,RangeX,RangeY,"#pi^{0}#rightarrow#gamma#gamma#rightarrow e^{+}e^{-}e^{+}e^{-}");
	
	RangeX[0] = -2.0; RangeX[1] = 2.0;
	RangeY[0] = 0.1; RangeY[1] = 2.0e5;
	if (fMode!=4 && fMode!=5)
	PlotSingle1DHistogram(histo_ESD_ConvGamma_Eta,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_ConvGamma_Eta.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_ConvGamma_Eta", "#eta_{#gamma}", "Events",suffix,0,0,RangeX,RangeY,"#pi^{0}#rightarrow#gamma#gamma#rightarrow e^{+}e^{-}e^{+}e^{-}");
	
	RangeX[0] = 0.; RangeX[1] = 24.9;
	RangeY[0] = 0.1; RangeY[1] = 8.0e7;
	PlotSingle1DHistogram(histo_PrimaryNegPions_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/PrimaryNegPions_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryNegPions_Pt", "p_{T,#pi^{-}} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 24.9;
	RangeY[0] = 0.1; RangeY[1] = 8.0e7;
	PlotSingle1DHistogram(histo_PrimaryPosPions_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/PrimaryPosPions_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryPosPions_Pt", "p_{T,#pi^{+}} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 2.*TMath::Pi();
	RangeY[0] = 3.0e6; RangeY[1] = 5.0e6;
	PlotSingle1DHistogram(histo_PrimaryNegPions_Phi,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/PrimaryNegPions_Phi.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryNegPions_Phi", "#varphi_{#pi^{-}}", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 2.*TMath::Pi();
	RangeY[0] = 3.0e6; RangeY[1] = 5.0e6;
	PlotSingle1DHistogram(histo_PrimaryPosPions_Phi,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/PrimaryPosPions_Phi.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryPosPions_Phi", "#varphi_{#pi^{+}}", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = -1.5; RangeX[1] = 1.5;
	RangeY[0] = 0.; RangeY[1] = 1.e8;
	PlotSingle1DHistogram(histo_PrimaryNegPions_Eta,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/PrimaryNegPions_Eta.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryNegPions_Eta", "#eta_{#pi^{-}}", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = -1.5; RangeX[1] = 1.5;
	RangeY[0] = 0.; RangeY[1] = 1.e8;
	PlotSingle1DHistogram(histo_PrimaryPosPions_Eta,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/PrimaryPosPions_Eta.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryPosPions_Eta", "#eta_{#pi^{+}}", "Events",suffix,0,0,RangeX,RangeY,"");

	//2D histos export
	RangeX[0] = 0.; RangeX[1] = 1.;
	RangeY[0] = 0.001; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryNegPions_ClsTPC,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PrimaryNegPions_ClsTPC.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryNegPions_ClsTPC", "TPC/findable clusters", "p_{T,#pi^{-}} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#pi^{-} clusters");

	RangeX[0] = 0.; RangeX[1] = 1.;
	RangeY[0] = 0.001; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPosPions_ClsTPC,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PrimaryPosPions_ClsTPC.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPosPions_ClsTPC", "TPC/findable clusters", "p_{T,#pi^{+}} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#pi^{+} clusters");
	
	RangeX[0] = -4.; RangeX[1] = 4.;
	RangeY[0] = 0.; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPions_DCAxy,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PrimaryPions_DCAxy.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPions_DCAxy", "DCA_{xy} [cm]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");

	RangeX[0] = -4.; RangeX[1] = 4.;
	RangeY[0] = 0.; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPions_DCAz,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PrimaryPions_DCAz.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPions_DCAz", "DCA_{z} [cm]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");

	RangeX[0] = 0.01; RangeX[1] = 100.;
	RangeY[0] = -10; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPions_TPCdEdx,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PrimaryPions_TPCdEdx.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPions_TPCdEdx", "p [GeV/c]", "#sigma_{TPC,#pi} [arb. u.]",suffix,1,0,1,RangeX,RangeY,"Pion candidates");

	RangeX[0] = 0.01; RangeX[1] = 100.;
	RangeY[0] = 0.; RangeY[1] = 170.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPions_TPCdEdxSignal,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PrimaryPions_TPCdEdxSignal.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPions_TPCdEdxSignal", "p [GeV/c]", "dE/dx [arb. u.]",suffix,1,0,1,RangeX,RangeY,"Pion candidates");

	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.; RangeY[1] = 20.;
	PlotSingle2DHistogram(histo_ESD_PiPlusPiNeg_InvMassPt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PiPlusPiNeg_InvMassPt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PiPlusPiNeg_InvMassPt", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#pi^{+}#pi^{-} pair inv. mass");

	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.; RangeY[1] = 20.;
	PlotSingle2DHistogram(histo_ESD_GammaGamma_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_GammaGamma_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_GammaGamma_InvMass_Pt", "M_{#gamma#gamma} [GeV/c^{2}]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#pi^{0}#rightarrow #gamma#gamma");

	RangeX[0] = 0.4; RangeX[1] = 0.9;
	RangeY[0] = 0.; RangeY[1] = 25.;
	PlotSingle2DHistogram(histo_ESD_Mother_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_Mother_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_Mother_InvMass_Pt", "M_{#pi^{+}#pi^{-}#pi^{0}} [GeV/c^{2}]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#splitline{Signal+Background}{pPb @ 5.02TeV}");

	RangeX[0] = 0.4; RangeX[1] = 0.9;
	RangeY[0] = 0.; RangeY[1] = 25.;
	PlotSingle2DHistogram(histo_ESD_Background_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_Background_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_Background_InvMass_Pt", "M_{#pi^{+}#pi^{-}#pi^{0}} [GeV/c^{2}]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#splitline{Background}{pPb @ 5.02TeV}");
	
	//slices
	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.; RangeY[1] = 6.0e6;
	//rho
	PlotSingle1DHistogram(histo_ESD_PiPlusPiNeg_InvMassPt->ProjectionX("ESD_PiPlusPiNeg_InvMassPtSlice1",4,4,""),Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PiPlusPiNeg_InvMassPtSlice1.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PiPlusPiNeg_InvMassPtSlice1", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.3,0.4) GeV/c");
	PlotSingle1DHistogram(histo_ESD_PiPlusPiNeg_InvMassPt->ProjectionX("ESD_PiPlusPiNeg_InvMassPtSlice2",5,5,""),Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PiPlusPiNeg_InvMassPtSlice2.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PiPlusPiNeg_InvMassPtSlice2", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.4,0.5) GeV/c");
	PlotSingle1DHistogram(histo_ESD_PiPlusPiNeg_InvMassPt->ProjectionX("ESD_PiPlusPiNeg_InvMassPtSlice3",6,6,""),Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PiPlusPiNeg_InvMassPtSlice3.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PiPlusPiNeg_InvMassPtSlice3", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.5,0.6) GeV/c");
	
	RangeX[0] = 0.0; RangeX[1] = 0.45;
	RangeY[0] = 0.; RangeY[1] = 3.0e3;
	//pi0
	PlotSingle1DHistogram(histo_ESD_GammaGamma_InvMass_Pt->ProjectionX("ESD_GammaGamma_InvMassPtSlice1",7,7,""),Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_GammaGamma_InvMassPtSlice1.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_GammaGamma_InvMassPtSlice1", "M_{#gamma#gamma} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.6,0.7) GeV/c",0.08,0.145);
	PlotSingle1DHistogram(histo_ESD_GammaGamma_InvMass_Pt->ProjectionX("ESD_GammaGamma_InvMassPtSlice2",8,8,""),Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_GammaGamma_InvMassPtSlice2.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_GammaGamma_InvMassPtSlice2", "M_{#gamma#gamma} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.7,0.8) GeV/c",0.08,0.145);
	PlotSingle1DHistogram(histo_ESD_GammaGamma_InvMass_Pt->ProjectionX("ESD_GammaGamma_InvMassPtSlice3",9,9,""),Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_GammaGamma_InvMassPtSlice3.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_GammaGamma_InvMassPtSlice3", "M_{#gamma#gamma} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.8,0.9) GeV/c",0.08,0.145);
	
	Double_t SlicesPt[17] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4., 6., 8., 10.0};
	Double_t RebinFactorsPtSlice[16] = {10, 8, 8, 5, 5, 5, 4, 4, 4, 4, 4, 5, 5, 8, 8, 10};
	Double_t MesonMassRange[2] = {0.0, 1.5};
	PlotMultipleSlicesOf2DHisto(histo_ESD_PiPlusPiNeg_InvMassPt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_PiPlusPiNeg_InvMassPt_Slices.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()), "ESD_PiPlusPiNeg_InvMassPt_Slices", "Pad", MesonMassRange, fdate.Data(), "omega", 4, 5, 1, 16, SlicesPt, RebinFactorsPtSlice, "#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}", "#pi^{+}#pi^{-}",  "#pi^{+}#pi^{-}", energy);
	}

	//loading from MCContainer and TrueConversionContainer
	TH1D* MC_AllGamma_Pt = (TH1D*)MCContainer->FindObject("MC_AllGamma_Pt");
	TH1D* MC_ConvGamma_Pt = (TH1D*)MCContainer->FindObject("MC_ConvGamma_Pt");
	TH1D* MC_AllPosPions_Pt = (TH1D*)MCContainer->FindObject("MC_AllPosPions_Pt");
	TH1D* MC_AllNegPions_Pt = (TH1D*)MCContainer->FindObject("MC_AllNegPions_Pt");
	TH1D* MC_GammaFromNeutralMeson_Pt = (TH1D*)MCContainer->FindObject("MC_GammaFromNeutralMeson_Pt");
	TH1D* MC_PosPionsFromNeutralMeson_Pt = (TH1D*)MCContainer->FindObject("MC_PosPionsFromNeutralMeson_Pt");
	TH1D* MC_NegPionsFromNeutralMeson_Pt = (TH1D*)MCContainer->FindObject("MC_NegPionsFromNeutralMeson_Pt");
	TH1D* MC_Eta_Pt = (TH1D*)MCContainer->FindObject("MC_Eta_Pt");
	TH1D* MC_EtaInAcc_Pt = (TH1D*)MCContainer->FindObject("MC_EtaInAcc_Pt");
	TH1D* MC_Omega_Pt = (TH1D*)MCContainer->FindObject("MC_Omega_Pt");
	TH1D* MC_OmegaInAcc_Pt = (TH1D*)MCContainer->FindObject("MC_OmegaInAcc_Pt");
	
	TH1D* ESD_TrueConvGamma_Pt = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueConvGamma_Pt");
	TH1D* ESD_TrueConvGammaFromNeutralMeson_Pt = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueConvGammaFromNeutralMeson_Pt");
	TH1D* ESD_TruePosPion_Pt = (TH1D*)TrueConversionContainer->FindObject("ESD_TruePosPion_Pt");
	TH1D* ESD_TrueNegPion_Pt = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueNegPion_Pt");
	TH1D* ESD_TruePosPionFromNeutralMeson_Pt = (TH1D*)TrueConversionContainer->FindObject("ESD_TruePosPionFromNeutralMeson_Pt");
	TH1D* ESD_TrueNegPionFromNeutralMeson_Pt = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueNegPionFromNeutralMeson_Pt");
	TH2D* ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt");
	TH2D* ESD_TrueMotherGG_InvMass_Pt = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueMotherGG_InvMass_Pt");
	TH2D* ESD_TrueMotherGGFromEta_InvMass_Pt = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueMotherGGFromEta_InvMass_Pt");
	TH2D* ESD_TrueMotherGGFromOmega_InvMass_Pt = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueMotherGGFromOmega_InvMass_Pt");
	TH2D* ESD_TruePiPlusPiNeg_InvMassPt = (TH2D*)TrueConversionContainer->FindObject("ESD_TruePiPlusPiNeg_InvMassPt");
	TH2D* ESD_TruePiPlusPiNegFromSameMother_InvMassPt = (TH2D*)TrueConversionContainer->FindObject("ESD_TruePiPlusPiNegFromSameMother_InvMassPt");
	TH2D* ESD_TruePiPlusPiNegFromEta_InvMassPt = (TH2D*)TrueConversionContainer->FindObject("ESD_TruePiPlusPiNegFromEta_InvMassPt");
	TH2D* ESD_TruePiPlusPiNegFromOmega_InvMassPt = (TH2D*)TrueConversionContainer->FindObject("ESD_TruePiPlusPiNegFromOmega_InvMassPt");
	
	//Purity calculation:
	{
	TH1D *histoTemp, *histoPtWidths;
	Double_t binWidth = 0.;
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.; RangeY[1] = 1.2;
	Int_t NumberOfPtSlices = 44;
	histoPtWidths = new TH1D("histoPtWidths", "histoPtWidths", NumberOfPtSlices, 0., 25.);
	Double_t PtSlicesRanges[45] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
									3.5, 4.0, 4.5, 5.0, 6., 7., 8., 9., 10.,11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25.};
	for(Int_t i = 1; i <= NumberOfPtSlices; i++)
	{
		histoPtWidths -> SetBinContent(i, PtSlicesRanges[i] - PtSlicesRanges[i-1]);
		histoPtWidths -> SetBinError(i, 0.);
	}
	
	if (fMode!=4 && fMode!=5)
	{
		histoTemp = (TH1D*)ESD_TrueConvGamma_Pt->Clone("histoTemp");
		histoTemp->Sumw2();
		histoTemp->Divide(ESD_TrueConvGamma_Pt, histo_ESD_ConvGamma_Pt, 1., 1., "");
		*histoTemp = *((TH1D*)histoTemp->Rebin(NumberOfPtSlices, "histoTemp1", PtSlicesRanges));
		histoTemp -> Divide(histoPtWidths);
		histoTemp -> Scale(ESD_TrueConvGamma_Pt->GetBinWidth(1));
		PlotSingle1DHistogram(histoTemp,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/Purity_Photons.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"Purity_Photons", "p_{T}", "Purity of #gamma's",suffix,0,0,RangeX,RangeY,"");
	}
	
	histoTemp = (TH1D*)ESD_TruePosPion_Pt->Clone("histoTemp");
	histoTemp->Sumw2();
	histoTemp->Divide(ESD_TruePosPion_Pt, histo_PrimaryPosPions_Pt, 1., 1., "");
	//*histoTemp = *((TH1D*)histoTemp->Rebin(NumberOfPtSlices, "histoTemp1", PtSlicesRanges));
	//histoTemp -> Divide(histoPtWidths);
	//histoTemp -> Scale(ESD_TruePosPion_Pt->GetBinWidth(1));
	PlotSingle1DHistogram(histoTemp,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/Purity_PosPions.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"Purity_PosPions", "p_{T}", "Purity of #pi^{+}'s",suffix,0,0,RangeX,RangeY,"");
	
	histoTemp = (TH1D*)ESD_TrueNegPion_Pt->Clone("histoTemp");
	histoTemp->Sumw2();
	histoTemp->Divide(ESD_TrueNegPion_Pt, histo_PrimaryNegPions_Pt, 1., 1., "");
	*histoTemp = *((TH1D*)histoTemp->Rebin(NumberOfPtSlices, "histoTemp1", PtSlicesRanges));
	histoTemp -> Divide(histoPtWidths);
	histoTemp -> Scale(ESD_TrueNegPion_Pt->GetBinWidth(1));
	PlotSingle1DHistogram(histoTemp,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/Purity_NegPions.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"Purity_NegPions", "p_{T}", "Purity of #pi^{-}'s",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 12.;
	RangeY[0] = 0.; RangeY[1] = 0.2;
	if (fMode!=4 && fMode!=5)
	{
		histoTemp = (TH1D*)ESD_TrueConvGammaFromNeutralMeson_Pt->Clone("histoTemp");
		histoTemp->Sumw2();
		histoTemp->Divide(ESD_TrueConvGammaFromNeutralMeson_Pt, ESD_TrueConvGamma_Pt, 1., 1., "");
		*histoTemp = *((TH1D*)histoTemp->Rebin(NumberOfPtSlices, "histoTemp1", PtSlicesRanges));
		histoTemp -> Divide(histoPtWidths);
		histoTemp -> Scale(ESD_TrueConvGammaFromNeutralMeson_Pt->GetBinWidth(1));
		PlotSingle1DHistogram(histoTemp,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/Purity_HowManyPhotonsFromNeutralMeson.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"Purity_HowManyPhotonsFromNeutralMeson", "p_{T}", "(true #gamma's from NM)/(true #gamma's)",suffix,0,0,RangeX,RangeY,"");
	}
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.; RangeY[1] = 1.2;
	Double_t RangeMass[2] = {0.08, 0.145};
	delete histoTemp;
	histoTemp = new TH1D("purityHisto", "purityHisto", ESD_TrueMotherGG_InvMass_Pt->GetNbinsY(), RangeX[0], RangeX[1]);
	GetPurityHistoFrom2DHisto(histoTemp, histo_ESD_GammaGamma_InvMass_Pt, ESD_TrueMotherGG_InvMass_Pt, RangeMass, RangeX);
	binWidth = histoTemp -> GetBinWidth(1);
	*histoTemp = *((TH1D*)histoTemp->Rebin(NumberOfPtSlices, "histoTemp1", PtSlicesRanges));
	histoTemp -> Divide(histoPtWidths);
	histoTemp -> Scale(binWidth);
	PlotSingle1DHistogram(histoTemp,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/Purity_NeutralPions.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"Purity_NeutralPions", "p_{T}", "Purity of #pi^{0}'s",suffix,0,0,RangeX,RangeY,"");
	delete histoTemp;
	
	RangeX[0] = 0.; RangeX[1] = 2.0;
	RangeY[0] = 0.; RangeY[1] = 0.02;
	histoTemp = new TH1D("purityHisto", "purityHisto", ESD_TrueMotherGG_InvMass_Pt->GetNbinsY(), RangeX[0], RangeX[1]);
	GetPurityHistoFrom2DHisto(histoTemp, ESD_TrueMotherGG_InvMass_Pt, ESD_TrueMotherGGFromEta_InvMass_Pt, RangeMass, RangeX);
	binWidth = histoTemp -> GetBinWidth(1);
	*histoTemp = *((TH1D*)histoTemp->Rebin(NumberOfPtSlices, "histoTemp1", PtSlicesRanges));
	histoTemp -> Divide(histoPtWidths);
	histoTemp -> Scale(binWidth);
	PlotSingle1DHistogram(histoTemp,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/Purity_HowManyNeutralPionsFromEta.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"Purity_HowManyNeutralPionsFromEta", "p_{T}", "(true #pi^{0}'s from #eta)/(true #pi^{0}'s)",suffix,0,0,RangeX,RangeY,"");
	delete histoTemp;
	
	RangeX[0] = 0.; RangeX[1] = 7.;
	RangeY[0] = 0.; RangeY[1] = 0.12;
	histoTemp = new TH1D("purityHisto", "purityHisto", ESD_TrueMotherGG_InvMass_Pt->GetNbinsY(), RangeX[0], RangeX[1]);
	GetPurityHistoFrom2DHisto(histoTemp, ESD_TrueMotherGG_InvMass_Pt, ESD_TrueMotherGGFromOmega_InvMass_Pt, RangeMass, RangeX);
	binWidth = histoTemp -> GetBinWidth(1);
	*histoTemp = *((TH1D*)histoTemp->Rebin(NumberOfPtSlices, "histoTemp1", PtSlicesRanges));
	histoTemp -> Divide(histoPtWidths);
	histoTemp -> Scale(binWidth);
	PlotSingle1DHistogram(histoTemp,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/Purity_HowManyNeutralPionsFromOmega.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"Purity_HowManyNeutralPionsFromOmega", "p_{T}", "(true #pi^{0}'s from #omega)/(true #pi^{0}'s)",suffix,0,0,RangeX,RangeY,"");
	delete histoTemp;
	}
	
	//now comes the ploting of MC folder
	{
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e9;
	PlotSingle1DHistogram(MC_AllGamma_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_AllGamma_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_AllGamma_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	if (fMode!=4 && fMode!=5)
	PlotSingle1DHistogram(MC_ConvGamma_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_ConvGamma_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_ConvGamma_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e9;
	PlotSingle1DHistogram(MC_AllPosPions_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_AllPosPions_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_AllPosPions_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e9;
	PlotSingle1DHistogram(MC_AllNegPions_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_AllNegPions_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_AllNegPions_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e7;
	PlotSingle1DHistogram(MC_GammaFromNeutralMeson_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_GammaFromNeutralMeson_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_GammaFromNeutralMeson_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	PlotSingle1DHistogram(MC_PosPionsFromNeutralMeson_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_PosPionsFromNeutralMeson_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_PosPionsFromNeutralMeson_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	PlotSingle1DHistogram(MC_NegPionsFromNeutralMeson_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_NegPionsFromNeutralMeson_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_NegPionsFromNeutralMeson_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e7;
	PlotSingle1DHistogram(MC_Eta_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_Eta_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_Eta_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e7;
	PlotSingle1DHistogram(MC_EtaInAcc_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_EtaInAcc_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_EtaInAcc_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	PlotSingle1DHistogram(MC_Omega_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_Omega_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_Omega_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	PlotSingle1DHistogram(MC_OmegaInAcc_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/MC_OmegaInAcc_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"MC_OmegaInAcc_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	}

	//now comes ploting of the True folder
	{
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e7;
	if (fMode!=4 && fMode!=5)
	PlotSingle1DHistogram(ESD_TrueConvGamma_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TrueConvGamma_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TrueConvGamma_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e7;
	if (fMode!=4 && fMode!=5)
	PlotSingle1DHistogram(ESD_TrueConvGammaFromNeutralMeson_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TrueConvGammaFromNeutralMeson_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TrueConvGammaFromNeutralMeson_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	PlotSingle1DHistogram(ESD_TruePosPion_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TruePosPion_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TruePosPion_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	PlotSingle1DHistogram(ESD_TrueNegPion_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TrueNegPion_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TrueNegPion_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	PlotSingle1DHistogram(ESD_TruePosPionFromNeutralMeson_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TruePosPionFromNeutralMeson_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TruePosPionFromNeutralMeson_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	PlotSingle1DHistogram(ESD_TrueNegPionFromNeutralMeson_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TrueNegPionFromNeutralMeson_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TrueNegPionFromNeutralMeson_Pt", "p_{T} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	//2D
	RangeX[0] = 0.4; RangeX[1] = 0.9;
	RangeY[0] = 0.0; RangeY[1] = 25.;
	PlotSingle2DHistogram(ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt", "M_{#pi^{+}#pi^{-}#pi^{0}} [GeV/c^{2}]", "p_{T,#pi^{+}#pi^{-}#pi^{0}} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.0; RangeX[1] = 0.45;
	RangeY[0] = 0.0; RangeY[1] = 25.;
	PlotSingle2DHistogram(ESD_TrueMotherGG_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TrueMotherGG_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TrueMotherGG_InvMass_Pt", "M_{#gamma#gamma} [GeV/c^{2}]", "p_{T,#gamma#gamma} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.0; RangeX[1] = 0.45;
	RangeY[0] = 0.0; RangeY[1] = 25.;
	PlotSingle2DHistogram(ESD_TrueMotherGGFromEta_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TrueMotherGGFromEta_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TrueMotherGGFromEta_InvMass_Pt", "M_{#gamma#gamma} [GeV/c^{2}]", "p_{T,#gamma#gamma} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.0; RangeX[1] = 0.45;
	RangeY[0] = 0.0; RangeY[1] = 25.;
	PlotSingle2DHistogram(ESD_TrueMotherGGFromOmega_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TrueMotherGGFromOmega_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TrueMotherGGFromOmega_InvMass_Pt", "M_{#gamma#gamma} [GeV/c^{2}]", "p_{T,#gamma#gamma} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.0; RangeY[1] = 20.;
	PlotSingle2DHistogram(ESD_TruePiPlusPiNeg_InvMassPt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TruePiPlusPiNeg_InvMassPt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TruePiPlusPiNeg_InvMassPt", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "p_{T,#pi^{+}#pi^{-}} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.0; RangeY[1] = 20.;
	PlotSingle2DHistogram(ESD_TruePiPlusPiNegFromSameMother_InvMassPt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TruePiPlusPiNegFromSameMother_InvMassPt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TruePiPlusPiNegFromSameMother_InvMassPt", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "p_{T,#pi^{+}#pi^{-}} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.0; RangeY[1] = 20.;
	PlotSingle2DHistogram(ESD_TruePiPlusPiNegFromEta_InvMassPt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TruePiPlusPiNegFromEta_InvMassPt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TruePiPlusPiNegFromEta_InvMassPt", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "p_{T,#pi^{+}#pi^{-}} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.0; RangeY[1] = 20.;
	PlotSingle2DHistogram(ESD_TruePiPlusPiNegFromOmega_InvMassPt,Form("%s/%s/%s/ExtractSignal/Plots_MC_QA/ESD_TruePiPlusPiNegFromOmega_InvMassPt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_TruePiPlusPiNegFromOmega_InvMassPt", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "p_{T,#pi^{+}#pi^{-}} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");
	
	}
}

void CreateQAHistosData(TString cutSelection, TString suffix, TString energy, TList* ESDContainer, Int_t fMode)
{
	TString fdate = ReturnDateString();
	
	cout<<"CREATING QA HISTOS DATA"<<endl;
	
	//1D histos
	TH1D* histo_NEvents = (TH1D*)ESDContainer->FindObject("NEvents");
	TH1D* histo_GoodESDTracks = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
	TH1D* histo_ESD_ConvGamma_Pt = (TH1D*)ESDContainer->FindObject("ESD_ConvGamma_Pt");
	TH1D* histo_ESD_ConvGamma_Eta = (TH1D*)ESDContainer->FindObject("ESD_ConvGamma_Eta");
	TH1D* histo_PrimaryNegPions_Pt = (TH1D*)ESDContainer->FindObject("ESD_PrimaryNegPions_Pt");
	TH1D* histo_PrimaryPosPions_Pt = (TH1D*)ESDContainer->FindObject("ESD_PrimaryPosPions_Pt");
	TH1D* histo_PrimaryNegPions_Phi = (TH1D*)ESDContainer->FindObject("ESD_PrimaryNegPions_Phi");
	TH1D* histo_PrimaryPosPions_Phi = (TH1D*)ESDContainer->FindObject("ESD_PrimaryPosPions_Phi");
	TH1D* histo_PrimaryNegPions_Eta = (TH1D*)ESDContainer->FindObject("ESD_PrimaryNegPions_Eta");
	TH1D* histo_PrimaryPosPions_Eta = (TH1D*)ESDContainer->FindObject("ESD_PrimaryPosPions_Eta");
	//2D histos
	TH2D* histo_ESD_PrimaryNegPions_ClsTPC = (TH2D*)ESDContainer->FindObject("ESD_PrimaryNegPions_ClsTPC");
	TH2D* histo_ESD_PrimaryPosPions_ClsTPC = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPosPions_ClsTPC");
	TH2D* histo_ESD_PrimaryPions_DCAxy = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPions_DCAxy");
	TH2D* histo_ESD_PrimaryPions_DCAz = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPions_DCAz");
	TH2D* histo_ESD_PrimaryPions_TPCdEdx = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPions_TPCdEdx");
	TH2D* histo_ESD_PrimaryPions_TPCdEdxSignal = (TH2D*)ESDContainer->FindObject("ESD_PrimaryPions_TPCdEdxSignal");
	TH2D* histo_ESD_PiPlusPiNeg_InvMassPt = (TH2D*)ESDContainer->FindObject("ESD_PiPlusPiNeg_InvMassPt");
	TH2D* histo_ESD_GammaGamma_InvMass_Pt = (TH2D*)ESDContainer->FindObject("ESD_GammaGamma_InvMass_Pt");
	TH2D* histo_ESD_Mother_InvMass_Pt = (TH2D*)ESDContainer->FindObject("ESD_Mother_InvMass_Pt");
	TH2D* histo_ESD_Background_InvMass_Pt = (TH2D*)ESDContainer->FindObject("ESD_Background_InvMass_Pt");
	
	//cout<<"histo_ESDConvGamma_Pt NbinsX = "<<histo_ESDConvGamma_Pt->GetNbinsX()<<endl;
	//cout<<"histo_ESDConvGamma_Pt Content in bin 6 = "<<histo_ESDConvGamma_Pt->GetBinContent(6)<<endl;
	
	Double_t RangeX[2],RangeY[2];
	
	//1D histos export
	RangeX[0] = 0.; RangeX[1] = 8.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	PlotSingle1DHistogram(histo_NEvents,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/NEvents.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"NEvents", " ", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 110.;
	RangeY[0] = 0.1; RangeY[1] = 4e6;
	PlotSingle1DHistogram(histo_GoodESDTracks,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/GoodESDTracks.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"GoodESDTracks", "Number of tracks", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 25.;
	RangeY[0] = 0.1; RangeY[1] = 1e8;
	if (fMode!=4 && fMode!=5)
	PlotSingle1DHistogram(histo_ESD_ConvGamma_Pt,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_ConvGamma_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_ConvGamma_Pt", "p_{T,#gamma}", "Events",suffix,0,1,RangeX,RangeY,"#pi^{0}#rightarrow#gamma#gamma#rightarrow e^{+}e^{-}e^{+}e^{-}");
	
	RangeX[0] = -2.0; RangeX[1] = 2.0;
	RangeY[0] = 0.1; RangeY[1] = 2.0e5;
	if (fMode!=4 && fMode!=5)
	PlotSingle1DHistogram(histo_ESD_ConvGamma_Eta,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_ConvGamma_Eta.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_ConvGamma_Eta", "#eta_{#gamma}", "Events",suffix,0,0,RangeX,RangeY,"#pi^{0}#rightarrow#gamma#gamma#rightarrow e^{+}e^{-}e^{+}e^{-}");
	
	RangeX[0] = 0.; RangeX[1] = 24.9;
	RangeY[0] = 0.1; RangeY[1] = 8.0e7;
	PlotSingle1DHistogram(histo_PrimaryNegPions_Pt,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/PrimaryNegPions_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryNegPions_Pt", "p_{T,#pi^{-}} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 24.9;
	RangeY[0] = 0.1; RangeY[1] = 8.0e7;
	PlotSingle1DHistogram(histo_PrimaryPosPions_Pt,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/PrimaryPosPions_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryPosPions_Pt", "p_{T,#pi^{+}} [GeV/c]", "Events",suffix,0,1,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 2.*TMath::Pi();
	RangeY[0] = 3.0e6; RangeY[1] = 4.5e6;
	PlotSingle1DHistogram(histo_PrimaryNegPions_Phi,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/PrimaryNegPions_Phi.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryNegPions_Phi", "#varphi_{#pi^{-}}", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = 0.; RangeX[1] = 2.*TMath::Pi();
	RangeY[0] = 3.0e6; RangeY[1] = 4.5e6;
	PlotSingle1DHistogram(histo_PrimaryPosPions_Phi,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/PrimaryPosPions_Phi.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryPosPions_Phi", "#varphi_{#pi^{+}}", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = -1.5; RangeX[1] = 1.5;
	RangeY[0] = 0.; RangeY[1] = 1.3e6;
	PlotSingle1DHistogram(histo_PrimaryNegPions_Eta,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/PrimaryNegPions_Eta.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryNegPions_Eta", "#eta_{#pi^{-}}", "Events",suffix,0,0,RangeX,RangeY,"");
	
	RangeX[0] = -1.5; RangeX[1] = 1.5;
	RangeY[0] = 0.; RangeY[1] = 1.3e6;
	PlotSingle1DHistogram(histo_PrimaryPosPions_Eta,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/PrimaryPosPions_Eta.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"PrimaryPosPions_Eta", "#eta_{#pi^{+}}", "Events",suffix,0,0,RangeX,RangeY,"");

	//2D histos export
	RangeX[0] = 0.; RangeX[1] = 1.;
	RangeY[0] = 0.001; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryNegPions_ClsTPC,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PrimaryNegPions_ClsTPC.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryNegPions_ClsTPC", "TPC/findable clusters", "p_{T,#pi^{-}} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#pi^{-} clusters");

	RangeX[0] = 0.; RangeX[1] = 1.;
	RangeY[0] = 0.001; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPosPions_ClsTPC,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PrimaryPosPions_ClsTPC.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPosPions_ClsTPC", "TPC/findable clusters", "p_{T,#pi^{+}} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#pi^{+} clusters");
	
	RangeX[0] = -4.; RangeX[1] = 4.;
	RangeY[0] = 0.; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPions_DCAxy,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PrimaryPions_DCAxy.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPions_DCAxy", "DCA_{xy} [cm]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");

	RangeX[0] = -4.; RangeX[1] = 4.;
	RangeY[0] = 0.; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPions_DCAz,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PrimaryPions_DCAz.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPions_DCAz", "DCA_{z} [cm]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"");

	RangeX[0] = 0.01; RangeX[1] = 100.;
	RangeY[0] = -10; RangeY[1] = 10.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPions_TPCdEdx,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PrimaryPions_TPCdEdx.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPions_TPCdEdx", "p [GeV/c]", "#sigma_{TPC,#pi} [arb. u.]",suffix,1,0,1,RangeX,RangeY,"Pion candidates");

	RangeX[0] = 0.01; RangeX[1] = 100.;
	RangeY[0] = 0.; RangeY[1] = 170.;
	PlotSingle2DHistogram(histo_ESD_PrimaryPions_TPCdEdxSignal,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PrimaryPions_TPCdEdxSignal.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PrimaryPions_TPCdEdxSignal", "p [GeV/c]", "dE/dx [arb. u.]",suffix,1,0,1,RangeX,RangeY,"Pion candidates");

	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.; RangeY[1] = 20.;
	PlotSingle2DHistogram(histo_ESD_PiPlusPiNeg_InvMassPt,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PiPlusPiNeg_InvMassPt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PiPlusPiNeg_InvMassPt", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#pi^{+}#pi^{-} pair inv. mass");

	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.; RangeY[1] = 20.;
	PlotSingle2DHistogram(histo_ESD_GammaGamma_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_GammaGamma_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_GammaGamma_InvMass_Pt", "M_{#gamma#gamma} [GeV/c^{2}]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#pi^{0}#rightarrow #gamma#gamma",0.1,0.145);

	RangeX[0] = 0.4; RangeX[1] = 0.9;
	RangeY[0] = 0.; RangeY[1] = 25.;
	PlotSingle2DHistogram(histo_ESD_Mother_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_Mother_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_Mother_InvMass_Pt", "M_{#pi^{+}#pi^{-}#pi^{0}} [GeV/c^{2}]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#splitline{Signal+Background}{pPb @ 5.02TeV}");

	RangeX[0] = 0.4; RangeX[1] = 0.9;
	RangeY[0] = 0.; RangeY[1] = 25.;
	PlotSingle2DHistogram(histo_ESD_Background_InvMass_Pt,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_Background_InvMass_Pt.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_Background_InvMass_Pt", "M_{#pi^{+}#pi^{-}#pi^{0}} [GeV/c^{2}]", "p_{T} [GeV/c]",suffix,0,0,1,RangeX,RangeY,"#splitline{Background}{pPb @ 5.02TeV}");
	
	//slices
	RangeX[0] = 0.0; RangeX[1] = 2.;
	RangeY[0] = 0.; RangeY[1] = 4.0e4;
	//rho
	PlotSingle1DHistogram(histo_ESD_PiPlusPiNeg_InvMassPt->ProjectionX("ESD_PiPlusPiNeg_InvMassPtSlice1",4,4,""),Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PiPlusPiNeg_InvMassPtSlice1.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PiPlusPiNeg_InvMassPtSlice1", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.3,0.4) GeV/c");
	PlotSingle1DHistogram(histo_ESD_PiPlusPiNeg_InvMassPt->ProjectionX("ESD_PiPlusPiNeg_InvMassPtSlice2",5,5,""),Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PiPlusPiNeg_InvMassPtSlice2.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PiPlusPiNeg_InvMassPtSlice2", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.4,0.5) GeV/c");
	PlotSingle1DHistogram(histo_ESD_PiPlusPiNeg_InvMassPt->ProjectionX("ESD_PiPlusPiNeg_InvMassPtSlice3",6,6,""),Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PiPlusPiNeg_InvMassPtSlice3.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_PiPlusPiNeg_InvMassPtSlice3", "M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.5,0.6) GeV/c");
	
	RangeX[0] = 0.0; RangeX[1] = 0.45;
	RangeY[0] = 0.; RangeY[1] = 3.0e3;
	//pi0
	PlotSingle1DHistogram(histo_ESD_GammaGamma_InvMass_Pt->ProjectionX("ESD_GammaGamma_InvMassPtSlice1",7,7,""),Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_GammaGamma_InvMassPtSlice1.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_GammaGamma_InvMassPtSlice1", "M_{#gamma#gamma} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.6,0.7) GeV/c",0.08,0.145);
	PlotSingle1DHistogram(histo_ESD_GammaGamma_InvMass_Pt->ProjectionX("ESD_GammaGamma_InvMassPtSlice2",8,8,""),Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_GammaGamma_InvMassPtSlice2.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_GammaGamma_InvMassPtSlice2", "M_{#gamma#gamma} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.7,0.8) GeV/c",0.08,0.145);
	PlotSingle1DHistogram(histo_ESD_GammaGamma_InvMass_Pt->ProjectionX("ESD_GammaGamma_InvMassPtSlice3",9,9,""),Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_GammaGamma_InvMassPtSlice3.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()),"ESD_GammaGamma_InvMassPtSlice3", "M_{#gamma#gamma} [GeV/c^{2}]", "Events",suffix,0,0,RangeX,RangeY,"p_{T}: (0.8,0.9) GeV/c",0.08,0.145);
	
	Double_t SlicesPt[17] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4., 6., 8., 10.0};
	Double_t RebinFactorsPtSlice[16] = {10, 8, 8, 5, 5, 5, 4, 4, 4, 4, 4, 5, 5, 8, 8, 10};
	Double_t MesonMassRange[2] = {0.0, 1.5};
	PlotMultipleSlicesOf2DHisto(histo_ESD_PiPlusPiNeg_InvMassPt,Form("%s/%s/%s/ExtractSignal/Plots_Data_QA/ESD_PiPlusPiNeg_InvMassPt_Slices.%s",cutSelection.Data(),energy.Data(),suffix.Data(),suffix.Data()), "ESD_PiPlusPiNeg_InvMassPt_Slices", "Pad", MesonMassRange, fdate.Data(), "omega", 4, 5, 1, 16, SlicesPt, RebinFactorsPtSlice, "#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}", "#pi^{+}#pi^{-}",  "#pi^{+}#pi^{-}", energy);
	
}

void ProcessEMLeftRight(TH1D* fFgr, TH1D* fBck, Double_t * fBGFitRangeEMLeft, Double_t * fBGFitRangeEMRight)
{
   //cout<<"{"<<endl<<"ProcessEM started"<<endl;
   for (Int_t binx= 0; binx < fFgr->GetNbinsX()+1; binx++){
      if(fFgr->GetBinContent(binx) == 0){
         fFgr->SetBinError(binx,1.);
         fFgr->SetBinContent(binx,0.);
      }
   }
   
   fBckNorm = (TH1D*)fBck->Clone("fBckNorm");
   fFgr->Sumw2();
   fBck->Sumw2();
   fBckNorm->Sumw2();
   
   Double_t 	norm = 1;
   
   Double_t 	r= fFgr->Integral(fFgr->GetXaxis()->FindBin(fBGFitRangeEMLeft[0]),fFgr->GetXaxis()->FindBin(fBGFitRangeEMLeft[1]));
   Double_t 	b= fBck->Integral(fBck->GetXaxis()->FindBin(fBGFitRangeEMLeft[0]),fBck->GetXaxis()->FindBin(fBGFitRangeEMLeft[1]));
   
   r+= fFgr->Integral(fFgr->GetXaxis()->FindBin(fBGFitRangeEMRight[0]),fFgr->GetXaxis()->FindBin(fBGFitRangeEMRight[1]));
   b+= fBck->Integral(fBck->GetXaxis()->FindBin(fBGFitRangeEMRight[0]),fBck->GetXaxis()->FindBin(fBGFitRangeEMRight[1]));
   
   if(b != 0) norm = r/b;
   
   fBckNorm->Scale(norm);
   //    cout<<"r="<<r<<" b="<<b<<" r/b="<<r/b<< " " << endl;
   
   Int_t numberOfZeros = 0;
   for (Int_t i = 1; i < fBckNorm->GetNbinsX()+1; i++){
      if (fBckNorm->GetBinContent(i) == 0){
         numberOfZeros++;
         if (norm > 1.){
            fBckNorm->SetBinError(i,1.);
            fBckNorm->SetBinContent(i,0.);
         }
      }
   }
//    cout << "counted " << numberOfZeros << " in the normalized BG : "<< (Double_t)numberOfZeros/fBck->GetNbinsX() << endl;
   
   fSignal = (TH1D*)fFgr->Clone("fSignal");
   fSignal->Sumw2();
 
   if ((Double_t)numberOfZeros/fBck->GetNbinsX()< 0.25) fSignal->Add(fBckNorm,-1.);
   //cout<<"ProcessEM finished"<<endl<<"}"<<endl;
}
                                                               
// Main Function
void ExtractSignalPiPlPiMiPiZeroTemplate(TString meson="", TString subDir="", TString file="", TString fileMC="", TString cutSelection="", TString Suffix="", TString optionMC="", TString optionEnergy="", TString optionCrystalBall="", TString directphotonPlots="", TString optionUseMinBiasEff="", TString optionPeriod="", TString optionAdvancedMesonQA="",Int_t numberOfBins = 30, Bool_t addSig = kFALSE, Bool_t makeQAplots = kTRUE){
	gROOT->Reset();

	// mode:	0 // new output PCM-PCM
	//			2 // new output PCM-EMCAL
	//			3 // new output PCM-PHOS
	//			4 // new output EMCAL-EMCAL
	//			5 // new output PHOS-PHOS
	//			9 // old output PCM-PCM
	
	if(directphotonPlots){}
   
	if(optionAdvancedMesonQA.Contains("AdvancedMesonQA")){fAdvancedMesonQA = kTRUE;}
	
	TString fTypeCutSelection="";
	TString fEventCutSelection="";
	TString fGammaCutSelection="";
	TString fClusterCutSelection="";
	TString fPionCutSelection="";
	TString fNeutralPionCutSelection="";
	TString fMesonCutSelection="";
	
	fCutSelection = cutSelection;
	TString fCutSelectionRead = cutSelection;
	
	//Int_t fMode = -1;
	fMode = ReturnSeparatedCutNumberPiPlPiMiPiZero(cutSelection,fTypeCutSelection, fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fPionCutSelection, fNeutralPionCutSelection, fMesonCutSelection);
	Int_t mode = fMode;
	
	cout << "\t MODE = " << mode << endl;
	
	TString fEventCutSelectionRead = fEventCutSelection.Data();
	TString fGammaCutSelectionRead = fGammaCutSelection.Data();
	TString fMesonCutSelectionRead = fMesonCutSelection.Data();
	if (addSig) {
		cout << "running added Signal" << endl;
		cout << fEventCutSelection.Data() << endl;
		fEventCutSelection.Replace(6,1,"2");
		cout << fEventCutSelection.Data() << endl;
		fEventCutSelectionRead = fEventCutSelection;
		fGammaCutSelectionRead = fGammaCutSelection;
		fMesonCutSelectionRead = fMesonCutSelection;
		if (mode==9)fCutSelectionRead = Form("%s%s_%s", fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
		else if (mode==0) fCutSelectionRead = Form("%s_%s_%s",fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
		cout << fCutSelectionRead.Data() << endl;
	}
	if(optionUseMinBiasEff.CompareTo("MinBiasEffOnly")==0 && optionMC.CompareTo("kTRUE") == 0){
		cout << "calculating MinBias Eff" << endl;
		cout << fEventCutSelection.Data() << endl;
		fEventCutSelection.Replace(1,2,"00");
		fEventCutSelectionRead = fEventCutSelection;
		fGammaCutSelectionRead = fGammaCutSelection;
		fMesonCutSelectionRead = fMesonCutSelection;      
		cout << fGammaCutSelection.Data() << endl;
		if (mode==9)fCutSelectionRead = Form("%s%s_%s", fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
		else if (mode==0)fCutSelectionRead = Form("%s_%s_%s", fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
		cout << fCutSelectionRead.Data() << endl;
	}
	
	
	StyleSettingsThesis();
	SetPlotStyle();
		
	fEnergyFlag = optionEnergy;
	fPrefix=meson;

	fPeriodFlag = optionPeriod;
	fdirectphoton = directphotonPlots;

	TString outputDir = Form("%s/%s/%s/ExtractSignal",cutSelection.Data(),optionEnergy.Data(),Suffix.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	gSystem->Exec("mkdir -p "+outputDir+"/Plots_Data_QA");
	gSystem->Exec("mkdir -p "+outputDir+"/Plots_MC_QA");
	
	cout<<"Pictures are saved as "<< Suffix.Data()<< endl;
	fdate = ReturnDateString();
	
	//****************************** Specification of collision system ************************************************
	TString textProcess = ReturnMesonString (fPrefix);
	if(textProcess.CompareTo("") == 0 ){
		cout << "Meson unknown" << endl;
		return ;
	}
	
	fTextMeasurement = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
	if (meson.CompareTo("Omega")==0) fTextMeasurement = "#omega #rightarrow #pi^{+} #pi^{-} #pi^{0}";
	if (meson.CompareTo("Eta")==0) fTextMeasurement = "#eta #rightarrow #pi^{+} #pi^{-} #pi^{0}";
	
	fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
	if (fCollisionSystem.CompareTo("") == 0){
		cout << "No correct collision system specification, has been given" << endl;
		return;
	}
	fDetectionProcess = ReturnFullTextReconstructionProcess(mode);
	
	cout << "Detection process is: " << fDetectionProcess.Data() << endl;
	
	//****************************** Choice of Fitting procedure ******************************************************
	if(optionCrystalBall.CompareTo("CrystalBall") == 0){// means we want to plot values for the pi0
		fCrysFitting=1;
		cout << "CrystalBall fit chosen ..." << endl;
	} else	{
		fCrysFitting=0;
		cout << "Gaussian fit chosen ..." << endl;
	}
	
	if(cutSelection.Length() == 0){
		cout<<"ERROR: Cut selection is not set, please do!"<<endl;
		return;
	}

	//***************************** Specification Data/MC ************************************************************
	if(optionMC.CompareTo("kTRUE") == 0){
		fIsMC = 1;
		fPrefix2 = "MC";
	} else {
		fIsMC = 0;
		fPrefix2 = "data";
	}
	
	//**************************** Determine Centrality *************************************************************
	centralityString = GetCentralityString(fEventCutSelection);
	if (centralityString.CompareTo("pp")==0){
		fTextCent = "MinBias";  
	} else {
		fTextCent = Form("%s central", centralityString.Data());
	}
	cout << "here" << endl;
	
	//***************************** Initialization of variables according to meson type ******************************
	//if(meson.CompareTo("Pi0") == 0){
	//	Initialize("Pi0",numberOfBins);
	/*} else*/ if (meson.CompareTo("Eta") == 0) {
		Initialize("Eta", numberOfBins);
	//} else if (meson.CompareTo("EtaPrim") == 0) {
	//	Initialize("EtaPrim",numberOfBins);
	} else if (meson.CompareTo("Omega") == 0) {
		Initialize("Omega", numberOfBins);
	//} else if(meson.CompareTo("Pi0EtaBinning") == 0) {
		//Initialize("Pi0EtaBinning",numberOfBins);
	} else	{
		cout<<"ERROR: First argument in the ExtractSignal(....) has to be either Eta or Omega"<<endl;
		return;
	}

	cout << "here" << endl;
	
	//************************* Start of Main routine ***************************************************************
	const char* fFileErrLogDatname = Form("%s/%s/%s_%s_FileErrLog%s_%s.dat",cutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelectionRead.Data());
	fFileErrLog.open(fFileErrLogDatname, ios::out);

	TFile f(file.Data());
	TFile fMC(fileMC.Data());
	
	TString nameMainDir = "";
	//if (mode == 9 || mode == 0) {nameMainDir = "GammaConvNeutralMesonPiPlPiMiPiZero_0_9";}
	//else if (mode == 2 || mode == 3) nameMainDir = "GammaConvCalo";
	
	nameMainDir = subDir;
	
	//if (mode == 6) nameMainDir = "GammaConvNeutralMesonPiPlPiMiPiZero_1";
	
	TList *TopDir =(TList*)f.Get(nameMainDir.Data());
	if(TopDir == NULL){
		cout<<"ERROR: TopDir not Found"<<endl;
		return;
	}
	TList *TopDirMC =(TList*)fMC.Get(nameMainDir.Data());
	if(TopDirMC == NULL){
		cout<<"ERROR: TopDirMC not Found"<<endl;
		return;
	}
	
	//data file
	TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
	if(HistosGammaConversion == NULL){
		cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in DATA File"<<endl;
		return;
	}
	TList *ESDContainer = (TList*) HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
	TList *BackgroundContainer = (TList*) HistosGammaConversion->FindObject(Form("%s Back histograms",fCutSelectionRead.Data()));
	TList *MotherContainer = (TList*) HistosGammaConversion->FindObject(Form("%s Mother histograms",fCutSelectionRead.Data()));
	
	//MC file
	TList *HistosGammaConversionMC = (TList*)TopDirMC->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
	if(HistosGammaConversionMC == NULL){
		cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in MC File"<<endl;
		return;
	}
	TList *ESDContainerMC = (TList*) HistosGammaConversionMC->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
	//TList *BackgroundContainerMC = (TList*) HistosGammaConversionMC->FindObject(Form("%s Back histograms",fCutSelectionRead.Data()));
	//TList *MotherContainerMC = (TList*) HistosGammaConversionMC->FindObject(Form("%s Mother histograms",fCutSelectionRead.Data()));
	TList *TrueContainerMC = (TList*) HistosGammaConversionMC->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
	
	cout << fMesonCutSelectionRead.Data() << endl;
	cout << fGammaCutSelectionRead.Data() << endl;   
	fNumberOfGoodESDTracks = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
	fEventQuality = (TH1D*)ESDContainer->FindObject("NEvents");
		
	TString rapidityRange;
	fYMaxMeson =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
	fBackgroundMultNumber = ReturnBackgroundMult(fMesonCutSelection);
		
	TString ObjectNameESD				= "ESD_Mother_InvMass_Pt";
	TString ObjectNameBck				= "ESD_Background_InvMass_Pt";
	TString ObjectNameTrueMother		= "ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt";
	
	SetCorrectMCHistogrammNames(optionPeriod.Data(), optionEnergy.Data());
	
	//TH2D *fRatioOfESDAndFromMC_Mother, *fRatioOfESDAndFromMC_Background;
	
	fGammaGammaInvMassVSPt = (TH2D*)ESDContainer->FindObject(ObjectNameESD.Data());
	fGammaGammaInvMassVSPt -> Sumw2();
	fBckInvMassVSPt = (TH2D*)ESDContainer->FindObject(ObjectNameBck.Data());
	fBckInvMassVSPt -> Sumw2();
	//fhistoMCTB = (TH2D*)ESDContainerMC->FindObject(ObjectNameESD.Data());
	//fhistoMCTB -> Sumw2();
	
	//this is new stuff
	fTrueInvMassVSPtMC = (TH2D*)TrueContainerMC->FindObject(ObjectNameTrueMother.Data());
	fTrueInvMassVSPtMC -> Sumw2();
	fSigBckInvMassVSPtMC = (TH2D*)ESDContainerMC->FindObject(ObjectNameESD.Data());
	fSigBckInvMassVSPtMC -> Sumw2();
	
	
	if (fTrueInvMassVSPtMC == NULL) {cout << "MC file does not contain ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt histogram\n...terminating..." << endl; return;}
	
	
	if (makeQAplots)
	{
		if (!fIsMC) CreateQAHistosData(cutSelection,Suffix,optionEnergy,ESDContainer, fMode); //MY ADDITION DATA
	}
	
	const char* FileDataLogname = Form("%s/%s/%s_%s_EffiCheck_RAWDATA%s_%s.dat",cutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelectionRead.Data());
	fFileDataLog.open(FileDataLogname, ios::out);

	ProduceBckProperWeighting(BackgroundContainer,MotherContainer);
	
	//now we overwrite background histo-slices
	TH1D *fHistoDataSignalPlusBackground[fNBinsPt], *fHistoMCSignalPlusBackground[fNBinsPt],
			*fHistoScaleBetweenSigAndBck[fNBinsPt], *fHistoTrueSignal[fNBinsPt], *fHistoScaleBetweenSigAndBckBeforeFit[fNBinsPt];
	TF1 *fFunctionScaleBetweenSigAndBck[fNBinsPt], *fFunctionScaleBetweenSigAndBckFullMass[fNBinsPt];
	
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++)
	{
		Int_t startBin = fBckInvMassVSPt -> GetYaxis() -> FindBin(fBinsPt[iPt]);
		Int_t endBin = fBckInvMassVSPt -> GetYaxis() -> FindBin(fBinsPt[iPt+1]);
		fHistoMappingBackInvMassPtBin[iPt] = fBckInvMassVSPt -> ProjectionX(Form("fHistoMappingBackInvMassPtBin_%i",iPt), startBin, endBin);
		fHistoMappingBackInvMassPtBin[iPt] -> Rebin(fNRebin[iPt]);
		
		//now we create array of pt slices for Data and MC sig+bck
		/*TString nameHisto = "";
		nameHisto = Form("fHistoDataSignalPlusBackground_%i",iPt);
		fHistoDataSignalPlusBackground[iPt] = new TH1D(nameHisto.Data(),nameHisto.Data(),fGammaGammaInvMassVSPt->GetNbinsX(),0.4,0.9);
		fHistoMCSignalPlusBackground[iPt] = new TH1D(nameHisto.Data(),nameHisto.Data(),fGammaGammaInvMassVSPt->GetNbinsX(),0.4,0.9);*/
		fHistoDataSignalPlusBackground[iPt] = fGammaGammaInvMassVSPt -> ProjectionX(Form("fHistoDataSignalPlusBackground_%i",iPt), startBin, endBin);
		fHistoMCSignalPlusBackground[iPt] = fSigBckInvMassVSPtMC -> ProjectionX(Form("fHistoMCSignalPlusBackground_%i",iPt), startBin, endBin);
		fHistoTrueSignal[iPt] = fTrueInvMassVSPtMC -> ProjectionX(Form("fHistoTrueSignal_%i",iPt), startBin, endBin);
		
		fHistoDataSignalPlusBackground[iPt] -> Sumw2();
		fHistoMCSignalPlusBackground[iPt] -> Sumw2();
		fHistoTrueSignal[iPt] -> Sumw2();
		fHistoDataSignalPlusBackground[iPt] -> Rebin(fNRebin[iPt]);
		fHistoMCSignalPlusBackground[iPt] -> Rebin(fNRebin[iPt]);
		fHistoTrueSignal[iPt] -> Rebin(fNRebin[iPt]);
		
		fHistoMCSignalPlusBackground[iPt] -> Add(fHistoTrueSignal[iPt], -1.);
		
		//this is just so fHistoScaleBetweenSigAndBck is defined
		fHistoScaleBetweenSigAndBck[iPt] = (TH1D*)fHistoDataSignalPlusBackground[iPt] -> Clone(Form("fHistoScaleBetweenSigAndBck_%i",iPt));
		fHistoScaleBetweenSigAndBck[iPt] -> Sumw2();
		
		if (fIsMC)
		fHistoScaleBetweenSigAndBck[iPt] -> Divide(fHistoDataSignalPlusBackground[iPt], fHistoMCSignalPlusBackground[iPt], 1., 1., "B");
		else
		fHistoScaleBetweenSigAndBck[iPt] -> Divide(fHistoDataSignalPlusBackground[iPt], fHistoMCSignalPlusBackground[iPt], 1., 1., "");
		
		fHistoScaleBetweenSigAndBckBeforeFit[iPt] = (TH1D*)fHistoScaleBetweenSigAndBck[iPt] -> Clone(Form("fHistoScaleBetweenSigAndBckBeforeFit_%i",iPt));
		
		fFunctionScaleBetweenSigAndBck[iPt] = new TF1(Form("fFunctionScaleBetweenSigAndBck_%i",iPt),TemplateBGCorrectionFitFunctionWOPeakRegion,fMesonMassRange[0], fMesonMassRange[1], 5);
		Double_t Fit_params[5] = {1., 1., -1., 1., -1.}, *Fit_errors;
		Fit_errors = new Double_t[5];
		fFunctionScaleBetweenSigAndBck[iPt] -> SetParameters(Fit_params);
		fHistoScaleBetweenSigAndBck[iPt] -> Fit(fFunctionScaleBetweenSigAndBck[iPt]);
		fFunctionScaleBetweenSigAndBck[iPt] -> GetParameters(&Fit_params[0]);
		Fit_errors = fFunctionScaleBetweenSigAndBck[iPt] -> GetParErrors();
		
		fFunctionScaleBetweenSigAndBckFullMass[iPt] = new TF1(Form("fFunctionScaleBetweenSigAndBckFullMass_%i",iPt),TemplateBGCorrectionFitFunction,fMesonMassRange[0], fMesonMassRange[1], 5);
		fFunctionScaleBetweenSigAndBckFullMass[iPt] -> SetParameters(Fit_params);
		fFunctionScaleBetweenSigAndBckFullMass[iPt] -> SetParErrors(Fit_errors);
		//fFunctionScaleBetweenSigAndBck[iPt] -> SetRange(fMesonMassRange[0], fMesonMassRange[1]);
		Double_t binWidth, leftBinEdge, rightBinEdge;
		//cout << "------------------------------------------------------------------------------" << endl;
		for (Int_t iBin = 1; iBin <= fHistoScaleBetweenSigAndBck[iPt]->GetNbinsX(); iBin++)
		{
			binWidth = fHistoScaleBetweenSigAndBck[iPt]->GetBinWidth(iBin);
			leftBinEdge = fHistoScaleBetweenSigAndBck[iPt] -> GetBinLowEdge(iBin);
			rightBinEdge = leftBinEdge + binWidth;
			fHistoScaleBetweenSigAndBck[iPt] -> SetBinContent(iBin, fFunctionScaleBetweenSigAndBckFullMass[iPt]->Integral(leftBinEdge, rightBinEdge)/binWidth);
			fHistoScaleBetweenSigAndBck[iPt] -> SetBinError(iBin, fFunctionScaleBetweenSigAndBckFullMass[iPt]->IntegralError(leftBinEdge, rightBinEdge)/binWidth);
		}
		
		fHistoMCSignalPlusBackground[iPt] -> Multiply(fHistoScaleBetweenSigAndBck[iPt]);
		
		fHistoMappingBackInvMassPtBin[iPt] = fHistoMCSignalPlusBackground[iPt];
		
		//fHistoDataSignalPlusBackground[iPt] -> Add(fHistoMCSignalPlusBackground[iPt], -1.);
		//fHistoMCSignalPlusBackground[iPt] -> Add(fHistoMCSignalPlusBackground[iPt], -1.);
		//their plotting is dona together with plotting of all outputs... (search for PlotInvMassInPtBins or MCAndDataSigPlusBck)
	}
	
	if(fIsMC){
		cout<<"STARTED BLOCK fIsMC"<<endl;
		TList *MCContainer = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
		TList *TrueConversionContainer = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
		cout << "here" << endl;
		
		if(MCContainer == NULL)
		{
			cout<<"ERROR: " << Form("%s MC histograms",fCutSelectionRead.Data()) << " not Found in File"<<endl;
			return;
		}
		else
		cout<<"MC analysis: MCContainer successfully initialized..."<<endl;
		
		
		if(TrueConversionContainer == NULL){
		cout<<"ERROR: " << Form("%s True histograms",fCutSelectionRead.Data()) << " not Found in File"<<endl;
		return;
		}
		else
		cout<<"MC analysis: TrueConversionContainer successfully initialized..."<<endl;
		
		if (makeQAplots) CreateQAHistosMC(cutSelection, Suffix, optionEnergy, ESDContainer, MCContainer, TrueConversionContainer, fMode);
		
		/*if( fMesonId == 111){
	// 			TH2D* fHisto2DimMCMesonPtWithinAcceptance = 
	// 			TH2D* fHisto2DimMCMesonPt = 
			fHistoMCMesonPtWithinAcceptance = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0Acc.Data());
			fHistoMCMesonPt =(TH1D*)MCContainer->FindObject(ObjectNameMCPi0.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
			fHistoMCMesonPtWOWeights =(TH1D*)MCContainer->FindObject(ObjectNameMCPi0WOWeights.Data());
			
		}*/
		if( fMesonId == 221){
	// 			TH2D* fHisto2DimMCMesonPtWithinAcceptance = (TH2D*)MCContainer->FindObject(ObjectNameMCEtaAcc.Data());
	// 			TH2D* fHisto2DimMCMesonPt = (TH2D*)MCContainer->FindObject(ObjectNameMCEta.Data());
				fHistoMCMesonPtWithinAcceptance = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaAcc.Data());
			fHistoMCMesonPt = (TH1D*)MCContainer->FindObject(ObjectNameMCEta.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
			fHistoMCMesonPtWOWeights =(TH1D*)MCContainer->FindObject(ObjectNameMCEtaWOWeights.Data());
		}
		if( fMesonId == 223){
	// 			TH2D* fHisto2DimMCMesonPtWithinAcceptance = (TH2D*)MCContainer->FindObject(ObjectNameMCEtaAcc.Data());
	// 			TH2D* fHisto2DimMCMesonPt = (TH2D*)MCContainer->FindObject(ObjectNameMCEta.Data());
				fHistoMCMesonPtWithinAcceptance = (TH1D*)MCContainer->FindObject(ObjectNameMCOmegaAcc.Data());
			fHistoMCMesonPt = (TH1D*)MCContainer->FindObject(ObjectNameMCOmega.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
			fHistoMCMesonPtWOWeights =(TH1D*)MCContainer->FindObject(ObjectNameMCOmegaWOWeights.Data());
		}
		fHistoMCMesonPt->Sumw2();
		fHistoMCMesonPtWithinAcceptance->Sumw2();
		if (fHistoMCMesonPtWOWeights){
		fHistoMCMesonPtWeights = (TH1D*)fHistoMCMesonPtWOWeights->Clone("WeightsMeson");
		fHistoMCMesonPtWeights->Divide(fHistoMCMesonPt,fHistoMCMesonPtWOWeights, 1.,1.,"B");
			
		}   
		fHistoTrueMesonInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrue.Data());
		fHistoTrueMesonInvMassVSPtWOWeights = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueWOWeights.Data());
		if (fHistoTrueMesonInvMassVSPtWOWeights == NULL) fHistoTrueMesonInvMassVSPtWOWeights = (TH2D*)fHistoTrueMesonInvMassVSPt->Clone("WOWeights");
		
		fProfileTrueMesonInvMassVSPtWeights = (TProfile2D*)TrueConversionContainer->FindObject(ObjectNameProfileWeights.Data());
		if (fProfileTrueMesonInvMassVSPtWeights == NULL) 
		{
			fProfileTrueMesonInvMassVSPtWeights = new TProfile2D("name_TP2D","name_TP2D", fHistoTrueMesonInvMassVSPt->GetNbinsX(), 0., 2., fHistoTrueMesonInvMassVSPt->GetNbinsY(), 0., 25., "");
			for (Int_t i = 1; i<fProfileTrueMesonInvMassVSPtWeights->GetNbinsX(); i++)
				for (Int_t j = 1; j<fProfileTrueMesonInvMassVSPtWeights->GetNbinsY(); j++)
					fProfileTrueMesonInvMassVSPtWeights->SetBinContent(i,j,1);
		}
		
		fHistoTrueMesonInvMassVSPtReweighted = (TH2D*)fHistoTrueMesonInvMassVSPtWOWeights->Clone("Reweighted");

		fHistoTrueMesonInvMassVSPtReweighted->Multiply(fProfileTrueMesonInvMassVSPtWeights);
		FillMassMCTrueMesonHistosArray(fHistoTrueMesonInvMassVSPt);
		FillMassMCTrueReweightedMesonHistosArray(fHistoTrueMesonInvMassVSPtReweighted);
		cout << "here" << endl;
		if (fAdvancedMesonQA) {
			if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5){
				fHistoTrueMesonCaloPhotonInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloPhoton.Data());
				FillMassMCTrueMesonCaloPhotonHistosArray(fHistoTrueMesonCaloPhotonInvMassVSPt);
				fHistoTrueMesonCaloConvPhotonInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloConvPhoton.Data());
				FillMassMCTrueMesonCaloConvPhotonHistosArray(fHistoTrueMesonCaloConvPhotonInvMassVSPt);
				fHistoTrueMesonCaloElectronInvMassVSPt=(TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloElectron.Data());
				FillMassMCTrueMesonCaloElectronHistosArray(fHistoTrueMesonCaloElectronInvMassVSPt);
				fHistoTrueMesonMergedClusterInvMassVSPt= (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloMerged.Data());
				FillMassMCTrueMesonCaloMergedClusterHistosArray(fHistoTrueMesonMergedClusterInvMassVSPt);
				fHistoTrueMesonMergedClusterPartConvInvMassVSPt=(TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloMergedPartConv.Data());
				FillMassMCTrueMesonCaloMergedClusterPartConvHistosArray(fHistoTrueMesonMergedClusterPartConvInvMassVSPt);
				fHistoTrueMesonCaloEMNonLeading=(TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloEMNonLeading.Data());
				FillMassMCTrueMesonCaloEMNonLeadingHistosArray(fHistoTrueMesonCaloEMNonLeading);
			} else {
				fHistoTrueContBckInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueContBck.Data());
				FillMassMCTrueContBckHistosArray(fHistoTrueContBckInvMassVSPt);
				fHistoTrueGGBckInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueGGBck.Data());
				FillMassMCTrueGGBckHistosArray(fHistoTrueGGBckInvMassVSPt);
				fHistoTrueAllBckInvMassVSPt =(TH2D*)fHistoTrueGGBckInvMassVSPt->Clone(ObjectNameTrueAllBck.Data());
				fHistoTrueAllBckInvMassVSPt->Sumw2();
				fHistoTrueAllBckInvMassVSPt->Add(fHistoTrueContBckInvMassVSPt);
				//cout<<"ADDING 1"<<endl;
				FillMassMCTrueAllBckHistosArray(fHistoTrueAllBckInvMassVSPt);
			}	
		}
		
		//fHistoTrueSecMesonInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueSec.Data());
		//FillMassMCTrueSecMesonHistosArray(fHistoTrueSecMesonInvMassVSPt);
		//fHistoTrueSecFromK0SMesonInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueSecFromK0S.Data());
		//FillMassMCTrueSecFromK0SMesonHistosArray(fHistoTrueSecFromK0SMesonInvMassVSPt); 
		
		cout<<"ENDED BLOCK fIsMC"<<endl; 
		
	}
	
	fMesonMassExpect = TDatabasePDG::Instance()->GetParticle(fMesonId)->Mass();
	if (fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("pPb") ){
		fNEvents = fEventQuality->GetBinContent(1);
	} else {
		fNEvents =  GetNEvents(fEventQuality);
	}
	
	TH1D *fBck = (TH1D*)fBckInvMassVSPt->ProjectionX("ESD_Background_InvMass");
	TH1D *fGammaGamma = (TH1D*)fGammaGammaInvMassVSPt->ProjectionX("ESD_Mother_InvMass");
	fPeakPosAlpha01 = (TH2D*)ESDContainer->FindObject("ESD_Mother_InvMass_vs_E_alpha");

	cout<< "The mass of the meson is: "<< fMesonMassExpect<< " Events analysed: "<< fNEvents<< endl;
	cout << "here" << endl;
	// Process the 1D invariant mass histos
	fGammaGamma->SetTitle(Form("%s %s",fGammaGamma->GetTitle(),fCutSelection.Data()));
	cout << "here" << endl;
	fGammaGamma->Rebin(fNRebinGlobal);
	//fGammaGamma->Scale(1./fNRebinGlobal);
	fBck->Rebin(fNRebinGlobal);
	//fBck->Scale(1./fNRebinGlobal);
	//cout<<"1. CALL OF ProcessEM(...)"<<endl;
	ProcessEM( fGammaGamma , fBck, fBGFitRange);
	fHistoMappingBackNormInvMass = fBckNorm;
	fHistoMappingSignalInvMass = fSignal;
	
	fGammaGamma->DrawCopy();
	fHistoMappingBackNormInvMass->DrawCopy("same");
	fHistoMappingSignalInvMass->DrawCopy("same");

	// Function to Project the 2D histos InvariantMass VS Pt into Invariant Mass spectrum
	FillMassHistosArray(fGammaGammaInvMassVSPt,fPeakPosAlpha01);
	
	//THIS IS NEW!!! 
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++)
	fHistoMappingGGInvMassPtBin[iPt] = fHistoDataSignalPlusBackground[iPt];
	
	cout << "here" << endl;
	//cout<<"2. CALL OF ProcessEM(...)"<<endl;
	ProcessEM( fMesonFullPtSignal, fMesonFullPtBackground, fBGFitRange);
	fMesonFullPtBackNorm = fBckNorm;
	
	//cout<<"3. CALL OF ProcessEM(...)"<<endl;
	ProcessEM( fFittingHistMidPtSignal, fFittingHistMidPtBackground, fBGFitRange);
	fFittingHistMidPtSignalSub = fSignal;
	if(fCrysFitting==0){
		fFileErrLog << "Using exp fit"<<endl;
		FitSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntRange,200,kTRUE);
		fFitSignalInvMassMidPt = fFitReco;
	} else {
		fFileErrLog << "Using Crystal Ball function"<<endl;
		FitCBSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntRange,200,kTRUE,"SinglefitfunctionMidPt");
		fFitSignalInvMassMidPt = fFitReco;
	}
	
		if (fIsMC){
			TH1D* fHistoMappingTrueMesonInvMassPtMidPt= NULL;
			fHistoMappingTrueMesonInvMassPtMidPt= new TH1D("TrueMassMidPt","TrueMassMidPt",fHistoTrueMesonInvMassVSPtReweighted->GetNbinsX(),0.,fHistoTrueMesonInvMassVSPtReweighted->GetXaxis()->GetBinUpEdge(fHistoTrueMesonInvMassVSPtReweighted->GetNbinsX()));
			fHistoMappingTrueMesonInvMassPtMidPt->Sumw2();
			Int_t startBin = fHistoTrueMesonInvMassVSPtReweighted->GetYaxis()->FindBin(fMidPt[0]+0.001);
			Int_t endBin = fHistoTrueMesonInvMassVSPtReweighted->GetYaxis()->FindBin(fMidPt[1]-0.001);

			//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
			//       cout<< "bin values::"<< fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
			//           << fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

			fHistoTrueMesonInvMassVSPtReweighted->ProjectionX("TrueMassMidPt",startBin,endBin,"e");
				
			fHistoMappingTrueMesonInvMassPtMidPt=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
			fHistoMappingTrueMesonInvMassPtMidPt->Rebin(fNRebin[5]);
			fHistoMappingTrueMesonInvMassPtMidPt->SetLineWidth(1);
			fHistoMappingTrueMesonInvMassPtMidPt->SetLineColor(2);
			FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtMidPt, fMesonIntRange,50,kTRUE);
			
			
		}	   
	
	//    TString nameMesonFittingMidPt= Form("%s/%s_%s_MesonSubtractedFittingMidPt%s_%s_%02d.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(), fCutSelection.Data(), 200, Suffix.Data());
	//    TString nameCanvasFittingMidPt= "MesonCanvasSubtractedFittingMidPt";
	//    TString namePadFittingMidPt= "MesonPadSubtractedFittingMidPt";
	//    
	TString fDecayChannel = "#gamma#gamma";
	if (meson.CompareTo("Omega")==0) fDecayChannel = "#pi^{+} #pi^{-} #pi^{0}";
	if (meson.CompareTo("Eta")==0) fDecayChannel = "#pi^{+} #pi^{-} #pi^{0}";
	// 
	//    PlotWithFitSubtractedInvMassSinglePtBin2( fFittingHistMidPtSignalSub, fFitSignalInvMassMidPt, nameMesonFittingMidPt, nameCanvasFittingMidPt, fMesonMassRange, fIsMC,fDecayChannel);

	delete fMidPt;

	if(fEnergyFlag.CompareTo("7TeV") == 0 ){
		if(fPrefix.CompareTo("Pi0") == 0){
			for (Int_t iPt=fStartPtBin;iPt<fNBinsPeakPt;iPt++){
				fFitPeakPosPtBin[iPt] =0x00;
				if(fCrysFitting==0){
				fFileErrLog << "Using exp fit"<<endl;
				FitPeakPosInvMassInPtBins(fHistoMappingPeakPosInvMassPtBin[iPt],iPt,kTRUE);
				fFitPeakPosPtBin[iPt] = fFitReco;
				} else {
				fFileErrLog << "Using Crystal Ball function"<<endl;
				FitCBSubtractedInvMassInPtBins(fHistoMappingPeakPosInvMassPtBin[iPt], fMesonIntRange,iPt,kFALSE,Form("CBFitFuncPeakPos%02d",iPt));
				fFitPeakPosPtBin[iPt] = fFitReco;
				}
				if (fFitPeakPosPtBin[iPt] !=0x00){
				fMesonMassPeakPos[iPt] = fFitPeakPosPtBin[iPt]->GetParameter(1);
				fMesonMassPeakPosError[iPt] = fFitPeakPosPtBin[iPt]->GetParError(1);
				CalculateFWHM( fFitPeakPosPtBin[iPt]);
				fMesonFWHMAlpha01[iPt] = fFWHMFunc;
				fMesonFWHMAlpha01Error[iPt] = fFWHMFuncError;
				
				} else {
				fMesonMassPeakPos[iPt] = 0.;
				fMesonMassPeakPosError[iPt] = 0.;
				fMesonFWHMAlpha01[iPt] = 0.;
				fMesonFWHMAlpha01Error[iPt] = 0.;
				}
			}
		}
	}
	
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){ // BEGIN ANALYSIS for each Pt bin

			cout << "---------------------------------------------------------------------------------" << endl;
			cout << "Begin Analysis Pt Bin " << iPt <<endl;
			cout << "---------------------------------------------------------------------------------" << endl;
		// Function to subtract GG minus Bck
			fFileDataLog << "---------------------------------------------------------------------------------" << endl;
			fFileDataLog << "----------------------------------new pT bin ------------------------------------" << endl;
			fFileDataLog << "---------------------------------------------------------------------------------" << endl;

		ProcessBckFitSubtraction(fHistoMappingGGInvMassPtBin[iPt],iPt,fPeakRange,fFitRange,optionEnergy,Suffix,cutSelection,meson);
		
		//cout<<"4. CALL OF ProcessEM(...)"<<endl;
		ProcessEM( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRange);
		fHistoMappingSignalInvMassPtBin[iPt] = fSignal;
		fHistoMappingBackNormInvMassPtBin[iPt] = fBckNorm;

		FitWithPol2ForBG(fHistoMappingGGInvMassPtBin[iPt], fMesonIntRange,iPt,kFALSE);
		fFitWithPol2ForBG[iPt] = fFitReco;

		/*		ProcessRatioSignalBackground(fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackNormInvMassPtBin[iPt]);
							fHistoMappingRatioSBInvMassPtBin[iPt] = fRatioSB;*/

		// Integrate the 2g histo

		//       cout<< "iPt"<< iPt<< " "<< "standard range"<<endl;
		// Fitting the subtracted spectra
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/right normalization" << endl;

		fFitSignalInvMassPtBin[iPt]=0x00;
		fFitSignalInvMassBackFitPtBin[iPt]=0x00;
		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			fFileDataLog << "Subtracted mixed event" << endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntRange,iPt,kFALSE);
			fFitSignalInvMassPtBin[iPt] = fFitReco;
			fFitSignalPeakPosInvMassPtBin[iPt] = fFitGausExp;
			fFitBckInvMassPtBin[iPt] = fFitLinearBck;
			fMesonYieldsResidualBckFunc[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncError[iPt] = fIntLinearBckError;
			fFileDataLog << "Subtracted polinomial fit" << endl;
			FitSubtractedInvMassInPtBins(fHistoMappingGGInvMassBackFitPtBin[iPt], fMesonIntRange,iPt,kFALSE);
			fFitSignalInvMassBackFitPtBin[iPt] = fFitReco;
			fFitSignalPeakPosInvMassBackFitPtBin[iPt] = fFitGausExp;
			fFitBckInvMassBackFitPtBin[iPt] = fFitLinearBck;
			fMesonYieldsResidualBckFuncBackFit[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncBackFitError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntRange,iPt,kFALSE,Form("CBFitFuncNormalBin%02d",iPt));
			fFitSignalInvMassPtBin[iPt] = fFitReco;
			fFitSignalPeakPosInvMassPtBin[iPt] = fFitGausExp;
			fFitBckInvMassPtBin[iPt] = fFitLinearBck;
			fMesonYieldsResidualBckFunc[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncError[iPt] = fIntLinearBckError;

		}

				//GetFWHM
		CalculateFWHM( fFitSignalInvMassPtBin[iPt]);
		fMesonFWHM[iPt] = fFWHMFunc;
		fMesonFWHMError[iPt] = fFWHMFuncError;

		
		if (fFitSignalInvMassPtBin[iPt] !=0x00){
			fMesonMass[iPt] = fFitSignalInvMassPtBin[iPt]->GetParameter(1);
			fMesonMassError[iPt] = fFitSignalInvMassPtBin[iPt]->GetParError(1);
			fMesonWidth[iPt] = fFitSignalInvMassPtBin[iPt]->GetParameter(2);
			fMesonWidthError[iPt] = fFitSignalInvMassPtBin[iPt]->GetParError(2);
	//          if (fPrefix.Contains("Pi0")) {
				fMesonCurIntRange[0] = fMesonIntRange[0] - (fMesonMassExpect-fMesonMass[iPt]);
				fMesonCurIntRangeWide[0] = fMesonIntRangeWide[0] - (fMesonMassExpect-fMesonMass[iPt]);
				fMesonCurIntRangeNarrow[0] = fMesonIntRangeNarrow[0] - (fMesonMassExpect-fMesonMass[iPt]);
				fMesonCurIntRange[1] = fMesonIntRange[1] - (fMesonMassExpect-fMesonMass[iPt]);
				fMesonCurIntRangeWide[1] = fMesonIntRangeWide[1] - (fMesonMassExpect-fMesonMass[iPt]);
				fMesonCurIntRangeNarrow[1] = fMesonIntRangeNarrow[1] - (fMesonMassExpect-fMesonMass[iPt]);
	//          } else {
	//             fMesonCurIntRange[0] = fMesonMass[iPt] - 8*fMesonFWHM[iPt]/2.36;
	//             fMesonCurIntRangeWide[0] = fMesonMass[iPt] - 11*fMesonFWHM[iPt]/2.36;
	//             fMesonCurIntRangeNarrow[0] = fMesonMass[iPt] - 5.5*fMesonFWHM[iPt]/2.36;
	//             fMesonCurIntRange[1] = fMesonMass[iPt] + 4*fMesonFWHM[iPt]/2.36;
	//             fMesonCurIntRangeWide[1] = fMesonMass[iPt] + 5.5*fMesonFWHM[iPt]/2.36;
	//             fMesonCurIntRangeNarrow[1] = fMesonMass[iPt] + 2.5*fMesonFWHM[iPt]/2.36;
	//          }   
		} else {
			fMesonMass[iPt] = 0.;
			fMesonMassError[iPt] = 0.;
			fMesonWidth[iPt] = 0.;
			fMesonWidthError[iPt] = 0.;
			fMesonCurIntRange[0] = fMesonIntRange[0];
			fMesonCurIntRangeWide[0] = fMesonIntRangeWide[0];
			fMesonCurIntRangeNarrow[0] = fMesonIntRangeNarrow[0];
			fMesonCurIntRange[1] = fMesonIntRange[1];
			fMesonCurIntRangeWide[1] = fMesonIntRangeWide[1];
			fMesonCurIntRangeNarrow[1] = fMesonIntRangeNarrow[1];
		}

		if (fFitSignalInvMassBackFitPtBin[iPt] !=0x00){
			cout<<fFitSignalInvMassBackFitPtBin[iPt]->GetParameter(1)<<endl;
			fMesonMassBackFit[iPt] = fFitSignalInvMassBackFitPtBin[iPt]->GetParameter(1);
			fMesonMassBackFitError[iPt] = fFitSignalInvMassBackFitPtBin[iPt]->GetParError(1);
			fMesonWidthBackFit[iPt] = fFitSignalInvMassBackFitPtBin[iPt]->GetParameter(2);
			fMesonWidthBackFitError[iPt] = fFitSignalInvMassBackFitPtBin[iPt]->GetParError(2);
			fMesonCurIntRangeBackFit[0] = fMesonIntRange[0] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
			fMesonCurIntRangeBackFit[1] = fMesonIntRange[1] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
		}
		else{
			fMesonMassBackFit[iPt] = 0;
			fMesonMassBackFitError[iPt] = 0;
			fMesonWidthBackFit[iPt] = 0;
			fMesonWidthBackFitError[iPt] = 0;
			fMesonCurIntRangeBackFit[0] = fMesonIntRange[0];
			fMesonCurIntRangeBackFit[1] = fMesonIntRange[1];
		}


		IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRange);
		fGGYields[iPt] = fYields;
		fGGYieldsError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRangeWide);
		fGGYieldsWide[iPt] = fYields;
		fGGYieldsWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
		fGGYieldsNarrow[iPt] = fYields;
		fGGYieldsNarrowError[iPt] = fYieldsError;

		// Integrate the bck histo
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRange);
		fBckYields[iPt] = fYields;
		fBckYieldsError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeWide);
		fBckYieldsWide[iPt] = fYields;
		fBckYieldsWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
		fBckYieldsNarrow[iPt] = fYields;
		fBckYieldsNarrowError[iPt] = fYieldsError;

		// Integrate the signal histo
		fFileDataLog<< endl <<"Signal histo normal range, right norm "<< fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
		IntegrateHistoInvMassStream( fHistoMappingGGInvMassBackFitPtBin[iPt], fMesonCurIntRangeBackFit); //cout<<"1. call of IntegrateHistoInvMassStream"<<endl;
		fMesonYieldsBackFit[iPt] = fYields;
		fMesonYieldsBackFitError[iPt] = fYieldsError;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange); //cout<<"2. call of IntegrateHistoInvMassStream"<<endl;
		fMesonYields[iPt] = fYields;
		fMesonYieldsError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo wide range, right norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRangeWide);	 //cout<<"3. call of IntegrateHistoInvMassStream"<<endl;
		fMesonYieldsWide[iPt] = fYields;
		fMesonYieldsWideError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo narrow range, right norm" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRangeNarrow);		// cout<<"4. call of IntegrateHistoInvMassStream"<<endl;
		fMesonYieldsNarrow[iPt] = fYields;
		fMesonYieldsNarrowError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		if(fIsMC){
			fFileDataLog<< endl <<"True histo normal range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			fFitTrueSignalInvMassPtBin[iPt]=0x00;
			if(fCrysFitting==0){
				fFileErrLog << "Using exp fit"<<endl;
				FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
			} else {
				fFileErrLog << "Using Crystal Ball function"<<endl;
				FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE,Form("CBFitFuncMCTrueBin%02d",iPt));
			}

			//	FitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
			if (fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() !=0){
				fFitTrueSignalInvMassPtBin[iPt] = fFitReco;
				if (fFitTrueSignalInvMassPtBin[iPt] != 0x00){
					fMesonTrueMass[iPt] = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(1);
					fMesonTrueMassError[iPt] = fFitTrueSignalInvMassPtBin[iPt]->GetParError(1);
					CalculateFWHM(fFitTrueSignalInvMassPtBin[iPt]);
					fMesonTrueFWHM[iPt] = fFWHMFunc;
					fMesonTrueFWHMError[iPt] = fFWHMFuncError;
					fFileDataLog << "TrueFWHM \t" << fMesonTrueFWHM[iPt] << "\t +-" << fMesonTrueFWHMError[iPt] << endl;
	//                if (fPrefix.Contains("Pi0")){
					fMesonTrueIntRange[0] = fMesonIntRange[0] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeWide[0] = fMesonIntRangeWide[0] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeNarrow[0] = fMesonIntRangeNarrow[0] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRange[1] = fMesonIntRange[1] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeWide[1] = fMesonIntRangeWide[1] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeNarrow[1] = fMesonIntRangeNarrow[1] - (fMesonMassExpect-fMesonTrueMass[iPt]);
	//                } else {
	//                   fMesonTrueIntRange[0] = fMesonTrueMass[iPt] - 8*fMesonTrueFWHM[iPt]/2.36;
	//                   fMesonTrueIntRangeWide[0] = fMesonTrueMass[iPt] - 11*fMesonTrueFWHM[iPt]/2.36;
	//                   fMesonTrueIntRangeNarrow[0] = fMesonTrueMass[iPt] - 5.5*fMesonTrueFWHM[iPt]/2.36;
	//                   fMesonTrueIntRange[1] = fMesonTrueMass[iPt] + 4*fMesonTrueFWHM[iPt]/2.36;
	//                   fMesonTrueIntRangeWide[1] = fMesonTrueMass[iPt] + 5.5*fMesonTrueFWHM[iPt]/2.36;
	//                   fMesonTrueIntRangeNarrow[1] = fMesonTrueMass[iPt] - 2.5*fMesonTrueFWHM[iPt]/2.36;
	//                }   
				} else {
					fMesonTrueMass[iPt] = 0.;
					fMesonTrueMassError[iPt] = 1.;
					fMesonTrueFWHM[iPt] = 0.;
					fMesonTrueFWHMError[iPt] = 0.;
					fMesonTrueIntRange[0] = fMesonIntRange[0];
					fMesonTrueIntRangeWide[0] = fMesonIntRangeWide[0];
					fMesonTrueIntRangeNarrow[0] = fMesonIntRangeNarrow[0];
					fMesonTrueIntRange[1] = fMesonIntRange[1];
					fMesonTrueIntRangeWide[1] = fMesonIntRangeWide[1];
					fMesonTrueIntRangeNarrow[1] = fMesonIntRangeNarrow[1];

				}
			}

			fFitTrueSignalInvMassPtReweightedBin[iPt]=0x00;
			cout << "here" << endl;
			if(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]){
				cout << "Using exp fit"<<endl;
				fFileErrLog << "Using exp fit"<<endl;
				FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonIntRange,iPt,kFALSE);
				if (fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->GetEntries() !=0){
					fFitTrueSignalInvMassPtReweightedBin[iPt] = fFitReco;
					if (fFitTrueSignalInvMassPtReweightedBin[iPt] != 0x00){
					fMesonTrueMassReweighted[iPt] = fFitTrueSignalInvMassPtReweightedBin[iPt]->GetParameter(1);
					fMesonTrueMassReweightedError[iPt] = fFitTrueSignalInvMassPtReweightedBin[iPt]->GetParError(1);
					CalculateFWHM(fFitTrueSignalInvMassPtReweightedBin[iPt]);
					fMesonTrueFWHMReweighted[iPt] = fFWHMFunc;
					fMesonTrueFWHMReweightedError[iPt] = fFWHMFuncError;
		//                if (fPrefix.Contains("Pi0")){
						fMesonTrueIntReweightedRange[0] = fMesonIntRange[0] - (fMesonMassExpect-fMesonTrueMassReweighted[iPt]);
						fMesonTrueIntReweightedRangeWide[0] = fMesonIntRangeWide[0] - (fMesonMassExpect-fMesonTrueMassReweighted[iPt]);
						fMesonTrueIntReweightedRangeNarrow[0] = fMesonIntRangeNarrow[0] - (fMesonMassExpect-fMesonTrueMassReweighted[iPt]);
						fMesonTrueIntReweightedRange[1] = fMesonIntRange[1] - (fMesonMassExpect-fMesonTrueMassReweighted[iPt]);
						fMesonTrueIntReweightedRangeWide[1] = fMesonIntRangeWide[1] - (fMesonMassExpect-fMesonTrueMassReweighted[iPt]);
						fMesonTrueIntReweightedRangeNarrow[1] = fMesonIntRangeNarrow[1] - (fMesonMassExpect-fMesonTrueMassReweighted[iPt]);
		//                } else {
		//                   fMesonTrueIntReweightedRange[0] = fMesonTrueMassReweighted[iPt] - 8*fMesonTrueFWHMReweighted[iPt]/2.36;
		//                   fMesonTrueIntReweightedRangeWide[0] = fMesonTrueMassReweighted[iPt] - 11*fMesonTrueFWHMReweighted[iPt]/2.36;
		//                   fMesonTrueIntReweightedRangeNarrow[0] = fMesonTrueMassReweighted[iPt] - 5.5*fMesonTrueFWHMReweighted[iPt]/2.36;
		//                   fMesonTrueIntReweightedRange[1] = fMesonTrueMassReweighted[iPt] + 4*fMesonTrueFWHMReweighted[iPt]/2.36;
		//                   fMesonTrueIntReweightedRangeWide[1] = fMesonTrueMassReweighted[iPt] + 5.5*fMesonTrueFWHMReweighted[iPt]/2.36;
		//                   fMesonTrueIntReweightedRangeNarrow[1] = fMesonTrueMassReweighted[iPt] - 2.5*fMesonTrueFWHMReweighted[iPt]/2.36;
		//                }   
					} else {
						fMesonTrueMassReweighted[iPt] = 0.;
						fMesonTrueMassReweightedError[iPt] = 1.;
						fMesonTrueFWHMReweighted[iPt] = 0.;
						fMesonTrueFWHMReweightedError[iPt] = 0.;
						fMesonTrueIntReweightedRange[0] = fMesonIntRange[0];
						fMesonTrueIntReweightedRangeWide[0] = fMesonIntRangeWide[0];
						fMesonTrueIntReweightedRangeNarrow[0] = fMesonIntRangeNarrow[0];
						fMesonTrueIntReweightedRange[1] = fMesonIntRange[1];
						fMesonTrueIntReweightedRangeWide[1] = fMesonIntRangeWide[1];
						fMesonTrueIntReweightedRangeNarrow[1] = fMesonIntRangeNarrow[1];

					}
				}
			}

			if (fAdvancedMesonQA && (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5)){
				if (fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->GetEntries() !=0){
					FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
					fFitTrueSignalCaloPhotonInvMassPtBin[iPt] = fFitReco;
					if (fFitTrueSignalCaloPhotonInvMassPtBin[iPt] != 0x00){
						fMesonTrueMassCaloPhoton[iPt] = fFitTrueSignalCaloPhotonInvMassPtBin[iPt]->GetParameter(1);
						fMesonTrueMassErrorCaloPhoton[iPt] = fFitTrueSignalCaloPhotonInvMassPtBin[iPt]->GetParError(1);
						CalculateFWHM(fFitTrueSignalCaloPhotonInvMassPtBin[iPt]);
						fMesonTrueFWHMCaloPhoton[iPt] = fFWHMFunc;
						fMesonTrueFWHMErrorCaloPhoton[iPt] = fFWHMFuncError;
					} else {
						fMesonTrueMassCaloPhoton[iPt] = 0.;
						fMesonTrueMassErrorCaloPhoton[iPt] = 1.;
						fMesonTrueFWHMCaloPhoton[iPt] = 0.;
						fMesonTrueFWHMErrorCaloPhoton[iPt] = 0.;
					}
				}	
				if (fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->GetEntries() !=0){
					FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
					fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt] = fFitReco;
					if (fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt] != 0x00){
						fMesonTrueMassCaloConvPhoton[iPt] = fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt]->GetParameter(1);
						fMesonTrueMassErrorCaloConvPhoton[iPt] = fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt]->GetParError(1);
						CalculateFWHM(fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt]);
						fMesonTrueFWHMCaloConvPhoton[iPt] = fFWHMFunc;
						fMesonTrueFWHMErrorCaloConvPhoton[iPt] = fFWHMFuncError;
					} else {
						fMesonTrueMassCaloConvPhoton[iPt] = 0.;
						fMesonTrueMassErrorCaloConvPhoton[iPt] = 1.;
						fMesonTrueFWHMCaloConvPhoton[iPt] = 0.;
						fMesonTrueFWHMErrorCaloConvPhoton[iPt] = 0.;
					}
				}	
				if (fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->GetEntries() !=0){
					FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
					fFitTrueSignalCaloElectronInvMassPtBin[iPt] = fFitReco;
					if (fFitTrueSignalCaloElectronInvMassPtBin[iPt] != 0x00){
						fMesonTrueMassCaloElectron[iPt] = fFitTrueSignalCaloElectronInvMassPtBin[iPt]->GetParameter(1);
						fMesonTrueMassErrorCaloElectron[iPt] = fFitTrueSignalCaloElectronInvMassPtBin[iPt]->GetParError(1);
						CalculateFWHM(fFitTrueSignalCaloElectronInvMassPtBin[iPt]);
						fMesonTrueFWHMCaloElectron[iPt] = fFWHMFunc;
						fMesonTrueFWHMErrorCaloElectron[iPt] = fFWHMFuncError;
					} else {
						fMesonTrueMassCaloElectron[iPt] = 0.;
						fMesonTrueMassErrorCaloElectron[iPt] = 1.;
						fMesonTrueFWHMCaloElectron[iPt] = 0.;
						fMesonTrueFWHMErrorCaloElectron[iPt] = 0.;
					}
				}	
				if (fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt]->GetEntries() !=0){
					FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
					fFitTrueSignalCaloEMNonLeadingInvMassPtBin[iPt] = fFitReco;
					if (fFitTrueSignalCaloEMNonLeadingInvMassPtBin[iPt] != 0x00){
						fMesonTrueMassCaloEMNonLeading[iPt] = fFitTrueSignalCaloEMNonLeadingInvMassPtBin[iPt]->GetParameter(1);
						fMesonTrueMassErrorCaloEMNonLeading[iPt] = fFitTrueSignalCaloEMNonLeadingInvMassPtBin[iPt]->GetParError(1);
						CalculateFWHM(fFitTrueSignalCaloEMNonLeadingInvMassPtBin[iPt]);
						fMesonTrueFWHMCaloEMNonLeading[iPt] = fFWHMFunc;
						fMesonTrueFWHMErrorCaloEMNonLeading[iPt] = fFWHMFuncError;
					} else {
						fMesonTrueMassCaloEMNonLeading[iPt] = 0.;
						fMesonTrueMassErrorCaloEMNonLeading[iPt] = 1.;
						fMesonTrueFWHMCaloEMNonLeading[iPt] = 0.;
						fMesonTrueFWHMErrorCaloEMNonLeading[iPt] = 0.;
					}
				}	
				if (fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]->GetEntries() !=0){
					FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
					fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt] = fFitReco;
					if (fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt] != 0x00){
						fMesonTrueMassCaloMergedCluster[iPt] = fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt]->GetParameter(1);
						fMesonTrueMassErrorCaloMergedCluster[iPt] = fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt]->GetParError(1);
						CalculateFWHM(fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt]);
						fMesonTrueFWHMCaloMergedCluster[iPt] = fFWHMFunc;
						fMesonTrueFWHMErrorCaloMergedCluster[iPt] = fFWHMFuncError;
					} else {
						fMesonTrueMassCaloMergedCluster[iPt] = 0.;
						fMesonTrueMassErrorCaloMergedCluster[iPt] = 1.;
						fMesonTrueFWHMCaloMergedCluster[iPt] = 0.;
						fMesonTrueFWHMErrorCaloMergedCluster[iPt] = 0.;
					}
				}	
				if (fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]->GetEntries() !=0){
					FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
					fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt] = fFitReco;
					if (fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt] != 0x00){
						fMesonTrueMassCaloMergedClusterPartConv[iPt] = fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt]->GetParameter(1);
						fMesonTrueMassErrorCaloMergedClusterPartConv[iPt] = fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt]->GetParError(1);
						CalculateFWHM(fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt]);
						fMesonTrueFWHMCaloMergedClusterPartConv[iPt] = fFWHMFunc;
						fMesonTrueFWHMErrorCaloMergedClusterPartConv[iPt] = fFWHMFuncError;
					} else {
						fMesonTrueMassCaloMergedClusterPartConv[iPt] = 0.;
						fMesonTrueMassErrorCaloMergedClusterPartConv[iPt] = 1.;
						fMesonTrueFWHMCaloMergedClusterPartConv[iPt] = 0.;
						fMesonTrueFWHMErrorCaloMergedClusterPartConv[iPt] = 0.;
					}
				}	


				
			}
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRange);		 //cout<<"5. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueYields[iPt] = fYields;
			fMesonTrueYieldsError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"True histo wide range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);	 //cout<<"6. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueYieldsWide[iPt] = fYields;
			fMesonTrueYieldsWideError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"True histo narrow range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);		 //cout<<"7. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueYieldsNarrow[iPt] = fYields;
			fMesonTrueYieldsNarrowError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			cout<< "Analyse reweighted histos" << endl;
			
			fFileDataLog<< endl <<"True histo normal range reweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonTrueIntReweightedRange);		 //cout<<"8. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueYieldsReweighted[iPt] = fYields;
			fMesonTrueYieldsReweightedError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"True histo wide range reweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonTrueIntReweightedRangeWide);		 //cout<<"9. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueYieldsReweightedWide[iPt] = fYields;
			fMesonTrueYieldsReweightedWideError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"True histo narrow range reweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonTrueIntRangeNarrow);		 //cout<<"10. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueYieldsReweightedNarrow[iPt] = fYields;
			fMesonTrueYieldsReweightedNarrowError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			/*
			fFileDataLog<< endl <<"TrueSecTotal histo" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueSecMesonInvMassPtBins[iPt], fMesonTrueIntRange);		 cout<<"11. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueSecYields[iPt] = fYields;
			fMesonTrueSecYieldsError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
				
			fFileDataLog<< endl <<"TrueSecTotal histo wide range  " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueSecMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);		 cout<<"12. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueSecYieldsWide[iPt] = fYields;
			fMesonTrueSecYieldsWideError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
				
			fFileDataLog<< endl <<"TrueSecTotal histo narrow range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueSecMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);		 cout<<"13. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueSecYieldsNarrow[iPt] = fYields;
			fMesonTrueSecYieldsNarrowError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"TrueSecK0s histo " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt], fMesonTrueIntRange);		 cout<<"14. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueSecFromK0SYields[iPt] = fYields;
			fMesonTrueSecFromK0SYieldsError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"TrueSecK0s histo wide range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);		 cout<<"15. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueSecFromK0SYieldsWide[iPt] = fYields;
			fMesonTrueSecFromK0SYieldsWideError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"TrueSecK0s histo narrow range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);		 cout<<"16. call of IntegrateHistoInvMassStream"<<endl;
			fMesonTrueSecFromK0SYieldsNarrow[iPt] = fYields;
			fMesonTrueSecFromK0SYieldsNarrowError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
*/


			if( (fGGYields[iPt] - fMesonTrueYields[iPt]) > 0) {
				fMesonTrueSB[iPt]   = fMesonTrueYields[iPt] / ( fGGYields[iPt] - fMesonTrueYields[iPt] );
				fMesonTrueSign[iPt] = fMesonTrueYields[iPt] / pow( ( fGGYields[iPt] - fMesonTrueYields[iPt] ) , 0.5);
				fMesonTrueSBError[iPt] = 0;
				fMesonTrueSignError[iPt] = 0;
			}
			else {
				fMesonTrueSB[iPt] = 0.;
				fMesonTrueSign[iPt] = 0.;
				fMesonTrueSBError[iPt] = 0.;
				fMesonTrueSignError[iPt] = 0.;
			}
		}
		fFileDataLog<< "Residual Background leftover norm integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFunc[iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[iPt] << endl<< endl;
		fTotalBckYields[iPt] = fBckYields[iPt] + fMesonYieldsResidualBckFunc[iPt];
		fTotalBckYieldsError[iPt] = pow(fBckYieldsError[iPt]*fBckYieldsError[iPt] + fMesonYieldsResidualBckFuncError[iPt]*fMesonYieldsResidualBckFuncError[iPt],0.5);

		fMesonYieldsCorResidualBckFunc[iPt] = fMesonYields[iPt]- fMesonYieldsResidualBckFunc[iPt];
		fMesonYieldsCorResidualBckFuncError[iPt] =
			pow((fMesonYieldsError[iPt]*fMesonYieldsError[iPt]+
				fMesonYieldsResidualBckFuncError[iPt]*fMesonYieldsResidualBckFuncError[iPt]),0.5);
		fMesonYieldsPerEvent[iPt]= fMesonYieldsCorResidualBckFunc[iPt]/fNEvents;
		fMesonYieldsPerEventError[iPt]= fMesonYieldsCorResidualBckFuncError[iPt]/fNEvents;

		fMesonYieldsCorResidualBckFuncBackFit[iPt] = fMesonYieldsBackFit[iPt]- fMesonYieldsResidualBckFuncBackFit[iPt];
		fMesonYieldsCorResidualBckFuncBackFitError[iPt] =
			pow((fMesonYieldsBackFitError[iPt]*fMesonYieldsBackFitError[iPt]+
				fMesonYieldsResidualBckFuncBackFitError[iPt]*fMesonYieldsResidualBckFuncBackFitError[iPt]),0.5);
		fMesonYieldsPerEventBackFit[iPt]= fMesonYieldsCorResidualBckFuncBackFit[iPt]/fNEvents;
		fMesonYieldsPerEventBackFitError[iPt]= fMesonYieldsCorResidualBckFuncBackFitError[iPt]/fNEvents;


		//Integrate Fit Function
		IntegrateFitFunc( fFitSignalPeakPosInvMassPtBin[iPt], fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange);
		fMesonYieldsFunc[iPt]=fYieldsFunc;


			
		if( fFitBckInvMassPtBin[iPt]->Integral(fMesonMass[iPt]-fMesonFWHM[iPt], fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2))!=0){
			Double_t background = fFitBckInvMassPtBin[iPt]->Integral(fMesonMass[iPt]-fMesonFWHM[iPt], 
																	fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2));
			Double_t backgroundErr = fFitBckInvMassPtBin[iPt]->IntegralError(fMesonMass[iPt]-fMesonFWHM[iPt], 
																			fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2));
			Double_t signal = fFitSignalInvMassPtBin[iPt]->Integral(fMesonMass[iPt]-fMesonFWHM[iPt], 
																	fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2)) - background;
			Double_t signalErr =pow( pow(fFitSignalInvMassPtBin[iPt]->IntegralError(fMesonMass[iPt]-fMesonFWHM[iPt],
																					fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2)),2 )+ pow(backgroundErr,2),0.5);
			fMesonSB[iPt] = signal/ background;
			fMesonSBError[iPt] = pow( pow(signalErr/background,2.)+pow(signal/(background *background )*backgroundErr ,2.) ,0.5); 
			fMesonSign[iPt] = signal/ pow(background + signal,0.5);
			fMesonSignError[iPt] = pow(pow( (pow(background + signal ,0.5) - 0.5*signal*pow(background+signal,-0.5))/(background+signal) * signalErr ,2) + pow( 0.5*pow(signal+background,-1.5),2) ,0.5); 
		}else{
			fMesonSB[iPt] = 0.;
			fMesonSBError[iPt] = 0.;
			fMesonSign[iPt] = 0.;
			fMesonSignError[iPt] = 0.;
		}

		//       cout<< "iPt"<< iPt<< " "<< "FWHM done"<<endl;

		// Wide integration mass window
		//       cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "wide range/right normalization" << endl;

		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntRangeWide,iPt,kFALSE);
			fMesonYieldsResidualBckFuncWide[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncWideError[iPt] = fIntLinearBckError;
		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntRangeWide,iPt,kFALSE,Form("CBFitFuncNormalWideBin%02d",iPt));
			fMesonYieldsResidualBckFuncWide[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncWideError[iPt] = fIntLinearBckError;

		}

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt],fMesonIntRangeWide,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover wide integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncWide[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncWideError[iPt] <<endl<< endl;
		fMesonYieldsCorResidualBckFuncWide[iPt] = fMesonYieldsWide[iPt]- fMesonYieldsResidualBckFuncWide[iPt];

		fTotalBckYieldsWide[iPt] = fBckYieldsWide[iPt] + fMesonYieldsResidualBckFuncWide[iPt];
		fTotalBckYieldsWideError[iPt] = pow(fBckYieldsWideError[iPt]*fBckYieldsWideError[iPt] + fMesonYieldsResidualBckFuncWideError[iPt]*fMesonYieldsResidualBckFuncWideError[iPt],0.5);

		fMesonYieldsCorResidualBckFuncWideError[iPt] =
			pow((fMesonYieldsWideError[iPt]*fMesonYieldsWideError[iPt]+
				fMesonYieldsResidualBckFuncWideError[iPt]*fMesonYieldsResidualBckFuncWideError[iPt]),0.5);
		fMesonYieldsPerEventWide[iPt]= fMesonYieldsCorResidualBckFuncWide[iPt]/fNEvents;
		fMesonYieldsPerEventWideError[iPt]= fMesonYieldsCorResidualBckFuncWideError[iPt]/fNEvents;

		if( fTotalBckYieldsWide[iPt]!=0){
			fMesonSBWide[iPt] = fMesonYieldsCorResidualBckFuncWide[iPt]/fTotalBckYieldsWide[iPt];
			fMesonSBWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncWideError[iPt]/fTotalBckYieldsWide[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncWide[iPt]/(fTotalBckYieldsWide[iPt]*fTotalBckYieldsWide[iPt])*fTotalBckYieldsWideError[iPt],2.) ,0.5);
			fMesonSignWide[iPt] = fMesonYieldsCorResidualBckFuncWide[iPt]/pow(fTotalBckYieldsWide[iPt],0.5);
			fMesonSignWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncWideError[iPt]/pow(fTotalBckYieldsWide[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncWide[iPt]/pow(fTotalBckYieldsWide[iPt],1.5)*fTotalBckYieldsWideError[iPt],2.) ,0.5);

		}else{
			fMesonSBWide[iPt] = 0.;
			fMesonSBWideError[iPt] = 0.;
			fMesonSignWide[iPt] = 0.;
			fMesonSignWideError[iPt] = 0.;
		}


		// Narrow integration mass window
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "narrow range/right normalization" << endl; ;
		//       cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;

		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntRangeNarrow,iPt,kFALSE);
			fMesonYieldsResidualBckFuncNarrow[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncNarrowError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntRangeNarrow,iPt,kFALSE, Form("CBFitFuncNormalNarrowBin%02d",iPt));
			fMesonYieldsResidualBckFuncNarrow[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncNarrowError[iPt] = fIntLinearBckError;

		}

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt],fMesonIntRangeNarrow,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover narrow integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncNarrow[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncNarrowError[iPt]<<endl<< endl;
		fMesonYieldsCorResidualBckFuncNarrow[iPt] = fMesonYieldsNarrow[iPt]- fMesonYieldsResidualBckFuncNarrow[iPt];
		fMesonYieldsCorResidualBckFuncNarrowError[iPt] =
			pow((fMesonYieldsNarrowError[iPt]*fMesonYieldsNarrowError[iPt]+
				fMesonYieldsResidualBckFuncNarrowError[iPt]*fMesonYieldsResidualBckFuncNarrowError[iPt]),0.5);
		fMesonYieldsPerEventNarrow[iPt]= fMesonYieldsCorResidualBckFuncNarrow[iPt]/fNEvents;
		fMesonYieldsPerEventNarrowError[iPt]= fMesonYieldsCorResidualBckFuncNarrowError[iPt]/fNEvents;

		fTotalBckYieldsNarrow[iPt] = fBckYieldsNarrow[iPt] + fMesonYieldsResidualBckFuncNarrow[iPt];
		fTotalBckYieldsNarrowError[iPt] = pow(fBckYieldsNarrowError[iPt]*fBckYieldsNarrowError[iPt] + fMesonYieldsResidualBckFuncNarrowError[iPt]*fMesonYieldsResidualBckFuncNarrowError[iPt],0.5);


		if( fTotalBckYieldsNarrow[iPt]!=0){
			fMesonSBNarrow[iPt] = fMesonYieldsCorResidualBckFuncNarrow[iPt]/fTotalBckYieldsNarrow[iPt];
			fMesonSBNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncNarrowError[iPt]/fTotalBckYieldsNarrow[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncNarrow[iPt]/(fTotalBckYieldsNarrow[iPt]*fTotalBckYieldsNarrow[iPt])*fTotalBckYieldsNarrowError[iPt],2.) ,0.5);
			fMesonSignNarrow[iPt] = fMesonYieldsCorResidualBckFuncNarrow[iPt]/pow(fTotalBckYieldsNarrow[iPt],0.5);
			fMesonSignNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncNarrowError[iPt]/pow(fTotalBckYieldsNarrow[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncNarrow[iPt]/pow(fTotalBckYieldsNarrow[iPt],1.5)*fTotalBckYieldsNarrowError[iPt],2.) ,0.5);

		}else{
			fMesonSBNarrow[iPt] = 0.;
			fMesonSBNarrowError[iPt] = 0.;
			fMesonSignNarrow[iPt] = 0.;
			fMesonSignNarrowError[iPt] = 0.;
		}

		//////////////////////////////// Start Analysis with  Normalization at the left of the Meson Peak

		// Function to subtract GG minus Bck
		cout << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
		//cout<<"5. CALL OF ProcessEM(...)"<<endl;
		ProcessEMLeftRight( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRangeLeft, fBGFitRange);
		fHistoMappingSignalInvMassLeftPtBin[iPt] = fSignal;
		fHistoMappingBackNormInvMassLeftPtBin[iPt] = fBckNorm;


		//       cout<< "iPt"<< iPt<< " "<< "standard range"<<endl;
		// Fitting the subtracted spectra
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/left normalization" << endl;

		fFitInvMassLeftPtBin[iPt] =0x00;
		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRange,iPt,kFALSE);
			fFitInvMassLeftPtBin[iPt] = fFitReco;
			fFitSignalPeakPosInvMassLeftPtBin[iPt] = fFitGausExp;
			fFitBckInvMassLeftPtBin[iPt] = fFitLinearBck;
			fMesonYieldsResidualBckFuncLeft[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRange,iPt,kFALSE, Form("CBFitFuncLeftBin%02d",iPt));
			fFitInvMassLeftPtBin[iPt] = fFitReco;
			fFitSignalPeakPosInvMassLeftPtBin[iPt] = fFitGausExp;
			fFitBckInvMassLeftPtBin[iPt] = fFitLinearBck;
			fMesonYieldsResidualBckFuncLeft[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftError[iPt] = fIntLinearBckError;

		}
		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRange,iPt,kFALSE);
		CalculateFWHM(fFitInvMassLeftPtBin[iPt]);
		fMesonFWHMLeft[iPt] = fFWHMFunc;
		fMesonFWHMLeftError[iPt] = fFWHMFuncError;
		
		if (fFitInvMassLeftPtBin[iPt] !=0x00){
			fMesonMassLeft[iPt] = fFitInvMassLeftPtBin[iPt]->GetParameter(1);
			fMesonMassLeftError[iPt] = fFitInvMassLeftPtBin[iPt]->GetParError(1);
			fMesonWidthLeft[iPt] = fFitInvMassLeftPtBin[iPt]->GetParameter(2);
			fMesonWidthLeftError[iPt] = fFitInvMassLeftPtBin[iPt]->GetParError(2);
	//          if (fPrefix.Contains("Pi0")){
				fMesonCurLeftIntRange[0] = fMesonIntRange[0] - (fMesonMassExpect-fMesonMassLeft[iPt]);
				fMesonCurLeftIntRangeWide[0] = fMesonIntRangeWide[0] - (fMesonMassExpect-fMesonMassLeft[iPt]);
				fMesonCurLeftIntRangeNarrow[0] = fMesonIntRangeNarrow[0] - (fMesonMassExpect-fMesonMassLeft[iPt]);
				fMesonCurLeftIntRange[1] = fMesonIntRange[1] - (fMesonMassExpect-fMesonMassLeft[iPt]);
				fMesonCurLeftIntRangeWide[1] = fMesonIntRangeWide[1] - (fMesonMassExpect-fMesonMassLeft[iPt]);
				fMesonCurLeftIntRangeNarrow[1] = fMesonIntRangeNarrow[1] - (fMesonMassExpect-fMesonMassLeft[iPt]);
	//          } else {
	//             fMesonCurLeftIntRange[0] = fMesonMassLeft[iPt] - 8*fMesonFWHMLeft[iPt]/2.36;
	//             fMesonCurLeftIntRangeWide[0] = fMesonMassLeft[iPt] - 11*fMesonFWHMLeft[iPt]/2.36;
	//             fMesonCurLeftIntRangeNarrow[0] = fMesonMassLeft[iPt] - 5.5*fMesonFWHMLeft[iPt]/2.36;
	//             fMesonCurLeftIntRange[1] = fMesonMassLeft[iPt] + 4*fMesonFWHMLeft[iPt]/2.36;
	//             fMesonCurLeftIntRangeWide[1] = fMesonMassLeft[iPt] + 5.5*fMesonFWHMLeft[iPt]/2.36;
	//             fMesonCurLeftIntRangeNarrow[1] = fMesonMassLeft[iPt] - 2.5*fMesonFWHMLeft[iPt]/2.36;
	//          }   
		} else {
			fMesonMassLeft[iPt] = 0.;
			fMesonMassLeftError[iPt] = 0.;
			fMesonWidthLeft[iPt] = 0.;
			fMesonWidthLeftError[iPt] = 0.;
			fMesonCurLeftIntRange[0] = fMesonIntRange[0];
			fMesonCurLeftIntRangeWide[0] = fMesonIntRangeWide[0];
			fMesonCurLeftIntRangeNarrow[0] = fMesonIntRangeNarrow[0];
			fMesonCurLeftIntRange[1] = fMesonIntRange[1];
			fMesonCurLeftIntRangeWide[1] = fMesonIntRangeWide[1];
			fMesonCurLeftIntRangeNarrow[1] = fMesonIntRangeNarrow[1];
		}

		// Integrate the bck histo
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
		fBckYieldsLeft[iPt] = fYields;
		fBckYieldsLeftError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);
		fBckYieldsLeftWide[iPt] = fYields;
		fBckYieldsLeftWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);
		fBckYieldsLeftNarrow[iPt] = fYields;
		fBckYieldsLeftNarrowError[iPt] = fYieldsError;

		// Integrate the signal histo
		fFileDataLog<< endl <<"Signal histo normal range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);		 //cout<<"3. call of IntegrateHistoInvMassStream"<<endl;
		fMesonYieldsLeft[iPt] = fYields;
		fMesonYieldsLeftError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo wide range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);		 //cout<<"3. call of IntegrateHistoInvMassStream"<<endl;
		fMesonYieldsLeftWide[iPt] = fYields;
		fMesonYieldsLeftWideError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo narrow range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);		 //cout<<"3. call of IntegrateHistoInvMassStream"<<endl;
		fMesonYieldsLeftNarrow[iPt] = fYields;
		fMesonYieldsLeftNarrowError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;

		fFileDataLog<< "Residual Background leftover norm integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeft[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncLeftError[iPt]<<endl<< endl;
		fMesonYieldsCorResidualBckFuncLeft[iPt] = fMesonYieldsLeft[iPt]- fMesonYieldsResidualBckFuncLeft[iPt];
		fMesonYieldsCorResidualBckFuncLeftError[iPt] =
			pow((fMesonYieldsLeftError[iPt]*fMesonYieldsLeftError[iPt]+
				fMesonYieldsResidualBckFuncLeftError[iPt]*fMesonYieldsResidualBckFuncLeftError[iPt]),0.5);
		fMesonYieldsLeftPerEvent[iPt]= fMesonYieldsCorResidualBckFuncLeft[iPt]/fNEvents;
		fMesonYieldsLeftPerEventError[iPt]= fMesonYieldsCorResidualBckFuncLeftError[iPt]/fNEvents;

		fTotalBckYieldsLeft[iPt] = fBckYieldsLeft[iPt] + fMesonYieldsResidualBckFuncLeft[iPt];
		fTotalBckYieldsLeftError[iPt] = pow(fBckYieldsLeftError[iPt]*fBckYieldsLeftError[iPt] + fMesonYieldsResidualBckFuncLeftError[iPt]*fMesonYieldsResidualBckFuncLeftError[iPt],0.5);

		//Integrate Fit Function
		IntegrateFitFunc( fFitSignalPeakPosInvMassLeftPtBin[iPt], fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
		fMesonYieldsFuncLeft[iPt]=fYieldsFunc;

		//GetFWHM
		
		if( fFitBckInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2))!=0){
			Double_t background = fFitBckInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2));
			Double_t backgroundErr = fFitBckInvMassLeftPtBin[iPt]->IntegralError(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2));
			Double_t signal = fFitInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2)) - background;
			Double_t signalErr =pow( pow(fFitInvMassLeftPtBin[iPt]->IntegralError(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2)),2 )+ pow(backgroundErr,2),0.5);
			fMesonSBLeft[iPt] = signal/ background;
			fMesonSBLeftError[iPt] = pow( pow(signalErr/background,2.)+pow(signal/(background *background )*backgroundErr ,2.) ,0.5); 
			fMesonSignLeft[iPt] = signal/ pow(background + signal,0.5);
			fMesonSignLeftError[iPt] = pow(pow( (pow(background + signal ,0.5) - 0.5*signal*pow(background+signal,-0.5))/(background+signal) * signalErr ,2) + pow( 0.5*pow(signal+background,-1.5),2) ,0.5); 
		}else{
			fMesonSBLeft[iPt] = 0.;
			fMesonSBLeftError[iPt] = 0.;
			fMesonSignLeft[iPt] = 0.;
			fMesonSignLeftError[iPt] = 0.;
		}

		// Wide integration mass window
		//       cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "wide range/left normalization" << endl;

		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRangeWide,iPt,kFALSE);
			fMesonYieldsResidualBckFuncLeftWide[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftWideError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRangeWide,iPt,kFALSE, Form("CBFitFuncLeftWideBin%02d",iPt));
			fMesonYieldsResidualBckFuncLeftWide[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftWideError[iPt] = fIntLinearBckError;

		}

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt],fMesonIntRangeWide,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover wide integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeftWide[iPt]<<"\t +- \t" << fMesonYieldsResidualBckFuncLeftWideError[iPt] <<endl<< endl;
		fMesonYieldsCorResidualBckFuncLeftWide[iPt] = fMesonYieldsLeftWide[iPt]- fMesonYieldsResidualBckFuncLeftWide[iPt];
		fMesonYieldsCorResidualBckFuncLeftWideError[iPt] =
			pow((fMesonYieldsLeftWideError[iPt]*fMesonYieldsLeftWideError[iPt]+ fMesonYieldsResidualBckFuncLeftWideError[iPt]*fMesonYieldsResidualBckFuncLeftWideError[iPt]),0.5);
		fMesonYieldsLeftPerEventWide[iPt]= fMesonYieldsCorResidualBckFuncLeftWide[iPt]/fNEvents;
		fMesonYieldsLeftPerEventWideError[iPt]= fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/fNEvents;

		fTotalBckYieldsLeftWide[iPt] = fBckYieldsLeftWide[iPt] + fMesonYieldsResidualBckFuncLeftWide[iPt];
		fTotalBckYieldsLeftWideError[iPt] = pow(fBckYieldsLeftWideError[iPt]*fBckYieldsLeftWideError[iPt] + fMesonYieldsResidualBckFuncLeftWideError[iPt]*fMesonYieldsResidualBckFuncLeftWideError[iPt],0.5);


		if( fTotalBckYieldsLeftWide[iPt]!=0){
			fMesonSBLeftWide[iPt] = fMesonYieldsCorResidualBckFuncLeftWide[iPt]/fTotalBckYieldsLeftWide[iPt];
			fMesonSBLeftWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/fTotalBckYieldsLeftWide[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncLeftWide[iPt]/(fTotalBckYieldsLeftWide[iPt]*fTotalBckYieldsLeftWide[iPt])*fTotalBckYieldsLeftWideError[iPt],2.) ,0.5);
			fMesonSignLeftWide[iPt] = fMesonYieldsCorResidualBckFuncLeftWide[iPt]/pow(fTotalBckYieldsLeftWide[iPt],0.5);
			fMesonSignLeftWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/pow(fTotalBckYieldsLeftWide[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncLeftWide[iPt]/pow(fTotalBckYieldsLeftWide[iPt],1.5)*fTotalBckYieldsLeftWideError[iPt],2.) ,0.5);
		}else{
			fMesonSBLeftWide[iPt] = 0.;
			fMesonSBLeftWideError[iPt] = 0.;
			fMesonSignLeftWide[iPt] = 0.;
			fMesonSignLeftWideError[iPt] = 0.;
		}


		// Narrow integration mass window
		//       cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "narrow range/left normalization" << endl;

		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRangeNarrow,iPt,kFALSE);
			fMesonYieldsResidualBckFuncLeftNarrow[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftNarrowError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRangeNarrow,iPt,kFALSE, Form("CBFitFuncLeftNarrowBin%02d",iPt));
			fMesonYieldsResidualBckFuncLeftNarrow[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftNarrowError[iPt] = fIntLinearBckError;

		}

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt],fMesonIntRangeNarrow,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover narrow integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeftNarrow[iPt]<<"\t +- \t" << fMesonYieldsResidualBckFuncLeftNarrowError[iPt] <<endl<< endl;
		fMesonYieldsCorResidualBckFuncLeftNarrow[iPt] = fMesonYieldsLeftNarrow[iPt]- fMesonYieldsResidualBckFuncLeftNarrow[iPt];
		fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt] =
			pow((fMesonYieldsLeftNarrowError[iPt]*fMesonYieldsLeftNarrowError[iPt]+
				fMesonYieldsResidualBckFuncLeftNarrowError[iPt]*fMesonYieldsResidualBckFuncLeftNarrowError[iPt]),0.5);
		fMesonYieldsLeftPerEventNarrow[iPt]= fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/fNEvents;
		fMesonYieldsLeftPerEventNarrowError[iPt]= fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/fNEvents;

		fTotalBckYieldsLeftNarrow[iPt] = fBckYieldsLeftNarrow[iPt] + fMesonYieldsResidualBckFuncLeftNarrow[iPt];
		fTotalBckYieldsLeftNarrowError[iPt] = pow(fBckYieldsLeftNarrowError[iPt]*fBckYieldsLeftNarrowError[iPt] + fMesonYieldsResidualBckFuncLeftNarrowError[iPt]*fMesonYieldsResidualBckFuncLeftNarrowError[iPt],0.5);


		if( fTotalBckYieldsLeftNarrow[iPt]!=0){
			fMesonSBLeftNarrow[iPt] = fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/fTotalBckYieldsLeftNarrow[iPt];
			fMesonSBLeftNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/fTotalBckYieldsLeftNarrow[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/(fTotalBckYieldsLeftNarrow[iPt]*fTotalBckYieldsLeftNarrow[iPt])*fTotalBckYieldsLeftNarrowError[iPt],2.) ,0.5);
			fMesonSignLeftNarrow[iPt] = fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],0.5);
			fMesonSignLeftNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],1.5)*fTotalBckYieldsLeftNarrowError[iPt],2.) ,0.5);

		}else{
			fMesonSBLeftNarrow[iPt] = 0.;
			fMesonSBLeftNarrowError[iPt] = 0.;
			fMesonSignLeftNarrow[iPt] = 0.;
			fMesonSignLeftNarrowError[iPt] = 0.;
		}

	}



	//******************** Data OUTPUTFILE ***************************************************
	const char* fileNameSysErrDat = Form("%s/%s/%s_%s_SystematicErrorYieldExtraction_RAWDATA_%s.dat",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(), cutSelection.Data());
	fstream fileSysErrDat;
	fileSysErrDat.open(fileNameSysErrDat, ios::out);
	fileSysErrDat << "Calculation of the systematic error due to the yield extraction RAWDATA " << endl;
	fileSysErrDat <<  endl;
	fileSysErrDat << "fPiPlPiMiPiZeroYields" << endl;
	fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr" << endl;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << fGGYields[iPt] << "+-" << fGGYieldsError[iPt] << "\t" <<
			fGGYieldsWide[iPt] << "+-" << fGGYieldsWideError[iPt] << "\t" <<
			fGGYieldsNarrow[iPt] << "+-" << fGGYieldsNarrowError[iPt] << endl;

	}
	fileSysErrDat <<  endl;
	fileSysErrDat << "fTotalBckYields" << endl;
	fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
			fTotalBckYields[iPt] << "+-" << fTotalBckYieldsError[iPt] << "\t" <<
			fTotalBckYieldsWide[iPt] << "+-" << fTotalBckYieldsWideError[iPt] << "\t" <<
			fTotalBckYieldsNarrow[iPt] << "+-" << fTotalBckYieldsNarrowError[iPt] << "\t" <<
			fTotalBckYieldsLeft[iPt] << "+-" << fTotalBckYieldsLeftError[iPt]<< "\t" <<
			fTotalBckYieldsLeftWide[iPt]<< "+-" << fTotalBckYieldsLeftWideError[iPt]<< "\t" <<
			fTotalBckYieldsLeftNarrow[iPt]<< "+-" << fTotalBckYieldsLeftNarrowError[iPt] << endl;
	}
	fileSysErrDat <<  endl;
	fileSysErrDat << "fMesonYields" << endl;
	fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
			fMesonYieldsCorResidualBckFunc[iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[iPt] << "\t" <<
			fMesonYieldsCorResidualBckFuncWide[iPt] << "+-" << fMesonYieldsCorResidualBckFuncWideError[iPt] << "\t" <<
			fMesonYieldsCorResidualBckFuncNarrow[iPt] << "+-" << fMesonYieldsCorResidualBckFuncNarrowError[iPt] << "\t" <<
			fMesonYieldsCorResidualBckFuncLeft[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftError[iPt]<< "\t" <<
			fMesonYieldsCorResidualBckFuncLeftWide[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftWideError[iPt]<< "\t" <<
			fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt] << endl;
	}
	if(fIsMC){
		fileSysErrDat <<  endl;
		fileSysErrDat << "TrueYields" << endl;
		fileSysErrDat << "Bin \t True \t True Wide \t True Narr " << endl;
		for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
			fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
				fMesonTrueYields[iPt] << "\t" <<
				fMesonTrueYieldsWide[iPt] << "\t" <<
				fMesonTrueYieldsNarrow[iPt] << endl;
		}
	}
	fileSysErrDat.close();
	//******************************** OUTPUT END ******************************************************
	
	TString nameMeson1 = Form("%s/%s_%s_MCAndDataSigPlusBck%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	TString nameCanvas1 = "MCAndDataSigPlusBck";
	TString namePad1 = "MCAndDataSigPlusBck";
	PlotInvMassInPtBins( fHistoDataSignalPlusBackground, fHistoMCSignalPlusBackground, nameMeson1, nameCanvas1, namePad1,  fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);

	nameMeson1 = Form("%s/%s_%s_MCAndDataSigPlusBckScaleHistogram%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	nameCanvas1 = "MCAndDataSigPlusBckScaleHistogram";
	namePad1 = "MCAndDataSigPlusBckScaleHistogram";
	PlotInvMassInPtBins( fHistoScaleBetweenSigAndBck, fHistoScaleBetweenSigAndBckBeforeFit, nameMeson1, nameCanvas1, namePad1,  fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);

	
	TString nameMeson = Form("%s/%s_%s_MesonWithBck%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	TString nameCanvas = "MesonWithBckCanvas";
	TString namePad = "MesonWithBckPad";
	PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassPtBin, nameMeson, nameCanvas, namePad, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem, "template BG");
		
	TString nameMesonSub= Form("%s/%s_%s_MesonSubtracted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
	TString nameCanvasSub= "MesonCanvasSubtracted";
	TString namePadSub= "MesonPadSubtracted";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem, "", kTRUE, "Fit", "template background subtr. M_{#pi^{+}#pi^{-}#pi^{0}}");
	
	nameMesonSub= Form("%s/%s_%s_MesonSubtractedBackFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
	nameCanvasSub= "MesonCanvasSubtractedBackFit";
	namePadSub= "MesonPadSubtractedBackFit";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassBackFitPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassBackFitPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem, "", kTRUE, "Fit", "template background subtr. M_{#pi^{+}#pi^{-}#pi^{0}}");
	//loop through Pt slices and plot inv mass around peak (only in DATA, not MC)
	if (!fIsMC)
	for (Int_t bintoplot = fStartPtBin+1; bintoplot<numberOfBins-1;bintoplot++)
	{
		if (meson.CompareTo("Omega")==0)
		{
			Double_t ptFromPlot = fBinsOmegapPbPt[bintoplot];
			Double_t ptToPlot = fBinsOmegapPbPt[bintoplot+1];
			PlotWithFitSubtractedInvMassSinglePtBinOmega( fHistoMappingGGInvMassBackFitPtBin[bintoplot], fFitSignalInvMassBackFitPtBin[bintoplot], Form("%s/%s/%s/ExtractSignal/OmegaPlotSinglePtSignal_%i.%s",cutSelection.Data(),optionEnergy.Data(),Suffix.Data(),bintoplot,Suffix.Data()), Form("OmegaPlotSinglePtSignal_%i",bintoplot), fMesonMassRange, kFALSE, "#omega #rightarrow #pi^{+} #pi^{-} #pi^{0}",Form("p_{T}: (%.1lf,%.1lf) GeV/c",ptFromPlot,ptToPlot));
		}
		if (meson.CompareTo("Eta")==0)
		{
			Double_t ptFromPlot = fBinsEtapPbPt[bintoplot];
			Double_t ptToPlot = fBinsEtapPbPt[bintoplot+1];
			PlotWithFitSubtractedInvMassSinglePtBinOmega( fHistoMappingGGInvMassBackFitPtBin[bintoplot], fFitSignalInvMassBackFitPtBin[bintoplot], Form("%s/%s/%s/ExtractSignal/EtaPlotSinglePtSignal_%i.%s",cutSelection.Data(),optionEnergy.Data(),Suffix.Data(),bintoplot,Suffix.Data()), Form("EtaPlotSinglePtSignal_%i",bintoplot), fMesonMassRange, kFALSE, "#eta #rightarrow #pi^{+} #pi^{-} #pi^{0}",Form("p_{T}: (%.1lf,%.1lf) GeV/c",ptFromPlot,ptToPlot));
		}
		
	}
	
	nameMesonSub= Form("%s/%s_%s_MesonBackFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
	nameCanvasSub= "MesonCanvasBackFit";
	namePadSub= "MesonPadBackFit";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fBackgroundFitPol, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem, "", kTRUE, "BG fit", "M_{#pi^{+}#pi^{-}#pi^{0}}");

	nameMeson= Form("%s/%s_%s_MesonWithBckLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	nameCanvas = "MesonWithBckCanvasLeft";
	namePad = "MesonWithBckPadLeft";
	PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassLeftPtBin, nameMeson, nameCanvas, namePad,  fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem, "template BG (left+right norm.)");

	nameMesonSub= Form("%s/%s_%s_MesonWithBGAndFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	nameCanvasSub= "MesonCanvasWithBGAndFit";
	namePadSub= "MesonPadWithBGAndFit";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins,fFitWithPol2ForBG, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC , fDecayChannel, fDetectionProcess, fCollisionSystem, "", kTRUE, "BG and FG fit", "M_{#pi^{+}#pi^{-}#pi^{0}}");

	nameMesonSub= Form("%s/%s_%s_MesonSubtractedLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	nameCanvasSub= "MesonCanvasSubtractedLeft";
	namePadSub= "MesonPadSubtractedLeft";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem, "", kTRUE, "Fit", "left+right norm. template background subtr. M_{#pi^{+}#pi^{-}#pi^{0}}");

	PlotExampleInvMassBins(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin] , fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2, fThesis, fCollisionSystem, fBinsPt, fDecayChannel);
		
	// if(fEnergyFlag.CompareTo("PbPb_2.76TeV") != 0 && fEnergyFlag.CompareTo("900GeV") != 0 && fEnergyFlag.CompareTo("2.76TeV") != 0 ){
	// 	if(fPrefix.CompareTo("Pi0") == 0 ) {
	// 		nameMesonSub= Form("%s/%s_%s_MesonFittedPeakPos%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(),fCutSelection.Data(), Suffix.Data());
	// 		nameCanvasSub= "MesonFittedPeakPos";
	// 		namePadSub= "MesonFittedPeakPos";
	// 		PlotWithFitPeakPosInvMassInPtBins(fHistoMappingPeakPosInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins ,fFitPeakPosPtBin, nameMesonSub, nameCanvasSub, namePadSub,  fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPeakPt, fBinsPeakPt, fTextMeasurement, fIsMC, fDecayChannel );
	// 	}
	// }
	if(fIsMC){
		TString nameMesonTrue= Form("%s/%s_%s_TrueMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
		TString nameCanvasTrue= "TrueMesonCanvasFitted";
		TString namePadTrue= "TrueMesonPadFitted";
		PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueMesonInvMassPtBins, fFitTrueSignalInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);
			
		/*nameMesonTrue= Form("%s/%s_%s_TrueMesonSecondary%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
		nameCanvasTrue= "TrueMesonCanvasSec";
		namePadTrue= "TrueMesonPadSec";
		PlotInvMassSecondaryInPtBins( fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueSecMesonInvMassPtBins, fHistoMappingTrueSecFromK0SMesonInvMassPtBins, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);
		*/
		if (fAdvancedMesonQA && (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5)){
			nameMesonTrue= Form("%s/%s_%s_TrueMesonCaloPhoton%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
			nameCanvasTrue= "TrueMesonCaloPhotonCanvasFitted";
			namePadTrue= "TrueMesonCaloPhotonPadFitted";
			PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloPhotonInvMassPtBins, fHistoMappingTrueMesonCaloPhotonInvMassPtBins, fFitTrueSignalCaloPhotonInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);

			nameMesonTrue= Form("%s/%s_%s_TrueMesonCaloConvPhoton%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
			nameCanvasTrue= "TrueMesonCaloConvPhotonCanvasFitted";
			namePadTrue= "TrueMesonCaloConvPhotonPadFitted";
			PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins, fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins, fFitTrueSignalCaloConvPhotonInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);

			nameMesonTrue= Form("%s/%s_%s_TrueMesonCaloElectron%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
			nameCanvasTrue= "TrueMesonCaloElectronCanvasFitted";
			namePadTrue= "TrueMesonCaloElectronPadFitted";
			PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloElectronInvMassPtBins, fHistoMappingTrueMesonCaloElectronInvMassPtBins, fFitTrueSignalCaloElectronInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);

			nameMesonTrue= Form("%s/%s_%s_TrueMesonCaloMergedCluster%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
			nameCanvasTrue= "TrueMesonCaloMergedClusterCanvasFitted";
			namePadTrue= "TrueMesonCaloMergedClusterPadFitted";
			PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins, fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins, fFitTrueSignalCaloMergedClusterInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);

			nameMesonTrue= Form("%s/%s_%s_TrueMesonCaloMergedClusterPartConv%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
			nameCanvasTrue= "TrueMesonCaloMergedClusterPartConvCanvasFitted";
			namePadTrue= "TrueMesonCaloMergedClusterPartConvPadFitted";
			PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins, fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins, fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);

			nameMesonTrue= Form("%s/%s_%s_TrueMesonCaloEMNonLeading%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
			nameCanvasTrue= "TrueMesonCaloEMNonLeadingCanvasFitted";
			namePadTrue= "TrueMesonCaloEMNonLeadingPadFitted";
			PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins, fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins, fFitTrueSignalCaloEMNonLeadingInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);
			
		}
	}

	CreatePtHistos();
	FillPtHistos();

	if(fIsMC){
	
		FillHistosArrayMC(fHistoMCMesonPtWithinAcceptance,fHistoMCMesonPt,fDeltaPt,Form("%s/%s_%s_MCYieldMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));

		CalculateMesonAcceptance();
		//       cout << "Calculated MesonAcceptance" << endl;

		fNameHistoEffi="MesonEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMeson,fHistoYieldTrueSecMeson,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="MesonEffiBackFitPt";
		CalculateMesonEfficiency(fHistoYieldMesonBackFit,fHistoYieldTrueSecMeson,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonBackFitFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonEffiBackFitPt = fHistoMCMesonEffiPt;

		fNameHistoEffi="MesonWideEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonWide,fHistoYieldTrueSecMesonWide,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonWideFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonWideEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonWideEffiFitPt = fHistoMCMesonEffiFitPt;

		fNameHistoEffi="MesonNarrowEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonNarrow,fHistoYieldTrueSecMesonNarrow,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonNarrowFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonNarrowEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonNarrowEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="MesonLeftEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonLeft,fHistoYieldTrueSecMeson,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonLefFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonLeftEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonLeftEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="MesonLeftNarrowEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonLeftNarrow,fHistoYieldTrueSecMesonNarrow,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonLefNarrowFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonLeftNarrowEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonLeftNarrowEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="MesonLeftWideEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonLeftWide,fHistoYieldTrueSecMesonWide,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonLeftWideFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonLeftWideEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonLeftWideEffiFitPt = fHistoMCMesonEffiFitPt;
			
		// True Meson (only once case, because no normalization
		fNameHistoEffi="TrueMesonEffiPt";
		cout << fNameHistoEffi.Data() << endl;
			CalculateMesonEfficiency(fHistoYieldTrueMeson,NULL,fNameHistoEffi,Form("%s/%s_%s_TrueYieldMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
	//       CalculateMesonEfficiency(fHistoYieldTrueMeson,fHistoYieldTrueSecMeson,fNameHistoEffi,Form("%s/%s_%s_TrueYieldMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMCTrueMesonEffiPt = fHistoMCMesonEffiPt; 
		fHistoMCTrueMesonEffiFitPt = fHistoMCMesonEffiFitPt; 
			
		fNameHistoEffi="TrueMesonWideEffiPt";
		cout << fNameHistoEffi.Data() << endl;
			CalculateMesonEfficiency(fHistoYieldTrueMesonWide,NULL,fNameHistoEffi,Form("%s/%s_%s_TrueYieldMesonWideFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		
	//       CalculateMesonEfficiency(fHistoYieldTrueMesonWide,fHistoYieldTrueSecMesonWide,fNameHistoEffi,Form("%s/%s_%s_TrueYieldMesonWideFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMCTrueMesonWideEffiPt = fHistoMCMesonEffiPt;
		fHistoMCTrueMesonWideEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="TrueMesonNarrowEffiPt";
		cout << fNameHistoEffi.Data() << endl;
			CalculateMesonEfficiency(fHistoYieldTrueMesonNarrow,NULL,fNameHistoEffi,Form("%s/%s_%s_TrueYieldMesonNarrowFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
	//       CalculateMesonEfficiency(fHistoYieldTrueMesonNarrow,fHistoYieldTrueSecMesonNarrow,fNameHistoEffi,Form("%s/%s_%s_TrueYieldMesonNarrowFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMCTrueMesonNarrowEffiPt = fHistoMCMesonEffiPt;
		fHistoMCTrueMesonNarrowEffiFitPt = fHistoMCMesonEffiFitPt;

		fNameHistoEffi="TrueMesonEffiPtReweighted";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldTrueMesonReweighted,NULL,fNameHistoEffi,"");
		fHistoMCTrueMesonEffiPtReweighted = fHistoMCMesonEffiPt; 
		
		
		fNameHistoEffi="TrueMesonWideEffiPtReweighted";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldTrueMesonReweightedWide,NULL,fNameHistoEffi,"");
		fHistoMCTrueMesonWideEffiPtReweighted = fHistoMCMesonEffiPt;

		fNameHistoEffi="TrueMesonNarrowEffiPtReweighted";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldTrueMesonReweightedNarrow,NULL,fNameHistoEffi,"");
		fHistoMCTrueMesonNarrowEffiPtReweighted = fHistoMCMesonEffiPt;

		
		fNameHistoFrac="TrueSecFrac";
		TH1D* fHistoYieldTrueMesonSecPlusPrim = (TH1D*)fHistoYieldTrueMeson->Clone("fHistoYieldTrueMesonSecPlusPrim");
		fHistoYieldTrueMesonSecPlusPrim->Add(fHistoYieldTrueSecMeson);
		//cout<<"ADDING 2"<<endl;
		fHistoYieldTrueSecFracMeson= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrim, fHistoYieldTrueSecMeson, fNameHistoFrac);
		fNameHistoFrac="TrueSecFracFromK0S";
		fHistoYieldTrueSecFracFromK0SMeson= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrim, fHistoYieldTrueSecFromK0SMeson, fNameHistoFrac);
			
		fNameHistoFrac="TrueSecFracNarrow";
		TH1D* fHistoYieldTrueMesonSecPlusPrimNarrow = (TH1D*)fHistoYieldTrueMesonNarrow->Clone("fHistoYieldTrueMesonSecPlusPrimNarrow");
		fHistoYieldTrueMesonSecPlusPrimNarrow->Add(fHistoYieldTrueSecMesonNarrow);
		//cout<<"ADDING 3"<<endl;
		fHistoYieldTrueSecFracMesonNarrow= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimNarrow, fHistoYieldTrueSecMesonNarrow, fNameHistoFrac);
		fNameHistoFrac="TrueSecFracFromK0SNarrow";
		fHistoYieldTrueSecFracFromK0SMesonNarrow= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimNarrow, fHistoYieldTrueSecFromK0SMesonNarrow, fNameHistoFrac);
			
		fNameHistoFrac="TrueSecFracWide";
		TH1D* fHistoYieldTrueMesonSecPlusPrimWide = (TH1D*)fHistoYieldTrueMesonWide->Clone("fHistoYieldTrueMesonSecPlusPrimWide");
		fHistoYieldTrueMesonSecPlusPrimWide->Add(fHistoYieldTrueSecMesonWide);
		//cout<<"ADDING 4"<<endl;
		fHistoYieldTrueSecFracMesonWide= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimWide, fHistoYieldTrueSecMesonWide, fNameHistoFrac);
		fNameHistoFrac="TrueSecFracFromK0SWide";
		fHistoYieldTrueSecFracFromK0SMesonWide= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimWide, fHistoYieldTrueSecFromK0SMesonWide, fNameHistoFrac);

		SaveCorrectionHistos(fCutSelection, fPrefix2);
	}
	SaveHistos(fIsMC, fCutSelection, fPrefix2);
	cout << "here" << endl;
	fFileErrLog.close();
	cout << "here" << endl;
	fFileDataLog.close();
	cout << "here" << endl;
	Delete();
	cout << "here" << endl;
}

void ProcessBckFitSubtraction(TH1D *fGammaGamma, Int_t i, Double_t * fPeakRangeDummy, Double_t *fFitRangeDummy, TString energy, TString suffix, TString cutSelection, TString meson)
{
 
   fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i] = (TH1D*)fGammaGamma->Clone(Form("GG_WithoutSigal_%i",i));
   fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Sumw2();

   for (Int_t binx= 0; binx < fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()+1; binx++){
      if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinContent(binx) == 0){
         fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,1.);
         fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,0.);
      }
      if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) > fPeakRangeDummy[0] && fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) < fPeakRangeDummy[1]){
         fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,0.0);
         fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,0.0);
      }
   }
   
	TF1 *FitBackFunc;
		
   fBackgroundFitPol[i] = NULL;
   
   FitBackFunc = new TF1("BGfit",FunctionBGExclusion,fFitRangeDummy[0],fFitRangeDummy[1],5);
   fBackgroundFitPol[i] = new TF1("BGfit","pol4",fFitRangeDummy[0],fFitRangeDummy[1]); 
   fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Fit(FitBackFunc,"QMRE0","",fFitRangeDummy[0],fFitRangeDummy[1]);//QMRE0
   
   Double_t FitParams[5];
   FitBackFunc->GetParameters(&FitParams[0]);
   fBackgroundFitPol[i]->SetParameters(FitParams);
   
   fBackgroundFitPol[i]->SetRange(fMesonFitRange[0],fMesonFitRange[1]);

   for (Int_t binx= 0; binx < fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()+1; binx++){
      if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) > fFitRangeDummy[0] && fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) < fFitRangeDummy[1]){
         fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,fBackgroundFitPol[i]->Eval(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx)));
         fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,fGammaGamma->GetBinError(binx));
      }
   }
   fHistoMappingGGInvMassBackFitPtBin[i] = (TH1D*) fGammaGamma->Clone(Form("GG_SubtractedSignal_%i",i));
   fHistoMappingGGInvMassBackFitPtBin[i]->Sumw2();
   
   if (meson.CompareTo("Omega")==0)
   PlotWithFitSubtractedInvMassSinglePtBinOmega( fHistoMappingGGInvMassBackFitPtBin[i], fBackgroundFitPol[i], Form("%s/%s/%s/ExtractSignal/OmegaPlotSinglePt_%i.%s",cutSelection.Data(),energy.Data(),suffix.Data(),i,suffix.Data()), Form("OmegaPlotSinglePt_%i",i), fFitRangeDummy, kFALSE, "#omega #rightarrow #pi^{+} #pi^{-} #pi^{0}",Form("p_{T}: (%.1lf,%.1lf) GeV/c",fBinsOmegapPbPt[i],fBinsOmegapPbPt[i+1]));
   if (meson.CompareTo("Eta")==0)
   PlotWithFitSubtractedInvMassSinglePtBinOmega( fHistoMappingGGInvMassBackFitPtBin[i], fBackgroundFitPol[i], Form("%s/%s/%s/ExtractSignal/EtaPlotSinglePt_%i.%s",cutSelection.Data(),energy.Data(),suffix.Data(),i,suffix.Data()), Form("EtaPlotSinglePt_%i",i), fFitRangeDummy, kFALSE, "#eta #rightarrow #pi^{+} #pi^{-} #pi^{0}",Form("p_{T}: (%.1lf,%.1lf) GeV/c",fBinsEtapPbPt[i],fBinsEtapPbPt[i+1]));
   
   //cout<<"fHistoMappingGGInvMassBackFitPtBin Nbins = "<<fHistoMappingGGInvMassBackFitPtBin[i]->GetNbinsX()<<endl;
   //cout<<"fHistoMappingGGInvMassBackFitWithoutSignalPtBin Nbins = "<<fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()<<endl;
   fHistoMappingGGInvMassBackFitPtBin[i]->Add(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i],-1);
   
}

void ProcessEM(TH1D* fGammaGamma, TH1D* fBck, Double_t * fBGFitRangeEM)
{
   //cout<<"{"<<endl<<"ProcessEM started"<<endl;
   for (Int_t binx= 0; binx < fGammaGamma->GetNbinsX()+1; binx++){
      if(fGammaGamma->GetBinContent(binx) == 0){
         fGammaGamma->SetBinError(binx,1.);
         fGammaGamma->SetBinContent(binx,0.);
      }
   }
   
   fBckNorm = (TH1D*)fBck->Clone("fBckNorm");
   fGammaGamma->Sumw2();
   fBck->Sumw2();
   fBckNorm->Sumw2();
   
   Double_t 	r= fGammaGamma->Integral(fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[0]),fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[1]));
   Double_t 	b= fBck->Integral(fBck->GetXaxis()->FindBin(fBGFitRangeEM[0]),fBck->GetXaxis()->FindBin(fBGFitRangeEM[1]));
   Double_t 	norm = 1;
   
   if(b != 0) norm = r/b;
   fBckNorm->Scale(norm);
   //    cout<<"r="<<r<<" b="<<b<<" r/b="<<r/b<< " " << endl;
   
   Int_t numberOfZeros = 0;
   for (Int_t i = 1; i < fBckNorm->GetNbinsX()+1; i++){
      if (fBckNorm->GetBinContent(i) == 0){
         numberOfZeros++;
         if (norm > 1.){
            fBckNorm->SetBinError(i,1.);
            fBckNorm->SetBinContent(i,0.);
         }
      }
   }
//    cout << "counted " << numberOfZeros << " in the normalized BG : "<< (Double_t)numberOfZeros/fBck->GetNbinsX() << endl;
   
   fSignal = (TH1D*)fGammaGamma->Clone("fSignal");
   fSignal->Sumw2();
   /*
   //MY CODE
   cout<<"fSignal Name = "<<fSignal->GetTitle()<<endl;
   cout<<"fSignal NbinsX = "<<fSignal->GetNbinsX()<<endl;
   cout<<"fBckNorm Name = "<<fBckNorm->GetTitle()<<endl;
   cout<<"fBckNorm NbinsX = "<<fBckNorm->GetNbinsX()<<endl;
   //END OF MY CODE
   */
   if ((Double_t)numberOfZeros/fBck->GetNbinsX()< 0.25) fSignal->Add(fBckNorm,-1.);
   //cout<<"ProcessEM finished"<<endl<<"}"<<endl;
}


void ProcessRatioSignalBackground(TH1D* fGammaGamma, TH1D* fBck)
{
   for (Int_t binx= 0; binx < fGammaGamma->GetNbinsX()+1; binx++){
      if(fGammaGamma->GetBinContent(binx) == 0){
         fGammaGamma->SetBinError(binx,1.);
         fGammaGamma->SetBinContent(binx,0.);
      }
   }

   fGammaGamma->Sumw2();
   fBck->Sumw2();
   fRatioSB = (TH1D*)fGammaGamma->Clone("RatioSB");
   fRatioSB->Divide(fGammaGamma, fBck, 1.,1.,"B");
}

void ProduceBckProperWeighting(TList* fBackgroundContainer,TList* fMotherContainer){
   
   THnSparseF* fSparseMotherZM;
   THnSparseF* fSparseBckZM;
   fSparseMotherZM = (THnSparseF*)fMotherContainer->FindObject("Back_Mother_InvMass_Pt_z_m");
   fSparseBckZM = (THnSparseF*)fBackgroundContainer->FindObject("Back_Back_InvMass_Pt_z_m");
   
 
   
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fHistoWeightsBGZbinVsMbin[iPt] = new  TH2F("BGWeights", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(),fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
      fHistoWeightsBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
      fHistoWeightsBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
      fHistoWeightsBGZbinVsMbin[iPt]->Sumw2();
      fHistoFillPerEventBGZbinVsMbin[iPt] = new  TH2F("BGPoolsFillstatus", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(),fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
      fHistoFillPerEventBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
      fHistoFillPerEventBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
      fHistoFillPerEventBGZbinVsMbin[iPt]->Sumw2();
   
      for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
         for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++){ 
            // pt
            fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
            fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
            // z
            fSparseMotherZM->GetAxis(2)->SetRange(z, z);
            fSparseBckZM->GetAxis(2)->SetRange(z, z);
            // m
            fSparseMotherZM->GetAxis(3)->SetRange(m,m);
            fSparseBckZM->GetAxis(3)->SetRange(m,m);
            
            fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
            fHistoMotherZMProj->Sumw2();
            fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
            fHistoBckZMProj->Sumw2();
        
            fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
            if (m==0 && z ==0){
               if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
                  delete fHistoMappingBackInvMassPtBin[iPt];
                  fHistoMappingBackInvMassPtBin[iPt]=NULL;
               }
               fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
               fHistoMappingBackInvMassPtBin[iPt]= (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
               fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
               for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
                  fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
                  fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,0.);
               }
            }
            Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
            Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
            if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
               fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
//                cout << z << "\t" << m << "\t"  << fScalingFactorBck[z][m] << "\t" << 20./fBackgroundMultNumber << endl;
               if ( fScalingFactorBck[z][m]> (20./fBackgroundMultNumber) ){
//                   cout << "fail safe entered" << endl;
                  fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
//                   cout << "\t" <<  z << "\t" << m << "\t"  << fScalingFactorBck[z][m] << "\t" << 20./fBackgroundMultNumber << endl;
               }
            }
            fHistoMappingBackInvMassPtBin[iPt]->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
            //cout<<"ADDING 7 (iPt factor [z][m])"<<endl;
            fHistoWeightsBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fScalingFactorBck[z][m]);
            fHistoFillPerEventBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fHistoBckZMProj->GetEntries());
            fHistoMotherZMProj->Clear();
            fHistoBckZMProj->Clear();

         }
      }
      fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
      //fHistoMappingBackInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
      for (Int_t ii = 0; ii < fHistoMappingBackInvMassPtBin[iPt]->GetNbinsX()+1; ii++){
         if(fHistoMappingBackInvMassPtBin[iPt]->GetBinContent(ii) == 0){
            fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
            fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,1.);
         }
      }
      fFileDataLog << "Scaling Background factors for Pt bin " << iPt << " z m " << endl;
      // 		cout << "Scaling Background factors for Pt bin " << iPt << " z m " << endl;
      for (Int_t z=0; z < fSparseMotherZM->GetAxis(2)->GetNbins()-1; z++){
         fFileDataLog << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
         // 			cout << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
      }
   }
   
   for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
      for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++) {

         // pt
         fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fFullPt[1]-0.001));
         fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fFullPt[1]-0.001));
         // z
         fSparseMotherZM->GetAxis(2)->SetRange(z+1, z+1);
         fSparseBckZM->GetAxis(2)->SetRange(z+1, z+1);
         // m
         fSparseMotherZM->GetAxis(3)->SetRange(m+1,m+1);
         fSparseBckZM->GetAxis(3)->SetRange(m+1,m+1);
         
         fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
         fHistoMotherZMProj->Sumw2();
         fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
         fHistoBckZMProj->Sumw2();
         
         fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
         if (m==0 && z ==0){
            fNameHistoBack = "Mapping_Back_InvMass_FullPt";
            fMesonFullPtBackground = (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
            fMesonFullPtBackground->Sumw2();
            for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
               fMesonFullPtBackground->SetBinContent(ii,0.);
               fMesonFullPtBackground->SetBinError(ii,0.);
            }
         }
         Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
         Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
         if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
            fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
            if ( fScalingFactorBck[z][m]>20./fBackgroundMultNumber ){
               fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
            }
         }
         fMesonFullPtBackground->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
         //cout<<"ADDING 8"<<endl;

         fHistoMotherZMProj->Clear();
         fHistoBckZMProj->Clear();
      }
   }
   fMesonFullPtBackground->Rebin(fNRebin[4]);
   //fMesonFullPtBackground->Scale(1./fNRebin[4]);
   

   for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
      for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++) {
        
         // pt
         fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fMidPt[1]-0.001));
         fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fMidPt[1]-0.001));
         // z
         fSparseMotherZM->GetAxis(2)->SetRange(z+1, z+1);
         fSparseBckZM->GetAxis(2)->SetRange(z+1, z+1);
         // m
         fSparseMotherZM->GetAxis(3)->SetRange(m+1,m+1);
         fSparseBckZM->GetAxis(3)->SetRange(m+1,m+1);
         
         fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
         fHistoMotherZMProj->Sumw2();
         fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
         fHistoBckZMProj->Sumw2();
         
         fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
         if (m==0 && z ==0){
            fNameHistoBack = "Mapping_Back_InvMass_MidPt";
            fFittingHistMidPtBackground = (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
            fFittingHistMidPtBackground->Sumw2();
            for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
               fFittingHistMidPtBackground->SetBinContent(ii,0.);
               fFittingHistMidPtBackground->SetBinError(ii,0.);
            }
         }
         Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
         Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
         if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
            fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
            if ( fScalingFactorBck[z][m]>20./fBackgroundMultNumber ){
               fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
            }
         }
         fFittingHistMidPtBackground->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
         //cout<<"ADDING 9"<<endl;
         
         fHistoMotherZMProj->Clear();
         fHistoBckZMProj->Clear();
      }
   }
   fFittingHistMidPtBackground->Rebin(fNRebin[4]);
   //fFittingHistMidPtBackground->Scale(1./fNRebin[4]);
}

void ProduceBckWithoutWeighting(TH2D *fBckInvMassVSPtDummy){
   //calculation background for midPt without weighting
   Int_t startBinMidPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[0]+0.001);
   Int_t endBinMidPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[1]-0.001);
   fFittingHistMidPtBackground = new TH1D("Mapping_Back_InvMass_MidPt","Mapping_Back_InvMass_MidPt",fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
	fFittingHistMidPtBackground->Sumw2();
   fBckInvMassVSPtDummy->ProjectionX("Mapping_Back_InvMass_MidPt",startBinMidPt,endBinMidPt);
   fFittingHistMidPtBackground=(TH1D*)gDirectory->Get("Mapping_Back_InvMass_MidPt");
   fFittingHistMidPtBackground->Rebin(fNRebin[4]);

   //calulation background for fullPt without weighting
   fMesonFullPtBackground = new TH1D("Mapping_Back_InvMass_FullPt","Mapping_Back_InvMass_FullPt",fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
   fMesonFullPtBackground->Sumw2();
   Int_t startBinFullPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[0]+0.001);
   Int_t endBinFullPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[1]-0.001);
   fBckInvMassVSPtDummy->ProjectionX("Mapping_Back_InvMass_FullPt",startBinFullPt,endBinFullPt);
   fMesonFullPtBackground=(TH1D*)gDirectory->Get("Mapping_Back_InvMass_FullPt");
   fMesonFullPtBackground->Rebin(fNRebin[4]);

   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
      if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
         delete fHistoMappingBackInvMassPtBin[iPt];
         fHistoMappingBackInvMassPtBin[iPt]=NULL;
      }
      fHistoMappingBackInvMassPtBin[iPt]=new TH1D(fNameHistoBack.Data(),fNameHistoBack.Data(),fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
		fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
      Int_t startBin = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
      Int_t endBin = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

      fBckInvMassVSPtDummy->ProjectionX(fNameHistoBack.Data(),startBin,endBin);
      fHistoMappingBackInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoBack.Data());
      if(fNRebin[iPt]>1){
         fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
      }

   }
}


void FillMassHistosArray(TH2D* fGammaGammaInvMassVSPtDummy, TH2D *fAlphaPeakPosDummy)
{
   fFittingHistMidPtSignal = new TH1D("Mapping_GG_InvMass_MidPt","Mapping_GG_InvMass_MidPt",fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,(double)fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
   cout << fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()) << endl;
   Int_t startBinMidPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[0]+0.001);
   Int_t endBinMidPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[1]-0.001);
	fFittingHistMidPtSignal->Sumw2();
   fGammaGammaInvMassVSPtDummy->ProjectionX("Mapping_GG_InvMass_MidPt",startBinMidPt,endBinMidPt);

   fFittingHistMidPtSignal=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_MidPt");
   fFittingHistMidPtSignal->Rebin(fNRebin[4]);
   //fFittingHistMidPtSignal->Scale(1./fNRebin[4]);
      cout << "Mid pt geschrieben" << endl;

   fMesonFullPtSignal = new TH1D("Mapping_GG_InvMass_FullPt","Mapping_GG_InvMass_FullPt",fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
	fMesonFullPtSignal->Sumw2();
   Int_t startBinFullPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[0]+0.001);
   Int_t endBinFullPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[1]-0.001);

   fGammaGammaInvMassVSPtDummy->ProjectionX("Mapping_GG_InvMass_FullPt",startBinFullPt,endBinFullPt);

   fMesonFullPtSignal=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_FullPt");
   fMesonFullPtSignal->Rebin(fNRebin[4]);
   //fMesonFullPtSignal->Scale(1./fNRebin[4]);
      cout << "Full pt geschrieben" << endl;

   if(fEnergyFlag.CompareTo("PbPb_2.76TeV") != 0 && fEnergyFlag.CompareTo("900GeV") != 0 && fEnergyFlag.CompareTo("2.76TeV") != 0 && fEnergyFlag.CompareTo("pPb_5.023TeV") != 0){
      if(fPrefix.CompareTo("Pi0") == 0){
         for(Int_t iPt=fStartPtBin;iPt<fNBinsPeakPt;iPt++){
            fNameHistoPP = Form("Mapping_PP_InvMass_in_Pt_Bin%02d", iPt);
            if(fHistoMappingPeakPosInvMassPtBin[iPt]!= NULL){
               delete fHistoMappingPeakPosInvMassPtBin[iPt];
               fHistoMappingPeakPosInvMassPtBin[iPt]=NULL;
            }
            fHistoMappingPeakPosInvMassPtBin[iPt]=new TH1D(fNameHistoPP.Data(),fNameHistoPP.Data(),fAlphaPeakPosDummy->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
				fHistoMappingPeakPosInvMassPtBin[iPt]->Sumw2();
            Int_t startBin = fAlphaPeakPosDummy->GetYaxis()->FindBin(fBinsPeakPt[iPt]+0.001);
            Int_t endBin = fAlphaPeakPosDummy->GetYaxis()->FindBin(fBinsPeakPt[iPt+1]-0.001);

            fAlphaPeakPosDummy->ProjectionX(fNameHistoPP.Data(),startBin,endBin);
            fHistoMappingPeakPosInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoPP.Data());

            if(fBinsPeakPtRebin[iPt]>1){
               fHistoMappingPeakPosInvMassPtBin[iPt]->Rebin(fBinsPeakPtRebin[iPt]);
               //fHistoMappingPeakPosInvMassPtBin[iPt]->Scale(1./fBinsPeakPtRebin[iPt]);
            }
         }
      }
   }
      cout << "nach Peak Pos" << endl;
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fNameHistoGG = Form("Mapping_GG_InvMass_in_Pt_Bin%02d", iPt);

      if(fHistoMappingGGInvMassPtBin[iPt]!= NULL){
         delete fHistoMappingGGInvMassPtBin[iPt];
         fHistoMappingGGInvMassPtBin[iPt]=NULL;
      }
      fHistoMappingGGInvMassPtBin[iPt]=new TH1D(fNameHistoGG.Data(),fNameHistoGG.Data(),fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
		fHistoMappingGGInvMassPtBin[iPt]->Sumw2();

      Int_t startBin = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
      Int_t endBin = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

//             cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
//             cout<< "bin values::"<< fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(startBin)<< " "
//                 << fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(endBin)<< endl;

      fGammaGammaInvMassVSPtDummy->ProjectionX(fNameHistoGG.Data(),startBin,endBin);


      fHistoMappingGGInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoGG.Data());

      if(fNRebin[iPt]>1){
         fHistoMappingGGInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
         //fHistoMappingGGInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
      }
   }
   //    cout << "each pt written" << endl;

}

void FillMassMCTrueMesonHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill)
{
	fHistoTrueMesonInvMassVSPtFill->Sumw2();
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fNameHistoTrue = Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);
      if(fHistoMappingTrueMesonInvMassPtBins[iPt]!= NULL){
         delete fHistoMappingTrueMesonInvMassPtBins[iPt];
         fHistoMappingTrueMesonInvMassPtBins[iPt]=NULL;
      }
      fHistoMappingTrueMesonInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonInvMassVSPtFill->GetNbinsX()));
		fHistoMappingTrueMesonInvMassPtBins[iPt]->Sumw2();
      Int_t startBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
      Int_t endBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

      //       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
      //       cout<< "bin values::"<< fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
      //           << fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

      fHistoTrueMesonInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");
		
      fHistoMappingTrueMesonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
		cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() << endl;
      if(fNRebin[iPt]>1){
         fHistoMappingTrueMesonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
      }
      fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineWidth(1);
      fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineColor(2);
   }

}

void FillMassMCTrueMesonCaloPhotonHistosArray(TH2D* fHistoTrueMesonCaloPhotonInvMassVSPtFill)
{
	fHistoTrueMesonCaloPhotonInvMassVSPtFill->Sumw2();
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrue = Form("Mapping_TrueMesonCaloPhoton_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt];
			fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonCaloPhotonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonCaloPhotonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonCaloPhotonInvMassVSPtFill->GetNbinsX()));
			fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->Sumw2();
		Int_t startBin = fHistoTrueMesonCaloPhotonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueMesonCaloPhotonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fHistoTrueMesonCaloPhotonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fHistoTrueMesonCaloPhotonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueMesonCaloPhotonInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");
			
		fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
			cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->GetEntries() << endl;
		if(fNRebin[iPt]>1){
			fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
		}
		fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->SetLineColor(2);
	}
}

void FillMassMCTrueMesonCaloConvPhotonHistosArray(TH2D* fHistoTrueMesonCaloConvPhotonInvMassVSPtFill)
{
	fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->Sumw2();
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrue = Form("Mapping_TrueMesonCaloConvPhoton_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt];
			fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->GetNbinsX()));
			fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->Sumw2();
		Int_t startBin = fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");
			
		fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
			cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->GetEntries() << endl;
		if(fNRebin[iPt]>1){
			fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
		}
		fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->SetLineColor(2);
	}
}


void FillMassMCTrueMesonCaloElectronHistosArray(TH2D* fHistoTrueMesonCaloElectronInvMassVSPtFill)
{
	fHistoTrueMesonCaloElectronInvMassVSPtFill->Sumw2();
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrue = Form("Mapping_TrueMesonCaloElectron_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt];
			fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonCaloElectronInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonCaloElectronInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonCaloElectronInvMassVSPtFill->GetNbinsX()));
			fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->Sumw2();
		Int_t startBin = fHistoTrueMesonCaloElectronInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueMesonCaloElectronInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fHistoTrueMesonCaloElectronInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fHistoTrueMesonCaloElectronInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueMesonCaloElectronInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");
			
		fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
			cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->GetEntries() << endl;
		if(fNRebin[iPt]>1){
			fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
		}
		fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->SetLineColor(2);
	}
}

void FillMassMCTrueMesonCaloEMNonLeadingHistosArray(TH2D* fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill)
{
	fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill->Sumw2();
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrue = Form("Mapping_TrueMesonCaloEMNonLeading_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt];
			fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill->GetNbinsX()));
			fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt]->Sumw2();
		Int_t startBin = fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueMesonCaloEMNonLeadingInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");
			
		fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
			cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt]->GetEntries() << endl;
		if(fNRebin[iPt]>1){
			fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
		}
		fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[iPt]->SetLineColor(2);
	}
}

void FillMassMCTrueMesonCaloMergedClusterHistosArray(TH2D* fHistoTrueMesonCaloMergedClusterInvMassVSPtFill)
{
	fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->Sumw2();
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrue = Form("Mapping_TrueMesonCaloMergedCluster_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt];
			fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->GetNbinsX()));
			fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]->Sumw2();
		Int_t startBin = fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");
			
		fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
			cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]->GetEntries() << endl;
		if(fNRebin[iPt]>1){
			fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
		}
		fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]->SetLineColor(2);
	}
}

void FillMassMCTrueMesonCaloMergedClusterPartConvHistosArray(TH2D* fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill)
{
	fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->Sumw2();
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrue = Form("Mapping_TrueMesonCaloMergedClusterPartConv_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt];
			fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->GetNbinsX()));
			fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]->Sumw2();
		Int_t startBin = fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");
			
		fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
			cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]->GetEntries() << endl;
		if(fNRebin[iPt]>1){
			fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
		}
		fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]->SetLineColor(2);
	}
}

void FillMassMCTrueReweightedMesonHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill)
{
   fHistoTrueMesonInvMassVSPtFill->Sumw2();
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fNameHistoTrue = Form("Mapping_TrueMeson_InvMassReweighted_in_Pt_Bin%02d", iPt);
      if(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]!= NULL){
         delete fHistoMappingTrueMesonInvMassPtReweightedBins[iPt];
         fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]=NULL;
      }
      fHistoMappingTrueMesonInvMassPtReweightedBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonInvMassVSPtFill->GetNbinsX()));
      fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->Sumw2();
      Int_t startBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
      Int_t endBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

      //       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
      //       cout<< "bin values::"<< fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
      //           << fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

      fHistoTrueMesonInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");
      
      fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
      cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->GetEntries() << endl;
      if(fNRebin[iPt]>1){
         fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->Rebin(fNRebin[iPt]);
      }
      fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->SetLineWidth(1);
      fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->SetLineColor(2);
   }

}


void FillMassMCTrueGGBckHistosArray(TH2D* fHistoTrueGGBckInvMassVSPtFill)
{
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fNameHistoTrueGGBck = Form("Mapping_TrueGGBck_InvMass_in_Pt_Bin%02d", iPt);
      if(fHistoMappingTrueGGBckInvMassPtBins[iPt]!= NULL){
         delete fHistoMappingTrueGGBckInvMassPtBins[iPt];
         fHistoMappingTrueGGBckInvMassPtBins[iPt]=NULL;
      }
      fHistoMappingTrueGGBckInvMassPtBins[iPt] = new TH1D(fNameHistoTrueGGBck.Data(),fNameHistoTrueGGBck.Data(),fHistoTrueGGBckInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueGGBckInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueGGBckInvMassVSPtFill->GetNbinsX()));
		fHistoMappingTrueGGBckInvMassPtBins[iPt]->Sumw2();
      Int_t startBin = fHistoTrueGGBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
      Int_t endBin = fHistoTrueGGBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

      //       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
      //       cout<< "bin values::"<< fHistoTrueGGBckInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
      //           << fHistoTrueGGBckInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

      fHistoTrueGGBckInvMassVSPtFill->ProjectionX(fNameHistoTrueGGBck.Data(),startBin,endBin);
      fHistoMappingTrueGGBckInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueGGBck.Data());
      if(fNRebin[iPt]>1){
         fHistoMappingTrueGGBckInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
         //fHistoMappingTrueGGBckInvMassPtBins[iPt]->Scale(1./fNRebin[iPt]);

      }
      fHistoMappingTrueGGBckInvMassPtBins[iPt]->SetLineWidth(1);
      fHistoMappingTrueGGBckInvMassPtBins[iPt]->SetLineColor(3);
   }

}

void FillMassMCTrueContBckHistosArray(TH2D* fHistoTrueContBckInvMassVSPtFill)
{
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fNameHistoTrueContBck = Form("Mapping_TrueContBck_InvMass_in_Pt_Bin%02d", iPt);
      if(fHistoMappingTrueContBckInvMassPtBins[iPt]!= NULL){
         delete fHistoMappingTrueContBckInvMassPtBins[iPt];
         fHistoMappingTrueContBckInvMassPtBins[iPt]=NULL;
      }
      fHistoMappingTrueContBckInvMassPtBins[iPt] = new TH1D(fNameHistoTrueContBck.Data(),fNameHistoTrueContBck.Data(),fHistoTrueContBckInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueContBckInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueContBckInvMassVSPtFill->GetNbinsX()));
		fHistoMappingTrueContBckInvMassPtBins[iPt]->Sumw2();
      Int_t startBin = fHistoTrueContBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
      Int_t endBin = fHistoTrueContBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

      //       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
      //       cout<< "bin values::"<< fHistoTrueContBckInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
      //           << fHistoTrueContBckInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

      fHistoTrueContBckInvMassVSPtFill->ProjectionX(fNameHistoTrueContBck.Data(),startBin,endBin);
      fHistoMappingTrueContBckInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueContBck.Data());
      if(fNRebin[iPt]>1){
         fHistoMappingTrueContBckInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
         //fHistoMappingTrueContBckInvMassPtBins[iPt]->Scale(1./fNRebin[iPt]);

      }
      fHistoMappingTrueContBckInvMassPtBins[iPt]->SetLineWidth(1);
      fHistoMappingTrueContBckInvMassPtBins[iPt]->SetLineColor(5);
   }
}

void FillMassMCTrueAllBckHistosArray(TH2D* fHistoTrueAllBckInvMassVSPtFill)
{
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fNameHistoTrueAllBck = Form("Mapping_TrueAllBck_InvMass_in_Pt_Bin%02d", iPt);
      if(fHistoMappingTrueAllBckInvMassPtBins[iPt]!= NULL){
         delete fHistoMappingTrueAllBckInvMassPtBins[iPt];
         fHistoMappingTrueAllBckInvMassPtBins[iPt]=NULL;
      }
      fHistoMappingTrueAllBckInvMassPtBins[iPt] = new TH1D(fNameHistoTrueAllBck.Data(),fNameHistoTrueAllBck.Data(),fHistoTrueAllBckInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueAllBckInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueAllBckInvMassVSPtFill->GetNbinsX()));
		fHistoMappingTrueAllBckInvMassPtBins[iPt]->Sumw2();
      Int_t startBin = fHistoTrueAllBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
      Int_t endBin = fHistoTrueAllBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

      //       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
      //       cout<< "bin values::"<< fHistoTrueAllBckInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
      //           << fHistoTrueAllBckInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

      fHistoTrueAllBckInvMassVSPtFill->ProjectionX(fNameHistoTrueAllBck.Data(),startBin,endBin);
      fHistoMappingTrueAllBckInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueAllBck.Data());
      if(fNRebin[iPt]>1){
         fHistoMappingTrueAllBckInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
         //fHistoMappingTrueAllBckInvMassPtBins[iPt]->Scale(1./fNRebin[iPt]);

      }
      fHistoMappingTrueAllBckInvMassPtBins[iPt]->SetLineWidth(1);
      fHistoMappingTrueAllBckInvMassPtBins[iPt]->SetLineColor(2);
   }
}

void FillMassMCTrueSecMesonHistosArray(TH2D* fHistoTrueSecMesonInvMassVSPtFill)
{
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fNameHistoTrueSec = Form("Mapping_TrueSecMeson_InvMass_in_Pt_Bin%02d", iPt);
      if(fHistoMappingTrueSecMesonInvMassPtBins[iPt]!= NULL){
         delete fHistoMappingTrueSecMesonInvMassPtBins[iPt];
         fHistoMappingTrueSecMesonInvMassPtBins[iPt]=NULL;
      }
      fHistoMappingTrueSecMesonInvMassPtBins[iPt] = new TH1D(fNameHistoTrueSec.Data(),fNameHistoTrueSec.Data(),fHistoTrueSecMesonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueSecMesonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueSecMesonInvMassVSPtFill->GetNbinsX()));
		fHistoMappingTrueSecMesonInvMassPtBins[iPt]->Sumw2();
      Int_t startBin = fHistoTrueSecMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
      Int_t endBin = fHistoTrueSecMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

      //       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
      //       cout<< "bin values::"<< fHistoTrueSecMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
      //           << fHistoTrueSecMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

      fHistoTrueSecMesonInvMassVSPtFill->ProjectionX(fNameHistoTrueSec.Data(),startBin,endBin);
      fHistoMappingTrueSecMesonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueSec.Data());
      if(fNRebin[iPt]>1){
         fHistoMappingTrueSecMesonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
         //fHistoMappingTrueSecMesonInvMassPtBins[iPt]->Scale(1./fNRebin[iPt]);

      }
      fHistoMappingTrueSecMesonInvMassPtBins[iPt]->SetLineWidth(1);
      fHistoMappingTrueSecMesonInvMassPtBins[iPt]->SetLineColor(2);
   }
}

void FillMassMCTrueSecFromK0SMesonHistosArray(TH2D* fHistoTrueSecFromK0SMesonInvMassVSPtFill)
{
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
      fNameHistoTrueSecFromK0S = Form("Mapping_TrueSecFromK0SMeson_InvMass_in_Pt_Bin%02d", iPt);
      if(fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]!= NULL){
         delete fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt];
         fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]=NULL;
      }
      fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt] = new TH1D(fNameHistoTrueSecFromK0S.Data(),fNameHistoTrueSecFromK0S.Data(),fHistoTrueSecFromK0SMesonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueSecFromK0SMesonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueSecFromK0SMesonInvMassVSPtFill->GetNbinsX()));
		fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]->Sumw2();
      Int_t startBin = fHistoTrueSecFromK0SMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
      Int_t endBin = fHistoTrueSecFromK0SMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

      //       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
      //       cout<< "bin values::"<< fHistoTrueSecFromK0SMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
      //           << fHistoTrueSecFromK0SMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

      fHistoTrueSecFromK0SMesonInvMassVSPtFill->ProjectionX(fNameHistoTrueSecFromK0S.Data(),startBin,endBin);
      fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueSecFromK0S.Data());
      if(fNRebin[iPt]>1){
         fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
         //fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]->Scale(1./fNRebin[iPt]);

      }
      fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]->SetLineWidth(1);
      fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]->SetLineColor(2);
   }

}

TH1D* CalculateSecondaryFractions(TH1D* histoRawYield, TH1D* histoRawYieldSec, TString nameHistoFrac){
   TH1D* histoFracSec = (TH1D*)histoRawYieldSec->Clone(nameHistoFrac.Data());
   histoFracSec->Divide(histoFracSec,histoRawYield,1.,1.,"B");
   return histoFracSec;
}

void CreatePtHistos(){

	fDeltaPt =			 	new TH1D("deltaPt","",fNBinsPt,fBinsPt);

	fHistoYieldMeson = 			new TH1D("histoYieldMeson","",fNBinsPt,fBinsPt);
	fHistoYieldMesonBackFit = 			new TH1D("histoYieldMesonBackFit","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMeson = 		new TH1D("histoYieldTrueMeson","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMesonReweighted =     new TH1D("histoYieldTrueMesonReweighted","",fNBinsPt,fBinsPt);
	fHistoYieldTrueSecMeson = 		new TH1D("histoYieldTrueSecMeson","",fNBinsPt,fBinsPt);
	fHistoYieldTrueSecFromK0SMeson = 		new TH1D("histoYieldTrueSecFromK0SMeson","",fNBinsPt,fBinsPt);

	fHistoYieldMesonPerEvent = 	new TH1D("histoYieldMesonPerEvent","",fNBinsPt,fBinsPt);
	fHistoYieldMesonPerEventBackFit = 	new TH1D("histoYieldMesonPerEventBackFit","",fNBinsPt,fBinsPt);
	fHistoSignMeson = 			new TH1D("histoSignMeson","",fNBinsPt,fBinsPt);
	fHistoSBMeson = 			new TH1D("histoSBMeson","",fNBinsPt,fBinsPt);
	if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0){fHistoMassPosition = new TH1D("histoMassPosition","",fNBinsPeakPt,fBinsPeakPtHalf);}
	if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0){fHistoFWHMMesonAlpha01 = new TH1D("histoFWHMMesonAlpha01","",fNBinsPeakPt,fBinsPeakPtHalf); }
	fHistoMassMeson = 			new TH1D("histoMassMeson","",fNBinsPt,fBinsPt);
	fHistoWidthMeson = 			new TH1D("histoWidthMeson","",fNBinsPt,fBinsPt);
	fHistoFWHMMeson = 			new TH1D("histoFWHMMeson","",fNBinsPt,fBinsPt);
	fHistoTrueMassMeson = 		new TH1D("histoTrueMassMeson","",fNBinsPt,fBinsPt);
	fHistoTrueMassMesonReweighted =      new TH1D("histoTrueMassMesonReweighted","",fNBinsPt,fBinsPt);
	fHistoTrueFWHMMeson = 		new TH1D("histoTrueFWHMMeson","",fNBinsPt,fBinsPt);
	fHistoTrueFWHMMesonReweighted =      new TH1D("histoTrueFWHMMesonReweighted","",fNBinsPt,fBinsPt);
	fHistoTrueSignMeson = 		new TH1D("histoTrueSignMeson","",fNBinsPt,fBinsPt);
	fHistoTrueSBMeson = 		new TH1D("histoTrueSBMeson","",fNBinsPt,fBinsPt);
	if (fAdvancedMesonQA && (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5)){
		fHistoTrueMassMesonCaloConvPhoton = 		new TH1D("histoTrueMassMesonCaloConvPhoton","",fNBinsPt,fBinsPt);
		fHistoTrueMassMesonCaloElectron = 		new TH1D("histoTrueMassMesonCaloConvElectron","",fNBinsPt,fBinsPt);
		fHistoTrueMassMesonCaloEMNonLeading = 		new TH1D("histoTrueMassMesonCaloEMNonLeading","",fNBinsPt,fBinsPt);
		fHistoTrueMassMesonCaloMergedCluster = 		new TH1D("histoTrueMassMesonCaloMergedCluster","",fNBinsPt,fBinsPt);
		fHistoTrueMassMesonCaloMergedPartConvCluster = 		new TH1D("histoTrueMassMesonCaloMergedPartConvCluster","",fNBinsPt,fBinsPt);
		fHistoTrueMassMesonCaloPhoton = 		new TH1D("histoTrueMassMesonCaloPhoton","",fNBinsPt,fBinsPt);
		fHistoTrueFWHMMesonCaloConvPhoton = 		new TH1D("histoTrueFWHMMesonCaloConvPhoton","",fNBinsPt,fBinsPt);
		fHistoTrueFWHMMesonCaloElectron = 		new TH1D("histoTrueFWHMMesonCaloConvElectron","",fNBinsPt,fBinsPt);
		fHistoTrueFWHMMesonCaloEMNonLeading = 		new TH1D("histoTrueFWHMMesonCaloEMNonLeading","",fNBinsPt,fBinsPt);
		fHistoTrueFWHMMesonCaloMergedCluster = 		new TH1D("histoTrueFWHMMesonCaloMergedCluster","",fNBinsPt,fBinsPt);
		fHistoTrueFWHMMesonCaloMergedPartConvCluster = 		new TH1D("histoTrueFWHMMesonCaloMergedPartConvCluster","",fNBinsPt,fBinsPt);
		fHistoTrueFWHMMesonCaloPhoton = 		new TH1D("histoTrueFWHMMesonCaloPhoton","",fNBinsPt,fBinsPt);
	
	}
	
	fHistoYieldMesonNarrow = 	new TH1D("histoYieldMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMesonNarrow = 	new TH1D("histoYieldTrueMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMesonReweightedNarrow =  new TH1D("histoYieldTrueMesonNarrowReweighted","",fNBinsPt,fBinsPt);
	fHistoYieldTrueSecMesonNarrow = 		new TH1D("histoYieldTrueSecMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldTrueSecFromK0SMesonNarrow = 		new TH1D("histoYieldTrueSecFromK0SMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldMesonPerEventNarrow = new TH1D("histoYieldMesonPerEventNarrow","",fNBinsPt,fBinsPt);
	fHistoSignMesonNarrow = 		new TH1D("histoSignMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoSBMesonNarrow = 		new TH1D("histoSBMesonNarrow","",fNBinsPt,fBinsPt);

	fHistoYieldMesonWide = 		new TH1D("histoYieldMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMesonWide = 	new TH1D("histoYieldTrueMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMesonReweightedWide =    new TH1D("histoYieldTrueMesonWideReweighted","",fNBinsPt,fBinsPt);
	fHistoYieldTrueSecMesonWide = 		new TH1D("histoYieldTrueSecMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldTrueSecFromK0SMesonWide = 		new TH1D("histoYieldTrueSecFromK0SMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldMesonPerEventWide = new TH1D("histoYieldMesonPerEventWide","",fNBinsPt,fBinsPt);
	fHistoSignMesonWide = 		new TH1D("histoSignMesonWide","",fNBinsPt,fBinsPt);
	fHistoSBMesonWide = 		new TH1D("histoSBMesonWide","",fNBinsPt,fBinsPt);

	// Histos for normalization at the left of the peak

	fHistoYieldMesonLeft = 		new TH1D("histoYieldMesonLeft","",fNBinsPt,fBinsPt);
	fHistoYieldMesonLeftPerEvent = new TH1D("histoYieldMesonLeftPerEvent","",fNBinsPt,fBinsPt);
	fHistoSignMesonLeft = 		new TH1D("histoSignMesonLeft","",fNBinsPt,fBinsPt);
	fHistoSBMesonLeft = 		new TH1D("histoSBMesonLeft","",fNBinsPt,fBinsPt);
	fHistoMassMesonLeft = 		new TH1D("histoMassMesonLeft","",fNBinsPt,fBinsPt);
	fHistoWidthMesonLeft = 		new TH1D("histoWidthMesonLeft","",fNBinsPt,fBinsPt);
	fHistoFWHMMesonLeft = 		new TH1D("histoFWHMMesonLeft","",fNBinsPt,fBinsPt);


	fHistoYieldMesonLeftNarrow = 	new TH1D("histoYieldMesonLeftNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldMesonLeftPerEventNarrow = new TH1D("histoYieldMesonLeftPerEventNarrow","",fNBinsPt,fBinsPt);
	fHistoSignMesonLeftNarrow = 	new TH1D("histoSignMesonLeftNarrow","",fNBinsPt,fBinsPt);
	fHistoSBMesonLeftNarrow = 	new TH1D("histoSBMesonLeftNarrow","",fNBinsPt,fBinsPt);


	fHistoYieldMesonLeftWide = 	new TH1D("histoYieldMesonLeftWide","",fNBinsPt,fBinsPt);
	fHistoYieldMesonLeftPerEventWide = new TH1D("histoYieldMesonLeftPerEventWide","",fNBinsPt,fBinsPt);
	fHistoSignMesonLeftWide = 	new TH1D("histoSignMesonLeftWide","",fNBinsPt,fBinsPt);
	fHistoSBMesonLeftWide = 		new TH1D("histoSBMesonLeftWide","",fNBinsPt,fBinsPt);

}

void FillPtHistos()
{
	if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0){
		for (Int_t iPt=fStartPtBin+1;iPt<fNBinsPeakPt+1;iPt++){
			fHistoMassPosition->SetBinContent(iPt,fMesonMassPeakPos[iPt-1]);
			fHistoMassPosition->SetBinError(iPt,fMesonMassPeakPosError[iPt-1]);
			fHistoFWHMMesonAlpha01->SetBinContent(iPt,fMesonFWHMAlpha01[iPt-1]);
			fHistoFWHMMesonAlpha01->SetBinError(iPt,fMesonFWHMAlpha01Error[iPt-1]);
		}
	}
	for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){

		fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
		fDeltaPt->SetBinError(iPt,0);


		fHistoMassMeson->SetBinContent(iPt,fMesonMass[iPt-1]);
		fHistoMassMeson->SetBinError(iPt,fMesonMassError[iPt-1]);
		// fHistoWidthMeson->SetBinContent(iPt);
		fHistoFWHMMeson->SetBinContent(iPt,fMesonFWHM[iPt-1]);
		fHistoFWHMMeson->SetBinError(iPt,fMesonFWHMError[iPt-1]);

		if (fIsMC) {
			fHistoTrueMassMeson->SetBinContent(iPt,fMesonTrueMass[iPt-1]);
			fHistoTrueMassMeson->SetBinError(iPt,fMesonTrueMassError[iPt-1]);
			fHistoTrueMassMesonReweighted->SetBinContent(iPt,fMesonTrueMassReweighted[iPt-1]);
			fHistoTrueMassMesonReweighted->SetBinError(iPt,fMesonTrueMassReweightedError[iPt-1]);
			
			fHistoTrueFWHMMeson->SetBinContent(iPt,fMesonTrueFWHM[iPt-1]);
			fHistoTrueFWHMMeson->SetBinError(iPt,fMesonTrueFWHMError[iPt-1]);
			fHistoTrueFWHMMesonReweighted->SetBinContent(iPt,fMesonTrueFWHMReweighted[iPt-1]);
			fHistoTrueFWHMMesonReweighted->SetBinError(iPt,fMesonTrueFWHMReweightedError[iPt-1]);
			
			fHistoTrueSignMeson->SetBinContent(iPt,fMesonTrueSign[iPt-1]);
			fHistoTrueSignMeson->SetBinError(iPt,fMesonTrueSignError[iPt-1]);
			fHistoTrueSBMeson->SetBinContent(iPt,fMesonTrueSB[iPt-1]);
			fHistoTrueSBMeson->SetBinError(iPt,fMesonTrueSBError[iPt-1]);
			
			if (fAdvancedMesonQA && (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5)){
				fHistoTrueMassMesonCaloConvPhoton->SetBinContent(iPt,fMesonTrueMassCaloConvPhoton[iPt-1]);
				fHistoTrueMassMesonCaloConvPhoton->SetBinError(iPt,fMesonTrueMassErrorCaloConvPhoton[iPt-1]);
				fHistoTrueMassMesonCaloElectron->SetBinContent(iPt,fMesonTrueMassCaloElectron[iPt-1]);
				fHistoTrueMassMesonCaloElectron->SetBinError(iPt,fMesonTrueMassErrorCaloElectron[iPt-1]);
				fHistoTrueMassMesonCaloEMNonLeading->SetBinContent(iPt,fMesonTrueMassCaloEMNonLeading[iPt-1]);
				fHistoTrueMassMesonCaloEMNonLeading->SetBinError(iPt,fMesonTrueMassErrorCaloEMNonLeading[iPt-1]);
				fHistoTrueMassMesonCaloMergedCluster->SetBinContent(iPt,fMesonTrueMassCaloMergedCluster[iPt-1]);
				fHistoTrueMassMesonCaloMergedCluster->SetBinError(iPt,fMesonTrueMassErrorCaloMergedCluster[iPt-1]);
				fHistoTrueMassMesonCaloMergedPartConvCluster->SetBinContent(iPt,fMesonTrueMassCaloMergedClusterPartConv[iPt-1]);
				fHistoTrueMassMesonCaloMergedPartConvCluster->SetBinError(iPt,fMesonTrueMassErrorCaloMergedClusterPartConv[iPt-1]);
				fHistoTrueMassMesonCaloPhoton->SetBinContent(iPt,fMesonTrueMassCaloPhoton[iPt-1]);
				fHistoTrueMassMesonCaloPhoton->SetBinError(iPt,fMesonTrueMassErrorCaloPhoton[iPt-1]);
				fHistoTrueFWHMMesonCaloConvPhoton->SetBinContent(iPt,fMesonTrueFWHMCaloConvPhoton[iPt-1]);
				fHistoTrueFWHMMesonCaloConvPhoton->SetBinError(iPt,fMesonTrueFWHMErrorCaloConvPhoton[iPt-1]);
				fHistoTrueFWHMMesonCaloElectron->SetBinContent(iPt,fMesonTrueFWHMCaloElectron[iPt-1]);
				fHistoTrueFWHMMesonCaloElectron->SetBinError(iPt,fMesonTrueFWHMErrorCaloElectron[iPt-1]);
				fHistoTrueFWHMMesonCaloEMNonLeading->SetBinContent(iPt,fMesonTrueFWHMCaloEMNonLeading[iPt-1]);
				fHistoTrueFWHMMesonCaloEMNonLeading->SetBinError(iPt,fMesonTrueFWHMErrorCaloEMNonLeading[iPt-1]);
				fHistoTrueFWHMMesonCaloMergedCluster->SetBinContent(iPt,fMesonTrueFWHMCaloMergedCluster[iPt-1]);
				fHistoTrueFWHMMesonCaloMergedCluster->SetBinError(iPt,fMesonTrueFWHMErrorCaloMergedCluster[iPt-1]);
				fHistoTrueFWHMMesonCaloMergedPartConvCluster->SetBinContent(iPt,fMesonTrueFWHMCaloMergedClusterPartConv[iPt-1]);
				fHistoTrueFWHMMesonCaloMergedPartConvCluster->SetBinError(iPt,fMesonTrueFWHMErrorCaloMergedClusterPartConv[iPt-1]);
				fHistoTrueFWHMMesonCaloPhoton->SetBinContent(iPt,fMesonTrueFWHMCaloPhoton[iPt-1]);
				fHistoTrueFWHMMesonCaloPhoton->SetBinError(iPt,fMesonTrueFWHMErrorCaloPhoton[iPt-1]);

				
			}	 
		}

		fHistoSignMeson->SetBinContent(iPt,fMesonSign[iPt-1]);
		fHistoSignMeson->SetBinError(iPt,fMesonSignError[iPt-1]);
		fHistoSBMeson->SetBinContent(iPt,fMesonSB[iPt-1]);
		fHistoSBMeson->SetBinError(iPt,fMesonSBError[iPt-1]);

		fHistoYieldMeson->SetBinContent(iPt,fMesonYieldsCorResidualBckFunc[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMeson->SetBinError(iPt,fMesonYieldsCorResidualBckFuncError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		fHistoYieldMesonBackFit->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncBackFit[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonBackFit->SetBinError(iPt,fMesonYieldsCorResidualBckFuncBackFitError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));


		if (fIsMC) {
			fHistoYieldTrueMeson->SetBinContent(iPt,fMesonTrueYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMeson->SetBinError(iPt,fMesonTrueYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonReweighted->SetBinContent(iPt,fMesonTrueYieldsReweighted[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonReweighted->SetBinError(iPt,fMesonTrueYieldsReweightedError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecMeson->SetBinContent(iPt,fMesonTrueSecYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecMeson->SetBinError(iPt,fMesonTrueSecYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecFromK0SMeson->SetBinContent(iPt,fMesonTrueSecFromK0SYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecFromK0SMeson->SetBinError(iPt,fMesonTrueSecFromK0SYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		
		}
		fHistoYieldMesonPerEvent->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFunc[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEvent->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		fHistoYieldMesonPerEventBackFit->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncBackFit[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEventBackFit->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncBackFitError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Narrow integration window
		fHistoSignMesonNarrow->SetBinContent(iPt,fMesonSignNarrow[iPt-1]);
		fHistoSignMesonNarrow->SetBinError(iPt,fMesonSignNarrowError[iPt-1]);
		fHistoSBMesonNarrow->SetBinContent(iPt,fMesonSBNarrow[iPt-1]);
		fHistoSBMesonNarrow->SetBinError(iPt,fMesonSBNarrowError[iPt-1]);

		fHistoYieldMesonNarrow->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonNarrow->SetBinError(iPt,fMesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		if (fIsMC) {
			fHistoYieldTrueMesonNarrow->SetBinContent(iPt,fMesonTrueYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonNarrow->SetBinError(iPt,fMesonTrueYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonReweightedNarrow->SetBinContent(iPt,fMesonTrueYieldsReweightedNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonReweightedNarrow->SetBinError(iPt,fMesonTrueYieldsReweightedNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecMesonNarrow->SetBinContent(iPt,fMesonTrueSecYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecMesonNarrow->SetBinError(iPt,fMesonTrueSecYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecFromK0SMesonNarrow->SetBinContent(iPt,fMesonTrueSecFromK0SYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecFromK0SMesonNarrow->SetBinError(iPt,fMesonTrueSecFromK0SYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		
		}
		fHistoYieldMesonPerEventNarrow->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEventNarrow->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Wide integration window
		fHistoSignMesonWide->SetBinContent(iPt,fMesonSignWide[iPt-1]);
		fHistoSignMesonWide->SetBinError(iPt,fMesonSignWideError[iPt-1]);
		fHistoSBMesonWide->SetBinContent(iPt,fMesonSBWide[iPt-1]);
		fHistoSBMesonWide->SetBinError(iPt,fMesonSBWideError[iPt-1]);

		fHistoYieldMesonWide->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonWide->SetBinError(iPt,fMesonYieldsCorResidualBckFuncWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		if (fIsMC) {
			fHistoYieldTrueMesonWide->SetBinContent(iPt,fMesonTrueYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonWide->SetBinError(iPt,fMesonTrueYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonReweightedWide->SetBinContent(iPt,fMesonTrueYieldsReweightedWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonReweightedWide->SetBinError(iPt,fMesonTrueYieldsReweightedWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecMesonWide->SetBinContent(iPt,fMesonTrueSecYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecMesonWide->SetBinError(iPt,fMesonTrueSecYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecFromK0SMesonWide->SetBinContent(iPt,fMesonTrueSecFromK0SYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueSecFromK0SMesonWide->SetBinError(iPt,fMesonTrueSecFromK0SYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		
		}
		fHistoYieldMesonPerEventWide->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEventWide->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Histos for integration at the left of the peak
		fHistoMassMesonLeft->SetBinContent(iPt,fMesonMassLeft[iPt-1]);
		fHistoMassMesonLeft->SetBinError(iPt,fMesonMassLeftError[iPt-1]);
		fHistoFWHMMesonLeft->SetBinContent(iPt,fMesonFWHMLeft[iPt-1]);
		fHistoFWHMMesonLeft->SetBinError(iPt,fMesonFWHMLeftError[iPt-1]);

		fHistoSignMesonLeft->SetBinContent(iPt,fMesonSignLeft[iPt-1]);
		fHistoSignMesonLeft->SetBinError(iPt,fMesonSignLeftError[iPt-1]);
		fHistoSBMesonLeft->SetBinContent(iPt,fMesonSBLeft[iPt-1]);
		fHistoSBMesonLeft->SetBinError(iPt,fMesonSBLeftError[iPt-1]);

		fHistoYieldMesonLeft->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeft[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeft->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEvent->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeft[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEvent->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Narrow integration window
		fHistoSignMesonLeftNarrow->SetBinContent(iPt,fMesonSignLeftNarrow[iPt-1]);
		fHistoSignMesonLeftNarrow->SetBinError(iPt,fMesonSignLeftNarrowError[iPt-1]);
		fHistoSBMesonLeftNarrow->SetBinContent(iPt,fMesonSBLeftNarrow[iPt-1]);
		fHistoSBMesonLeftNarrow->SetBinError(iPt,fMesonSBLeftNarrowError[iPt-1]);

		fHistoYieldMesonLeftNarrow->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftNarrow->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventNarrow->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventNarrow->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Wide integration window
		fHistoSignMesonLeftWide->SetBinContent(iPt,fMesonSignLeftWide[iPt-1]);
		fHistoSignMesonLeftWide->SetBinError(iPt,fMesonSignLeftWideError[iPt-1]);
		fHistoSBMesonLeftWide->SetBinContent(iPt,fMesonSBLeftWide[iPt-1]);
		fHistoSBMesonLeftWide->SetBinError(iPt,fMesonSBLeftWideError[iPt-1]);

		fHistoYieldMesonLeftWide->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftWide->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventWide->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventWide->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
	}
}





void FitSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntRangeFit, Int_t ptBin, Bool_t vary)
{

	//    cout<<"Start Fitting spectra"<<endl;
	fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
	Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
	Double_t mesonAmplitudeMin;
	Double_t mesonAmplitudeMax;
	if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
		mesonAmplitudeMin = mesonAmplitude*98./100.;
		mesonAmplitudeMax = mesonAmplitude*115./100.;
		if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) mesonAmplitudeMin = mesonAmplitude*92./100.;
	} else {
		mesonAmplitudeMin = mesonAmplitude*50./100.;
		mesonAmplitudeMax = mesonAmplitude*120./100.;
		if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5){
			mesonAmplitudeMin = mesonAmplitude*10./100.;
		}	
	}
	fFitReco= NULL;
	fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",fMesonFitRange[0],fMesonFitRange[1]);

	fFitGausExp =NULL;
	fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

	fFitLinearBck = NULL;
	fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);


	fFitReco->SetParameter(0,mesonAmplitude);
	fFitReco->SetParameter(1,fMesonMassExpect);
	fFitReco->SetParameter(2,fMesonWidthExpect);
	//    if(vary){
		fFitReco->SetParameter(3,fMesonLambdaTail);
	//    } else {
	//       fFitReco->FixParameter(3,fMesonLambdaTail);
	//    }
	fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
	//    fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
	fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
	fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
	//    if(vary){
		fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
	// }

	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

	
	fFitReco->SetLineColor(3);
	fFitReco->SetLineWidth(1);
	fFitReco->SetLineStyle(1);

	if (vary && !fIsMC){
		fMesonLambdaTail = fFitReco->GetParameter(3);
		fMesonLambdaTailRange[0] = 0.9*fFitReco->GetParameter(3);
		fMesonLambdaTailRange[1] = 1.1*fFitReco->GetParameter(3);
		fMesonWidthExpect = fFitReco->GetParameter(2);
		fMesonWidthRange[0] = 0.5*fFitReco->GetParameter(2);
		fMesonWidthRange[1] = 1.5*fFitReco->GetParameter(2);
	}
	fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
	fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
	fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
	fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));

	fFitGausExp->SetParError(0,fFitReco->GetParError(0));
	fFitGausExp->SetParError(1,fFitReco->GetParError(1));
	fFitGausExp->SetParError(2,fFitReco->GetParError(2));
	fFitGausExp->SetParError(3,fFitReco->GetParError(3));

	fFitLinearBck->SetParameter(0,fFitReco->GetParameter(4));
	fFitLinearBck->SetParameter(1,fFitReco->GetParameter(5));

	fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
	fFitLinearBck->SetParError(1,fFitReco->GetParError(5));

	Int_t binCenterStart;
	Double_t startBinEdge;
	Int_t binCenterEnd;
	Double_t endBinEdge;

	TVirtualFitter * fitter = TVirtualFitter::GetFitter();

	fIntLinearBck = 0;
	fIntLinearBckError = 0;
	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
		binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[0]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[1]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

		Int_t nFreePar = fFitReco->GetNumberFreeParameters();
		double * covMatrix = fitter->GetCovarianceMatrix();

		Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
			0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);

		Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitReco->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitReco->GetParError(5),2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

		fFileDataLog << "Parameter for bin " << ptBin << endl;
		fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
		fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5)<< endl;

		fIntLinearBck = intLinearBack/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		fIntLinearBckError = errorLinearBck/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
	} else {
		fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
	}
	fFitReco->DrawCopy("same");
	
}

void FitTrueInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntRangeFit, Int_t ptBin, Bool_t vary)
{
	//    cout<<"Start Fitting spectra"<<endl;
	fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
	Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
	Double_t mesonAmplitudeMin;
	Double_t mesonAmplitudeMax;
	mesonAmplitudeMin = mesonAmplitude*95./100.;
	mesonAmplitudeMax = mesonAmplitude*130./100.;
	if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5){
		mesonAmplitudeMin = mesonAmplitude*20./100.;
	}
	fFitReco = NULL;
	fFitReco = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

	fFitReco->SetParameter(0,mesonAmplitude);
	fFitReco->SetParameter(1,fMesonMassExpect);
	fFitReco->SetParameter(2,fMesonWidthExpect);
	fFitReco->SetParameter(3,fMesonLambdaTailMC);

	fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
	fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
	fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
	fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
	// //    fFitReco->SetParLimits(1,fMesonMassExpect*0.8,fMesonMassExpect*1.3);
	// //    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
	// //    if(vary){
	// 	   fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
	
	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"RME0");

	//    cout << TString(gMinuit->fCstatu.Data()).Data() << endl;

	if (vary){
		fMesonLambdaTailMC = fFitReco->GetParameter(3);
	// 	   fMesonWidthExpectMC = fMesonMassExpect*0.03;
		}
	fFitReco->SetLineColor(3);
	fFitReco->SetLineWidth(1);
	fFitReco->SetLineStyle(1);

	//    Int_t binCenterStart = 0;
	//    Double_t startBinEdge = 0;;
	//    Int_t binCenterEnd = 0;
	//    Double_t endBinEdge = 0;

	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
	//       binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[0]-(fMesonMassExpect-fFitReco->GetParameter(1)));
	//       startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
	//       binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[1]-(fMesonMassExpect-fFitReco->GetParameter(1)));
	//       endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

		fFileDataLog << "Parameter for bin " << ptBin << endl;
		fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
	} else {
		fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
	}
	fFitReco->DrawCopy("same");
	if (fMesonIntRangeFit){}
}


void FitPeakPosInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Int_t ptBin, Bool_t vary)
{

   //    cout<<"Start Fitting spectra"<<endl;
   fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(0.05,0.3);
   Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
   Double_t mesonAmplitudeMin;
   Double_t mesonAmplitudeMax;
   if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
      mesonAmplitudeMin = mesonAmplitude*80./100.;
      mesonAmplitudeMax = mesonAmplitude*115./100.;
   } else {
      mesonAmplitudeMin = mesonAmplitude*80./110.;
      mesonAmplitudeMax = mesonAmplitude*115./100.;
   }
   fFitReco = NULL;
   fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x+[6]*x*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x+[6]*x*x)",0.05,0.3);

   fFitGausExp = NULL;
   fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",0.05,0.3);

   fFitLinearBck = NULL;
   fFitLinearBck = new TF1("Linear","[0]+[1]*x+[2]*x*x",0.05,0.3);


   fFitReco->SetParameter(0,mesonAmplitude);
   fFitReco->SetParameter(1,fMesonMassExpect);
   fFitReco->SetParameter(2,fMesonWidthExpect);
   if(vary){
      fFitReco->SetParameter(3,fMesonLambdaTail);
   } else {
      fFitReco->FixParameter(3,fMesonLambdaTail);
   }
   fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
   fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
   fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
   if(vary){fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);}

   fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
   fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

   fFitReco->SetLineColor(3);
   fFitReco->SetLineWidth(1);
   fFitReco->SetLineStyle(1);

   //if (vary) fMesonLambdaTail = fFitReco->GetParameter(3);

   fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
   fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
   fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
   fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));

   fFitGausExp->SetParError(0,fFitReco->GetParError(0));
   fFitGausExp->SetParError(1,fFitReco->GetParError(1));
   fFitGausExp->SetParError(2,fFitReco->GetParError(2));
   fFitGausExp->SetParError(3,fFitReco->GetParError(3));

   fFitLinearBck->SetParameter(0,fFitReco->GetParameter(4));
   fFitLinearBck->SetParameter(1,fFitReco->GetParameter(5));
   fFitLinearBck->SetParameter(1,fFitReco->GetParameter(6));

   fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
   fFitLinearBck->SetParError(1,fFitReco->GetParError(5));
   fFitLinearBck->SetParError(1,fFitReco->GetParError(6));

   if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
      fFileDataLog << "Parameter for bin " << ptBin << endl;
      fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
      fFileDataLog << "Quadratic: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5) << "\t "<<fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<< endl;
   } else {
      fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
   }
   fFitReco->DrawCopy("same");
}


void FitCBSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntRangeFit, Int_t ptBin,Bool_t vary ,TString functionname)
{


   fFileErrLog<<"Start Fitting spectra with CB fit"<<endl;
   fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
   Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
   Double_t mesonAmplitudeMin = mesonAmplitude*50./100.;
   Double_t mesonAmplitudeMax = mesonAmplitude*115./100.;

   fFitReco = NULL;
   fFitReco = new TF1(functionname,CrystalBallBck,fMesonFitRange[0],fMesonFitRange[1],7);

   fFitGausExp = NULL;
   fFitGausExp = new TF1("optionCrystalBall",CrystalBall,fMesonFitRange[0],fMesonFitRange[1],5);

   fFitLinearBck = NULL;
   fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);

   fFitReco->SetParameter(0,mesonAmplitude);
   fFitReco->SetParameter(1,fMesonMassExpect);
   fFitReco->SetParameter(2,fMesonWidthExpect);
   fFitReco->SetParameter(3,2.);
   fFitReco->SetParameter(4,0.7);
   fFitReco->SetParameter(5,0.);
   fFitReco->SetParameter(6,1.);

   if (!vary) {
      fFitReco->FixParameter(3,fCBn);
      fFitReco->FixParameter(4,fCBAlpha);
   }

   fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
   fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
   fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);

   fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRE0");
   fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRE0");

   fFitReco->SetLineColor(3);
   fFitReco->SetLineWidth(1);
   fFitReco->SetLineStyle(1);

   if (vary) {
      fCBAlpha = fFitReco->GetParameter(4);
      fCBn = fFitReco->GetParError(3);
   }

   fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
   fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
   fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
   fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));
   fFitGausExp->SetParameter(4,fFitReco->GetParameter(4));

   fFitGausExp->SetParError(0,fFitReco->GetParError(0));
   fFitGausExp->SetParError(1,fFitReco->GetParError(1));
   fFitGausExp->SetParError(2,fFitReco->GetParError(2));
   fFitGausExp->SetParError(3,fFitReco->GetParError(3));
   fFitGausExp->SetParError(4,fFitReco->GetParError(4));

   fFitLinearBck->SetParameter(0,fFitReco->GetParameter(5));
   fFitLinearBck->SetParameter(1,fFitReco->GetParameter(6));

   fFitLinearBck->SetParError(0,fFitReco->GetParError(5));
   fFitLinearBck->SetParError(1,fFitReco->GetParError(6));

   Int_t binCenterStart;
   Double_t startBinEdge;
   Int_t binCenterEnd;
   Double_t endBinEdge;

   TVirtualFitter * fitter = TVirtualFitter::GetFitter();

   if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
      binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[0]-(fMesonMassExpect-fFitReco->GetParameter(1)));
      startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
      binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[1]-(fMesonMassExpect-fFitReco->GetParameter(1)));
      endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

      Int_t nFreePar = fFitReco->GetNumberFreeParameters();
      double * covMatrix = fitter->GetCovarianceMatrix();

      Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
         0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);

      Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitReco->GetParError(5),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitReco->GetParError(6),2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

      fFileDataLog << "Parameter for bin " << ptBin << endl;
      fFileDataLog << "CrystalBall: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<< "\t "<< fFitReco->GetParameter(4) <<"+-" << fFitReco->GetParError(4)<<endl;
      fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(5)<<"+-" << fFitReco->GetParError(5) << "\t "<<fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<< endl;

      fIntLinearBck = intLinearBack/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
      fIntLinearBckError = errorLinearBck/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
   } else {
      fFileErrLog << "Fitting failed in " << ptBin << " with status::" << gMinuit->fCstatu.Data() <<"why failed?"<<endl << endl;
   }
   fFitReco->DrawCopy("same");

}

void FitWithPol2ForBG(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntRangeFit, Int_t ptBin, Bool_t vary)
{
   //    cout<<"Start Fitting spectra"<<endl;
   Int_t startBinSearch = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[0]);
   Int_t endBinSearch = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[1]);
   Double_t mesonAmplitude = 0;
   for (Int_t i = startBinSearch; i < endBinSearch+1; i++){
      if (mesonAmplitude < fHistoMappingSignalInvMassPtBinSingle->GetBinContent(i)){
         mesonAmplitude = fHistoMappingSignalInvMassPtBinSingle->GetBinContent(i);
      }
   }
   mesonAmplitude = mesonAmplitude-fHistoMappingSignalInvMassPtBinSingle->GetBinContent(endBinSearch);

   fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0] ,fMesonMassRange[1]);
   Double_t mesonAmplitudeMin;
   Double_t mesonAmplitudeMax;
   if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
      mesonAmplitudeMin = mesonAmplitude*80./100.;
      mesonAmplitudeMax = mesonAmplitude*115./100.;
   } else {
      mesonAmplitudeMin = mesonAmplitude*70./100.;
      mesonAmplitudeMax = mesonAmplitude*110./100.;
   }
   Double_t fitRange[2];
   if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
      fitRange[0] = 	fMesonMassRange[0];
      fitRange[1] = 	fMesonMassRange[1];
   } else {
      fitRange[0] = 	0.3;
      fitRange[1] = 	0.8;
   }
   
   if (fPrefix.CompareTo("Omega") ==0) {fitRange[0]=0.4; fitRange[1]=0.9;}

   fFitReco = NULL;
   fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x+[6]*x*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x+[6]*x*x)",fitRange[0] ,fitRange[1]);

   fFitGausExp = NULL;
   fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",0.05,0.3);

   fFitLinearBck = NULL;
   fFitLinearBck = new TF1("Linear","[0]+[1]*x+[2]*x*x",fitRange[0] ,fitRange[1]);

   fFitReco->SetParameter(0,mesonAmplitude);
   fFitReco->SetParameter(1,fMesonMassExpect);
   fFitReco->SetParameter(2,fMesonWidthExpect);
   if(vary){
      fFitReco->SetParameter(3,fMesonLambdaTail);
   } else {
      fFitReco->FixParameter(3,fMesonLambdaTail);
   }
   fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
   fFitReco->SetParLimits(1,fMesonMassExpect*95./100.,fMesonMassExpect*105./100.);
   fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
   if(vary){fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);}

   fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
   fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

   fFitReco->SetLineColor(3);
   fFitReco->SetLineWidth(1);
   fFitReco->SetLineStyle(1);

   //if (vary) fMesonLambdaTail = fFitReco->GetParameter(3);

   fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
   fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
   fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
   fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));

   fFitGausExp->SetParError(0,fFitReco->GetParError(0));
   fFitGausExp->SetParError(1,fFitReco->GetParError(1));
   fFitGausExp->SetParError(2,fFitReco->GetParError(2));
   fFitGausExp->SetParError(3,fFitReco->GetParError(3));

   fFitLinearBck->SetParameter(0,fFitReco->GetParameter(4));
   fFitLinearBck->SetParameter(1,fFitReco->GetParameter(5));
   fFitLinearBck->SetParameter(2,fFitReco->GetParameter(6));

   fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
   fFitLinearBck->SetParError(1,fFitReco->GetParError(5));
   fFitLinearBck->SetParError(2,fFitReco->GetParError(6));

   if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
      fFileDataLog << "Parameter for bin " << ptBin << endl;
      fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
      fFileDataLog << "Quadratic: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5) << "\t "<<fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<< endl;
   } else {
      fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
   }
   fFitReco->DrawCopy("same");
}




void IntegrateHistoInvMass(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt)
{
   Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
   Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
   fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
}


void IntegrateHistoInvMassStream(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt)
{
   Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
   Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
   fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
   for ( Int_t M = binLowMassMeson; M < binHighMassMeson+1; M++){
      fFileDataLog << M << "\t" << fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(M) <<"\t" <<fHistoMappingSignalInvMassPtBinSingle->GetBinContent(M)<< "+-"<< fHistoMappingSignalInvMassPtBinSingle->GetBinError(M)<< endl;
   }
}


void IntegrateFitFunc(TF1 * fFunc, TH1D *  fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntRangeInt)
{

   fYieldsFunc = fFunc->Integral(fMesonIntRangeInt[0],fMesonIntRangeInt[1])/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

}



void FillHistosArrayMC(TH1D* fHistoMCMesonPtWithinAcceptanceFill, TH1D * fHistoMCMesonPtFill, TH1D * fDeltaPtFill, TString nameCanvas)
{
//    Char_t nameHisto[100] = "fHistoMCMesonPtEtaWithinAcceptance";
	fHistoMCMesonPtWithinAcceptanceFill->Sumw2();
   fHistoMCMesonWithinAccepPt = (TH1D*)fHistoMCMesonPtWithinAcceptanceFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
   fHistoMCMesonWithinAccepPt->Divide(fDeltaPtFill);
	fHistoMCMesonPtFill->Sumw2();
   fHistoMCMesonPt1 = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
   fHistoMCMesonPt1->Divide(fDeltaPtFill);

   
//    TH1D* fHistoMCMesonPt1ScaledByPt = (TH1D*)fHistoMCMesonPtFill->Clone("ScaledByPt");
//    fHistoMCMesonPt1ScaledByPt->Rebin(1);
// // 	fHistoMCMesonPt1ScaledByPt = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
//    for (Int_t i = 1; i < fHistoMCMesonPt1ScaledByPt->GetNbinsX()+1; i++){
//       fHistoMCMesonPt1ScaledByPt->SetBinContent(i, fHistoMCMesonPt1ScaledByPt->GetBinContent(i)/fHistoMCMesonPt1ScaledByPt->GetBinWidth(i)/fNEvents/2/TMath::Pi()/fHistoMCMesonPt1ScaledByPt->GetBinCenter(i)/(2*fYMaxMeson));//
//       fHistoMCMesonPt1ScaledByPt->SetBinError(i, fHistoMCMesonPt1ScaledByPt->GetBinError(i)/fHistoMCMesonPt1ScaledByPt->GetBinWidth(i)/2/TMath::Pi()/fNEvents/fHistoMCMesonPt1ScaledByPt->GetBinCenter(i)/(2*fYMaxMeson));//fHistoMCMesonPt1ScaledByPt->GetBinCenter(i)/
// //       		fHistoMCMesonPt1ScaledByPt->SetBinContent(i, fHistoMCMesonPt1ScaledByPt->GetBinContent(i)/fNEvents/fHistoMCMesonPt1ScaledByPt->GetBinWidth(i));// 
// //        		fHistoMCMesonPt1ScaledByPt->SetBinError(i, fHistoMCMesonPt1ScaledByPt->GetBinError(i)/fNEvents/fHistoMCMesonPt1ScaledByPt->GetBinWidth(i));
//    }
//    TCanvas* canvasRecTrueFit = new TCanvas();
//    canvasRecTrueFit->SetLogy(1);
// // 	canvasRecTrueFit->SetLogx(1);
//    // 	Double_t paramTest[3] = {1.0e-2,5.,0.13};<<"-" << fBinsPt[iPt+1]
// 	Double_t maximumPt = fHistoMCMesonPt1ScaledByPt->GetXaxis()->GetBinUpEdge(fHistoMCMesonPt1ScaledByPt->GetNbinsX());
// 	for( Int_t i = fHistoMCMesonPt1ScaledByPt->GetNbinsX(); i > 1 ; i--){
// 		if (	fHistoMCMesonPt1ScaledByPt->GetBinContent(i) != 0){
// 			maximumPt = fHistoMCMesonPt1ScaledByPt->GetXaxis()->GetBinUpEdge(i);
// 			i = 0;
// 		}
// 	}	
// //    maximumPt= 10;
//    Double_t parameters[3] = {2.,  6.18,1.35};
// //    TF1* fitInputMC = FitObject("h","fitInputMC","Pi0",fHistoMCMesonPt1ScaledByPt,0.1 ,maximumPt,parameters, "NRME+");//,paramTest,"NRME+");
//    fitInputMC->FixParameter(1,6.18);
//    fitInputMC->FixParameter(2,1.35);
//    TF1* fitInputMC = FitObject("l","fitInputMC","Pi0",fHistoMCMesonPt1ScaledByPt,0.1 ,maximumPt,NULL, "NRME+");//,paramTest,"NRME+");
//    TFitResultPtr resultInputMC = fHistoMCMesonPt1ScaledByPt->Fit(fitInputMC,"SNRME+","",0.2 ,maximumPt);
//    DrawGammaSetMarkerTF1(fitInputMC, 1, 1.5, kBlue);
// 	
//    fHistoMCMesonPt1ScaledByPt->GetYaxis()->SetRangeUser(1e-9,1e5);
//    fHistoMCMesonPt1ScaledByPt->Draw("pe1");
//    cout << WriteParameterToFile(fitInputMC)<< endl;	
// 	cout << fitInputMC->Integral(0,25,resultInputMC->GetParams())<< endl;
// //    fitInputMC->SetRange(0.,25.);
//    fitInputMC->Draw("same,l");
//    canvasRecTrueFit->SaveAs(nameCanvas.Data());
	if (nameCanvas){}	
}



void CalculateMesonAcceptance()
{
   fHistoMCMesonAcceptPt = new TH1D("fMCMesonAccepPt","",fNBinsPt,fBinsPt);
   fHistoMCMesonAcceptPt->Sumw2();

   fHistoMCMesonAcceptPt->Divide(fHistoMCMesonWithinAccepPt,fHistoMCMesonPt1,1.,1.,"B");
   fHistoMCMesonAcceptPt->DrawCopy();
   fFileDataLog << endl << "Calculation of the Acceptance" << endl;
   for ( Int_t i = 1; i < fHistoMCMesonAcceptPt->GetNbinsX()+1 ; i++){
      fFileDataLog << "Bin " << i << "\t"<< fHistoMCMesonAcceptPt->GetBinCenter(i)<< "\t" << fHistoMCMesonAcceptPt->GetBinContent(i) << "\t" << fHistoMCMesonAcceptPt->GetBinError(i) <<endl;
   }
}

void CalculateMesonEfficiency(TH1D* fMC_fMesonYieldsPt, TH1D* fMC_SecondaryYieldPt, TString nameEfi, TString nameCanvas )
{
   fHistoMCMesonEffiPt = new TH1D(nameEfi.Data(),"",fNBinsPt,fBinsPt);
	
   fHistoMCMesonEffiPt->Sumw2();
   fHistoMCMesonEffiPt->Add(fMC_fMesonYieldsPt,1.);
   fHistoMCMesonEffiPt->Add(fMC_SecondaryYieldPt,-1.);
   //cout<<"ADDING 10"<<endl;
	
	
   TH1D* fHistoMCMesonEffiScaledByPt = (TH1D*)fHistoMCMesonEffiPt->Clone("ScaledByPt");
   for (Int_t i = 1; i < fHistoMCMesonEffiScaledByPt->GetNbinsX()+1; i++){
      fHistoMCMesonEffiScaledByPt->SetBinContent(i, fHistoMCMesonEffiScaledByPt->GetBinContent(i)/fNEvents);//fHistoMCMesonEffiScaledByPt->GetBinCenter(i)/2/TMath::Pi()
      fHistoMCMesonEffiScaledByPt->SetBinError(i, fHistoMCMesonEffiScaledByPt->GetBinError(i)/fNEvents	);//fHistoMCMesonEffiScaledByPt->GetBinCenter(i)/2/TMath::Pi()
   }
	
   fHistoMCMesonEffiPt->Divide(fHistoMCMesonEffiPt,fHistoMCMesonWithinAccepPt,1.,1.,"B");
   fFileDataLog << endl << "Calculation of the Efficiency" << nameEfi.Data()<< endl;
   for ( Int_t i = 1; i < fHistoMCMesonEffiPt->GetNbinsX()+1 ; i++){
      fFileDataLog << "Bin " << i << "\t" << fHistoMCMesonEffiPt->GetBinCenter(i)<< "\t"<< fHistoMCMesonEffiPt->GetBinContent(i) << "\t" << fHistoMCMesonEffiPt->GetBinError(i) <<endl;
   }
   fHistoMCMesonEffiFitPt = (TH1D*)fHistoMCMesonEffiPt->Clone(Form("%s_Fit",nameEfi.Data()));
//    if (fPrefix.CompareTo("Pi0") == 0 ||
//        (fPrefix.CompareTo("Eta") == 0 && (fEnergyFlag.CompareTo("7TeV") == 0 || fEnergyFlag.CompareTo("2.76TeV") == 0) )||
//        (fPrefix.CompareTo("Pi0EtaBinning") == 0 && (fEnergyFlag.CompareTo("7TeV") == 0 || fEnergyFlag.CompareTo("2.76TeV") == 0))){
//       TCanvas* canvasRecTrueFit = new TCanvas();
//       canvasRecTrueFit->SetLogy(1);
//       //   	Double_t paramTest[3] = {24,6.7,-6.5};
//       TF1* fitPtQCD = FitObject("qcd","fitPtQCD","Pi0",fHistoMCMesonEffiScaledByPt,fBinsPt[fStartPtBin] ,fBinsPt[fNBinsPt], NULL,"NRME+");//,paramTest,"NRME+");
//       cout << WriteParameterToFile(fitPtQCD)<< endl;
//       TFitResultPtr resultQCD = fHistoMCMesonEffiScaledByPt->Fit(fitPtQCD,"SNRME+","",fBinsPt[fStartPtBin] ,fBinsPt[fNBinsPt]);
//
//       DrawGammaSetMarkerTF1(fitPtQCD, 1, 1.5, kBlue);
//       fHistoMCMesonEffiScaledByPt->GetYaxis()->SetRangeUser(1e-9,1e2);
//       fHistoMCMesonEffiScaledByPt->Draw("pe1");
//       fitPtQCD->Draw("same,l");
//
//       TH1D* fHistoMCMesonWithinAccepPtScaledByPt = (TH1D*)fHistoMCMesonWithinAccepPt->Clone("ScaledByPtWithinAcceptance");
//       for (Int_t i = 1; i < fHistoMCMesonWithinAccepPtScaledByPt->GetNbinsX()+1; i++){
//          // 		fHistoMCMesonWithinAccepPtScaledByPt->SetBinContent(i, fHistoMCMesonWithinAccepPtScaledByPt->GetBinContent(i)/fHistoMCMesonWithinAccepPtScaledByPt->GetBinCenter(i)/fNEvents/2/TMath::Pi());
//          //fHistoMCMesonWithinAccepPtScaledByPt->SetBinError(i, fHistoMCMesonWithinAccepPtScaledByPt->GetBinError(i)/fHistoMCMesonWithinAccepPtScaledByPt->GetBinCenter(i)/2/TMath::Pi()/fNEvents	);
//          fHistoMCMesonWithinAccepPtScaledByPt->SetBinContent(i, fHistoMCMesonWithinAccepPtScaledByPt->GetBinContent(i)/fNEvents);
//          fHistoMCMesonWithinAccepPtScaledByPt->SetBinError(i, fHistoMCMesonWithinAccepPtScaledByPt->GetBinError(i)/fNEvents	);
//       }
//       Double_t paramTest2[3] = {1.0e5,7.,0.13};
//       TF1* fitPtMCWithinAcceptance = FitObject("qcd","fitPtMCWithinAcceptance","Pi0",fHistoMCMesonWithinAccepPtScaledByPt,fBinsPt[fStartPtBin] ,fBinsPt[fNBinsPt],NULL,"NRME+");
//       TFitResultPtr resultMCWithinAcceptance= fHistoMCMesonWithinAccepPtScaledByPt->Fit(fitPtMCWithinAcceptance,"SQNRME+","",fBinsPt[fStartPtBin] ,fBinsPt[fNBinsPt]);
//
//       DrawGammaSetMarkerTF1(fitPtMCWithinAcceptance, 1, 1.5, kGray);
//       fHistoMCMesonWithinAccepPtScaledByPt->SetLineColor(kGray);
//       fHistoMCMesonWithinAccepPtScaledByPt->SetMarkerColor(kGray);
//       // // 	fHistoMCMesonWithinAccepPtScaledByPt->GetYaxis()->SetRangeUser(1e-6,100);
//       fHistoMCMesonWithinAccepPtScaledByPt->Draw("pe1,same");
//       cout << WriteParameterToFile(fitPtMCWithinAcceptance)<< endl;
//       // 	fitPtMCWithinAcceptance->SetRange(0.,10.);
//       fitPtMCWithinAcceptance->Draw("same,l");
//       canvasRecTrueFit->SaveAs(nameCanvas.Data());
//
//       for (Int_t i = 2; i < fHistoMCMesonEffiFitPt->GetNbinsX()+1; i++){
//          Double_t ptStart = fHistoMCMesonEffiFitPt->GetXaxis()->GetBinLowEdge(i);
//          Double_t ptEnd = fHistoMCMesonEffiFitPt->GetXaxis()->GetBinUpEdge(i);
//          Double_t binWidth = ptEnd-ptStart;
//          Double_t intRec;
//          Double_t errorIntRec;
//          intRec = fitPtQCD->Integral(ptStart, ptEnd, resultQCD->GetParams());
//          errorIntRec = fitPtQCD->IntegralError(ptStart, ptEnd,resultQCD->GetParams(),resultQCD->GetCovarianceMatrix().GetMatrixArray() );
//
//          Double_t intProd;
//          Double_t errorIntProd ;
//          intProd = fitPtMCWithinAcceptance->Integral(ptStart, ptEnd, resultMCWithinAcceptance->GetParams());
//          errorIntProd = fitPtMCWithinAcceptance->IntegralError(ptStart, ptEnd, resultMCWithinAcceptance->GetParams(), resultMCWithinAcceptance->GetCovarianceMatrix().GetMatrixArray() );
//          fHistoMCMesonEffiFitPt->SetBinContent(i,intRec/ intProd);
//          Double_t errorRatio = TMath::Sqrt(TMath::Power(errorIntRec/intProd,2) +TMath::Power(errorIntProd/intProd,2))*intRec/intProd;
// 			if (errorRatio == 0){
// 				errorRatio = 0.04*intRec/ intProd;
// 			}
// // 			cout << "Int Rec: " <<intRec << "\t +- "<<errorIntRec/intRec*100 << " \t Int Prod: " << intProd << "\t+- " << errorIntProd/intProd*100 << "\t Effi: "<< intRec/intProd << "\t +-" << errorRatio/(intRec/intProd)*100<<endl;
//          fHistoMCMesonEffiFitPt->SetBinError(i, errorRatio);
//       }
//    }
   if (nameCanvas){}
}



void SaveHistos(Int_t optionMC, TString fCutID, TString fPrefix3)
{
	const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1WithoutCorrection%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix3.Data(),fPeriodFlag.Data(),fCutID.Data());
	fOutput1 = new TFile(nameOutput,"RECREATE");

	cout << "Begin writing Uncorrected File" << endl;
	
	fHistoYieldMeson->Write();
	fHistoYieldMesonPerEvent->Write();
	fHistoYieldMesonBackFit->Write();
	fHistoYieldMesonPerEventBackFit->Write();
	fHistoSignMeson->Write();
	fHistoSBMeson->Write();

	fHistoYieldMesonNarrow->Write();
	fHistoYieldMesonPerEventNarrow->Write();
	fHistoSignMesonNarrow->Write();
	fHistoSBMesonNarrow->Write();

	fHistoYieldMesonWide->Write();
	fHistoYieldMesonPerEventWide->Write();
	fHistoSignMesonWide->Write();
	fHistoSBMesonWide->Write();

	if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0 ){fHistoMassPosition->Write();}
	if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0 ){fHistoFWHMMesonAlpha01->Write();}
	fHistoMassMeson->Write();
	fHistoWidthMeson->Write();
	fHistoFWHMMeson->Write();
	fDeltaPt->Write();
	
	fHistoYieldMesonLeft->Write();
	fHistoYieldMesonLeftPerEvent->Write();
	fHistoSignMesonLeft->Write();
	fHistoSBMesonLeft->Write();

	fHistoYieldMesonLeftNarrow->Write();
	fHistoYieldMesonLeftPerEventNarrow->Write();
	fHistoSignMesonLeftNarrow->Write();
	fHistoSBMesonLeftNarrow->Write();

	fHistoYieldMesonLeftWide->Write();
	fHistoYieldMesonLeftPerEventWide->Write();
	fHistoSignMesonLeftWide->Write();
	fHistoSBMesonLeftWide->Write();

	fHistoMassMesonLeft->Write();
	fHistoWidthMesonLeft->Write();
	fHistoFWHMMesonLeft->Write();
	fMesonFullPtSignal->Write();
	fMesonFullPtBackground->Write();
	fMesonFullPtBackNorm->SetName("Mapping_BackNorm_InvMass_FullPt");
	fMesonFullPtBackNorm->Write();
	fNumberOfGoodESDTracks->Write();
	fEventQuality->Write();

	TString nameHistoSignal;
	TString nameHistoPeakPos;
	TString nameHistoBckNorm;
	TString fitnameSignal;
	TString nameHistoSignalPos;
	if (fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0) {
		for (Int_t ii=fStartPtBin;ii<fNBinsPeakPt;ii++){
			nameHistoPeakPos = Form("InvMassAlpha01_in_Pt_Bin%02d", ii);
			fHistoMappingPeakPosInvMassPtBin[ii]->Write(nameHistoPeakPos.Data());
			fNameFitSignalPos = Form("Signal_InvMassFitPos_in_Pt_Bin%02d", ii);
			if(fFitPeakPosPtBin[ii]!=0x00) fFitPeakPosPtBin[ii]->Write(fNameFitSignalPos.Data());
		}
	}
	for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
		fHistoWeightsBGZbinVsMbin[ii]->Write(Form("BGWeights_%02d", ii));
		fHistoFillPerEventBGZbinVsMbin[ii]->Scale(1./fNEvents);
		fHistoFillPerEventBGZbinVsMbin[ii]->Write(Form("BGPoolsFillstatus_%02d", ii));
	
		fHistoMappingGGInvMassPtBin[ii]->Write();
		nameHistoBckNorm = Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d", ii);
		fHistoMappingBackNormInvMassPtBin[ii]->Write(nameHistoBckNorm.Data());
		nameHistoSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", ii);
		fHistoMappingSignalInvMassPtBin[ii]->Write(nameHistoSignal.Data());
		fitnameSignal = Form("Signal_InvMassFit_in_Pt_Bin%02d", ii);
		if(fFitSignalInvMassPtBin[ii]!=0x00) fFitSignalInvMassPtBin[ii]->Write(fitnameSignal.Data());
	}

	if(optionMC){
		fHistoTrueSignMeson->Write();
		fHistoTrueSBMeson->Write();
		fHistoMCMesonPtWithinAcceptance->Write();
		fHistoMCMesonWithinAccepPt->Write(); // Proper bins in Pt
		fHistoMCMesonPt1->Write(); // Proper bins in Pt
		fHistoYieldTrueMeson->Write();
		fHistoYieldTrueMesonWide->Write();
		fHistoYieldTrueMesonNarrow->Write();
		fHistoYieldTrueMesonReweighted->Write();
		fHistoYieldTrueMesonReweightedWide->Write();
		fHistoYieldTrueMesonReweightedNarrow->Write();
		
		fHistoTrueMesonInvMassVSPt->Write();
		fHistoYieldTrueSecMeson->Write();
		fHistoYieldTrueSecFromK0SMeson->Write();
		fHistoYieldTrueSecMesonWide->Write();
		fHistoYieldTrueSecFromK0SMesonWide->Write();
		fHistoYieldTrueSecMesonNarrow->Write();
		fHistoYieldTrueSecFromK0SMesonNarrow->Write();
		for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
			fHistoMappingTrueMesonInvMassPtBins[ii]->Write();
			fHistoMappingTrueMesonInvMassPtReweightedBins[ii]->Write();
			if (fAdvancedMesonQA){
				if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5){
					if (fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[ii])fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[ii]->Write();
					if (fHistoMappingTrueMesonCaloElectronInvMassPtBins[ii])fHistoMappingTrueMesonCaloElectronInvMassPtBins[ii]->Write();
					if (fHistoMappingTrueMesonCaloPhotonInvMassPtBins[ii])fHistoMappingTrueMesonCaloPhotonInvMassPtBins[ii]->Write();
					if (fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[ii])fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[ii]->Write();
					if (fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[ii])fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[ii]->Write();
					if (fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[ii])fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[ii]->Write();
					if (fFitTrueSignalCaloConvPhotonInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloConvPhotonInvMassPtBin[ii]->Write();
					if (fFitTrueSignalCaloPhotonInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloPhotonInvMassPtBin[ii]->Write();
					if (fFitTrueSignalCaloElectronInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloElectronInvMassPtBin[ii]->Write();
					if (fFitTrueSignalCaloMergedClusterInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloMergedClusterInvMassPtBin[ii]->Write();
					if (fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[ii]->Write();
					if (fFitTrueSignalCaloEMNonLeadingInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloEMNonLeadingInvMassPtBin[ii]->Write();
				} else {
					fHistoMappingTrueGGBckInvMassPtBins[ii]->Write();
					fHistoMappingTrueContBckInvMassPtBins[ii]->Write();
					fHistoMappingTrueAllBckInvMassPtBins[ii]->Write();
				}	
			}
			//fHistoMappingTrueSecMesonInvMassPtBins[ii]->Write();
			//fHistoMappingTrueSecFromK0SMesonInvMassPtBins[ii]->Write();	
			if (fFitTrueSignalInvMassPtBin[ii]!=0x00) fFitTrueSignalInvMassPtBin[ii]->Write();
			if (fFitTrueSignalInvMassPtReweightedBin[ii]!=0x00) fFitTrueSignalInvMassPtReweightedBin[ii]->Write();
		}
	}

	cout << "End writing Uncorrected File" << endl;
	
	fOutput1->Write();
	fOutput1->Close();
}

void SaveCorrectionHistos(TString fCutID, TString fPrefix3)
{
	const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1CorrectionHistos%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix3.Data(),fPeriodFlag.Data(),fCutID.Data());
	fOutput2 = new TFile(nameOutput,"RECREATE");
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	cout<<nameOutput<<endl;
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;

	cout << "Begin writing Correction File" << endl;
	fHistoMCMesonAcceptPt->Write();
	fHistoTrueMesonInvMassVSPt->Write();
	fHistoMonteMesonEffiPt->Write();
	fHistoMonteMesonEffiBackFitPt->Write();
	fHistoMonteMesonNarrowEffiPt->Write();
	fHistoMonteMesonWideEffiPt->Write();
	fHistoMonteMesonLeftEffiPt->Write();
	fHistoMonteMesonLeftNarrowEffiPt->Write();
	fHistoMonteMesonLeftWideEffiPt->Write();
	fHistoMCTrueMesonEffiPt->Write();
	fHistoMCTrueMesonEffiPtReweighted->Write();
	fHistoMCTrueMesonNarrowEffiPt->Write();
	fHistoMCTrueMesonNarrowEffiPtReweighted->Write();
	fHistoMCTrueMesonWideEffiPt->Write();
	fHistoMCTrueMesonWideEffiPtReweighted->Write();
	fHistoMonteMesonEffiFitPt->Write();
	fHistoMonteMesonNarrowEffiFitPt->Write();
	fHistoMonteMesonWideEffiFitPt->Write();
	fHistoMonteMesonLeftEffiFitPt->Write();
	fHistoMonteMesonLeftNarrowEffiFitPt->Write();
	fHistoMonteMesonLeftWideEffiFitPt->Write();
	fHistoMCTrueMesonEffiFitPt->Write();
	fHistoMCTrueMesonNarrowEffiFitPt->Write();
	fHistoMCTrueMesonWideEffiFitPt->Write();
	
	fHistoTrueMassMeson->Write();
	fHistoTrueMassMesonReweighted->Write();
	fHistoTrueFWHMMeson->Write();
	fHistoTrueFWHMMesonReweighted->Write();
	if (fAdvancedMesonQA && (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5)){
		if (fHistoTrueMassMesonCaloPhoton) fHistoTrueMassMesonCaloPhoton->Write();
		if (fHistoTrueMassMesonCaloElectron) fHistoTrueMassMesonCaloElectron->Write();
		if (fHistoTrueMassMesonCaloConvPhoton) fHistoTrueMassMesonCaloConvPhoton->Write();
		if (fHistoTrueMassMesonCaloMergedCluster) fHistoTrueMassMesonCaloMergedCluster->Write();
		if (fHistoTrueMassMesonCaloMergedPartConvCluster) fHistoTrueMassMesonCaloMergedPartConvCluster->Write();
		if (fHistoTrueMassMesonCaloEMNonLeading) fHistoTrueMassMesonCaloEMNonLeading->Write();
		if (fHistoTrueFWHMMesonCaloPhoton) fHistoTrueFWHMMesonCaloPhoton->Write();
		if (fHistoTrueFWHMMesonCaloElectron) fHistoTrueFWHMMesonCaloElectron->Write();
		if (fHistoTrueFWHMMesonCaloConvPhoton) fHistoTrueFWHMMesonCaloConvPhoton->Write();
		if (fHistoTrueFWHMMesonCaloMergedCluster) fHistoTrueFWHMMesonCaloMergedCluster->Write();
		if (fHistoTrueFWHMMesonCaloMergedPartConvCluster) fHistoTrueFWHMMesonCaloMergedPartConvCluster->Write();
		if (fHistoTrueFWHMMesonCaloEMNonLeading) fHistoTrueFWHMMesonCaloEMNonLeading->Write();
	
	}
	fHistoMCMesonPt1->SetName("MC_Meson_genPt");
	fHistoMCMesonPt1->Write(); // Proper bins in Pt
	fHistoMCMesonPt->SetName("MC_Meson_genPt_oldBin");
	fHistoMCMesonPt->Scale(1./fHistoMCMesonPt->GetBinWidth(5));
	//    fHistoMCMesonPt->GetXaxis()->SetRangeUser(0.,25.);
	fHistoMCMesonPt->Write(); // Proper bins in Pt
	if (fHistoMCMesonPtWOWeights){
		fHistoMCMesonPtWOWeights->SetName("MC_Meson_genPt_WOWeights");
		fHistoMCMesonPtWOWeights->Scale(1./fHistoMCMesonPtWOWeights->GetBinWidth(5));
	//       fHistoMCMesonPtWOWeights->GetXaxis()->SetRangeUser(0.,25.);
		fHistoMCMesonPtWOWeights->Write(); // Proper bins in Pt
	}
	if (fHistoMCMesonPtWeights){
		fHistoMCMesonPtWeights->Write("MC_Meson_genPt_Weights"); // Proper bins in Pt
	}   
	fEventQuality->Write();
	fHistoYieldTrueSecMeson->Write();
	fHistoYieldTrueSecFromK0SMeson->Write();
	fHistoYieldTrueSecMesonWide->Write();
	fHistoYieldTrueSecFromK0SMesonWide->Write();
	fHistoYieldTrueSecMesonNarrow->Write();
	fHistoYieldTrueSecFromK0SMesonNarrow->Write();
	fHistoYieldTrueSecFracMeson->Write();
	fHistoYieldTrueSecFracFromK0SMeson->Write();
	fHistoYieldTrueSecFracMesonWide->Write();
	fHistoYieldTrueSecFracFromK0SMesonWide->Write();
	fHistoYieldTrueSecFracMesonNarrow->Write();
	fHistoYieldTrueSecFracFromK0SMesonNarrow->Write();

	fHistoYieldTrueMeson->Write();
	fHistoYieldTrueMesonWide->Write();
	fHistoYieldTrueMesonNarrow->Write();
	fHistoYieldTrueMesonReweighted->Write();
	fHistoYieldTrueMesonReweightedWide->Write();
	fHistoYieldTrueMesonReweightedNarrow->Write();
	
	cout << "end writing Correction File" << endl;
	
	fOutput2->Write();
	fOutput2->Close();
}

void Initialize(TString setPi0, Int_t numberOfBins){
	
	//cout << "MODE in INITIALIZE = " << fMode << endl;
	
	if (setPi0.CompareTo("Eta") == 0){
		fNBinsPt = 		numberOfBins;
		fBinsPt= 			new Double_t[20];
		fNRebin = 		new Int_t[19];

		if (fEnergyFlag.CompareTo("7TeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	5;
			fRow = 		3;
			
			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins for Eta, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt = 12;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] = fBinsEta7TeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta7TeVPtRebin[i];
			}
			fExampleBin = 6;
		} else if (fEnergyFlag.CompareTo("8TeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	5;
			fRow = 		3;
			
			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins for Eta, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt = 12;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] = fBinsEta8TeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVPtRebin[i];
			}
			fExampleBin = 6;
		} else if (fEnergyFlag.CompareTo("2.76TeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	3;
			fRow = 		3;

			if (fNBinsPt > 7) {
				cout << "You have chosen to have more than 7 bins for Eta, this is not possible, it will be reduced to 7" << endl;
				fNBinsPt = 7;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta2760GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPtRebin[i];
			}
			fExampleBin = 4;
		} else if (fEnergyFlag.CompareTo("900GeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	2;
			fRow = 		2;

			if (fNBinsPt > 3) {
				cout << "You have chosen to have more than 3 bins for Eta, this is not possible, it will be reduced to 3" << endl;
				fNBinsPt = 3;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta900GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta900GeVPtRebin[i];
			}
			fExampleBin = 2;
		} else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
	//          fStartPtBin=	1;
	//          fColumn = 	2;
	//          fRow = 		2;
	// 
				fStartPtBin=   2;
				fColumn =   4;
				fRow =      3;
				if (fNBinsPt > 12) {
					cout << "You have chosen to have more than 12 bins, this is not possible, it will be reduced to 12" << endl;
					fNBinsPt = 12;
				}
				for (Int_t i = 0; i < fNBinsPt+1; i++) {
					fBinsPt[i] = fBinsEtaHIPtLHC11h[i];
					if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaHIPtRebinLHC11hFinerBinning[i];
				}         
			fExampleBin = 4;
		} else if( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) { 
			fStartPtBin=   8;
			fColumn =   3;
			fRow =      2;

			if (fNBinsPt > 13) {
				cout << "You have chosen to have more than 13 bins, this is not possible, it will be reduced to 13" << endl;
				fNBinsPt = 13;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEtapPbPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEtapPbPtRebin[i];
			}

			fExampleBin = 9;
		}
		/*fPeakRange = 		new Double_t[2]; 	fPeakRange[0]=	0.48; 	fPeakRange[1]=	0.58; //eta 0.9
		fFitRange = 		new Double_t[2]; 	fFitRange[0]=	0.4; 	fFitRange[1]=	0.65; //eta 0.9
		fBGFitRange = 		new Double_t[2]; 	fBGFitRange[0]=	0.58; 	fBGFitRange[1]=	0.79; //eta 0.9
		fBGFitRangeLeft = 	new Double_t[2]; 	fBGFitRangeLeft[0]=	0.35; 	fBGFitRangeLeft[1]=	0.48;  // eta 09
		fMesonPlotRange = 	new Double_t[2]; 	fMesonPlotRange[0]=	0.53; 	fMesonPlotRange[1]=	0.560;
		fMesonIntRange = 	new Double_t[2]; 	fMesonIntRange[0] = 0.50; 	fMesonIntRange[1]=	0.57;
		fMesonIntRangeWide = new Double_t[2]; 	fMesonIntRangeWide[0]=0.48; 	fMesonIntRangeWide[1]=0.58;
		fMesonIntRangeNarrow = new Double_t[2]; fMesonIntRangeNarrow[0]=0.515; fMesonIntRangeNarrow[1]=0.56;
		fMesonMassRange = 	new Double_t[2]; 	fMesonMassRange[0]=	0.35; 	fMesonMassRange[1]=	0.79;
		fMesonFitRange = 	new Double_t[2]; 	fMesonFitRange[0]=	0.4; 	fMesonFitRange[1]=	0.7;*/
		
		fPeakRange = 		new Double_t[2]; 	fPeakRange[0]=	0.53; 	fPeakRange[1]=	0.57; //eta 0.9
		fFitRange = 		new Double_t[2]; 	fFitRange[0]=	0.4; 	fFitRange[1]=	0.7; //eta 0.9
		fBGFitRange = 		new Double_t[2]; 	fBGFitRange[0]=	0.6; 	fBGFitRange[1]=	0.65; //eta 0.9
		fBGFitRangeLeft = 	new Double_t[2]; 	fBGFitRangeLeft[0]=	0.4; 	fBGFitRangeLeft[1]=	0.48;  // eta 09
		fMesonPlotRange = 	new Double_t[2]; 	fMesonPlotRange[0]=	0.53; 	fMesonPlotRange[1]=	0.560;
		fMesonIntRange = 	new Double_t[2]; 	fMesonIntRange[0] = 0.50; 	fMesonIntRange[1]=	0.57;
		fMesonIntRangeWide = new Double_t[2]; 	fMesonIntRangeWide[0]=0.53; 	fMesonIntRangeWide[1]=0.58;
		fMesonIntRangeNarrow = new Double_t[2]; fMesonIntRangeNarrow[0]=0.515; fMesonIntRangeNarrow[1]=0.56;
		fMesonMassRange = 	new Double_t[2]; 	fMesonMassRange[0]=	0.4; 	fMesonMassRange[1]=	0.7;
		fMesonFitRange = 	new Double_t[2]; 	fMesonFitRange[0]=	0.5; 	fMesonFitRange[1]=	0.6; //0.4, 0.7
		
		//		fMesonFitRangeWithoutPeak = new Double_t[2]; fMesonFitRangeWithoutPeak[0] = 0.4; fMesonFitRangeWithoutPeak[0] = 0.7;
		fMesonId=			221;
		if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
			fMesonWidthExpect = 0.010;
		} else {
			fMesonWidthExpect = 0.01;
		}   
		fMesonLambdaTail = 	0.007;
		fMesonWidthRange = 	new Double_t[2]; 	 fMesonWidthRange[0]=0.005;     fMesonWidthRange[1]=0.02;
		fMesonLambdaTailRange = new Double_t[2]; fMesonLambdaTailRange[0]=0.004; fMesonLambdaTailRange[1]=0.01;
		fMidPt = 							new Double_t[2]; fMidPt[0] = 1.5; fMidPt[1] = 2.5;
		
		
		if ((fMode == 2) || (fMode == 3)){
			fPeakRange[0]=	0.52; 	fPeakRange[1]=	0.57; //eta 0.9
			fFitRange[0]=	0.4; 	fFitRange[1]=	0.7; //eta 0.9
			fBGFitRange[0]=	0.6; 	fBGFitRange[1]=	0.65; //eta 0.9
			fBGFitRangeLeft[0]=	0.4; 	fBGFitRangeLeft[1]=	0.48;  // eta 09
			fMesonFitRange[0]=	0.48; 	fMesonFitRange[1]=	0.62;
			fMesonLambdaTail = 	0.007;
			fMesonWidthRange[0]=0.005; 	fMesonWidthRange[1]=0.02;
			fMesonWidthExpect = 0.012;
			fMesonLambdaTailRange[0]=0.004; fMesonLambdaTailRange[1]=0.01;
		}
		
		if ((fMode == 4) || (fMode == 5)){
			fPeakRange[0]=	0.52; 	fPeakRange[1]=	0.57; //eta 0.9
			fFitRange[0]=	0.4; 	fFitRange[1]=	0.7; //eta 0.9
			fBGFitRange[0]=	0.62; 	fBGFitRange[1]=	0.67; //eta 0.9
			fBGFitRangeLeft[0]=	0.4; 	fBGFitRangeLeft[1]=	0.48;  // eta 09
			fMesonFitRange[0]=	0.48; 	fMesonFitRange[1]=	0.62;
			fMesonLambdaTail = 	0.007;
			fMesonWidthRange[0]=0.005; 	fMesonWidthRange[1]=0.02;
			fMesonWidthExpect = 0.012;
			fMesonLambdaTailRange[0]=0.004; fMesonLambdaTailRange[1]=0.01;
		}

		
		
	} else if (setPi0.CompareTo("Omega") == 0){
		fNBinsPt = 		numberOfBins;
		fBinsPt= 			new Double_t[20];
		fNRebin = 		new Int_t[19];

		if (fEnergyFlag.CompareTo("7TeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	5;
			fRow = 		3;
			
			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins for Eta, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt = 12;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] = fBinsEta7TeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta7TeVPtRebin[i];
			}
			fExampleBin = 6;
		} else if (fEnergyFlag.CompareTo("8TeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	5;
			fRow = 		3;
			
			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins for Eta, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt = 12;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] = fBinsEta8TeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVPtRebin[i];
			}
			fExampleBin = 6;
		} else if (fEnergyFlag.CompareTo("2.76TeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	3;
			fRow = 		3;

			if (fNBinsPt > 7) {
				cout << "You have chosen to have more than 7 bins for Eta, this is not possible, it will be reduced to 7" << endl;
				fNBinsPt = 7;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta2760GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPtRebin[i];
			}
			fExampleBin = 4;
		} else if (fEnergyFlag.CompareTo("900GeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	2;
			fRow = 		2;

			if (fNBinsPt > 3) {
				cout << "You have chosen to have more than 3 bins for Eta, this is not possible, it will be reduced to 3" << endl;
				fNBinsPt = 3;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta900GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta900GeVPtRebin[i];
			}
			fExampleBin = 2;
		} else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
	//          fStartPtBin=	1;
	//          fColumn = 	2;
	//          fRow = 		2;
	// 
				fStartPtBin=   2;
				fColumn =   4;
				fRow =      3;
				if (fNBinsPt > 12) {
					cout << "You have chosen to have more than 12 bins, this is not possible, it will be reduced to 12" << endl;
					fNBinsPt = 12;
				}
				for (Int_t i = 0; i < fNBinsPt+1; i++) {
					fBinsPt[i] = fBinsEtaHIPtLHC11h[i];
					if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaHIPtRebinLHC11hFinerBinning[i];
				}         
			fExampleBin = 4;
		} else if( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) { 
			fStartPtBin=   1;
			fColumn =   4;
			fRow =      2;

			if (fNBinsPt > 7) {
				cout << "You have chosen to have more than 7 bins, this is not possible, it will be reduced to 7" << endl;
				fNBinsPt = 7;
			}
			
			if (fMode == 4)
			{
				fStartPtBin=   3;
				fColumn =   3;
				fRow =      2;
				//fNBinsPt = TMath::Min(fNBinsPt, 7);
			}
			
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsOmegapPbPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsOmegapPbPtRebin[i];
			}

			fExampleBin = 3;
		}
		/*fPeakRange = 		new Double_t[2]; 	fPeakRange[0]=	0.68; 	fPeakRange[1]=	0.81; //eta 0.9
		fFitRange = 		new Double_t[2]; 	fFitRange[0]=	0.55; 	fFitRange[1]=	0.83; //eta 0.9
		fBGFitRange = 		new Double_t[2]; 	fBGFitRange[0]=	0.81; 	fBGFitRange[1]=	0.89; //eta 0.9
		fBGFitRangeLeft = 	new Double_t[2]; 	fBGFitRangeLeft[0]=	0.4; 	fBGFitRangeLeft[1]=	0.68;  // eta 09
		fMesonPlotRange = 	new Double_t[2]; 	fMesonPlotRange[0]=	0.75; 	fMesonPlotRange[1]=	0.79;
		fMesonIntRange = 	new Double_t[2]; 	fMesonIntRange[0] = 0.7; 	fMesonIntRange[1]=	0.8;
		fMesonIntRangeWide = new Double_t[2]; 	fMesonIntRangeWide[0]=0.68; 	fMesonIntRangeWide[1]=0.81;
		fMesonIntRangeNarrow = new Double_t[2]; fMesonIntRangeNarrow[0]=0.715; fMesonIntRangeNarrow[1]=0.8;
		fMesonMassRange = 	new Double_t[2]; 	fMesonMassRange[0]=	0.4; 	fMesonMassRange[1]=	0.89;
		fMesonFitRange = 	new Double_t[2]; 	fMesonFitRange[0]=	0.55; 	fMesonFitRange[1]=	0.85;
		*/
		
		fPeakRange = 		new Double_t[2]; 	fPeakRange[0]=	0.75; 	fPeakRange[1]=	0.81; //eta 0.9
		fFitRange = 		new Double_t[2]; 	fFitRange[0]=	0.55; 	fFitRange[1]=	0.89; //eta 0.9
		fBGFitRange = 		new Double_t[2]; 	fBGFitRange[0]=	0.83; 	fBGFitRange[1]=	0.89; //eta 0.9
		fBGFitRangeLeft = 	new Double_t[2]; 	fBGFitRangeLeft[0]=	0.6; 	fBGFitRangeLeft[1]=	0.69;  // eta 09
		fMesonPlotRange = 	new Double_t[2]; 	fMesonPlotRange[0]=	0.75; 	fMesonPlotRange[1]=	0.79;
		fMesonIntRange = 	new Double_t[2]; 	fMesonIntRange[0] = 0.7; 	fMesonIntRange[1]=	0.8;
		fMesonIntRangeWide = new Double_t[2]; 	fMesonIntRangeWide[0]=0.7; 	fMesonIntRangeWide[1]=0.81;
		fMesonIntRangeNarrow = new Double_t[2]; fMesonIntRangeNarrow[0]=0.715; fMesonIntRangeNarrow[1]=0.8;
		fMesonMassRange = 	new Double_t[2]; 	fMesonMassRange[0]=	0.55; 	fMesonMassRange[1]=	0.89;
		fMesonFitRange = 	new Double_t[2]; 	fMesonFitRange[0]=	0.65; 	fMesonFitRange[1]=	0.85; //0.55, 0.9
		
		//		fMesonFitRangeWithoutPeak = new Double_t[2]; fMesonFitRangeWithoutPeak[0] = 0.4; fMesonFitRangeWithoutPeak[0] = 0.7;
		fMesonId=			223;
		if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
			fMesonWidthExpect = 0.015;
		} else {
			fMesonWidthExpect = 0.01;
		}   
		
		fMesonLambdaTail = 	0.007;
		fMesonWidthRange = 	new Double_t[2]; 	 fMesonWidthRange[0]=0.005;     fMesonWidthRange[1]=0.02;
		fMesonLambdaTailRange = new Double_t[2]; fMesonLambdaTailRange[0]=0.004; fMesonLambdaTailRange[1]=0.01;
		fMidPt = 							new Double_t[2]; fMidPt[0] = 1.5; fMidPt[1] = 2.5;
		if ((fMode == 2) || (fMode == 3)){
			fPeakRange[0]=	0.75; 	fPeakRange[1]=	0.81; //eta 0.9
			fFitRange[0]=	0.55; 	fFitRange[1]=	0.89; //eta 0.9
			fBGFitRange[0]=	0.83; 	fBGFitRange[1]=	0.89; //eta 0.9
			fBGFitRangeLeft[0]=	0.6; 	fBGFitRangeLeft[1]=	0.67;  // eta 09
			fMesonFitRange[0]=	0.6; 	fMesonFitRange[1]=	0.85;
			fMesonLambdaTail = 	0.007;
			fMesonWidthRange[0]=0.005; 	fMesonWidthRange[1]=0.02;
			fMesonWidthExpect = 0.01;
			fMesonLambdaTailRange[0]=0.004; fMesonLambdaTailRange[1]=0.01;
		}
		
		if ((fMode == 4) || (fMode == 5)){
			fPeakRange[0]=	0.75; 	fPeakRange[1]=	0.81; //eta 0.9
			fFitRange[0]=	0.55; 	fFitRange[1]=	0.89; //eta 0.9
			fBGFitRange[0]=	0.83; 	fBGFitRange[1]=	0.89; //eta 0.9
			fBGFitRangeLeft[0]=	0.6; 	fBGFitRangeLeft[1]=	0.66;  // eta 09
			fMesonFitRange[0]=	0.6; 	fMesonFitRange[1]=	0.85;
			fMesonLambdaTail = 	0.007;
			fMesonWidthRange[0]=0.005; 	fMesonWidthRange[1]=0.02;
			fMesonWidthExpect = 0.01;
			fMesonLambdaTailRange[0]=0.004; fMesonLambdaTailRange[1]=0.01;
		}

	} 
	
	fMesonLambdaTailMC = 	fMesonLambdaTail;
	fMesonWidthExpectMC = 	fMesonWidthExpect;
	fMesonWidthRangeMC = 	new Double_t[2]; 	fMesonWidthRangeMC[0]=fMesonWidthRange[0]; 	fMesonWidthRangeMC[1]=fMesonWidthRange[1];
	fFullPt = 				new Double_t[2]; 	fFullPt[0] = 0.4; 							fFullPt[1] = 15.;

	fMesonCurIntRange = 									new Double_t[2];
	fMesonCurIntRangeBackFit = 								new Double_t[2];
	fMesonCurIntRangeWide = 								new Double_t[2];
	fMesonCurIntRangeNarrow = 								new Double_t[2];
	fMesonCurLeftIntRange = 								new Double_t[2];
	fMesonCurLeftIntRangeWide = 							new Double_t[2];
	fMesonCurLeftIntRangeNarrow = 							new Double_t[2];
	fMesonTrueIntRange = 									new Double_t[2];
	fMesonTrueIntRangeWide = 								new Double_t[2];
	fMesonTrueIntRangeNarrow = 								new Double_t[2];
	fMesonTrueIntReweightedRange =							new Double_t[2];
	fMesonTrueIntReweightedRangeWide =						new Double_t[2];
	fMesonTrueIntReweightedRangeNarrow =					new Double_t[2];

   
	fGGYields = 											new Double_t[fNBinsPt];
	fMesonYieldsBackFit =									new Double_t[fNBinsPt];
	fBckYields = 											new Double_t[fNBinsPt];
	fTotalBckYields = 										new Double_t[fNBinsPt];
	fMesonYields = 											new Double_t[fNBinsPt];
	fMesonTrueYields = 										new Double_t[fNBinsPt];
	fMesonTrueYieldsReweighted =							new Double_t[fNBinsPt];
	fMesonTrueSecYields = 									new Double_t[fNBinsPt];
	fMesonTrueSecFromK0SYields = 							new Double_t[fNBinsPt];
	fMesonYieldsFunc = 										new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFunc = 							new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncBackFit = 					new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFunc = 						new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncBackFit = 				new Double_t[fNBinsPt];
	fMesonYieldsPerEvent = 									new Double_t[fNBinsPt];
	fMesonYieldsPerEventBackFit = 							new Double_t[fNBinsPt];
	fMesonMass = 											new Double_t[fNBinsPt];
	fMesonMassBackFit = 									new Double_t[fNBinsPt];
	fMesonWidth = 											new Double_t[fNBinsPt];
	fMesonWidthBackFit = 									new Double_t[fNBinsPt];
	fMesonSB = 												new Double_t[fNBinsPt];
	fMesonSign = 											new Double_t[fNBinsPt];
	fMesonFWHM = 											new Double_t[fNBinsPt];
	fMesonFWHMAlpha01 =										new Double_t[fNBinsPt];
	fMesonTrueMass = 										new Double_t[fNBinsPt];
	fMesonTrueMassReweighted = 								new Double_t[fNBinsPt];
	fMesonTrueFWHM = 										new Double_t[fNBinsPt];
	fMesonTrueFWHMReweighted =								new Double_t[fNBinsPt];
	fMesonTrueSB = 											new Double_t[fNBinsPt];
	fMesonTrueSign = 										new Double_t[fNBinsPt];
	
	fMesonTrueMassCaloPhoton = 								new Double_t[fNBinsPt];
	fMesonTrueMassCaloElectron = 							new Double_t[fNBinsPt];
	fMesonTrueMassCaloConvPhoton = 							new Double_t[fNBinsPt];
	fMesonTrueMassCaloMergedCluster = 						new Double_t[fNBinsPt];
	fMesonTrueMassCaloMergedClusterPartConv = 				new Double_t[fNBinsPt];
	fMesonTrueMassCaloEMNonLeading = 						new Double_t[fNBinsPt];
	fMesonTrueFWHMCaloPhoton = 								new Double_t[fNBinsPt];
	fMesonTrueFWHMCaloElectron = 							new Double_t[fNBinsPt];
	fMesonTrueFWHMCaloConvPhoton = 							new Double_t[fNBinsPt];
	fMesonTrueFWHMCaloMergedCluster = 						new Double_t[fNBinsPt];
	fMesonTrueFWHMCaloMergedClusterPartConv = 				new Double_t[fNBinsPt];
	fMesonTrueFWHMCaloEMNonLeading = 						new Double_t[fNBinsPt];
	
	// Normalization at the left of the peak
	fGGYieldsLeft = 										new Double_t[fNBinsPt];
	fBckYieldsLeft = 										new Double_t[fNBinsPt];
	fTotalBckYieldsLeft = 									new Double_t[fNBinsPt];
	fMesonYieldsLeft = 										new Double_t[fNBinsPt];
	fMesonYieldsFuncLeft = 									new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeft = 						new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeft = 					new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEvent = 								new Double_t[fNBinsPt];
	fMesonMassLeft = 										new Double_t[fNBinsPt];
	fMesonWidthLeft = 										new Double_t[fNBinsPt];
	fMesonSBLeft = 											new Double_t[fNBinsPt];
	fMesonSignLeft = 										new Double_t[fNBinsPt];
	fMesonFWHMLeft = 										new Double_t[fNBinsPt];

	// Narrow Integration Window
	fGGYieldsNarrow = 										new Double_t[fNBinsPt];
	fBckYieldsNarrow = 										new Double_t[fNBinsPt];
	fTotalBckYieldsNarrow =	 								new Double_t[fNBinsPt];
	fMesonYieldsNarrow = 									new Double_t[fNBinsPt];
	fMesonTrueYieldsNarrow = 								new Double_t[fNBinsPt];
	fMesonTrueYieldsReweightedNarrow =						new Double_t[fNBinsPt];
	fMesonTrueSecYieldsNarrow = 							new Double_t[fNBinsPt];
	fMesonTrueSecFromK0SYieldsNarrow = 						new Double_t[fNBinsPt];
	fMesonYieldsFuncNarrow = 								new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncNarrow = 					new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncNarrow = 					new Double_t[fNBinsPt];
	fMesonYieldsPerEventNarrow = 							new Double_t[fNBinsPt];
	fMesonSBNarrow = 										new Double_t[fNBinsPt];
	fMesonSignNarrow = 										new Double_t[fNBinsPt];

	fGGYieldsLeftNarrow = 									new Double_t[fNBinsPt];
	fBckYieldsLeftNarrow = 									new Double_t[fNBinsPt];
	fTotalBckYieldsLeftNarrow = 							new Double_t[fNBinsPt];
	fMesonYieldsLeftNarrow = 								new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftNarrow = 							new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftNarrow = 				new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftNarrow = 				new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventNarrow = 						new Double_t[fNBinsPt];
	fMesonSBLeftNarrow = 									new Double_t[fNBinsPt];
	fMesonSignLeftNarrow = 									new Double_t[fNBinsPt];

	// Wide Integration Window
	fGGYieldsWide = 										new Double_t[fNBinsPt];
	fBckYieldsWide = 										new Double_t[fNBinsPt];
	fTotalBckYieldsWide =	 								new Double_t[fNBinsPt];
	fMesonYieldsWide = 										new Double_t[fNBinsPt];
	fMesonTrueYieldsWide = 									new Double_t[fNBinsPt];
	fMesonTrueYieldsReweightedWide =						new Double_t[fNBinsPt];
	fMesonTrueSecYieldsWide = 								new Double_t[fNBinsPt];
	fMesonTrueSecFromK0SYieldsWide = 						new Double_t[fNBinsPt];
	fMesonYieldsFuncWide = 									new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncWide = 						new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncWide =					new Double_t[fNBinsPt];
	fMesonYieldsPerEventWide = 								new Double_t[fNBinsPt];
	fMesonSBWide = 											new Double_t[fNBinsPt];
	fMesonSignWide = 										new Double_t[fNBinsPt];

	fGGYieldsLeftWide = 									new Double_t[fNBinsPt];
	fBckYieldsLeftWide = 									new Double_t[fNBinsPt];
	fTotalBckYieldsLeftWide = 								new Double_t[fNBinsPt];
	fMesonYieldsLeftWide = 									new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftWide = 								new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftWide =					new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftWide = 				new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventWide = 							new Double_t[fNBinsPt];
	fMesonSBLeftWide = 										new Double_t[fNBinsPt];
	fMesonSignLeftWide = 									new Double_t[fNBinsPt];

	fGGYieldsError = 										new Double_t[fNBinsPt];
	fMesonYieldsBackFitError =								new Double_t[fNBinsPt];
	fBckYieldsError = 										new Double_t[fNBinsPt];
	fTotalBckYieldsError = 									new Double_t[fNBinsPt];
	fMesonYieldsError = 									new Double_t[fNBinsPt];
	fMesonTrueYieldsError = 								new Double_t[fNBinsPt];
	fMesonTrueYieldsReweightedError =						new Double_t[fNBinsPt];
	fMesonTrueSecYieldsError = 								new Double_t[fNBinsPt];
	fMesonTrueSecFromK0SYieldsError = 						new Double_t[fNBinsPt];
	fMesonYieldsFuncError = 								new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncError = 						new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncBackFitError = 				new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncError =					new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncBackFitError =			new Double_t[fNBinsPt];
	fMesonYieldsPerEventError = 							new Double_t[fNBinsPt];
	fMesonYieldsPerEventBackFitError = 						new Double_t[fNBinsPt];
	fMesonMassError = 										new Double_t[fNBinsPt];
	fMesonMassBackFitError = 								new Double_t[fNBinsPt];
	fMesonWidthError = 										new Double_t[fNBinsPt];
	fMesonWidthBackFitError = 								new Double_t[fNBinsPt];
	fMesonSBError = 										new Double_t[fNBinsPt];
	fMesonSignError = 										new Double_t[fNBinsPt];
	fMesonFWHMError = 										new Double_t[fNBinsPt];
	fMesonFWHMAlpha01Error = 								new Double_t[fNBinsPt];
	fMesonTrueMassError = 									new Double_t[fNBinsPt];
	fMesonTrueMassReweightedError =		 					new Double_t[fNBinsPt];
	fMesonTrueFWHMError = 									new Double_t[fNBinsPt];
	fMesonTrueFWHMReweightedError = 						new Double_t[fNBinsPt];
	fMesonTrueSBError = 									new Double_t[fNBinsPt];
	fMesonTrueSignError = 									new Double_t[fNBinsPt];
	fMesonTrueMassErrorCaloPhoton = 						new Double_t[fNBinsPt];
	fMesonTrueMassErrorCaloElectron = 						new Double_t[fNBinsPt];
	fMesonTrueMassErrorCaloConvPhoton = 					new Double_t[fNBinsPt];
	fMesonTrueMassErrorCaloMergedCluster = 					new Double_t[fNBinsPt];
	fMesonTrueMassErrorCaloMergedClusterPartConv = 			new Double_t[fNBinsPt];
	fMesonTrueMassErrorCaloEMNonLeading = 					new Double_t[fNBinsPt];
	fMesonTrueFWHMErrorCaloPhoton = 						new Double_t[fNBinsPt];
	fMesonTrueFWHMErrorCaloElectron = 						new Double_t[fNBinsPt];
	fMesonTrueFWHMErrorCaloConvPhoton = 					new Double_t[fNBinsPt];
	fMesonTrueFWHMErrorCaloMergedCluster = 					new Double_t[fNBinsPt];
	fMesonTrueFWHMErrorCaloMergedClusterPartConv = 			new Double_t[fNBinsPt];
	fMesonTrueFWHMErrorCaloEMNonLeading = 					new Double_t[fNBinsPt];
	
	fGGYieldsLeftError = 									new Double_t[fNBinsPt];
	fBckYieldsLeftError =	 								new Double_t[fNBinsPt];
	fTotalBckYieldsLeftError = 								new Double_t[fNBinsPt];
	fMesonYieldsLeftError = 								new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftError = 							new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftError = 					new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftError = 				new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventError = 						new Double_t[fNBinsPt];
	fMesonMassLeftError = 									new Double_t[fNBinsPt];
	fMesonWidthLeftError = 									new Double_t[fNBinsPt];
	fMesonSBLeftError = 									new Double_t[fNBinsPt];
	fMesonSignLeftError = 									new Double_t[fNBinsPt];
	fMesonFWHMLeftError = 									new Double_t[fNBinsPt];

	// Narrow integration Window
	fGGYieldsNarrowError = 									new Double_t[fNBinsPt];
	fBckYieldsNarrowError = 								new Double_t[fNBinsPt];
	fTotalBckYieldsNarrowError = 							new Double_t[fNBinsPt];
	fMesonYieldsNarrowError = 								new Double_t[fNBinsPt];
	fMesonTrueYieldsNarrowError = 							new Double_t[fNBinsPt];
	fMesonTrueYieldsReweightedNarrowError =					new Double_t[fNBinsPt];
	fMesonTrueSecYieldsNarrowError = 						new Double_t[fNBinsPt];
	fMesonTrueSecFromK0SYieldsNarrowError = 				new Double_t[fNBinsPt];
	fMesonYieldsFuncNarrowError = 							new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncNarrowError =				new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncNarrowError = 			new Double_t[fNBinsPt];
	fMesonYieldsPerEventNarrowError = 						new Double_t[fNBinsPt];
	fMesonSBNarrowError = 									new Double_t[fNBinsPt];
	fMesonSignNarrowError = 								new Double_t[fNBinsPt];

	fGGYieldsLeftNarrowError = 								new Double_t[fNBinsPt];
	fBckYieldsLeftNarrowError = 							new Double_t[fNBinsPt];
	fTotalBckYieldsLeftNarrowError = 						new Double_t[fNBinsPt];
	fMesonYieldsLeftNarrowError = 							new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftNarrowError = 						new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftNarrowError =			new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftNarrowError = 		new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventNarrowError = 					new Double_t[fNBinsPt];
	fMesonSBLeftNarrowError = 								new Double_t[fNBinsPt];
	fMesonSignLeftNarrowError = 							new Double_t[fNBinsPt];

	// Wide integration Window
	fGGYieldsWideError = 									new Double_t[fNBinsPt];
	fBckYieldsWideError = 									new Double_t[fNBinsPt];
	fTotalBckYieldsWideError = 								new Double_t[fNBinsPt];
	fMesonYieldsWideError = 								new Double_t[fNBinsPt];
	fMesonTrueYieldsWideError = 							new Double_t[fNBinsPt];
	fMesonTrueYieldsReweightedWideError =					new Double_t[fNBinsPt];
	fMesonTrueSecYieldsWideError = 							new Double_t[fNBinsPt];
	fMesonTrueSecFromK0SYieldsWideError = 					new Double_t[fNBinsPt];	
	fMesonYieldsFuncWideError = 							new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncWideError = 					new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncWideError = 				new Double_t[fNBinsPt];
	fMesonYieldsPerEventWideError = 						new Double_t[fNBinsPt];
	fMesonSBWideError = 									new Double_t[fNBinsPt];
	fMesonSignWideError = 									new Double_t[fNBinsPt];

	fGGYieldsLeftWideError = 								new Double_t[fNBinsPt];
	fBckYieldsLeftWideError = 								new Double_t[fNBinsPt];
	fTotalBckYieldsLeftWideError = 							new Double_t[fNBinsPt];
	fMesonYieldsLeftWideError = 							new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftWideError = 						new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftWideError = 				new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftWideError = 			new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventWideError = 					new Double_t[fNBinsPt];
	fMesonSBLeftWideError = 								new Double_t[fNBinsPt];
	fMesonSignLeftWideError = 								new Double_t[fNBinsPt];

	fHistoMappingTrueMesonInvMassPtBins = 					new TH1D*[fNBinsPt];    
	fHistoMappingTrueMesonInvMassPtReweightedBins =  		new TH1D*[fNBinsPt];    
	fHistoMappingTrueGGBckInvMassPtBins = 					new TH1D*[fNBinsPt];    
	fHistoMappingTrueContBckInvMassPtBins = 				new TH1D*[fNBinsPt];    
	fHistoMappingTrueAllBckInvMassPtBins = 					new TH1D*[fNBinsPt];    
	fHistoMappingTrueSecMesonInvMassPtBins = 				new TH1D*[fNBinsPt];    
	fHistoMappingTrueSecFromK0SMesonInvMassPtBins = 		new TH1D*[fNBinsPt];    
	fHistoMappingTrueMesonCaloPhotonInvMassPtBins =			new TH1D*[fNBinsPt];    
	fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins =		new TH1D*[fNBinsPt];    
	fHistoMappingTrueMesonCaloElectronInvMassPtBins =		new TH1D*[fNBinsPt];    
	fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins =	new TH1D*[fNBinsPt];    
	fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins =	new TH1D*[fNBinsPt];    
	fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins =	new TH1D*[fNBinsPt];    

	fHistoWeightsBGZbinVsMbin = 							new TH2F*[fNBinsPt];    
	fHistoFillPerEventBGZbinVsMbin = 						new TH2F*[fNBinsPt];    
	
	fHistoMappingGGInvMassPtBin = 							new TH1D*[fNBinsPt];    
	fHistoMappingGGInvMassBackFitPtBin =					new TH1D*[fNBinsPt];
	fHistoMappingGGInvMassBackFitWithoutSignalPtBin =		new TH1D*[fNBinsPt];
	fHistoMappingBackInvMassPtBin = 						new TH1D*[fNBinsPt];
	fHistoMappingBackNormInvMassPtBin = 					new TH1D*[fNBinsPt];
	fHistoMappingSignalInvMassPtBin =						new TH1D*[fNBinsPt];
	fHistoMappingRatioSBInvMassPtBin= 						new TH1D*[fNBinsPt];

	fBackgroundFitPol = 									new TF1*[fNBinsPt];
	fFitSignalInvMassPtBin = 								new TF1*[fNBinsPt];
	fFitSignalInvMassBackFitPtBin = 						new TF1*[fNBinsPt];
	fFitTrueSignalInvMassPtBin = 							new TF1*[fNBinsPt];
	fFitTrueSignalInvMassPtReweightedBin =					new TF1*[fNBinsPt];
	fFitTrueSignalCaloConvPhotonInvMassPtBin =				new TF1*[fNBinsPt];
	fFitTrueSignalCaloElectronInvMassPtBin =				new TF1*[fNBinsPt];
	fFitTrueSignalCaloEMNonLeadingInvMassPtBin =			new TF1*[fNBinsPt];
	fFitTrueSignalCaloMergedClusterInvMassPtBin = 			new TF1*[fNBinsPt];
	fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin =	new TF1*[fNBinsPt];
	fFitTrueSignalCaloPhotonInvMassPtBin = 					new TF1*[fNBinsPt];
	
	fFitSignalPeakPosInvMassPtBin = 						new TF1*[fNBinsPt];
	fFitSignalPeakPosInvMassBackFitPtBin = 					new TF1*[fNBinsPt];
	fFitBckInvMassPtBin = 									new TF1*[fNBinsPt];
	fFitBckInvMassBackFitPtBin =							new TF1*[fNBinsPt];
	fFitRatioInvMassPtBin = 								new TF1*[fNBinsPt];
	// Histograms for normalization on the left of the peak
	fHistoMappingBackNormInvMassLeftPtBin = 				new TH1D*[fNBinsPt];
	fHistoMappingSignalInvMassLeftPtBin = 					new TH1D*[fNBinsPt];
	fHistoMappingPeakPosInvMassPtBin = 						new TH1D*[fNBinsPt];

	fFitInvMassLeftPtBin = 									new TF1*[fNBinsPt];
	fFitSignalPeakPosInvMassLeftPtBin = 					new TF1*[fNBinsPt];
	fFitBckInvMassLeftPtBin = 								new TF1*[fNBinsPt];
	fFitWithPol2ForBG = 									new TF1*[fNBinsPt];
	fFitPeakPosPtBin = 										new TF1*[fNBinsPt];
	fMesonMassPeakPos =										new Double_t[fNBinsPt];
	fMesonMassPeakPosError = 								new Double_t[fNBinsPt];

	for(Int_t i = 0;i<fNBinsPt; i++){
		fHistoMappingTrueMesonInvMassPtBins[i] = 					NULL;
		fHistoMappingTrueMesonInvMassPtReweightedBins[i] =			NULL;// array of histos for pt slices
		fHistoMappingTrueGGBckInvMassPtBins[i] =	 				NULL;    // array of histos for pt slices
		fHistoMappingTrueContBckInvMassPtBins[i] = 					NULL;    // array of histos for pt slices
		fHistoMappingTrueAllBckInvMassPtBins[i] = 					NULL;    // array of histos for pt slices
		fHistoMappingTrueSecMesonInvMassPtBins[i] =	 				NULL;
		fHistoMappingTrueSecFromK0SMesonInvMassPtBins[i] = 			NULL;
		fHistoMappingTrueMesonCaloPhotonInvMassPtBins[i] =			NULL;    
		fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[i] =		NULL;    
		fHistoMappingTrueMesonCaloElectronInvMassPtBins[i] =		NULL;    
		fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[i] =	NULL;    
		fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[i] =	NULL;    
		fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins[i] =	NULL;    
		
		fHistoMappingGGInvMassPtBin[i] = 							NULL;
		fHistoMappingGGInvMassBackFitPtBin[i] = 					NULL;
		fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i] = 		NULL;
		fBackgroundFitPol[i] = 										NULL;
		fHistoMappingBackInvMassPtBin[i] = 							NULL;
		fHistoMappingBackNormInvMassPtBin[i] = 						NULL;
		fHistoMappingSignalInvMassPtBin[i] = 						NULL;
		fHistoMappingRatioSBInvMassPtBin[i] = 						NULL;

		fFitSignalInvMassPtBin[i] = 								NULL;
		fFitSignalInvMassBackFitPtBin[i] = 							NULL;
		fFitTrueSignalInvMassPtBin[i] = 							NULL;
		fFitTrueSignalInvMassPtReweightedBin[i] = 					NULL;
		fFitTrueSignalCaloConvPhotonInvMassPtBin[i] = 				NULL;
		fFitTrueSignalCaloElectronInvMassPtBin[i] = 				NULL;
		fFitTrueSignalCaloEMNonLeadingInvMassPtBin[i] = 			NULL;
		fFitTrueSignalCaloMergedClusterInvMassPtBin[i] = 			NULL;
		fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[i] = 	NULL;
		fFitTrueSignalCaloPhotonInvMassPtBin[i] = 					NULL;

		fFitSignalPeakPosInvMassPtBin[i] = 							NULL;
		fFitSignalPeakPosInvMassBackFitPtBin[i] = 					NULL;
		fFitBckInvMassPtBin[i] = 									NULL;
		fFitBckInvMassBackFitPtBin[i] = 							NULL;
		fFitRatioInvMassPtBin[i] = 									NULL;
		// Histograms for normalization on the left of the peak
		fHistoMappingBackNormInvMassLeftPtBin[i] = 					NULL;
		fHistoMappingSignalInvMassLeftPtBin[i] = 					NULL;
		fHistoMappingPeakPosInvMassPtBin[i] = 						NULL;
		
		fFitInvMassLeftPtBin[i] = 									NULL;
		fFitSignalPeakPosInvMassLeftPtBin[i] = 						NULL;
		fFitBckInvMassLeftPtBin[i] = 								NULL;
		fFitWithPol2ForBG[i] = 										NULL;
		fFitPeakPosPtBin[i] = 										NULL;
	}
}

void CalculateFWHM(TF1 * fFunc)
{
   // Default function
	TF1* fFunc_def;
	fFunc_def = new TF1("fFunc_def","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);
	fFunc_def->SetParameter(0,fFunc->GetParameter(0));
	fFunc_def->SetParameter(1,fFunc->GetParameter(1));
	fFunc_def->SetParameter(2,fFunc->GetParameter(2));
	fFunc_def->SetParameter(3,fFunc->GetParameter(3));



	//FWHM
	fFWHMFunc = fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,fFunc_def->GetParameter(1), fMesonFitRange[1]) - fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_def->GetParameter(1));

	//FWHM error +
	TF1* fFunc_plus;
	//	fFunc_plus = fFunc;
	fFunc_plus = new TF1("fFunc_plus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

	fFunc_plus->SetParameter(0,fFunc->GetParameter(0) + fFunc->GetParError(0));
	fFunc_plus->SetParameter(1,fFunc->GetParameter(1) + fFunc->GetParError(1));
	fFunc_plus->SetParameter(2,fFunc->GetParameter(2) + fFunc->GetParError(2));
	fFunc_plus->SetParameter(3,fFunc->GetParameter(3) + fFunc->GetParError(3));
	Double_t FWHM_plus = fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fFunc_plus->GetParameter(1), fMesonFitRange[1]) - fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_plus->GetParameter(1));

	//FWHM error -
	TF1* fFunc_minus;
	//	fFunc_minus = fFunc;
	fFunc_minus = new TF1("fFunc_minus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);
	fFunc_minus->SetParameter(0,fFunc->GetParameter(0) - fFunc->GetParError(0));
	fFunc_minus->SetParameter(1,fFunc->GetParameter(1) - fFunc->GetParError(1));
	fFunc_minus->SetParameter(2,fFunc->GetParameter(2) - fFunc->GetParError(2));
	fFunc_minus->SetParameter(3,fFunc->GetParameter(3) - fFunc->GetParError(3));
	Double_t FWHM_minus =  fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fFunc_minus->GetParameter(1), fMesonFitRange[1]) -fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_minus->GetParameter(1));

	Double_t Error1 = TMath::Abs(fFWHMFunc-FWHM_plus);
	Double_t Error2 = TMath::Abs(fFWHMFunc-FWHM_minus);

	if(Error1>=Error2) fFWHMFuncError = Error1;
	if(Error1<Error2) fFWHMFuncError = Error2;
	}

	//Crystal ball function for signal +linear background, parameters are 0:normalization,1:mean,2:sigma,3:n,4:alpha;
	Double_t CrystalBallBck(Double_t *x,Double_t *par) {

	Double_t t = (x[0]-par[1])/par[2];
	if (par[4] < 0) t = -t;

	Double_t absAlpha = fabs((Double_t)par[4]);

	if (t >= -absAlpha) {
		return par[0]*exp(-0.5*t*t)+par[5]+par[6]*x[0];
	}
	else {
		Double_t a =  TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
		Double_t b= par[3]/absAlpha - absAlpha;

		return par[0]*(a/TMath::Power(b - t, par[3]))+par[5]+par[6]*x[0];
	}
}


//Crystal ball function for signal, parameters are 0:normalization,1:mean,2:sigma,3:n,4:alpha;
Double_t CrystalBall(Double_t *x,Double_t *par) {

	Double_t t = (x[0]-par[1])/par[2];
	if (par[4] < 0) t = -t;

	Double_t absAlpha = fabs((Double_t)par[4]);

	if (t >= -absAlpha) {
		return par[0]*exp(-0.5*t*t);
	}
	else {
		Double_t a =  TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
		Double_t b= par[3]/absAlpha - absAlpha;

		return par[0]*(a/TMath::Power(b - t, par[3]));
	}
}


void Delete(){
	if (fBinsPt) delete fBinsPt;
	if (fPeakRange) delete fPeakRange;
	if (fFitRange) delete fFitRange;
	if (fBGFitRange) delete fBGFitRange;
	if (fBGFitRangeLeft) delete fBGFitRangeLeft;
	if (fMesonPlotRange) delete fMesonPlotRange;
	if (fMesonIntRange) delete fMesonIntRange;
	if (fMesonIntRangeWide) delete fMesonIntRangeWide;
	if (fMesonIntRangeNarrow) delete fMesonIntRangeNarrow;
	if (fMesonMassRange) delete fMesonMassRange;
	if (fMesonFitRange) delete fMesonFitRange;
	if (fMesonWidthRange) delete fMesonWidthRange;
	if (fMesonLambdaTailRange) delete fMesonLambdaTailRange;
	if (fNRebin) delete fNRebin;
	if (fGGYields) delete fGGYields;
	if (fMesonYieldsBackFit) delete fMesonYieldsBackFit;
	if (fBckYields) delete fBckYields;
	if (fMesonYields) delete fMesonYields;
	if (fMesonTrueYields) delete fMesonTrueYields;
	if (fMesonTrueYieldsReweighted) delete fMesonTrueYieldsReweighted;
	if (fMesonYieldsFunc) delete fMesonYieldsFunc;
	if (fMesonYieldsResidualBckFunc) delete fMesonYieldsResidualBckFunc;
	if (fMesonYieldsResidualBckFuncBackFit) delete fMesonYieldsResidualBckFuncBackFit;
	if (fMesonYieldsCorResidualBckFunc) delete fMesonYieldsCorResidualBckFunc;
	if (fMesonYieldsCorResidualBckFuncBackFit) delete fMesonYieldsCorResidualBckFuncBackFit;
	if (fMesonYieldsPerEvent) delete fMesonYieldsPerEvent;
	if (fMesonYieldsPerEventBackFit) delete fMesonYieldsPerEventBackFit;
	if (fMesonMass) delete fMesonMass;
	if (fMesonWidth) delete fMesonWidth;
	if (fMesonSB) delete fMesonSB;
	if (fMesonSign) delete fMesonSign;
	if (fMesonTrueSB) delete fMesonTrueSB;
	if (fMesonTrueSign) delete fMesonTrueSign;
	if (fMesonFWHM) delete fMesonFWHM;
	if (fGGYieldsLeft) delete fGGYieldsLeft;
	if (fBckYieldsLeft) delete fBckYieldsLeft;
	if (fMesonYieldsLeft) delete fMesonYieldsLeft;
	if (fMesonYieldsFuncLeft) delete fMesonYieldsFuncLeft;
	if (fMesonYieldsResidualBckFuncLeft) delete fMesonYieldsResidualBckFuncLeft;
	if (fMesonYieldsCorResidualBckFuncLeft) delete fMesonYieldsCorResidualBckFuncLeft;
	if (fMesonYieldsLeftPerEvent) delete fMesonYieldsLeftPerEvent;
	if (fMesonMassLeft) delete fMesonMassLeft;
	if (fMesonWidthLeft) delete fMesonWidthLeft;
	if (fMesonSBLeft) delete fMesonSBLeft;
	if (fMesonSignLeft) delete fMesonSignLeft;
	if (fMesonFWHMLeft) delete fMesonFWHMLeft;
	if (fGGYieldsNarrow) delete fGGYieldsNarrow;
	if (fBckYieldsNarrow) delete fBckYieldsNarrow;
	if (fMesonYieldsNarrow) delete fMesonYieldsNarrow;
	if (fMesonTrueYieldsNarrow) delete fMesonTrueYieldsNarrow;
	if (fMesonTrueYieldsReweightedNarrow) delete fMesonTrueYieldsReweightedNarrow;
	if (fMesonYieldsFuncNarrow) delete fMesonYieldsFuncNarrow;
	if (fMesonYieldsResidualBckFuncNarrow) delete fMesonYieldsResidualBckFuncNarrow;
	if (fMesonYieldsCorResidualBckFuncNarrow) delete fMesonYieldsCorResidualBckFuncNarrow;
	if (fMesonYieldsPerEventNarrow) delete fMesonYieldsPerEventNarrow;
	if (fMesonSBNarrow) delete fMesonSBNarrow;
	if (fMesonSignNarrow) delete fMesonSignNarrow;
	if (fGGYieldsLeftNarrow) delete fGGYieldsLeftNarrow;
	if (fBckYieldsLeftNarrow) delete fBckYieldsLeftNarrow;
	if (fMesonYieldsLeftNarrow) delete fMesonYieldsLeftNarrow;
	if (fMesonYieldsFuncLeftNarrow) delete fMesonYieldsFuncLeftNarrow;
	if (fMesonYieldsResidualBckFuncLeftNarrow) delete fMesonYieldsResidualBckFuncLeftNarrow;
	if (fMesonYieldsCorResidualBckFuncLeftNarrow) delete fMesonYieldsCorResidualBckFuncLeftNarrow;
	if (fMesonYieldsLeftPerEventNarrow) delete fMesonYieldsLeftPerEventNarrow;
	if (fMesonSBLeftNarrow) delete fMesonSBLeftNarrow;
	if (fMesonSignLeftNarrow) delete fMesonSignLeftNarrow;
	if (fGGYieldsWide) delete fGGYieldsWide;
	if (fBckYieldsWide) delete fBckYieldsWide;
	if (fMesonYieldsWide) delete fMesonYieldsWide;
	if (fMesonTrueYieldsWide) delete fMesonTrueYieldsWide;
	if (fMesonTrueYieldsReweightedWide) delete fMesonTrueYieldsReweightedWide;
	if (fMesonYieldsFuncWide) delete fMesonYieldsFuncWide;
	if (fMesonYieldsResidualBckFuncWide) delete fMesonYieldsResidualBckFuncWide;
	if (fMesonYieldsCorResidualBckFuncWide) delete fMesonYieldsCorResidualBckFuncWide;
	if (fMesonYieldsPerEventWide) delete fMesonYieldsPerEventWide;
	if (fMesonSBWide) delete fMesonSBWide;
	if (fMesonSignWide) delete fMesonSignWide;
	if (fGGYieldsLeftWide) delete fGGYieldsLeftWide;
	if (fBckYieldsLeftWide) delete fBckYieldsLeftWide;
	if (fMesonYieldsLeftWide) delete fMesonYieldsLeftWide;
	if (fMesonYieldsFuncLeftWide) delete fMesonYieldsFuncLeftWide;
	if (fMesonYieldsResidualBckFuncLeftWide) delete fMesonYieldsResidualBckFuncLeftWide;
	if (fMesonYieldsCorResidualBckFuncLeftWide) delete fMesonYieldsCorResidualBckFuncLeftWide;
	if (fMesonYieldsLeftPerEventWide) delete fMesonYieldsLeftPerEventWide;
	if (fMesonSBLeftWide) delete fMesonSBLeftWide;
	if (fMesonSignLeftWide) delete fMesonSignLeftWide;
	if (fMesonYieldsBackFitError) delete fMesonYieldsBackFitError;
	if (fBckYieldsError) delete fBckYieldsError;
	if (fMesonYieldsError) delete fMesonYieldsError;
	if (fMesonYieldsFuncError) delete fMesonYieldsFuncError;
	if (fMesonYieldsResidualBckFuncError) delete fMesonYieldsResidualBckFuncError;
	if (fMesonYieldsResidualBckFuncBackFitError) delete fMesonYieldsResidualBckFuncBackFitError;
	if (fMesonYieldsCorResidualBckFuncError) delete fMesonYieldsCorResidualBckFuncError;
	if (fMesonYieldsCorResidualBckFuncBackFitError) delete fMesonYieldsCorResidualBckFuncBackFitError;
	if (fMesonYieldsPerEventError) delete fMesonYieldsPerEventError;
	if (fMesonYieldsPerEventBackFitError) delete fMesonYieldsPerEventBackFitError;
	if (fMesonMassError) delete fMesonMassError;
	if (fMesonWidthError) delete fMesonWidthError;
	if (fMesonSBError) delete fMesonSBError;
	if (fMesonSignError) delete fMesonSignError;
	if (fMesonTrueSBError) delete fMesonTrueSBError;
	if (fMesonTrueSignError) delete fMesonTrueSignError;
	if (fMesonFWHMError) delete fMesonFWHMError;
	if (fGGYieldsLeftError) delete fGGYieldsLeftError;
	if (fBckYieldsLeftError) delete fBckYieldsLeftError;
	if (fMesonYieldsLeftError) delete fMesonYieldsLeftError;
	if (fMesonYieldsFuncLeftError) delete fMesonYieldsFuncLeftError;
	if (fMesonYieldsResidualBckFuncLeftError) delete fMesonYieldsResidualBckFuncLeftError;
	if (fMesonYieldsCorResidualBckFuncLeftError) delete fMesonYieldsCorResidualBckFuncLeftError;
	if (fMesonYieldsLeftPerEventError) delete fMesonYieldsLeftPerEventError;
	if (fMesonMassLeftError) delete fMesonMassLeftError;
	if (fMesonWidthLeftError) delete fMesonWidthLeftError;
	if (fMesonSBLeftError) delete fMesonSBLeftError;
	if (fMesonSignLeftError) delete fMesonSignLeftError;
	if (fMesonFWHMLeftError) delete fMesonFWHMLeftError;
	if (fGGYieldsNarrowError) delete fGGYieldsNarrowError;
	if (fBckYieldsNarrowError) delete fBckYieldsNarrowError;
	if (fMesonYieldsNarrowError) delete fMesonYieldsNarrowError;
	if (fMesonYieldsFuncNarrowError) delete fMesonYieldsFuncNarrowError;
	if (fMesonYieldsResidualBckFuncNarrowError) delete fMesonYieldsResidualBckFuncNarrowError;
	if (fMesonYieldsCorResidualBckFuncNarrowError) delete fMesonYieldsCorResidualBckFuncNarrowError;
	if (fMesonYieldsPerEventNarrowError) delete fMesonYieldsPerEventNarrowError;
	if (fGGYieldsLeftNarrowError) delete fGGYieldsLeftNarrowError;
	if (fBckYieldsLeftNarrowError) delete fBckYieldsLeftNarrowError;
	if (fMesonYieldsLeftNarrowError) delete fMesonYieldsLeftNarrowError;
	if (fMesonYieldsFuncLeftNarrowError) delete fMesonYieldsFuncLeftNarrowError;
	if (fMesonYieldsResidualBckFuncLeftNarrowError) delete fMesonYieldsResidualBckFuncLeftNarrowError;
	if (fMesonYieldsCorResidualBckFuncLeftNarrowError) delete fMesonYieldsCorResidualBckFuncLeftNarrowError;
	if (fMesonYieldsLeftPerEventNarrowError) delete fMesonYieldsLeftPerEventNarrowError;
	if (fMesonSBLeftNarrowError) delete fMesonSBLeftNarrowError;
	if (fMesonSignLeftNarrowError) delete fMesonSignLeftNarrowError;
	if (fGGYieldsWideError) delete fGGYieldsWideError;
	if (fBckYieldsWideError) delete fBckYieldsWideError;
	if (fMesonYieldsWideError) delete fMesonYieldsWideError;
	if (fMesonYieldsFuncWideError) delete fMesonYieldsFuncWideError;
	if (fMesonYieldsResidualBckFuncWideError) delete fMesonYieldsResidualBckFuncWideError;
	if (fMesonYieldsCorResidualBckFuncWideError) delete fMesonYieldsCorResidualBckFuncWideError;
	if (fMesonYieldsPerEventWideError) delete fMesonYieldsPerEventWideError;
	if (fMesonSBWideError) delete fMesonSBWideError;
	if (fMesonSignWideError) delete fMesonSignWideError;
	if (fGGYieldsLeftWideError) delete fGGYieldsLeftWideError;
	if (fBckYieldsLeftWideError) delete fBckYieldsLeftWideError;
	if (fMesonYieldsLeftWideError) delete fMesonYieldsLeftWideError;
	if (fMesonYieldsFuncLeftWideError) delete fMesonYieldsFuncLeftWideError;
	if (fMesonYieldsResidualBckFuncLeftWideError) delete fMesonYieldsResidualBckFuncLeftWideError;
	if (fMesonYieldsCorResidualBckFuncLeftWideError) delete fMesonYieldsCorResidualBckFuncLeftWideError;
	if (fMesonYieldsLeftPerEventWideError) delete fMesonYieldsLeftPerEventWideError;
	if (fMesonSBLeftWideError) delete fMesonSBLeftWideError;
	if (fMesonSignLeftWideError) delete fMesonSignLeftWideError;
	if (fHistoMappingTrueMesonInvMassPtBins) delete fHistoMappingTrueMesonInvMassPtBins;
	if (fHistoMappingTrueMesonInvMassPtReweightedBins) delete fHistoMappingTrueMesonInvMassPtReweightedBins;
	if (fHistoMappingTrueGGBckInvMassPtBins) delete fHistoMappingTrueGGBckInvMassPtBins;
	if (fHistoMappingTrueContBckInvMassPtBins) delete fHistoMappingTrueContBckInvMassPtBins;
	if (fHistoMappingTrueAllBckInvMassPtBins) delete fHistoMappingTrueAllBckInvMassPtBins;
	if (fHistoMappingGGInvMassPtBin) delete fHistoMappingGGInvMassPtBin;
	if (fHistoMappingBackInvMassPtBin) delete fHistoMappingBackInvMassPtBin;
	if (fHistoMappingBackNormInvMassPtBin) delete fHistoMappingBackNormInvMassPtBin;
	if (fHistoMappingSignalInvMassPtBin) delete fHistoMappingSignalInvMassPtBin;
	if (fHistoMappingRatioSBInvMassPtBin) delete fHistoMappingRatioSBInvMassPtBin;
	if (fFitSignalInvMassPtBin) delete fFitSignalInvMassPtBin;
	if (fFitSignalPeakPosInvMassPtBin) delete fFitSignalPeakPosInvMassPtBin;
	if (fFitBckInvMassPtBin) delete fFitBckInvMassPtBin;
	if (fHistoMappingBackNormInvMassLeftPtBin) delete fHistoMappingBackNormInvMassLeftPtBin;
	if (fHistoMappingSignalInvMassLeftPtBin) delete fHistoMappingSignalInvMassLeftPtBin;
	if (fFitInvMassLeftPtBin) delete fFitInvMassLeftPtBin;
	if (fFitSignalPeakPosInvMassLeftPtBin) delete fFitSignalPeakPosInvMassLeftPtBin;
	if (fFitBckInvMassLeftPtBin) delete fFitBckInvMassLeftPtBin;
	if (fFitRatioInvMassPtBin) delete fFitRatioInvMassPtBin;
	if (fFitWithPol2ForBG) delete fFitWithPol2ForBG;
	if (fMesonFWHMAlpha01) delete fMesonFWHMAlpha01;
	if (fMesonFWHMAlpha01Error) delete fMesonFWHMAlpha01Error;
	if (fHistoWeightsBGZbinVsMbin) delete fHistoWeightsBGZbinVsMbin;
	if (fHistoFillPerEventBGZbinVsMbin) delete fHistoFillPerEventBGZbinVsMbin;
}


void SetCorrectMCHistogrammNames(TString optionEnergyPeriods, TString optionEnergyEnergy){
	cout << "standard MC chosen" << endl;
	ObjectNameTrue 						= "ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt";
	ObjectNameTrueWOWeights 			= "ESD_TrueMotherPiPlPiMiPiZeroWOWeights_InvMass_Pt";	//  n/a
	ObjectNameProfileWeights 			= "ESD_TruePrimaryMotherWeights_InvMass_Pt";			//  n/a
	ObjectNameTrueSec 					= "ESD_TrueSecondaryMother_InvMass_Pt";					//  n/a
	ObjectNameTrueSecFromK0S 			= "ESD_TrueSecondaryMotherFromK0s_InvMass_Pt";			//  n/a
	ObjectNameMCEtaAcc 					= "MC_EtaInAcc_Pt";
	ObjectNameMCOmegaAcc 				= "MC_OmegaInAcc_Pt";
	ObjectNameMCEta 					= "MC_Eta_Pt";
	ObjectNameMCEtaWOWeights 			= "MC_Eta_WOWeights_Pt";								//  n/a
	ObjectNameMCOmega 					= "MC_Omega_Pt";
	ObjectNameMCOmegaWOWeights 			= "MC_Omega_WOWeights_Pt";								//  n/a
	ObjectNameTrueGGBck     			= "ESD_TrueBckGG_InvMass_Pt";							//  n/a even in previous
	ObjectNameTrueContBck 				= "ESD_TrueBckCont_InvMass_Pt";							//  n/a even in previous
	ObjectNameTrueAllBck    			= "ESD_TrueAllCont_InvMass_Pt";							//  n/a even in previous
	ObjectNameTrueCaloPhoton 			= "ESD_TrueMotherCaloPhoton_InvMass_Pt";
	ObjectNameTrueCaloConvPhoton 		= "ESD_TrueMotherCaloConvertedPhoton_InvMass_Pt";
	ObjectNameTrueCaloElectron 			= "ESD_TrueMotherCaloElectron_InvMass_Pt";
	ObjectNameTrueCaloMerged 			= "ESD_TrueMotherCaloMergedCluster_InvMass_Pt";
	ObjectNameTrueCaloEMNonLeading	 	= "ESD_TrueMotherCaloEMNonLeading_InvMass_Pt";
	ObjectNameTrueCaloMergedPartConv	= "ESD_TrueMotherCaloMergedClusterPartConv_InvMass_Pt";
	if (optionEnergyPeriods || optionEnergyEnergy){}
}

/*
void SetCorrectMCHistogrammNames(TString optionEnergyPeriods, TString optionEnergyEnergy){
	cout << "standard MC chosen" << endl;
	ObjectNameTrue 						= "ESD_TruePrimaryMother_InvMass_Pt";
	ObjectNameTrueWOWeights 			= "ESD_TruePrimaryMotherW0Weights_InvMass_Pt";
	ObjectNameProfileWeights 			= "ESD_TruePrimaryMotherWeights_InvMass_Pt";
	ObjectNameTrueSec 					= "ESD_TrueSecondaryMother_InvMass_Pt";
	ObjectNameTrueSecFromK0S 			= "ESD_TrueSecondaryMotherFromK0s_InvMass_Pt";
	ObjectNameMCPi0Acc 					= "MC_Pi0InAcc_Pt";
	ObjectNameMCEtaAcc 					= "MC_EtaInAcc_Pt";
	ObjectNameMCOmegaAcc 				= "MC_OmegaInAcc_Pt";
	ObjectNameMCPi0 					= "MC_Pi0_Pt";
	ObjectNameMCPi0WOWeights 			= "MC_Pi0_WOWeights_Pt";
	ObjectNameMCEta 					= "MC_Eta_Pt";
	ObjectNameMCEtaWOWeights 			= "MC_Eta_WOWeights_Pt";
	ObjectNameMCOmega 					= "MC_Omega_Pt";
	ObjectNameMCOmegaWOWeights 			= "MC_Omega_WOWeights_Pt";
	ObjectNameTrueGGBck     			= "ESD_TrueBckGG_InvMass_Pt";
	ObjectNameTrueContBck 				= "ESD_TrueBckCont_InvMass_Pt";
	ObjectNameTrueAllBck    			= "ESD_TrueAllCont_InvMass_Pt";
	ObjectNameTrueCaloPhoton 			= "ESD_TrueMotherCaloPhoton_InvMass_Pt";
	ObjectNameTrueCaloConvPhoton 		= "ESD_TrueMotherCaloConvertedPhoton_InvMass_Pt";
	ObjectNameTrueCaloElectron 			= "ESD_TrueMotherCaloElectron_InvMass_Pt";
	ObjectNameTrueCaloMerged 			= "ESD_TrueMotherCaloMergedCluster_InvMass_Pt";
	ObjectNameTrueCaloEMNonLeading	 	= "ESD_TrueMotherCaloEMNonLeading_InvMass_Pt";
	ObjectNameTrueCaloMergedPartConv	= "ESD_TrueMotherCaloMergedClusterPartConv_InvMass_Pt";
	if (optionEnergyPeriods || optionEnergyEnergy){}
}
*/
