#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdlib.h>
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
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CalculateGammaToPi0V2.h"
#include "../CommonHeaders/FittingGammaConversion.h"
//#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "TMath.h"
#include "TSpline.h"

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;


void  CalculateGammaToPi0V2(    TString nameFileGamma   = "",
                                TString nameFilePi0     = "",
                                TString cutSel          = "",
                                TString suffix          = "gif",
                                TString nameMeson       = "",
                                TString isMC            = "",
                                TString option          = "",
                                TString fEstimatePileup = "",
                                Int_t mode              = 9
                            ){

	gROOT->Reset();

	StyleSettingsThesis();  
	SetPlotStyle();

	Color_t colorCocktailAllDecay               = kBlack;
	Color_t colorCocktailPi0                    = kRed+2;
	Color_t colorCocktailEta                    = kBlue+1;
	Color_t colorCocktailEtaP                   = kOrange+1;
	Color_t colorCocktailOmega                  = kYellow+2;
	Color_t colorCocktailPhi                    = kViolet;
	Color_t colorCocktailRho0                   = kAzure-2;
	
	cout << Form("%s/%s/%s/GammaToPi0",cutSel.Data(),option.Data(),suffix.Data()) << endl;
	outputDir                                   = Form("%s/%s/%s/GammaToPi0",cutSel.Data(),option.Data(),suffix.Data());
	gSystem->Exec("mkdir "+outputDir);
    
    // old cutnumber seperation
	TString fEventCutSelection                  = "";
	TString fGammaCutSelection                  = "";
	TString fClusterCutSelection                = "";
	TString fElectronCutSelection               = "";
	TString fMesonCutSelection                  = "";
	TString fCutSelection                       = cutSel;
	if (mode == 9){
		ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
		fEventCutSelection                      = fGammaCutSelection(0,7);
		fGammaCutSelection                      = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
		cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
	} else {
		ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
	}
    
	TString centralityCutNumberLow              = fEventCutSelection(GetEventCentralityMinCutPosition(),1);
    //TString centralityCutNumberHigh             = fEventCutSelection(GetEventCentralityMaxCutPosition(),1);
	Int_t centCutNumberI                        = centralityCutNumberLow.Atoi();
    //Int_t centCutNumberBI                       = centralityCutNumberHigh.Atoi();

	TString extraSignals(fEventCutSelection(GetEventRejectExtraSignalsCutPosition(),1));
	TString centrality                          = GetCentralityString(fEventCutSelection);
	FindCocktailFile(fEventCutSelection, fGammaCutSelection, centrality,option);
	TString centralityAdd                       = GetCentralityStringWoPer(fEventCutSelection);

	textPi0New                                  = Form("Gamma_%s",nameMeson.Data());

	if (isMC.CompareTo("kTRUE") == 0){
		textPrefix2                             = "recMC";
		pictDrawingOptions[1]                   = kTRUE;
	} else {
		textPrefix2                             = "data";
		pictDrawingOptions[1]                   = kFALSE;
	}

	TString mesonType;
	if ((nameMeson.CompareTo("Pi0") == 0) || (nameMeson.CompareTo("Pi0EtaBinning") == 0)){
        pictDrawingOptions[3]                   = kTRUE;
		Meson_text                              = "#pi^{0}";
		mesonType                               = "Pi0";
	} else {
		pictDrawingOptions[3]                   = kFALSE;
		Meson_text                              = "#eta";
		mesonType                               = "Eta";
	}
    
    Bool_t kDoPileup                            = kFALSE;
    if (fEstimatePileup.CompareTo("EstimateTrainPileUp") == 0)
        kDoPileup                               = kTRUE;

	One                                         = new TF1("One","1",0,25);
	One->SetLineWidth(1.2);
	One->SetLineColor(1);

	TString nameFinalResDat                     = Form("%s/%s/Gamma_%s_FinalExtraction_%s.dat",cutSel.Data(),option.Data(), textPrefix2.Data(), cutSel.Data());
	fstream fileFinalResults;
	fileFinalResults.open(nameFinalResDat.Data(), ios::out);
	
    // ***************** gamma file **********************************************************
	fileGamma                                   = new TFile(nameFileGamma);
	if(option.CompareTo("PbPb_2.76TeV") == 0){
		histoGammaSpecCorrPurity                = (TH1D*)fileGamma->Get("GammaUnfold");
	} else if (option.CompareTo("pPb_5.023TeV") == 0){
		histoGammaSpecCorrPurity                = (TH1D*)fileGamma->Get("GammaUnfold");
    } else if (option.CompareTo("13TeV") == 0){
        if (kDoPileup)
            histoGammaSpecCorrPurity            = (TH1D*)fileGamma->Get("GammaCorrUnfoldPileUp_Pt");
        else
            histoGammaSpecCorrPurity            = (TH1D*)fileGamma->Get("GammaCorrUnfold_Pt");
    } else {
        histoGammaSpecCorrPurity                = (TH1D*)fileGamma->Get("GammaUnfold");
    }
    if (!histoGammaSpecCorrPurity) {
        cout << "ERROR: GammaCorrUnfold_Pt not in gamma file" << endl;
        return;
    }
    
	histoGammaSpecMCAll                         = (TH1D*)fileGamma->Get("GammaSpecMC");                     // MC input
    histoGammaSpecCorrESDMC                     = (TH1D*)fileGamma->Get("GammaSpecCorrESDMC");              // fully corrected MC rec
    
	histoMCDecaySumGammaPt                      = (TH1D*)fileGamma->Get("MC_DecayGammaAll_Pt");
	histoMCDecayPi0GammaPt                      = (TH1D*)fileGamma->Get("MC_DecayGammaPi0_Pt");
	histoMCDecayEtaGammaPt                      = (TH1D*)fileGamma->Get("MC_DecayGammaEta_Pt");
	histoMCDecayEtapGammaPt                     = (TH1D*)fileGamma->Get("MC_DecayGammaEtap_Pt");
	histoMCDecayOmegaGammaPt                    = (TH1D*)fileGamma->Get("MC_DecayGammaOmega_Pt");
	histoMCDecayRho0GammaPt                     = (TH1D*)fileGamma->Get("MC_DecayGammaRho_Pt");
	histoMCDecayPhiGammaPt                      = (TH1D*)fileGamma->Get("MC_DecayGammaPhi_Pt");
	histMCAllMinusDecay                         = (TH1D*)fileGamma->Get("MC_Direct_or_All_Minus_Decay");

    // ***************** pi0 file ************************************************************
	filePi0                                     = new TFile(nameFilePi0);
    histoCorrectedPi0YieldNormalEff             = (TH1D*)filePi0->Get("CorrectedYieldNormEff");
	histoCorrectedPi0Yield                      = (TH1D*)filePi0->Get("CorrectedYieldTrueEff");
	histoCorrectedPi0YieldWide                  = (TH1D*)filePi0->Get("CorrectedYieldTrueEffWide");
	histoCorrectedPi0YieldNarrow                = (TH1D*)filePi0->Get("CorrectedYieldTrueEffNarrow");
	histoMCYieldMeson                           = (TH1D*)filePi0->Get("MCYield_Meson");
	histoMCYieldMesonOldBin                     = (TH1D*)filePi0->Get("MCYield_Meson_oldBin");

	histoIncRatioPurity                         = (TH1D*) histoGammaSpecCorrPurity->Clone("Rec_IncRatioPurity");
	histoIncRatioPurity->Divide(histoIncRatioPurity,histoCorrectedPi0YieldNormalEff,1,1,"");

	histoIncRatioPurityTrueEff                  = (TH1D*) histoGammaSpecCorrPurity->Clone("IncRatioPurity_trueEff");
	histoIncRatioPurityTrueEff->Divide(histoIncRatioPurityTrueEff,histoCorrectedPi0Yield,1,1,"");

	histoIncRatioPurityTrueEffWide              = (TH1D*) histoGammaSpecCorrPurity->Clone("IncRatioPurity_trueEffWide");
	histoIncRatioPurityTrueEffWide->Divide(histoIncRatioPurityTrueEffWide,histoCorrectedPi0YieldWide,1,1,"");
	
	histoIncRatioPurityTrueEffNarrow            = (TH1D*) histoGammaSpecCorrPurity->Clone("IncRatioPurity_trueEffNarrow");
	histoIncRatioPurityTrueEffNarrow->Divide(histoIncRatioPurityTrueEffNarrow,histoCorrectedPi0YieldNarrow,1,1,"");
	
	histoMCIncRatio                             = (TH1D*) histoGammaSpecMCAll->Clone("MC_IncRatio");
	histoMCIncRatio->Divide(histoGammaSpecMCAll,histoMCYieldMeson,1,1,"");

    histoMCesdIncRatio                          = (TH1D*) histoGammaSpecCorrESDMC->Clone("MCesd_IncRatio");
	histoMCesdIncRatio->Divide(histoGammaSpecCorrESDMC,histoMCYieldMeson,1,1,"");
	
    histoIncRatioGammaMC                        = (TH1D*) histoGammaSpecMCAll->Clone("MC_MCGamma_DataPi0_IncRatio");
	histoIncRatioGammaMC->Divide(histoIncRatioGammaMC,histoCorrectedPi0YieldNormalEff,1,1,"");
	
    cout << __LINE__ << endl;

	// MC Data Ratios
	histoTruePi0MCData                          = (TH1D*) histoMCYieldMeson->Clone("TruePi0MCData");
	histoTruePi0MCData->Divide(histoCorrectedPi0Yield);

	histoNormalPi0MCData                        = (TH1D*) histoMCYieldMeson->Clone("NormalPi0MCData");
	histoNormalPi0MCData->Divide(histoCorrectedPi0YieldNormalEff);

	histoPurityGammaMCData                      = (TH1D*) histoGammaSpecMCAll->Clone("PurityGammaMCData");
	histoPurityGammaMCData->Divide(histoGammaSpecCorrPurity);

    cout << __LINE__ << endl;

	//MC Decay Ratios
	histoDecayRatioSumGamma                     = (TH1D*) histoMCDecaySumGammaPt->Clone("MC_DecayRatio_SumGamma");
    histoDecayRatioSumGamma->Divide(histoMCDecaySumGammaPt,histoMCYieldMesonOldBin,1,1,"");

	histoDecayRatioPi0Gamma                     = (TH1D*) histoMCDecayPi0GammaPt->Clone("MC_DecayRatio_Pi0Gamma");
	histoDecayRatioPi0Gamma->Divide(histoMCDecayPi0GammaPt,histoMCYieldMesonOldBin,1,1,"");

    histoDecayRatioEtaGamma                     = (TH1D*) histoMCDecayEtaGammaPt->Clone("MC_DecayRatio_EtaGamma");
    histoDecayRatioEtaGamma->Divide(histoMCDecayEtaGammaPt,histoMCYieldMesonOldBin,1,1,"");

	histoDecayRatioEtapGamma                    = (TH1D*) histoMCDecayEtapGammaPt->Clone("MC_DecayRatio_EtapGamma");
	histoDecayRatioEtapGamma->Divide(histoMCDecayEtapGammaPt,histoMCYieldMesonOldBin,1,1,"");

	histoDecayRatioOmegaGamma                   = (TH1D*) histoMCDecayOmegaGammaPt->Clone("MC_DecayRatio_OmegaGamma");
	histoDecayRatioOmegaGamma->Divide(histoMCDecayOmegaGammaPt,histoMCYieldMesonOldBin,1,1,"");

	histoDecayRatioRho0Gamma                    = (TH1D*) histoMCDecayRho0GammaPt->Clone("MC_DecayRatio_Rho0Gamma");
	histoDecayRatioRho0Gamma->Divide(histoMCDecayRho0GammaPt,histoMCYieldMesonOldBin,1,1,"");

	histoDecayRatioPhiGamma                     = (TH1D*) histoMCDecayPhiGammaPt->Clone("MC_DecayRatio_PhiGamma");
	histoDecayRatioPhiGamma->Divide(histoMCDecayPhiGammaPt,histoMCYieldMesonOldBin,1,1,"");
	
    cout << __LINE__ << endl;

	//**********************************************************************************
	//******************** NLO Calculatins Ratio ***************************************
	//**********************************************************************************

	Bool_t doNLOComparison                      = kTRUE;
	
	// Get the NLO Results
    // xSection also depends on trigger -> implement decision according to trigger choice
	if(option.CompareTo("7TeV") == 0){
		xSection                                = xSection7TeV;
		fileNameNLOPhotonHalf                   = "ExternalInput/Theory/ALICENLOcalcDirectPhoton7TeVmuhalf.dat";
		fileNameNLOPhotonOne                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton7TeVmu.dat";
		fileNameNLOPhotonTwo                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton7TeVtwomu.dat";
	}
	if(option.CompareTo("8TeV") == 0){
		xSection                                = xSection7TeV;
		fileNameNLOPhotonHalf                   = "ExternalInput/Theory/ALICENLOcalcDirectPhoton8TeVmuhalf.dat";
		fileNameNLOPhotonOne                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton8TeVmu.dat";
		fileNameNLOPhotonTwo                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton8TeVtwomu.dat";
	}
    if(option.CompareTo("13TeV") == 0){
        xSection                                = xSection7TeV;
        fileNameNLOPhotonHalf                   = "ExternalInput/Theory/ALICENLOcalcDirectPhoton7TeVmuhalf.dat";
        fileNameNLOPhotonOne                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton7TeVmu.dat";
        fileNameNLOPhotonTwo                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton7TeVtwomu.dat";
    }
    if(option.CompareTo("PbPb_2.76TeV") == 0 || option.CompareTo("2.76TeV") == 0){
		xSection                                = xSection2760GeV;
		fileNameNLOPhotonHalf                   = "ExternalInput/Theory/ALICENLOcalcDirectPhoton2760GeVmuhalf.dat";
		fileNameNLOPhotonOne                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton2760GeVmu.dat";
		fileNameNLOPhotonTwo                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton2760GeVtwomu.dat";
	}
	if(option.CompareTo("900GeV") == 0){
		xSection                                = xSection900GeV;
		fileNameNLOPhotonHalf                   = "ExternalInput/Theory/ALICENLOcalcDirectPhoton2760GeVmu.dat";
		fileNameNLOPhotonOne                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton2760GeVmu.dat";
		fileNameNLOPhotonTwo                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton2760GeVmu.dat";
	}
	if(option.CompareTo("pPb_5.023TeV") == 0){
		xSection                                = xSection5023GeVppINEL;
		fileNameNLOPhotonHalf                   = "ExternalInput/Theory/ALICENLOcalcDirectPhoton5023GeVHalfMu.dat";
		fileNameNLOPhotonOne                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton5023GeVMu.dat";
		fileNameNLOPhotonTwo                    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton5023GeVTwoMu.dat";
		doNLOComparison                         = kTRUE;
	}

	Double_t etascaling                         = 1./1.6;
	etascaling                                  = 1.;           // why is this done?

	if (doNLOComparison){
        inHalf.open(fileNameNLOPhotonHalf,ios_base::in);
        cout << fileNameNLOPhotonHalf << endl;
		
        while(!inHalf.eof()){
			nlinesNLOHalf++;

			//               Pt                              DirectPhoton           FragmentationPhoton         SumPhoton
			inHalf >> ptNLOPhotonHalf[nlinesNLOHalf] >> muHalfD[nlinesNLOHalf] >> muHalfF[nlinesNLOHalf] >> muHalfDF[nlinesNLOHalf];
			muHalfF[nlinesNLOHalf]              = muHalfF[nlinesNLOHalf]/xSection/recalcBarn*fcmult*(etascaling);
			muHalfD[nlinesNLOHalf]              = muHalfD[nlinesNLOHalf]/xSection/recalcBarn*fcmult*(etascaling);
            muHalfDF[nlinesNLOHalf]             = muHalfDF[nlinesNLOHalf]/xSection/recalcBarn*fcmult*(etascaling);
		}
		inHalf.close();

		inOne.open(fileNameNLOPhotonOne,ios_base::in);
		cout << fileNameNLOPhotonOne << endl;

		while(!inOne.eof()){
			nlinesNLOOne++;
			inOne >> ptNLOPhotonOne[nlinesNLOOne] >> muOneD[nlinesNLOOne] >> muOneF[nlinesNLOOne] >> muOneDF[nlinesNLOOne];
			muOneF[nlinesNLOOne]                = muOneF[nlinesNLOOne]/xSection/recalcBarn*fcmult*(etascaling);
			muOneD[nlinesNLOOne]                = muOneD[nlinesNLOOne]/xSection/recalcBarn*fcmult*(etascaling);
            muOneDF[nlinesNLOOne]               = muOneDF[nlinesNLOOne]/xSection/recalcBarn*fcmult*(etascaling);
		}
		inOne.close();

		inTwo.open(fileNameNLOPhotonTwo,ios_base::in);
		cout << fileNameNLOPhotonTwo << endl;

		while(!inTwo.eof()){
			nlinesNLOTwo++;
			inTwo >> ptNLOPhotonTwo[nlinesNLOTwo] >> muTwoD[nlinesNLOTwo] >> muTwoF[nlinesNLOTwo] >> muTwoDF[nlinesNLOTwo];
			muTwoF[nlinesNLOTwo]                = muTwoF[nlinesNLOTwo]/xSection/recalcBarn*fcmult*(etascaling);
			muTwoD[nlinesNLOTwo]                = muTwoD[nlinesNLOTwo]/xSection/recalcBarn*fcmult*(etascaling);
            muTwoDF[nlinesNLOTwo]               = muTwoDF[nlinesNLOTwo]/xSection/recalcBarn*fcmult*(etascaling);
        }
		inTwo.close();
        
        cout << __LINE__ << endl;

		TCanvas *NLOCalculationcanvas           = GetAndSetCanvas("NLOCalculationcanvas");
		DrawGammaCanvasSettings( NLOCalculationcanvas, 0.1, 0.015, 0.01, 0.09);
		TPad* padNLOHistos                      = new TPad("padNLOHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
		DrawGammaPadSettings( padNLOHistos, 0.1, 0.015, 0.01, 0.);
		padNLOHistos->Draw();
		TPad* padNLORatios                      = new TPad("padNLORatios", "", 0., 0., 1., 0.25,-1, -1, -2);
		DrawGammaPadSettings( padNLORatios, 0.1, 0.015, 0.0, 0.35);
		padNLORatios->Draw();

		padNLOHistos->cd();
		padNLOHistos->SetLogy();

		SetHistogramm(histoMCDecaySumGammaPt,"p_{T} (GeV/c)","#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)",5e-12,250);
		
		graphNLOCalcMuHalf                      = new TGraphErrors(nlinesNLOHalf,ptNLOPhotonHalf,muHalfDF,0,0);
		graphNLOCalcMuHalf->SetName("graphNLOCalcMuHalf");
		graphNLOCalcMuHalf->SetTitle("graphNLOCalcMuHalf");
		graphNLOCalcMuHalf->GetYaxis()->SetRangeUser(0.1e-12,1);
		DrawGammaSetMarkerTGraph(graphNLOCalcMuHalf, 24, 1., kRed+2, kRed+2);
        
        cout << __LINE__ << ": start fitting to NLO calculations" << endl;
        
		//Double_t paramGraphNLOMuHalf[3]       = {0.1,3.4,0.5};
		fitNLOMuHalf                            = FitObject("m","fitNLOMuHalf","Pi0",graphNLOCalcMuHalf,2.,25.);//,paramGraphNLOMuHalf);        // mod power law
		DrawGammaSetMarkerTF1( fitNLOMuHalf, 8, 2.0, kBlue+2);
        
		graphNLOCalcMuOne                       = new TGraphErrors(nlinesNLOOne,ptNLOPhotonOne,muOneDF,0,0);
		graphNLOCalcMuOne->SetName("graphNLOCalcMuOne");
		DrawGammaSetMarkerTGraph(graphNLOCalcMuOne, 27, 1., kRed, kRed);
        
		//Double_t paramGraphNLOMuOne[3]        = {0.1,0.37,5.1};
		//Double_t paramGraphNLOMuOne[3]        = {0.1,3.5,0.5};
		fitNLOMuOne                             = FitObject("m","fitNLOMuOne","Pi0",graphNLOCalcMuOne,2.,25.);//,paramGraphNLOMuOne);
		DrawGammaSetMarkerTF1( fitNLOMuOne, 8, 2.0, kBlue);

		graphNLOCalcMuTwo                       = new TGraphErrors(nlinesNLOTwo,ptNLOPhotonTwo,muTwoDF,0,0);
		graphNLOCalcMuTwo->SetName("graphNLOCalcMuTwo");
		DrawGammaSetMarkerTGraph(graphNLOCalcMuTwo, 30, 1., kRed-2, kRed-2);
        
		//Double_t paramGraphNLOMuTwo[3]        = {0.1,0.37,5.1};
		//Double_t paramGraphNLOMuTwo[3]        = {0.1,3.5,0.5};
		fitNLOMuTwo                             = FitObject("m","fitNLOMuTwo","Pi0",graphNLOCalcMuTwo,2.,25.);//,paramGraphNLOMuTwo);
		DrawGammaSetMarkerTF1( fitNLOMuTwo, 8, 2.0, kBlue-2);
        
        histoMCDecaySumGammaPt->Draw("");
        histoGammaSpecCorrPurity->Draw("same");
        fitNLOMuHalf->Draw("same");
        graphNLOCalcMuHalf->Draw("p,same");
        fitNLOMuOne->Draw("same");
        graphNLOCalcMuOne->Draw("p,same");
        fitNLOMuTwo->Draw("same");
        graphNLOCalcMuTwo->Draw("p,same");
        //histMCAllMinusDecay->Draw("same");

		TLegend* leg_NLOCalculationcanvas       = GetAndSetLegend(0.5,0.5,8);
        leg_NLOCalculationcanvas->AddEntry(graphNLOCalcMuHalf,"Direct Photon NLO Calc #mu = p_{T}/2","lep");
        leg_NLOCalculationcanvas->AddEntry(fitNLOMuHalf,"Fit to Direct Photon NLO Calc #mu = p_{T}/2");
		leg_NLOCalculationcanvas->AddEntry(graphNLOCalcMuOne,"Direct Photon NLO Calc #mu = p_{T}","lep");
        leg_NLOCalculationcanvas->AddEntry(fitNLOMuOne,"Fit to Direct Photon NLO Calc #mu = p_{T}");
		leg_NLOCalculationcanvas->AddEntry(graphNLOCalcMuTwo,"Direct Photon NLO Calc #mu = 2p_{T}","lep");
        leg_NLOCalculationcanvas->AddEntry(fitNLOMuTwo,"Fit to Direct Photon NLO Calc #mu = 2p_{T}");
		//leg_NLOCalculationcanvas->AddEntry(histMCAllMinusDecay,"Not Decay Photons from MC");
		leg_NLOCalculationcanvas->AddEntry(histoMCDecaySumGammaPt,"Decay Photons from MC");
		leg_NLOCalculationcanvas->AddEntry(histoGammaSpecCorrPurity,"Inclusive Photons from data");
        leg_NLOCalculationcanvas->Draw();

		padNLORatios->cd();
		padNLORatios->SetLogy();

		histRatioNLOMuHalf                      = (TH1D*) histoMCDecaySumGammaPt->Clone("histRatioNLOMuHalf");
		histRatioNLOMuOne                       = (TH1D*) histoMCDecaySumGammaPt->Clone("histRatioNLOMuOne");
		histRatioNLOMuTwo                       = (TH1D*) histoMCDecaySumGammaPt->Clone("histRatioNLOMuTwo");

		histRatioNLOMuHalf                      = CalculateHistoRatioToFitNLO(histRatioNLOMuHalf,fitNLOMuHalf,2.);
		histRatioNLOMuOne                       = CalculateHistoRatioToFitNLO(histRatioNLOMuOne,fitNLOMuOne,2.);
		histRatioNLOMuTwo                       = CalculateHistoRatioToFitNLO(histRatioNLOMuTwo,fitNLOMuTwo,2.);

        SetHistogramm(histRatioNLOMuHalf,"p_{T} (GeV/c)","#frac{Decay #gamma}{NLO Direct #gamma}",5,3e2);
        SetHistogramm(histRatioNLOMuOne,"p_{T} (GeV/c)","#frac{Decay #gamma}{NLO Direct #gamma}",5,3e2);
        SetHistogramm(histRatioNLOMuTwo,"p_{T} (GeV/c)","#frac{Decay #gamma}{NLO Direct #gamma}",5,3e2);

		histRatioNLOMuHalf->GetYaxis()->CenterTitle(kTRUE);
		histRatioNLOMuOne->GetYaxis()->CenterTitle(kTRUE);
		histRatioNLOMuTwo->GetYaxis()->CenterTitle(kTRUE);
		histRatioNLOMuHalf->GetYaxis()->SetTitleSize(0.1);
		histRatioNLOMuHalf->GetYaxis()->SetTitleOffset(0.4);
		histRatioNLOMuHalf->GetYaxis()->CenterTitle(kTRUE);

		DrawGammaSetMarker(histRatioNLOMuHalf, 24, 0.5, kRed+2, kRed+2);
		DrawGammaSetMarker(histRatioNLOMuOne, 27, 0.5, kRed, kRed);
		DrawGammaSetMarker(histRatioNLOMuTwo, 30, 0.5, kRed-2, kRed-2);

        histRatioNLOMuHalf->Draw();
        histRatioNLOMuOne->Draw("e1,same");
        histRatioNLOMuTwo->Draw("e1,same");

		NLOCalculationcanvas->Print(Form("%s/%s_NLOCalculations_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

		graphNLOCalcMuHalfRebin                 = RebinNLOGraph(graphNLOCalcMuHalf);
		graphNLOCalcMuOneRebin                  = RebinNLOGraph(graphNLOCalcMuOne);
		graphNLOCalcMuTwoRebin                  = RebinNLOGraph(graphNLOCalcMuTwo);

		histNLOCalcMuTwoRebin                   = GraphToHist(graphNLOCalcMuTwoRebin,25,"histNLOCalcMuTwoRebin");
		histNLOCalcMuOneRebin                   = GraphToHist(graphNLOCalcMuOneRebin,25,"histNLOCalcMuOneRebin");
		histNLOCalcMuHalfRebin                  = GraphToHist(graphNLOCalcMuHalfRebin,25,"histNLOCalcMuHalfRebin");
	}
    
    cout << __LINE__ << endl;

	//**********************************************************************************
	//******************** Inclusive Ratio *********************************************
	//**********************************************************************************
	TCanvas* canvasIncRatioPurity               = GetAndSetCanvas("canvasIncRatioPurity");
	SetHistogramm(histoIncRatioPurity, "p_{T} (GeV/c)", "Ratio Inclusive #gamma/#pi^{0}",0,2);
	DrawGammaSetMarker(histoIncRatioPurity, 20, 2.0, 1, 1);
	histoIncRatioPurity->Draw();
	if(option.CompareTo("PbPb_2.76TeV") == 0) DrawCentrality(0.8,0.9,centrality);
	canvasIncRatioPurity->Print(Form("%s/%s_IncRatioPurity_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

	TCanvas* canvasIncRatioPurityTrueEff        = GetAndSetCanvas("canvasIncRatioPurityTrueEff");
	SetHistogramm(histoIncRatioPurityTrueEff, "p_{T} (GeV/c)", "Ratio Inclusive #gamma/#pi^{0}",0,2);
	DrawGammaSetMarker(histoIncRatioPurityTrueEff, 20, 2.0, 4, 4);
	histoIncRatioPurityTrueEff->Draw();
	canvasIncRatioPurityTrueEff->Print(Form("%s/%s_IncRatioPurity_trueEff_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

	TCanvas* canvasMCIncRatio                   = GetAndSetCanvas("canvasMCIncRatio");
	SetHistogramm(histoMCIncRatio, "p_{T} (GeV/c)", "Ratio Inclusive #gamma/#pi^{0}",0,2);
	SetHistogramm(histoMCesdIncRatio, "p_{T} (GeV/c)", "Ratio Inclusive #gamma/#pi^{0}",0,2);
    DrawGammaSetMarker(histoMCIncRatio, 24, 2.0, 2, 2);
    DrawGammaSetMarker(histoMCesdIncRatio, 20, 2.0, 3, 3);
	histoMCesdIncRatio->Draw();
	histoMCIncRatio->Draw("same");
	if(option.CompareTo("PbPb_2.76TeV") == 0) DrawCentrality(0.8,0.9,centrality);
	canvasMCIncRatio->Print(Form("%s/%s_MC_IncRatio_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

	TCanvas* canvasIncRatioAll                  = GetAndSetCanvas("canvasIncRatioAll");
	histoIncRatioPurity->Draw("e1");
	histoIncRatioPurityTrueEff->Draw("e1,same");
	histoMCIncRatio->Draw("e1,same");
	TLegend* leg_canvasIncRatioAll              = GetAndSetLegend(0.18,0.7,5);
	leg_canvasIncRatioAll->AddEntry(histoIncRatioPurity,"Ratio with normal Eff Purity");
	leg_canvasIncRatioAll->AddEntry(histoIncRatioPurityTrueEff,"Ratio with true Eff Purity");
	leg_canvasIncRatioAll->AddEntry(histoMCIncRatio,"MC Ratio");
	leg_canvasIncRatioAll->Draw();
	if(option.CompareTo("PbPb_2.76TeV") == 0) DrawCentrality(0.8,0.9,centrality);
	canvasIncRatioAll->Print(Form("%s/%s_IncRatio_all_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    
	//**********************************************************************************
	//********************** Charged Pions *********************************************
	//**********************************************************************************
	TFile* fileChargedPionInput;
	Bool_t doChargedPionComp                    = kFALSE;
	TH1D* histoChargedPionSpecStatErr           = NULL;
	TH1D* histoChargedPionSpecSystErr           = NULL;
	
	Double_t *oHagparameters                    = new Double_t[5]; // for pp
	oHagparameters[0]                           = 36;
	oHagparameters[1]                           = 0.378866;
	oHagparameters[2]                           = 0.0784619;
	oHagparameters[3]                           = 0.687722;
	oHagparameters[4]                           = 6.18396;

	TString forOutput                           = "";
	if (option.CompareTo("pPb_5.023TeV")==0){
		fileChargedPionInput                    = new TFile("ExternalInputpPb/ChargedPionSpectrapPb_4_Apr_2014.root");
		histoChargedPionSpecStatErr             = (TH1D*)fileChargedPionInput->Get("histoChargedPionSpecFullPtStatpPb");
		histoChargedPionSpecSystErr             = (TH1D*)fileChargedPionInput->Get("histoChargedPionSpecFullPtSyspPb");

		doChargedPionComp                       = kTRUE;
	}
    
	TSpline3 *splineChargedHigh                 = NULL;
	TF1* fitChargedPions                        = NULL;
	if (doChargedPionComp){
		fitChargedPions                         = FitObject("oHag","chargedPionFit","Pi0",NULL);
		cout << "here" << endl;
// 		fitChargedPions->SetParLimits(0, 0.1, 100);
// 		fitChargedPions->SetParLimits(1, 0.1, 100);
// 		fitChargedPions->SetRange(0.5,15.);
// 		fitChargedPions->SetParLimits(1, par1min3, par1max3);
// 		fitChargedPions->SetParLimits(3, 4, 10);
// 		fitChargedPions->SetParLimits(3, par3min3, par3max3);
// 		fitChargedPions->SetParLimits(4, par4min3, par4max3);

// 		fitChargedPions->FixParameter(1,oHagparameters[1]);
// 		fitChargedPions->FixParameter(2,oHagparameters[2]);
// 		fitChargedPions->FixParameter(3,oHagparameters[3]);
// 		fitChargedPions->FixParameter(4,oHagparameters[4]);
		histoChargedPionSpecStatErr->Fit(fitChargedPions,"INRME+","",0.4,15.);
		cout << "here" << endl;
		
		DrawGammaSetMarkerTF1(fitChargedPions, 1, 1.0, kBlue-2);
		fileFinalResults<<"fitChargedPions "<<endl;
		forOutput= WriteParameterToFile(fitChargedPions);
		fileFinalResults << forOutput.Data() << endl;

// 		splineChargedHigh = new TSpline3(histoChargedPionSpecStatErr);
// 		splineChargedHigh->SetLineColor(kRed);
// 		splineChargedHigh->SetLineStyle(1);
        
		TCanvas *ChargedCanvas                  = GetAndSetCanvas("ChargedCanvas");
		ChargedCanvas->SetLogy();
		ChargedCanvas->SetLogx();
		
		TH2F * histoDummyPlotSpectra            = new TH2F("histoDummyPlotSpectra","histoDummyPlotSpectra",1000,0.09,30.,1000,1e-8,1e5);
		SetStyleHistoTH2ForGraphs(histoDummyPlotSpectra, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#i{y}} (GeV^{-2}#it{c})", 0.032,0.04, 0.04,0.04, 1,1.55);
		histoDummyPlotSpectra->GetYaxis()->SetRangeUser(1e-8,1e3);
		histoDummyPlotSpectra->GetXaxis()->SetRangeUser(0.09,20);
		histoDummyPlotSpectra->DrawCopy(); 
		
		DrawGammaSetMarker(histoChargedPionSpecStatErr, 20, 1.5, kBlue-2, kBlue-2);
		DrawGammaSetMarker(histoCorrectedPi0Yield, 24, 2.5, 1, 1);
		
		histoCorrectedPi0Yield->DrawCopy("same");
		histoChargedPionSpecStatErr->DrawCopy("same");
		
		TLegend* leg_ChargedCanvas              = GetAndSetLegend(0.13,0.12,6);
		leg_ChargedCanvas->AddEntry(histoCorrectedPi0Yield,"Neutral pions","pl");
		leg_ChargedCanvas->AddEntry(histoChargedPionSpecStatErr,"charged pions","pl");
		leg_ChargedCanvas->AddEntry(fitChargedPions,"Fit to charged pions","l");
		leg_ChargedCanvas->Draw();

		histoCorrectedPi0Yield->DrawCopy("same");
		histoChargedPionSpecStatErr->DrawCopy("same");
// 		splineChargedHigh->Draw("same");	
		fitChargedPions->Draw("same");	
		histoCorrectedPi0Yield->DrawCopy("same");
		ChargedCanvas->Print(Form("%s/%s_Spectra_Charged_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
		doChargedPionComp=kFALSE;

        TCanvas* canvasRatioCharged             = new TCanvas("canvasRatioCharged","",200,10,1200,700);  // gives the page size
        DrawGammaCanvasSettings(canvasRatioCharged,  0.08, 0.01, 0.015, 0.13);
        canvasRatioCharged->SetLogx();
        
        TH2F * ratioChargedToFit;
        ratioChargedToFit                       = new TH2F("ratioChargedToFit","ratioChargedToFit",1000,0.23,30.,1000,0.5,1.85);
        SetStyleHistoTH2ForGraphs(ratioChargedToFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.064, 0.05,0.064, 0.8,0.6, 512, 505);
        ratioChargedToFit->DrawCopy();
	
        TH1D* histoRatioChargedToFit            = (TH1D*)histoChargedPionSpecStatErr->Clone("ratiotofit");
        histoRatioChargedToFit                  = CalculateHistoRatioToFit (histoRatioChargedToFit, fitChargedPions);
	
        DrawGammaSetMarker(histoRatioChargedToFit, 20, 1.5, kBlue-2, kBlue-2);
        histoRatioChargedToFit->Draw("same,e1,x0");
	
        DrawGammaLines(0., 30.,1., 1.,0.1);
		
        canvasRatioCharged->Update();
        canvasRatioCharged->Print(Form("%s/%s_Ratio_ChargedToFit_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
	}
    
    cout << __LINE__ << endl;

	//***************************************************************************************************************
	//**************************************** Fitting photon spectrum **********************************************
	//***************************************************************************************************************
	TString fitGammaA                           = "h";      // Hagedorn
	TString fitGammaB                           = "l";      // Levy
	Double_t fitMinPt                           = 0.3;
    Double_t fitMaxPt                           = 16;

	if(option.CompareTo("PbPb_2.76TeV") == 0){
		if(centCutNumberI<4){
			fitGammaA                           = "QCD";
			fitGammaB                           = "oHag";
			fitMaxPt                            = 14;
		} else{
			fitGammaA                           = "QCD";
			fitGammaB                           = "oHag";
			fitMaxPt                            = 14;
		}
	}
	
    ConversionGammaFitA                         = FitObject(fitGammaA,"ConversionGammaFitA","Pi0",histoGammaSpecCorrPurity,fitMinPt,fitMaxPt);
	DrawGammaSetMarkerTF1(ConversionGammaFitA, 1, 2.0, kBlue-2);
	fileFinalResults << "CorrectedYieldTrueEff " << fitGammaA << endl;
	forOutput                                   = WriteParameterToFile(ConversionGammaFitA);
	fileFinalResults << forOutput.Data() << endl;
	
	ConversionGammaFitB                         = FitObject(fitGammaB,"ConversionGammaFitB","Pi0",histoGammaSpecCorrPurity,fitMinPt,fitMaxPt);
	DrawGammaSetMarkerTF1(ConversionGammaFitB, 2, 2.0, kRed-3);
	fileFinalResults << "CorrectedYieldTrueEff " << fitGammaB << endl;
	forOutput                                   = WriteParameterToFile(ConversionGammaFitB);
	fileFinalResults << forOutput.Data() << endl;
    
	ConversionGammaFitB->SetLineColor(2);
	ConversionGammaFitB->SetLineStyle(2);
	
	Int_t canvasWidth                           = 1400;
	Int_t canvasHeight                          = 1500;
	Int_t textSizeLabelsPixel                   = 40;
	Double_t relativeMarginsX                   = 0.12;
	Int_t absmargin                             = relativeMarginsX*canvasWidth;
	
	Double_t scalefacMargin                     = 0.21;
	TCanvas* ConversionGammaCanvas              = GetAndSetCanvas("ConversionGammaSpeccanvas", 0., 0., canvasWidth ,canvasHeight);
	DrawGammaCanvasSettings( ConversionGammaCanvas, relativeMarginsX, 0.015, 0.01, 0.09);
	TPad* padConversionGamma                    = new TPad("padConversionGammaHistos", "", 0., 0.33, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padConversionGamma, relativeMarginsX, 0.015, 0.01, 0.);
	padConversionGamma->Draw();
	TPad* padConversionGammaRatio               = new TPad("padConversionGammaRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
	DrawGammaPadSettings( padConversionGammaRatio, relativeMarginsX, 0.015, 0.0, 0.2);
	padConversionGammaRatio->Draw();
    
    cout << __LINE__ << endl;
	
	//***********************************************************************************************************************************************
	//********************************************** Calculate textsizes ****************************************************************************
	//***********************************************************************************************************************************************
	Double_t textsizeLabelsUpper                = 0;
	Double_t textsizeFacUpper                   = 0;
	Double_t marginOffsetUpper                  = 0;
	if (padConversionGamma->XtoPixel(padConversionGamma->GetX2()) < padConversionGamma->YtoPixel(padConversionGamma->GetY1())){
		textsizeLabelsUpper                     = (Double_t)textSizeLabelsPixel/padConversionGamma->XtoPixel(padConversionGamma->GetX2()) ;
		textsizeFacUpper                        = (Double_t)1./padConversionGamma->XtoPixel(padConversionGamma->GetX2()) ;
		marginOffsetUpper                       = scalefacMargin/(textsizeFacUpper*absmargin);
	} else {
		textsizeLabelsUpper                     = (Double_t)textSizeLabelsPixel/padConversionGamma->YtoPixel(padConversionGamma->GetY1());
		textsizeFacUpper                        = (Double_t)1./padConversionGamma->YtoPixel(padConversionGamma->GetY1());
		marginOffsetUpper                       = scalefacMargin/(textsizeFacUpper*absmargin);
	}
	cout << __LINE__ << ": marginOffsetUpper = " << marginOffsetUpper << endl;

	Double_t textsizeLabelsLower                = 0;
	Double_t textsizeFacLower                   = 0;
	Double_t marginOffsetLower                  = 0;
	if (padConversionGammaRatio->XtoPixel(padConversionGammaRatio->GetX2()) < padConversionGammaRatio->YtoPixel(padConversionGammaRatio->GetY1())){
		textsizeLabelsLower                     = (Double_t)textSizeLabelsPixel/padConversionGammaRatio->XtoPixel(padConversionGammaRatio->GetX2()) ;
		textsizeFacLower                        = (Double_t)1./padConversionGammaRatio->XtoPixel(padConversionGammaRatio->GetX2()) ;
		marginOffsetLower                       = scalefacMargin/(textsizeFacLower*absmargin);
	} else {
		textsizeLabelsLower                     = (Double_t)textSizeLabelsPixel/padConversionGammaRatio->YtoPixel(padConversionGammaRatio->GetY1());
		textsizeFacLower                        = (Double_t)1./padConversionGammaRatio->YtoPixel(padConversionGammaRatio->GetY1());
		marginOffsetLower                       = scalefacMargin/(textsizeFacLower*absmargin);
	}
	//***********************************************************************************************************************************************
	
	padConversionGamma->cd();
	padConversionGamma->SetLogy();
	
	SetStyleHistoTH1ForGraphs(histoGammaSpecCorrPurity, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",textsizeLabelsUpper,textsizeLabelsUpper, textsizeLabelsUpper,textsizeLabelsUpper, 1.,marginOffsetUpper);
	DrawGammaSetMarker(histoGammaSpecCorrPurity, 24, 2.0, kBlack, kBlack);
	
	histoGammaSpecCorrPurity->Draw("e1");
	ConversionGammaFitA->Draw("same");
	ConversionGammaFitB->Draw("Csame");

	TLegend* leg_ConversionGammaSpeccanvas      = GetAndSetLegend(0.5,0.83,2);
	leg_ConversionGammaSpeccanvas->AddEntry(ConversionGammaFitA,fitGammaA);
	leg_ConversionGammaSpeccanvas->AddEntry(ConversionGammaFitB,fitGammaB);
	leg_ConversionGammaSpeccanvas->Draw();

	padConversionGammaRatio->cd();

	TH2F * histo2DSpectraRatio;
	histo2DSpectraRatio                         = new TH2F("histo2DSpectraRatio","histo2DSpectraRatio",1000,0.0,20.,1000,0.5,1.9);
	SetStyleHistoTH2ForGraphs(histo2DSpectraRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Yield}{Fit}", textsizeLabelsLower,textsizeLabelsLower, textsizeLabelsLower,textsizeLabelsLower, 1.,marginOffsetLower);
	histo2DSpectraRatio->GetXaxis()->SetRangeUser(0.,histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
	histo2DSpectraRatio->GetXaxis()->SetTitleFont(62);
	histo2DSpectraRatio->GetYaxis()->SetTitleFont(62);
	histo2DSpectraRatio->GetYaxis()->CenterTitle(kTRUE);
	histo2DSpectraRatio->DrawCopy(); 

	histRatioConversionGammaA                   = (TH1D*) histoGammaSpecCorrPurity->Clone("histRatioConversionGammaA");
    histRatioConversionGammaA                   = CalculateHistoRatioToFit(histRatioConversionGammaA,ConversionGammaFitA);
    histRatioConversionGammaB                   = (TH1D*) histoGammaSpecCorrPurity->Clone("histRatioConversionGammaB");
	histRatioConversionGammaB                   = CalculateHistoRatioToFit(histRatioConversionGammaB,ConversionGammaFitB);

	DrawGammaSetMarker(histRatioConversionGammaA, 20, 1., kBlue-2, kBlue-2);
	DrawGammaSetMarker(histRatioConversionGammaB, 24, 1., kRed-3, kRed-3);
    
    One->Draw();
	histRatioConversionGammaA->Draw("same");
	histRatioConversionGammaB->Draw("same");

	ConversionGammaCanvas->Print(Form("%s/%s_Spectra_ConversionGamma_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    
	//***********************************************************************************************************************************************
	//****************************************** Fitting pi0 and eta if possible ********************************************************************
	//***********************************************************************************************************************************************
	TString fitPi0A                             = "oHag";       // modified Hagedorn
	TString fitPi0B                             = "qcd";        // [0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))
	TString fitPi0C                             = "rad";        // Radoslav

	Double_t *qcdparameters                     = new Double_t[5]; // for pp
	qcdparameters[0]                            = 0.0978884;
	qcdparameters[1]                            = 5.7545;
	qcdparameters[2]                            = -4.71103;
	qcdparameters[3]                            = 0.876807;
	qcdparameters[4]                            = 1.26592;

	fitMinPt                                    = 0.5;
	fitMaxPt                                    = 14;

	if(option.CompareTo("PbPb_2.76TeV") == 0){
		oHagparameters                          = NULL;
		qcdparameters                           = NULL;
		if(centCutNumberI<4){
            fitPi0A                             = "oHag";
            fitPi0B                             = "qcd";
            fitPi0C                             = "rad";
		} else{
            fitPi0A                             = "oHag";
            fitPi0B                             = "qcd";
            fitPi0C                             = "rad";
		}
	} else{
		oHagparameters                          = NULL;
		qcdparameters                           = NULL;
		fitPi0A                                 = "h";
		fitPi0B                                 = "l";
		fitPi0C                                 = "oHag";
	}

	cout<<"-----------------------------------------------------------------"<<endl;
	cout<<"---------------------- Begin Fitting Pi0 ------------------------"<<endl;
	cout<<"-----------------------------------------------------------------"<<endl;

	fitPi0YieldA                                = FitObject(fitPi0A,"fitPi0YieldA","Pi0",histoCorrectedPi0Yield,fitMinPt,fitMaxPt,NULL,"IQNRME+");
	fitPi0YieldA->SetRange(0,16);
    DrawGammaSetMarkerTF1(fitPi0YieldA, 1, 2.0, kBlue-2);
    fileFinalResults << "CorrectedYieldTrueEff hagedorn" << endl;
    forOutput                                   = WriteParameterToFile(fitPi0YieldA);
    fileFinalResults << forOutput.Data() << endl;
    
    cout << __LINE__ << endl;
    
	fitPi0YieldAHighpt                          = FitObject("rad","fitPi0YieldAHighPt","Pi0",histoCorrectedPi0Yield,fitMinPt,fitMaxPt,NULL,"IQNRME+");
	fitPi0YieldAHighpt->SetRange(0,16);
	DrawGammaSetMarkerTF1(fitPi0YieldAHighpt, 1, 2.0, kGreen+2);
	fileFinalResults << "CorrectedYieldTrueEff rad " << endl;
	forOutput                                   = WriteParameterToFile(fitPi0YieldAHighpt);
	fileFinalResults << forOutput.Data() << endl;

    cout << __LINE__ << endl;

	fitPi0YieldB                                = FitObject(fitPi0B,"fitPi0YieldB","Pi0",histoCorrectedPi0Yield,fitMinPt,fitMaxPt,oHagparameters,"IQNRME+");
	DrawGammaSetMarkerTF1(fitPi0YieldB, 2, 2.0, kRed-3);
	fitPi0YieldB->SetRange(0,16);
	fileFinalResults << "CorrectedYieldTrueEff Levy" << endl;
	forOutput                                   = WriteParameterToFile(fitPi0YieldB);
	fileFinalResults << forOutput.Data() << endl;
    
    cout << __LINE__ << endl;

	fitPi0YieldC                                = FitObject(fitPi0C,"fitPi0YieldC","Pi0",histoCorrectedPi0Yield,fitMinPt,fitMaxPt,oHagparameters,"IQNRME+");//,NULL,"IQNRME+");
	DrawGammaSetMarkerTF1(fitPi0YieldC, 3, 2.0, kYellow+2);
	fitPi0YieldC->SetRange(0,16);
	TFitResultPtr FitresultC;
	TString resultC                             = "";
	Int_t iterations                            = 11;
	for(Int_t i = 0; i<iterations;i++){
		cout << "Fit Iteration " << i << endl;
		FitresultC                              = histoCorrectedPi0Yield->Fit(fitPi0YieldC, "SIQNRME+", "", fitMinPt, fitMaxPt);
        resultC                                 = gMinuit->fCstatu;
        if(!resultC.CompareTo("SUCCESSFUL")) break;
		else iterations++;
		if(iterations>600) break;
	}
	fileFinalResults << "CorrectedYieldTrueEff modified Hagedorn" << endl;
	forOutput                                   = WriteParameterToFile(fitPi0YieldC);
	fileFinalResults << forOutput.Data() << endl;
    
	fitPi0YieldCWide                            = FitObject(fitPi0C,"fitPi0YieldCWide","Pi0",histoCorrectedPi0YieldWide,fitMinPt,fitMaxPt,oHagparameters,"IQNRME+");//,NULL,"IQNRME+");
	DrawGammaSetMarkerTF1(fitPi0YieldCWide, 3, 2.0, kYellow+2);
	fitPi0YieldCWide->SetRange(0,16);
	TFitResultPtr FitresultCWide;
	TString resultCWide                         = "";
	iterations                                  = 11;
	for(Int_t i = 0; i<iterations;i++){
		cout << "Fit Iteration " << i << endl;
		FitresultCWide                          = histoCorrectedPi0YieldWide->Fit(fitPi0YieldCWide, "SIQNRME+", "", fitMinPt, fitMaxPt);
		resultCWide                             = gMinuit->fCstatu;
		if(!resultCWide.CompareTo("SUCCESSFUL")) break;
		else iterations++;
		if(iterations>600) break;
	}
	fileFinalResults << "CorrectedYieldTrueEff modified Hagedorn wide" << endl;
	forOutput                                   = WriteParameterToFile(fitPi0YieldCWide);
	fileFinalResults << forOutput.Data() << endl;

    cout << __LINE__ << endl;
    
	fitPi0YieldCNarrow                          = FitObject(fitPi0C,"fitPi0YieldCNarrow","Pi0",histoCorrectedPi0YieldNarrow,fitMinPt,fitMaxPt,oHagparameters,"IQNRME+");//,NULL,"IQNRME+");
	DrawGammaSetMarkerTF1(fitPi0YieldCNarrow, 3, 2.0, kYellow+2);
	fitPi0YieldCNarrow->SetRange(0,16);
	TFitResultPtr FitresultCNarrow;
	TString resultCNarrow                       = "";
	iterations                                  = 11;
	for(Int_t i = 0; i<iterations;i++){
		cout << "Fit Iteration " << i << endl;
		FitresultCNarrow                        = histoCorrectedPi0YieldNarrow->Fit(fitPi0YieldCNarrow, "SIQNRME+", "", fitMinPt, fitMaxPt);
		resultCNarrow                           = gMinuit->fCstatu;
		if(!resultCNarrow.CompareTo("SUCCESSFUL")) break;
		else iterations++;
		if(iterations>600) break;
	}
	fileFinalResults << "CorrectedYieldTrueEff modified Hagedorn wide" << endl;
	forOutput                                   = WriteParameterToFile(fitPi0YieldCNarrow);
	fileFinalResults << forOutput.Data() << endl;
	
	cout<<"------------------------------------------------------------------"<<endl;
	cout<<"----------------------- End Fitting Pi0 --------------------------"<<endl;
	cout<<"------------------------------------------------------------------"<<endl;

    TString cutSelEta                           = fGammaCutSelection.Data();
    TString fitEtaA                             = "h";
    TString fitEtaB                             = "l";
    TString fitEtaC                             = "l";
	
	cout << "THIS HAS STILL TO BE FIXED!!!" << endl;
	//return;
	
	fileFinalResults << cutSelEta.Data() << endl;       // why is there a different cutselection chosen for eta?
	cutSelEta.Replace(39,1,"0");                        // using cutSel instead
	fileFinalResults << cutSelEta.Data() << endl;

    cout << "cutSel             " << cutSel.Data() << endl;
    cout << "fGammaCutSelection " << fGammaCutSelection.Data() << endl;
    cout << "cutSelEta          " << cutSelEta.Data() << endl;
	
	//fileEta = new TFile(Form("%s/%s/Eta_data_GammaConvV1Correction_%s.root",cutSelEta.Data(),option.Data(),cutSelEta.Data()));
    fileEta                                     = new TFile(Form("%s/%s/Eta_data_GammaConvV1Correction_%s.root",cutSel.Data(),option.Data(),cutSel.Data()));

    if (fileEta->IsZombie()) fileEta            = NULL;

    if (fileEta) {
		histoCorrectedEtaYieldNormalEff         = (TH1D*)fileEta->Get("CorrectedYieldNormEff");
		histoCorrectedEtaYield                  = (TH1D*)fileEta->Get("CorrectedYieldTrueEff");
        
        //Double_t binningPi0[50];
        //binningPi0[0] = 0;
		Double_t binningEta[20];
        binningEta[0]                           = 0;
        Double_t additionalBins[5]              = {7,8,10,12,14};
		
		for (Int_t i = 1; i < histoCorrectedPi0Yield->GetNbinsX()+1;i++){
			fileFinalResults << "Pi0 \t"<<i << "\t" <<histoCorrectedPi0Yield->GetXaxis()->GetBinUpEdge(i) << endl;
            //binningPi0[i] = histoCorrectedPi0Yield->GetXaxis()->GetBinUpEdge(i);
		}	
		for (Int_t i = 1; i < histoCorrectedEtaYield->GetNbinsX()+1;i++){
			fileFinalResults << "Eta \t"<< i << "\t" <<histoCorrectedEtaYield->GetXaxis()->GetBinUpEdge(i) << endl;
			binningEta[i]                       = histoCorrectedEtaYield->GetXaxis()->GetBinUpEdge(i);
		}
		
		fileFinalResults << "bin with 6 GeV Eta: " <<  histoCorrectedEtaYield->GetXaxis()->FindBin(6) << endl;
		fileFinalResults << "bin with 6 GeV Pi0: " <<  histoCorrectedPi0Yield->GetXaxis()->FindBin(6) << endl;
		fileFinalResults << "maximum number of bins Pi0: " <<  histoCorrectedPi0Yield->GetNbinsX() << endl;
		fileFinalResults << "additional bins in Pi0: " << histoCorrectedPi0Yield->GetNbinsX()-histoCorrectedPi0Yield->GetXaxis()->FindBin(6)+1 << endl;
        
		Int_t b                                 = 0;
		for (Int_t i = histoCorrectedEtaYield->GetXaxis()->FindBin(6); i < histoCorrectedEtaYield->GetXaxis()->FindBin(6)+5;i++){ 
			binningEta[i]                       = additionalBins[b];
			b++;
		}
		for (Int_t i = 0; i < histoCorrectedEtaYield->GetXaxis()->FindBin(6)+5; i++){
			fileFinalResults << i << "\t"<<binningEta[i] << endl;	
		}
		
		histoCorrectedEtaYieldPatched           = new TH1D("histoCorrectedEtaYieldPatched","",histoCorrectedEtaYield->GetXaxis()->FindBin(6)+4,binningEta);
		histoCorrectedEtaYieldPatched->Sumw2();
		for (Int_t i = 1; i < histoCorrectedEtaYield->GetXaxis()->FindBin(6); i++){
			histoCorrectedEtaYieldPatched->SetBinContent(i,histoCorrectedEtaYield->GetBinContent(i));
			histoCorrectedEtaYieldPatched->SetBinError(i,histoCorrectedEtaYield->GetBinError(i));
		}	
		b                                       = histoCorrectedPi0Yield->GetXaxis()->FindBin(6);
		for (Int_t i = histoCorrectedEtaYield->GetXaxis()->FindBin(6); i < histoCorrectedEtaYieldPatched->GetNbinsX()+1; i++){
			histoCorrectedEtaYieldPatched->SetBinContent(i,histoCorrectedPi0Yield->GetBinContent(b)*0.47);
			histoCorrectedEtaYieldPatched->SetBinError(i,histoCorrectedPi0Yield->GetBinError(b)*0.47);
			b++;
		}
        
        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"---------------------- Begin Fitting Eta ------------------------"<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;
        
		fitEtaYieldA                            = FitObject(fitEtaA,"fitEtaYieldA","Eta",histoCorrectedEtaYield,.0,8.0,NULL,"IQNRME+");
		fitEtaYieldA->SetRange(0,16);
		fileFinalResults << "CorrectedEtaTrueEff hagedorn" << endl;
		forOutput                               = WriteParameterToFile(fitEtaYieldA);
		fileFinalResults << forOutput.Data() << endl;

		fitEtaYieldB                            = FitObject(fitEtaB,"fitEtaYieldB","Eta",histoCorrectedEtaYield,.0,8.0,NULL,"IQNRME+");
		fitEtaYieldB->SetRange(0,16);
		fileFinalResults << "CorrectedEta Levy" << endl;
		forOutput                               = WriteParameterToFile(fitEtaYieldB);
		fileFinalResults << forOutput.Data() << endl;
				
		fitEtaYieldC                            = FitObject(fitEtaC,"fitEtaYieldC","Eta",histoCorrectedEtaYieldPatched,.0,14.0,NULL,"IQNRME+");
		fitEtaYieldC->SetRange(0,16);
		fileFinalResults << "CorrectedEta patched modified Hagedorn" << endl;
		fitEtaYieldB->SetLineColor(2);
		forOutput                               = WriteParameterToFile(fitEtaYieldC);
		fileFinalResults << forOutput.Data() << endl;
		
        fitEtaYieldA->SetLineColor(kRed+2);
        fitEtaYieldA->SetLineStyle(5);
		fitEtaYieldB->SetLineColor(kGreen-7);
		fitEtaYieldB->SetLineStyle(5);
		fitEtaYieldC->SetLineColor(kBlue-7);
		fitEtaYieldC->SetLineStyle(6);
        
        cout<<"------------------------------------------------------------------"<<endl;
        cout<<"----------------------- End Fitting Eta --------------------------"<<endl;
        cout<<"------------------------------------------------------------------"<<endl;
	}
	//***********************************************************************************************************************************************

	filePhi                                     = new TFile("ExternalInputpPb/PhiSpectraMinimumBias_13042014.root");
	if (filePhi->IsZombie())
		filePhi                                 = NULL;
	if (filePhi){
		histoCorrectedPhiYieldStat              = (TH1D*) filePhi->Get("stat_cen0100_phi");
		histoCorrectedPhiYieldSyst              = (TH1D*) filePhi->Get("sys_cen0100_phi");
		fitPhiYield                             = FitObject("l","fitPhiYield","Phi",histoCorrectedPhiYieldStat,0.4,13,NULL,"IQNRME+");
        fitPhiYield->SetParameter(0,0.2);
		fitPhiYield->SetParameter(1,7.67588);
        //fitPhiYield->SetParLimits(1,5,9);
		fitPhiYield->SetParameter(2,0.400);
		histoCorrectedPhiYieldStat->Fit(fitPhiYield, "IQNRME+", "", 0.4, 13);
		DrawGammaSetMarkerTF1(fitPhiYield, 2, 2.0, kViolet-2);
		fitPhiYield->SetRange(0,16);
		fileFinalResults << "Phi Yield Levy" << endl;
		forOutput                               = WriteParameterToFile(fitPhiYield);
		fileFinalResults << forOutput.Data() << endl;
	}	
    
    cout << __LINE__ << endl;
	
	TCanvas *ConversionSpeccanvas               = GetAndSetCanvas("ConversionSpeccanvas",0., 0., canvasWidth ,canvasHeight);
	DrawGammaCanvasSettings( ConversionSpeccanvas, relativeMarginsX, 0.015, 0.01, 0.09);
	TPad* padConversionHistos                   = new TPad("padConversionHistos", "", 0., 0.33, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padConversionHistos, relativeMarginsX, 0.015, 0.01, 0.);
	padConversionHistos->Draw();
	TPad* padConversionRatios                   = new TPad("padConversionRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
	DrawGammaPadSettings( padConversionRatios, relativeMarginsX, 0.015, 0.0, 0.2);
	padConversionRatios->Draw();

	padConversionHistos->cd();
	padConversionHistos->SetLogy();
	
	TH2F * histo2DSpectraPi0;
	histo2DSpectraPi0                           = new TH2F("histo2DSpectraPi0","histo2DSpectraPi0",1000,0.0,20.,1000,2e-9,1000);
	SetStyleHistoTH2ForGraphs(histo2DSpectraPi0, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",textsizeLabelsUpper,textsizeLabelsUpper, textsizeLabelsUpper,textsizeLabelsUpper, 1.,marginOffsetUpper);
	histo2DSpectraPi0->GetXaxis()->SetRangeUser(0.,histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
	histo2DSpectraPi0->GetXaxis()->SetTitleFont(62);
	histo2DSpectraPi0->GetYaxis()->SetTitleFont(62);
	//histo2DSpectraPi0->GetYaxis()->CenterTitle(kTRUE);
	histo2DSpectraPi0->DrawCopy(); 

	DrawGammaSetMarker(histoCorrectedPi0Yield, 20, 2.0, kBlack, kBlack);
	histoCorrectedPi0Yield->Draw("e1,same");
	//fitPi0YieldA->Draw("same");
	//fitPi0YieldAHighpt->Draw("same");
	fitPi0YieldB->Draw("Csame");
	fitPi0YieldC->Draw("same");

	if(fileEta){
		DrawGammaSetMarker(histoCorrectedEtaYield, 20, 2.0, kBlue+2, kBlue+2);
		DrawGammaSetMarker(histoCorrectedEtaYieldPatched, 24, 2.0, kBlue-7, kBlue-7);
        histoCorrectedEtaYield->Draw("e1,same");
        histoCorrectedEtaYieldPatched->Draw("e1,same");
        //fitEtaYieldA->Draw("same");
		fitEtaYieldB->Draw("same");
		fitEtaYieldC->Draw("same");
	}

    TLegend* leg_ConversionSpeccanvas;
    if (fileEta) {
        leg_ConversionSpeccanvas                = GetAndSetLegend(0.45,0.63,7);
    } else {
        leg_ConversionSpeccanvas                = GetAndSetLegend(0.45,0.63,3);
    }
    leg_ConversionSpeccanvas->AddEntry(histoCorrectedPi0Yield,"#pi^{0}","pl");
	//    leg_ConversionSpeccanvas->AddEntry(fitPi0YieldA,Form("Parametrisation A %s Chi2 %.2f",fitPi0A.Data(),fitPi0YieldA->GetChisquare()),"l");
	leg_ConversionSpeccanvas->AddEntry(fitPi0YieldB,Form("Parametrisation B %s #chi^{2}/ndf %.2f",fitPi0B.Data(),fitPi0YieldB->GetChisquare()/fitPi0YieldB->GetNDF()),"l");
	leg_ConversionSpeccanvas->AddEntry(fitPi0YieldC,Form("Parametrisation C %s #chi^{2}/ndf %.2f",fitPi0C.Data(),fitPi0YieldC->GetChisquare()/fitPi0YieldC->GetNDF()),"l");
	//    leg_ConversionSpeccanvas->AddEntry(fitPi0YieldAHighpt,Form("Parametrisation D rad Chi2 %.2f",fitPi0YieldAHighpt->GetChisquare()),"l");
	if(fileEta){
        leg_ConversionSpeccanvas->AddEntry(histoCorrectedEtaYield,"#eta","pl");
        leg_ConversionSpeccanvas->AddEntry(histoCorrectedEtaYield,"#eta patched","pl");
        leg_ConversionSpeccanvas->AddEntry(fitEtaYieldB,Form("Parametrisation Eta B %s #chi^{2}/ndf %.2f",fitEtaB.Data(),fitEtaYieldB->GetChisquare()/fitEtaYieldB->GetNDF()),"l");
		leg_ConversionSpeccanvas->AddEntry(fitEtaYieldC,Form("Parametrisation Eta C %s #chi^{2}/ndf %.2f",fitEtaC.Data(),fitEtaYieldC->GetChisquare()/fitEtaYieldC->GetNDF()),"l");
	}
	leg_ConversionSpeccanvas->Draw();
	
	padConversionRatios->cd();
	histo2DSpectraRatio->DrawCopy();
	
	histRatioConversionPi0A                     = (TH1D*) histoCorrectedPi0Yield->Clone("histRatioConversionPi0A");
	histRatioConversionPi0AHighpt               = (TH1D*) histoCorrectedPi0Yield->Clone("histRatioConversionPi0A");
	histRatioConversionPi0B                     = (TH1D*) histoCorrectedPi0Yield->Clone("histRatioConversionPi0B");
	histRatioConversionPi0C                     = (TH1D*) histoCorrectedPi0Yield->Clone("histRatioConversionPi0C");

    cout << __LINE__ << " division with cocktailPi0 not possible, not known yet" << endl;
    
	histRatioConversionPi0A                     = CalculateHistoRatioToFit(histRatioConversionPi0A,fitPi0YieldA);
	histRatioConversionPi0AHighpt->Divide(histRatioConversionPi0AHighpt,cocktailPi0);
	histRatioConversionPi0B                     = CalculateHistoRatioToFit(histRatioConversionPi0B,fitPi0YieldB);
	histRatioConversionPi0C                     = CalculateHistoRatioToFit(histRatioConversionPi0C,fitPi0YieldC);

    cout << __LINE__ << endl;
    
	DrawGammaSetMarker(histRatioConversionPi0A, 20, 2.,kBlue-2,kBlue-2);
	DrawGammaSetMarker(histRatioConversionPi0AHighpt, 20, 2.,kGreen+2,kGreen+2);
	DrawGammaSetMarker(histRatioConversionPi0B, 24, 2.,kRed-3,kRed-3);
	DrawGammaSetMarker(histRatioConversionPi0C, 25, 2.,kYellow+2,kYellow+2);
	
	histRatioConversionPi0B->SetLineStyle(2);
	histRatioConversionPi0C->SetLineStyle(3);
	
	One->Draw("same");
	//histRatioConversionPi0A->Draw("same");
	//histRatioConversionPi0AHighpt->Draw("same");
	histRatioConversionPi0B->Draw("same");
	histRatioConversionPi0C->Draw("same");

	if(fileEta){
		histRatioConversionEtaA                 = (TH1D*) histoCorrectedEtaYield->Clone("histRatioConversionEtaA");
		histRatioConversionEtaB                 = (TH1D*) histoCorrectedEtaYield->Clone("histRatioConversionEtaB");
		histRatioConversionEtaC                 = (TH1D*) histoCorrectedEtaYield->Clone("histRatioConversionEtaC");

		histRatioConversionEtaA                 = CalculateHistoRatioToFit(histRatioConversionEtaA,fitEtaYieldA);
		histRatioConversionEtaB                 = CalculateHistoRatioToFit(histRatioConversionEtaB,fitEtaYieldB);
		histRatioConversionEtaC                 = CalculateHistoRatioToFit(histRatioConversionEtaC,fitEtaYieldC);
		
        //DrawGammaSetMarker(histRatioConversionEtaA, 20, 2.,kRed+2,kRed+2);
        //histRatioConversionEtaA->SetLineStyle(5);
		DrawGammaSetMarker(histRatioConversionEtaB, 20, 2.,kGreen-7,kGreen-7);
		histRatioConversionEtaB->SetLineStyle(5);
		DrawGammaSetMarker(histRatioConversionEtaC, 24, 2.,kBlue-7,kBlue-7);
		histRatioConversionEtaB->SetLineStyle(6);
		
        //histRatioConversionEtaA->Draw("same");
		histRatioConversionEtaB->Draw("same");
		histRatioConversionEtaC->Draw("same");
	}

	ConversionSpeccanvas->Print(Form("%s/%s_Spectra_ConversionPi0_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    
	// Fitted Double Ratios

	// // ---------------------------------------- Random Error Fit ----------------------------------

	TRandom3 *random                            = new TRandom3(0);
	Int_t nIter                                 = 1;
    //Int_t nIterPi0High                          = 0;
    //Int_t nIterPi0Low                           = 0;

	TH1D *LowErrorPi0Yield                      = (TH1D*) histoCorrectedPi0Yield->Clone("LowErrorPi0Yield");
	TH1D *HighErrorPi0Yield                     = (TH1D*) histoCorrectedPi0Yield->Clone("HighErrorPi0Yield");

	Double_t highParameterPi0A                  = 0;
	Double_t highParameterPi0B                  = 0;
	Double_t highParameterPi0C                  = 0;
	Double_t highParameterPi0D                  = 0;
	Double_t highParameterPi0E                  = 0;

	Double_t lowParameterPi0A                   = 0;
	Double_t lowParameterPi0B                   = 0;
	Double_t lowParameterPi0C                   = 0;
	Double_t lowParameterPi0D                   = 0;
	Double_t lowParameterPi0E                   = 0;

	for(Int_t it = 1; it<nIter+1; it++){
		for(Int_t i = 0; i<histoCorrectedPi0Yield->GetNbinsX(); i++){
			//Double_t randomDouble               = random->Gaus(0,histoCorrectedPi0Yield->GetBinError(i+1));
			Double_t randomDouble               = random->Rndm(0)*0.25;
			//if(TMath::Abs(randomDouble)>histoCorrectedPi0Yield->GetBinError(i+1)*0.5) randomDouble = 0;
			LowErrorPi0Yield->SetBinContent(i+1, LowErrorPi0Yield->GetBinContent(i+1) - (randomDouble)*LowErrorPi0Yield->GetBinError(i+1) );
			HighErrorPi0Yield->SetBinContent(i+1, HighErrorPi0Yield->GetBinContent(i+1) + (randomDouble)*HighErrorPi0Yield->GetBinError(i+1) );
		}
        
		FitRandomLow                            = FitObject(fitPi0C,"FitRandomLow","Pi0",LowErrorPi0Yield,fitMinPt,16);
		FitRandomHigh                           = FitObject(fitPi0C,"FitRandomHigh","Pi0",HighErrorPi0Yield,fitMinPt,16);

		lowParameterPi0A                        = lowParameterPi0A + FitRandomLow->GetParameter(0);
		lowParameterPi0B                        = lowParameterPi0B + FitRandomLow->GetParameter(1);
		lowParameterPi0C                        = lowParameterPi0C + FitRandomLow->GetParameter(2);
		lowParameterPi0D                        = lowParameterPi0D + FitRandomLow->GetParameter(3);
		lowParameterPi0E                        = lowParameterPi0E + FitRandomLow->GetParameter(4);

		highParameterPi0A                       = highParameterPi0A + FitRandomHigh->GetParameter(0);
		highParameterPi0B                       = highParameterPi0B + FitRandomHigh->GetParameter(1);
		highParameterPi0C                       = highParameterPi0C + FitRandomHigh->GetParameter(2);
		highParameterPi0D                       = highParameterPi0D + FitRandomHigh->GetParameter(3);
		highParameterPi0E                       = highParameterPi0E + FitRandomHigh->GetParameter(4);
	}

    cout << __LINE__ << endl;

	FitRandomLow->SetParameter(0,lowParameterPi0A/nIter);
	FitRandomLow->SetParameter(1,lowParameterPi0B/nIter);
	FitRandomLow->SetParameter(2,lowParameterPi0C/nIter);
	FitRandomLow->SetParameter(3,lowParameterPi0D/nIter);
	FitRandomLow->SetParameter(4,lowParameterPi0E/nIter);

	FitRandomHigh->SetParameter(0,highParameterPi0A/nIter);
	FitRandomHigh->SetParameter(1,highParameterPi0B/nIter);
	FitRandomHigh->SetParameter(2,highParameterPi0C/nIter);
	FitRandomHigh->SetParameter(3,highParameterPi0D/nIter);
	FitRandomHigh->SetParameter(4,highParameterPi0E/nIter);

	histoIncRatioFitPurity                      = (TH1D*) histoIncRatioPurityTrueEff->Clone("histoIncRatioFitPurity");
	histoIncRatioFitPurityWide                  = (TH1D*) histoIncRatioPurityTrueEffWide->Clone("histoIncRatioFitPurityWide");
	histoIncRatioFitPurityNarrow                = (TH1D*) histoIncRatioPurityTrueEffNarrow->Clone("histoIncRatioFitPurityNarrow");

	histoIncRatioHighFitPurity                  = (TH1D*) histoIncRatioPurityTrueEff->Clone("histoIncRatioHighFitPurity");
	histoIncRatioLowFitPurity                   = (TH1D*) histoIncRatioPurityTrueEff->Clone("histoIncRatioLowFitPurity");
	histoCorrectedPi0YieldFit                   = (TH1D*) histoIncRatioPurityTrueEff->Clone("CorrectedYieldTrueEffPi0Fit");
	histoCorrectedPi0YieldFitWide               = (TH1D*) histoIncRatioPurityTrueEff->Clone("CorrectedYieldTrueEffPi0FitWide");
	histoCorrectedPi0YieldFitNarrow             = (TH1D*) histoIncRatioPurityTrueEff->Clone("CorrectedYieldTrueEffPi0FitNarrow");

	for(Int_t bin = 1; bin<histoIncRatioFitPurity->GetNbinsX()+1; bin++){

		Double_t ptStart                        = histoIncRatioPurityTrueEff->GetBinLowEdge(bin);
		Double_t ptEnd                          = histoIncRatioPurityTrueEff->GetBinLowEdge(bin+1);
		Double_t binWidth                       = histoIncRatioPurityTrueEff->GetBinWidth(bin);

		histoCorrectedPi0YieldFit->SetBinContent(bin,fitPi0YieldC->Eval(histoIncRatioPurityTrueEff->GetBinCenter(bin)));
		histoCorrectedPi0YieldFit->SetBinError(bin,(fitPi0YieldC->IntegralError(ptStart, ptEnd, FitresultC->GetParams(), FitresultC->GetCovarianceMatrix().GetMatrixArray() ) / binWidth));
		histoCorrectedPi0YieldFitWide->SetBinContent(bin,fitPi0YieldCWide->Eval(histoIncRatioPurityTrueEffWide->GetBinCenter(bin)));
		histoCorrectedPi0YieldFitWide->SetBinError(bin,(fitPi0YieldCWide->IntegralError(ptStart, ptEnd, FitresultCWide->GetParams(), FitresultCWide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth));
		histoCorrectedPi0YieldFitNarrow->SetBinContent(bin,fitPi0YieldCNarrow->Eval(histoIncRatioPurityTrueEffNarrow->GetBinCenter(bin)));
		histoCorrectedPi0YieldFitNarrow->SetBinError(bin,(fitPi0YieldCNarrow->IntegralError(ptStart, ptEnd, FitresultCNarrow->GetParams(), FitresultCNarrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth));
		
		histoIncRatioFitPurity->SetBinContent(bin,fitPi0YieldC->Eval(histoIncRatioPurityTrueEff->GetBinCenter(bin)));
		histoIncRatioFitPurity->SetBinError(bin,histoCorrectedPi0Yield->GetBinError(bin));
		histoIncRatioFitPurityWide->SetBinContent(bin,fitPi0YieldCWide->Eval(histoIncRatioPurityTrueEffWide->GetBinCenter(bin)));
		histoIncRatioFitPurityWide->SetBinError(bin,histoCorrectedPi0YieldWide->GetBinError(bin));
		histoIncRatioFitPurityNarrow->SetBinContent(bin,fitPi0YieldCNarrow->Eval(histoIncRatioPurityTrueEffNarrow->GetBinCenter(bin)));
		histoIncRatioFitPurityNarrow->SetBinError(bin,histoCorrectedPi0YieldNarrow->GetBinError(bin));

		histoIncRatioLowFitPurity->SetBinContent(bin,FitRandomLow->Eval(histoIncRatioPurityTrueEff->GetBinCenter(bin)));
		histoIncRatioLowFitPurity->SetBinError(bin,histoIncRatioPurityTrueEff->GetBinError(bin));
		histoIncRatioHighFitPurity->SetBinContent(bin,FitRandomHigh->Eval(histoIncRatioPurityTrueEff->GetBinCenter(bin)));
		histoIncRatioHighFitPurity->SetBinError(bin,histoIncRatioPurityTrueEff->GetBinError(bin));
	}
	
	DrawGammaSetMarker(histoIncRatioPurityTrueEff, 20, 2.0, 1, 1);
	
	histoIncRatioFitPurity->Divide(histoGammaSpecCorrPurity,histoIncRatioFitPurity);
	histoIncRatioFitPurityWide->Divide(histoGammaSpecCorrPurity,histoIncRatioFitPurityWide);
	histoIncRatioFitPurityNarrow->Divide(histoGammaSpecCorrPurity,histoIncRatioFitPurityNarrow);
	histoIncRatioLowFitPurity->Divide(histoGammaSpecCorrPurity,histoIncRatioLowFitPurity);
	histoIncRatioHighFitPurity->Divide(histoGammaSpecCorrPurity,histoIncRatioHighFitPurity);

	SetHistogramm(histoIncRatioFitPurity,"p_{T} (GeV/c)", "Ratio Inclusive #gamma/#pi^{0}",0.0,3.0);
	
	TCanvas *canvasIncRatioFit                  = GetAndSetCanvas("IncRatioFit");
	DrawGammaSetMarker(histoIncRatioFitPurity, 20, 2.0, 7, 7);

	histoIncRatioPurityTrueEff->DrawCopy("e1");
	histoIncRatioFitPurity->DrawCopy("same,e1");
	
    TLegend* leg_canvasIncRatioFit              = GetAndSetLegend(0.4,0.8,3);
	leg_canvasIncRatioFit->AddEntry(histoIncRatioPurityTrueEff,"Inclusive Ratio");
	leg_canvasIncRatioFit->AddEntry(histoIncRatioFitPurity,"Ratio #pi^{0} Fit");
	leg_canvasIncRatioFit->Draw();

	if(option.CompareTo("PbPb_2.76TeV") == 0) DrawCentrality(0.8,0.9,centrality);

	canvasIncRatioFit->Print(Form("%s/%s_IncRatio_Fit_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;

	//**********************************************************************************
	//******************** Comparison MC / Data ****************************************
	//**********************************************************************************

	TCanvas* canvasComparisonMCData             = GetAndSetCanvas("canvasComparisonMCData");
	DrawGammaSetMarker(histoTruePi0MCData, 22, 2.0, kOrange-2, kOrange-2);
	DrawGammaSetMarker(histoNormalPi0MCData, 22, 2.0, kOrange+2, kOrange+2);
	DrawGammaSetMarker(histoPurityGammaMCData, 22, 2.0, kBlue+2, kBlue+2);

	SetHistogramm( histoTruePi0MCData, "p_{T} (GeV/c)", "Comparison MC / Data",0.0,7.0);
	SetHistogramm( histoNormalPi0MCData, "p_{T} (GeV/c)", "Comparison MC / Data",0.0,2.0);
	SetHistogramm( histoPurityGammaMCData, "p_{T} (GeV/c)", "Comparison MC / Data",0.0,7.0);
	histoTruePi0MCData->Draw("e1");
	histoNormalPi0MCData->Draw("e1,same");
	histoPurityGammaMCData->Draw("e1,same");

	TLegend* leg_ComparisonMCData               = GetAndSetLegend(0.12,0.7,4);
	leg_ComparisonMCData->AddEntry(histoTruePi0MCData,"Data #pi^{0} True Eff/ MC #pi^{0}");
	leg_ComparisonMCData->AddEntry(histoNormalPi0MCData,"Data #pi^{0} Normal Eff/ MC #pi^{0}");
	leg_ComparisonMCData->AddEntry(histoPurityGammaMCData,"Data #gamma Purity / MC #gamma");
	leg_ComparisonMCData->Draw();

	canvasComparisonMCData->Print(Form("%s/%s_ComparisonMCData_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;
	
	//**********************************************************************************
	//******************** MC Decay Ratios *********************************************
	//**********************************************************************************

	TCanvas* canvasDecayGammaRatioMC            = GetAndSetCanvas("canvasDecayGammaRatioMC");
	canvasDecayGammaRatioMC->SetLogy();

    SetHistogramm(histoDecayRatioSumGamma,"p_{T} (GeV/c)","Ratio #gamma from Decay / #pi^{0}",1e-5, 1e2);
    SetHistogramm(histoDecayRatioPi0Gamma,"p_{T} (GeV/c)","Ratio #gamma from Decay / #pi^{0}");
    SetHistogramm(histoDecayRatioEtaGamma,"p_{T} (GeV/c)","Ratio #gamma from Decay / #pi^{0}");
    SetHistogramm(histoDecayRatioEtapGamma,"p_{T} (GeV/c)","Ratio #gamma from Decay / #pi^{0}");
    SetHistogramm(histoDecayRatioOmegaGamma,"p_{T} (GeV/c)","Ratio #gamma from Decay / #pi^{0}");
    SetHistogramm(histoDecayRatioPhiGamma,"p_{T} (GeV/c)","Ratio #gamma from Decay / #pi^{0}");
    SetHistogramm(histoDecayRatioRho0Gamma,"p_{T} (GeV/c)","Ratio #gamma from Decay / #pi^{0}");

    histoDecayRatioSumGamma->GetXaxis()->SetRangeUser(0.0,(histoGammaSpecCorrPurity->GetXaxis())->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));

	DrawGammaSetMarker(histoDecayRatioSumGamma, 22, 2.0, colorCocktailAllDecay, colorCocktailAllDecay);
	DrawGammaSetMarker(histoDecayRatioPi0Gamma, 22, 2.0,colorCocktailPi0,colorCocktailPi0);
	DrawGammaSetMarker(histoDecayRatioEtaGamma, 22, 2.0, colorCocktailEta, colorCocktailEta);
	DrawGammaSetMarker(histoDecayRatioEtapGamma, 22, 2.0, colorCocktailEtaP,colorCocktailEtaP);
	DrawGammaSetMarker(histoDecayRatioOmegaGamma, 22, 2.0, colorCocktailOmega,colorCocktailOmega);
	DrawGammaSetMarker(histoDecayRatioPhiGamma, 22, 2.0,colorCocktailPhi,colorCocktailPhi);
	DrawGammaSetMarker(histoDecayRatioRho0Gamma, 22, 2.0, colorCocktailRho0,colorCocktailRho0);
	
	histoDecayRatioSumGamma->Draw("chist");
	histoDecayRatioPi0Gamma->Draw("chistsame");
	histoDecayRatioEtaGamma->Draw("chistsame");
	histoDecayRatioEtapGamma->Draw("chistsame");
	histoDecayRatioOmegaGamma->Draw("chistsame");
	histoDecayRatioRho0Gamma->Draw("chistsame");
	histoDecayRatioPhiGamma->Draw("chistsame");

	cocktail                                    = new TLatex(0.45,0.9,"all decay #gamma");
	SetStyleTLatex( cocktail, 0.04,4,colorCocktailAllDecay,42);
	cocktail->Draw();
	tpi                                         = new TLatex(0.45,0.855,"#pi^{0} #rightarrow #gamma#gamma (e^{+}e^{-}#gamma)");
	SetStyleTLatex( tpi, 0.04,4,colorCocktailPi0,42);
	tpi->Draw();   
	teta                                        = new TLatex(0.45,0.81,"#eta #rightarrow #gamma#gamma (#pi^{+}#pi^{-}#gamma,e^{+}e^{-}#gamma,#pi^{0}#gamma#gamma)");
	SetStyleTLatex( teta, 0.04,4,colorCocktailEta,42);
	teta->Draw();
	tomega                                      = new TLatex(0.45,0.765,"#omega #rightarrow #pi^{0}#gamma (#eta#gamma)");
	SetStyleTLatex( tomega, 0.04,4,colorCocktailOmega,42);
	tomega->Draw();
	tetaprime                                   = new TLatex(0.75,0.855,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)");
	SetStyleTLatex( tetaprime, 0.04,4,colorCocktailEtaP,42);
	tetaprime->Draw();
	tphi                                        = new TLatex(0.75,0.81,"#phi #rightarrow #eta#gamma (#pi^{0}#gamma, #omega#gamma)");
	SetStyleTLatex( tphi, 0.04,4,colorCocktailPhi,42);
	tphi->Draw();
	trho                                        = new TLatex(0.75,0.765,"#rho #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)");
	SetStyleTLatex( trho, 0.04,4,colorCocktailRho0,42);
	trho->Draw();

	canvasDecayGammaRatioMC->Print(Form("%s/%s_MC_DecayRatio_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    
	//**********************************************************************************
	//*********************** Cocktail *************************************************
	//**********************************************************************************
	cocktailFile                                = new TFile(Form("CocktailInput/%s.root",fcocktailFunc.Data()));
    cocktailDir                                 = (TDirectoryFile*) cocktailFile->Get("Combined_7TeV_Hagedorn");
    
    if (fcocktailFuncMtScaledEta.CompareTo("") != 0) {
        cocktailFileMtScaledEta                 = new TFile(Form("CocktailInput/%s.root",fcocktailFuncMtScaledEta.Data()));
        cocktailDirMtScaledEta                  = (TDirectoryFile*) cocktailFileMtScaledEta->Get(fcocktailFuncMtScaledEta);
    } else {
        cout << "fcocktailFuncMtScaledEta not defined" << endl;
        cocktailFileMtScaledEta                 = NULL;
    }
    
    if (fcocktailFuncChargedPions.CompareTo("") != 0) {
        cocktailFileChargedPions                = new TFile(Form("CocktailInput/%s.root",fcocktailFuncChargedPions.Data()));
        cocktailDirChargedPions                 = (TDirectoryFile*) cocktailFileChargedPions->Get(fcocktailFuncChargedPions);
    } else {
        cout << "fcocktailFuncChargedPions not defined" << endl;
        cocktailFileChargedPions                = NULL;
    }
    
	cout<<"**********************************************************************************"<<endl;
	cout<<"**********************************************************************************"<<endl;
	cout<<fcocktailFunc<<endl;
	cout<<"**********************************************************************************"<<endl;
	cout<<"**********************************************************************************"<<endl;

	cocktailPi0                                 = (TH1D* )cocktailDir->Get("ptPion");
	cocktailEta                                 = (TH1D* )cocktailDir->Get("ptEta");
	cocktailPhi                                 = (TH1D* )cocktailDir->Get("ptPhi");

    if (cocktailPi0) cout << "found cocktailPi0" << endl;
    if (cocktailEta) cout << "found cocktailEta" << endl;
    if (cocktailPhi) cout << "found cocktailPhi" << endl;
    
    TCanvas *cocktailCanvasMesonSpec            = GetAndSetCanvas("cocktailCanvasMesonSpec");
	cocktailCanvasMesonSpec->SetLogy();
	//cocktailCanvasMesonSpec->SetGridx();
	//cocktailCanvasMesonSpec->SetGridy();
	
	//if(option.CompareTo("PbPb_2.76TeV") == 0) SetHistogramm(cocktailPi0,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)",5e-8,1e3);
	//else SetHistogramm(cocktailPi0,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)",5e-10,10);

	SetHistogramm(cocktailEta,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
	SetHistogramm(cocktailPhi,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");

    DrawGammaSetMarker(cocktailEta, 20, 2.0,colorCocktailEta,colorCocktailEta);
    DrawGammaSetMarker(cocktailPhi, 20, 2.0,colorCocktailPhi,colorCocktailPhi);
    
	DrawAutoGammaMesonHistos( cocktailPi0, 
							"", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)", 
							kTRUE, 2.,5e-10, kFALSE,
							kFALSE, -0.004, 0.020, 
							kTRUE, 0.,(histoCorrectedPi0Yield->GetXaxis())->GetBinUpEdge(histoCorrectedPi0Yield->GetNbinsX()));
	DrawGammaSetMarker(cocktailPi0, 20, 2.0,colorCocktailPi0,colorCocktailPi0);
	
    cout << __LINE__ << endl;
	
	cocktailPi0->Draw("chist");
	cocktailPi0                                 = RebinTH1D(cocktailPi0,histoCorrectedPi0Yield);
	cocktailPi0->DrawCopy("same");
	cocktailEta->DrawCopy("csamehist");
	cocktailPhi->DrawCopy("csamehist");
	fitPi0YieldC->DrawCopy("same");
    
	// explanation of points
	DrawGammaSetMarker(histoCorrectedPi0Yield, 4, 1., 1, 1);
	histoCorrectedPi0Yield->Draw("e1,same");
	if (fileEta)histoCorrectedEtaYield->Draw("e1,same");
	if (filePhi){
		if (fitPhiYield) fitPhiYield->Draw("same");
		DrawGammaSetMarker(histoCorrectedPhiYieldStat, 24, 2.0,colorCocktailPhi,colorCocktailPhi);
		histoCorrectedPhiYieldStat->Draw("e1,same");
	}	
	TLatex* textPi0Measured                     = new TLatex(4.0,0.002,"#pi^{0} from ALICE Data");
	textPi0Measured->SetTextColor(kBlack);
	textPi0Measured->SetTextSize(0.035);
	textPi0Measured->SetTextFont(42);
	textPi0Measured->Draw();

	TLatex* textEtaMeasured                     = new TLatex(1.0,0.0001,"#eta from ALICE Data");
	textEtaMeasured->SetTextFont(42);
	textEtaMeasured->SetTextColor(kBlue+2);
	textEtaMeasured->SetTextSize(0.035);
	textEtaMeasured->Draw();
	
	TLatex* textPi0Cocktail                     = new TLatex(1.,1.5,"#pi^{0}");
	textPi0Cocktail->SetTextColor(colorCocktailPi0);
	textPi0Cocktail->SetTextSize(0.04);
	textPi0Cocktail->SetTextFont(42);
	textPi0Cocktail->Draw();

	TLatex* textEtaCocktail                     = new TLatex(.7,0.005,"#eta");
	textEtaCocktail->SetTextColor(colorCocktailEta);
	textEtaCocktail->SetTextSize(0.04);
	textEtaCocktail->SetTextFont(42);
	textEtaCocktail->Draw();
	
    TLegend* leg_MesonSpectra;
    if (cocktailEta && cocktailPhi) {
        leg_MesonSpectra                        = GetAndSetLegend(0.7,0.75,5);
    } else {
        leg_MesonSpectra                        = GetAndSetLegend(0.7,0.75,3);
    }
	leg_MesonSpectra->AddEntry(cocktailPi0,"Cocktail #pi^{0}");
    leg_MesonSpectra->AddEntry(histoCorrectedPi0Yield,"#pi^{0} data");
    leg_MesonSpectra->AddEntry(fitPi0YieldC,"fit to #pi^{0} data","l");
    if (cocktailEta) {
        leg_MesonSpectra->AddEntry(cocktailEta,"Cocktail #eta");
        leg_MesonSpectra->AddEntry(histoCorrectedEtaYield,"#eta data");
    }
    if (cocktailPhi) {
        leg_MesonSpectra->AddEntry(histoCorrectedPhiYieldStat,"#phi data");
        if (fitPhiYield) leg_MesonSpectra->AddEntry(fitPhiYield,"fit to #phi data", "l");
    }
	leg_MesonSpectra->Draw();

	cocktailCanvasMesonSpec->Print(Form("%s/%s_Cocktail_MesonSpectra_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    
    //******************************************************************************************************************
	//**********************************Plotting Cocktail Spectra ******************************************************
	//******************************************************************************************************************
	cocktailAllGamma                            = (TH1D*) cocktailDir->Get("ptg");
	cocktailPi0Gamma                            = (TH1D*) cocktailDir->Get("ptgPion");
	cocktailEtaGamma                            = (TH1D*) cocktailDir->Get("ptgEta");
	cocktailEtapGamma                           = (TH1D*) cocktailDir->Get("ptgEtaprime");
	cocktailOmegaGamma                          = (TH1D*) cocktailDir->Get("ptgOmega");
	cocktailPhiGamma                            = (TH1D*) cocktailDir->Get("ptgPhi");
	cocktailRhoGamma                            = (TH1D*) cocktailDir->Get("ptgRho");
	//cocktailSigmaGamma                        = (TH1D*) cocktailDir->Get("ptgSigma");
    
	TCanvas *cocktailCanvasSpec                 = GetAndSetCanvas("cocktailCanvasSpec");
	cocktailCanvasSpec->SetLogy();
	SetHistogramm(cocktailAllGamma,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)",5e-13,1e3);
    SetHistogramm(cocktailPi0Gamma,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
    SetHistogramm(cocktailEtaGamma,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
    SetHistogramm(cocktailEtapGamma,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
    SetHistogramm(cocktailOmegaGamma,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
    SetHistogramm(cocktailPhiGamma,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
    SetHistogramm(cocktailRhoGamma,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
    //SetHistogramm(cocktailSigmaGamma,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");

    DrawGammaSetMarker(cocktailAllGamma, 20, 2.0,colorCocktailAllDecay,colorCocktailAllDecay);
	DrawGammaSetMarker(cocktailPi0Gamma, 20, 2.0,colorCocktailPi0,colorCocktailPi0);
	DrawGammaSetMarker(cocktailEtaGamma, 20, 2.0,colorCocktailEta,colorCocktailEta);
	DrawGammaSetMarker(cocktailEtapGamma, 20, 2.0,colorCocktailEtaP,colorCocktailEtaP);
	DrawGammaSetMarker(cocktailOmegaGamma, 20, 2.0,colorCocktailOmega,colorCocktailOmega);
	DrawGammaSetMarker(cocktailPhiGamma, 20, 2.0,colorCocktailPhi,colorCocktailPhi);
	DrawGammaSetMarker(cocktailRhoGamma, 20, 2.0,colorCocktailRho0,colorCocktailRho0);
	//DrawGammaSetMarker(cocktailSigmaGamma, 20, 2.0,kGray,kGray);

	cocktailAllGamma->GetXaxis()->SetRangeUser(0,(histoGammaSpecCorrPurity->GetXaxis())->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
	cocktailAllGamma->GetYaxis()->SetTitleSize(0.035);
	cocktailAllGamma->GetYaxis()->SetTitleOffset(1.3);
	cocktailAllGamma->Draw("chist");
	cocktailPi0Gamma->Draw("csamehist");
	cocktailEtaGamma->Draw("csamehist");
	cocktailEtapGamma->Draw("csamehist");
	cocktailOmegaGamma->Draw("csamehist");
	cocktailPhiGamma->Draw("csamehist");
	cocktailRhoGamma->Draw("csamehist");
	//cocktailSigmaGamma->Draw("csamehist");
	
	cocktail                                    = new TLatex(0.65,0.9,"all decay #gamma");
	SetStyleTLatex( cocktail, 0.04,4,colorCocktailAllDecay,42);
	cocktail->Draw();
	
	tpi                                         = new TLatex(0.65,0.855,"#pi^{0} #rightarrow #gamma#gamma (e^{+}e^{-}#gamma)");
	SetStyleTLatex( tpi, 0.04,4,colorCocktailPi0,42);
	tpi->Draw();
	
	teta                                        = new TLatex(0.65,0.81,"#eta #rightarrow #gamma#gamma (#pi^{+}#pi^{-}#gamma,e^{+}e^{-}#gamma,#pi^{0}#gamma#gamma)");
	SetStyleTLatex( teta, 0.04,4,colorCocktailEta,42);
	teta->Draw();
	
	tomega                                      = new TLatex(0.65,0.765,"#omega #rightarrow #pi^{0}#gamma (#eta#gamma)");
	SetStyleTLatex( tomega, 0.04,4,colorCocktailOmega,42);
	tomega->Draw();
	
	tetaprime                                   = new TLatex(0.65,0.725,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)");
	SetStyleTLatex( tetaprime, 0.04,4,colorCocktailEtaP,42);
	tetaprime->Draw();
	
	tphi                                        = new TLatex(0.65,0.68,"#phi #rightarrow #eta#gamma (#pi^{0}#gamma, #omega#gamma)");
	SetStyleTLatex( tphi, 0.04,4,colorCocktailPhi,42);
	tphi->Draw();
	
	trho                                        = new TLatex(0.65,0.635,"#rho #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)");
	SetStyleTLatex( trho, 0.04,4,colorCocktailRho0,42);
	trho->Draw();

	cocktailCanvasSpec->Print(Form("%s/%s_Cocktail_Spectra_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    
    if (cocktailFileMtScaledEta){
		cocktailAllGammaMtScaledEta             = (TH1D*) cocktailDirMtScaledEta->Get("ptg");
		cocktailPi0GammaMtScaledEta             = (TH1D*) cocktailDirMtScaledEta->Get("ptgPion");
		cocktailEtaGammaMtScaledEta             = (TH1D*) cocktailDirMtScaledEta->Get("ptgEta");
		cocktailEtapGammaMtScaledEta            = (TH1D*) cocktailDirMtScaledEta->Get("ptgEtaprime");
		cocktailOmegaGammaMtScaledEta           = (TH1D*) cocktailDirMtScaledEta->Get("ptgOmega");
		cocktailPhiGammaMtScaledEta             = (TH1D*) cocktailDirMtScaledEta->Get("ptgPhi");
		cocktailRhoGammaMtScaledEta             = (TH1D*) cocktailDirMtScaledEta->Get("ptgRho");
		//cocktailSigmaGammaMtScaledEta         = (TH1D*) cocktailDirMtScaledEta->Get("ptgSigma");
		cocktailPi0MtScaledEta                  = (TH1D*) cocktailDirMtScaledEta->Get("ptPion");
		cocktailEtaMtScaledEta                  = (TH1D*) cocktailDirMtScaledEta->Get("ptEta");

		cocktailCanvasSpec->SetLogy();
		SetHistogramm(cocktailAllGammaMtScaledEta,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)",5e-13,1e3);
		DrawGammaSetMarker(cocktailAllGammaMtScaledEta, 20, 2.0,colorCocktailAllDecay,colorCocktailAllDecay);
		SetHistogramm(cocktailPi0GammaMtScaledEta,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
		DrawGammaSetMarker(cocktailPi0GammaMtScaledEta, 20, 2.0,colorCocktailPi0,colorCocktailPi0);
		SetHistogramm(cocktailEtaGammaMtScaledEta,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
		DrawGammaSetMarker(cocktailEtaGammaMtScaledEta, 20, 2.0,colorCocktailEta,colorCocktailEta);
		SetHistogramm(cocktailEtapGammaMtScaledEta,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
		DrawGammaSetMarker(cocktailEtapGammaMtScaledEta, 20, 2.0,colorCocktailEtaP,colorCocktailEtaP);
		SetHistogramm(cocktailOmegaGammaMtScaledEta,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
		DrawGammaSetMarker(cocktailOmegaGammaMtScaledEta, 20, 2.0,colorCocktailOmega,colorCocktailOmega);
		SetHistogramm(cocktailPhiGammaMtScaledEta,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
		DrawGammaSetMarker(cocktailPhiGammaMtScaledEta, 20, 2.0,colorCocktailPhi,colorCocktailPhi);
		SetHistogramm(cocktailRhoGammaMtScaledEta,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
		DrawGammaSetMarker(cocktailRhoGammaMtScaledEta, 20, 2.0,colorCocktailRho0,colorCocktailRho0);
		// SetHistogramm(cocktailSigmaGammaMtScaledEta,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
		// DrawGammaSetMarker(cocktailSigmaGammaMtScaledEta, 20, 2.0,kGray,kGray);

		cocktailAllGammaMtScaledEta->GetXaxis()->SetRangeUser(0,(histoGammaSpecCorrPurity->GetXaxis())->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
		cocktailAllGammaMtScaledEta->GetYaxis()->SetTitleSize(0.035);
		cocktailAllGammaMtScaledEta->GetYaxis()->SetTitleOffset(1.3);
		cocktailAllGammaMtScaledEta->Draw("chist");
		cocktailPi0GammaMtScaledEta->Draw("csamehist");
		cocktailEtaGammaMtScaledEta->Draw("csamehist");
		cocktailEtapGammaMtScaledEta->Draw("csamehist");
		cocktailOmegaGammaMtScaledEta->Draw("csamehist");
		cocktailPhiGammaMtScaledEta->Draw("csamehist");
		cocktailRhoGammaMtScaledEta->Draw("csamehist");
		//cocktailSigmaGammaMtScaledEta->Draw("csamehist");
		
		cocktail->Draw();
		tpi->Draw();
		teta->Draw();
		tomega->Draw();
		tetaprime->Draw();
		tphi->Draw();
		trho->Draw();

		cocktailCanvasSpec->Print(Form("%s/%s_Cocktail_Spectra_mtScaledEta_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));
	}
    
    cout << __LINE__ << endl;

	TH1D* cocktailChargedPionsCorrectBinning    = (TH1D*) histoCorrectedPi0Yield->Clone("cocktailChargedPionsCorrectBinning");
	cout << "here" <<endl;
	if (cocktailFileChargedPions){
		cocktailAllGammaChargedPions            = (TH1D*) cocktailDirChargedPions->Get("ptg");
		cocktailPi0GammaChargedPions            = (TH1D*) cocktailDirChargedPions->Get("ptgPi0");
		cocktailEtaGammaChargedPions            = (TH1D*) cocktailDirChargedPions->Get("ptgEta");
		cocktailEtapGammaChargedPions           = (TH1D*) cocktailDirChargedPions->Get("ptgEtaprime");
		cocktailOmegaGammaChargedPions          = (TH1D*) cocktailDirChargedPions->Get("ptgOmega");
		cocktailPhiGammaChargedPions            = (TH1D*) cocktailDirChargedPions->Get("ptgPhi");
		cocktailRhoGammaChargedPions            = (TH1D*) cocktailDirChargedPions->Get("ptgRho");
		//cocktailSigmaGammaChargedPions        = (TH1D*) cocktailDirChargedPions->Get("ptgSigma");
		cocktailPi0ChargedPions                 = (TH1D*) cocktailDirChargedPions->Get("ptPi0");
		cocktailEtaChargedPions                 = (TH1D*) cocktailDirChargedPions->Get("ptEta");

		cocktailCanvasSpec->SetLogy();
        
		SetHistogramm(cocktailAllGammaChargedPions,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)",5e-13,1e3);
        SetHistogramm(cocktailPi0GammaChargedPions,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
        SetHistogramm(cocktailEtaGammaChargedPions,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
        SetHistogramm(cocktailEtapGammaChargedPions,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
        SetHistogramm(cocktailOmegaGammaChargedPions,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
        SetHistogramm(cocktailPhiGammaChargedPions,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
        SetHistogramm(cocktailRhoGammaChargedPions,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
        //SetHistogramm(cocktailSigmaGammaChargedPions,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");

        DrawGammaSetMarker(cocktailAllGammaChargedPions, 20, 2.0,colorCocktailAllDecay,colorCocktailAllDecay);
		DrawGammaSetMarker(cocktailPi0GammaChargedPions, 20, 2.0,colorCocktailPi0,colorCocktailPi0);
		DrawGammaSetMarker(cocktailEtaGammaChargedPions, 20, 2.0,colorCocktailEta,colorCocktailEta);
		DrawGammaSetMarker(cocktailEtapGammaChargedPions, 20, 2.0,colorCocktailEtaP,colorCocktailEtaP);
		DrawGammaSetMarker(cocktailOmegaGammaChargedPions, 20, 2.0,colorCocktailOmega,colorCocktailOmega);
		DrawGammaSetMarker(cocktailPhiGammaChargedPions, 20, 2.0,colorCocktailPhi,colorCocktailPhi);
		DrawGammaSetMarker(cocktailRhoGammaChargedPions, 20, 2.0,colorCocktailRho0,colorCocktailRho0);
		//DrawGammaSetMarker(cocktailSigmaGammaChargedPions, 20, 2.0,kGray,kGray);

		cocktailAllGammaChargedPions->GetXaxis()->SetRangeUser(0,(histoGammaSpecCorrPurity->GetXaxis())->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
		cocktailAllGammaChargedPions->GetYaxis()->SetTitleSize(0.035);
		cocktailAllGammaChargedPions->GetYaxis()->SetTitleOffset(1.3);
		cocktailAllGammaChargedPions->Draw("chist");
		cocktailPi0GammaChargedPions->Draw("csamehist");
		cocktailEtaGammaChargedPions->Draw("csamehist");
		cocktailEtapGammaChargedPions->Draw("csamehist");
		cocktailOmegaGammaChargedPions->Draw("csamehist");
		cocktailPhiGammaChargedPions->Draw("csamehist");
		cocktailRhoGammaChargedPions->Draw("csamehist");
		//cocktailSigmaGammaChargedPions->Draw("csamehist");
		
		cocktail->Draw();
		tpi->Draw();
		teta->Draw();
		tomega->Draw();
		tetaprime->Draw();
		tphi->Draw();
		trho->Draw();

		cocktailCanvasSpec->Print(Form("%s/%s_Cocktail_Spectra_chargedPions_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));
		
		cout << "here" <<endl;
		cocktailCanvasMesonSpec->cd();
		
		SetHistogramm(cocktailEtaChargedPions,"p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)");
		DrawGammaSetMarker(cocktailEtaChargedPions, 20, 2.0,kMagenta,kMagenta);

		DrawAutoGammaMesonHistos( cocktailPi0ChargedPions, 
								"", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c)", 
								kTRUE, 2.,5e-10, kFALSE,
								kFALSE, -0.004, 0.020, 
								kTRUE, 0.,(histoCorrectedPi0Yield->GetXaxis())->GetBinUpEdge(histoCorrectedPi0Yield->GetNbinsX()));
		DrawGammaSetMarker(cocktailPi0ChargedPions, 20, 2.0,kRed,kRed);
		
		cocktailPi0ChargedPions->Draw("chist");
		TF1* fitChargedPionsCocktail        = FitObject("oHag","chargedPionFit","Pi0",NULL);
		cout << "here" << endl;
        //fitChargedPions->SetParLimits(0, 0.1, 100);
        //fitChargedPions->SetParLimits(1, 0.1, 100);
        //fitChargedPions->SetRange(0.5,15.);
        //fitChargedPions->SetParLimits(1, par1min3, par1max3);
        //fitChargedPions->SetParLimits(3, 4, 10);
        //fitChargedPions->SetParLimits(3, par3min3, par3max3);
        //fitChargedPions->SetParLimits(4, par4min3, par4max3);

        //fitChargedPions->FixParameter(1,oHagparameters[1]);
        //fitChargedPions->FixParameter(2,oHagparameters[2]);
        //fitChargedPions->FixParameter(3,oHagparameters[3]);
        //fitChargedPions->FixParameter(4,oHagparameters[4]);
		cocktailPi0ChargedPions->Fit(fitChargedPionsCocktail,"INRME+","",0.4,15.);
		DrawGammaSetMarkerTF1(fitChargedPionsCocktail, 1, 1.0, kBlue-2);
		fitChargedPionsCocktail->SetRange(0.1,15.);
        //fitChargedPionsCocktail->Draw("same");
        //DrawGammaSetMarker(cocktailPi0ChargedPions, 24, 1., kRed+2, kRed+2);
        //cocktailPi0ChargedPions->DrawCopy("same");
		
		TFitResultPtr FitresultCocktail;
		TString resultCocktail                      = "";
		Int_t iterations                            = 11;
		for(Int_t i = 0; i<iterations;i++){
			cout<<"Fit Iteration "<<i<<endl;
			FitresultCocktail                       = cocktailPi0ChargedPions->Fit(fitChargedPionsCocktail, "SIQNRME+", "", 0.4, 15);
			resultCocktail                          = gMinuit->fCstatu;
			if(!resultC.CompareTo("SUCCESSFUL")) break;
			else iterations++;
			if(iterations>600) break;
		}
		
		cout << "here" <<endl;

		fitChargedPionsCocktail->Draw("same");
		cocktailEtaChargedPions->DrawCopy("csamehist");

		// explanation of points
		DrawGammaSetMarker(histoCorrectedPi0Yield, 24, 1., kGreen+2, kGreen+2);
		histoCorrectedPi0Yield->Draw("e1,same");		
		DrawGammaSetMarker(histoChargedPionSpecStatErr, 4, 1., 1, 1);
		histoChargedPionSpecStatErr->Draw("e1,same");

		for( Int_t t = 1; t < histoChargedPionSpecStatErr->GetNbinsX()+1; t++){
			cout << t << "\t" << histoChargedPionSpecStatErr->GetXaxis()->GetBinLowEdge(t) << "\t" << histoChargedPionSpecStatErr->GetXaxis()->GetBinUpEdge(t) << endl;
		}	
		for(Int_t bin = 1; bin<cocktailChargedPionsCorrectBinning->GetNbinsX()+1; bin++){

			Double_t ptStart                        = cocktailChargedPionsCorrectBinning->GetBinLowEdge(bin);
			Double_t ptEnd                          = cocktailChargedPionsCorrectBinning->GetBinLowEdge(bin+1);
			Double_t binWidth                       = cocktailChargedPionsCorrectBinning->GetBinWidth(bin);

			cocktailChargedPionsCorrectBinning->SetBinContent(bin,fitChargedPionsCocktail->Eval(cocktailChargedPionsCorrectBinning->GetBinCenter(bin)));
			cocktailChargedPionsCorrectBinning->SetBinError(bin,(fitChargedPionsCocktail->IntegralError(ptStart, ptEnd, FitresultCocktail->GetParams(), FitresultCocktail->GetCovarianceMatrix().GetMatrixArray() ) / binWidth));
		}

		cout << "here" <<endl;	

		TLatex* textChargedPion                     = new TLatex(3.0,0.02,"#pi^{#pm} from ALICE Data");
		TLatex* textNeutralPion                     = new TLatex(4.0,0.002,"#pi^{0} from ALICE Data");
		TLatex* textCocktailChargedPion             = new TLatex(1.,1.5,"#pi^{#pm}");
		TLatex* textCocktailScaledEtaChargedPions   = new TLatex(.7,0.005,"#eta");

		textNeutralPion->SetTextColor(kGreen+2);
		textNeutralPion->SetTextSize(0.035);
		textNeutralPion->SetTextFont(42);
		textNeutralPion->Draw();

		textChargedPion->SetTextColor(kBlack);
		textChargedPion->SetTextSize(0.035);
		textChargedPion->SetTextFont(42);
		textChargedPion->Draw();
		textCocktailChargedPion->SetTextColor(kRed);
		textCocktailChargedPion->SetTextSize(0.04);
		textCocktailChargedPion->SetTextFont(42);
		textCocktailChargedPion->Draw();
		textCocktailScaledEtaChargedPions->SetTextColor(kMagenta);
		textCocktailScaledEtaChargedPions->SetTextSize(0.04);
		textCocktailScaledEtaChargedPions->SetTextFont(42);
		textCocktailScaledEtaChargedPions->Draw();

		TLegend* leg_MesonSpectraChargedPions       = GetAndSetLegend(0.7,0.75,2);
		leg_MesonSpectraChargedPions->AddEntry(cocktailPi0ChargedPions,"Cocktail #pi^{0}");
		if(cocktailEta) leg_MesonSpectraChargedPions->AddEntry(cocktailEtaChargedPions,"Cocktail #eta","l");
		leg_MesonSpectraChargedPions->Draw();

		cocktailCanvasMesonSpec->Print(Form("%s/%s_Cocktail_MesonSpectraChargedPion_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));
	}

    cout << __LINE__ << endl;
	
	//******************************************************************************************************************
	//**********************************Plotting Cocktail Ratios *******************************************************
	//******************************************************************************************************************

	TCanvas *cocktailCanvasRatio                    = GetAndSetCanvas("cocktailCanvasRatio");
	cocktailCanvasRatio->SetLogy();

	cocktailAllGammaPi0                             = (TH1D*) cocktailDir->Get("sumgammapi0");
	cocktailPi0GammaPi0                             = (TH1D*) cocktailDir->Get("pi0gammapi0");
	cocktailEtaGammaPi0                             = (TH1D*) cocktailDir->Get("etagammapi0");
	cocktailOmegaGammaPi0                           = (TH1D*) cocktailDir->Get("omegagammapi0");
	cocktailEtapGammaPi0                            = (TH1D*) cocktailDir->Get("etapgammapi0");
	cocktailPhiGammaPi0                             = (TH1D*) cocktailDir->Get("phigammapi0");
	cocktailRhoGammaPi0                             = (TH1D*) cocktailDir->Get("rhogammapi0");
	//cocktailSigmaGammaPi0                         = (TH1D*) cocktailDir->Get("sigmagammapi0");

	SetHistogramm(cocktailAllGammaPi0,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
	SetHistogramm(cocktailPi0GammaPi0,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
    SetHistogramm(cocktailEtaGammaPi0,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
    SetHistogramm(cocktailEtapGammaPi0,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
    SetHistogramm(cocktailOmegaGammaPi0,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
    SetHistogramm(cocktailPhiGammaPi0,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
    SetHistogramm(cocktailRhoGammaPi0,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
    //SetHistogramm(cocktailSigmaGammaPi0,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");

    DrawGammaSetMarker(cocktailAllGammaPi0, 20, 2.0,colorCocktailAllDecay,colorCocktailAllDecay);
	DrawGammaSetMarker(cocktailPi0GammaPi0, 20, 2.0,colorCocktailPi0,colorCocktailPi0);
	DrawGammaSetMarker(cocktailEtaGammaPi0, 20, 2.0,colorCocktailEta,colorCocktailEta);
	DrawGammaSetMarker(cocktailEtapGammaPi0, 20, 2.0,colorCocktailEtaP,colorCocktailEtaP);
	DrawGammaSetMarker(cocktailOmegaGammaPi0, 20, 2.0,colorCocktailOmega,colorCocktailOmega);
	DrawGammaSetMarker(cocktailPhiGammaPi0, 20, 2.0,colorCocktailPhi,colorCocktailPhi);
	DrawGammaSetMarker(cocktailRhoGammaPi0, 20, 2.0,colorCocktailRho0,colorCocktailRho0);
	//DrawGammaSetMarker(cocktailSigmaGammaPi0, 20, 2.0,kGray,kGray);

	cocktailAllGammaPi0->GetXaxis()->SetRangeUser(0.,(histoGammaSpecCorrPurity->GetXaxis())->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
	cocktailAllGammaPi0->GetYaxis()->SetRangeUser(1e-5,1e2);
	cocktailAllGammaPi0->Draw("chist");
	cocktailPi0GammaPi0->Draw("csamehist");
	cocktailEtaGammaPi0->Draw("csamehist");
	cocktailEtapGammaPi0->Draw("csamehist");
	cocktailOmegaGammaPi0->Draw("csamehist");
	cocktailPhiGammaPi0->Draw("csamehist");
	cocktailRhoGammaPi0->Draw("csamehist");
	//cocktailSigmaGammaPi0->Draw("csamehist");
	
	cout << "here" <<endl;
	cocktail                                        = new TLatex(0.45,0.9,"all decay #gamma");
	SetStyleTLatex( cocktail, 0.04,4,colorCocktailAllDecay,42);
	cocktail->Draw();
	tpi                                             = new TLatex(0.45,0.855,"#pi^{0} #rightarrow #gamma#gamma (e^{+}e^{-}#gamma)");
	SetStyleTLatex( tpi, 0.04,4,colorCocktailPi0,42);
	tpi->Draw();
	teta                                            = new TLatex(0.45,0.81,"#eta #rightarrow #gamma#gamma (#pi^{+}#pi^{-}#gamma,e^{+}e^{-}#gamma,#pi^{0}#gamma#gamma)");
	SetStyleTLatex( teta, 0.04,4,colorCocktailEta,42);
	teta->Draw();
	tomega                                          = new TLatex(0.45,0.765,"#omega #rightarrow #pi^{0}#gamma (#eta#gamma)");
	SetStyleTLatex( tomega, 0.04,4,colorCocktailOmega,42);
	tomega->Draw();
	tetaprime                                       = new TLatex(0.75,0.855,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)");
	SetStyleTLatex( tetaprime, 0.04,4,colorCocktailEtaP,42);
	tetaprime->Draw();
	tphi                                            = new TLatex(0.75,0.81,"#phi #rightarrow #eta#gamma (#pi^{0}#gamma, #omega#gamma)");
	SetStyleTLatex( tphi, 0.04,4,colorCocktailPhi,42);
	tphi->Draw();
	trho                                            = new TLatex(0.75,0.765,"#rho #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)");
	SetStyleTLatex( trho, 0.04,4,colorCocktailRho0,42);
	trho->Draw();
	
	cocktailCanvasRatio->Print(Form("%s/%s_Cocktail_Ratios_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    
	if (cocktailFileChargedPions){
		cocktailCanvasRatio->SetLogy();

		cocktailAllGammaPi0ChargedPions             = (TH1D*) cocktailDirChargedPions->Get("sumgammapi0");
		cocktailPi0GammaPi0ChargedPions             = (TH1D*) cocktailDirChargedPions->Get("pi0gammapi0");
		cocktailEtaGammaPi0ChargedPions             = (TH1D*) cocktailDirChargedPions->Get("etagammapi0");
		cocktailOmegaGammaPi0ChargedPions           = (TH1D*) cocktailDirChargedPions->Get("omegagammapi0");
		cocktailEtapGammaPi0ChargedPions            = (TH1D*) cocktailDirChargedPions->Get("etapgammapi0");
		cocktailPhiGammaPi0ChargedPions             = (TH1D*) cocktailDirChargedPions->Get("phigammapi0");
		cocktailRhoGammaPi0ChargedPions             = (TH1D*) cocktailDirChargedPions->Get("rhogammapi0");
		//cocktailSigmaGammaPi0ChargedPions         = (TH1D*) cocktailDir->Get("sigmagammapi0");
        
		SetHistogramm(cocktailAllGammaPi0ChargedPions,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailPi0GammaPi0ChargedPions,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailEtaGammaPi0ChargedPions,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailEtapGammaPi0ChargedPions,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailOmegaGammaPi0ChargedPions,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailPhiGammaPi0ChargedPions,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailRhoGammaPi0ChargedPions,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        //SetHistogramm(cocktailSigmaGammaPi0ChargedPions,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");

        DrawGammaSetMarker(cocktailAllGammaPi0ChargedPions, 20, 2.0,colorCocktailAllDecay,colorCocktailAllDecay);
		DrawGammaSetMarker(cocktailPi0GammaPi0ChargedPions, 20, 2.0,colorCocktailPi0,colorCocktailPi0);
		DrawGammaSetMarker(cocktailEtaGammaPi0ChargedPions, 20, 2.0,colorCocktailEta,colorCocktailEta);
		DrawGammaSetMarker(cocktailEtapGammaPi0ChargedPions, 20, 2.0,colorCocktailEtaP,colorCocktailEtaP);
		DrawGammaSetMarker(cocktailOmegaGammaPi0ChargedPions, 20, 2.0,colorCocktailOmega,colorCocktailOmega);
		DrawGammaSetMarker(cocktailPhiGammaPi0ChargedPions, 20, 2.0,colorCocktailPhi,colorCocktailPhi);
		DrawGammaSetMarker(cocktailRhoGammaPi0ChargedPions, 20, 2.0,colorCocktailRho0,colorCocktailRho0);
		//DrawGammaSetMarker(cocktailSigmaGammaPi0ChargedPions, 20, 2.0,kGray,kGray);

		cocktailAllGammaPi0ChargedPions->GetXaxis()->SetRangeUser(0.,(histoGammaSpecCorrPurity->GetXaxis())->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
		cocktailAllGammaPi0ChargedPions->GetYaxis()->SetRangeUser(1e-5,1e2);
		cocktailAllGammaPi0ChargedPions->Draw("chist");
		cocktailPi0GammaPi0ChargedPions->Draw("csamehist");
		cocktailEtaGammaPi0ChargedPions->Draw("csamehist");
		cocktailEtapGammaPi0ChargedPions->Draw("csamehist");
		cocktailOmegaGammaPi0ChargedPions->Draw("csamehist");
		cocktailPhiGammaPi0ChargedPions->Draw("csamehist");
		cocktailRhoGammaPi0ChargedPions->Draw("csamehist");
		//cocktailSigmaGammaPi0ChargedPions->Draw("csamehist");
		
		cocktail->Draw();
		tpi->Draw();
		teta->Draw();
		tomega->Draw();
		tetaprime->Draw();
		tphi->Draw();
		trho->Draw();
		
		cocktailCanvasRatio->Print(Form("%s/%s_Cocktail_Ratios_ChargedPions_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));
	}
    
    cout << __LINE__ << endl;

	if (cocktailFileMtScaledEta){
		cocktailCanvasRatio->SetLogy();

		cocktailAllGammaPi0MtScaledEta              = (TH1D*) cocktailDirMtScaledEta->Get("sumgammapi0");
		cocktailPi0GammaPi0MtScaledEta              = (TH1D*) cocktailDirMtScaledEta->Get("pi0gammapi0");
		cocktailEtaGammaPi0MtScaledEta              = (TH1D*) cocktailDirMtScaledEta->Get("etagammapi0");
		cocktailOmegaGammaPi0MtScaledEta            = (TH1D*) cocktailDirMtScaledEta->Get("omegagammapi0");
		cocktailEtapGammaPi0MtScaledEta             = (TH1D*) cocktailDirMtScaledEta->Get("etapgammapi0");
		cocktailPhiGammaPi0MtScaledEta              = (TH1D*) cocktailDirMtScaledEta->Get("phigammapi0");
		cocktailRhoGammaPi0MtScaledEta              = (TH1D*) cocktailDirMtScaledEta->Get("rhogammapi0");
		//cocktailSigmaGammaPi0MtScaledEta          = (TH1D*) cocktailDir->Get("sigmagammapi0");

		SetHistogramm(cocktailAllGammaPi0MtScaledEta,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailPi0GammaPi0MtScaledEta,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailEtaGammaPi0MtScaledEta,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailEtapGammaPi0MtScaledEta,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailOmegaGammaPi0MtScaledEta,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailPhiGammaPi0MtScaledEta,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        SetHistogramm(cocktailRhoGammaPi0MtScaledEta,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");
        //SetHistogramm(cocktailSigmaGammaPi0MtScaledEta,"p_{T} (GeV/c)","Ratio #gamma / #pi^{0} from different sources");

        DrawGammaSetMarker(cocktailAllGammaPi0MtScaledEta, 20, 2.0,colorCocktailAllDecay,colorCocktailAllDecay);
		DrawGammaSetMarker(cocktailPi0GammaPi0MtScaledEta, 20, 2.0,colorCocktailPi0,colorCocktailPi0);
		DrawGammaSetMarker(cocktailEtaGammaPi0MtScaledEta, 20, 2.0,colorCocktailEta,colorCocktailEta);
		DrawGammaSetMarker(cocktailEtapGammaPi0MtScaledEta, 20, 2.0,colorCocktailEtaP,colorCocktailEtaP);
		DrawGammaSetMarker(cocktailOmegaGammaPi0MtScaledEta, 20, 2.0,colorCocktailOmega,colorCocktailOmega);
		DrawGammaSetMarker(cocktailPhiGammaPi0MtScaledEta, 20, 2.0,colorCocktailPhi,colorCocktailPhi);
		DrawGammaSetMarker(cocktailRhoGammaPi0MtScaledEta, 20, 2.0,colorCocktailRho0,colorCocktailRho0);
		//DrawGammaSetMarker(cocktailSigmaGammaPi0MtScaledEta, 20, 2.0,kGray,kGray);

		cocktailAllGammaPi0MtScaledEta->GetXaxis()->SetRangeUser(0.,(histoGammaSpecCorrPurity->GetXaxis())->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
		cocktailAllGammaPi0MtScaledEta->GetYaxis()->SetRangeUser(1e-5,1e2);
		cocktailAllGammaPi0MtScaledEta->Draw("chist");
		cocktailPi0GammaPi0MtScaledEta->Draw("csamehist");
		cocktailEtaGammaPi0MtScaledEta->Draw("csamehist");
		cocktailEtapGammaPi0MtScaledEta->Draw("csamehist");
		cocktailOmegaGammaPi0MtScaledEta->Draw("csamehist");
		cocktailPhiGammaPi0MtScaledEta->Draw("csamehist");
		cocktailRhoGammaPi0MtScaledEta->Draw("csamehist");
		//cocktailSigmaGammaPi0MtScaledEta->Draw("csamehist");
		
		cocktail->Draw();
		tpi->Draw();
		teta->Draw();
		tomega->Draw();
		tetaprime->Draw();
		tphi->Draw();
		trho->Draw();
		
		cocktailCanvasRatio->Print(Form("%s/%s_Cocktail_Ratios_MtScaledEta_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));
	}

    cout << __LINE__ << endl;
	
	//    return;
	TH1D *cocktailAllGammaPi0original               = (TH1D*)cocktailAllGammaPi0->Clone("cocktailAllGammaPi0original");

	//**********************************************************************************
	//******************** NLO Direct Photon Ratio **************************************
	//**********************************************************************************
	if (doNLOComparison){
	
		TF1* QGPin7Tev                              = new TF1("QGPin7Tev","1+[0]*exp(-x/[1])+[2]*exp(-x/[3])",0,100);
		QGPin7Tev->SetParameters(4.24707e+01, 2.85766e+00, 1.63353e+02, 4.53670e-01);

		TH1D *cocktailAllGammaNLO                   = (TH1D*) cocktailAllGamma->Clone("cocktailAllGammaNLO");
		graphNLOCalcMuHalfcopyA                     = (TGraphErrors*) graphNLOCalcMuHalf->Clone("graphNLOCalcMuHalfcopyA");
		graphNLOCalcMuHalfcopyB                     = (TGraphErrors*) graphNLOCalcMuOne->Clone("graphNLOCalcMuHalfcopyB");
		graphNLOCalcMuHalfcopyC                     = (TGraphErrors*) graphNLOCalcMuTwo->Clone("graphNLOCalcMuHalfcopyC");

		xHalf                                       = graphNLOCalcMuHalfcopyA->GetX();
		yHalf                                       = graphNLOCalcMuHalfcopyA->GetY();
		eyHalf                                      = graphNLOCalcMuHalfcopyA->GetEY();

		xOne                                        = graphNLOCalcMuHalfcopyB->GetX();
		yOne                                        = graphNLOCalcMuHalfcopyB->GetY();
		eyOne                                       = graphNLOCalcMuHalfcopyB->GetEY();

		xTwo                                        = graphNLOCalcMuHalfcopyC->GetX();
		yTwo                                        = graphNLOCalcMuHalfcopyC->GetY();
		eyTwo=                                      graphNLOCalcMuHalfcopyC->GetEY();

		TString cocktailFit                         = "xqcd";
		//if(fcocktailFunc.Contains("QCD") || fcocktailFunc.Contains("qcd"))
		//    cocktailFit = "qcd";
		//if(fcocktailFunc.Contains("oHag") || fcocktailFunc.Contains("ohag"))
		//    cocktailFit = "oHag";

		TF1 *cocktailFitAllGammaForNLO              = (TF1*) FitObject(cocktailFit,"cocktailFit","Pi0",cocktailAllGammaNLO,2.0,16);

		TCanvas *canvasCocktailFitForNLO            = GetAndSetCanvas("canvasCocktailFitForNLO");
		canvasCocktailFitForNLO->SetLogy();
		cocktailAllGammaNLO->DrawCopy();
		cocktailFitAllGammaForNLO->Draw("same");

		canvasCocktailFitForNLO->Print(Form("%s/%scocktailAllGammaNLO_cocktailFitForNLO_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

        cout << __LINE__ << endl;

		for(Int_t bin = 0; bin < graphNLOCalcMuHalf->GetN(); bin++){// Scale with N_collisions
			yHalf[bin]                              = (1 +( yHalf[bin] / (cocktailFitAllGammaForNLO->Eval(xHalf[bin]))));
			yOne[bin]                               = (1 +( yOne[bin] / (cocktailFitAllGammaForNLO->Eval(xOne[bin]))));
			yTwo[bin]                               = (1 +( yTwo[bin] / (cocktailFitAllGammaForNLO->Eval(xTwo[bin]))));
		}

		half                                        = new TGraphErrors(graphNLOCalcMuHalf->GetN(),xHalf,yHalf,eyHalf);
		one                                         = new TGraphErrors(graphNLOCalcMuOne->GetN(),xOne,yOne,eyOne);
		two                                         = new TGraphErrors(graphNLOCalcMuTwo->GetN(),xTwo,yTwo,eyTwo);
	}

    cout << __LINE__ << endl;

	//**********************************************************************************
	//******************** Double Ratios ***********************************************
	//**********************************************************************************

	cocktailAllGammaPi0                             = RebinTH1D(cocktailAllGammaPi0,histoIncRatioPurity);
	TH1D *cocktailAllGammaRebined                   = RebinTH1D(cocktailAllGamma,histoIncRatioPurity);
	TH1D *cocktailPi0Rebined                        = RebinTH1D(cocktailPi0,histoIncRatioPurity);

	//cocktailPi0Gamma                              = RebinTH1D(cocktailPi0Gamma,histoIncRatioPurity);
	cocktailAllGammaPi0->Divide(cocktailAllGammaRebined,cocktailPi0Rebined);
	//cocktailAllGammaPi0->Divide(cocktailAllGammaRebined,cocktailPi0Gamma);
	if (cocktailFileMtScaledEta){
		cocktailAllGammaPi0MtScaledEta              = RebinTH1D(cocktailAllGammaPi0MtScaledEta,histoIncRatioPurity);
		TH1D *cocktailAllGammaRebinedMtScaledEta    = RebinTH1D(cocktailAllGammaMtScaledEta,histoIncRatioPurity);
		TH1D *cocktailPi0RebinedMtScaledEta         = RebinTH1D(cocktailPi0MtScaledEta,histoIncRatioPurity);
		cocktailAllGammaPi0MtScaledEta->Divide(cocktailAllGammaRebinedMtScaledEta,cocktailPi0RebinedMtScaledEta);
	}
	if (cocktailFileChargedPions){
		cocktailAllGammaPi0ChargedPions             = RebinTH1D(cocktailAllGammaPi0ChargedPions,histoIncRatioPurity);
		TH1D *cocktailAllGammaRebinedChargedPions   = RebinTH1D(cocktailAllGammaChargedPions,histoIncRatioPurity);
		TH1D *cocktailPi0RebinedChargedPions        = RebinTH1D(cocktailPi0ChargedPions,histoIncRatioPurity);
		cocktailAllGammaPi0ChargedPions->Divide(cocktailAllGammaRebinedChargedPions,cocktailPi0Rebined);
	}
    
    cout << __LINE__ << endl;
	
	/////////////////////////////
	TCanvas* canvasIncRatioPurityTrueEffAll         = GetAndSetCanvas(" canvasIncRatioPurityTrueEffAll");
	SetHistogramm(histoIncRatioPurityTrueEff, "p_{T} (GeV/c)", "Ratio Inclusive #gamma/#pi^{0}",0,2);
	DrawGammaSetMarker(histoIncRatioPurityTrueEff, 20, 2.0, 4, 4);
	histoIncRatioPurityTrueEff->Draw();
	cocktailAllGammaPi0->Draw("chistsame");
	
	canvasIncRatioPurityTrueEffAll->Print(Form("%s/%s_IncRatioPurity_trueEff_all_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
	///////////////////////////////

    cout << __LINE__ << endl;

	TCanvas* canvasMCDoubleRatioSum                 = GetAndSetCanvas("canvasMCDoubleRatioSum");
	canvasMCDoubleRatioSum->SetLogy(0);
	histoMCDoubleRatioSum                           = (TH1D*)histoMCIncRatio->Clone("MC_DoubleRatio_Sum");
    histoMCesdDoubleRatioSum                        = (TH1D*)histoMCesdIncRatio->Clone("MCesd_DoubleRatio_Sum");
    histoDecayRatioSumGamma                         = RebinTH1D(histoDecayRatioSumGamma,histoMCIncRatio);
	histoMCDoubleRatioSum->Divide(histoMCDoubleRatioSum,histoDecayRatioSumGamma,1,1,"");
	histoMCesdDoubleRatioSum->Divide(histoMCesdDoubleRatioSum,histoDecayRatioSumGamma,1,1,"");

	DrawGammaSetMarker(histoMCDoubleRatioSum, 22, 2.0, 1, 1);
	DrawGammaSetMarker(histoMCesdDoubleRatioSum, 23, 2.0, 2, 2);
    
	SetHistogramm(histoMCDoubleRatioSum,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.5,2.0);
	SetHistogramm(histoMCesdDoubleRatioSum,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.5,2.0);
    
    histoMCesdDoubleRatioSum->Draw("");
	histoMCDoubleRatioSum->Draw("same");

	canvasMCDoubleRatioSum->Print(Form("%s/%s_DoubleRatioMC_Sum_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

	histoDoubleRatioConversionNormalPurity          = (TH1D*) histoIncRatioPurity->Clone("DoubleRatioConversionNormalPurity");
	histoDoubleRatioConversionNormalPurity->Divide(cocktailAllGammaPi0);
	histoDoubleRatioConversionTrueEffPurity         = (TH1D*) histoIncRatioPurityTrueEff->Clone("DoubleRatioConversionTrueEffPurity");
	histoDoubleRatioConversionTrueEffPurity->Divide(cocktailAllGammaPi0);
	histoDoubleRatioConversionTrueEffPurityWide     = (TH1D*) histoIncRatioPurityTrueEffWide->Clone("DoubleRatioConversionTrueEffPurityWide");
	histoDoubleRatioConversionTrueEffPurityWide->Divide(cocktailAllGammaPi0);
	histoDoubleRatioConversionTrueEffPurityNarrow   = (TH1D*) histoIncRatioPurityTrueEffNarrow->Clone("DoubleRatioConversionTrueEffPurityNarrow");
	histoDoubleRatioConversionTrueEffPurityNarrow->Divide(cocktailAllGammaPi0);
	histoDoubleRatioConversionOnlyGamma             = (TH1D*) histoGammaSpecCorrPurity->Clone("DoubleRatioConversionOnlyGamma");
	histoDoubleRatioConversionOnlyGamma->Divide(cocktailAllGammaRebined);

    cout << __LINE__ << endl;

	//if(option.CompareTo("7TeV") == 0){
	//	histoDoubleRatioCombinedPurity                      = (TH1D*) histoIncRatioCombinedPurity->Clone("DoubleRatioCombinedPurity");  // doesn't exist
	//	histoDoubleRatioCombinedPurity->Divide(cocktailAllGammaPi0);
	//}
    
    cout << __LINE__ << endl;

	histoDoubleRatioFitPi0YieldPurity                       = (TH1D*) histoIncRatioFitPurity->Clone("DoubleRatioConversionFitPurity");
	histoDoubleRatioFitPi0YieldPurity->Divide(cocktailAllGammaPi0);
	histoDoubleRatioFitPi0YieldPurityWide                   = (TH1D*) histoIncRatioFitPurityWide->Clone("DoubleRatioConversionFitPurityWide");
	histoDoubleRatioFitPi0YieldPurityWide->Divide(cocktailAllGammaPi0);
	histoDoubleRatioFitPi0YieldPurityNarrow                 = (TH1D*) histoIncRatioFitPurityNarrow->Clone("DoubleRatioConversionFitPurityNarrow");
	histoDoubleRatioFitPi0YieldPurityNarrow->Divide(cocktailAllGammaPi0);
 
    cout << __LINE__ << endl;

	if (cocktailFileMtScaledEta){
		histoDoubleRatioConversionTrueEffPurityMtScaledEta  = (TH1D*) histoIncRatioPurityTrueEff->Clone("DoubleRatioConversionTrueEffPurity");
		histoDoubleRatioConversionTrueEffPurityMtScaledEta->Divide(cocktailAllGammaPi0MtScaledEta);
		SetHistogramm(histoDoubleRatioConversionTrueEffPurityMtScaledEta,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,2.0);
	
		histoDoubleRatioFitPi0YieldPurityMtScaledEta        = (TH1D*) histoIncRatioFitPurity->Clone("DoubleRatioConversionFitPurity");
		histoDoubleRatioFitPi0YieldPurityMtScaledEta->Divide(cocktailAllGammaPi0MtScaledEta);
		SetHistogramm(histoDoubleRatioFitPi0YieldPurityMtScaledEta,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,2.0);
	}
    
    cout << __LINE__ << endl;
    
	TH1D* histoIncRatioPurityTrueEffChargedPions            = NULL;
	TH1D* histoCorrectedChargedPionYield                    = NULL;
	if (cocktailFileChargedPions){
		
		histoCorrectedChargedPionYield                      = (TH1D*)histoCorrectedPi0Yield->Clone("histoCorrectedChargedPionYield");
		if (doChargedPionComp){
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(0.35),(histoChargedPionSpecStatErr->GetBinContent(9)*histoChargedPionSpecStatErr->GetBinCenter(9) *histoChargedPionSpecStatErr->GetBinWidth(9) +histoChargedPionSpecStatErr->GetBinContent(10) *histoChargedPionSpecStatErr->GetBinCenter(10)*histoChargedPionSpecStatErr->GetBinWidth(10))/0.35/0.1);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(0.45),(histoChargedPionSpecStatErr->GetBinContent(11)*histoChargedPionSpecStatErr->GetBinCenter(11) *histoChargedPionSpecStatErr->GetBinWidth(11) +histoChargedPionSpecStatErr->GetBinContent(12)*histoChargedPionSpecStatErr->GetBinCenter(12)*histoChargedPionSpecStatErr->GetBinWidth(12))/0.45/0.1);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(0.55),(histoChargedPionSpecStatErr->GetBinContent(13)*histoChargedPionSpecStatErr->GetBinCenter(13) *histoChargedPionSpecStatErr->GetBinWidth(13) +histoChargedPionSpecStatErr->GetBinContent(14)*histoChargedPionSpecStatErr->GetBinCenter(14)*histoChargedPionSpecStatErr->GetBinWidth(14))/0.55/0.1);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(0.65),(histoChargedPionSpecStatErr->GetBinContent(15)*histoChargedPionSpecStatErr->GetBinCenter(15) *histoChargedPionSpecStatErr->GetBinWidth(15) +histoChargedPionSpecStatErr->GetBinContent(16)*histoChargedPionSpecStatErr->GetBinCenter(16)*histoChargedPionSpecStatErr->GetBinWidth(16))/0.65/0.1);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(0.75),(histoChargedPionSpecStatErr->GetBinContent(17)*histoChargedPionSpecStatErr->GetBinCenter(17) *histoChargedPionSpecStatErr->GetBinWidth(17) +histoChargedPionSpecStatErr->GetBinContent(18)*histoChargedPionSpecStatErr->GetBinCenter(18)*histoChargedPionSpecStatErr->GetBinWidth(18))/0.75/0.1);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(0.9),(histoChargedPionSpecStatErr->GetBinContent(19)*histoChargedPionSpecStatErr->GetBinCenter(19) *histoChargedPionSpecStatErr->GetBinWidth(19) +histoChargedPionSpecStatErr->GetBinContent(20)*histoChargedPionSpecStatErr->GetBinCenter(20)*histoChargedPionSpecStatErr->GetBinWidth(20) +histoChargedPionSpecStatErr->GetBinContent(21) *histoChargedPionSpecStatErr->GetBinCenter(21)*histoChargedPionSpecStatErr->GetBinWidth(21) +histoChargedPionSpecStatErr->GetBinContent(22)*histoChargedPionSpecStatErr->GetBinCenter(22) *histoChargedPionSpecStatErr->GetBinWidth(22))/0.9/0.2);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(1.1),(histoChargedPionSpecStatErr->GetBinContent(23)*histoChargedPionSpecStatErr->GetBinCenter(23) *histoChargedPionSpecStatErr->GetBinWidth(23) +histoChargedPionSpecStatErr->GetBinContent(24)*histoChargedPionSpecStatErr->GetBinCenter(24)*histoChargedPionSpecStatErr->GetBinWidth(24))/1.1/0.2);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(1.3),(histoChargedPionSpecStatErr->GetBinContent(25)*histoChargedPionSpecStatErr->GetBinCenter(25) *histoChargedPionSpecStatErr->GetBinWidth(25) +histoChargedPionSpecStatErr->GetBinContent(26)*histoChargedPionSpecStatErr->GetBinCenter(26)*histoChargedPionSpecStatErr->GetBinWidth(26))/1.3/0.2);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(1.5),(histoChargedPionSpecStatErr->GetBinContent(27)*histoChargedPionSpecStatErr->GetBinCenter(27) *histoChargedPionSpecStatErr->GetBinWidth(27) +histoChargedPionSpecStatErr->GetBinContent(28)*histoChargedPionSpecStatErr->GetBinCenter(28)*histoChargedPionSpecStatErr->GetBinWidth(28))/1.5/0.2);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(1.7),(histoChargedPionSpecStatErr->GetBinContent(29)*histoChargedPionSpecStatErr->GetBinCenter(29) *histoChargedPionSpecStatErr->GetBinWidth(29) +histoChargedPionSpecStatErr->GetBinContent(30)*histoChargedPionSpecStatErr->GetBinCenter(30)*histoChargedPionSpecStatErr->GetBinWidth(30))/1.7/0.2);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(1.9),(histoChargedPionSpecStatErr->GetBinContent(31)*histoChargedPionSpecStatErr->GetBinCenter(31) *histoChargedPionSpecStatErr->GetBinWidth(31) +histoChargedPionSpecStatErr->GetBinContent(32)*histoChargedPionSpecStatErr->GetBinCenter(32)*histoChargedPionSpecStatErr->GetBinWidth(32))/1.9/0.2);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(2.1),(histoChargedPionSpecStatErr->GetBinContent(33)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(2.3),(histoChargedPionSpecStatErr->GetBinContent(34)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(2.5),(histoChargedPionSpecStatErr->GetBinContent(35)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(2.7),(histoChargedPionSpecStatErr->GetBinContent(36)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(2.9),(histoChargedPionSpecStatErr->GetBinContent(37)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(3.1),(histoChargedPionSpecStatErr->GetBinContent(38)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(3.3),(histoChargedPionSpecStatErr->GetBinContent(39)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(3.5),(histoChargedPionSpecStatErr->GetBinContent(40)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(3.7),(histoChargedPionSpecStatErr->GetBinContent(41)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(3.9),(histoChargedPionSpecStatErr->GetBinContent(42)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(4.25),(histoChargedPionSpecStatErr->GetBinContent(43)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(4.75),(histoChargedPionSpecStatErr->GetBinContent(44)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(5.25),(histoChargedPionSpecStatErr->GetBinContent(45)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(5.75),(histoChargedPionSpecStatErr->GetBinContent(46)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(6.5),(histoChargedPionSpecStatErr->GetBinContent(47)*histoChargedPionSpecStatErr->GetBinCenter(47) *histoChargedPionSpecStatErr->GetBinWidth(47) +histoChargedPionSpecStatErr->GetBinContent(48)*histoChargedPionSpecStatErr->GetBinCenter(48)*histoChargedPionSpecStatErr->GetBinWidth(48))/6.5/1.);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(7.5),(histoChargedPionSpecStatErr->GetBinContent(49)));
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(9.0),(histoChargedPionSpecStatErr->GetBinContent(50)*histoChargedPionSpecStatErr->GetBinCenter(50) *histoChargedPionSpecStatErr->GetBinWidth(50) +histoChargedPionSpecStatErr->GetBinContent(51)*histoChargedPionSpecStatErr->GetBinCenter(51)*histoChargedPionSpecStatErr->GetBinWidth(51))/9/2.);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(11.0),(histoChargedPionSpecStatErr->GetBinContent(52)*histoChargedPionSpecStatErr->GetBinCenter(52) *histoChargedPionSpecStatErr->GetBinWidth(52) +histoChargedPionSpecStatErr->GetBinContent(53)*histoChargedPionSpecStatErr->GetBinCenter(53)*histoChargedPionSpecStatErr->GetBinWidth(53))/11./2.);
			histoCorrectedChargedPionYield->SetBinContent(histoCorrectedChargedPionYield->GetXaxis()->FindBin(13.0),(histoChargedPionSpecStatErr->GetBinContent(54)*histoChargedPionSpecStatErr->GetBinCenter(54) *histoChargedPionSpecStatErr->GetBinWidth(54) +histoChargedPionSpecStatErr->GetBinContent(55)*histoChargedPionSpecStatErr->GetBinCenter(55)*histoChargedPionSpecStatErr->GetBinWidth(55))/13./2.);
		}
		
		histoIncRatioPurityTrueEffChargedPions              = (TH1D*) histoGammaSpecCorrPurity->Clone("IncRatioPurity_trueEff");
		histoIncRatioPurityTrueEffChargedPions->Divide(histoIncRatioPurityTrueEffChargedPions,histoCorrectedChargedPionYield,1,1,"");
		histoIncRatioPurityTrueEffChargedPions->Draw();

		histoDoubleRatioConversionTrueEffPurityChargedPions = (TH1D*) histoIncRatioPurityTrueEffChargedPions->Clone("DoubleRatioConversionTrueEffPurity");
		histoDoubleRatioConversionTrueEffPurityChargedPions->Divide(cocktailAllGammaPi0ChargedPions);
		SetHistogramm(histoDoubleRatioConversionTrueEffPurityChargedPions,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,2.0);
	
		histoDoubleRatioFitPi0YieldPurityChargedPions       = (TH1D*) histoIncRatioPurityTrueEffChargedPions->Clone("DoubleRatioConversionFitPurity");
		histoDoubleRatioFitPi0YieldPurityChargedPions->Divide(cocktailAllGammaPi0ChargedPions);
		SetHistogramm(histoDoubleRatioFitPi0YieldPurityChargedPions,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,2.0);
	}
    
    cout << __LINE__ << endl;

	histoDoubleRatioConversionLowFitPurity                  = (TH1D*) histoIncRatioLowFitPurity->Clone("DoubleRatioConversionLowFitPurity");
	histoDoubleRatioConversionLowFitPurity->Divide(cocktailAllGammaPi0);
	histoDoubleRatioConversionHighFitPurity                 = (TH1D*) histoIncRatioHighFitPurity->Clone("DoubleRatioConversionHighFitPurity");
	histoDoubleRatioConversionHighFitPurity->Divide(cocktailAllGammaPi0);

	SetHistogramm(histoDoubleRatioFitPi0YieldPurity,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,2.0);
	SetHistogramm(histoDoubleRatioConversionLowFitPurity,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,3.5);
	SetHistogramm(histoDoubleRatioConversionHighFitPurity,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,3.5);
	SetHistogramm(histoDoubleRatioFitPi0YieldPurityWide,"#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,2.0);
	SetHistogramm(histoDoubleRatioFitPi0YieldPurityNarrow,"#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,2.0);
	SetHistogramm(histoDoubleRatioConversionTrueEffPurityWide,"#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,1.8);
	SetHistogramm(histoDoubleRatioConversionTrueEffPurityNarrow,"#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,1.8);
    
	DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurityWide, 27, 0.8, kMagenta+1,kMagenta+1);
	DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurityNarrow, 27, 0.8, kMagenta+1,kMagenta+1);
	DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurityWide, 20, 0.8, 1, 1);
	DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurityNarrow, 20, 0.8, 1, 1);

    cout << __LINE__ << endl;
    
	TCanvas *canvasConversionFitDoubleRatioSum              = new TCanvas("canvasConversionFitDoubleRatioSum","",0,0,1350,1000);  // gives the page size
	DrawGammaCanvasSettings( canvasConversionFitDoubleRatioSum, 0.09, 0.02, 0.02, 0.09);
	canvasConversionFitDoubleRatioSum->cd();
		
		TH2F * histo2DDoubleRatioPlotting;
		histo2DDoubleRatioPlotting                          = new TH2F("histo2DDoubleRatioPlotting","histo2DDoubleRatioPlotting",1000,0.0,20.,1000,0.85,1.65);
		SetStyleHistoTH2ForGraphs(histo2DDoubleRatioPlotting, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.04,0.04, 0.04,0.04, 1.,1.);
		histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.,histoDoubleRatioConversionTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioConversionTrueEffPurity->GetNbinsX()));
		histo2DDoubleRatioPlotting->GetXaxis()->SetTitleFont(62);
		histo2DDoubleRatioPlotting->GetYaxis()->SetTitleFont(62);
        //histo2DDoubleRatioPlotting->GetYaxis()->CenterTitle(kTRUE);
		histo2DDoubleRatioPlotting->DrawCopy(); 

		DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurity, 20, 2.0, kBlack, kBlack);   
		DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurity, 20, 2.0, kBlue+2, kBlue+2);   
		
		One->Draw("same");

        //histoDoubleRatioConversionTrueEffPurity->DrawCopy("same");
		histoDoubleRatioFitPi0YieldPurity->DrawCopy("same");
		if (cocktailFileMtScaledEta){
			DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurityMtScaledEta, 21, 2.0, kRed-7, kRed-7);   
			DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurityMtScaledEta, 25, 2.0, kRed+2, kRed+2);   
            //histoDoubleRatioConversionTrueEffPurityMtScaledEta->DrawCopy("same");
			histoDoubleRatioFitPi0YieldPurityMtScaledEta->DrawCopy("same");
		}

		if (cocktailFileChargedPions){
			DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurityChargedPions, 21, 2.0, kGreen-7, kGreen-7);   
			DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurityChargedPions, 25, 2.0, kGreen+2, kGreen+2);   
            //histoDoubleRatioConversionTrueEffPurityChargedPions->DrawCopy("same");
			histoDoubleRatioFitPi0YieldPurityChargedPions->DrawCopy("same");
		}
		
		TLegend* legendDoubleConversionFit                  = GetAndSetLegend(0.14,0.7,4,1,"Direct Photon Signal via Conversions");
        //legendDoubleConversionFit->AddEntry(histoDoubleRatioConversionTrueEffPurity,"Measured Direct Photon Signal","p");
		legendDoubleConversionFit->AddEntry(histoDoubleRatioFitPi0YieldPurity,"Measured Direct Photon Signal fitted #pi^{0}, measured #eta","p");
		if (cocktailFileMtScaledEta){
            //legendDoubleConversionFit->AddEntry(histoDoubleRatioConversionTrueEffPurityMtScaledEta,"Measured Direct Photon Signal, mt-scaled #eta","p");
			legendDoubleConversionFit->AddEntry(histoDoubleRatioFitPi0YieldPurityMtScaledEta,"Measured Direct Photon Signal fitted #pi^{0}, mt-scaled #eta","p");
		}
		if (cocktailFileChargedPions){
            //legendDoubleConversionFit->AddEntry(histoDoubleRatioConversionTrueEffPurityChargedPions,"Measured Direct Photon Signal, mt-scaled #eta","p");
			legendDoubleConversionFit->AddEntry(histoDoubleRatioFitPi0YieldPurityChargedPions,"Measured Direct Photon Signal exchanged #pi^{0} by #pi^{#pm},","p");
			legendDoubleConversionFit->AddEntry((TObject*)0,"scaled #eta accordingly","");
		}

		legendDoubleConversionFit->Draw();

	canvasConversionFitDoubleRatioSum->Print(Form("%s/%s_DoubleRatioFit_Comparisons_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    
	if (doNLOComparison){
		// ------------------------------ NLO Calculations -----------------------------
		TGraphErrors *DirectPhotonDoubleNLOhalf             = (TGraphErrors*) half->Clone("doubleRatioNLOhalf");
		TGraphErrors *DirectPhotonDoubleNLOone              = (TGraphErrors*) one->Clone("doubleRatioNLOone");
		TGraphErrors *DirectPhotonDoubleNLOtwo              = (TGraphErrors*) two->Clone("doubleRatioNLOtwo");

		TGraphErrors *DirectPhotonNLOhalf                   = (TGraphErrors*) graphNLOCalcMuHalf->Clone("graphNLOCalcMuHalf");
		TGraphErrors *DirectPhotonNLOone                    = (TGraphErrors*) graphNLOCalcMuOne->Clone("graphNLOCalcMuOne");
		TGraphErrors *DirectPhotonNLOtwo                    = (TGraphErrors*) graphNLOCalcMuTwo->Clone("graphNLOCalcMuTwo");
		
		Double_t *errorup                                   = new Double_t[DirectPhotonDoubleNLOhalf->GetN()];
		Double_t *errorlow                                  = new Double_t[DirectPhotonDoubleNLOtwo->GetN()];

		Double_t *errorSpecup                               = new Double_t[DirectPhotonDoubleNLOhalf->GetN()];
		Double_t *errorSpeclow                              = new Double_t[DirectPhotonDoubleNLOtwo->GetN()];

		yHalf                                               = DirectPhotonDoubleNLOhalf->GetY();
		yOne                                                = DirectPhotonDoubleNLOone->GetY();
		yTwo                                                = DirectPhotonDoubleNLOtwo->GetY();

		Double_t *ySpecHalf                                 = DirectPhotonNLOhalf->GetY();
		Double_t *ySpecOne                                  = DirectPhotonNLOone->GetY();
		Double_t *ySpecTwo                                  = DirectPhotonNLOtwo->GetY();

		for(Int_t i = 0;i<DirectPhotonDoubleNLOhalf->GetN(); i++){
			errorup[i]                                      = yHalf[i]-yOne[i];
			errorlow[i]                                     = -yTwo[i]+yOne[i];
			errorSpecup[i]                                  = ySpecHalf[i]-ySpecOne[i];
			errorSpeclow[i]                                 = -ySpecTwo[i]+ySpecOne[i];
		}

		TGraphAsymmErrors *NLODoubleRatio                   = new TGraphAsymmErrors( DirectPhotonDoubleNLOhalf->GetN(), DirectPhotonDoubleNLOone->GetX(), DirectPhotonDoubleNLOone->GetY(), DirectPhotonDoubleNLOone->GetEX(), DirectPhotonDoubleNLOone->GetEX(), errorlow,errorup );
		TGraphAsymmErrors *NLO                              = new TGraphAsymmErrors( DirectPhotonDoubleNLOhalf->GetN(), DirectPhotonNLOone->GetX(), DirectPhotonNLOone->GetY(), DirectPhotonNLOone->GetEX(), DirectPhotonNLOone->GetEX(),errorSpeclow ,errorSpecup );
        
		NLODoubleRatio->SetLineColor(kAzure);
		NLODoubleRatio->SetFillColor(kAzure);
		NLODoubleRatio->SetLineWidth(3.0);
		NLODoubleRatio->SetMarkerSize(0);
		
		NLODoubleRatio->Print();
		NLO->SetLineColor(kAzure);
		NLO->SetFillColor(kAzure);
		NLO->SetLineWidth(3.0);
		NLO->SetMarkerSize(0);
		NLO->RemovePoint(0);
        
		TF1 *NLOdoubleRatioFit                              = new TF1("NLOdoubleRatioFit","(x<=2)*(1.0+[0]*x)+(x>2)*([1]+[2]*x+[3]*x*x)",0.,4.0);
		NLODoubleRatio->Fit(NLOdoubleRatioFit,"NRME+","",0,4.0);
		NLOdoubleRatioFit->SetLineColor(2);
        
        cout << __LINE__ << endl;

		// ------------------------------ NLO Calculations -----------------------------
		histo2DDoubleRatioPlotting->DrawCopy();
		NLODoubleRatio->RemovePoint(0);
		One->Draw("same");
		NLODoubleRatio->Draw("lp3");

		DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurity, 20, 2.0, kBlack, kBlack);   
		DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurity, 20, 2.0, kBlue+2, kBlue+2);   
		histoDoubleRatioConversionTrueEffPurity->DrawCopy("same");
		histoDoubleRatioFitPi0YieldPurity->DrawCopy("same");

		TLegend* legendDoubleConversionFit2             = GetAndSetLegend(0.14,0.7,4,1,"Direct Photon Signal via Conversions");
		legendDoubleConversionFit2->AddEntry(histoDoubleRatioConversionTrueEffPurity,"Measured Direct Photon Signal","p");
		legendDoubleConversionFit2->AddEntry(histoDoubleRatioFitPi0YieldPurity,"Measured Direct Photon Signal, fitted #pi^{0}","p");
		if(option.CompareTo("PbPb_2.76TeV") == 0) legendDoubleConversionFit2->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 2.76TeV scaled N_{coll}","l");
		else if(option.CompareTo("pPb_5.023TeV") == 0) legendDoubleConversionFit2->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 5.023TeV scaled N_{coll}","l");
		else legendDoubleConversionFit2->AddEntry(NLODoubleRatio,"pp NLO Direct Photon","l");
		legendDoubleConversionFit2->Draw();
        //legendDoubleConversionFit2->AddEntry(NLOdoubleRatioFit,"Fit to NLO Double Ratio","l");

		canvasConversionFitDoubleRatioSum->Print(Form("%s/%s_DoubleRatioFit_Sum_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));
	
        /*
        TCanvas *canvasConversionDoubleRatioSum         = GetAndSetCanvas("canvasConversionDoubleRatioSum");
        
 		Double_t binningDoubleRatioHIPions[14]          = {1.6,1.8,2.0,2.2,2.4,2.6,3.0,3.5,4.0,5.0,6.0,8.0,11.0,14.0};
  		TH1D *histoBinningDoubleRatioPions              = new TH1D("","",13,binningDoubleRatioHIPions);
 
 		SetHistogramm(histoDoubleRatioConversionTrueEffPurity,"p_{T} (GeV/c)","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})",0.7,1.8);
 		DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurity, 20, 2.0, 1, 1);
 	
        if(!option.CompareTo("PbPb_2.76TeV")){
 			histoDoubleRatioConversionTrueEffPurity->GetYaxis()->SetRangeUser(0.5,3.5);
 		}
 		else{
 			histoDoubleRatioConversionTrueEffPurity->GetYaxis()->SetRangeUser(0.7,2.5);
 		}
 
 		histoDoubleRatioConversionTrueEffPurity->Draw("");
 		NLODoubleRatio->Draw("lp3");
 		histoDoubleRatioConversionTrueEffPurity->Draw("same");
 		One->Draw("same");
 		DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurity, 20, 2.0, 1, 1);
 		histoDoubleRatioConversionTrueEffPurity->DrawCopy("same,e1");
 		
 		TLegend* legendDoubleConversion                 = GetAndSetLegend(0.16,0.80, 3);
 		legendDoubleConversion->AddEntry(histoDoubleRatioConversionTrueEffPurity,"Direct photons signal with neutral pions");
 		if(option.CompareTo("PbPb_2.76TeV") == 0) legendDoubleConversion->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 2.76TeV scaled N_{coll}","l");
 		else if(option.CompareTo("pPb_5.023TeV") == 0) legendDoubleConversion->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 5.023TeV scaled N_{coll}","l");
 		else legendDoubleConversion->AddEntry(NLODoubleRatio,"pp NLO Direct Photon","l");
 
 		legendDoubleConversion->Draw();
 
 		canvasConversionDoubleRatioSum->Print(Form("%s/%s_DoubleRatio_Sum_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));
        */

		// --------------------- Direct Photon Spectrum -------------------------
		TH1D *histoDirectPhotonSpectrum                 = (TH1D*)histoGammaSpecCorrPurity->Clone("histoDirectPhotonSpectrum");
		TH1D *histoThermalPhotonSpectrum                = (TH1D*)histoGammaSpecCorrPurity->Clone("histoThermalPhotonSpectrum");
		TH1D *histoPromptPhotonSpectrum                 = (TH1D*)histoGammaSpecCorrPurity->Clone("histoPromptPhotonSpectrum");
		
        histoDirectPhotonSpectrum->GetYaxis()->SetRangeUser(1e-7,200);
		histoThermalPhotonSpectrum->GetYaxis()->SetRangeUser(1e-7,200);
        histoPromptPhotonSpectrum->GetYaxis()->SetRangeUser(1e-7,200);

        DrawGammaSetMarker(histoDirectPhotonSpectrum, 24, 1.1, 1,1);
		DrawGammaSetMarker(histoThermalPhotonSpectrum, 21, 1.1, kRed,kRed);
		DrawGammaSetMarker(histoPromptPhotonSpectrum, 24, 1.1, kBlue,kBlue);
		
        TCanvas* DirectPhotons                          = GetAndSetCanvas("DirectPhotons",0.15,0.1,1000 ,1500);
		DirectPhotons->SetLogy();

		for(Int_t i = 1; i < histoDirectPhotonSpectrum->GetNbinsX(); i++){
			histoDirectPhotonSpectrum->SetBinContent(i+1,0);
			histoThermalPhotonSpectrum->SetBinContent(i+1,0);
		}
		for(Int_t i = 1; i < histoDirectPhotonSpectrum->GetNbinsX(); i++){
			Double_t binContent                         = -1;
			Double_t binError                           = -1;
			binContent                                  = 1-(1./( histoDoubleRatioConversionTrueEffPurity->GetBinContent(i+1)));
			binContent                                  = binContent*histoGammaSpecCorrPurity->GetBinContent(i+1);
			binError                                    = 1-(1./( histoDoubleRatioConversionTrueEffPurity->GetBinContent(i+1) + histoDoubleRatioConversionTrueEffPurity->GetBinError(i+1)));
			binError                                    = binError*histoGammaSpecCorrPurity->GetBinContent(i+1);
			binError                                    = binError-binContent;
			histoDirectPhotonSpectrum->SetBinContent(i+1,binContent);
			histoDirectPhotonSpectrum->SetBinError(i+1,binError);

			if(histoDirectPhotonSpectrum->GetBinCenter(i+1) < 4.5){
				binContent                              = 1-(1./(1+ histoDoubleRatioConversionTrueEffPurity->GetBinContent(i+1) - NLOdoubleRatioFit->Eval(histoDoubleRatioConversionTrueEffPurity->GetBinCenter(i+1))));
				binContent                              = binContent*histoGammaSpecCorrPurity->GetBinContent(i+1);
				histoThermalPhotonSpectrum->SetBinContent(i+1,binContent);
				histoThermalPhotonSpectrum->SetBinError(i+1,binError);
			}
			else{
				histoThermalPhotonSpectrum->SetBinContent(i+1,0);
				histoThermalPhotonSpectrum->SetBinError(i+1,0);
			}
			if(histoDirectPhotonSpectrum->GetBinCenter(i+1) < 4.5){
				binContent                              = 1-(1./(NLOdoubleRatioFit->Eval(histoDoubleRatioConversionTrueEffPurity->GetBinCenter(i+1))));
				binContent                              = binContent*histoGammaSpecCorrPurity->GetBinContent(i+1);
				histoPromptPhotonSpectrum->SetBinContent(i+1,binContent);
				histoPromptPhotonSpectrum->SetBinError(i+1,0.0000001);
			}
			else{
				histoPromptPhotonSpectrum->SetBinContent(i+1,0);
				histoPromptPhotonSpectrum->SetBinError(i+1,0);
			}
		}

		histo2DSpectraPi0->GetYaxis()->SetLabelSize(0.035);
		histo2DSpectraPi0->GetXaxis()->SetLabelSize(0.035);
		histo2DSpectraPi0->GetYaxis()->SetTitleSize(0.035);
		histo2DSpectraPi0->GetXaxis()->SetTitleSize(0.035);
		histo2DSpectraPi0->GetYaxis()->SetTitleOffset(1.8);
		histo2DSpectraPi0->GetXaxis()->SetTitleOffset(1.4);

		histo2DSpectraPi0->GetXaxis()->SetRangeUser(0.6,5);
		histo2DSpectraPi0->GetYaxis()->SetRangeUser(1e-4,10);

		histo2DSpectraPi0->DrawCopy();

		histoDirectPhotonSpectrum->Draw("same");
		histoThermalPhotonSpectrum->Draw("same");

		histoPromptPhotonSpectrum->GetXaxis()->SetRangeUser(1.2,14);
		histoPromptPhotonSpectrum->Draw("lpsame");

		TF1 *exponetialA                                = new TF1("Exponential Fit A","[0]*exp(-x/[1])",0.8,2.2);
		TF1 *exponetialB                                = new TF1("Exponential Fit B","[0]*exp(-x/[1])",0.8,2.2);
		TF1 *exponetialC                                = new TF1("Exponential Fit B","[0]*exp(-x/[1])",0.8,3.5);
		TF1 *exponetialD                                = new TF1("Exponential Fit B","[0]*exp(-x/[1])",0.8,3.5);
        
		exponetialA->SetParameters(1.,220);
		exponetialB->SetParameters(1.,0.3);
		exponetialC->SetParameters(1.,0.3);
		exponetialD->SetParameters(1.,0.3);

		histoDirectPhotonSpectrum->Fit(exponetialA,"QNRME+","",0.8,2.2);
		histoThermalPhotonSpectrum->Fit(exponetialB,"QNRME+","",0.8,2.2);
		histoThermalPhotonSpectrum->Fit(exponetialC,"QNRME+","",0.8,3.5);
		histoDirectPhotonSpectrum->Fit(exponetialD,"QNRME+","",0.8,3.5);

		exponetialA->SetLineColor(kGreen-2);
		exponetialA->SetLineStyle(1);
		exponetialB->SetLineColor(kOrange-3);
		exponetialB->SetLineStyle(7);
		exponetialC->SetLineColor(kMagenta-2);
		exponetialC->SetLineStyle(3);
		exponetialD->SetLineColor(kYellow-1);
		exponetialD->SetLineStyle(5);

		exponetialA->Draw("same");
		exponetialD->Draw("same");
		exponetialC->Draw("same");
		exponetialB->Draw("same");

		NLO->SetMarkerSize(0);
		NLO->Draw("same");
		histoThermalPhotonSpectrum->Draw("same");
		histoDirectPhotonSpectrum->Draw("same");
		TLegend* legendDirectPhoton                     = GetAndSetLegend(0.37,0.665, 6);
		legendDirectPhoton->SetTextSize(0.031);
		legendDirectPhoton->SetTextFont(42);
		legendDirectPhoton->AddEntry(histoDirectPhotonSpectrum,"Direct Photon Spectrum","p");
		legendDirectPhoton->AddEntry(histoThermalPhotonSpectrum,"Direct - NLO double ratio fit","p");
		legendDirectPhoton->AddEntry(NLO,"NLO Photons (PDF: CTEQ6M5 FF: DSS)","l");
		legendDirectPhoton->AddEntry(histoPromptPhotonSpectrum,"NLO Photons from doube ratio fit","p");
		legendDirectPhoton->AddEntry((TObject*)0,"Exponential fits:","");
		legendDirectPhoton->AddEntry(exponetialA,Form("Direct Photons, T = %.0f #pm %.0f^{stat} MeV",exponetialA->GetParameter(1)*1000,exponetialA->GetParError(1)*1000),"l");
		legendDirectPhoton->AddEntry(exponetialD,Form("Up to 3.5 GeV/c, T = %.0f #pm %.0f^{stat} MeV",exponetialD->GetParameter(1)*1000,exponetialD->GetParError(1)*1000),"l");
		legendDirectPhoton->AddEntry(exponetialB,Form("Subtracted NLO, T = %.0f #pm %.0f^{stat} MeV",exponetialB->GetParameter(1)*1000,exponetialB->GetParError(1)*1000),"l");
		legendDirectPhoton->AddEntry(exponetialC,Form("Up to 3.5 GeV/c, T = %.0f #pm %.0f^{stat} MeV",exponetialC->GetParameter(1)*1000,exponetialC->GetParError(1)*1000),"l");
		legendDirectPhoton->Draw();

		DirectPhotons->Print(Form("%s/%s_DirectPhotons_Sum_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

		histo2DSpectraPi0->GetXaxis()->SetRangeUser(0.6,15);
		histo2DSpectraPi0->GetYaxis()->SetRangeUser(1e-7,10);

        cout << __LINE__ << endl;

		TCanvas* DirectPhotons2                         = GetAndSetCanvas("DirectPhotons",0.15,0.1);
		DirectPhotons2->SetLogy();
        
		histo2DSpectraPi0->DrawCopy();

		histoDirectPhotonSpectrum->DrawCopy("same");
		NLO->Draw("same");

		//legendDirectPhoton2->Draw();

		DirectPhotons2->Print(Form("%s/%s_DirectPhotons2_Sum_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));

        cout << __LINE__ << endl;

		// Inclusive Photon RAA
        /*
        if(option.CompareTo("2.76TeV") == 0){
         TCanvas *canvasIncPhotonRaa                     = GetAndSetCanvas("canvasIncPhotonRaa");
         TFile *fileGammapp                              = new TFile("pp_2.76TeV_a/0000010002093250003800000_01631031009/2.76TeV/Gamma_Pi0_data_GammaConvV1Correction_0000010002093250003800000_01631031009.root");

         TFile *RAAFile                                  = new TFile("RAAFile.root");

         TH1D *Raa020                                    = (TH1D*) RAAFile->Get("Raa_0020");
         TH1D *Raa2040                                   = (TH1D*) RAAFile->Get("Raa_0240");
         TH1D *Raa4080                                   = (TH1D*) RAAFile->Get("Raa_0480");

         TH1D *histoGammaSpecCorrPuritypp                = (TH1D*) fileGammapp->Get("GammaUnfold");
         TH1D *IncPhotonRaa                              = (TH1D*) histoGammaSpecMCAll->Clone("IncPhotonRaa");

         histoGammaSpecCorrPuritypp->Scale(fcmult);
         IncPhotonRaa->Divide(histoGammaSpecCorrPurity,histoGammaSpecCorrPuritypp,1,1);

         SetHistogramm(Raa020,"p_{T} (GeV/c)","Inc Photon R_{AA}",0.0,1.05);
         SetHistogramm(Raa2040,"p_{T} (GeV/c)","Inc Photon R_{AA}",0.0,1.05);
         SetHistogramm(Raa4080,"p_{T} (GeV/c)","Inc Photon R_{AA}",0.0,1.05);

         DrawGammaSetMarker(Raa020, 20, 1.7, 1,1);
         DrawGammaSetMarker(Raa2040, 20, 1.7, 4,4);
         DrawGammaSetMarker(Raa4080, 20, 1.7, 2,2);

         Raa020->Draw("e1");
         Raa2040->Draw("e1,same");
         Raa4080->Draw("e1,same");
         One->Draw("same");

         TLegend* legendRAA                              = GetAndSetLegend(0.15,0.6,3,"Inclusive photon R_{AA} 2.76TeV");
         legendRAA->AddEntry(Raa020,"0-20%","p");
         legendRAA->AddEntry(Raa2040,"20-40%","p");
         legendRAA->AddEntry(Raa4080,"40-80%","p");
         legendRAA->Draw();
         
         canvasIncPhotonRaa->Print(Form("%s/%s_IncPhotonRaa_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));
        }
        */

        fileFinalResults.close();

        cout << __LINE__ << endl;
        
		// for systematics of cocktail:
        if ( (fcocktailFuncModA.CompareTo("") != 0) && (fcocktailFuncModB.CompareTo("") != 0) ) {
            cocktailFileModA                            = new TFile(Form("CocktailInput/%s.root",fcocktailFuncModA.Data()));
            cocktailFileModB                            = new TFile(Form("CocktailInput/%s.root",fcocktailFuncModB.Data()));
            cocktailDirModA                             = (TDirectoryFile*) cocktailFileModA->Get(fcocktailFuncModA);
            cocktailDirModB                             = (TDirectoryFile*) cocktailFileModB->Get(fcocktailFuncModB);
            
            cocktailAllGammaPi0ModA                     = (TH1D*) cocktailDirModA->Get("sumgammapi0");
            cocktailAllGammaModA                        = (TH1D*) cocktailDirModA->Get("ptg");
            cocktailPi0ModA                             = (TH1D*) cocktailDirModA->Get("ptPi0");
            cocktailAllGammaPi0ModA                     = RebinTH1D(cocktailAllGammaPi0ModA,histoIncRatioPurity);
            TH1D *cocktailAllGammaRebinedModA           = RebinTH1D(cocktailAllGammaModA,histoIncRatioPurity);
            TH1D *cocktailPi0RebinedModA                = RebinTH1D(cocktailPi0ModA,histoIncRatioPurity);
            cocktailAllGammaPi0ModA->Divide(cocktailAllGammaRebinedModA,cocktailPi0RebinedModA);
            histoDoubleRatioFitPi0YieldPurityModA       = (TH1D*) histoIncRatioFitPurity->Clone("DoubleRatioConversionFitPurityModA");
            histoDoubleRatioFitPi0YieldPurityModA->Divide(cocktailAllGammaPi0ModA);
            histoDoubleRatioConversionTrueEffPurityModA = (TH1D*) histoIncRatioPurityTrueEff->Clone("DoubleRatioConversionTrueEffPurityModA");
            histoDoubleRatioConversionTrueEffPurityModA->Divide(cocktailAllGammaPi0ModA);
            
            cocktailAllGammaPi0ModB                     = (TH1D*) cocktailDirModB->Get("sumgammapi0");
            cocktailAllGammaModB                        = (TH1D*) cocktailDirModB->Get("ptg");
            cocktailPi0ModB                             = (TH1D*) cocktailDirModB->Get("ptPi0");
            cocktailAllGammaPi0ModB                     = RebinTH1D(cocktailAllGammaPi0ModB,histoIncRatioPurity);
            TH1D *cocktailAllGammaRebinedModB           = RebinTH1D(cocktailAllGammaModB,histoIncRatioPurity);
            TH1D *cocktailPi0RebinedModB                = RebinTH1D(cocktailPi0ModB,histoIncRatioPurity);
            cocktailAllGammaPi0ModB->Divide(cocktailAllGammaRebinedModB,cocktailPi0RebinedModB);
            histoDoubleRatioFitPi0YieldPurityModB       = (TH1D*) histoIncRatioFitPurity->Clone("DoubleRatioConversionFitPurityModB");
            histoDoubleRatioFitPi0YieldPurityModB->Divide(cocktailAllGammaPi0ModB);
            histoDoubleRatioConversionTrueEffPurityModB = (TH1D*) histoIncRatioPurityTrueEff->Clone("DoubleRatioConversionTrueEffPurityModB");
            histoDoubleRatioConversionTrueEffPurityModB->Divide(cocktailAllGammaPi0ModB);
        }
        
        cout << __LINE__ << endl;
		
		TString nameOutputFile                          = Form("%s/%s/%s_%s_GammaConvV1_InclusiveRatio_%s.root",cutSel.Data(),option.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data());
		fileCorrectedOutput                             = new TFile(nameOutputFile,"RECREATE");
        cout << __LINE__ << endl;
        
			// pi0 quantities
			histoTruePi0MCData->Write();
			histoNormalPi0MCData->Write();
			fitPi0YieldA->Write();
			fitPi0YieldB->Write();
			fitPi0YieldC->Write();
			histoCorrectedPi0Yield->Write();
			histoCorrectedPi0YieldWide->Write();
			histoCorrectedPi0YieldNarrow->Write(); 
			histoCorrectedPi0YieldFit->Write();
			histoCorrectedPi0YieldFitWide->Write();
			histoCorrectedPi0YieldFitNarrow->Write();
			
        cout << __LINE__ << endl;

			// gamma quantities
			histoPurityGammaMCData->Write();
			ConversionGammaFitA->Write();
			histoDirectPhotonSpectrum->Write();
			//histoGammaSpecCorrPurity->Write("histoGammaSpecCorrPurity");
			histoGammaSpecCorrPurity->Write("histoGammaSpecCorrPurity");
				
        cout << __LINE__ << endl;

			// Double ratio
			if(histoDoubleRatioFitPi0YieldPurity) histoDoubleRatioFitPi0YieldPurity->Write();
			if(histoDoubleRatioFitPi0YieldPurityWide) histoDoubleRatioFitPi0YieldPurityWide->Write();
			if(histoDoubleRatioFitPi0YieldPurityNarrow) histoDoubleRatioFitPi0YieldPurityNarrow->Write();
			if(histoDoubleRatioConversionHighFitPurity) histoDoubleRatioConversionHighFitPurity->Write();
			if(histoDoubleRatioConversionLowFitPurity) histoDoubleRatioConversionLowFitPurity->Write();
			if(histoDoubleRatioFitPi0YieldPurityModA) histoDoubleRatioFitPi0YieldPurityModA->Write();
			if(histoDoubleRatioFitPi0YieldPurityModB) histoDoubleRatioFitPi0YieldPurityModB->Write();
			//histoDoubleRatioCombinedPurity->Write();
            if(histoDoubleRatioConversionTrueEffPurity) histoDoubleRatioConversionTrueEffPurity->Write();
			if(histoDoubleRatioConversionTrueEffPurityWide) histoDoubleRatioConversionTrueEffPurityWide->Write();
			if(histoDoubleRatioConversionTrueEffPurityNarrow) histoDoubleRatioConversionTrueEffPurityNarrow->Write();
			if(histoDoubleRatioConversionTrueEffPurityModA) histoDoubleRatioConversionTrueEffPurityModA->Write();
			if(histoDoubleRatioConversionTrueEffPurityModB) histoDoubleRatioConversionTrueEffPurityModB->Write();
			if(histoDoubleRatioConversionOnlyGamma) histoDoubleRatioConversionOnlyGamma->Write();
			
        cout << __LINE__ << endl;

			// inclusive ratio
			histoIncRatioPurity->Write();
			histoIncRatioPurityTrueEff->Write();
			histoIncRatioPurityTrueEffWide->Write();
			histoIncRatioPurityTrueEffNarrow->Write();
			histoMCIncRatio->Write();
			histoIncRatioGammaMC->Write();
			histoIncRatioFitPurity->Write();
			histoIncRatioFitPurityWide->Write();
			histoIncRatioFitPurityNarrow->Write();
			histoIncRatioLowFitPurity->Write();
			histoIncRatioHighFitPurity->Write();
			//histoMCesdIncRatio->Write();
			
        cout << __LINE__ << endl;

			// cocktail
			cocktailAllGammaPi0original->Write("sumgammapi0");
			cocktailPi0GammaPi0->Write();
			cocktailEtaGammaPi0->Write();
			cocktailOmegaGammaPi0->Write();
			cocktailEtapGammaPi0->Write();
			cocktailPhiGammaPi0->Write();
			cocktailRhoGammaPi0->Write();
			cocktailAllGamma->Write();
			cocktailPi0Gamma->Write();
			cocktailEtaGamma->Write();
			cocktailEtapGamma->Write();
			cocktailOmegaGamma->Write();
			cocktailPhiGamma->Write();
			cocktailRhoGamma->Write();
			cocktailPi0->Write();

        cout << __LINE__ << endl;

			// Mc cocktail
			histoDecayRatioSumGamma->Write();
			histoDecayRatioPi0Gamma->Write();
			histoDecayRatioEtaGamma->Write();
			histoDecayRatioEtapGamma->Write();
			histoDecayRatioOmegaGamma->Write();
			histoDecayRatioRho0Gamma->Write();
			histoDecayRatioPhiGamma->Write();
			histoMCDoubleRatioSum->Write();
			
        cout << __LINE__ << endl;

			//NLO
			half->Write("doubleRatioNLOhalf");
			one->Write("doubleRatioNLOone");
			two->Write("doubleRatioNLOtwo");
			graphNLOCalcMuTwo->Write("graphNLOCalcMuTwo");
			graphNLOCalcMuOne->Write("graphNLOCalcMuOne");
			graphNLOCalcMuHalf->Write("graphNLOCalcMuHalf");
			NLODoubleRatio->Write("graphNLODoubleRatio");
			NLO->Write("graphNLODirGamma");
		fileCorrectedOutput->Close();

        cout << __LINE__ << endl;
        
		cout<<"CorrectedYieldTrueEff "<<fitPi0A<<endl;
		cout<<fitPi0YieldA->GetChisquare()/fitPi0YieldA->GetNDF()<<endl;
		cout<<fitPi0YieldA->GetParameter(0)<<endl;
		cout<<fitPi0YieldA->GetParameter(1)<<endl;
		cout<<fitPi0YieldA->GetParameter(2)<<endl;
		cout<<fitPi0YieldA->GetParameter(3)<<endl;
		cout<<fitPi0YieldA->GetParameter(4)<<endl;
		cout<<fitPi0YieldA->GetParameter(5)<<endl;

		cout<<"CorrectedYieldTrueEff "<<fitPi0B<<endl;
		cout<<fitPi0YieldB->GetChisquare()/fitPi0YieldB->GetNDF()<<endl;
		cout<<fitPi0YieldB->GetParameter(0)<<endl;
		cout<<fitPi0YieldB->GetParameter(1)<<endl;
		cout<<fitPi0YieldB->GetParameter(2)<<endl;
		cout<<fitPi0YieldB->GetParameter(3)<<endl;
		cout<<fitPi0YieldB->GetParameter(4)<<endl;

		cout<<"CorrectedYieldTrueEff "<<fitPi0C<<endl;
		cout<<fitPi0YieldC->GetChisquare()/fitPi0YieldC->GetNDF()<<endl;
		cout<<fitPi0YieldC->GetParameter(0)<<endl;
		cout<<fitPi0YieldC->GetParameter(1)<<endl;
		cout<<fitPi0YieldC->GetParameter(2)<<endl;
		cout<<fitPi0YieldC->GetParameter(3)<<endl;
		cout<<fitPi0YieldC->GetParameter(4)<<endl;
	}
}