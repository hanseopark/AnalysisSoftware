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
#include "CalculateGammaToPi0V3.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "TMath.h"
#include "TSpline.h"

extern TRandom    *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem    *gSystem;


void  CalculateGammaToPi0V3(    TString nameFileGamma   = "",
                                TString nameFilePi0     = "",
                                TString nameFileCocktail= "",
                                TString cutSel          = "",
                                TString suffix          = "pdf",
                                TString nameMeson       = "",
                                TString isMC            = "",
                                TString option          = "",
                                TString fEstimatePileup = "",
                                Int_t mode              = 0
                            ){
    // switch systematics on/off
    Bool_t doSysErr                             = kTRUE;
    if (!option.CompareTo("7TeV") && mode == 4)
        doSysErr                                = kFALSE;
    
// Setting the general style
   gROOT->Reset();
   StyleSettingsThesis();
   SetPlotStyle();

// Separating cutnumber and retrieving centrality and number of events
   SeparateCutnumberString(cutSel,mode);

// Create strings for naming
   CreateNamingStrings(nameMeson,isMC);
   TString collisionSystem                      = ReturnFullCollisionsSystem(option);
   TString detectionProcess                     = ReturnFullTextReconstructionProcess(mode);
   
// Creating output directory and output file
   cout << "Output directory with plots:" << endl;
   cout << Form("%s/%s/%s/GammaToPi0",cutSel.Data(),option.Data(),suffix.Data()) << endl;
   TString outputDir                            = Form("%s/%s/%s/GammaToPi0",cutSel.Data(),option.Data(),suffix.Data());
   gSystem->Exec("mkdir "+outputDir);

   TString nameFinalResDat                      = Form("%s/%s/Gamma_%s_FinalExtraction_%s.dat",cutSel.Data(),option.Data(), textPrefix2.Data(), cutSel.Data());
   fstream fileFinalResults;
   fileFinalResults.open(nameFinalResDat.Data(), ios::out);

   if (!fEstimatePileup.CompareTo("EstimateTrainPileUp"))
      kDoPileup                                 = kTRUE;

// Opening gamma file and loading data gamma spectrum
   fileGamma                                    = new TFile(nameFileGamma);
   if (!option.CompareTo("7TeV") || !option.CompareTo("13TeV")){
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
    histoMCDecaySumGammaPt                      = (TH1D*)fileGamma->Get("MC_DecayGammaAll_Pt");
    histoGammaSpecMCAll                         = (TH1D*)fileGamma->Get("GammaSpecMC");

// Opening pion file and loading spectra
    filePi0                                     = new TFile(nameFilePi0);
    histoCorrectedPi0YieldNormalEff             = (TH1D*)filePi0->Get("CorrectedYieldNormEff");
    histoCorrectedPi0Yield                      = (TH1D*)filePi0->Get("CorrectedYieldTrueEff");
    histoCorrectedPi0YieldWide                  = (TH1D*)filePi0->Get("CorrectedYieldTrueEffWide");
    histoCorrectedPi0YieldNarrow                = (TH1D*)filePi0->Get("CorrectedYieldTrueEffNarrow");
    histoMCYieldMeson                           = (TH1D*)filePi0->Get("MCYield_Meson");
    histoMCYieldMesonOldBin                     = (TH1D*)filePi0->Get("MCYield_Meson_oldBin");

// Creating inclusive gamma ratio
    histoIncRatioPurityTrueEff                  = (TH1D*) histoGammaSpecCorrPurity->Clone("IncRatioPurity_trueEff");
    histoIncRatioPurityTrueEff->Divide(histoIncRatioPurityTrueEff,histoCorrectedPi0Yield,1,1,"");

    histoMCIncRatio                             = (TH1D*) histoGammaSpecMCAll->Clone("MC_IncRatio");
    histoMCIncRatio->Divide(histoGammaSpecMCAll,histoMCYieldMeson,1,1,"");

// Opening systematics file
   TString fileNameSysErrGamma                  ="GammaSystematicErrorsCalculated/SystematicErrorAveraged_Gamma_7TeV_2016_11_30.dat"; // default
   TString fileNameSysErrInclRatio              ="GammaSystematicErrorsCalculated/SystematicErrorAveraged_IncRatio_7TeV_2016_11_30.dat"; // default
   TString fileNameSysErrDoubleRatio            ="GammaSystematicErrorsCalculated/SystematicErrorAveraged_DoubleRatio_7TeV_2016_11_30.dat"; // default
   if(option.CompareTo("7TeV") == 0){
        fileNameSysErrGamma                     = "GammaSystematicErrorsCalculated/SystematicErrorAveraged_Gamma_7TeV_2016_12_15.dat";
        fileNameSysErrInclRatio                 = "GammaSystematicErrorsCalculated/SystematicErrorAveraged_IncRatio_7TeV_2016_12_15.dat";
        fileNameSysErrDoubleRatio               = "GammaSystematicErrorsCalculated/SystematicErrorAveraged_DoubleRatio_7TeV_2016_12_15.dat";
   }
   fileSysErrGamma.open(fileNameSysErrGamma,ios_base::in);
   cout << fileNameSysErrGamma << endl;

   while(!fileSysErrGamma.eof() && nPointsGamma < 100){
       fileSysErrGamma >> relSystErrorGammaDown[nPointsGamma] >> relSystErrorGammaUp[nPointsGamma]>>    relSystErrorWOMaterialGammaDown[nPointsGamma] >> relSystErrorWOMaterialGammaUp[nPointsGamma];
        cout << nPointsGamma << "\t"  << relSystErrorGammaDown[nPointsGamma] << "\t"  <<relSystErrorGammaUp[nPointsGamma] << "\t" << relSystErrorWOMaterialGammaDown[nPointsGamma] << "\t"  <<relSystErrorWOMaterialGammaUp[nPointsGamma] << endl;;
        nPointsGamma++;
   }
   fileSysErrGamma.close();
   nPointsGamma = nPointsGamma-1;

   fileSysErrInclRatio.open(fileNameSysErrInclRatio,ios_base::in);
   cout << fileNameSysErrInclRatio << endl;

   while(!fileSysErrInclRatio.eof() && nPointsInclRatio < 100){
       fileSysErrInclRatio >> relSystErrorInclRatioDown[nPointsInclRatio] >> relSystErrorInclRatioUp[nPointsInclRatio]>>    relSystErrorWOMaterialInclRatioDown[nPointsInclRatio] >> relSystErrorWOMaterialInclRatioUp[nPointsInclRatio];
        cout << nPointsInclRatio << "\t"  << relSystErrorInclRatioDown[nPointsInclRatio] << "\t"  <<relSystErrorInclRatioUp[nPointsInclRatio] << "\t" << relSystErrorWOMaterialInclRatioDown[nPointsInclRatio] << "\t"  <<relSystErrorWOMaterialInclRatioUp[nPointsInclRatio] << endl;;
        nPointsInclRatio++;
   }
   fileSysErrInclRatio.close();
   nPointsInclRatio = nPointsInclRatio-1;

   fileSysErrDoubleRatio.open(fileNameSysErrDoubleRatio,ios_base::in);
   cout << fileNameSysErrDoubleRatio << endl;

   while(!fileSysErrDoubleRatio.eof() && nPointsDoubleRatio < 100){
       fileSysErrDoubleRatio >> relSystErrorDoubleRatioDown[nPointsDoubleRatio] >> relSystErrorDoubleRatioUp[nPointsDoubleRatio]>>    relSystErrorWOMaterialDoubleRatioDown[nPointsDoubleRatio] >> relSystErrorWOMaterialDoubleRatioUp[nPointsDoubleRatio];
        cout << nPointsDoubleRatio << "\t"  << relSystErrorDoubleRatioDown[nPointsDoubleRatio] << "\t"  <<relSystErrorDoubleRatioUp[nPointsDoubleRatio] << "\t" << relSystErrorWOMaterialDoubleRatioDown[nPointsDoubleRatio] << "\t"  <<relSystErrorWOMaterialDoubleRatioUp[nPointsDoubleRatio] << endl;;
        nPointsDoubleRatio++;
   }
   fileSysErrDoubleRatio.close();
   nPointsDoubleRatio = nPointsDoubleRatio-1;
   
    if (doSysErr) {
        graphGammaYieldSysErr   = CalculateSysErrFromRelSysHisto( histoGammaSpecCorrPurity, "Pi0SystError",relSystErrorGammaDown , relSystErrorGammaUp, 2, nPointsGamma);
        graphInclRatioSysErr    = CalculateSysErrFromRelSysHisto( histoIncRatioPurityTrueEff, "Pi0SystErrorA",relSystErrorInclRatioDown , relSystErrorInclRatioUp, 2, nPointsInclRatio);
    } else {
        graphGammaYieldSysErr   = NULL;
        graphInclRatioSysErr    = NULL;
    }
   
//**********************************************************************************
//***                      NLO Calculatins Ratio                                 ***
//**********************************************************************************
    Bool_t doNLOComparison                      = kTRUE;
    if (doNLOComparison){
        fileTheoryCompilation                   = new TFile("ExternalInput/Theory/TheoryCompilationPP.root");
        directoryGamma                          = (TDirectoryFile*)fileTheoryCompilation->Get("DirectPhoton");
   // Load NLO input from TheoryCompilationPP.root
        if (!option.CompareTo("900GeV")) {
            cout << "ERROR: No calculations in theory compilation file yet!" << endl;
            return;
        } else if (!option.CompareTo("2.76TeV") || !option.CompareTo("PbPb_2.76TeV")) {
            if (!triggerCutNumber.CompareTo("1") && !subTriggerCutNumber.CompareTo("0")) {   // kINT7
                graphDirectPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT7_2760GeV");
                graphPromptPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT7_2760GeV");
                graphFragmentationPhotonNLO     = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT7_2760GeV");
            } else {
                graphDirectPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT1_2760GeV");
                graphPromptPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT1_2760GeV");
                graphFragmentationPhotonNLO     = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT1_2760GeV");
            }
        } else if (!option.CompareTo("5TeV") || !option.CompareTo("pPb_5.023TeV")) {
            cout << "ERROR: No inv. yield calculations in compilation file yet!" << endl;
            return;
        } else if (!option.CompareTo("7TeV")) {
            if (!triggerCutNumber.CompareTo("1") && !subTriggerCutNumber.CompareTo("0")) {   // kINT7
                graphDirectPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT7_7TeV");
                graphPromptPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT7_7TeV");
                graphFragmentationPhotonNLO     = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT7_7TeV");
            } else {
                graphDirectPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT1_7TeV");
                graphPromptPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT1_7TeV");
                graphFragmentationPhotonNLO     = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT1_7TeV");
            }
        } else if (!option.CompareTo("8TeV")) {
            graphDirectPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYield_8TeV");
            graphPromptPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYield_8TeV");
            graphFragmentationPhotonNLO         = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYield_8TeV");
        } else if (!option.CompareTo("13TeV")) {
            graphDirectPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT7_7TeV");
            graphPromptPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT7_7TeV");
            graphFragmentationPhotonNLO         = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT7_7TeV");
            cout << "WARNING: No 13TeV inv. yield calculations in compilation file yet! (using 7TeV calculations)" << endl;
        } else {
            cout << "ERROR: Energy or collision system not known!" << endl;
            return;
        }

   // scale with Ncoll (for pPb and PbPb)
        graphDirectPhotonNLO                    = (TGraphAsymmErrors*)ScaleGraph(graphDirectPhotonNLO, fNcoll);
        graphPromptPhotonNLO                    = (TGraphAsymmErrors*)ScaleGraph(graphPromptPhotonNLO, fNcoll);
        graphFragmentationPhotonNLO             = (TGraphAsymmErrors*)ScaleGraph(graphFragmentationPhotonNLO, fNcoll);

   // Create canvas and pads
        TCanvas *NLOCalculationcanvas           = GetAndSetCanvas("NLOCalculationcanvas",0.,0.,1000,1350);
        TPad* padNLOHistos                      = new TPad("padNLOHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padNLOHistos, 0.14, 0.017, 0.01, 0.);
        padNLOHistos->Draw();
        TPad* padNLORatios                      = new TPad("padNLORatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padNLORatios, 0.14, 0.017, 0.0, 0.22);
        padNLORatios->Draw();
        padNLOHistos->cd();
        padNLOHistos->SetLogy();

   // Set axis range and labels
        TH2F * histoDummy                       = new TH2F("histoDummy","histoDummy",1000,0., 25,10000,5e-12, 250);
        SetStyleHistoTH2ForGraphs(histoDummy, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                0.85*textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.77,1.7);
        histoDummy->Draw();

        graphDirectPhotonNLO->SetTitle("graphNLOCalcDirectPhoton");
        graphDirectPhotonNLO->GetYaxis()->SetRangeUser(0.1e-12,1);
        DrawGammaSetMarkerTGraphAsym(graphDirectPhotonNLO, 24, 1., kRed-1, kRed-1);

        graphPromptPhotonNLO->SetTitle("graphNLOCalcPromptPhoton");
        graphPromptPhotonNLO->GetYaxis()->SetRangeUser(0.1e-12,1);
        DrawGammaSetMarkerTGraphAsym(graphPromptPhotonNLO, 24, 1., kBlue+2, kBlue+2);

        graphFragmentationPhotonNLO->SetTitle("graphNLOCalcFragmentationPhoton");
        graphFragmentationPhotonNLO->GetYaxis()->SetRangeUser(0.1e-12,1);
        DrawGammaSetMarkerTGraphAsym(graphFragmentationPhotonNLO, 24, 1., kCyan-2, kCyan-2);

        cout << __LINE__ << ": start fitting to NLO calculations" << endl;
        fitNLODirectPhoton                      = FitObject("m","fitNLODirectPhoton","Pi0",graphDirectPhotonNLO,2.,25.);                    // mod power law
        DrawGammaSetMarkerTF1( fitNLODirectPhoton, 8, 2.0, kRed-1);
        fitNLOPromptPhoton                      = FitObject("m","fitNLOPromptPhoton","Pi0",graphPromptPhotonNLO,2.,25.);                    // mod power law
        DrawGammaSetMarkerTF1( fitNLOPromptPhoton, 8, 2.0, kBlue+2);
        fitNLOFragmentationPhoton               = FitObject("m","fitNLOFragmentationPhoton","Pi0",graphFragmentationPhotonNLO,2.,25.);      // mod power law
        DrawGammaSetMarkerTF1( fitNLOFragmentationPhoton, 8, 2.0, kCyan-2);

        histoMCDecaySumGammaPt->Draw("same");
        histoGammaSpecCorrPurity->Draw("same");
        fitNLODirectPhoton->Draw("same");
        graphDirectPhotonNLO->Draw("p,same");
        fitNLOPromptPhoton->Draw("same");
        graphPromptPhotonNLO->Draw("p,same");
        fitNLOFragmentationPhoton->Draw("same");
        graphFragmentationPhotonNLO->Draw("p,same");
        
   // Create legend
        TLegend* leg_NLOCalculationcanvas       = GetAndSetLegend(0.3,0.6,8);
        leg_NLOCalculationcanvas->AddEntry(graphDirectPhotonNLO,       "Direct Photon NLO Calc", "lep");
        leg_NLOCalculationcanvas->AddEntry(fitNLODirectPhoton,         "Fit to Direct Photon NLO Calc");
        leg_NLOCalculationcanvas->AddEntry(graphPromptPhotonNLO,       "Prompt Photon NLO Calc", "lep");
        leg_NLOCalculationcanvas->AddEntry(fitNLOPromptPhoton,         "Fit to Prompt Photon NLO Calc");
        leg_NLOCalculationcanvas->AddEntry(graphFragmentationPhotonNLO,"Fragmentation Photon NLO Calc", "lep");
        leg_NLOCalculationcanvas->AddEntry(fitNLOFragmentationPhoton,  "Fit to Fragmentation Photon NLO Calc");
        leg_NLOCalculationcanvas->AddEntry(histoMCDecaySumGammaPt,     "Decay Photons from MC");
        leg_NLOCalculationcanvas->AddEntry(histoGammaSpecCorrPurity,   "Inclusive Photons from data");
        leg_NLOCalculationcanvas->SetTextSize(0.035);
        leg_NLOCalculationcanvas->Draw();

   // Calculating ratio
        padNLORatios->cd();
        padNLORatios->SetLogy();

        histRatioNLODirectPhoton                = (TH1D*) histoMCDecaySumGammaPt->Clone("histRatioNLODirectPhoton");
        histRatioNLODirectPhoton                = CalculateHistoRatioToFitNLO(histRatioNLODirectPhoton,fitNLODirectPhoton,2.);
        textSizeSpectra                         = 0.1;
        SetStyleHistoTH1ForGraphs(histRatioNLODirectPhoton, "#it{p}_{T} (GeV/#it{c})","#frac{Decay #gamma}{NLO Direct #gamma}",textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.95,0.6);

        histRatioNLODirectPhoton->GetYaxis()->CenterTitle(kTRUE);
        DrawGammaSetMarker(histRatioNLODirectPhoton, 24, 0.5, kRed+2, kRed+2);

        histRatioNLODirectPhoton->Draw();
        NLOCalculationcanvas->Print(Form("%s/%s_NLOCalculations_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
        cout << "NLO calculations done" << endl;
    }


//**********************************************************************************
//***                      Inclusive Ratio                                       ***
//**********************************************************************************
    Bool_t doInclusiveRatios                    = kTRUE;
    if (doInclusiveRatios){
   // Create ratio with true efficiency
        TCanvas* canvasIncRatioPurityTrueEff    = GetAndSetCanvas("canvasIncRatioPurityTrueEff",0.095,0.09,1000,815);
        SetHistogramm(histoIncRatioPurityTrueEff, "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
        DrawGammaSetMarker(histoIncRatioPurityTrueEff, 20, 2.0, 4, 4);
        histoIncRatioPurityTrueEff->Draw();
        if (graphInclRatioSysErr) {
            DrawGammaSetMarkerTGraphAsym(graphInclRatioSysErr,26,0,kBlue-8,kBlue-8,2,kTRUE);
            graphInclRatioSysErr->Draw("p,2,same");
        }
        PlotLatexLegend(0.66, 0.75, 0.045,collisionSystem,detectionProcess);
        canvasIncRatioPurityTrueEff->Print(Form("%s/%s_IncRatioPurity_trueEff_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

        // Create ratio for MC
        TCanvas* canvasMCIncRatio               = GetAndSetCanvas("canvasMCIncRatio",0.095,0.09,1000,815);
        SetHistogramm(histoMCIncRatio, "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
        DrawGammaSetMarker(histoMCIncRatio, 24, 2.0, 2, 2);
        histoMCIncRatio->Draw("same");
        PlotLatexLegend(0.66, 0.75, 0.045,collisionSystem,detectionProcess);
        if(option.CompareTo("PbPb_2.76TeV") == 0) DrawCentrality(0.8,0.9,centrality);
        canvasMCIncRatio->Print(Form("%s/%s_MC_IncRatio_%s.%s",outputDir.Data(),textPi0New.Data(),centralityAdd.Data(),suffix.Data()));
        
   // Create plot with both ratios
        TCanvas* canvasIncRatioAll              = GetAndSetCanvas("canvasIncRatioAll",0.095,0.09,1000,815);
        histoIncRatioPurityTrueEff->Draw("e1");
        histoMCIncRatio->Draw("e1,same");
        TLegend* leg_canvasIncRatioAll          = GetAndSetLegend(0.18,0.7,5);
        leg_canvasIncRatioAll->AddEntry(histoIncRatioPurityTrueEff,"Ratio with true Eff Purity");
        leg_canvasIncRatioAll->AddEntry(histoMCIncRatio,"MC Ratio");
        leg_canvasIncRatioAll->Draw();
        PlotLatexLegend(0.66, 0.75, 0.045,collisionSystem,detectionProcess);
        if(option.CompareTo("PbPb_2.76TeV") == 0) DrawCentrality(0.8,0.9,centrality);
        canvasIncRatioAll->Print(Form("%s/%s_IncRatio_all_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
        cout << "Inclusive ratios plotted" << endl;
    }
  
//**********************************************************************************
//***                      Photon spectra data                                   ***
//**********************************************************************************
    TCanvas* CanvasGammaSpecSingle              = GetAndSetCanvas("CanvasGammaSpecSingle", 0.12, 0.1, 1000 ,1350);
    DrawGammaCanvasSettings( CanvasGammaSpecSingle, 0.16, 0.02, 0.015, 0.07);
    CanvasGammaSpecSingle->SetLogy();
    CanvasGammaSpecSingle->SetLogx();
    
    Double_t        minPt                       = 0.2;
    Double_t        maxPt                       = 20;
    
    histoGammaSpecCorrPurity->GetXaxis()->SetLabelOffset(-1e-2);
    DrawAutoGammaMesonHistos( histoGammaSpecCorrPurity,
                             "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                             kTRUE, 20.,5e-10, kFALSE,
                             kFALSE, -0.004, 0.020,
                             kTRUE, minPt,maxPt,62,0.04,62,0.03,0.7,1.7);
    DrawGammaSetMarker(histoGammaSpecCorrPurity, 20, 1.5,kBlue+2,kBlue+2);
    histoGammaSpecCorrPurity->GetYaxis()->SetRangeUser(3e-9, 10);
    
    histoGammaSpecCorrPurity->DrawCopy("e1,same");
    if (graphGammaYieldSysErr) {
        DrawGammaSetMarkerTGraphAsym(graphGammaYieldSysErr,26,0,kBlue-8,kBlue-8,2,kTRUE);
        graphGammaYieldSysErr->Draw("p,2,same");
    }
    
    PlotLatexLegend(0.64, 0.78, 0.045,collisionSystem,detectionProcess,2);
    drawLatex("#gamma's from ALICE Data", 1.7, 0.000000001, kBlack,0.035);
    
    TLegend* leg_GammaSpectra;
    leg_GammaSpectra                            = GetAndSetLegend(0.2,0.2,2);
    leg_GammaSpectra->AddEntry(histoGammaSpecCorrPurity,"Inclusive photons");
    if(graphGammaYieldSysErr)
        leg_GammaSpectra->AddEntry(graphGammaYieldSysErr,"Systematic uncertainty");
    
    leg_GammaSpectra->Draw();
    
    CanvasGammaSpecSingle->Print(Form("%s/%s_GammaSpectrum_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

  
//**********************************************************************************
//***                      Fitting photon spectrum                               ***
//**********************************************************************************
    Bool_t doPhotonSpectra                      = kTRUE;
    if (doInclusiveRatios){
       TString fitGammaA                        = "h";
       TString fitGammaB                        = "l";
        
   // Get minimum and maximum pT for fits from histogram
        fitMinPt                                = 0.4;
        for (Int_t i=1; i<=histoGammaSpecCorrPurity->GetNbinsX(); i++) {
            if (histoGammaSpecCorrPurity->GetBinContent(i) == 0) {
                continue;
            } else {
                fitMinPt                        = histoGammaSpecCorrPurity->GetXaxis()->GetBinLowEdge(i);
                break;
            }
        }
        fitMaxPt                                = histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX());

   // Special PbPb setting for fits
        if(option.CompareTo("PbPb_2.76TeV") == 0){
            if(centCutNumberI<4){
                fitGammaA                       = "QCD";
                fitGammaB                       = "oHag";
                fitMaxPt                        = 14;
            } else{
                fitGammaA                       = "QCD";
                fitGammaB                       = "oHag";
                fitMaxPt                        = 14;
            }
        }
       
        TString fitOptions                      = "IQNRME+";
        if (!option.CompareTo("7TeV") && mode == 4) {
            fitOptions                          = "QNRME+";
        }

   // Start fitting hagedorn and tsallis
       ConversionGammaFitA                      = FitObject(fitGammaA,"ConversionGammaFitA","Pi0",histoGammaSpecCorrPurity,fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(ConversionGammaFitA, 1, 2.0, kBlue-2);
        fileFinalResults << "CorrectedYieldTrueEff " << fitGammaA << endl;
        forOutput                               = WriteParameterToFile(ConversionGammaFitA);
        fileFinalResults << forOutput.Data() << endl;

        ConversionGammaFitB                     = FitObject(fitGammaB,"ConversionGammaFitB","Pi0",histoGammaSpecCorrPurity,fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(ConversionGammaFitB, 2, 2.0, kRed-3);
        fileFinalResults << "CorrectedYieldTrueEff " << fitGammaB << endl;
        forOutput                               = WriteParameterToFile(ConversionGammaFitB);
        fileFinalResults << forOutput.Data() << endl;

        ConversionGammaFitB->SetLineColor(2);
        ConversionGammaFitB->SetLineStyle(2);
        
   // Create canvas and pads with 3:1 ratio
        TCanvas* ConversionGammaCanvas          = GetAndSetCanvas("ConversionGammaSpeccanvas", 0., 0., 1000 ,1350);
        DrawGammaCanvasSettings( ConversionGammaCanvas, 0.12, 0.015, 0.01, 0.09);
        TPad* padConversionGamma                = new TPad("padConversionGammaHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padConversionGamma, 0.14, 0.017, 0.01, 0.);
        padConversionGamma->Draw();
        TPad* padConversionGammaRatio           = new TPad("padConversionGammaRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padConversionGammaRatio, 0.14, 0.017, 0.0, 0.22);
        padConversionGammaRatio->Draw();
        padConversionGamma->cd();
        padConversionGamma->SetLogy();
        padConversionGamma->SetLogx();
        
   // Upper part of plot with histograms
        textSizeSpectra=0.035;
        Double_t        ptMin                   = 0.3;
        if (mode==4)    ptMin                   = 1.5;
        Double_t        ptMax                   = 40.;
        TH2F * histoDummy2                      = new TH2F("histoDummy2","histoDummy2",1000,ptMin,ptMax,10000,5e-9, 9);
        SetStyleHistoTH2ForGraphs(histoDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                 0.85*textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.77,1.7);
        histoDummy2->GetXaxis()->SetRangeUser(ptMin,ptMax);
        histoDummy2->Draw();
        DrawGammaSetMarker(histoGammaSpecCorrPurity, 24, 2.0, kBlack, kBlack);
        ConversionGammaFitA->SetLineColor(kBlue-2);
        ConversionGammaFitB->SetLineColor(kRed-3);

        histoGammaSpecCorrPurity->Draw("same,e1");
        ConversionGammaFitA->Draw("same");
        ConversionGammaFitB->Draw("Csame");
 
        PlotLatexLegend(0.66, 0.75, 0.045,collisionSystem,detectionProcess,2);
        TLegend* leg_ConversionGammaSpeccanvas  = GetAndSetLegend(0.66,0.65,3);
        leg_ConversionGammaSpeccanvas->SetTextSize(0.04);
        leg_ConversionGammaSpeccanvas->AddEntry(histoGammaSpecCorrPurity,"#gamma data", "lp");
        leg_ConversionGammaSpeccanvas->AddEntry(ConversionGammaFitA,"Hagedorn", "l");
        leg_ConversionGammaSpeccanvas->AddEntry(ConversionGammaFitB,"Levy", "l");
        leg_ConversionGammaSpeccanvas->Draw();
        
   // Lower part of plot with ratio
        padConversionGammaRatio->cd();
        padConversionGammaRatio->SetLogx();

        histRatioConversionGammaA               = (TH1D*) histoGammaSpecCorrPurity->Clone("histRatioConversionGammaA");
        histRatioConversionGammaA               = CalculateHistoRatioToFit(histRatioConversionGammaA,ConversionGammaFitA);
        histRatioConversionGammaB               = (TH1D*) histoGammaSpecCorrPurity->Clone("histRatioConversionGammaB");
        histRatioConversionGammaB               = CalculateHistoRatioToFit(histRatioConversionGammaB,ConversionGammaFitB);
        textSizeSpectra                         = 0.1;

        TH1D* dummy                             = new TH1D("dummy", "dummy",1000, ptMin, ptMax);
        SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})", "data/fit", textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.95,0.6);
        dummy->GetYaxis()->SetRangeUser(0.5, 1.55);
        dummy->GetXaxis()->SetLabelOffset(-1e-2);

        DrawGammaSetMarker(histRatioConversionGammaA, 20, 1.5, kBlue-2, kBlue-2);
        DrawGammaSetMarker(histRatioConversionGammaB, 20, 1.5, kRed-3, kRed-3);

        dummy->Draw();
        histRatioConversionGammaA->Draw("e1,same");
        DrawGammaLines(ptMin, ptMax, 1.0, 1.0,1.0, kGray+2 ,7);
        histRatioConversionGammaB->Draw("e1,same");
        
   // Save plot
        ConversionGammaCanvas->Print(Form("%s/%s_Spectra_ConversionGamma_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
        cout << "Photon spectra plotted" << endl;
    }


//**********************************************************************************
//***                      Fitting pi0 and eta if possible                       ***
//**********************************************************************************
    Bool_t doPionFitting                        = kTRUE;
    if (doPionFitting){
        TString fitPi0A                         = "oHag";
        TString fitPi0B                         = "qcd";
        TString fitPi0C                         = "rad";
        
   // Get min and max pT from histogram
        fitMinPt                                = 0.5;
        for (Int_t i=1; i<=histoCorrectedPi0Yield->GetNbinsX(); i++) {
            if (histoCorrectedPi0Yield->GetBinContent(i) == 0) {
                continue;
            } else {
                fitMinPt                        = histoCorrectedPi0Yield->GetXaxis()->GetBinLowEdge(i);
                break;
            }
        }
        fitMaxPt                                = histoCorrectedPi0Yield->GetXaxis()->GetBinUpEdge(histoCorrectedPi0Yield->GetNbinsX());
        
   // Special fitting options for PbPb
        if(option.CompareTo("PbPb_2.76TeV") == 0){
            if(centCutNumberI<4){
                fitPi0A                         = "oHag";
                fitPi0B                         = "qcd";
                fitPi0C                         = "rad";
            } else{
                fitPi0A                         = "oHag";
                fitPi0B                         = "qcd";
                fitPi0C                         = "rad";
            }
        } else{
            fitPi0A                             = "h";
            fitPi0B                             = "l";
            fitPi0C                             = "oHag";
        }
       
       TString fitOptions                       = "IQNRME+";
       if (!option.CompareTo("7TeV") && mode == 4) {
           fitOptions                           = "QNRME+";
       }

        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"---------------------- Begin Fitting Pi0 ------------------------"<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;
        
   // Fitting Hagedorn
        fitPi0YieldA                            = FitObject(fitPi0A,"fitPi0YieldA","Pi0",histoCorrectedPi0Yield,fitMinPt,fitMaxPt,NULL,fitOptions);
        fitPi0YieldA->SetRange(fitMinPt,fitMaxPt);
        DrawGammaSetMarkerTF1(fitPi0YieldA, 1, 2.0, kBlue+1);
        fileFinalResults << "CorrectedYieldTrueEff hagedorn" << endl;
        forOutput                               = WriteParameterToFile(fitPi0YieldA);
        fileFinalResults << forOutput.Data() << endl;

   // Fitting Levy-Tsallis
        fitPi0YieldB                            = FitObject(fitPi0B,"fitPi0YieldB","Pi0",histoCorrectedPi0Yield,fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(fitPi0YieldB, 2, 2.0, kRed+1);
        fitPi0YieldB->SetRange(fitMinPt,fitMaxPt);
        fileFinalResults << "CorrectedYieldTrueEff Levy" << endl;
        forOutput                               = WriteParameterToFile(fitPi0YieldB);
        fileFinalResults << forOutput.Data() << endl;

   // Fitting modified Hagedorn
        fitPi0YieldC                            = FitObject(fitPi0C,"fitPi0YieldC","Pi0",histoCorrectedPi0Yield,fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(fitPi0YieldC, 3, 2.0, kGreen+1);
        fitPi0YieldC->SetRange(fitMinPt,fitMaxPt);
        fileFinalResults << "CorrectedYieldTrueEff modified Hagedorn" << endl;
        forOutput                               = WriteParameterToFile(fitPi0YieldC);
        fileFinalResults << forOutput.Data() << endl;

        cout<<"------------------------------------------------------------------"<<endl;
        cout<<"----------------------- End Fitting Pi0 --------------------------"<<endl;
        cout<<"------------------------------------------------------------------"<<endl;

        TCanvas *ConversionSpeccanvas           = GetAndSetCanvas("ConversionSpeccanvas",0., 0., 1000 ,1350);
        DrawGammaCanvasSettings( ConversionSpeccanvas, 0.12, 0.015, 0.01, 0.09);
        TPad* padConversionHistos               = new TPad("padConversionHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padConversionHistos, 0.14, 0.017, 0.01, 0.);
        padConversionHistos->Draw();
        TPad* padConversionRatios               = new TPad("padConversionRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padConversionRatios, 0.14, 0.017, 0.0, 0.22);
        padConversionRatios->Draw();

        padConversionHistos->cd();
        padConversionHistos->SetLogy();
        padConversionHistos->SetLogx();

   // Plotting spectrum and fits in upper pad
        textSizeSpectra                         = 0.035;
        Double_t        ptMin                   = 0.2;
        if (mode==4)    ptMin                   = 1.5;
        Double_t        ptMax                   = 40.;
        TH2F * histoDummy3                      = new TH2F("histoDummy3","histoDummy3",1000,ptMin,1.5*histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()),1000,2e-9,999);
        SetStyleHistoTH2ForGraphs(histoDummy3, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                  0.85*textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.77,1.7);
        histoDummy3->DrawCopy();

        DrawGammaSetMarker(histoCorrectedPi0Yield, 20, 1.5, kBlack, kBlack);
        fitPi0YieldA->SetLineColor(kBlue-3);
        fitPi0YieldB->SetLineColor(kRed-3);
        fitPi0YieldC->SetLineColor(kGreen-2);
        
        histoCorrectedPi0Yield->Draw("e1,same");
        fitPi0YieldA->Draw("Csame");
        fitPi0YieldB->Draw("Csame");
        fitPi0YieldC->Draw("same");
        PlotLatexLegend(0.2, 0.06, 0.045,collisionSystem,detectionProcess,3);

        TLegend* leg_ConversionSpeccanvas       = GetAndSetLegend(0.3,0.76,4);
        leg_ConversionSpeccanvas->SetTextSize(0.04);
        leg_ConversionSpeccanvas->AddEntry(histoCorrectedPi0Yield,"#pi^{0} corr. yield","pl");
        leg_ConversionSpeccanvas->AddEntry(fitPi0YieldA,Form("Hagedorn: %s #chi^{2}/ndf %.2f",fitPi0B.Data(),fitPi0YieldB->GetChisquare()/fitPi0YieldB->GetNDF()),"l");
        leg_ConversionSpeccanvas->AddEntry(fitPi0YieldB,Form("Levy-Tsallis: %s #chi^{2}/ndf %.2f",fitPi0B.Data(),fitPi0YieldB->GetChisquare()/fitPi0YieldB->GetNDF()),"l");
        leg_ConversionSpeccanvas->AddEntry(fitPi0YieldC,Form("Mod. Hagedorn: %s #chi^{2}/ndf %.2f",fitPi0C.Data(),fitPi0YieldC->GetChisquare()/fitPi0YieldC->GetNDF()),"l");
        leg_ConversionSpeccanvas->Draw();

   // Plotting ratios of data and fits in lower pad
        padConversionRatios->cd();
        padConversionRatios->SetLogx();

        histRatioConversionPi0A                 = (TH1D*) histoCorrectedPi0Yield->Clone("histRatioConversionPi0A");
        histRatioConversionPi0A                 = CalculateHistoRatioToFit(histRatioConversionPi0A,fitPi0YieldA);
        histRatioConversionPi0B                 = (TH1D*) histoCorrectedPi0Yield->Clone("histRatioConversionPi0B");
        histRatioConversionPi0B                 = CalculateHistoRatioToFit(histRatioConversionPi0B,fitPi0YieldB);
        histRatioConversionPi0C                 = (TH1D*) histoCorrectedPi0Yield->Clone("histRatioConversionPi0C");
        histRatioConversionPi0C                 = CalculateHistoRatioToFit(histRatioConversionPi0C,fitPi0YieldC);

        textSizeSpectra                         = 0.1;
        TH2F * histoDummy31                     = new TH2F("histoDummy31","histoDummy31",1000,ptMin,1.5*histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()),1000,0.41, 1.59);
        SetStyleHistoTH2ForGraphs(histoDummy31, "#it{p}_{T} (GeV/#it{c})", "data/fit",textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.95,0.6);
        histoDummy31->GetYaxis()->SetLabelOffset(0.01);
        histoDummy31->GetXaxis()->SetLabelOffset(-1e-2);
        histoDummy31->DrawCopy();

        DrawGammaSetMarker(histRatioConversionPi0A, 20, 1.4,kBlue-3,kBlue-3);
        DrawGammaSetMarker(histRatioConversionPi0B, 24, 1.4,kRed-3,kRed-3);
        DrawGammaSetMarker(histRatioConversionPi0C, 25, 1.4,kGreen-2,kGreen-2);

        DrawGammaLines(ptMin,1.5*histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()), 1.0, 1.0,1.0, kGray+2 ,7);
        histRatioConversionPi0A->Draw("e1,same");
        histRatioConversionPi0B->Draw("e1,same");
        histRatioConversionPi0C->Draw("e1,same");
        
   // Save plot
        ConversionSpeccanvas->Print(Form("%s/%s_Spectra_ConversionPi0_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
    }

//**********************************************************************************
//***                      Inclusive ratio with & w/o fit                        ***
//**********************************************************************************
    Bool_t doInclusiveFitRatio                  = kTRUE;
    if (doInclusiveFitRatio){
        histoIncRatioFitPurity                  = (TH1D*) histoIncRatioPurityTrueEff->Clone("histoIncRatioFitPurity");
        
   // Set datapoints in histoIncRatioFitPurity to the fit values from fitPi0YieldB (Tsallis fit)
   for(Int_t bin = 1; bin<histoIncRatioFitPurity->GetNbinsX()+1; bin++){
      histoIncRatioFitPurity->SetBinContent(bin,fitPi0YieldC->Eval(histoIncRatioPurityTrueEff->GetBinCenter(bin))); //nschmidt2016 changed to hagedorn
      histoIncRatioFitPurity->SetBinError(bin,histoCorrectedPi0Yield->GetBinError(bin));
   }
   DrawGammaSetMarker(histoIncRatioPurityTrueEff, 20, 2.0, 1, 1); 
   histoIncRatioFitPurity->Divide(histoGammaSpecCorrPurity,histoIncRatioFitPurity);

        TCanvas* canvasIncRatioFit              = GetAndSetCanvas("canvasIncRatioFit",0.095,0.09,1000,815);
        SetHistogramm(histoIncRatioFitPurity, "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
        DrawGammaSetMarker(histoIncRatioFitPurity, 20, 2.0, kBlue-3, kBlue-3);

        histoIncRatioPurityTrueEff->DrawCopy("e1");
        histoIncRatioFitPurity->DrawCopy("same,e1");

        PlotLatexLegend(0.6, 0.75, 0.045,collisionSystem,detectionProcess,3);

        TLegend* leg_canvasIncRatioFit          = GetAndSetLegend(0.6,0.65,2);
        leg_canvasIncRatioFit->AddEntry(histoIncRatioPurityTrueEff,"Inclusive Ratio");
        leg_canvasIncRatioFit->AddEntry(histoIncRatioFitPurity,"Ratio #pi^{0} Fit");
        leg_canvasIncRatioFit->Draw();

        if(option.CompareTo("PbPb_2.76TeV") == 0) DrawCentrality(0.8,0.9,centrality);
        canvasIncRatioFit->Print(Form("%s/%s_IncRatio_Fit_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
   }

//**********************************************************************************
//***                      Cocktail                                              ***
//**********************************************************************************
    Bool_t doCocktailLoading                    = kTRUE;
    if (doCocktailLoading){
        cocktailFile                            = new TFile(nameFileCocktail);

        cout<<"loading cocktail file: "<<nameFileCocktail<<endl;

        if (!option.CompareTo("7TeV") || !option.CompareTo("13TeV")){
            cocktailPi0                         = (TH1D* )cocktailFile->Get("Pi0_Pt");
            cocktailAllGamma                    = (TH1D* )cocktailFile->Get("Gamma_Pt");
            cocktailAllGammaPi0                 = (TH1D* )cocktailFile->Get("Gamma_From_Pi0_Pt");
        } else {
            cocktailPi0                         = (TH1D* )cocktailDir->Get("ptPion");
            cocktailAllGamma                    = (TH1D* )cocktailDir->Get("sumgammapi0");
            cocktailAllGammaPi0                 = (TH1D* )cocktailDir->Get("sumgammapi0");
        }
        if (cocktailAllGamma)      cout << "found cocktailGamma"    << endl;
        if (cocktailAllGammaPi0)   cout << "found cocktailGammaPi0" << endl;
        if (cocktailPi0)           cout << "found cocktailPi0"      << endl;
    }
    
//**********************************************************************************
//***                      Meson spectra data - cocktail                         ***
//**********************************************************************************
    Bool_t doMesonSpectra                       = kTRUE;
    if (doMesonSpectra){
        TCanvas* cocktailCanvasMesonSpec        = GetAndSetCanvas("cocktailCanvasMesonSpec", 0.12, 0.1, 1000 ,1350);
        DrawGammaCanvasSettings( cocktailCanvasMesonSpec, 0.16, 0.02, 0.015, 0.07);
        cocktailCanvasMesonSpec->SetLogy();
        cocktailCanvasMesonSpec->SetLogx();

        Double_t        minPt                   = 0.2;
        if (mode==4)    minPt                   = 1.5;
        
        cocktailPi0->GetXaxis()->SetLabelOffset(-1e-2);
        DrawAutoGammaMesonHistos( cocktailPi0,
                                 "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                 kTRUE, 20.,5e-10, kFALSE,
                                 kFALSE, -0.004, 0.020,
                                 kTRUE, minPt,(histoCorrectedPi0Yield->GetXaxis())->GetBinUpEdge(histoCorrectedPi0Yield->GetNbinsX()),62,0.04,62,0.03,0.7,1.7);
        DrawGammaSetMarker(cocktailPi0, 20, 2.0,colorCocktailPi0,colorCocktailPi0);
        
        //cocktailPi0->Draw("chist");
        cocktailPi0->DrawCopy("e1,same");
        
        fitPi0YieldB->SetLineColor(1);
        fitPi0YieldB->DrawCopy("same");

        DrawGammaSetMarker(histoCorrectedPi0Yield, 4, 2.0, 1, 1);
        histoCorrectedPi0Yield->Draw("e1,same");

        PlotLatexLegend(0.65, 0.8, 0.045,collisionSystem,detectionProcess,3);
        //drawLatex("#pi^{0} from ALICE Data", 4.0, 0.002, kBlack,0.035);
        drawLatex("#pi^{0} from ALICE Data", 1.7, 0.000000001, kBlack,0.035);
        
        TLegend* leg_MesonSpectra;
        leg_MesonSpectra                        = GetAndSetLegend(0.65,0.7,3);
        leg_MesonSpectra->AddEntry(cocktailPi0,"Cocktail #pi^{0}");
        leg_MesonSpectra->AddEntry(histoCorrectedPi0Yield,"#pi^{0} data");
        leg_MesonSpectra->AddEntry(fitPi0YieldB,"fit to #pi^{0} data","l");
        leg_MesonSpectra->Draw();

        cocktailCanvasMesonSpec->Print(Form("%s/%s_Cocktail_MesonSpectra_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
    }

//**********************************************************************************
//***                      NLO Direct Photon Ratio                               ***
//**********************************************************************************
    if (doNLOComparison){
        TF1* QGPin7Tev                          = new TF1("QGPin7Tev","1+[0]*exp(-x/[1])+[2]*exp(-x/[3])",0,100);
        QGPin7Tev->SetParameters(4.24707e+01, 2.85766e+00, 1.63353e+02, 4.53670e-01);

        cocktailAllGammaNLO                     = (TH1D*) cocktailAllGamma->Clone("cocktailAllGammaNLO");
        graphDirectPhotonNLOCopy                = (TGraphAsymmErrors*)graphDirectPhotonNLO->Clone("graphNLOCalcCopy");

        xVal                                    = graphDirectPhotonNLOCopy->GetX();
        xErr                                    = graphDirectPhotonNLOCopy->GetEX();
        yVal                                    = graphDirectPhotonNLOCopy->GetY();
        yErrUp                                  = graphDirectPhotonNLOCopy->GetEYhigh();
        yErrDown                                = graphDirectPhotonNLOCopy->GetEYlow();

        TString cocktailFit                     = "xqcd";
        TString fitOptions                      = "QNRME+";//"IQNRME+";
        if (!option.CompareTo("7TeV") && mode == 4) {
            fitOptions                          = "QNRME+";
        }

        cocktailFitAllGammaForNLO               = (TF1*)FitObject(cocktailFit,"cocktailFit","Pi0",cocktailAllGammaNLO,2.0,16,NULL,fitOptions);

        for (Int_t bin=0; bin<graphDirectPhotonNLOCopy->GetN(); bin++) {
            yVal[bin]                           = (1 + ( yVal[bin] / (cocktailFitAllGammaForNLO->Eval(xVal[bin]))));
        }
        
        // ------------------------------ NLO Calculations -----------------------------
        NLODoubleRatio                          = new TGraphAsymmErrors(graphDirectPhotonNLOCopy->GetN(), xVal, yVal, xErr, xErr, yErrDown, yErrUp);
        NLODoubleRatio->SetName("NLODoubleRatio");
        NLO                                     = (TGraphAsymmErrors*)graphDirectPhotonNLO->Clone("NLO");
    }


//**********************************************************************************
//***                      Double Ratios                                         ***
//**********************************************************************************
   Bool_t doDoubleRatios                        = kTRUE;
   if (doDoubleRatios){
       cocktailAllGammaPi0                      = RebinTH1D(cocktailAllGammaPi0,histoIncRatioPurityTrueEff);
       cocktailAllGammaRebinned                 = RebinTH1D(cocktailAllGamma,histoIncRatioPurityTrueEff);
       cocktailPi0Rebinned                      = RebinTH1D(cocktailPi0,histoIncRatioPurityTrueEff);
       cocktailAllGammaPi0->Divide(cocktailAllGammaRebinned,cocktailPi0Rebinned);

       histoDoubleRatioConversionTrueEffPurity  = (TH1D*) histoIncRatioPurityTrueEff->Clone("DoubleRatioConversionTrueEffPurity");
       histoDoubleRatioConversionTrueEffPurity->Divide(cocktailAllGammaPi0);
       if (doSysErr)    graphDoubleRatioSysErr  = CalculateSysErrFromRelSysHisto( histoDoubleRatioConversionTrueEffPurity, "DoubleRatioSystError",relSystErrorDoubleRatioDown , relSystErrorDoubleRatioUp, 2, nPointsDoubleRatio);
       else             graphDoubleRatioSysErr  = NULL;
           
       
       histoDoubleRatioFitPi0YieldPurity        = (TH1D*) histoIncRatioFitPurity->Clone("DoubleRatioConversionFitPurity");
       histoDoubleRatioFitPi0YieldPurity->Divide(cocktailAllGammaPi0);
       if (doSysErr)graphDoubleRatioFitSysErr   = CalculateSysErrFromRelSysHisto( histoDoubleRatioFitPi0YieldPurity, "DoubleRatioFitError",relSystErrorDoubleRatioDown , relSystErrorDoubleRatioUp, 2, nPointsDoubleRatio);
       else         graphDoubleRatioFitSysErr   = NULL;

       // double ratio combined
       TCanvas *canvasConversionFitDoubleRatioSum = new TCanvas("canvasConversionFitDoubleRatioSum","",0.095,0.09,1000,815);
       DrawGammaCanvasSettings( canvasConversionFitDoubleRatioSum, 0.09, 0.02, 0.02, 0.09);
       canvasConversionFitDoubleRatioSum->cd();
       canvasConversionFitDoubleRatioSum->SetLogx();

       Double_t minPt                           = 0.3;
       Double_t maxPt                           = 40.;
       Double_t minY                            = 0.85;
       Double_t maxY                            = 1.65;
       if (mode==4) {
           minPt                                = 1.5;
           minY                                 = 0.5;
       }

       TString      method                      = "conversions";
       if (mode==4) method                      = "EMCal";

       TH2F * histo2DDoubleRatioPlotting       = new TH2F("histo2DDoubleRatioPlotting","histo2DDoubleRatioPlotting",1000,minPt,maxPt,1000,minY,maxY);
       SetStyleHistoTH2ForGraphs(histo2DDoubleRatioPlotting, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.04,0.04, 0.04,0.04, 1.,1.);
       histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.,histoDoubleRatioConversionTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioConversionTrueEffPurity->GetNbinsX())*1.5);
       histo2DDoubleRatioPlotting->GetXaxis()->SetTitleFont(62);
       histo2DDoubleRatioPlotting->GetYaxis()->SetTitleFont(62);
       histo2DDoubleRatioPlotting->GetXaxis()->SetLabelOffset(-1e-2);
       histo2DDoubleRatioPlotting->DrawCopy();

       DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurity, 20, 2.0, kBlack, kBlack);
       DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurity, 20, 2.0, kBlue+2, kBlue+2);
       if (graphDoubleRatioSysErr)      DrawGammaSetMarkerTGraphAsym(graphDoubleRatioSysErr,26,0,kGray+2,kGray+2,2,kTRUE);
       if (graphDoubleRatioFitSysErr)   DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitSysErr,26,0,kBlue-8,kBlue-8,2,kTRUE);

       histoDoubleRatioFitPi0YieldPurity->DrawCopy("same,E1");
       DrawGammaLines(0., histoDoubleRatioConversionTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioConversionTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
       histoDoubleRatioConversionTrueEffPurity->DrawCopy("same,E1");
       if (graphDoubleRatioSysErr)      graphDoubleRatioSysErr->Draw("p,2,same");
       if (graphDoubleRatioFitSysErr)   graphDoubleRatioFitSysErr->Draw("p,2,same");

       TLegend* legendDoubleConversionFit       = GetAndSetLegend(0.14,0.76,3,1,Form("direct photon signal via %s", method.Data()));
       legendDoubleConversionFit->AddEntry(histoDoubleRatioConversionTrueEffPurity,"measured direct photon signal","p");
       legendDoubleConversionFit->AddEntry(histoDoubleRatioFitPi0YieldPurity,"Measured Direct Photon Signal fitted #pi^{0}","p");
       legendDoubleConversionFit->Draw();

       canvasConversionFitDoubleRatioSum->Print(Form("%s/%s_DoubleRatioComparison_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
       
       // double ratio fit
       TCanvas *canvasConversionFitDoubleRatioSum2 = new TCanvas("canvasConversionFitDoubleRatioSum2","",0.095,0.09,1000,815);
       DrawGammaCanvasSettings( canvasConversionFitDoubleRatioSum2, 0.09, 0.02, 0.02, 0.09);
       canvasConversionFitDoubleRatioSum2->cd();
       canvasConversionFitDoubleRatioSum2->SetLogx();
       
       histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.2, 20);
       histo2DDoubleRatioPlotting->DrawCopy();
       
       DrawGammaLines(0., 20. , 1.0, 1.0,2.0, kGray+2 ,7);
       histoDoubleRatioFitPi0YieldPurity->DrawCopy("same,E1");
       if (graphDoubleRatioFitSysErr)   graphDoubleRatioFitSysErr->Draw("p,2,same");
       
       TLegend* legendDoubleConversionFit2       = GetAndSetLegend(0.14,0.8,2,1,Form("direct photon signal via %s", method.Data()));
       legendDoubleConversionFit2->AddEntry(histoDoubleRatioFitPi0YieldPurity,"Measured Direct Photon Signal fitted #pi^{0}","p");
       legendDoubleConversionFit2->Draw();
       
       canvasConversionFitDoubleRatioSum2->Print(Form("%s/%s_DoubleRatioFit_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

       // double ratio no fit
       TCanvas *canvasConversionFitDoubleRatioSum3 = new TCanvas("canvasConversionFitDoubleRatioSum","",0.095,0.09,1000,815);
       DrawGammaCanvasSettings( canvasConversionFitDoubleRatioSum3, 0.09, 0.02, 0.02, 0.09);
       canvasConversionFitDoubleRatioSum3->cd();
       canvasConversionFitDoubleRatioSum3->SetLogx();

       histo2DDoubleRatioPlotting->DrawCopy();
       
       DrawGammaLines(0., histoDoubleRatioConversionTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioConversionTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
       histoDoubleRatioConversionTrueEffPurity->DrawCopy("same,E1");
       if (graphDoubleRatioSysErr) graphDoubleRatioSysErr->Draw("p,2,same");
       
       TLegend* legendDoubleConversionFit3       = GetAndSetLegend(0.14,0.8,2,1,Form("direct photon signal via %s", method.Data()));
       legendDoubleConversionFit3->AddEntry(histoDoubleRatioConversionTrueEffPurity,"measured direct photon signal","p");
       legendDoubleConversionFit3->Draw();
       
       canvasConversionFitDoubleRatioSum3->Print(Form("%s/%s_DoubleRatio_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
       
       // do theory comparisons
       if (doNLOComparison){
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
         
           NLOdoubleRatioFit                   = new TF1("NLOdoubleRatioFit","(x<=2)*(1.0+[0]*x)+(x>2)*([1]+[2]*x+[3]*x*x)",0.,4.0);
           NLODoubleRatio->Fit(NLOdoubleRatioFit,"NRME+","",0,4.0);
           NLOdoubleRatioFit->SetLineColor(2);
           
           // ------------------------------ NLO Calculations -----------------------------
           // combined
           TCanvas *canvasConversionFitDoubleRatioNLO = new TCanvas("canvasConversionFitDoubleRatioNLO","",0.095,0.09,1000,815);
           DrawGammaCanvasSettings( canvasConversionFitDoubleRatioNLO, 0.09, 0.02, 0.02, 0.09);
           canvasConversionFitDoubleRatioNLO->cd();
           canvasConversionFitDoubleRatioNLO->SetLogx();

           histo2DDoubleRatioPlotting->DrawCopy();
           NLODoubleRatio->RemovePoint(0);
           DrawGammaLines(0., histoDoubleRatioConversionTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioConversionTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
           NLODoubleRatio->Draw("lp3");

           DrawGammaSetMarker(histoDoubleRatioConversionTrueEffPurity, 20, 2.0, kBlack, kBlack);
           DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurity, 20, 2.0, kBlue+2, kBlue+2);
           histoDoubleRatioConversionTrueEffPurity->DrawCopy("same");
           histoDoubleRatioFitPi0YieldPurity->DrawCopy("same");
           if (graphDoubleRatioSysErr)     graphDoubleRatioSysErr->Draw("p,2,same");
           if (graphDoubleRatioFitSysErr)  graphDoubleRatioFitSysErr->Draw("p,2,same");
        
           TLegend* legendDoubleConversionFitNLO = GetAndSetLegend(0.14,0.72,4,1,Form("direct photon signal via %s", method.Data()));
           legendDoubleConversionFitNLO->AddEntry(histoDoubleRatioConversionTrueEffPurity,"measured direct photon signal","p");
           legendDoubleConversionFitNLO->AddEntry(histoDoubleRatioFitPi0YieldPurity,"measured direct photon signal, fitted #pi^{0}","p");
           if(option.CompareTo("PbPb_2.76TeV") == 0)       legendDoubleConversionFitNLO->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 2.76TeV scaled N_{coll}","l");
           else if(option.CompareTo("pPb_5.023TeV") == 0)  legendDoubleConversionFitNLO->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 5.023TeV scaled N_{coll}","l");
           else                                            legendDoubleConversionFitNLO->AddEntry(NLODoubleRatio,"pp NLO Direct Photon","l");
           legendDoubleConversionFitNLO->Draw();

           canvasConversionFitDoubleRatioNLO->Print(Form("%s/%s_DoubleRatioComparison_NLO_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));

           // fit
           TCanvas *canvasConversionFitDoubleRatioNLO2 = new TCanvas("canvasConversionFitDoubleRatioNLO2","",0.095,0.09,1000,815);
           DrawGammaCanvasSettings( canvasConversionFitDoubleRatioNLO2, 0.09, 0.02, 0.02, 0.09);
           canvasConversionFitDoubleRatioNLO2->cd();
           canvasConversionFitDoubleRatioNLO2->SetLogx();
           
           histo2DDoubleRatioPlotting->DrawCopy();
           NLODoubleRatio->RemovePoint(0);
           DrawGammaLines(0., histoDoubleRatioConversionTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioConversionTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
           NLODoubleRatio->Draw("lp3");
           
           histoDoubleRatioFitPi0YieldPurity->DrawCopy("same");
           if (graphDoubleRatioFitSysErr) graphDoubleRatioFitSysErr->Draw("p,2,same");
           
           TLegend* legendDoubleConversionFitNLO2 = GetAndSetLegend(0.14,0.76,3,1,Form("direct photon signal via %s", method.Data()));
           legendDoubleConversionFitNLO2->AddEntry(histoDoubleRatioFitPi0YieldPurity,"measured direct photon signal, fitted #pi^{0}","p");
           if(option.CompareTo("PbPb_2.76TeV") == 0)       legendDoubleConversionFitNLO2->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 2.76TeV scaled N_{coll}","l");
           else if(option.CompareTo("pPb_5.023TeV") == 0)  legendDoubleConversionFitNLO2->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 5.023TeV scaled N_{coll}","l");
           else                                            legendDoubleConversionFitNLO2->AddEntry(NLODoubleRatio,"pp NLO Direct Photon","l");
           legendDoubleConversionFitNLO2->Draw();
           
           canvasConversionFitDoubleRatioNLO2->Print(Form("%s/%s_DoubleRatioFit_NLO_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
           
           // measured
           TCanvas *canvasConversionFitDoubleRatioNLO3 = new TCanvas("canvasConversionFitDoubleRatioNLO","",0.095,0.09,1000,815);
           DrawGammaCanvasSettings( canvasConversionFitDoubleRatioNLO3, 0.09, 0.02, 0.02, 0.09);
           canvasConversionFitDoubleRatioNLO3->cd();
           canvasConversionFitDoubleRatioNLO3->SetLogx();
           
           histo2DDoubleRatioPlotting->DrawCopy();
           NLODoubleRatio->RemovePoint(0);
           DrawGammaLines(0., histoDoubleRatioConversionTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioConversionTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
           NLODoubleRatio->Draw("lp3");
           
           histoDoubleRatioConversionTrueEffPurity->DrawCopy("same");
           if (graphDoubleRatioSysErr) graphDoubleRatioSysErr->Draw("p,2,same");
           
           TLegend* legendDoubleConversionFitNLO3 = GetAndSetLegend(0.14,0.76,3,1,Form("direct photon signal via %s", method.Data()));
           legendDoubleConversionFitNLO3->AddEntry(histoDoubleRatioConversionTrueEffPurity,"measured direct photon signal","p");
           if(option.CompareTo("PbPb_2.76TeV") == 0)       legendDoubleConversionFitNLO3->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 2.76TeV scaled N_{coll}","l");
           else if(option.CompareTo("pPb_5.023TeV") == 0)  legendDoubleConversionFitNLO3->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 5.023TeV scaled N_{coll}","l");
           else                                            legendDoubleConversionFitNLO3->AddEntry(NLODoubleRatio,"pp NLO Direct Photon","l");
           legendDoubleConversionFitNLO3->Draw();
           
           canvasConversionFitDoubleRatioNLO3->Print(Form("%s/%s_DoubleRatio_NLO_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
           
           // --------------------- Direct Photon Spectrum -------------------------
           if (graphDoubleRatioSysErr) {
               // calculate direct photon spectrum
               histoDirectPhotonSpectrum                = (TH1D*)histoGammaSpecCorrPurity->Clone("histoDirectPhotonSpectrum");
               histoThermalPhotonSpectrum               = (TH1D*)histoGammaSpecCorrPurity->Clone("histoThermalPhotonSpectrum");
               histoPromptPhotonSpectrum                = (TH1D*)histoGammaSpecCorrPurity->Clone("histoPromptPhotonSpectrum");
               
               Bool_t doUpperLimits                     = kTRUE;
               if (!doUpperLimits) {
                   for(Int_t i = 1; i < histoDirectPhotonSpectrum->GetNbinsX(); i++){
                       histoDirectPhotonSpectrum->SetBinContent(    i+1,0);
                       histoThermalPhotonSpectrum->SetBinContent(   i+1,0);
                   }
                   for(Int_t i = 1; i < histoDirectPhotonSpectrum->GetNbinsX(); i++){
                       Double_t binContent              = -1;
                       Double_t binError                = -1;
                       binContent                       = 1-(1./( histoDoubleRatioConversionTrueEffPurity->GetBinContent(i+1)));
                       binContent                       = binContent*histoGammaSpecCorrPurity->GetBinContent(i+1);
                       binError                         = 1-(1./( histoDoubleRatioConversionTrueEffPurity->GetBinContent(i+1) + histoDoubleRatioConversionTrueEffPurity->GetBinError(i+1)));
                       binError                         = binError*histoGammaSpecCorrPurity->GetBinContent(i+1);
                       binError                         = binError-binContent;
                       histoDirectPhotonSpectrum->SetBinContent(i+1,binContent);
                       histoDirectPhotonSpectrum->SetBinError(  i+1,binError);
                       
                       if(histoDirectPhotonSpectrum->GetBinCenter(i+1) < 4.5){
                           binContent                   = 1-(1./(1+ histoDoubleRatioConversionTrueEffPurity->GetBinContent(i+1) - NLOdoubleRatioFit->Eval(histoDoubleRatioConversionTrueEffPurity->GetBinCenter(i+1))));
                           binContent                   = binContent*histoGammaSpecCorrPurity->GetBinContent(i+1);
                           histoThermalPhotonSpectrum->SetBinContent(   i+1,binContent);
                           histoThermalPhotonSpectrum->SetBinError(     i+1,binError);
                       } else {
                           histoThermalPhotonSpectrum->SetBinContent(   i+1,0);
                           histoThermalPhotonSpectrum->SetBinError(     i+1,0);
                       }
                       if(histoDirectPhotonSpectrum->GetBinCenter(i+1) < 4.5){
                           binContent                   = 1-(1./(NLOdoubleRatioFit->Eval(histoDoubleRatioConversionTrueEffPurity->GetBinCenter(i+1))));
                           binContent                   = binContent*histoGammaSpecCorrPurity->GetBinContent(i+1);
                           histoPromptPhotonSpectrum->SetBinContent(    i+1,binContent);
                           histoPromptPhotonSpectrum->SetBinError(      i+1,0.0000001);
                       } else {
                           histoPromptPhotonSpectrum->SetBinContent(    i+1,0);
                           histoPromptPhotonSpectrum->SetBinError(      i+1,0);
                       }
                   }
               } else {
                   histoDoubleRatioUpperLimits          = GetUpperLimitsHisto(histoDoubleRatioConversionTrueEffPurity,  graphDoubleRatioSysErr, 0.95, 0.004, 1e2);
                   Double_t binContent                  = -1;
                   Double_t binError                    = -1;
                   for (Int_t i=1; i<histoDirectPhotonSpectrum->GetNbinsX()+1; i++) {
                       if (!histoDirectPhotonSpectrum->GetBinContent(i)) continue;
                       binContent                       = histoDirectPhotonSpectrum->GetBinContent(i) * (1 - 1/histoDoubleRatioUpperLimits->GetBinContent(i));
                       binError                         = 0.;
                       histoDirectPhotonSpectrum->SetBinContent(i, binContent);
                       histoDirectPhotonSpectrum->SetBinError(  i, binError);
                   }
               }
               
               // plot direct photon spectrum
               TCanvas* CanvasDirGamma              = GetAndSetCanvas("CanvasDirGamma", 0.12, 0.1, 1000 ,1350);
               DrawGammaCanvasSettings( CanvasDirGamma, 0.16, 0.02, 0.015, 0.07);
               CanvasDirGamma->SetLogy();
               CanvasDirGamma->SetLogx();
               
               histoDirectPhotonSpectrum->SetLineColor(kBlack);
               histoDirectPhotonSpectrum->SetLineWidth(2);
               
               TH2F * dummyDir       = new TH2F("dummyDir","dummyDir",1000,0.2,20.0,1000,5e-10,2);
               SetStyleHistoTH2ForGraphs(dummyDir, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.04,0.04, 0.04,0.04, 1.,1.);
               dummyDir->GetXaxis()->SetRangeUser(0.,histoDoubleRatioConversionTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioConversionTrueEffPurity->GetNbinsX())*1.5);
               dummyDir->GetXaxis()->SetTitleFont(62);
               dummyDir->GetYaxis()->SetTitleFont(62);
               dummyDir->GetXaxis()->SetLabelOffset(-2e-2);
               dummyDir->GetXaxis()->SetTitleOffset(0.8);
               dummyDir->GetYaxis()->SetTitleOffset(1.5);
               dummyDir->DrawCopy();
               
               PlotLatexLegend(0.66, 0.78, 0.045,collisionSystem,detectionProcess,2);
               drawLatex("#gamma's from ALICE Data", 1.7, 0.000000001, kBlack,0.035);
               TArrow *ar2[50];
               for(Int_t i=1;i<histoDirectPhotonSpectrum->GetNbinsX()+1;i++)
               {
                   if (!histoDirectPhotonSpectrum->GetBinContent(i)) continue;
                   ar2[i] = new TArrow(histoDirectPhotonSpectrum->GetBinCenter(i),histoDirectPhotonSpectrum->GetBinContent(i),histoDirectPhotonSpectrum->GetBinCenter(i),histoDirectPhotonSpectrum->GetBinContent(i)*0.1,0.02,"|->");
                   ar2[i]->SetAngle(40);
                   ar2[i]->SetLineWidth(2);
                   ar2[i]->Draw();
               }
              
               NLO->GetXaxis()->SetRangeUser(NLO->GetXaxis()->GetXmin(), histoDirectPhotonSpectrum->GetXaxis()->GetXmax()); 
               NLO->SetLineWidth(2);
               NLO->Draw("lp3");
               
               TLegend* leg_GammaSpectra;
               leg_GammaSpectra                            = GetAndSetLegend(0.2,0.2,2);
               leg_GammaSpectra->AddEntry(histoDirectPhotonSpectrum,"direct photon upper limits",   "l");
               leg_GammaSpectra->AddEntry(NLO,                      "pp NLO direct photon",         "l");
               leg_GammaSpectra->Draw();
               
               CanvasDirGamma->Print(Form("%s/%s_DirectPhotonSpectrum_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));
           } else {
               cout << "No systematic uncertainties on double ratio given, skipping direct photon calculation." << endl;
           }
        }
    }

    fileFinalResults.close();
    
//**********************************************************************************
//***                 Inclusive and decay ratio                                  ***
//**********************************************************************************
    TCanvas* canvasIncRatioDecRatio         = GetAndSetCanvas("canvasIncRatioDecRatio",0.095,0.09,1000,815);
    canvasIncRatioDecRatio->SetLogx();

    SetHistogramm(cocktailAllGammaPi0,          "#it{p}_{T} (GeV/#it{c})", "#gamma/#pi^{0}",0,1.8);
    SetHistogramm(histoIncRatioPurityTrueEff,   "#it{p}_{T} (GeV/#it{c})", "#gamma/#pi^{0}",0,1.8);

    DrawGammaSetMarker(histoIncRatioPurityTrueEff,  20, 2.0, kBlack,  kBlack);
    DrawGammaSetMarker(cocktailAllGammaPi0,         33, 2.0, kBlue-2,  kBlue-2);

    histoIncRatioPurityTrueEff->GetXaxis()->SetRangeUser(0.2, 20);
    histoIncRatioPurityTrueEff->Draw();
    cocktailAllGammaPi0->Draw("same");

    PlotLatexLegend(0.66, 0.75, 0.045,collisionSystem,detectionProcess);

    TLegend* legendIncRatioDecRatio         = GetAndSetLegend(0.65,0.6,3);
    legendIncRatioDecRatio->AddEntry(histoIncRatioPurityTrueEff,    "(#gamma_{incl.} / #pi^{0})_{meas.}", "lp");
    legendIncRatioDecRatio->AddEntry(cocktailAllGammaPi0,           "(#gamma_{dec.} / #pi^{0})_{sim.}",  "lp");
    legendIncRatioDecRatio->Draw();

    canvasIncRatioDecRatio->Print(Form("%s/%s_IncRatio_DecRatio_trueEff_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data(),suffix.Data()));


//**********************************************************************************
//***                      Save Histograms                                       ***
//**********************************************************************************   
    Bool_t doSaving                             = kTRUE;
    if (doSaving){
        TString nameOutputFile                  = Form("%s/%s/%s_%s_GammaConvV1_InclusiveRatio_%s.root",cutSel.Data(),option.Data(),textPi0New.Data(),textPrefix2.Data(),centralityAdd.Data());
        fileCorrectedOutput                     = new TFile(nameOutputFile,"RECREATE");

   // pi0 quantities
        fitPi0YieldA->Write(            fitPi0YieldA->GetName(),            TObject::kOverwrite);
        fitPi0YieldB->Write(            fitPi0YieldB->GetName(),            TObject::kOverwrite);
        fitPi0YieldC->Write(            fitPi0YieldC->GetName(),            TObject::kOverwrite);
        histoCorrectedPi0Yield->Write(  histoCorrectedPi0Yield->GetName(),  TObject::kOverwrite);

   // gamma quantities
        if (ConversionGammaFitA)        ConversionGammaFitA->Write(         ConversionGammaFitA->GetName(),         TObject::kOverwrite);
        if (histoDirectPhotonSpectrum)  histoDirectPhotonSpectrum->Write(   histoDirectPhotonSpectrum->GetName(),   TObject::kOverwrite);
        if (histoGammaSpecCorrPurity)   histoGammaSpecCorrPurity->Write(    "histoGammaSpecCorrPurity",             TObject::kOverwrite);

   // Double ratio
        if(histoDoubleRatioFitPi0YieldPurity)           histoDoubleRatioFitPi0YieldPurity->Write(           histoDoubleRatioFitPi0YieldPurity->GetName(),               TObject::kOverwrite);
        if(histoDoubleRatioConversionHighFitPurity)     histoDoubleRatioConversionHighFitPurity->Write(     histoDoubleRatioConversionHighFitPurity->GetName(),         TObject::kOverwrite);
        if(histoDoubleRatioConversionLowFitPurity)      histoDoubleRatioConversionLowFitPurity->Write(      histoDoubleRatioConversionLowFitPurity->GetName(),          TObject::kOverwrite);
        if(histoDoubleRatioFitPi0YieldPurityModA)       histoDoubleRatioFitPi0YieldPurityModA->Write(       histoDoubleRatioFitPi0YieldPurityModA->GetName(),           TObject::kOverwrite);
        if(histoDoubleRatioFitPi0YieldPurityModB)       histoDoubleRatioFitPi0YieldPurityModB->Write(       histoDoubleRatioFitPi0YieldPurityModB->GetName(),           TObject::kOverwrite);
        if(histoDoubleRatioConversionTrueEffPurity)     histoDoubleRatioConversionTrueEffPurity->Write(     histoDoubleRatioConversionTrueEffPurity->GetName(),         TObject::kOverwrite);
        if(histoDoubleRatioConversionTrueEffPurityModA) histoDoubleRatioConversionTrueEffPurityModA->Write( histoDoubleRatioConversionTrueEffPurityModA->GetName(),     TObject::kOverwrite);
        if(histoDoubleRatioConversionTrueEffPurityModB) histoDoubleRatioConversionTrueEffPurityModB->Write( histoDoubleRatioConversionTrueEffPurityModB->GetName(),     TObject::kOverwrite);
        if(histoDoubleRatioUpperLimits)                 histoDoubleRatioUpperLimits->Write(Form("%s_UpperLimits", histoDoubleRatioConversionTrueEffPurity->GetName()),  TObject::kOverwrite);

   // systematics
        if (graphGammaYieldSysErr)      graphGammaYieldSysErr->Write(       Form("%s_SystErr", histoGammaSpecCorrPurity->GetName()),                TObject::kOverwrite);
        if (graphInclRatioSysErr)       graphInclRatioSysErr->Write(        Form("%s_SystErr", histoIncRatioPurityTrueEff->GetName()),              TObject::kOverwrite);
        if (graphDoubleRatioSysErr)     graphDoubleRatioSysErr->Write(      Form("%s_SystErr", histoDoubleRatioConversionTrueEffPurity->GetName()), TObject::kOverwrite);
        if (graphDoubleRatioFitSysErr)  graphDoubleRatioFitSysErr->Write(   Form("%s_SystErr", histoDoubleRatioFitPi0YieldPurity->GetName()),       TObject::kOverwrite);
        
   // inclusive ratio
        histoIncRatioPurityTrueEff->Write(  histoIncRatioPurityTrueEff->GetName(),  TObject::kOverwrite);
        histoMCIncRatio->Write(             histoMCIncRatio->GetName(),             TObject::kOverwrite);
        histoIncRatioFitPurity->Write(      histoIncRatioFitPurity->GetName(),      TObject::kOverwrite);

   // cocktail
        cocktailAllGamma->Write(cocktailAllGamma->GetName(),    TObject::kOverwrite);
        cocktailPi0->Write(     cocktailPi0->GetName(),         TObject::kOverwrite);
        
   //NLO
        fitNLODirectPhoton->Write(          "fitNLODirectPhoton",           TObject::kOverwrite);
        graphPromptPhotonNLO->Write(        "graphPromptPhotonNLO",         TObject::kOverwrite);
        fitNLOPromptPhoton->Write(          "fitNLOPromptPhoton",           TObject::kOverwrite);
        graphFragmentationPhotonNLO->Write( "graphFragmentationPhotonNLO",  TObject::kOverwrite);
        fitNLOFragmentationPhoton->Write(   "fitNLOFragmentationPhoton",    TObject::kOverwrite);

        histRatioNLODirectPhoton->Write(    "histRatioNLODirectPhoton",     TObject::kOverwrite);
        NLODoubleRatio->Write(              "graphNLODoubleRatio",          TObject::kOverwrite);
        NLO->Write(                         "graphNLODirGamma",             TObject::kOverwrite);
        
        fileCorrectedOutput->Close();
    }
}