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

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
Double_t FindLargestBin1DHist(TH1* hist ){
    Double_t largestContent     = 0;
    for (Int_t i= 0; i < hist->GetNbinsX(); i++){
        if (largestContent < hist->GetBinContent(i)){
            largestContent = hist->GetBinContent(i);
        }    
    }    
    return largestContent;
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
Double_t FindSmallestBin1DHist(TH1* hist, Double_t maxStart = 1e6 ){
    Double_t smallesContent     = maxStart;
    for (Int_t i= 0; i < hist->GetNbinsX(); i++){
        if (hist->GetBinContent(i) != 0 && smallesContent > hist->GetBinContent(i)){
            smallesContent = hist->GetBinContent(i);
        }    
    }    
    return smallesContent;
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
void  CalculateGammaToPi0V3(    TString nameFileGamma   = "",
                                TString nameFilePi0     = "",
                                TString nameFileCocktail= "",
                                TString cutSel          = "",
                                TString suffix          = "pdf",
                                TString nameMeson       = "",
                                TString isMC            = "",
                                TString fEnergy          = "",
                                TString fEstimatePileup = "",
                                Int_t mode              = 0
                            ){
    // switch systematics on/off
    Bool_t doSysErr                             = kFALSE;
    if (!fEnergy.CompareTo("7TeV") && mode == 0)
        doSysErr                                = kTRUE;
    
    // Setting the general style
    gROOT->Reset();
    StyleSettingsThesis();
    SetPlotStyle();

    // Separating cutnumber and retrieving centrality and number of events
    SeparateCutnumberString(cutSel,mode);

    // Create strings for naming
    CreateNamingStrings(nameMeson,isMC);
    TString collisionSystem                      = ReturnFullCollisionsSystem(fEnergy);
    TString detectionProcess                     = ReturnFullTextReconstructionProcess(mode);
    if(fEnergy.Contains("PbPb")){
        collisionSystem                          = Form("%s %s", centrality.Data(), collisionSystem.Data());
    } else if(fEnergy.Contains("pPb") && !centrality.Contains("0-100")){
        collisionSystem                          = Form("%s %s", centrality.Data(), collisionSystem.Data());
    }

    // Creating output directory and output file
    cout << "Output directory with plots:" << endl;
    cout << Form("%s/%s/%s/GammaToPi0",cutSel.Data(),fEnergy.Data(),suffix.Data()) << endl;
    TString outputDir                            = Form("%s/%s/%s/GammaToPi0",cutSel.Data(),fEnergy.Data(),suffix.Data());
    gSystem->Exec("mkdir "+outputDir);

    TString nameFinalResDat                      = Form("%s/%s/Gamma_%s_FinalExtraction_%s.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data(), cutSel.Data());
    fstream fileFinalResults;
    fileFinalResults.open(nameFinalResDat.Data(), ios::out);

    // Opening gamma file and loading data gamma spectrum
    fileGamma                                   = new TFile(nameFileGamma);
    histoGammaSpecCorrPurity                    = (TH1D*)fileGamma->Get("GammaCorrUnfold_Pt");
    if (!histoGammaSpecCorrPurity) {
        cout << "ERROR: GammaCorrUnfold_Pt not in gamma file" << endl;
        return;
    }    
    histoMCDecaySumGammaPt                      = (TH1D*)fileGamma->Get("MC_DecayGammaAll_Pt");
    histoGammaSpecMCAll                         = (TH1D*)fileGamma->Get("GammaSpecMC");

    // Opening pion file and loading spectra
    filePi0                                     = new TFile(nameFilePi0);
    histoCorrectedPi0Yield                      = (TH1D*)filePi0->Get("CorrectedYieldTrueEff");
    histoMCYieldMeson                           = (TH1D*)filePi0->Get("MCYield_Meson");
    histoMCYieldMesonOldBin                     = (TH1D*)filePi0->Get("MCYield_Meson_oldBin");

    Bool_t haveEta                              = kFALSE;
    TString nameFileEta                         = nameFilePi0.ReplaceAll("Pi0","Eta");
    fileEta                                     = new TFile(nameFileEta);
    if (!fileEta->IsZombie()){
        histoCorrectedEtaYield                  = (TH1D*)fileEta->Get("CorrectedYieldTrueEff");
        if (histoCorrectedEtaYield)
            haveEta                             = kTRUE;
    }
    
    // Creating inclusive gamma ratio
    histoIncRatioPurityTrueEff                  = (TH1D*) histoGammaSpecCorrPurity->Clone("IncRatioPurity_trueEff");
    histoIncRatioPurityTrueEff->Divide(histoIncRatioPurityTrueEff,histoCorrectedPi0Yield,1,1,"");

    histoMCIncRatio                             = (TH1D*) histoGammaSpecMCAll->Clone("MC_IncRatio");
    histoMCIncRatio->Divide(histoGammaSpecMCAll,histoMCYieldMeson,1,1,"");

    // Opening systematics file
    TString fileNameSysErrGamma                  =""; // default
    TString fileNameSysErrInclRatio              =""; // default
    TString fileNameSysErrDoubleRatio            =""; // default
    if(fEnergy.CompareTo("7TeV") == 0){
        fileNameSysErrGamma                     = "GammaSystematicErrorsCalculated/SystematicErrorAveraged_Gamma_7TeV_2016_12_15.dat";
        fileNameSysErrInclRatio                 = "GammaSystematicErrorsCalculated/SystematicErrorAveraged_IncRatio_7TeV_2016_12_15.dat";
        fileNameSysErrDoubleRatio               = "GammaSystematicErrorsCalculated/SystematicErrorAveraged_DoubleRatio_7TeV_2016_12_15.dat";
    }

    if (doSysErr){
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
        
        graphGammaYieldSysErr   = CalculateSysErrFromRelSysHisto( histoGammaSpecCorrPurity, "Pi0SystError",relSystErrorGammaDown , relSystErrorGammaUp, 2, nPointsGamma);
        graphInclRatioSysErr    = CalculateSysErrFromRelSysHisto( histoIncRatioPurityTrueEff, "Pi0SystErrorA",relSystErrorInclRatioDown , relSystErrorInclRatioUp, 2, nPointsInclRatio);
    } else {
        graphGammaYieldSysErr   = NULL;
        graphInclRatioSysErr    = NULL;
    }
   
    //**********************************************************************************
    //***                  Load NLO Calculatins Ratio                                 ***
    //**********************************************************************************
    Bool_t doNLOComparison                      = kTRUE;
    if (doNLOComparison){
        fileTheoryCompilation                   = new TFile("ExternalInput/Theory/TheoryCompilationPP.root");
        directoryGamma                          = (TDirectoryFile*)fileTheoryCompilation->Get("DirectPhoton");
        // Load NLO input from TheoryCompilationPP.root
        if (!fEnergy.CompareTo("900GeV")) {
            cout << "WARNING: USING 7TEV NLO CALCULATIONS FOR 900GeV!" << endl;
            graphDirectPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT7_7TeV");
            graphPromptPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT7_7TeV");
            graphFragmentationPhotonNLO         = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT7_7TeV");
        } else if (!fEnergy.CompareTo("2.76TeV") || !fEnergy.CompareTo("PbPb_2.76TeV")) {
            if (!triggerCutNumber.CompareTo("1") && !subTriggerCutNumber.CompareTo("0")) {   // kINT7
                graphDirectPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT7_2760GeV");
                graphPromptPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT7_2760GeV");
                graphFragmentationPhotonNLO     = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT7_2760GeV");
            } else {
                graphDirectPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT1_2760GeV");
                graphPromptPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT1_2760GeV");
                graphFragmentationPhotonNLO     = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT1_2760GeV");
            }
        } else if (!fEnergy.CompareTo("5TeV") || !fEnergy.CompareTo("pPb_5.023TeV")) {
            cout << "WARNING: No inv. yield calculations in compilation file yet!" << endl;
            graphDirectPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT7_7TeV");
            graphPromptPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT7_7TeV");
            graphFragmentationPhotonNLO     = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT7_7TeV");
        } else if (!fEnergy.CompareTo("7TeV")) {
            if (!triggerCutNumber.CompareTo("1") && !subTriggerCutNumber.CompareTo("0")) {   // kINT7
                graphDirectPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT7_7TeV");
                graphPromptPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT7_7TeV");
                graphFragmentationPhotonNLO     = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT7_7TeV");
            } else {
                graphDirectPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT1_7TeV");
                graphPromptPhotonNLO            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT1_7TeV");
                graphFragmentationPhotonNLO     = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT1_7TeV");
            }
        } else if (!fEnergy.CompareTo("8TeV")) {
            graphDirectPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYield_8TeV");
            graphPromptPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYield_8TeV");
            graphFragmentationPhotonNLO         = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYield_8TeV");
        } else if (!fEnergy.CompareTo("13TeV")) {
            graphDirectPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT7_7TeV");
            graphPromptPhotonNLO                = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT7_7TeV");
            graphFragmentationPhotonNLO         = (TGraphAsymmErrors*)directoryGamma->Get("graphFragmentationPhotonNLOVogelsangInvYieldINT7_7TeV");
            cout << "WARNING: No 13TeV inv. yield calculations in compilation file yet! (using 7TeV calculations)" << endl;
        } else {
            cout << "ERROR: Energy or collision system not known!" << endl;
            doNLOComparison                     = kFALSE;
        }
    }
    
    //****************************************************************************************
    //***                  Ony proceed with NLO calculations if available ********************
    //****************************************************************************************
    if (doNLOComparison){
        // scale with Ncoll (for pPb and PbPb)
        graphDirectPhotonNLO                    = (TGraphAsymmErrors*)ScaleGraph(graphDirectPhotonNLO, fNcoll);
        graphPromptPhotonNLO                    = (TGraphAsymmErrors*)ScaleGraph(graphPromptPhotonNLO, fNcoll);
        graphFragmentationPhotonNLO             = (TGraphAsymmErrors*)ScaleGraph(graphFragmentationPhotonNLO, fNcoll);

        // Create canvas and pads
        TCanvas *canvasNLOCalculations          = GetAndSetCanvas("canvasNLOCalculations",0.,0.,1000,1350);
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
        TLegend* legendNLOCalculations          = GetAndSetLegend(0.3,0.6,8);
        legendNLOCalculations->AddEntry(graphDirectPhotonNLO,       "Direct Photon NLO Calc", "lep");
        legendNLOCalculations->AddEntry(fitNLODirectPhoton,         "Fit to Direct Photon NLO Calc");
        legendNLOCalculations->AddEntry(graphPromptPhotonNLO,       "Prompt Photon NLO Calc", "lep");
        legendNLOCalculations->AddEntry(fitNLOPromptPhoton,         "Fit to Prompt Photon NLO Calc");
        legendNLOCalculations->AddEntry(graphFragmentationPhotonNLO,"Fragmentation Photon NLO Calc", "lep");
        legendNLOCalculations->AddEntry(fitNLOFragmentationPhoton,  "Fit to Fragmentation Photon NLO Calc");
        legendNLOCalculations->AddEntry(histoMCDecaySumGammaPt,     "Decay Photons from MC");
        legendNLOCalculations->AddEntry(histoGammaSpecCorrPurity,   "Inclusive Photons from data");
        legendNLOCalculations->SetTextSize(0.035);
        legendNLOCalculations->Draw();

        // Calculating ratio
        padNLORatios->cd();
        padNLORatios->SetLogy();

        histoRatioNLODirectPhoton               = (TH1D*) histoMCDecaySumGammaPt->Clone("histoRatioNLODirectPhoton");
        histoRatioNLODirectPhoton               = CalculateHistoRatioToFitNLO(histoRatioNLODirectPhoton,fitNLODirectPhoton,2.);
        textSizeSpectra                         = 0.1;
        SetStyleHistoTH1ForGraphs(histoRatioNLODirectPhoton, "#it{p}_{T} (GeV/#it{c})","#frac{Decay #gamma}{NLO Direct #gamma}",textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.95,0.6);

        histoRatioNLODirectPhoton->GetYaxis()->CenterTitle(kTRUE);
        DrawGammaSetMarker(histoRatioNLODirectPhoton, 24, 0.5, kRed+2, kRed+2);

        histoRatioNLODirectPhoton->Draw();
        canvasNLOCalculations->Print(Form("%s/%s_%s_NLOCalculations.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        cout << "NLO calculations done" << endl;
        delete histoDummy;
        delete padNLOHistos;
        delete padNLORatios;
        delete canvasNLOCalculations;
    }


    //**********************************************************************************
    //***                      Inclusive Ratio                                       ***
    //**********************************************************************************
    Bool_t doInclusiveRatios                    = kTRUE;
    if (doInclusiveRatios){
        // Create ratio with true efficiency
        TCanvas* canvasIncRatio    = GetAndSetCanvas("canvasIncRatio",0.095,0.09,1000,815);
        SetHistogramm(histoIncRatioPurityTrueEff, "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
        DrawGammaSetMarker(histoIncRatioPurityTrueEff, 20, 2.0, 4, 4);
        histoIncRatioPurityTrueEff->Draw();
        if (graphInclRatioSysErr) {
            DrawGammaSetMarkerTGraphAsym(graphInclRatioSysErr,26,0,kBlue-8,kBlue-8,2,kTRUE);
            graphInclRatioSysErr->Draw("p,2,same");
        }
        PlotLatexLegend(0.95, 0.93-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
        canvasIncRatio->Print(Form("%s/%s_%s_IncRatioPurity_trueEff.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));

        // Create ratio for MC
        canvasIncRatio->cd();
        SetHistogramm(histoMCIncRatio, "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
        DrawGammaSetMarker(histoMCIncRatio, 24, 2.0, 2, 2);
        histoMCIncRatio->Draw();
        PlotLatexLegend(0.95, 0.93-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
        canvasIncRatio->Print(Form("%s/%s_MC_IncRatio.%s",outputDir.Data(),nameOutputLabel.Data(),suffix.Data()));
        
        // Create plot with both ratios
        canvasIncRatio->cd();
        histoIncRatioPurityTrueEff->Draw("e1");
        histoMCIncRatio->Draw("e1,same");
        TLegend* legendIncRatioAll              = GetAndSetLegend2(0.18,0.93-2*0.045,0.5,0.93,0.045,1,"",42,0.15);
        legendIncRatioAll->AddEntry(histoIncRatioPurityTrueEff,"Ratio with true Eff Purity");
        legendIncRatioAll->AddEntry(histoMCIncRatio,"MC Ratio");
        legendIncRatioAll->Draw();
        PlotLatexLegend(0.95, 0.93-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
        canvasIncRatio->Print(Form("%s/%s_%s_IncRatio_all.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        delete canvasIncRatio;
        cout << "Inclusive ratios plotted" << endl;
    }
  
    //**********************************************************************************
    //***                      Photon spectra data                                   ***
    //**********************************************************************************
    TCanvas* canvasGammaSpectraSingle              = GetAndSetCanvas("canvasGammaSpectraSingle", 0.12, 0.1, 1000 ,1350);
    DrawGammaCanvasSettings( canvasGammaSpectraSingle, 0.16, 0.02, 0.015, 0.07);
    canvasGammaSpectraSingle->SetLogy();
    canvasGammaSpectraSingle->SetLogx();
    
    textSizeSpectra=0.035;
    Double_t        minPt                       = 0.2;
    Double_t        maxPt                       = 20;

    TH2F * histoDummy2                      = new TH2F("histoDummy2","histoDummy2",10000,0.12,20,10000,5e-9, 90);
    SetStyleHistoTH2ForGraphs(histoDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                0.85*textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.77,1.7);
    histoDummy2->GetXaxis()->SetLabelOffset(-0.01);
    histoDummy2->Draw();
        DrawGammaSetMarker(histoGammaSpecCorrPurity, 20, 1.5,kBlue+2,kBlue+2);
        histoGammaSpecCorrPurity->GetYaxis()->SetRangeUser(3e-9, 10);
        
        histoGammaSpecCorrPurity->DrawCopy("e1,same");
        if (graphGammaYieldSysErr) {
            DrawGammaSetMarkerTGraphAsym(graphGammaYieldSysErr,26,0,kBlue-8,kBlue-8,2,kTRUE);
            graphGammaYieldSysErr->Draw("p,2,same");
        }

        Int_t nEntriesGammaSpec     = 1;
        if(graphGammaYieldSysErr)
            nEntriesGammaSpec++;
        TLegend* legendGammaDirSpectra              = GetAndSetLegend2(0.2,0.13,0.5,0.13+0.045*nEntriesGammaSpec,0.045,1,"",42,0.15);
        legendGammaDirSpectra->AddEntry(histoGammaSpecCorrPurity,"Inclusive photons");
        if(graphGammaYieldSysErr)
            legendGammaDirSpectra->AddEntry(graphGammaYieldSysErr,"Systematic uncertainty");
        legendGammaDirSpectra->Draw();
        
        PlotLatexLegend(0.93, 0.93-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
        
    canvasGammaSpectraSingle->Print(Form("%s/%s_%s_GammaSpectrum.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
    delete canvasGammaSpectraSingle;
  
    //**********************************************************************************
    //***                      Fitting photon spectrum                               ***
    //**********************************************************************************
    Bool_t doPhotonSpectra                      = kTRUE;
    if (doInclusiveRatios){
       TString fitNameGammaA                        = "h";
       TString fitNameGammaB                        = "l";
        
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
        if(fEnergy.CompareTo("PbPb_2.76TeV") == 0){
            if(centCutNumberI<4){
                fitNameGammaA                   = "QCD";
                fitNameGammaB                   = "oHag";
                fitMaxPt                        = 14;
            } else{
                fitNameGammaA                   = "QCD";
                fitNameGammaB                   = "oHag";
                fitMaxPt                        = 14;
            }
        }
       
        TString fitOptions                      = "QNRME+";
        if (!fEnergy.CompareTo("7TeV") && mode == 4) {
            fitOptions                          = "QNRME+";
        }

        // Start fitting hagedorn and tsallis
        fitGammaA                               = FitObject(fitNameGammaA,"fitGammaA","Pi0",histoGammaSpecCorrPurity,fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(fitGammaA, 1, 2.0, kBlue-2);
        fileFinalResults << "CorrectedYieldTrueEff " << fitNameGammaA << endl;
        forOutput                               = WriteParameterToFile(fitGammaA);
        fileFinalResults << forOutput.Data() << endl;

        fitGammaB                               = FitObject(fitNameGammaB,"fitGammaB","Pi0",histoGammaSpecCorrPurity,fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(fitGammaB, 2, 2.0, kRed-3);
        fileFinalResults << "CorrectedYieldTrueEff " << fitNameGammaB << endl;
        forOutput                               = WriteParameterToFile(fitGammaB);
        fileFinalResults << forOutput.Data() << endl;

        fitGammaB->SetLineColor(2);
        fitGammaB->SetLineStyle(2);
        
        // Create canvas and pads with 3:1 ratio
        TCanvas* canvasGammaFits          = GetAndSetCanvas("ConversionGammaSpeccanvas", 0., 0., 1000 ,1350);
        DrawGammaCanvasSettings( canvasGammaFits, 0.12, 0.015, 0.01, 0.09);
        TPad* padGammaFits                = new TPad("padGammaFitsHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padGammaFits, 0.14, 0.017, 0.01, 0.);
        padGammaFits->Draw();
        TPad* padGammaFitsRatio           = new TPad("padGammaFitsRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padGammaFitsRatio, 0.14, 0.017, 0.0, 0.22);
        padGammaFitsRatio->Draw();
        padGammaFits->cd();
        padGammaFits->SetLogy();
        padGammaFits->SetLogx();
        
        // Upper part of plot with histograms
        Double_t        ptMin                   = 0.3;
        if (mode==4)    ptMin                   = 1.5;
        Double_t        ptMax                   = 40.;
        histoDummy2->Draw();
            DrawGammaSetMarker(histoGammaSpecCorrPurity, 24, 2.0, kBlack, kBlack);
            fitGammaA->SetLineColor(kBlue-2);
            fitGammaB->SetLineColor(kRed-3);

            histoGammaSpecCorrPurity->Draw("same,e1");
            fitGammaA->Draw("same");
            fitGammaB->Draw("Csame");
    
            PlotLatexLegend(0.93, 0.95-0.045*2, 0.045,collisionSystem,detectionProcess,2,31);
            TLegend* legendGammaSpectraFits  = GetAndSetLegend2(0.7, 0.93-0.045*4, 0.9, 0.93-0.045*2, 0.045,1,"",42,0.2);
            legendGammaSpectraFits->AddEntry(histoGammaSpecCorrPurity,"#gamma data", "lp");
            legendGammaSpectraFits->AddEntry(fitGammaA,"Hagedorn", "l");
            legendGammaSpectraFits->AddEntry(fitGammaB,"Levy", "l");
            legendGammaSpectraFits->Draw();
            
        // Lower part of plot with ratio
        padGammaFitsRatio->cd();
        padGammaFitsRatio->SetLogx();

            histoRatioFitGammaA                     = (TH1D*) histoGammaSpecCorrPurity->Clone("histoRatioFitGammaA");
            histoRatioFitGammaA                     = CalculateHistoRatioToFit(histoRatioFitGammaA,fitGammaA);
            histoRatioFitGammaB                     = (TH1D*) histoGammaSpecCorrPurity->Clone("histoRatioFitGammaB");
            histoRatioFitGammaB                     = CalculateHistoRatioToFit(histoRatioFitGammaB,fitGammaB);
            textSizeSpectra                         = 0.1;

            TH1D* dummy                             = new TH1D("dummy", "dummy",1000, 0.12, 20);
            SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})", "data/fit", textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.95,0.6);
            dummy->GetYaxis()->SetRangeUser(0.5, 1.55);
            dummy->GetXaxis()->SetLabelOffset(-0.02);

            DrawGammaSetMarker(histoRatioFitGammaA, 20, 1.5, kBlue-2, kBlue-2);
            DrawGammaSetMarker(histoRatioFitGammaB, 20, 1.5, kRed-3, kRed-3);

            dummy->Draw();
            histoRatioFitGammaA->Draw("e1,same");
            DrawGammaLines(0.12, 20, 1.0, 1.0,1.0, kGray+2 ,7);
            histoRatioFitGammaB->Draw("e1,same");
        
        // Save plot
        canvasGammaFits->Print(Form("%s/%s_%s_Spectra_Gamma.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        cout << "Photon spectra plotted" << endl;
        
        delete dummy;
        delete padGammaFits;
        delete padGammaFitsRatio;
        delete canvasGammaFits;
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
        
        // Special fitting fEnergys for PbPb
        if(fEnergy.CompareTo("PbPb_2.76TeV") == 0){
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
       
       TString fitOptions                       = "QNRME+";
       if (!fEnergy.CompareTo("7TeV") && mode == 4) {
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

        TCanvas *canvasGammaSpectra             = GetAndSetCanvas("canvasGammaSpectra",0., 0., 1000 ,1350);
        DrawGammaCanvasSettings( canvasGammaSpectra, 0.12, 0.015, 0.01, 0.09);
        TPad* padGammaHistos                    = new TPad("padGammaHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padGammaHistos, 0.14, 0.017, 0.01, 0.);
        padGammaHistos->Draw();
        TPad* padGammaRatios                    = new TPad("padGammaRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padGammaRatios, 0.14, 0.017, 0.0, 0.22);
        padGammaRatios->Draw();

        padGammaHistos->cd();
        padGammaHistos->SetLogy();
        padGammaHistos->SetLogx();

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
        
            PlotLatexLegend(0.93, 0.95-0.045*2, 0.045,collisionSystem,detectionProcess,2,31);

            TLegend* legendGammaSpectra       = GetAndSetLegend2(0.18, 0.05, 0.5, 0.05+0.045*4, 0.04,1,"",42,0.15);
            legendGammaSpectra->AddEntry(histoCorrectedPi0Yield,"#pi^{0} corr. yield","pl");
            legendGammaSpectra->AddEntry(fitPi0YieldA,Form("Hagedorn: %s #chi^{2}/ndf %.2f",fitPi0B.Data(),fitPi0YieldB->GetChisquare()/fitPi0YieldB->GetNDF()),"l");
            legendGammaSpectra->AddEntry(fitPi0YieldB,Form("Levy-Tsallis: %s #chi^{2}/ndf %.2f",fitPi0B.Data(),fitPi0YieldB->GetChisquare()/fitPi0YieldB->GetNDF()),"l");
            legendGammaSpectra->AddEntry(fitPi0YieldC,Form("Mod. Hagedorn: %s #chi^{2}/ndf %.2f",fitPi0C.Data(),fitPi0YieldC->GetChisquare()/fitPi0YieldC->GetNDF()),"l");
            legendGammaSpectra->Draw();

        // Plotting ratios of data and fits in lower pad
        padGammaRatios->cd();
        padGammaRatios->SetLogx();

            histoRatioFitPi0A                       = (TH1D*) histoCorrectedPi0Yield->Clone("histoRatioFitPi0A");
            histoRatioFitPi0A                       = CalculateHistoRatioToFit(histoRatioFitPi0A,fitPi0YieldA);
            histoRatioFitPi0B                       = (TH1D*) histoCorrectedPi0Yield->Clone("histoRatioFitPi0B");
            histoRatioFitPi0B                       = CalculateHistoRatioToFit(histoRatioFitPi0B,fitPi0YieldB);
            histoRatioFitPi0C                       = (TH1D*) histoCorrectedPi0Yield->Clone("histoRatioFitPi0C");
            histoRatioFitPi0C                       = CalculateHistoRatioToFit(histoRatioFitPi0C,fitPi0YieldC);

            textSizeSpectra                         = 0.1;
            TH2F * histoDummy31                     = new TH2F("histoDummy31","histoDummy31",1000,ptMin,1.5*histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()),1000,0.41, 1.59);
            SetStyleHistoTH2ForGraphs(histoDummy31, "#it{p}_{T} (GeV/#it{c})", "data/fit",textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.95,0.6);
            histoDummy31->GetYaxis()->SetLabelOffset(0.01);
            histoDummy31->GetXaxis()->SetLabelOffset(-1e-2);
            histoDummy31->DrawCopy();

            DrawGammaSetMarker(histoRatioFitPi0A, 20, 1.4,kBlue-3,kBlue-3);
            DrawGammaSetMarker(histoRatioFitPi0B, 24, 1.4,kRed-3,kRed-3);
            DrawGammaSetMarker(histoRatioFitPi0C, 25, 1.4,kGreen-2,kGreen-2);

            DrawGammaLines(ptMin,1.5*histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()), 1.0, 1.0,1.0, kGray+2 ,7);
            histoRatioFitPi0A->Draw("e1,same");
            histoRatioFitPi0B->Draw("e1,same");
            histoRatioFitPi0C->Draw("e1,same");
            
        // Save plot
        canvasGammaSpectra->Print(Form("%s/%s_%s_Spectra_Pi0.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        delete histoDummy3;
        delete histoDummy31;
        delete padGammaHistos;
        delete padGammaRatios;
        delete canvasGammaSpectra;
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
            
            PlotLatexLegend(0.95, 0.95-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
            TLegend* legendIncRatioFit          = GetAndSetLegend2(0.7,0.93-4*0.045,0.93,0.93-2*0.045,0.045,1,"",42,0.15);
            legendIncRatioFit->AddEntry(histoIncRatioPurityTrueEff,"Inclusive Ratio");
            legendIncRatioFit->AddEntry(histoIncRatioFitPurity,"Ratio #pi^{0} Fit");
            legendIncRatioFit->Draw();
        canvasIncRatioFit->Print(Form("%s/%s_%s_IncRatio_Fit.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        delete canvasIncRatioFit;
    }

    //**********************************************************************************
    //***                      Cocktail                                              ***
    //**********************************************************************************
    Bool_t doCocktailLoading                    = kTRUE;
    if (doCocktailLoading){
        cocktailFile                            = new TFile(nameFileCocktail);

        cout<<"loading cocktail file: "<<nameFileCocktail<<endl;

        if (!fEnergy.CompareTo("900GeV") || !fEnergy.CompareTo("7TeV") || !fEnergy.CompareTo("8TeV") || !fEnergy.CompareTo("13TeV") || !fEnergy.CompareTo("pPb_5.023TeV") || !fEnergy.CompareTo("PbPb_2.76TeV")){
            cocktailPi0                         = (TH1D* )cocktailFile->Get("Pi0_Pt");
            cocktailEta                         = (TH1D* )cocktailFile->Get("Eta_Pt");
            cocktailAllGamma                    = (TH1D* )cocktailFile->Get("Gamma_Pt");
            cocktailAllGammaPi0                 = (TH1D* )cocktailFile->Get("Gamma_From_Pi0_Pt");
        } else {
            cocktailPi0                         = (TH1D* )cocktailDir->Get("ptPion");
            cocktailEta                         = (TH1D* )cocktailDir->Get("ptEta");
            cocktailAllGamma                    = (TH1D* )cocktailDir->Get("sumgammapi0");
            cocktailAllGammaPi0                 = (TH1D* )cocktailDir->Get("sumgammapi0");
        }
        if (cocktailAllGamma)      cout << "found cocktailGamma"    << endl;
        if (cocktailAllGammaPi0)   cout << "found cocktailGammaPi0" << endl;
        if (cocktailPi0)           cout << "found cocktailPi0"      << endl;
        if (cocktailEta)           cout << "found cocktailEta"      << endl;
    }
    
    //**********************************************************************************
    //***                      Meson spectra data - cocktail                         ***
    //**********************************************************************************
    Bool_t doMesonSpectra                       = kTRUE;
    if (doMesonSpectra){
        TCanvas* canvasMesonSpectra             = GetAndSetCanvas("canvasMesonSpectra", 0.12, 0.1, 1000 ,1350);
        DrawGammaCanvasSettings( canvasMesonSpectra, 0.16, 0.02, 0.015, 0.07);
        canvasMesonSpectra->SetLogy();
        canvasMesonSpectra->SetLogx();

        Double_t        minPt                   = 0.2;
        if (mode==4)    minPt                   = 1.5;
        textSizeSpectra                         = 0.04;

        TH1F * histoDummy4                      = new TH1F("histoDummy4","histoDummy4",1000,minPt,1.2*histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
        SetStyleHistoTH1ForGraphs(histoDummy4, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                  0.85*textSizeSpectra, textSizeSpectra, 0.85*textSizeSpectra, textSizeSpectra, 0.77,1.7);
        histoDummy4->GetXaxis()->SetLabelOffset(-0.02);
        histoDummy4->GetYaxis()->SetRangeUser(FindSmallestBin1DHist(histoCorrectedPi0Yield)/100.,FindLargestBin1DHist(histoCorrectedPi0Yield)*50.);
        histoDummy4->DrawCopy();
        
            DrawGammaSetMarker(cocktailPi0, 20, 2.0,kRed+2,kRed+2);
            cocktailPi0->DrawCopy("e1,same");
            
            fitPi0YieldB->SetLineColor(1);
            fitPi0YieldB->DrawCopy("same");

            DrawGammaSetMarker(histoCorrectedPi0Yield, 4, 2.0, 1, 1);
            histoCorrectedPi0Yield->Draw("e1,same");
            
            Int_t nCocktails                        = 0;
            if (haveEta){
                if(cocktailEta ){
                    DrawGammaSetMarker(cocktailEta, 20, 2.0,kBlue+1,kBlue+1);
                    cocktailEta->DrawCopy("e1,same");
                    nCocktails++;
                }    
                DrawGammaSetMarker(histoCorrectedEtaYield, 4, 2.0, kGray+2, kGray+2);
                histoCorrectedEtaYield->Draw("e1,same");
                nCocktails++;
            }

            PlotLatexLegend(0.93, 0.93-0.045*2, 0.045,collisionSystem,detectionProcess,3,31);
            
            TLegend* legendMesonSpectra                 = GetAndSetLegend2(0.19, 0.11, 0.45, 0.11+0.045*(3+nCocktails), 0.045,1,"",42,0.22); 
            legendMesonSpectra->AddEntry(cocktailPi0,"Cocktail #pi^{0}");
            legendMesonSpectra->AddEntry(histoCorrectedPi0Yield,"#pi^{0} data");
            legendMesonSpectra->AddEntry(fitPi0YieldB,"fit to #pi^{0} data","l");
            if (haveEta){
                if(cocktailEta)legendMesonSpectra->AddEntry(cocktailEta,"Cocktail #eta");
                legendMesonSpectra->AddEntry(histoCorrectedEtaYield,"#eta data");
            }    
            legendMesonSpectra->Draw();

        canvasMesonSpectra->Print(Form("%s/%s_%s_Cocktail_MesonSpectra.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        delete canvasMesonSpectra;
    }

    //**********************************************************************************
    //***                      NLO Direct Photon Ratio                               ***
    //**********************************************************************************
    if (doNLOComparison){
        cocktailAllGammaNLO                     = (TH1D*) cocktailAllGamma->Clone("cocktailAllGammaNLO");
        graphDirectPhotonNLOCopy                = (TGraphAsymmErrors*)graphDirectPhotonNLO->Clone("graphNLOCalcCopy");

        xVal                                    = graphDirectPhotonNLOCopy->GetX();
        xErr                                    = graphDirectPhotonNLOCopy->GetEX();
        yVal                                    = graphDirectPhotonNLOCopy->GetY();
        yErrUp                                  = graphDirectPhotonNLOCopy->GetEYhigh();
        yErrDown                                = graphDirectPhotonNLOCopy->GetEYlow();

        TString cocktailFit                     = "xqcd";
        TString fitOptions                      = "QNRME+";//"IQNRME+";
        if (!fEnergy.CompareTo("7TeV") && mode == 4) {
            fitOptions                          = "QNRME+";
        }

        fitCocktailAllGammaForNLO               = (TF1*)FitObject(cocktailFit,"cocktailFit","Pi0",cocktailAllGammaNLO,2.0,16,NULL,fitOptions);

        for (Int_t bin=0; bin<graphDirectPhotonNLOCopy->GetN(); bin++) {
            yVal[bin]                           = (1 + ( yVal[bin] / (fitCocktailAllGammaForNLO->Eval(xVal[bin]))));
        }
        
        // ------------------------------ NLO Calculations -----------------------------
        graphNLODoubleRatio                     = new TGraphAsymmErrors(graphDirectPhotonNLOCopy->GetN(), xVal, yVal, xErr, xErr, yErrDown, yErrUp);
        graphNLODoubleRatio->SetName("graphNLODoubleRatio");
        graphNLODirGammaSpectra                 = (TGraphAsymmErrors*)graphDirectPhotonNLO->Clone("graphNLODirGammaSpectra");
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

        histoDoubleRatioTrueEffPurity  = (TH1D*) histoIncRatioPurityTrueEff->Clone("DoubleRatioConversionTrueEffPurity");
        histoDoubleRatioTrueEffPurity->Divide(cocktailAllGammaPi0);
        if (doSysErr)    graphDoubleRatioSysErr  = CalculateSysErrFromRelSysHisto( histoDoubleRatioTrueEffPurity, "DoubleRatioSystError",relSystErrorDoubleRatioDown , relSystErrorDoubleRatioUp, 2, nPointsDoubleRatio);
        else             graphDoubleRatioSysErr  = NULL;
            
        
        histoDoubleRatioFitPi0YieldPurity        = (TH1D*) histoIncRatioFitPurity->Clone("DoubleRatioConversionFitPurity");
        histoDoubleRatioFitPi0YieldPurity->Divide(cocktailAllGammaPi0);
        if (doSysErr)graphDoubleRatioFitSysErr   = CalculateSysErrFromRelSysHisto( histoDoubleRatioFitPi0YieldPurity, "DoubleRatioFitError",relSystErrorDoubleRatioDown , relSystErrorDoubleRatioUp, 2, nPointsDoubleRatio);
        else         graphDoubleRatioFitSysErr   = NULL;

        // double ratio combined
        TCanvas *canvasDoubleRatio = new TCanvas("canvasDoubleRatio","",0.095,0.09,1000,815);
        DrawGammaCanvasSettings( canvasDoubleRatio, 0.09, 0.02, 0.02, 0.09);
        canvasDoubleRatio->cd();
        canvasDoubleRatio->SetLogx();

        Double_t minPt                           = 0.3;
        Double_t maxPt                           = 40.;
        Double_t minY                            = 0.85;
        Double_t maxY                            = 1.65;
        if (mode==4) {
            minPt                                = 1.5;
            minY                                 = 0.5;
        }

        TString      method                      = "PCM";
        if (mode==4) method                      = "EMCal";

        TH2F * histo2DDoubleRatioPlotting       = new TH2F("histo2DDoubleRatioPlotting","histo2DDoubleRatioPlotting",1000,minPt,maxPt,1000,minY,maxY);
        SetStyleHistoTH2ForGraphs(histo2DDoubleRatioPlotting, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.04,0.04, 0.04,0.04, 1.,1.);
        histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.,histoDoubleRatioTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity->GetNbinsX())*1.5);
        histo2DDoubleRatioPlotting->GetXaxis()->SetTitleFont(62);
        histo2DDoubleRatioPlotting->GetYaxis()->SetTitleFont(62);
        histo2DDoubleRatioPlotting->GetXaxis()->SetLabelOffset(-1e-2);
        histo2DDoubleRatioPlotting->DrawCopy();

            DrawGammaSetMarker(histoDoubleRatioTrueEffPurity, 20, 2.0, kBlack, kBlack);
            DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurity, 20, 2.0, kBlue+2, kBlue+2);
            if (graphDoubleRatioSysErr)      DrawGammaSetMarkerTGraphAsym(graphDoubleRatioSysErr,26,0,kGray+2,kGray+2,2,kTRUE);
            if (graphDoubleRatioFitSysErr)   DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitSysErr,26,0,kBlue-8,kBlue-8,2,kTRUE);

            histoDoubleRatioFitPi0YieldPurity->DrawCopy("same,E1");
            DrawGammaLines(0., histoDoubleRatioTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
            histoDoubleRatioTrueEffPurity->DrawCopy("same,E1");
            if (graphDoubleRatioSysErr)      graphDoubleRatioSysErr->Draw("p,2,same");
            if (graphDoubleRatioFitSysErr)   graphDoubleRatioFitSysErr->Draw("p,2,same");

            TLegend* legendDoubleConversionFit       = GetAndSetLegend2(0.14,0.93-2*0.045,0.5,0.93,0.045,1,"",42,0.2); 
            legendDoubleConversionFit->AddEntry(histoDoubleRatioTrueEffPurity,"Data","p");
            legendDoubleConversionFit->AddEntry(histoDoubleRatioFitPi0YieldPurity,"Data, fitted #pi^{0}","p");
            legendDoubleConversionFit->Draw();

            PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
            
        canvasDoubleRatio->Print(Form("%s/%s_%s_DoubleRatioComparison.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        
        // double ratio fit
        canvasDoubleRatio->cd();        
        histo2DDoubleRatioPlotting->DrawCopy();
        
            DrawGammaLines(0., 20. , 1.0, 1.0,2.0, kGray+2 ,7);
            histoDoubleRatioFitPi0YieldPurity->DrawCopy("same,E1");
            if (graphDoubleRatioFitSysErr)   graphDoubleRatioFitSysErr->Draw("p,2,same");
            
            TLegend* legendDoubleConversionFit2       = GetAndSetLegend2(0.14,0.93-1*0.045,0.5,0.93,0.045,1,"",42,0.2); 
            legendDoubleConversionFit2->AddEntry(histoDoubleRatioFitPi0YieldPurity,"Data, fitted #pi^{0}","p");
            legendDoubleConversionFit2->Draw();
            PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
        
        canvasDoubleRatio->Print(Form("%s/%s_%s_DoubleRatioFit.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));

        // double ratio no fit
        canvasDoubleRatio->cd();
        histo2DDoubleRatioPlotting->DrawCopy();
        
            DrawGammaLines(0., histoDoubleRatioTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
            histoDoubleRatioTrueEffPurity->DrawCopy("same,E1");
            if (graphDoubleRatioSysErr) graphDoubleRatioSysErr->Draw("p,2,same");
            PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
            
            TLegend* legendDoubleConversionFit3       = GetAndSetLegend2(0.14,0.93-1*0.045,0.5,0.93,0.045,1,"",42,0.2); 
            legendDoubleConversionFit3->AddEntry(histoDoubleRatioTrueEffPurity,"Data","p");
            legendDoubleConversionFit3->Draw();
        
        canvasDoubleRatio->Print(Form("%s/%s_%s_DoubleRatio.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        delete canvasDoubleRatio;
        
        // do theory comparisons
        if (doNLOComparison){
            graphNLODoubleRatio->SetLineColor(kAzure);
            graphNLODoubleRatio->SetFillColor(kAzure);
            graphNLODoubleRatio->SetLineWidth(3.0);
            graphNLODoubleRatio->SetMarkerSize(0);
            
            graphNLODoubleRatio->Print();
            graphNLODirGammaSpectra->SetLineColor(kAzure);
            graphNLODirGammaSpectra->SetFillColor(kAzure);
            graphNLODirGammaSpectra->SetLineWidth(3.0);
            graphNLODirGammaSpectra->SetMarkerSize(0);
            graphNLODirGammaSpectra->RemovePoint(0);
            
            fitNLODoublRatio                   = new TF1("fitNLODoublRatio","(x<=2)*(1.0+[0]*x)+(x>2)*([1]+[2]*x+[3]*x*x)",0.,4.0);
            graphNLODoubleRatio->Fit(fitNLODoublRatio,"NRME+","",0,4.0);
            fitNLODoublRatio->SetLineColor(2);
            
            Int_t nLinesNLOLegends  = 2;
            if (fEnergy.CompareTo("PbPb_2.76TeV") == 0 || fEnergy.CompareTo("pPb_5.023TeV") == 0) 
                nLinesNLOLegends    = 3;
            // ------------------------------ NLO Calculations -----------------------------
            // combined
            TCanvas *canvasDoublRatioNLO = new TCanvas("canvasDoublRatioNLO","",0.095,0.09,1000,815);
            DrawGammaCanvasSettings( canvasDoublRatioNLO, 0.09, 0.02, 0.02, 0.09);
            canvasDoublRatioNLO->cd();
            canvasDoublRatioNLO->SetLogx();
            
            histo2DDoubleRatioPlotting->DrawCopy();
                graphNLODoubleRatio->RemovePoint(0);
                DrawGammaLines(0., histoDoubleRatioTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
                graphNLODoubleRatio->Draw("lp3");

                DrawGammaSetMarker(histoDoubleRatioTrueEffPurity, 20, 2.0, kBlack, kBlack);
                DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurity, 20, 2.0, kBlue+2, kBlue+2);
                histoDoubleRatioTrueEffPurity->DrawCopy("same");
                histoDoubleRatioFitPi0YieldPurity->DrawCopy("same");
                if (graphDoubleRatioSysErr)     graphDoubleRatioSysErr->Draw("p,2,same");
                if (graphDoubleRatioFitSysErr)  graphDoubleRatioFitSysErr->Draw("p,2,same");
                
                TLegend* legendDoubleConversionFitNLO = GetAndSetLegend2(0.14,0.93-(nLinesNLOLegends+1+0.5)*0.045,0.5,0.93,0.045,1,"",42,0.2); 
                legendDoubleConversionFitNLO->AddEntry(histoDoubleRatioTrueEffPurity,"Data","p");
                legendDoubleConversionFitNLO->AddEntry(histoDoubleRatioFitPi0YieldPurity,"Data, fitted #pi^{0}","p");
                if(fEnergy.CompareTo("PbPb_2.76TeV") == 0)       legendDoubleConversionFitNLO->AddEntry(graphNLODoubleRatio,"pp #gamma_{dir} NLO #sqrt{s} = 2.76 TeV","l");
                else if(fEnergy.CompareTo("pPb_5.023TeV") == 0)  legendDoubleConversionFitNLO->AddEntry(graphNLODoubleRatio,"pp #gamma_{dir} NLO #sqrt{s} = 5.02 TeV ","l");
                else                                             legendDoubleConversionFitNLO->AddEntry(graphNLODoubleRatio,"pp NLO Direct Photon","l");
                if (fEnergy.CompareTo("PbPb_2.76TeV") == 0 || fEnergy.CompareTo("pPb_5.023TeV") == 0) legendDoubleConversionFitNLO->AddEntry((TObject*)0,"scaled by N_{coll}","");

                legendDoubleConversionFitNLO->Draw();

                PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
            
            canvasDoublRatioNLO->Print(Form("%s/%s_%s_DoubleRatioComparison_NLO.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));

            // fit
            canvasDoublRatioNLO->cd();           
            histo2DDoubleRatioPlotting->DrawCopy();
                DrawGammaLines(0., histoDoubleRatioTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
                graphNLODoubleRatio->Draw("lp3");
                
                histoDoubleRatioFitPi0YieldPurity->DrawCopy("same");
                if (graphDoubleRatioFitSysErr) graphDoubleRatioFitSysErr->Draw("p,2,same");
                
                TLegend* legendDoubleConversionFitNLO2 = GetAndSetLegend2(0.14,0.93-(nLinesNLOLegends+0.5)*0.045,0.5,0.93,0.045,1,"",42,0.2); 
                legendDoubleConversionFitNLO2->AddEntry(histoDoubleRatioFitPi0YieldPurity,"Data, fitted #pi^{0}","p");
                if(fEnergy.CompareTo("PbPb_2.76TeV") == 0)       legendDoubleConversionFitNLO2->AddEntry(graphNLODoubleRatio,"pp #gamma_{dir} NLO #sqrt{s} = 2.76 TeV","l");
                else if(fEnergy.CompareTo("pPb_5.023TeV") == 0)  legendDoubleConversionFitNLO2->AddEntry(graphNLODoubleRatio,"pp #gamma_{dir} NLO #sqrt{s} = 5.02 TeV ","l");
                else                                             legendDoubleConversionFitNLO2->AddEntry(graphNLODoubleRatio,"pp NLO Direct Photon","l");
                if (fEnergy.CompareTo("PbPb_2.76TeV") == 0 || fEnergy.CompareTo("pPb_5.023TeV") == 0) legendDoubleConversionFitNLO2->AddEntry((TObject*)0,"scaled by N_{coll}","");
                legendDoubleConversionFitNLO2->Draw();
                PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
                
            canvasDoublRatioNLO->Print(Form("%s/%s_%s_DoubleRatioFit_NLO.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
            
            // measured
            canvasDoublRatioNLO->cd();                       
            histo2DDoubleRatioPlotting->DrawCopy();
                DrawGammaLines(0., histoDoubleRatioTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
                graphNLODoubleRatio->Draw("lp3");
                
                histoDoubleRatioTrueEffPurity->DrawCopy("same");
                if (graphDoubleRatioSysErr) graphDoubleRatioSysErr->Draw("p,2,same");
                
                TLegend* legendDoubleConversionFitNLO3 = GetAndSetLegend2(0.14,0.93-(nLinesNLOLegends+0.5)*0.045,0.5,0.93,0.045,1,"",42,0.2); 
                legendDoubleConversionFitNLO3->AddEntry(histoDoubleRatioTrueEffPurity,"Data","p");
                if(fEnergy.CompareTo("PbPb_2.76TeV") == 0)       legendDoubleConversionFitNLO3->AddEntry(graphNLODoubleRatio,"pp #gamma_{dir} NLO #sqrt{s} = 2.76 TeV","l");
                else if(fEnergy.CompareTo("pPb_5.023TeV") == 0)  legendDoubleConversionFitNLO3->AddEntry(graphNLODoubleRatio,"pp #gamma_{dir} NLO #sqrt{s} = 5.02 TeV ","l");
                else                                             legendDoubleConversionFitNLO3->AddEntry(graphNLODoubleRatio,"pp NLO Direct Photon","l");
                if (fEnergy.CompareTo("PbPb_2.76TeV") == 0 || fEnergy.CompareTo("pPb_5.023TeV") == 0) legendDoubleConversionFitNLO3->AddEntry((TObject*)0,"scaled by N_{coll}","");
                legendDoubleConversionFitNLO3->Draw();
                PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
                
            canvasDoublRatioNLO->Print(Form("%s/%s_%s_DoubleRatio_NLO.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
            delete canvasDoublRatioNLO;
        }
       
        // --------------------- Direct Photon Spectrum -------------------------
        if (graphDoubleRatioSysErr) {
            // calculate direct photon spectrum
            histoDirectPhotonSpectrum                = (TH1D*)histoGammaSpecCorrPurity->Clone("histoDirectPhotonSpectrum");
            histoThermalPhotonSpectrum               = (TH1D*)histoGammaSpecCorrPurity->Clone("histoThermalPhotonSpectrum");
            histoPromptPhotonSpectrum                = (TH1D*)histoGammaSpecCorrPurity->Clone("histoPromptPhotonSpectrum");
            
            Bool_t doUpperLimits                     = kFALSE;
            if (!doUpperLimits) {
                for(Int_t i = 1; i < histoDirectPhotonSpectrum->GetNbinsX(); i++){
                    histoDirectPhotonSpectrum->SetBinContent(    i+1,0);
                    histoThermalPhotonSpectrum->SetBinContent(   i+1,0);
                }
                for(Int_t i = 1; i < histoDirectPhotonSpectrum->GetNbinsX(); i++){
                    Double_t binContent              = -1;
                    Double_t binError                = -1;
                    binContent                       = 1-(1./( histoDoubleRatioTrueEffPurity->GetBinContent(i+1)));
                    binContent                       = binContent*histoGammaSpecCorrPurity->GetBinContent(i+1);
                    binError                         = 1-(1./( histoDoubleRatioTrueEffPurity->GetBinContent(i+1) + histoDoubleRatioTrueEffPurity->GetBinError(i+1)));
                    binError                         = binError*histoGammaSpecCorrPurity->GetBinContent(i+1);
                    binError                         = binError-binContent;
                    histoDirectPhotonSpectrum->SetBinContent(i+1,binContent);
                    histoDirectPhotonSpectrum->SetBinError(  i+1,binError);
                    
                    if (doNLOComparison){
                        if(histoDirectPhotonSpectrum->GetBinCenter(i+1) < 4.5){
                            binContent                   = 1-(1./(1+ histoDoubleRatioTrueEffPurity->GetBinContent(i+1) - fitNLODoublRatio->Eval(histoDoubleRatioTrueEffPurity->GetBinCenter(i+1))));
                            binContent                   = binContent*histoGammaSpecCorrPurity->GetBinContent(i+1);
                            histoThermalPhotonSpectrum->SetBinContent(   i+1,binContent);
                            histoThermalPhotonSpectrum->SetBinError(     i+1,binError);
                        } else {
                            histoThermalPhotonSpectrum->SetBinContent(   i+1,0);
                            histoThermalPhotonSpectrum->SetBinError(     i+1,0);
                        }
                        if(histoDirectPhotonSpectrum->GetBinCenter(i+1) < 4.5){
                            binContent                   = 1-(1./(fitNLODoublRatio->Eval(histoDoubleRatioTrueEffPurity->GetBinCenter(i+1))));
                            binContent                   = binContent*histoGammaSpecCorrPurity->GetBinContent(i+1);
                            histoPromptPhotonSpectrum->SetBinContent(    i+1,binContent);
                            histoPromptPhotonSpectrum->SetBinError(      i+1,0.0000001);
                        } else {
                            histoPromptPhotonSpectrum->SetBinContent(    i+1,0);
                            histoPromptPhotonSpectrum->SetBinError(      i+1,0);
                        }
                    }    
                }
            } else {
                histoDoubleRatioUpperLimits          = GetUpperLimitsHisto(histoDoubleRatioTrueEffPurity,  graphDoubleRatioSysErr, 0.95, 1e-6, 1e4);
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
            SetStyleHistoTH2ForGraphs(dummyDir, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})", 0.04,0.04, 0.04,0.04, 1.,1.1);
            dummyDir->GetXaxis()->SetRangeUser(0.,histoDoubleRatioTrueEffPurity->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity->GetNbinsX())*1.5);
            dummyDir->GetXaxis()->SetTitleFont(62);
            dummyDir->GetYaxis()->SetTitleFont(62);
            dummyDir->GetXaxis()->SetLabelOffset(-2e-2);
            dummyDir->GetXaxis()->SetTitleOffset(0.8);
            dummyDir->GetYaxis()->SetTitleOffset(1.5);
            dummyDir->DrawCopy();
            
            PlotLatexLegend(0.93, 0.95-0.045*2, 0.045,collisionSystem,detectionProcess,2,31);
            TArrow *ar2[50];
            for(Int_t i=1;i<histoDirectPhotonSpectrum->GetNbinsX()+1;i++)
            {
                if (!histoDirectPhotonSpectrum->GetBinContent(i) > 5e-10) continue;
                ar2[i] = new TArrow(histoDirectPhotonSpectrum->GetBinCenter(i),histoDirectPhotonSpectrum->GetBinContent(i),histoDirectPhotonSpectrum->GetBinCenter(i),histoDirectPhotonSpectrum->GetBinContent(i)*0.1,0.02,"|->");
                ar2[i]->SetAngle(40);
                ar2[i]->SetLineWidth(2);
                ar2[i]->Draw();
            }
            
            if (doNLOComparison){
                graphNLODirGammaSpectra->GetXaxis()->SetRangeUser(graphNLODirGammaSpectra->GetXaxis()->GetXmin(), histoDirectPhotonSpectrum->GetXaxis()->GetXmax()); 
                graphNLODirGammaSpectra->SetLineWidth(2);
                graphNLODirGammaSpectra->Draw("lp3");
            }
            
            TLegend* legendGammaDirSpectra               = GetAndSetLegend(0.2,0.2,2);
            legendGammaDirSpectra->AddEntry(histoDirectPhotonSpectrum,"direct photon upper limits",   "l");
            if (doNLOComparison) legendGammaDirSpectra->AddEntry(graphNLODirGammaSpectra,  "pp NLO direct photon",         "l");
            legendGammaDirSpectra->Draw();
            
            CanvasDirGamma->Print(Form("%s/%s_%s_DirectPhotonSpectrum.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        } else {
            cout << "No systematic uncertainties on double ratio given, skipping direct photon calculation." << endl;
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

    PlotLatexLegend(0.95, 0.93-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
    TLegend* legendIncRatioDecRatio          = GetAndSetLegend2(0.7,0.93-6*0.045,0.93,0.93-3*0.045,0.045,1,"",42,0.15);
    legendIncRatioDecRatio->AddEntry(histoIncRatioPurityTrueEff,    "(#gamma_{incl.} / #pi^{0})_{meas.}", "lp");
    legendIncRatioDecRatio->AddEntry(cocktailAllGammaPi0,           "(#gamma_{dec.} / #pi^{0})_{sim.}",  "lp");
    legendIncRatioDecRatio->Draw();

    canvasIncRatioDecRatio->Print(Form("%s/%s_%s_IncRatio_DecRatio_trueEff.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
    delete canvasIncRatioDecRatio;

    //**********************************************************************************
    //***                      Save Histograms                                       ***
    //**********************************************************************************   
    Bool_t doSaving                             = kTRUE;
    if (doSaving){
        TString nameOutputFile                  = Form("%s/%s/%s_%s_GammaConvV1_InclusiveRatio.root",cutSel.Data(),fEnergy.Data(),nameOutputLabel.Data(),nameRec.Data());
        fileCorrectedOutput                     = new TFile(nameOutputFile,"RECREATE");

        // pi0 quantities
        fitPi0YieldA->Write(            fitPi0YieldA->GetName(),            TObject::kOverwrite);
        fitPi0YieldB->Write(            fitPi0YieldB->GetName(),            TObject::kOverwrite);
        fitPi0YieldC->Write(            fitPi0YieldC->GetName(),            TObject::kOverwrite);
        histoCorrectedPi0Yield->Write(  histoCorrectedPi0Yield->GetName(),  TObject::kOverwrite);

        // gamma quantities
        if (fitGammaA)        fitGammaA->Write(         fitGammaA->GetName(),         TObject::kOverwrite);
        if (histoDirectPhotonSpectrum)  histoDirectPhotonSpectrum->Write(   histoDirectPhotonSpectrum->GetName(),   TObject::kOverwrite);
        if (histoGammaSpecCorrPurity)   histoGammaSpecCorrPurity->Write(    "histoGammaSpecCorrPurity",             TObject::kOverwrite);

        // Double ratio
        if(histoDoubleRatioFitPi0YieldPurity)           histoDoubleRatioFitPi0YieldPurity->Write(           histoDoubleRatioFitPi0YieldPurity->GetName(),               TObject::kOverwrite);
        if(histoDoubleRatioTrueEffPurity)               histoDoubleRatioTrueEffPurity->Write(     histoDoubleRatioTrueEffPurity->GetName(),         TObject::kOverwrite);
        if(histoDoubleRatioUpperLimits)                 histoDoubleRatioUpperLimits->Write(Form("%s_UpperLimits", histoDoubleRatioTrueEffPurity->GetName()),  TObject::kOverwrite);

        // systematics
        if (graphGammaYieldSysErr)      graphGammaYieldSysErr->Write(       Form("%s_SystErr", histoGammaSpecCorrPurity->GetName()),                TObject::kOverwrite);
        if (graphInclRatioSysErr)       graphInclRatioSysErr->Write(        Form("%s_SystErr", histoIncRatioPurityTrueEff->GetName()),              TObject::kOverwrite);
        if (graphDoubleRatioSysErr)     graphDoubleRatioSysErr->Write(      Form("%s_SystErr", histoDoubleRatioTrueEffPurity->GetName()), TObject::kOverwrite);
        if (graphDoubleRatioFitSysErr)  graphDoubleRatioFitSysErr->Write(   Form("%s_SystErr", histoDoubleRatioFitPi0YieldPurity->GetName()),       TObject::kOverwrite);
        
        // inclusive ratio
        histoIncRatioPurityTrueEff->Write(  histoIncRatioPurityTrueEff->GetName(),  TObject::kOverwrite);
        histoMCIncRatio->Write(             histoMCIncRatio->GetName(),             TObject::kOverwrite);
        histoIncRatioFitPurity->Write(      histoIncRatioFitPurity->GetName(),      TObject::kOverwrite);

        // cocktail
        cocktailAllGamma->Write(cocktailAllGamma->GetName(),    TObject::kOverwrite);
        cocktailPi0->Write(     cocktailPi0->GetName(),         TObject::kOverwrite);
        
        //NLO
        if (doNLOComparison){
            fitNLODirectPhoton->Write(          "fitNLODirectPhoton",           TObject::kOverwrite);
            graphPromptPhotonNLO->Write(        "graphPromptPhotonNLO",         TObject::kOverwrite);
            fitNLOPromptPhoton->Write(          "fitNLOPromptPhoton",           TObject::kOverwrite);
            graphFragmentationPhotonNLO->Write( "graphFragmentationPhotonNLO",  TObject::kOverwrite);
            fitNLOFragmentationPhoton->Write(   "fitNLOFragmentationPhoton",    TObject::kOverwrite);

            histoRatioNLODirectPhoton->Write(    "histoRatioNLODirectPhoton",     TObject::kOverwrite);
            graphNLODoubleRatio->Write(         "graphgraphNLODoubleRatio",     TObject::kOverwrite);
            graphNLODirGammaSpectra->Write(     "graphNLODirGamma",             TObject::kOverwrite);
        }
        
        fileCorrectedOutput->Close();
        delete fileCorrectedOutput;
    }
}