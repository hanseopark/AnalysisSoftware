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
#include "CalculateGammaToPi0V4.h"
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
void  CalculateGammaToPi0V4(    TString nameFileGamma       = "",
                                TString nameFilePi0         = "",
                                TString nameFileCocktail    = "",
                                TString cutSel              = "",
                                TString suffix              = "pdf",
                                TString nameMeson           = "",
                                TString isMC                = "",
                                TString fEnergy             = "",
                                Int_t mode                  = 0,
                                TString nameFileFitsShift   = ""

                            ){

    // Setting the general style
    gROOT->Reset();
    StyleSettingsThesis();
    SetPlotStyle();

    // Separating cutnumber and retrieving centrality and number of events
    SeparateCutnumberString(cutSel, mode, fEnergy);

    // Create strings for naming
    CreateNamingStrings(nameMeson,isMC);
    TString collisionSystem                      = ReturnFullCollisionsSystem(fEnergy);
    TString detectionProcess                     = ReturnFullTextReconstructionProcess(mode);
    if(fEnergy.Contains("PbPb")){
        collisionSystem                          = Form("%s %s", centrality.Data(), collisionSystem.Data());
    } else if(fEnergy.Contains("pPb") && !centrality.Contains("0-100")){
        collisionSystem                          = Form("%s %s", centrality.Data(), collisionSystem.Data());
    }
    TString centralityString2                    = GetCentralityString(fEventCutSelection);

    // check if fit file for binshifting has to be adjusted for every energy
    TF1* fitBinShiftPi0                          = 0x0;
    TF1* fitBinShiftGamma                        = 0x0;
    Bool_t doBinShiftForDR                       = kFALSE;
    TString addNameBinshift                      = "";
    if (nameFileFitsShift.CompareTo("") != 0 ){
        doBinShiftForDR                          = kTRUE;
        addNameBinshift                          = "YShifted";
    }

    if (doBinShiftForDR){
        TFile *fileFitsBinShift                = new TFile(nameFileFitsShift);
        fitBinShiftPi0                            = (TF1*)fileFitsBinShift->Get("TsallisFitPi0");
        if(!fitBinShiftPi0)
          fitBinShiftPi0                          = (TF1*)fileFitsBinShift->Get("TwoComponentModelFitPi0");
        if(!fitBinShiftPi0)
          fitBinShiftPi0                          = (TF1*)fileFitsBinShift->Get(Form("Pi0%s/TsallisFitPi0",fEnergy.Data()));
        if(!fitBinShiftPi0)
          fitBinShiftPi0                          = (TF1*)fileFitsBinShift->Get(Form("Pi0%s/Fits/TsallisFitPi0",fEnergy.Data()));
        if(!fitBinShiftPi0)
          fitBinShiftPi0                          = (TF1*)fileFitsBinShift->Get(Form("Pi0%s/Fits/TwoComponentModelFitPi0",fEnergy.Data()));
        if(!fitBinShiftPi0){
          cout << "missing pi0 bin shift fit... will not do binshifting!" << endl;
          doBinShiftForDR                         = kFALSE;
        }

        fitBinShiftGamma                          = (TF1*)fileFitsBinShift->Get(Form("TsallisFitGamma"));
        if(!fitBinShiftGamma)
          fitBinShiftGamma                        = (TF1*)fileFitsBinShift->Get(Form("TwoComponentModelFitGamma"));
        if(!fitBinShiftGamma)
          fitBinShiftGamma                        = (TF1*)fileFitsBinShift->Get(Form("ModHagedornFitGamma"));
        if(!fitBinShiftGamma){
          cout << "missing gamma bin shift fit... will not do binshifting!" << endl;
          doBinShiftForDR                         = kFALSE;
        }
        if(fitBinShiftGamma && fitBinShiftPi0){
          cout << fitBinShiftPi0 << " - " << fitBinShiftGamma << endl;
          cout << "fits for shifting found " << endl;
        }
    }
    // Creating output directory and output file
    cout << "Output directory with plots:" << endl;
    cout << Form("%s/%s/%s/GammaToPi0V4",cutSel.Data(),fEnergy.Data(),suffix.Data()) << endl;
    TString outputDir                            = Form("%s/%s/%s/GammaToPi0V4",cutSel.Data(),fEnergy.Data(),suffix.Data());
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
    if(doBinShiftForDR){
      histoGammaSpecCorrPurityBinShift          = (TH1D*)histoGammaSpecCorrPurity->Clone("CorrectedYieldGammaBinShift");
      histoGammaSpecCorrPurityBinShift          = ApplyYshiftIndividualSpectra( histoGammaSpecCorrPurityBinShift, fitBinShiftGamma);
      cout << "bin shifted gamma" << endl;
    }
    histoMCDecaySumGammaPt                      = (TH1D*)fileGamma->Get("MC_DecayGammaAll_Pt");
    histoGammaSpecMCAll                         = (TH1D*)fileGamma->Get("GammaSpecMC");

    // Opening pion file and loading spectra
    filePi0                                     = new TFile(nameFilePi0);
    for (Int_t k = 0; k < 6; k++){
        if(mode == 4 || mode == 5) histoCorrectedPi0Yield[k]               = (TH1D*)filePi0->Get(Form("CorrectedYieldNormEff%s",nameIntRanges[k].Data()));
        else                       histoCorrectedPi0Yield[k]               = (TH1D*)filePi0->Get(Form("CorrectedYieldTrueEff%s",nameIntRanges[k].Data()));
    }
    if(doBinShiftForDR){
      histoCorrectedPi0YieldBinShift            = (TH1D*)histoCorrectedPi0Yield[0]->Clone("CorrectedYieldPi0BinShift");
      histoCorrectedPi0YieldBinShift            = ApplyYshiftIndividualSpectra( histoCorrectedPi0YieldBinShift, fitBinShiftPi0);
      cout << "bin shifted pi0" << endl;
    }
    graphSysPi0PileUpOptions                    = (TGraphAsymmErrors*)filePi0->Get(Form("Pi0_SystErrorRel_BGEstimate_%s",centralityString2.Data()));
    graphSysPi0PileUpIterations                 = (TGraphAsymmErrors*)filePi0->Get(Form("Pi0_SystErrorRel_BGEstimateIterations_%s",centralityString2.Data()));


    histoMCYieldMeson                           = (TH1D*)filePi0->Get("MCYield_Meson");
    histoMCYieldMesonOldBin                     = (TH1D*)filePi0->Get("MCYield_Meson_oldBin");

    Bool_t haveEta                              = kFALSE;
    TString nameFileEta                         = nameFilePi0.ReplaceAll("Pi0","Eta");
    fileEta                                     = new TFile(nameFileEta);
    if (!fileEta->IsZombie()){
        if(mode == 4 || mode == 5) histoCorrectedEtaYield                  = (TH1D*)fileEta->Get("CorrectedYieldNormEff");
        else                       histoCorrectedEtaYield                  = (TH1D*)fileEta->Get("CorrectedYieldTrueEff");
        if (histoCorrectedEtaYield)
            haveEta                             = kTRUE;
    }

    // Creating inclusive gamma ratio
    for (Int_t k = 0; k < 6; k++){
        histoIncRatioPurityTrueEff[k]           = (TH1D*) histoGammaSpecCorrPurity->Clone(Form("IncRatioPurity_trueEff%s",nameIntRanges[k].Data()));
        histoIncRatioPurityTrueEff[k]->Divide(histoIncRatioPurityTrueEff[k],histoCorrectedPi0Yield[k],1,1,"");
    }
    if(doBinShiftForDR){
      histoIncRatioPurityTrueEffBinShift           = (TH1D*) histoGammaSpecCorrPurityBinShift->Clone(Form("IncRatioPurity_trueEff_%s",addNameBinshift.Data()));
      histoIncRatioPurityTrueEffBinShift->Divide(histoIncRatioPurityTrueEffBinShift,histoCorrectedPi0YieldBinShift,1,1,"");
    }
    TString nameStatErrorCheckDat               = Form("%s/%s/Gamma_%s_StatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
    fstream fileStatErrorCheck;
    fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
    fileStatErrorCheck << "# pT \t yield \t absErr \t relErr(%)" << endl;
    for (Int_t i=1; i<=histoGammaSpecCorrPurity->GetNbinsX(); i++){
        if (histoGammaSpecCorrPurity->GetBinContent(i) != 0)
            fileStatErrorCheck << histoGammaSpecCorrPurity->GetBinCenter(i) << "\t" << histoGammaSpecCorrPurity->GetBinContent(i) << "\t" << histoGammaSpecCorrPurity->GetBinError(i) << "\t" << histoGammaSpecCorrPurity->GetBinError(i)/histoGammaSpecCorrPurity->GetBinContent(i)*100 << endl;
    }
    fileStatErrorCheck.close();
    nameStatErrorCheckDat                       = Form("%s/%s/Gamma_%s_RelStatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
    fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
    fileStatErrorCheck << "pT \t GammaStat" << endl;
    for (Int_t i=1; i<=histoIncRatioPurityTrueEff[0]->GetNbinsX(); i++){
        if (histoIncRatioPurityTrueEff[0]->GetBinContent(i) != 0){
            Double_t relGammaErr    = histoGammaSpecCorrPurity->GetBinError(i)/histoGammaSpecCorrPurity->GetBinContent(i)*100 ;
            Double_t relTotErr      = histoGammaSpecCorrPurity->GetBinError(i)/histoGammaSpecCorrPurity->GetBinContent(i)*100 ;
            fileStatErrorCheck << histoIncRatioPurityTrueEff[0]->GetBinCenter(i) << "\t" << relGammaErr << "\t" << relTotErr << endl;
        }
    }
    nameStatErrorCheckDat                       = Form("%s/%s/Pi0_StatError_%s_ForDirGammaEst.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
    fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
    fileStatErrorCheck << "# pT \t yield \t absErr \t relErr(%)" << endl;
    for (Int_t i=1; i<=histoCorrectedPi0Yield[0]->GetNbinsX(); i++){
        if (histoCorrectedPi0Yield[0]->GetBinContent(i) != 0)
            fileStatErrorCheck << histoCorrectedPi0Yield[0]->GetBinCenter(i) << "\t" << histoCorrectedPi0Yield[0]->GetBinContent(i)<< "\t" << histoCorrectedPi0Yield[0]->GetBinError(i) << "\t" << histoCorrectedPi0Yield[0]->GetBinError(i)/histoCorrectedPi0Yield[0]->GetBinContent(i)*100 << endl;
    }
    fileStatErrorCheck.close();
    nameStatErrorCheckDat                       = Form("%s/%s/IncGammaToPi0_%s_StatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
    fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
    fileStatErrorCheck << "# pT \t ratio \t absErr \t relErr(%)" << endl;
    for (Int_t i=1; i<=histoIncRatioPurityTrueEff[0]->GetNbinsX(); i++){
        if (histoIncRatioPurityTrueEff[0]->GetBinContent(i) != 0)
        fileStatErrorCheck << histoIncRatioPurityTrueEff[0]->GetBinCenter(i) << "\t" << histoIncRatioPurityTrueEff[0]->GetBinContent(i)<< "\t" << histoIncRatioPurityTrueEff[0]->GetBinError(i) << "\t" << histoIncRatioPurityTrueEff[0]->GetBinError(i)/histoIncRatioPurityTrueEff[0]->GetBinContent(i)*100 << endl;
    }
    fileStatErrorCheck.close();
    nameStatErrorCheckDat                       = Form("%s/%s/IncGammaToPi0_%s_RelStatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
    fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
    fileStatErrorCheck << "pT \t GammaStat \t Pi0Stat" << endl;
    for (Int_t i=1; i<=histoIncRatioPurityTrueEff[0]->GetNbinsX(); i++){
        if (histoIncRatioPurityTrueEff[0]->GetBinContent(i) != 0){
            Double_t relGammaErr    = histoGammaSpecCorrPurity->GetBinError(i)/histoGammaSpecCorrPurity->GetBinContent(i)*100 ;
            Double_t relPi0Err      = histoCorrectedPi0Yield[0]->GetBinError(i)/histoCorrectedPi0Yield[0]->GetBinContent(i)*100;
            Double_t relTotErr      = histoIncRatioPurityTrueEff[0]->GetBinError(i)/histoIncRatioPurityTrueEff[0]->GetBinContent(i)*100;
            fileStatErrorCheck << histoIncRatioPurityTrueEff[0]->GetBinCenter(i) << "\t" << relGammaErr << "\t" << relPi0Err << "\t" << relTotErr << endl;
        }
    }
    fileStatErrorCheck.close();

    histoMCIncRatio                             = (TH1D*) histoGammaSpecMCAll->Clone("MC_IncRatio");
    histoMCIncRatio->Divide(histoGammaSpecMCAll,histoMCYieldMeson,1,1,"");

    //**********************************************************************************
    //***                      Bin Shifting                                          ***
    //**********************************************************************************
    if(doBinShiftForDR){
      TCanvas* canvasBinShifting    = GetAndSetCanvas("canvasBinShifting",0.095,0.09,1000,815);
      canvasBinShifting->cd();
      histoGammaSpecCorrPurityBinShiftCorr = (TH1D*)histoGammaSpecCorrPurityBinShift->Clone("InclusiveRatioBinShiftCorrection");
      histoGammaSpecCorrPurityBinShiftCorr->Divide(histoGammaSpecCorrPurityBinShiftCorr,histoGammaSpecCorrPurity,1,1,"");
      SetHistogramm(histoGammaSpecCorrPurityBinShiftCorr, "#it{p}_{T} (GeV/#it{c})", "Shifted / Standard",0.8,1.2);
      DrawGammaSetMarker(histoGammaSpecCorrPurityBinShiftCorr, 20, 2.0, kGreen+2, kGreen+2);
      histoGammaSpecCorrPurityBinShiftCorr->Draw();
      PlotLatexLegend(0.95, 0.96-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
      canvasBinShifting->Print(Form("%s/%s_Gamma_BinShiftCorrection.%s",outputDir.Data(),nameRec.Data(),suffix.Data()));

      canvasBinShifting->cd();
      histoCorrectedPi0YieldBinShiftCorr = (TH1D*)histoCorrectedPi0YieldBinShift->Clone("InclusiveRatioBinShiftCorrection");
      histoCorrectedPi0YieldBinShiftCorr->Divide(histoCorrectedPi0YieldBinShiftCorr,histoCorrectedPi0Yield[0],1,1,"");
      SetHistogramm(histoCorrectedPi0YieldBinShiftCorr, "#it{p}_{T} (GeV/#it{c})", "Shifted / Standard",0.8,1.2);
      DrawGammaSetMarker(histoCorrectedPi0YieldBinShiftCorr, 20, 2.0, kGreen+2, kGreen+2);
      histoCorrectedPi0YieldBinShiftCorr->Draw();
      PlotLatexLegend(0.95, 0.96-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
      canvasBinShifting->Print(Form("%s/%s_Pi0_BinShiftCorrection.%s",outputDir.Data(),nameRec.Data(),suffix.Data()));
    }

    //**********************************************************************************
    //***                      Cocktail                                              ***
    //**********************************************************************************
    Bool_t doCocktailLoading                    = kTRUE;
    if (doCocktailLoading){
        cocktailFile                            = new TFile(nameFileCocktail);

        cout<<"loading cocktail file: "<<nameFileCocktail<<endl;

        if (!fEnergy.CompareTo("900GeV") || !fEnergy.CompareTo("2.76TeV")|| !fEnergy.CompareTo("7TeV") || !fEnergy.CompareTo("8TeV") || !fEnergy.CompareTo("13TeV") ||
            !fEnergy.CompareTo("pPb_5.023TeV") ||
            !fEnergy.CompareTo("PbPb_2.76TeV")){
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
    //***                      Inclusive Ratio                                       ***
    //**********************************************************************************
    Bool_t doInclusiveRatios                    = kTRUE;
    if (doInclusiveRatios){
        // Create ratio with true efficiency
        TCanvas* canvasIncRatio    = GetAndSetCanvas("canvasIncRatio",0.095,0.09,1000,815);
        SetHistogramm(histoIncRatioPurityTrueEff[0], "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
        DrawGammaSetMarker(histoIncRatioPurityTrueEff[0], 20, 2.0, 4, 4);
        histoIncRatioPurityTrueEff[0]->Draw();
        PlotLatexLegend(0.95, 0.96-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
        canvasIncRatio->Print(Form("%s/%s_%s_IncRatioPurity_trueEff.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));

        if(doBinShiftForDR){
          canvasIncRatio->cd();
          SetHistogramm(histoIncRatioPurityTrueEffBinShift, "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
          DrawGammaSetMarker(histoIncRatioPurityTrueEffBinShift, 20, 2.0, kGreen+2, kGreen+2);
          histoIncRatioPurityTrueEffBinShift->Draw();
          PlotLatexLegend(0.95, 0.96-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
          canvasIncRatio->Print(Form("%s/%s_%s_IncRatioPurity_trueEff_%s.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),addNameBinshift.Data(),suffix.Data()));

          // canvasIncRatio->cd();
          // SetHistogramm(histoIncRatioPurityTrueEffBinShift, "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
          // histoIncRatioPurityTrueEff[0]->Draw();
          // histoIncRatioPurityTrueEffBinShift->Draw("same");
          // PlotLatexLegend(0.95, 0.96-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
          // canvasIncRatio->Print(Form("%s/%s_%s_IncRatioPurity_trueEff_BinShiftComparison.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));

          canvasIncRatio->cd();
          histoIncRatioPurityTrueEffBinShiftCorr = (TH1D*)histoIncRatioPurityTrueEffBinShift->Clone("InclusiveRatioBinShiftCorrection");
          histoIncRatioPurityTrueEffBinShiftCorr->Divide(histoIncRatioPurityTrueEffBinShiftCorr,histoIncRatioPurityTrueEff[0],1,1,"");
          SetHistogramm(histoIncRatioPurityTrueEffBinShiftCorr, "#it{p}_{T} (GeV/#it{c})", "Shifted / Standard",0.8,1.2);
          DrawGammaSetMarker(histoIncRatioPurityTrueEffBinShiftCorr, 20, 2.0, kGreen+2, kGreen+2);
          histoIncRatioPurityTrueEffBinShiftCorr->Draw();
          PlotLatexLegend(0.95, 0.96-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
          canvasIncRatio->Print(Form("%s/%s_%s_IncRatioPurity_trueEff_BinShiftCorrection.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        }

        // Create ratio for MC
        canvasIncRatio->cd();
        SetHistogramm(histoMCIncRatio, "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
        DrawGammaSetMarker(histoMCIncRatio, 24, 2.0, 2, 2);
        histoMCIncRatio->Draw();
        PlotLatexLegend(0.95, 0.96-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
        canvasIncRatio->Print(Form("%s/%s_MC_IncRatio.%s",outputDir.Data(),nameOutputLabel.Data(),suffix.Data()));

        // Create plot with both ratios
        canvasIncRatio->cd();
        histoIncRatioPurityTrueEff[0]->Draw("e1");
        histoMCIncRatio->Draw("e1,same");
        TLegend* legendIncRatioAll              = GetAndSetLegend2(0.18,0.93-2*0.045,0.5,0.93,0.045,1,"",42,0.15);
        legendIncRatioAll->AddEntry(histoIncRatioPurityTrueEff[0],"Ratio with true Eff Purity");
        legendIncRatioAll->AddEntry(histoMCIncRatio,"MC Ratio");
        legendIncRatioAll->Draw();
        PlotLatexLegend(0.95, 0.96-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
        canvasIncRatio->Print(Form("%s/%s_%s_IncRatio_all.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        delete canvasIncRatio;
        cout << "Inclusive ratios plotted" << endl;
    }

    //**********************************************************************************
    //***                      Photon spectra data                                   ***
    //**********************************************************************************
    TCanvas* canvasGammaSpectraSingle           = GetAndSetCanvas("canvasGammaSpectraSingle", 0.12, 0.1, 1000 ,1350);
    DrawGammaCanvasSettings( canvasGammaSpectraSingle, 0.16, 0.02, 0.015, 0.07);
    canvasGammaSpectraSingle->SetLogy();
    canvasGammaSpectraSingle->SetLogx();

    textSizeSpectra=0.04;
    Double_t minPt                              = 0.2;
    if(mode==2) minPt                           = 0.5;
    if(mode==4) minPt                           = 0.8;
    Double_t maxPt                              = 50;
    Int_t nEntriesGammaSpec                     = 1;

    TH2F * histoDummy2                          = new TH2F("histoDummy2","histoDummy2",10000,minPt,maxPt,10000,5e-9, 90);
    SetStyleHistoTH2ForGraphs(histoDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                0.85*textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.77,1.78);
    histoDummy2->GetXaxis()->SetLabelOffset(-0.01);
    histoDummy2->Draw();
        DrawGammaSetMarker(histoGammaSpecCorrPurity, 20, 1.5,kBlue+2,kBlue+2);
        histoGammaSpecCorrPurity->GetYaxis()->SetRangeUser(3e-9, 10);

        if (cocktailAllGamma){
            nEntriesGammaSpec++;
            DrawGammaSetMarker(cocktailAllGamma, 2, 1.5,kBlack,kBlack);
            cocktailAllGamma->Draw("same,hist,l");
        }
        histoGammaSpecCorrPurity->DrawCopy("e1,same");

        TLegend* legendGammaDirSpectra          = GetAndSetLegend2(0.2,0.095,0.5,0.095+0.045*nEntriesGammaSpec,0.045,1,"",42,0.15);
        legendGammaDirSpectra->AddEntry(histoGammaSpecCorrPurity,"Inclusive photons");
        if (cocktailAllGamma) legendGammaDirSpectra->AddEntry(cocktailAllGamma, "decay cocktail","l");
        legendGammaDirSpectra->Draw();

        PlotLatexLegend(0.93, 0.98-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);

    canvasGammaSpectraSingle->Print(Form("%s/%s_%s_GammaSpectrum.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
    if(doBinShiftForDR){
      histoDummy2->Draw();
          DrawGammaSetMarker(histoGammaSpecCorrPurityBinShift, 20, 1.5,kGreen+2,kGreen+2);
          histoGammaSpecCorrPurityBinShift->GetYaxis()->SetRangeUser(3e-9, 10);

          if (cocktailAllGamma){
              nEntriesGammaSpec++;
              DrawGammaSetMarker(cocktailAllGamma, 2, 1.5,kBlack,kBlack);
              cocktailAllGamma->Draw("same,hist,l");
          }
          histoGammaSpecCorrPurityBinShift->DrawCopy("e1,same");
          legendGammaDirSpectra->Draw();
          PlotLatexLegend(0.93, 0.98-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);

      canvasGammaSpectraSingle->Print(Form("%s/%s_%s_GammaSpectrum_%s.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),addNameBinshift.Data(),suffix.Data()));
    }
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
            if(centrality.CompareTo("0-10%")==0){
                fitNameGammaA                   = "xqcd";
                fitNameGammaB                   = "rad";
                fitMaxPt                        = 14;
            } else{
                fitNameGammaA                   = "rad";
                fitNameGammaB                   = "oHag";
                fitMaxPt                        = 14;
            }
        }

        TString fitOptions                      = "QNRME+";
        if (!fEnergy.CompareTo("7TeV") && mode == 4) {
            fitOptions                          = "QNRME+";
        }

        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"---------------------- Begin Fitting Gamma ----------------------"<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;

        fitGammaA                               = FitObject(fitNameGammaA,"fitGammaA","Pi0",histoGammaSpecCorrPurity,fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(fitGammaA, 1, 2.0, kBlue-2);
        fileFinalResults << "IncRatioPurity_trueEff " << fitNameGammaA << endl;
        forOutput                               = WriteParameterToFile(fitGammaA);
        fileFinalResults << forOutput.Data() << endl;

        fitGammaB                               = FitObject(fitNameGammaB,"fitGammaB","Pi0",histoGammaSpecCorrPurity,fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(fitGammaB, 2, 2.0, kRed-3);
        fileFinalResults << "IncRatioPurity_trueEff " << fitNameGammaB << endl;
        forOutput                               = WriteParameterToFile(fitGammaB);
        fileFinalResults << forOutput.Data() << endl;

        cout<<"-----------------------------------------------------------------"<<endl;
        cout<<"------------------------ End Fitting Gamma ----------------------"<<endl;
        cout<<"-----------------------------------------------------------------"<<endl;

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
        if (mode==4)    ptMin                   = 0.8;
        Double_t        ptMax                   = 40.;
        histoDummy2->Draw();
            DrawGammaSetMarker(histoGammaSpecCorrPurity, 24, 2.0, kBlack, kBlack);
            fitGammaA->SetLineColor(kBlue-2);
            fitGammaB->SetLineColor(kRed-3);

            histoGammaSpecCorrPurity->Draw("same,e1");
            fitGammaA->Draw("same");
            fitGammaB->Draw("Csame");

            PlotLatexLegend(0.93, 0.95-0.045*2, 0.045,collisionSystem,detectionProcess,2,31);
            TLegend* legendGammaSpectraFits  = GetAndSetLegend2(0.7, 0.93-0.045*5, 0.9, 0.93-0.045*2, 0.045,1,"",42,0.2);
            legendGammaSpectraFits->AddEntry(histoGammaSpecCorrPurity,"#gamma data", "lp");
            legendGammaSpectraFits->AddEntry(fitGammaA,fitNameGammaA, "l");
            legendGammaSpectraFits->AddEntry(fitGammaB,fitNameGammaB, "l");
            legendGammaSpectraFits->Draw();

        // Lower part of plot with ratio
        padGammaFitsRatio->cd();
        padGammaFitsRatio->SetLogx();

            histoRatioFitGammaA                     = (TH1D*) histoGammaSpecCorrPurity->Clone("histoRatioFitGammaA");
            histoRatioFitGammaA                     = CalculateHistoRatioToFit(histoRatioFitGammaA,fitGammaA);
            histoRatioFitGammaB                     = (TH1D*) histoGammaSpecCorrPurity->Clone("histoRatioFitGammaB");
            histoRatioFitGammaB                     = CalculateHistoRatioToFit(histoRatioFitGammaB,fitGammaB);
            textSizeSpectra                         = 0.1;

            Double_t dummyMinX = 0.2;
            if(mode==2) dummyMinX = 0.5;
            if(mode==4) dummyMinX = 0.8;
            TH1D* dummy                             = new TH1D("dummy", "dummy",1000, dummyMinX, 20);
            SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})", "data/fit", textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.95,0.6);
            dummy->GetYaxis()->SetRangeUser(0.5, 1.55);
            dummy->GetXaxis()->SetLabelOffset(-0.02);

            DrawGammaSetMarker(histoRatioFitGammaA, 20, 1.5, kBlue-2, kBlue-2);
            DrawGammaSetMarker(histoRatioFitGammaB, 20, 1.5, kRed-3, kRed-3);

            dummy->Draw();
            histoRatioFitGammaA->Draw("e1,same");
            DrawGammaLines(dummyMinX, 20, 1.0, 1.0,1.0, kGray+2 ,7);
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
        for (Int_t i=1; i<=histoCorrectedPi0Yield[0]->GetNbinsX(); i++) {
            if (histoCorrectedPi0Yield[0]->GetBinContent(i) == 0) {
                continue;
            } else {
                fitMinPt                        = histoCorrectedPi0Yield[0]->GetXaxis()->GetBinLowEdge(i);
                break;
            }
        }
        fitMaxPt                                = histoCorrectedPi0Yield[0]->GetXaxis()->GetBinUpEdge(histoCorrectedPi0Yield[0]->GetNbinsX());

        // Special fitting fEnergys for PbPb
        if(fEnergy.CompareTo("PbPb_2.76TeV") == 0){
            if(centrality.CompareTo("0-10%")==0){
                fitPi0A                         = "oHag";
                fitPi0B                         = "rad";
                fitPi0C                         = "xqcd";
            } else{
                fitPi0A                         = "oHag";
                fitPi0B                         = "rad";
                fitPi0C                         = "qcd";
            }
        } else if (fEnergy.CompareTo("pPb_5.023TeV") == 0){
            fitPi0A                             = "h";
            fitPi0B                             = "l";
            fitPi0C                             = "oHag";
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
        fitPi0YieldA                            = FitObject(fitPi0A,"fitPi0YieldA","Pi0",histoCorrectedPi0Yield[0],fitMinPt,fitMaxPt,NULL,fitOptions);
        fitPi0YieldA->SetRange(fitMinPt,fitMaxPt);
        DrawGammaSetMarkerTF1(fitPi0YieldA, 1, 2.0, kBlue+1);
        fileFinalResults << "CorrectedYieldTrueEff " << fitPi0A << endl;
        forOutput                               = WriteParameterToFile(fitPi0YieldA);
        fileFinalResults << forOutput.Data() << endl;

        // Fitting Levy-Tsallis
        fitPi0YieldB                            = FitObject(fitPi0B,"fitPi0YieldB","Pi0",histoCorrectedPi0Yield[0],fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(fitPi0YieldB, 2, 2.0, kRed+1);
        fitPi0YieldB->SetRange(fitMinPt,fitMaxPt);
        fileFinalResults << "CorrectedYieldTrueEff " << fitPi0B << endl;
        forOutput                               = WriteParameterToFile(fitPi0YieldB);
        fileFinalResults << forOutput.Data() << endl;

        // Fitting modified Hagedorn
        fitPi0YieldC                            = FitObject(fitPi0C,"fitPi0YieldC","Pi0",histoCorrectedPi0Yield[0],fitMinPt,fitMaxPt,NULL,fitOptions);
        DrawGammaSetMarkerTF1(fitPi0YieldC, 3, 2.0, kGreen+1);
        fitPi0YieldC->SetRange(fitMinPt,fitMaxPt);
        fileFinalResults << "CorrectedYieldTrueEff " << fitPi0C << endl;
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
        if (mode==2)    ptMin                   = 0.5;
        if (mode==4)    ptMin                   = 0.8;
        Double_t        ptMax                   = 40.;
        TH2F * histoDummy3                      = new TH2F("histoDummy3","histoDummy3",1000,ptMin,1.5*histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()),1000,2e-9,999);
        SetStyleHistoTH2ForGraphs(histoDummy3, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                  0.85*textSizeSpectra,textSizeSpectra, textSizeSpectra,textSizeSpectra, 0.77,1.7);
        histoDummy3->DrawCopy();

            DrawGammaSetMarker(histoCorrectedPi0Yield[0], 20, 1.5, kBlack, kBlack);
            fitPi0YieldA->SetLineColor(kBlue-3);
            fitPi0YieldB->SetLineColor(kRed-3);
            fitPi0YieldC->SetLineColor(kGreen-2);

            histoCorrectedPi0Yield[0]->Draw("e1,same");
            fitPi0YieldA->Draw("Csame");
            fitPi0YieldB->Draw("Csame");
            fitPi0YieldC->Draw("same");

            PlotLatexLegend(0.93, 0.95-0.045*2, 0.045,collisionSystem,detectionProcess,2,31);

            TLegend* legendGammaSpectra       = GetAndSetLegend2(0.18, 0.05, 0.5, 0.05+0.045*4, 0.04,1,"",42,0.15);
            legendGammaSpectra->AddEntry(histoCorrectedPi0Yield[0],"#pi^{0} corr. yield","pl");
            legendGammaSpectra->AddEntry(fitPi0YieldA,Form("%s #chi^{2}/ndf %.2f",fitPi0A.Data(),fitPi0YieldA->GetChisquare()/fitPi0YieldA->GetNDF()),"l");
            legendGammaSpectra->AddEntry(fitPi0YieldB,Form("%s #chi^{2}/ndf %.2f",fitPi0B.Data(),fitPi0YieldB->GetChisquare()/fitPi0YieldB->GetNDF()),"l");
            legendGammaSpectra->AddEntry(fitPi0YieldC,Form("%s #chi^{2}/ndf %.2f",fitPi0C.Data(),fitPi0YieldC->GetChisquare()/fitPi0YieldC->GetNDF()),"l");
            legendGammaSpectra->Draw();

        // Plotting ratios of data and fits in lower pad
        padGammaRatios->cd();
        padGammaRatios->SetLogx();

            histoRatioFitPi0A                       = (TH1D*) histoCorrectedPi0Yield[0]->Clone("histoRatioFitPi0A");
            histoRatioFitPi0A                       = CalculateHistoRatioToFit(histoRatioFitPi0A,fitPi0YieldA);
            histoRatioFitPi0B                       = (TH1D*) histoCorrectedPi0Yield[0]->Clone("histoRatioFitPi0B");
            histoRatioFitPi0B                       = CalculateHistoRatioToFit(histoRatioFitPi0B,fitPi0YieldB);
            histoRatioFitPi0C                       = (TH1D*) histoCorrectedPi0Yield[0]->Clone("histoRatioFitPi0C");
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
        for (Int_t k = 0; k < 6; k++){
            histoIncRatioFitPurity[k]           = (TH1D*) histoIncRatioPurityTrueEff[k]->Clone(Form("histoIncRatioFitPurity%s",nameIntRanges[k].Data()));
            histoPi0Fitted[k]                   = (TH1D*) histoIncRatioPurityTrueEff[k]->Clone(Form("histoPi0Fitted%s",nameIntRanges[k].Data()));
            // Set datapoints in histoIncRatioFitPurity to the fit values from fitPi0YieldC (Tsallis fit)
            for(Int_t bin = 1; bin<histoIncRatioFitPurity[k]->GetNbinsX()+1; bin++){
                if(fEnergy.CompareTo("PbPb_2.76TeV") == 0){
                    if(histoIncRatioPurityTrueEff[k]->GetBinCenter(bin) <=3.5){
                        histoIncRatioFitPurity[k]->SetBinContent(bin,fitPi0YieldB->Eval(histoIncRatioPurityTrueEff[k]->GetBinCenter(bin)));
                        histoIncRatioFitPurity[k]->SetBinError(bin,histoCorrectedPi0Yield[k]->GetBinError(bin));
                        histoPi0Fitted[k]->SetBinContent(bin,fitPi0YieldB->Eval(histoPi0Fitted[k]->GetBinCenter(bin)));
                        histoPi0Fitted[k]->SetBinError(bin,histoCorrectedPi0Yield[k]->GetBinError(bin));
                    } else {
                        histoIncRatioFitPurity[k]->SetBinContent(bin,fitPi0YieldC->Eval(histoIncRatioPurityTrueEff[k]->GetBinCenter(bin)));
                        histoIncRatioFitPurity[k]->SetBinError(bin,histoCorrectedPi0Yield[k]->GetBinError(bin));
                        histoPi0Fitted[k]->SetBinContent(bin,fitPi0YieldC->Eval(histoPi0Fitted[k]->GetBinCenter(bin)));
                        histoPi0Fitted[k]->SetBinError(bin,histoCorrectedPi0Yield[k]->GetBinError(bin));
                    }
                } else {
                    histoIncRatioFitPurity[k]->SetBinContent(bin,fitPi0YieldC->Eval(histoIncRatioPurityTrueEff[k]->GetBinCenter(bin))); //nschmidt2016 changed to mod hagedorn
                    histoIncRatioFitPurity[k]->SetBinError(bin,histoCorrectedPi0Yield[k]->GetBinError(bin));
                    histoPi0Fitted[k]->SetBinContent(bin,fitPi0YieldC->Eval(histoPi0Fitted[k]->GetBinCenter(bin))); //nschmidt2016 changed to mod hagedorn
                    histoPi0Fitted[k]->SetBinError(bin,histoCorrectedPi0Yield[k]->GetBinError(bin));
                }
            }
            histoIncRatioFitPurity[k]->Divide( histoGammaSpecCorrPurity, histoIncRatioFitPurity[k]);
        }

        DrawGammaSetMarker(histoIncRatioPurityTrueEff[0], 20, 2.0, 1, 1);

        TCanvas* canvasIncRatioFit              = GetAndSetCanvas("canvasIncRatioFit",0.095,0.09,1000,815);
        SetHistogramm(histoIncRatioFitPurity[0], "#it{p}_{T} (GeV/#it{c})", "Ratio Inclusive #gamma/#pi^{0}",0,2);
        DrawGammaSetMarker(histoIncRatioFitPurity[0], 20, 2.0, kBlue-3, kBlue-3);

            histoIncRatioPurityTrueEff[0]->DrawCopy("e1");
            histoIncRatioFitPurity[0]->DrawCopy("same,e1");

            PlotLatexLegend(0.95, 0.95-2*0.045, 0.045,collisionSystem,detectionProcess,2,31);
            TLegend* legendIncRatioFit          = GetAndSetLegend2(0.7,0.93-4*0.045,0.93,0.93-2*0.045,0.045,1,"",42,0.15);
            legendIncRatioFit->AddEntry(histoIncRatioPurityTrueEff[0],"Inclusive Ratio");
            legendIncRatioFit->AddEntry(histoIncRatioFitPurity[0],"Ratio #pi^{0} Fit");
            legendIncRatioFit->Draw();
        canvasIncRatioFit->Print(Form("%s/%s_%s_IncRatio_Fit.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        delete canvasIncRatioFit;

        nameStatErrorCheckDat                       = Form("%s/%s/IncGammaToPi0Fit_%s_StatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
        fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
        fileStatErrorCheck << "# pT \t ratio \t absErr \t relErr(%)" << endl;
        for (Int_t i=1; i<=histoIncRatioFitPurity[0]->GetNbinsX(); i++){
            if (histoIncRatioFitPurity[0]->GetBinContent(i) != 0)
                fileStatErrorCheck << histoIncRatioFitPurity[0]->GetBinCenter(i) << "\t" << histoIncRatioFitPurity[0]->GetBinContent(i)<< "\t" << histoIncRatioFitPurity[0]->GetBinError(i) << "\t" << histoIncRatioFitPurity[0]->GetBinError(i)/histoIncRatioFitPurity[0]->GetBinContent(i)*100 << endl;
        }
        fileStatErrorCheck.close();

        nameStatErrorCheckDat                       = Form("%s/%s/IncGammaToPi0Fit_%s_RelStatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
        fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
        fileStatErrorCheck << "pT \t GammaStat \t Pi0Stat" << endl;
        for (Int_t i=1; i<=histoIncRatioFitPurity[0]->GetNbinsX(); i++){
            if (histoIncRatioFitPurity[0]->GetBinContent(i) != 0){
                Double_t relGammaErr    = histoGammaSpecCorrPurity->GetBinError(i)/histoGammaSpecCorrPurity->GetBinContent(i)*100 ;
                Double_t relPi0Err      = histoPi0Fitted[0]->GetBinError(i)/histoPi0Fitted[0]->GetBinContent(i)*100;
                Double_t relTotErr      = histoIncRatioFitPurity[0]->GetBinError(i)/histoIncRatioFitPurity[0]->GetBinContent(i)*100;
                fileStatErrorCheck << histoIncRatioFitPurity[0]->GetBinCenter(i) << "\t" << relGammaErr << "\t" << relPi0Err << "\t" << relTotErr << endl;
            }
        }
        fileStatErrorCheck.close();
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
        if (mode==2)    minPt                   = 0.5;
        if (mode==4)    minPt                   = 0.8;
        textSizeSpectra                         = 0.04;

        TH1F * histoDummy4                      = new TH1F("histoDummy4","histoDummy4",1000,minPt,1.2*histoGammaSpecCorrPurity->GetXaxis()->GetBinUpEdge(histoGammaSpecCorrPurity->GetNbinsX()));
        SetStyleHistoTH1ForGraphs(histoDummy4, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                  0.85*textSizeSpectra, textSizeSpectra, 0.85*textSizeSpectra, textSizeSpectra, 0.77,1.7);
        histoDummy4->GetXaxis()->SetLabelOffset(-0.02);
        histoDummy4->GetYaxis()->SetRangeUser(FindSmallestBin1DHist(histoCorrectedPi0Yield[0])/100.,FindLargestBin1DHist(histoCorrectedPi0Yield[0])*50.);
        histoDummy4->DrawCopy();

            DrawGammaSetMarker(cocktailPi0, 20, 2.0,kRed+2,kRed+2);
            cocktailPi0->DrawCopy("e1,same");

            fitPi0YieldC->SetLineColor(1);
            fitPi0YieldC->DrawCopy("same");

            DrawGammaSetMarker(histoCorrectedPi0Yield[0], 4, 2.0, 1, 1);
            histoCorrectedPi0Yield[0]->Draw("e1,same");

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

            PlotLatexLegend(0.93, 0.98-0.045*2, 0.045,collisionSystem,detectionProcess,2,31);

            TLegend* legendMesonSpectra                 = GetAndSetLegend2(0.19, 0.11, 0.45, 0.11+0.045*(3+nCocktails), 0.045,1,"",42,0.22);
            legendMesonSpectra->AddEntry(cocktailPi0,"Cocktail #pi^{0}");
            legendMesonSpectra->AddEntry(histoCorrectedPi0Yield[0],"#pi^{0} data");
            legendMesonSpectra->AddEntry(fitPi0YieldC,"fit to #pi^{0} data","l");
            if (haveEta){
                if(cocktailEta)legendMesonSpectra->AddEntry(cocktailEta,"Cocktail #eta");
                legendMesonSpectra->AddEntry(histoCorrectedEtaYield,"#eta data");
            }
            legendMesonSpectra->Draw();

        canvasMesonSpectra->Print(Form("%s/%s_%s_Cocktail_MesonSpectra.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        delete canvasMesonSpectra;
    }

    nameStatErrorCheckDat                       = Form("%s/%s/Pi0Cocktail_%s_StatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
    fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
    fileStatErrorCheck << "# pT \t cocktail \t absErr \t relErr(%)" << endl;
    for (Int_t i=1; i<=cocktailPi0->GetNbinsX(); i++){
        if (cocktailPi0->GetBinContent(i) != 0)
            fileStatErrorCheck << cocktailPi0->GetBinCenter(i) << "\t" << cocktailPi0->GetBinContent(i)<< "\t" << cocktailPi0->GetBinError(i) << "\t" << cocktailPi0->GetBinError(i)/cocktailPi0->GetBinContent(i)*100 << endl;
    }
    fileStatErrorCheck.close();

    //**********************************************************************************
    //***                      Double Ratios                                         ***
    //**********************************************************************************
    Bool_t doDoubleRatios                        = kTRUE;
    if (doDoubleRatios){
        cocktailAllGammaPi0                      = RebinTH1D(cocktailAllGammaPi0,histoIncRatioPurityTrueEff[0]);
        cocktailAllGammaRebinned                 = RebinTH1D(cocktailAllGamma,histoIncRatioPurityTrueEff[0]);
        cocktailPi0Rebinned                      = RebinTH1D(cocktailPi0,histoIncRatioPurityTrueEff[0]);
        cocktailAllGammaPi0->Divide(cocktailAllGammaRebinned,cocktailPi0Rebinned);

        for (Int_t k = 0; k < 6; k++){
            histoDoubleRatioTrueEffPurity[k]     = (TH1D*) histoIncRatioPurityTrueEff[k]->Clone(Form("DoubleRatioTrueEffPurity%s", nameIntRanges[k].Data()));
            histoDoubleRatioTrueEffPurity[k]->Divide(cocktailAllGammaPi0);
            if(doBinShiftForDR){
              histoDoubleRatioTrueEffPurityBinShift= (TH1D*) histoIncRatioPurityTrueEffBinShift->Clone(Form("DoubleRatioTrueEffPurity_%s", addNameBinshift.Data()));
              histoDoubleRatioTrueEffPurityBinShift->Divide(cocktailAllGammaPi0);
            }
            histoDoubleRatioFitPi0YieldPurity[k] = (TH1D*) histoIncRatioFitPurity[k]->Clone(Form("DoubleRatioFitPurity%s", nameIntRanges[k].Data()));
            histoDoubleRatioFitPi0YieldPurity[k]->Divide(cocktailAllGammaPi0);
        }

        // double ratio combined
        TCanvas *canvasDoubleRatio = new TCanvas("canvasDoubleRatio","",0.095,0.09,1000,815);
        DrawGammaCanvasSettings( canvasDoubleRatio, 0.09, 0.02, 0.02, 0.09);
        canvasDoubleRatio->cd();
        canvasDoubleRatio->SetLogx();

        Double_t minPt                           = 0.3;
        Double_t maxPtLog                        = 50.;
        Double_t minY                            = 0.85;
        Double_t maxY                            = 1.65;
        if (mode==2){
            minPt                                = 0.5;
            minY                                 = 0.75;
        }
        if (mode==4) {
            minPt                                = 0.8;
            minY                                 = 0.5;
            if(fEnergy.CompareTo("8TeV")==0 || fEnergy.CompareTo("7TeV")==0 || fEnergy.CompareTo("900GeV")==0){
              minPt = 1.0;
              minY = 0.7;
            }
        }

        TString      method                      = "PCM";
        if (mode==4) method                      = "EMCal";

        TH2F * histo2DDoubleRatioPlotting       = new TH2F("histo2DDoubleRatioPlotting","histo2DDoubleRatioPlotting",1000,minPt,maxPtLog,1000,minY,maxY);
        SetStyleHistoTH2ForGraphs(histo2DDoubleRatioPlotting, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.04,0.04, 0.04,0.04, 1.,1.);
        histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.,histoDoubleRatioTrueEffPurity[0]->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity[0]->GetNbinsX())*1.5);
        histo2DDoubleRatioPlotting->GetXaxis()->SetTitleFont(62);
        histo2DDoubleRatioPlotting->GetYaxis()->SetTitleFont(62);
        histo2DDoubleRatioPlotting->GetXaxis()->SetLabelOffset(-1e-2);
        histo2DDoubleRatioPlotting->DrawCopy();

            DrawGammaSetMarker(histoDoubleRatioTrueEffPurity[0], 20, 2.0, kBlack, kBlack);
            DrawGammaSetMarker(histoDoubleRatioFitPi0YieldPurity[0], 20, 2.0, kBlue+2, kBlue+2);

            histoDoubleRatioFitPi0YieldPurity[0]->DrawCopy("same,E1");
            DrawGammaLines(0., histoDoubleRatioTrueEffPurity[0]->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity[0]->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
            histoDoubleRatioTrueEffPurity[0]->DrawCopy("same,E1");

            TLegend* legendDoubleConversionFit       = GetAndSetLegend2(0.14,0.93-2*0.045,0.5,0.93,0.045,1,"",42,0.2);
            legendDoubleConversionFit->AddEntry(histoDoubleRatioTrueEffPurity[0],"Data","p");
            legendDoubleConversionFit->AddEntry(histoDoubleRatioFitPi0YieldPurity[0],"Data, fitted #pi^{0}","p");
            legendDoubleConversionFit->Draw();

            PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);

        canvasDoubleRatio->Print(Form("%s/%s_%s_DoubleRatioComparison.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));

        // double ratio bin shift comparison
        if(doBinShiftForDR){
          canvasDoubleRatio->cd();
          histo2DDoubleRatioPlotting->DrawCopy();
              DrawGammaSetMarker(histoDoubleRatioTrueEffPurityBinShift, 20, 2.0, kGreen+2, kGreen+2);

              DrawGammaLines(0., 20. , 1.0, 1.0,2.0, kGray+2 ,7);
              histoDoubleRatioTrueEffPurity[0]->DrawCopy("same,E1");
              histoDoubleRatioTrueEffPurityBinShift->DrawCopy("same,E1");

              legendDoubleConversionFit->Clear();
              legendDoubleConversionFit->AddEntry(histoDoubleRatioTrueEffPurity[0],"Data","p");
              legendDoubleConversionFit->AddEntry(histoDoubleRatioTrueEffPurityBinShift,"Data Bin Shifted","p");
              legendDoubleConversionFit->Draw();
              PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);

          canvasDoubleRatio->Print(Form("%s/%s_%s_DoubleRatioBinShiftComparison.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        }

        // double ratio fit
        canvasDoubleRatio->cd();
        histo2DDoubleRatioPlotting->DrawCopy();

            DrawGammaLines(0., 20. , 1.0, 1.0,2.0, kGray+2 ,7);
            histoDoubleRatioFitPi0YieldPurity[0]->DrawCopy("same,E1");

            TLegend* legendDoubleConversionFit2       = GetAndSetLegend2(0.14,0.93-1*0.045,0.5,0.93,0.045,1,"",42,0.2);
            legendDoubleConversionFit2->AddEntry(histoDoubleRatioFitPi0YieldPurity[0],"Data, fitted #pi^{0}","p");
            legendDoubleConversionFit2->Draw();
            PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);

        canvasDoubleRatio->Print(Form("%s/%s_%s_DoubleRatioFit.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));

        // double ratio no fit
        canvasDoubleRatio->cd();
        histo2DDoubleRatioPlotting->DrawCopy();

            DrawGammaLines(0., histoDoubleRatioTrueEffPurity[0]->GetXaxis()->GetBinUpEdge(histoDoubleRatioTrueEffPurity[0]->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
            histoDoubleRatioTrueEffPurity[0]->DrawCopy("same,E1");
            PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);

            TF1* fitRGamma          = new TF1("fitRGamma","[0]",1,3);
            TF1* fitRGammaFitPi0    = new TF1("fitRGammaFitPi0","[0]",1,3);
            histoDoubleRatioFitPi0YieldPurity[0]->Fit(fitRGammaFitPi0,"QRME0","",1,3);
            histoDoubleRatioTrueEffPurity[0]->Fit(fitRGamma,"QRME0","",1,3);
            fileFinalResults << "std.: " << fitRGamma->GetParameter(0) << "+-" << fitRGamma->GetParError(0) << endl;
            fileFinalResults << "fitted.: " << fitRGammaFitPi0->GetParameter(0) << "+-" << fitRGammaFitPi0->GetParError(0) << endl;

            TLegend* legendDoubleConversionFit3       = GetAndSetLegend2(0.14,0.93-1*0.045,0.5,0.93,0.045,1,"",42,0.2);
            legendDoubleConversionFit3->AddEntry(histoDoubleRatioTrueEffPurity[0],"Data","p");
            legendDoubleConversionFit3->Draw();

        canvasDoubleRatio->Print(Form("%s/%s_%s_DoubleRatio.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
        delete canvasDoubleRatio;

        nameStatErrorCheckDat                       = Form("%s/%s/RGamma_%s_StatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
        fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
        fileStatErrorCheck << "# pT \t DR \t absErr \t relErr(%)" << endl;
        for (Int_t i=1; i<=histoDoubleRatioTrueEffPurity[0]->GetNbinsX(); i++){
            if (histoDoubleRatioTrueEffPurity[0]->GetBinContent(i) != 0)
                fileStatErrorCheck << histoDoubleRatioTrueEffPurity[0]->GetBinCenter(i) << "\t" << histoDoubleRatioTrueEffPurity[0]->GetBinContent(i)<< "\t" << histoDoubleRatioTrueEffPurity[0]->GetBinError(i) << "\t" << histoDoubleRatioTrueEffPurity[0]->GetBinError(i)/histoDoubleRatioTrueEffPurity[0]->GetBinContent(i)*100 << endl;
        }
        fileStatErrorCheck.close();

        nameStatErrorCheckDat                       = Form("%s/%s/RGammaFitPi0_%s_StatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
        fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
        fileStatErrorCheck << "# pT \t DR \t absErr \t relErr(%)" << endl;
        for (Int_t i=1; i<=histoDoubleRatioFitPi0YieldPurity[0]->GetNbinsX(); i++){
            if (histoDoubleRatioFitPi0YieldPurity[0]->GetBinContent(i) != 0)
                fileStatErrorCheck << histoDoubleRatioFitPi0YieldPurity[0]->GetBinCenter(i) << "\t" << histoDoubleRatioFitPi0YieldPurity[0]->GetBinContent(i)<< "\t" << histoDoubleRatioFitPi0YieldPurity[0]->GetBinError(i) << "\t" << histoDoubleRatioFitPi0YieldPurity[0]->GetBinError(i)/histoDoubleRatioFitPi0YieldPurity[0]->GetBinContent(i)*100 << endl;
        }
        fileStatErrorCheck.close();

        nameStatErrorCheckDat                       = Form("%s/%s/RGamma_%s_RelStatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
        fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
        fileStatErrorCheck << "pT \t GammaStat \t Pi0Stat \t CocktailStat" << endl;
        for (Int_t i=1; i<=histoDoubleRatioTrueEffPurity[0]->GetNbinsX(); i++){
            if (histoDoubleRatioTrueEffPurity[0]->GetBinContent(i) != 0){
                Double_t relGammaErr    = histoGammaSpecCorrPurity->GetBinError(i)/histoGammaSpecCorrPurity->GetBinContent(i)*100 ;
                Double_t relPi0Err      = histoCorrectedPi0Yield[0]->GetBinError(i)/histoCorrectedPi0Yield[0]->GetBinContent(i)*100;
                Double_t relCocktErr    = cocktailPi0->GetBinError(i)/cocktailPi0->GetBinContent(i)*100;
                Double_t relTotErr      = histoDoubleRatioTrueEffPurity[0]->GetBinError(i)/histoDoubleRatioTrueEffPurity[0]->GetBinContent(i)*100;
                fileStatErrorCheck << histoDoubleRatioTrueEffPurity[0]->GetBinCenter(i) << "\t" << relGammaErr << "\t" << relPi0Err << "\t" << relCocktErr << "\t"<< relTotErr << endl;
            }
        }
        fileStatErrorCheck.close();

        nameStatErrorCheckDat                       = Form("%s/%s/RGammaFitPi0_%s_RelStatError.dat",cutSel.Data(),fEnergy.Data(), nameRec.Data());
        fileStatErrorCheck.open(nameStatErrorCheckDat.Data(), ios::out);
        fileStatErrorCheck << "pT \t GammaStat \t Pi0Stat \t CocktailStat" << endl;
        for (Int_t i=1; i<=histoDoubleRatioTrueEffPurity[0]->GetNbinsX(); i++){
            if (histoDoubleRatioFitPi0YieldPurity[0]->GetBinContent(i) != 0){
                Double_t relGammaErr    = histoGammaSpecCorrPurity->GetBinError(i)/histoGammaSpecCorrPurity->GetBinContent(i)*100 ;
                Double_t relPi0Err      = histoPi0Fitted[0]->GetBinError(i)/histoPi0Fitted[0]->GetBinContent(i)*100;
                Double_t relCocktErr    = cocktailPi0->GetBinError(i)/cocktailPi0->GetBinContent(i)*100;
                Double_t relTotErr      = histoDoubleRatioFitPi0YieldPurity[0]->GetBinError(i)/histoDoubleRatioFitPi0YieldPurity[0]->GetBinContent(i)*100;
                fileStatErrorCheck << histoDoubleRatioFitPi0YieldPurity[0]->GetBinCenter(i) << "\t" << relGammaErr << "\t" << relPi0Err << "\t" << relCocktErr << "\t"<< relTotErr << endl;
            }
        }
        fileStatErrorCheck.close();

    }

    fileFinalResults.close();

    //**********************************************************************************
    //***                 Inclusive and decay ratio                                  ***
    //**********************************************************************************
    TCanvas* canvasIncRatioDecRatio         = GetAndSetCanvas("canvasIncRatioDecRatio",0.095,0.09,1000,815);
    canvasIncRatioDecRatio->SetLogx();

    SetHistogramm(cocktailAllGammaPi0,          "#it{p}_{T} (GeV/#it{c})", "#gamma/#pi^{0}",0,1.8);
    SetHistogramm(histoIncRatioPurityTrueEff[0],   "#it{p}_{T} (GeV/#it{c})", "#gamma/#pi^{0}",0,1.8);

    DrawGammaSetMarker(histoIncRatioPurityTrueEff[0],  20, 2.0, kBlack,  kBlack);
    DrawGammaSetMarker(cocktailAllGammaPi0,         33, 2.0, kBlue-2,  kBlue-2);

    histoIncRatioPurityTrueEff[0]->GetXaxis()->SetRangeUser(0.2, 20);
    if(mode==2) histoIncRatioPurityTrueEff[0]->GetXaxis()->SetRangeUser(0.5, 20);
    if(mode==4) histoIncRatioPurityTrueEff[0]->GetXaxis()->SetRangeUser(0.8, 20);
    histoIncRatioPurityTrueEff[0]->Draw();
    cocktailAllGammaPi0->Draw("same");

    PlotLatexLegend(0.95, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
    TLegend* legendIncRatioDecRatio          = GetAndSetLegend2(0.7,0.93-6*0.045,0.93,0.93-3*0.045,0.045,1,"",42,0.15);
    legendIncRatioDecRatio->AddEntry(histoIncRatioPurityTrueEff[0],    "(#gamma_{incl.} / #pi^{0})_{meas.}", "lp");
    legendIncRatioDecRatio->AddEntry(cocktailAllGammaPi0,           "(#gamma_{dec.} / #pi^{0})_{sim.}",  "lp");
    legendIncRatioDecRatio->Draw();

    canvasIncRatioDecRatio->Print(Form("%s/%s_%s_IncRatio_DecRatio_trueEff.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),suffix.Data()));
    if(doBinShiftForDR){
      canvasIncRatioDecRatio->cd();

      SetHistogramm(histoIncRatioPurityTrueEffBinShift,   "#it{p}_{T} (GeV/#it{c})", "#gamma/#pi^{0}",0,1.8);

      DrawGammaSetMarker(histoIncRatioPurityTrueEffBinShift,  20, 2.0, kGreen+2,  kGreen+2);

      histoIncRatioPurityTrueEffBinShift->GetXaxis()->SetRangeUser(0.2, 20);
      if(mode==2) histoIncRatioPurityTrueEffBinShift->GetXaxis()->SetRangeUser(0.5, 20);
      if(mode==4) histoIncRatioPurityTrueEffBinShift->GetXaxis()->SetRangeUser(0.8, 20);
      histoIncRatioPurityTrueEffBinShift->Draw();
      cocktailAllGammaPi0->Draw("same");

      PlotLatexLegend(0.95, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
      TLegend* legendIncRatioDecRatioBS          = GetAndSetLegend2(0.7,0.93-6*0.045,0.93,0.93-3*0.045,0.045,1,"",42,0.15);
      legendIncRatioDecRatioBS->AddEntry(histoIncRatioPurityTrueEffBinShift,    "(#gamma_{incl.} / #pi^{0})_{meas.}^{bin. shift.}", "lp");
      legendIncRatioDecRatioBS->AddEntry(cocktailAllGammaPi0,           "(#gamma_{dec.} / #pi^{0})_{sim.}",  "lp");
      legendIncRatioDecRatioBS->Draw();

      canvasIncRatioDecRatio->Print(Form("%s/%s_%s_IncRatio_DecRatio_trueEff_%s.%s",outputDir.Data(),nameOutputLabel.Data(),nameRec.Data(),addNameBinshift.Data(),suffix.Data()));
    }
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
        for (Int_t k = 0; k< 6; k++){
            histoCorrectedPi0Yield[k]->Write( histoCorrectedPi0Yield[k]->GetName(),  TObject::kOverwrite);
        }
        if(doBinShiftForDR && histoCorrectedPi0YieldBinShift)histoCorrectedPi0YieldBinShift->Write( histoCorrectedPi0YieldBinShift->GetName(),  TObject::kOverwrite);
        // pi0 syst related to OOB pileup
        graphSysPi0PileUpOptions->Write("Pi0_SystErrorRel_OOBPileup_Options",  TObject::kOverwrite);
        graphSysPi0PileUpIterations->Write("Pi0_SystErrorRel_OOBPileup_Iterations",  TObject::kOverwrite);

        // gamma quantities
        if (fitGammaA)                  fitGammaA->Write(                   fitGammaA->GetName(),                   TObject::kOverwrite);
        if (histoGammaSpecCorrPurity)   histoGammaSpecCorrPurity->Write(    "histoGammaSpecCorrPurity",             TObject::kOverwrite);
        if (doBinShiftForDR){
          if (histoGammaSpecCorrPurityBinShift)   histoGammaSpecCorrPurityBinShift->Write(    Form("histoGammaSpecCorrPurity_%s",addNameBinshift.Data()), TObject::kOverwrite);
          if (histoGammaSpecCorrPurityBinShiftCorr)   histoGammaSpecCorrPurityBinShiftCorr->Write(    "histoGammaSpecCorrPurityBinShiftCorr",             TObject::kOverwrite);
        }

        // Double ratio
        if(histoDoubleRatioUpperLimits)             histoDoubleRatioUpperLimits->Write(             Form("%s_UpperLimits", histoDoubleRatioTrueEffPurity[0]->GetName()),  TObject::kOverwrite);
        for (Int_t k = 0; k< 6; k++){
            if(histoDoubleRatioTrueEffPurity[k])    histoDoubleRatioTrueEffPurity[k]->Write(        histoDoubleRatioTrueEffPurity[k]->GetName(),        TObject::kOverwrite);
            if(histoDoubleRatioFitPi0YieldPurity[k])histoDoubleRatioFitPi0YieldPurity[k]->Write(    histoDoubleRatioFitPi0YieldPurity[k]->GetName(),    TObject::kOverwrite);
        }
        if(doBinShiftForDR && histoDoubleRatioTrueEffPurityBinShift)    histoDoubleRatioTrueEffPurityBinShift->Write(        histoDoubleRatioTrueEffPurityBinShift->GetName(),        TObject::kOverwrite);
        // inclusive ratio
        for (Int_t k = 0; k< 6; k++){
            histoIncRatioPurityTrueEff[k]->Write(  histoIncRatioPurityTrueEff[k]->GetName(),  TObject::kOverwrite);
            histoIncRatioFitPurity[k]->Write(      histoIncRatioFitPurity[k]->GetName(),      TObject::kOverwrite);
        }
        if(doBinShiftForDR){
          if (histoIncRatioPurityTrueEffBinShift) histoIncRatioPurityTrueEffBinShift->Write(  histoIncRatioPurityTrueEffBinShift->GetName(),  TObject::kOverwrite);
          if (histoIncRatioPurityTrueEffBinShiftCorr) histoIncRatioPurityTrueEffBinShiftCorr->Write(  histoIncRatioPurityTrueEffBinShiftCorr->GetName(),  TObject::kOverwrite);
        }
        histoMCIncRatio->Write(             histoMCIncRatio->GetName(),             TObject::kOverwrite);

        // cocktail
        cocktailAllGamma->Write(cocktailAllGamma->GetName(),    TObject::kOverwrite);

        TString particlesInCocktail[14]   = {"Pi0","Eta","EtaPrime","omega","rho0","rho+","rho-","phi","Delta0","Delta+","Sigma0","K0s","K0l","Lambda"};
        TH1D* dummyCocktailHist;
        TH1D* dummyGammaCocktailHist;
        for(Int_t i=0; i<14;i++){
            dummyCocktailHist             = NULL;
            dummyGammaCocktailHist        = NULL;
            dummyCocktailHist             = (TH1D* )cocktailFile->Get(Form("%s_Pt",particlesInCocktail[i].Data()));
            dummyGammaCocktailHist        = (TH1D* )cocktailFile->Get(Form("Gamma_From_%s_Pt",particlesInCocktail[i].Data()));
            if(dummyCocktailHist) dummyCocktailHist->Write(     dummyCocktailHist->GetName(),         TObject::kOverwrite);
            if(dummyGammaCocktailHist) dummyGammaCocktailHist->Write(     dummyGammaCocktailHist->GetName(),         TObject::kOverwrite);
        }

        fileCorrectedOutput->Close();
        delete fileCorrectedOutput;
    }
}
