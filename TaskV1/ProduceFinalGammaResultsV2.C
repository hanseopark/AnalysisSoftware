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
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

void ProduceFinalGammaResultsV2(    TString configurationFileName   = "configurationFileName.txt",
                                    TString optionEnergy            = "",
                                    TString suffix                  = "eps",
                                    Int_t mode                      = 0
                                ){

    StyleSettingsThesis();
    SetPlotStyle();


    ifstream in(configurationFileName.Data());
    cout<<"Available Triggers:"<<endl;
    TString fCutNumber                  = "";
    TString nameSysFiles[3]             = {""};
    Int_t doSys                         = 0;
    Bool_t doSysGamma                   = 0;
    Bool_t doSysIncGammaToPi0           = 0;
    Bool_t doSysDR                      = 0;
    Bool_t doTheo                       = kFALSE;
    TString nameTheoryFile              = "";
    TString nameTheoryGraphsErr         = "";
    TString nameTheoryGraphsCenter      = "";
    Bool_t doDirGamma                   = kFALSE;
    Int_t lines                         = 0;
    while(!in.eof() && lines < 1){
        TString doTheoString            = "";
        TString doDirGammaString        = "";
        in >> fCutNumber >> nameSysFiles[0] >>nameSysFiles[1] >> nameSysFiles[2]    >> doTheoString    >> nameTheoryFile >> nameTheoryGraphsErr >> nameTheoryGraphsCenter >> doDirGammaString;
        cout << fCutNumber.Data() << endl;
        if ( nameSysFiles[0].CompareTo("bla") != 0 ){
            doSys                       = doSys+100;
            doSysGamma                  = kTRUE;
        }
        if ( nameSysFiles[1].CompareTo("bla") != 0 ){
            doSys                       = doSys+10;
            doSysIncGammaToPi0          = kTRUE;
        }
        if ( nameSysFiles[2].CompareTo("bla") != 0 ){
            doSys                       = doSys+1;
            doSysDR                     = kTRUE;
        }
        if (doSys == 111)
            cout << "will do full sys: \t" << nameSysFiles[0].Data() <<  "\t" << nameSysFiles[1].Data() << "\t" << nameSysFiles[2].Data() << endl;
        else if (doSys > 0)
            cout << "will do part of sys: \t" << nameSysFiles[0].Data() <<  "\t" << nameSysFiles[1].Data() << "\t" << nameSysFiles[2].Data() << endl;
        else
            cout << "no sys has been defined" <<endl;
        cout<< "theo: " << doTheoString << endl;
        if (doTheoString.CompareTo("kTRUE") == 0)
            doTheo                      = kTRUE;
        if (doTheo){
            cout << "trying to read theo" << endl;
            cout << "path: " << nameTheoryFile.Data() << endl;
            cout << "path: " << nameTheoryGraphsErr.Data() << endl;
            cout << "path: " << nameTheoryGraphsCenter.Data() << endl;
            if (nameTheoryFile.CompareTo("bla") == 0){
                doTheo      = kFALSE;
                cout << "theory file not correctly specified" << endl;
            }
        }
        if (doDirGammaString.CompareTo("kTRUE") == 0)
            doDirGamma                  = kTRUE;
        if (doDirGamma)
            cout << "will calculate dir gamma spectra" << endl;
        lines++;
    }

    TString dateForOutput           = ReturnDateStringForOutput();
    TString collisionSystem         = ReturnFullCollisionsSystem(optionEnergy);
    TString collisionSystemOutput   = ReturnCollisionEnergyOutputString(optionEnergy);
    TString centrality              = GetCentralityString(fCutNumber);
    TString centralityW0Per         = GetCentralityStringWoPer(fCutNumber);
    if (centrality.CompareTo("pp") == 0 || centrality.CompareTo("0-100%") == 0 ) {
        centrality                  = "";
        centralityW0Per             = "";
    }
    TString detectionProcess        = ReturnFullTextReconstructionProcess(mode);

    TString system                  = "PCM";
    if (mode == 2) system           = "PCM-EMC";
    if (mode == 3) system           = "PCM-PHOS";
    if (mode == 4) system           = "EMC";
    if (mode == 5) system           = "PHOS";
    if (mode == 10) system          = "mEMC";
    if (mode == 11) system          = "mPHOS";

    // defining output directory
    TString outputDir       = Form("%s/%s/FinalGammaResults_%s_%s%s", suffix.Data(), dateForOutput.Data(), system.Data(),centralityW0Per.Data(), collisionSystemOutput.Data() );
    gSystem->Exec("mkdir -p "+outputDir);



    //*************************************************************************************************
    //******************** read from input files ******************************************************
    //*************************************************************************************************
    TString inputFileName               = Form("%s/%s/Gamma_Pi0_data_GammaConvV1_InclusiveRatio.root", fCutNumber.Data(), optionEnergy.Data());
    cout << "trying to read: " << inputFileName.Data() << endl;
    TFile *fileInput                    = new TFile(inputFileName.Data());
    if (fileInput->IsZombie()) {
        cout << "file couldn't be read, aborting....";
        return;
    }

    // read stat error hists
    TH1D* histoIncGamma                     = (TH1D*) fileInput->Get("histoGammaSpecCorrPurity");
    TH1D* histoPi0Spectrum                  = (TH1D*) fileInput->Get("CorrectedYieldTrueEff");
    TH1D* histoIncRatio                     = (TH1D*) fileInput->Get("IncRatioPurity_trueEff");
    TH1D* histoIncRatioPi0Fit               = (TH1D*) fileInput->Get("histoIncRatioFitPurity");
    TH1D* histoDR                           = (TH1D*) fileInput->Get("DoubleRatioTrueEffPurity");
    TH1D* histoDRFit                        = (TH1D*) fileInput->Get("DoubleRatioFitPurity");
    TH1D* histococktailAllGamma             = (TH1D*) fileInput->Get("Gamma_Pt");
    
    TString inputFileNameAdditional               = Form("%s/%s/Gamma_Pi0_data_GammaConvV1Correction_%s.root", fCutNumber.Data(), optionEnergy.Data(), fCutNumber.Data());
    cout << "trying to read: " << inputFileNameAdditional.Data() << endl;
    TFile *fileInputAdditional                    = new TFile(inputFileNameAdditional.Data());
    if (fileInputAdditional->IsZombie()) {
        cout << "file couldn't be read, aborting....";
        return;
    }

    // read stat error hists
    TH1D* histoPileupCorrection             = (TH1D*) fileInputAdditional->Get("PileUpCorrectionFactor");
    TH1D* histoGammaRawYields               = (TH1D*) fileInputAdditional->Get("GammaRaw_Pt");
    
    
    TString inputFileNameAdditional2               = Form("%s/%s/Pi0_MC_GammaConvV1CorrectionHistos_%s.root", fCutNumber.Data(), optionEnergy.Data(), fCutNumber.Data());
    cout << "trying to read: " << inputFileNameAdditional2.Data() << endl;
    TFile *fileInputAdditional2                    = new TFile(inputFileNameAdditional2.Data());
    if (fileInputAdditional2->IsZombie()) {
        cout << "file couldn't be read, aborting....";
        return;
    }

    // read stat error hists
    TH1D* histoGammaPurity                  = (TH1D*) fileInputAdditional2->Get("GammaTruePurity_Pt");
    TH1D* histoGammaConvProb                = (TH1D*) fileInputAdditional2->Get("GammaConvProb_MCPt");
    TH1D* histoGammaRecoEff                 = (TH1D*) fileInputAdditional2->Get("GammaRecoEff_Pt");


    // calculate sys error graphs
    TGraphAsymmErrors* graphIncGammaSysErr      =  NULL;
    TGraphAsymmErrors* graphIncRatioSysErr      =  NULL;
    TGraphAsymmErrors* graphIncRatioPi0FitSysErr=  NULL;
    TGraphAsymmErrors* graphDRSysErr            =  NULL;
    TGraphAsymmErrors* graphDRPi0FitSysErr      =  NULL;

    if (doSysGamma){
        ifstream fileSysErrGamma;
        Int_t nPointsGamma                                          = 0;
        Double_t ptSysGamma[50];
        Double_t relSystErrorGammaUp[50];
        Double_t relSystErrorGammaDown[50];
        Double_t relSystErrorWOMaterialGammaUp[50];
        Double_t relSystErrorWOMaterialGammaDown[50];
        fileSysErrGamma.open(nameSysFiles[0].Data(),ios_base::in);
        cout << nameSysFiles[0].Data() << endl;

        while(!fileSysErrGamma.eof() && nPointsGamma < 100){
            fileSysErrGamma >> ptSysGamma[nPointsGamma] >> relSystErrorGammaDown[nPointsGamma] >> relSystErrorGammaUp[nPointsGamma]>>    relSystErrorWOMaterialGammaDown[nPointsGamma] >> relSystErrorWOMaterialGammaUp[nPointsGamma];
            cout << nPointsGamma << "\t"  << relSystErrorGammaDown[nPointsGamma] << "\t"  <<relSystErrorGammaUp[nPointsGamma] << "\t" << relSystErrorWOMaterialGammaDown[nPointsGamma] << "\t"  <<relSystErrorWOMaterialGammaUp[nPointsGamma] << endl;;
            nPointsGamma++;
        }
        fileSysErrGamma.close();
        nPointsGamma                = nPointsGamma-1;
        graphIncGammaSysErr         = CalculateSysErrFromRelSysHistoWithPtBins( histoIncGamma, "GammaSystError",
                                                                            relSystErrorGammaDown, relSystErrorGammaUp, ptSysGamma, nPointsGamma);
    }

    if (doSysIncGammaToPi0) {
        ifstream fileSysErrInclRatio;
        Int_t nPointsInclRatio                                      = 0;
        Double_t ptSysInclRatioFit[50];
        Double_t ptSysInclRatio[50];
        Double_t relSystErrorInclRatioUp[50];
        Double_t relSystErrorInclRatioDown[50];
        Double_t relSystErrorWOMaterialInclRatioUp[50];
        Double_t relSystErrorWOMaterialInclRatioDown[50];

        fileSysErrInclRatio.open(nameSysFiles[1].Data(),ios_base::in);
        cout << nameSysFiles[1].Data() << endl;

        while(!fileSysErrInclRatio.eof() && nPointsInclRatio < 100){
            fileSysErrInclRatio >> ptSysInclRatio[nPointsInclRatio] >> relSystErrorInclRatioDown[nPointsInclRatio] >> relSystErrorInclRatioUp[nPointsInclRatio]>>    relSystErrorWOMaterialInclRatioDown[nPointsInclRatio] >> relSystErrorWOMaterialInclRatioUp[nPointsInclRatio];
            cout << nPointsInclRatio << "\t"  << relSystErrorInclRatioDown[nPointsInclRatio] << "\t"  <<relSystErrorInclRatioUp[nPointsInclRatio] << "\t" << relSystErrorWOMaterialInclRatioDown[nPointsInclRatio] << "\t"  <<relSystErrorWOMaterialInclRatioUp[nPointsInclRatio] << endl;;
            nPointsInclRatio++;
        }
        fileSysErrInclRatio.close();
        nPointsInclRatio            = nPointsInclRatio-1;
        graphIncRatioSysErr         = CalculateSysErrFromRelSysHistoWithPtBins( histoIncRatio, "InclRatioSysError",
                                                                                relSystErrorInclRatioDown, relSystErrorInclRatioUp, ptSysInclRatio, nPointsInclRatio);
        graphIncRatioPi0FitSysErr   = CalculateSysErrFromRelSysHistoWithPtBins( histoIncRatioPi0Fit, "InclRatioSysError",
                                                                                relSystErrorInclRatioDown, relSystErrorInclRatioUp, ptSysInclRatio, nPointsInclRatio);
    }

    if (doSysDR){
        ifstream fileSysErrDoubleRatio;
        Int_t nPointsDoubleRatio                                    = 0;
        Double_t ptSysDoubleRatio[50];
        Double_t relSystErrorDoubleRatioUp[50];
        Double_t relSystErrorDoubleRatioDown[50];
        Double_t relSystErrorWOMaterialDoubleRatioUp[50];
        Double_t relSystErrorWOMaterialDoubleRatioDown[50];

        fileSysErrDoubleRatio.open(nameSysFiles[2].Data(),ios_base::in);
        cout << nameSysFiles[2].Data() << endl;

        while(!fileSysErrDoubleRatio.eof() && nPointsDoubleRatio < 100){
            fileSysErrDoubleRatio >> ptSysDoubleRatio[nPointsDoubleRatio] >> relSystErrorDoubleRatioDown[nPointsDoubleRatio] >> relSystErrorDoubleRatioUp[nPointsDoubleRatio]>>    relSystErrorWOMaterialDoubleRatioDown[nPointsDoubleRatio] >> relSystErrorWOMaterialDoubleRatioUp[nPointsDoubleRatio];
            cout << nPointsDoubleRatio << "\t"  << relSystErrorDoubleRatioDown[nPointsDoubleRatio] << "\t"  <<relSystErrorDoubleRatioUp[nPointsDoubleRatio] << "\t" << relSystErrorWOMaterialDoubleRatioDown[nPointsDoubleRatio] << "\t"  <<relSystErrorWOMaterialDoubleRatioUp[nPointsDoubleRatio] << endl;;
            nPointsDoubleRatio++;
        }
        fileSysErrDoubleRatio.close();
        nPointsDoubleRatio          = nPointsDoubleRatio-1;

        graphDRSysErr               = CalculateSysErrFromRelSysHistoWithPtBins( histoDR, "DRSysError",
                                                                                relSystErrorDoubleRatioDown, relSystErrorDoubleRatioUp, ptSysDoubleRatio, nPointsDoubleRatio);
        graphDRPi0FitSysErr         = CalculateSysErrFromRelSysHistoWithPtBins( histoDRFit, "DRSysError",
                                                                                relSystErrorDoubleRatioDown, relSystErrorDoubleRatioUp, ptSysDoubleRatio, nPointsDoubleRatio);

    }

    // read theory graphs
    TGraphAsymmErrors* graphNLODR           = NULL;
    TGraph* graphNLODRCenter                = NULL;
    if (doTheo ){
        TFile* fileTheo                     = new TFile(nameTheoryFile.Data());
        if (fileTheo->IsZombie()) {
            cout << "couldn't read theo file, jumping....";
            doTheo                          = kFALSE;
        } else {
            if (nameTheoryGraphsErr.CompareTo("bla") != 0){
                graphNLODR                  = (TGraphAsymmErrors*)fileTheo->Get(nameTheoryGraphsErr.Data());
            }
            if (nameTheoryGraphsCenter.CompareTo("bla") != 0){
                graphNLODRCenter            = (TGraph*)fileTheo->Get(nameTheoryGraphsCenter.Data());
            }
            if (!graphNLODR && !graphNLODRCenter)
                doTheo                      = kFALSE;

        }

    }
    //***************************************************************************************************
    //***************************************************************************************************
    //***************************************************************************************************
    Int_t nLinesNLOLegends  = 2;
    if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0 || optionEnergy.CompareTo("pPb_5.023TeV") == 0)
        nLinesNLOLegends    = 3;
    Double_t minPt                              = 0.2;
    Double_t maxPtLog                        = 40.;
    Double_t minY                            = 0.75;
    Double_t maxY                            = 1.65;
    if(mode == 0 && optionEnergy.CompareTo("PbPb_2.76TeV") == 0){
        maxY                                 = 2.;
        minY                                 = 0.75;
    } else if (mode==4) {
        minPt                                = 1.5;
        minY                                 = 0.5;
    }

    TCanvas *canvasDoublRatioNLO = new TCanvas("canvasDoublRatioNLO","",0.095,0.09,1000,815);
    DrawGammaCanvasSettings( canvasDoublRatioNLO, 0.09, 0.02, 0.02, 0.09);
    canvasDoublRatioNLO->cd();
    canvasDoublRatioNLO->SetLogx();

    TH2F * histo2DDoubleRatioPlotting       = new TH2F("histo2DDoubleRatioPlotting","histo2DDoubleRatioPlotting",1000,minPt,maxPtLog,1000,minY,maxY);
    SetStyleHistoTH2ForGraphs(histo2DDoubleRatioPlotting, "#it{p}_{T} (GeV/#it{c})","R_{#gamma}", 0.04,0.04, 0.04,0.04, 1.,1.);
    if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0)
        histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.4,histoDR->GetXaxis()->GetBinUpEdge(histoDR->GetNbinsX())*1.5);
    else
        histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.,histoDR->GetXaxis()->GetBinUpEdge(histoDR->GetNbinsX())*1.5);
    histo2DDoubleRatioPlotting->GetXaxis()->SetTitleFont(62);
    histo2DDoubleRatioPlotting->GetYaxis()->SetTitleFont(62);
    histo2DDoubleRatioPlotting->GetXaxis()->SetLabelOffset(-1e-2);
    histo2DDoubleRatioPlotting->DrawCopy();

        histo2DDoubleRatioPlotting->DrawCopy();
            if (doTheo && graphNLODR) {
                DrawGammaSetMarkerTGraphAsym(graphNLODR, 0, 0, kAzure-9, kAzure-9, 0.2, kTRUE, kAzure-9);
                graphNLODR->Draw("3,same");
            }
            if (doTheo && graphNLODRCenter){
                DrawGammaNLOTGraph( graphNLODRCenter, 2, 7, kAzure+2);
                graphNLODRCenter->Draw("lc,same");
            }
            if (graphDRPi0FitSysErr){
                DrawGammaSetMarkerTGraphAsym(graphDRPi0FitSysErr,26,0,kBlue-8,kBlue-8,2,kTRUE);
                graphDRPi0FitSysErr->Draw("p,2,same");
            }
            if (graphDRSysErr){
                DrawGammaSetMarkerTGraphAsym(graphDRSysErr,26,0,kGray+2,kGray+2,2,kTRUE);
                graphDRSysErr->Draw("p,2,same");
            }
            DrawGammaSetMarker(histoDR, 20, 2.0, kBlack, kBlack);
            DrawGammaSetMarker(histoDRFit, 20, 2.0, kBlue+2, kBlue+2);
            DrawGammaLines(0., histoDR->GetXaxis()->GetBinUpEdge(histoDR->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);

            histoDR->DrawCopy("same");
            histoDRFit->DrawCopy("same");

            TLegend* legendDoubleConversionFitNLO = GetAndSetLegend2(0.14,0.93-(nLinesNLOLegends+1+0.5)*0.045,0.5,0.93,0.045,1,"",42,0.2);
            legendDoubleConversionFitNLO->AddEntry(histoDR,"Data","p");
            legendDoubleConversionFitNLO->AddEntry(histoDRFit,"Data, fitted #pi^{0}","p");
            if (doTheo && graphNLODR){
                if(optionEnergy.CompareTo("PbPb_2.76TeV") == 0)       legendDoubleConversionFitNLO->AddEntry(graphNLODR,"pp #gamma_{dir} NLO #sqrt{s} = 2.76 TeV","f");
                else if(optionEnergy.CompareTo("pPb_5.023TeV") == 0)  legendDoubleConversionFitNLO->AddEntry(graphNLODR,"pp #gamma_{dir} NLO #sqrt{s} = 5.02 TeV ","f");
                else                                             legendDoubleConversionFitNLO->AddEntry(graphNLODR,"pp NLO Direct Photon","f");
                if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0 || optionEnergy.CompareTo("pPb_5.023TeV") == 0) legendDoubleConversionFitNLO->AddEntry((TObject*)0,"scaled by N_{coll}","");
            }
            legendDoubleConversionFitNLO->Draw();

        histo2DDoubleRatioPlotting->Draw("same,axis");
    canvasDoublRatioNLO->Print(Form("%s/DoubleRatioComparison_NLO.%s",outputDir.Data(),suffix.Data()));

    //*************************************************************************************************
    // calculate the direct photon spectrum
    //*************************************************************************************************
    TH1D* histoDirectPhotonSpectrum = NULL;
    TH1D* histoDoubleRatioUpperLimits= NULL;
    TH1D* histoDirectPhotonSpectrumFit= NULL;
    TH1D* histoDoubleRatioUpperLimitsFit= NULL;
    if(doDirGamma){
      cout << __LINE__ << endl;
        // calculation for non-fitted pi0 double ratio
        if(graphDRSysErr&&histoIncGamma&&histoDR){
          cout << __LINE__ << endl;
            histoDirectPhotonSpectrum            = (TH1D*)histoIncGamma->Clone("histoDirectPhotonSpectrum");
            histoDoubleRatioUpperLimits          = GetUpperLimitsHisto(histoDR,  graphDRSysErr, 0.95, 1e-6, 1e4);
            Double_t binContent                  = -1;
            Double_t binError                    = -1;
            for (Int_t i=1; i<histoDirectPhotonSpectrum->GetNbinsX()+1; i++) {
                if (!histoDirectPhotonSpectrum->GetBinContent(i)) continue;
                binContent                       = histoDirectPhotonSpectrum->GetBinContent(i) * (1 - 1/histoDoubleRatioUpperLimits->GetBinContent(i));
                binError                         = 0.;
                histoDirectPhotonSpectrum->SetBinContent(i, binContent);
                histoDirectPhotonSpectrum->SetBinError(  i, binError);
            }
            cout << __LINE__ << endl;
        }
        // calculation for fitted pi0 double ratio
        if(graphDRPi0FitSysErr&&histoIncGamma&&histoDRFit){
          cout << __LINE__ << endl;
            histoDirectPhotonSpectrumFit        = (TH1D*)histoIncGamma->Clone("histoDirectPhotonSpectrumFit");
            histoDoubleRatioUpperLimitsFit      = GetUpperLimitsHisto(histoDRFit,  graphDRPi0FitSysErr, 0.95, 1e-6, 1e4);
            Double_t binContent                  = -1;
            Double_t binError                    = -1;
            for (Int_t i=1; i<histoDirectPhotonSpectrumFit->GetNbinsX()+1; i++) {
                if (!histoDirectPhotonSpectrumFit->GetBinContent(i)) continue;
                binContent                      = histoDirectPhotonSpectrumFit->GetBinContent(i) * (1 - 1/histoDoubleRatioUpperLimitsFit->GetBinContent(i));
                binError                        = 0.;
                histoDirectPhotonSpectrumFit->SetBinContent(i, binContent);
                histoDirectPhotonSpectrumFit->SetBinError(  i, binError);
            }
            cout << __LINE__ << endl;
        }
    }

    //*************************************************************************************************
    // put everything in common output per system
    //*************************************************************************************************
    TString optionOutput                = "pp";
    TString fileNameOutputComp          = Form("%s_%sResultsFullCorrection_PP.root","data",system.Data());
    if (optionEnergy.Contains("pPb")){
        fileNameOutputComp              = Form("%s_%sResultsFullCorrection_pPb.root","data",system.Data());
        optionOutput                    = "pPb";
    } else if (optionEnergy.Contains("PbPb")){
        fileNameOutputComp              = Form("%s_%sResultsFullCorrection_PbPb.root","data",system.Data());
        optionOutput                    = "PbPb";
    }
    TFile* fileGammaFinal               = new TFile(fileNameOutputComp,"UPDATE");

        // create subdirectory for respective energy
        TDirectoryFile* directoryGamma      = (TDirectoryFile*)fileGammaFinal->Get(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() ));
        if (!directoryGamma){
            fileGammaFinal->mkdir(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() ));
            directoryGamma      = (TDirectoryFile*)fileGammaFinal->Get(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() ));
            fileGammaFinal->cd(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() ));
        } else {
            fileGammaFinal->cd(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() ));
        }

            // writing double ratio quantities
            if (histoDR){
                SetHistogramm(histoDR,"#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})");
                histoDR->Write("DoubleRatioStatError",TObject::kOverwrite);
            }
            if (graphDRSysErr) graphDRSysErr->Write("DoubleRatioSystError",TObject::kOverwrite);

            if (histoDRFit){
                SetHistogramm(histoDRFit,"#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})");
                histoDRFit->Write("DoubleRatioPi0FitStatError",TObject::kOverwrite);
            }
            if(graphDRPi0FitSysErr) graphDRPi0FitSysErr->Write("DoubleRatioPi0FitSystError",TObject::kOverwrite);

            // writing inclusive ratio quantities
            if (histoIncRatio){
                SetHistogramm(histoIncRatio,"#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
                histoIncRatio->Write("IncRatioStatError",TObject::kOverwrite);
            }
            if (graphIncRatioSysErr) graphIncRatioSysErr->Write("IncRatioSystError",TObject::kOverwrite);

            if (histoIncRatioPi0Fit){
                SetHistogramm(histoIncRatioPi0Fit,"#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
                histoIncRatioPi0Fit->Write("IncRatioPi0FitStatError",TObject::kOverwrite);
            }
            if(graphIncRatioPi0FitSysErr) graphIncRatioPi0FitSysErr->Write("IncRatioPi0FitSystError",TObject::kOverwrite);

            // writing inclusive gamma spectrum
            if (histoIncGamma){
                SetHistogramm(histoIncGamma,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#gamma}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
                histoIncGamma->Write("IncGammaStatError",TObject::kOverwrite);
            }
            if (graphIncGammaSysErr) graphIncGammaSysErr->Write("IncGammaSystError",TObject::kOverwrite);

            // writing pi0 used pi0 spectrum
            if (histoPi0Spectrum){
                SetHistogramm(histoPi0Spectrum,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
                histoPi0Spectrum->Write("Pi0StatError",TObject::kOverwrite);
            }
            // // writing direct photon spectrum
            if(histoDirectPhotonSpectrum)            histoDirectPhotonSpectrum->Write("DirectPhotonSpectrum",TObject::kOverwrite);
            if(histoDirectPhotonSpectrumFit)            histoDirectPhotonSpectrumFit->Write("DirectPhotonSpectrumFit",TObject::kOverwrite);


          if(histoPileupCorrection)            histoPileupCorrection->Write("PileUpCorrectionFactor",TObject::kOverwrite);
          if(histoGammaRawYields)            histoGammaRawYields->Write("GammaRawYields",TObject::kOverwrite);
          if(histoGammaPurity)            histoGammaPurity->Write("GammaTruePurity",TObject::kOverwrite);
          if(histoGammaConvProb)            histoGammaConvProb->Write("GammaConversionProbability",TObject::kOverwrite);
          if(histoGammaRecoEff)            histoGammaRecoEff->Write("GammaRecoEfficiency",TObject::kOverwrite);

      // Cocktail histograms
      directoryGamma->mkdir("Cocktail");
      TDirectoryFile* directoryCocktail = (TDirectoryFile*)directoryGamma->Get("Cocktail");
      directoryGamma->cd("Cocktail");
      if(histococktailAllGamma) histococktailAllGamma->Write(histococktailAllGamma->GetName(),TObject::kOverwrite);
      TString particlesInCocktail[14]   = {"Pi0","Eta","EtaPrim","omega","rho0","rho+","rho-","phi","Delta0","Delta+","Sigma0","K0s","K0l","Lambda"};
      TH1D* dummyCocktailHist;
      TH1D* dummyGammaCocktailHist;
      TH1D* dummyGammaRatioCocktailHist;
      for(Int_t i=0; i<14;i++){
          dummyCocktailHist             = NULL;
          dummyGammaCocktailHist        = NULL;
          dummyGammaRatioCocktailHist   = NULL;
          dummyCocktailHist             = (TH1D* )fileInput->Get(Form("%s_Pt",particlesInCocktail[i].Data()));
          dummyGammaCocktailHist        = (TH1D* )fileInput->Get(Form("Gamma_From_%s_Pt",particlesInCocktail[i].Data()));
          if(dummyCocktailHist) dummyCocktailHist->Write(dummyCocktailHist->GetName(),TObject::kOverwrite);
          if(dummyGammaCocktailHist){
              dummyGammaCocktailHist->Write(dummyGammaCocktailHist->GetName(),TObject::kOverwrite);
              if(histococktailAllGamma){
                  dummyGammaRatioCocktailHist = (TH1D* )fileInput->Get(Form("Gamma_From_%s_Pt",particlesInCocktail[i].Data()));
                  dummyGammaRatioCocktailHist->Divide(histococktailAllGamma);
                  dummyGammaRatioCocktailHist->Write(Form("Gamma_From_%s_Ratio_To_All",particlesInCocktail[i].Data()), TObject::kOverwrite);
              }
          }
      }

    fileGammaFinal->Write();
    fileGammaFinal->Close();

}

