/****************************************************************************************************************************
******      Friederike Bock, friederike.bock@cern.ch                                      *****
******      Hikari Murakami, Pedro Gonzalez
******      Efficiency study 
******      LHC12f1a,LHC12f1,bLHC12i3, LHC15g1a (ancLHC11a-WSDD)
******
*****************************************************************************************************************************/

#include <Riostream.h>
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
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

extern TRandom* gRandom;
extern TBenchmark* gBenchmark;
extern TSystem* gSystem;
extern TMinuit* gMinuit;

//**********************************************************************************************************************
//***************************** Main function **************************************************************************
//**********************************************************************************************************************
void ExtractInputForWeights(    TString suffix                  = "pdf", 
                                TString nameFileData            = "data_PCMResults_pp.root", 
                                Bool_t combinedSpectrum         = kFALSE,
                                Int_t nGenerators               = 0,
                                TString configFile              = "",
                                TString fOptEnergy              = "2.76TeV",
                                Bool_t runDrawReweighted        = kTRUE 
                            ){

    //**********************************************************************************************************************
    //*************************************************** general style settings *******************************************
    //**********************************************************************************************************************
    gROOT->Reset();   
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();  
    SetPlotStyle();
    
    TString dateForOutput           = ReturnDateStringForOutput();
    TString outputDir               = Form("%s/%s/ExtractInputForWeights",suffix.Data(),dateForOutput.Data());
        
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputFileData.root ",nameFileData.Data(),outputDir.Data() ));
        
    ifstream in(configFile.Data());
    if (nGenerators > 5) 
        nGenerators                 = 5;
    Int_t nrOfGenerators            = 0;
    TString generatorName[5]        = {"", "", "", "", ""};
    TString fileNameMCPi0[5]        = {"", "", "", "", ""};
    TString fileNameMCEta[5]        = {"", "", "", "", ""};
    TString outputString[5]         = {"", "", "", "", ""};
    
    while(!in.eof() && nrOfGenerators<nGenerators ){
        // file format: generatorName "\t" filename for pi0 input file as output from CorrectSignal.C "\t" filename for eta input file as output from CorrectSignal.C "\t" generalized string to be written to the weighting file
        in >> generatorName[nrOfGenerators] >> fileNameMCPi0[nrOfGenerators] >> fileNameMCEta[nrOfGenerators] >> outputString[nrOfGenerators];
        cout<< generatorName[nrOfGenerators] << "\t"<<fileNameMCPi0[nrOfGenerators]<< "\t" << fileNameMCEta[nrOfGenerators] << "\t"<< outputString[nrOfGenerators] << endl;
        nrOfGenerators++;
    }

    TString collisionSystem                         = ReturnFullCollisionsSystem(fOptEnergy);
    TString collisionSystemForWriting               = ReturnCollisionEnergyStringForTheory(fOptEnergy);
    Color_t colorData                               = GetColorDefaultColor( fOptEnergy, "", "", kFALSE);
    Style_t markerStyleSpectrum                     = GetDefaultMarkerStyle( fOptEnergy, "", "");
    Size_t  markerSizeSpectrum                      = GetDefaultMarkerSize( fOptEnergy, "", "");
    
    Color_t colorMC[5]                              = { 1, 1, 1, 1, 1 };
    for (Int_t i = 0; i< nGenerators; i++){
        colorMC[i]                                  = GetColorDefaultColor( fOptEnergy, generatorName[i], "", kFALSE);
    }    
    
    Style_t lineStyleMC[5]                          = { 1, 2, 3, 4, 5} ;    
    Style_t markerStyleMC[5]                        = { 20, 25, 30, 29, 32};

    TString nameHistoPi0                            = "CorrectedYieldPi0";
    TString nameHistoEta                            = "CorrectedYieldEta";
    TString nameFitPi0                              = "TwoComponentModelFitYieldPi0";
    TString nameFitEta                              = "TwoComponentModelFitYieldEta";
    if (combinedSpectrum){
        nameHistoPi0                                = "graphInvYieldINELPi0Comb2760GeVATotErr";
        nameHistoEta                                = "graphInvYieldINELEtaComb2760GeVATotErr";
    }
                
    //**********************************************************************************************************************
    //************************************** read reconstructed data *******************************************************
    //**********************************************************************************************************************
    TFile* fileDataInput                            = new TFile(nameFileData);
    
    TDirectory* directoryPi0                        = (TDirectory*)fileDataInput->Get(Form("Pi0%s",fOptEnergy.Data())); 
    TDirectory* directoryEta                        = (TDirectory*)fileDataInput->Get(Form("Eta%s",fOptEnergy.Data()));
    
    TGraphAsymmErrors* graphYieldPi0                = NULL;
    TGraphAsymmErrors* graphYieldEta                = NULL;
    TF1* fitPi0Yield                                = NULL;
    TF1* fitEtaYield                                = NULL;
    if (!combinedSpectrum){
        TH1D* histoYieldPi0                         = (TH1D*)directoryPi0->Get(nameHistoPi0.Data());
        TH1D* histoYieldEta                         = (TH1D*)directoryEta->Get(nameHistoEta.Data());
        graphYieldPi0                               = new TGraphAsymmErrors(histoYieldPi0);
        graphYieldEta                               = new TGraphAsymmErrors(histoYieldEta);
        // remove empty bins of spectrum
        while (graphYieldPi0->GetY()[0] == 0) graphYieldPi0->RemovePoint(0);
        while (graphYieldEta->GetY()[0] == 0) graphYieldEta->RemovePoint(0);
    } else {
        graphYieldPi0                               = (TGraphAsymmErrors*)directoryPi0->Get(nameHistoPi0.Data());
        graphYieldEta                               = (TGraphAsymmErrors*)directoryEta->Get(nameHistoEta.Data());
        fitPi0Yield                                 = (TF1*)directoryPi0->Get(nameFitPi0.Data());
        fitEtaYield                                 = (TF1*)directoryEta->Get(nameFitEta.Data());
    }
    
    //**********************************************************************************************************************
    //****************************************Fit Pi0 Spectra **************************************************************
    //**********************************************************************************************************************    
    if (!fitPi0Yield){
        fitPi0Yield     = FitObject("l","fitPi0Yield","Pi0");
        graphYieldPi0->Fit(fitPi0Yield,"QNRMEI+","",graphYieldPi0->GetX()[0]-graphYieldPi0->GetErrorXlow(0), graphYieldPi0->GetX()[graphYieldPi0->GetN()-1]-graphYieldPi0->GetErrorXlow(graphYieldPi0->GetN()-1));
    }    
    // print fit result to shell
    cout << WriteParameterToFile(fitPi0Yield)<< endl;	

    //**********************************************************************************************************************
    //****************************************Fit Eta Spectra **************************************************************
    //**********************************************************************************************************************    
    if (!fitEtaYield){
        fitEtaYield     = FitObject("l","fitEtaYield","Eta");
        graphYieldEta->Fit(fitEtaYield,"QNRMEI+","",graphYieldEta->GetX()[0]-graphYieldEta->GetErrorXlow(0), graphYieldEta->GetX()[graphYieldEta->GetN()-1]-graphYieldEta->GetErrorXlow(graphYieldEta->GetN()-1));
    }    
    // print fit result to shell
    cout << WriteParameterToFile(fitEtaYield)<< endl;	
    
    // *********************************************************************************************************************
    // ******************************* Read inputs from different MC's *****************************************************
    // *********************************************************************************************************************
    TFile* filePi0MCInput[5]                        = {NULL, NULL, NULL, NULL, NULL};
    TFile* fileEtaMCInput[5]                        = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoPi0InputMCWOWeights[5]               = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoPi0InputMCWWeights[5]                = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoPi0Efficiency[5]                     = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaInputMCWOWeights[5]               = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaInputMCWWeights[5]                = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaEfficiency[5]                     = {NULL, NULL, NULL, NULL, NULL};
    
    for (Int_t i = 0; i < nGenerators; i++){
        cout << fileNameMCPi0[i].Data() << "\t" << fileNameMCEta[i].Data() << endl;
        filePi0MCInput[i]                           = new TFile(fileNameMCPi0[i]);
        histoPi0InputMCWOWeights[i]                 = (TH1D*)filePi0MCInput[i]->Get("MCYield_Meson_oldBinWOWeights");
        histoPi0InputMCWWeights[i]                 = (TH1D*)filePi0MCInput[i]->Get("MCYield_Meson_oldBin");
        
        histoPi0Efficiency[i]                       = (TH1D*)filePi0MCInput[i]->Get("TrueMesonEffiPt");
        fileEtaMCInput[i]                           = new TFile(fileNameMCEta[i]);
        histoEtaInputMCWOWeights[i]                 = (TH1D*)fileEtaMCInput[i]->Get("MCYield_Meson_oldBinWOWeights");
        histoEtaInputMCWWeights[i]                  = (TH1D*)fileEtaMCInput[i]->Get("MCYield_Meson_oldBin");
        histoEtaEfficiency[i]                       = (TH1D*)fileEtaMCInput[i]->Get("TrueMesonEffiPt");
    }    
    
    //**********************************************************************************************************************
    //****************************************Pi0 Spectra compared to MC****************************************************
    //**********************************************************************************************************************
    TCanvas* canvasSpectra  = new TCanvas("canvasSpectra", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasSpectra,  0.13, 0.01, 0.015, 0.08);
    canvasSpectra->SetLogy();
    canvasSpectra->SetLogx();
    
    TH1F * histo1DSpectra;
    histo1DSpectra          = new TH1F("histo1DSpectra", "histo1DSpectra",1000, 0.23, graphYieldPi0->GetX()[graphYieldPi0->GetN()-1]*2);
    SetStyleHistoTH1ForGraphs( histo1DSpectra, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
                            0.03, 0.035, 0.03, 0.035, 0.83, 1.7);
    histo1DSpectra->GetYaxis()->SetRangeUser(graphYieldPi0->GetY()[graphYieldPi0->GetN()-1]/10, graphYieldPi0->GetX()[0]*10);
    histo1DSpectra->GetXaxis()->SetLabelOffset(-0.01);
    histo1DSpectra->GetYaxis()->SetLabelOffset(0.01);
    histo1DSpectra->DrawCopy(); 

        TLegend* legendSpectraPi0            = GetAndSetLegend2(0.16, 0.12, 0.73, 0.12+0.035*(1+nGenerators)*1.15, 0.035, 1, "", 42, 0.15);
        DrawGammaSetMarkerTGraphAsym(graphYieldPi0, markerStyleSpectrum, markerSizeSpectrum, colorData , colorData);
        graphYieldPi0->Draw("p,same,e1");
        legendSpectraPi0->AddEntry(graphYieldPi0,collisionSystem.Data(),"pe");
        
        for (Int_t i = 0; i < nGenerators; i++){
            SetStyleHisto(histoPi0InputMCWOWeights[i], 3., lineStyleMC[i], colorMC[i]);  
            histoPi0InputMCWOWeights[i]->Draw("same,hist,c");
            legendSpectraPi0->AddEntry(histoPi0InputMCWOWeights[i],generatorName[i],"l");
        }    
        graphYieldPi0->Draw("p,same,e1");
        
        TLatex *labelSpectraPi0Label        = new TLatex(0.95,0.92,"#pi^{0} #rightarrow #gamma #gamma");
        SetStyleTLatex( labelSpectraPi0Label, 0.035, 4, 1, 42, kTRUE, 31);
        labelSpectraPi0Label->Draw();
        
        legendSpectraPi0->Draw();
    
    canvasSpectra->Update();
    canvasSpectra->Print(Form("%s/Pi0_Spectra_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));
    
    if (histoPi0InputMCWWeights[0]){
        histo1DSpectra->DrawCopy(); 

            graphYieldPi0->Draw("p,same,e1");
            
            for (Int_t i = 0; i < nGenerators; i++){
                if (histoPi0InputMCWWeights[i]){
                    SetStyleHisto(histoPi0InputMCWWeights[i], 3., lineStyleMC[i], colorMC[i]);  
                    histoPi0InputMCWWeights[i]->Draw("same,hist,c");
                }    
            }    
            graphYieldPi0->Draw("p,same,e1");
            labelSpectraPi0Label->Draw();
            legendSpectraPi0->Draw();
        
        canvasSpectra->Update();
        canvasSpectra->Print(Form("%s/Pi0_Spectra_MCWeighted_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));
    
    }
    
    // plot result with input
    canvasSpectra->cd();
    histo1DSpectra->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphYieldPi0, markerStyleSpectrum, markerSizeSpectrum, colorData , colorData);
        graphYieldPi0->Draw("p,same,e1");

        fitPi0Yield->SetLineColor(colorData);
        fitPi0Yield->Draw("same");
    
        labelSpectraPi0Label->Draw();
    
    canvasSpectra->Print(Form("%s/Pi0_Spectra_WithFit_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ******************************* Ratio of data to fit and MC input to fit for pi0 2.76 TeV ****************************
    // **********************************************************************************************************************
    // definition of canvas for ratios
    TCanvas* canvasRatioToFit = new TCanvas("canvasRatioToFit","",1550,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToFit,  0.08, 0.015, 0.015, 0.07);
    canvasRatioToFit->SetGridx(0);
    canvasRatioToFit->SetGridy(0);
    canvasRatioToFit->SetLogx(1);

    // Calculation of ratio histograms
    TGraphAsymmErrors* graphRatioPi0DatatoFit   = CalculateGraphErrRatioToFit (graphYieldPi0, fitPi0Yield);
    TH1D* histoPi0RatioMCInputsToFit[5]         = {NULL, NULL, NULL, NULL, NULL};
    for (Int_t i = 0; i < nGenerators; i++){
        histoPi0RatioMCInputsToFit[i]           = CalculateHistoRatioToFit (histoPi0InputMCWOWeights[i], fitPi0Yield,kTRUE);
        SetStyleHisto(histoPi0RatioMCInputsToFit[i], 2, lineStyleMC[i], colorMC[i]);
    }

    
    // drawing
    canvasRatioToFit->cd();    
    TH1F * histo1DRatio;
    histo1DRatio          = new TH1F("histo1DRatio", "histo1DRatio",1000, 0.23, graphYieldPi0->GetX()[graphYieldPi0->GetN()-1]*2);
    SetStyleHistoTH1ForGraphs( histo1DRatio, "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
                            0.03, 0.035, 0.03, 0.035, 0.83, 1.2);
    histo1DRatio->GetYaxis()->SetRangeUser(0, 3);
    histo1DRatio->GetXaxis()->SetLabelOffset(-0.01);
    histo1DRatio->GetYaxis()->SetLabelOffset(0.01);
    histo1DRatio->DrawCopy(); 
    
        TLegend* legendRatioPi0 =  GetAndSetLegend2(0.11, 0.12, 0.4, 0.12+0.035*(1+nGenerators)*1.15, 0.035, 1, "", 42, 0.2); 
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0DatatoFit, markerStyleSpectrum, markerSizeSpectrum, colorData , colorData);
        graphRatioPi0DatatoFit->Draw("e,p");  
        legendRatioPi0->AddEntry(graphRatioPi0DatatoFit,"Data/fit to Data","pe");
        
        for (Int_t i = 0; i < nGenerators; i++){
            histoPi0RatioMCInputsToFit[i]->Draw("same,hist,l");
            legendRatioPi0->AddEntry(histoPi0RatioMCInputsToFit[i],Form("%s/ fit to Data", generatorName[i].Data()),"l");
        }
        graphRatioPi0DatatoFit->Draw("e,p");  
        legendRatioPi0->Draw();

        TLatex *labelEnergy2760GeVRatio     = new TLatex(0.12,0.92,collisionSystem.Data());
        SetStyleTLatex( labelEnergy2760GeVRatio, 0.035 , 4, 1, 42, kTRUE, 11);
        labelEnergy2760GeVRatio->Draw();    
        TLatex *labelSpectraPi0LabelRatio   = new TLatex(0.12,0.88,"#pi^{0} #rightarrow #gamma #gamma");
        SetStyleTLatex( labelSpectraPi0LabelRatio, 0.035, 4, 1, 42, kTRUE, 11);
        labelSpectraPi0LabelRatio->Draw();
    
        DrawGammaLines(0.23, graphYieldPi0->GetX()[graphYieldPi0->GetN()-1]*2 ,1., 1.,0.1);
    canvasRatioToFit->Update();
    canvasRatioToFit->SaveAs(Form("%s/Pi0_RatioToDataFit_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));


    //**********************************************************************************************************************
    //**************************************** Eta Spectra compared to MC****************************************************
    //**********************************************************************************************************************    
    canvasSpectra->cd();
    histo1DSpectra->DrawCopy(); 

        TLegend* legendSpectraEta            = GetAndSetLegend2(0.16, 0.12, 0.73, 0.12+0.035*(1+nGenerators)*1.15, 0.035, 1, "", 42, 0.15);
        DrawGammaSetMarkerTGraphAsym(graphYieldEta, markerStyleSpectrum, markerSizeSpectrum, colorData , colorData);
        graphYieldEta->Draw("p,same,e1");
        legendSpectraEta->AddEntry(graphYieldEta,collisionSystem.Data(),"pe");
        
        for (Int_t i = 0; i < nGenerators; i++){
            SetStyleHisto(histoEtaInputMCWOWeights[i], 3., lineStyleMC[i], colorMC[i]);  
            histoEtaInputMCWOWeights[i]->Draw("same,hist,c");
            legendSpectraEta->AddEntry(histoEtaInputMCWOWeights[i],generatorName[i],"l");
        }    
        graphYieldEta->Draw("p,same,e1");
        
        TLatex *labelSpectraEtaLabel        = new TLatex(0.95,0.92,"#eta #rightarrow #gamma #gamma");
        SetStyleTLatex( labelSpectraEtaLabel, 0.035, 4, 1, 42, kTRUE, 31);
        labelSpectraEtaLabel->Draw();
        
        legendSpectraEta->Draw();
    
    canvasSpectra->Update();
    canvasSpectra->Print(Form("%s/Eta_Spectra_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));

    if (histoEtaInputMCWWeights[0]){
        histo1DSpectra->DrawCopy(); 

            graphYieldEta->Draw("p,same,e1");
            legendSpectraEta->AddEntry(graphYieldEta,collisionSystem.Data(),"pe");
            
            for (Int_t i = 0; i < nGenerators; i++){
                if (histoEtaInputMCWWeights[i]){
                    SetStyleHisto(histoEtaInputMCWWeights[i], 3., lineStyleMC[i], colorMC[i]);  
                    histoEtaInputMCWWeights[i]->Draw("same,hist,c");
                }
            }    
            graphYieldEta->Draw("p,same,e1");
            labelSpectraEtaLabel->Draw();
            
            legendSpectraEta->Draw();
    
        canvasSpectra->Update();
        canvasSpectra->Print(Form("%s/Eta_Spectra_MCWeighted_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));
    
    }
    
    // plot result with input
    canvasSpectra->cd();
    histo1DSpectra->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphYieldEta, markerStyleSpectrum, markerSizeSpectrum, colorData , colorData);
        graphYieldEta->Draw("p,same,e1");

        fitEtaYield->SetLineColor(colorData);
        fitEtaYield->Draw("same");
    
        labelSpectraEtaLabel->Draw();
    
    canvasSpectra->Print(Form("%s/Eta_Spectra_WithFit_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************* Ratio of data to fit and MC input to fit for eta *************************************
    // **********************************************************************************************************************
    // Calculation of ratio histograms
    TGraphAsymmErrors* graphRatioEtaDatatoFit   = CalculateGraphErrRatioToFit (graphYieldEta, fitEtaYield);
    TH1D* histoEtaRatioMCInputsToFit[5]         = {NULL, NULL, NULL, NULL, NULL};
    for (Int_t i = 0; i < nGenerators; i++){
        histoEtaRatioMCInputsToFit[i]           = CalculateHistoRatioToFit (histoEtaInputMCWOWeights[i], fitEtaYield,kTRUE);
        SetStyleHisto(histoEtaRatioMCInputsToFit[i], 2, lineStyleMC[i], colorMC[i]);
    }

    // drawing
    canvasRatioToFit->cd();    
    histo1DRatio->DrawCopy(); 
    
        TLegend* legendRatioEta =  GetAndSetLegend2(0.11, 0.12, 0.4, 0.12+0.035*(1+nGenerators)*1.15, 0.035, 1, "", 42, 0.2); 
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaDatatoFit, markerStyleSpectrum, markerSizeSpectrum, colorData , colorData);
        graphRatioEtaDatatoFit->Draw("e,p");  
        legendRatioEta->AddEntry(graphRatioEtaDatatoFit,"Data/fit to Data","pe");
        
        for (Int_t i = 0; i < nGenerators; i++){
            histoEtaRatioMCInputsToFit[i]->Draw("same,hist,l");
            legendRatioEta->AddEntry(histoEtaRatioMCInputsToFit[i],Form("%s/ fit to Data", generatorName[i].Data()),"l");
        }
        graphRatioEtaDatatoFit->Draw("e,p");  
        legendRatioEta->Draw();

        labelEnergy2760GeVRatio->Draw();    
        TLatex *labelSpectraEtaLabelRatio   = new TLatex(0.12,0.88,"#eta #rightarrow #gamma #gamma");
        SetStyleTLatex( labelSpectraEtaLabelRatio, 0.035, 4, 1, 42, kTRUE, 11);
        labelSpectraEtaLabelRatio->Draw();
    
        DrawGammaLines(0.23, graphYieldPi0->GetX()[graphYieldPi0->GetN()-1]*2 ,1., 1.,0.1);
    canvasRatioToFit->Update();
    canvasRatioToFit->SaveAs(Form("%s/Eta_RatioToDataFit_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));
    
    
    //	**********************************************************************************************************************
    //	******************************Compare Efficiencies for pi0 and eta for different MC *******************************
    //	**********************************************************************************************************************
    TCanvas* canvasEffi     = new TCanvas("canvasEffi", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEffi,  0.09, 0.01, 0.015, 0.08);
    canvasEffi->SetLogy();
    
        TH2F * histo2DPi0Effi;
        histo2DPi0Effi                           = new TH2F("histo2DPi0Effi", "histo2DPi0Effi",1000, 0., 10.5, 1000, 1e-5, 1e-2 );
        SetStyleHistoTH2ForGraphs( histo2DPi0Effi, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#pi^{0}}", 
                                0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
        histo2DPi0Effi->GetYaxis()->SetLabelOffset(0.01);
        histo2DPi0Effi->DrawCopy(); 

        TLegend* legendEfficiencyPi0 =  GetAndSetLegend2(0.56, 0.12, 0.93, 0.12+0.035*(nGenerators)*1.15, 0.035, 1, "", 42, 0.2); 
        for (Int_t i = 0; i < nGenerators; i++){
            DrawGammaSetMarker(histoPi0Efficiency[i], markerStyleMC[i], markerSizeSpectrum, colorMC[i] , colorMC[i]);
            histoPi0Efficiency[i]->Draw("p,same,e1");
            legendEfficiencyPi0->AddEntry(histoPi0Efficiency[i],generatorName[i].Data(),"p");
        }    
        legendEfficiencyPi0->Draw();
    
    canvasEffi->Update();
    canvasEffi->Print(Form("%s/Pi0_Efficiency_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));

        TH2F * histo2DEtaEffi;
        histo2DEtaEffi                           = new TH2F("histo2DEtaEffi", "histo2DEtaEffi",1000, 0., 10.5, 1000, 1e-5, 1e-2 );
        SetStyleHistoTH2ForGraphs( histo2DEtaEffi, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#eta}", 
                                0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
        histo2DEtaEffi->GetYaxis()->SetLabelOffset(0.01);
        histo2DEtaEffi->DrawCopy(); 

        TLegend* legendEfficiencyEta =  GetAndSetLegend2(0.56, 0.12, 0.93, 0.12+0.035*(nGenerators)*1.15, 0.035, 1, "", 42, 0.2); 
        for (Int_t i = 0; i < nGenerators; i++){
            DrawGammaSetMarker(histoEtaEfficiency[i], markerStyleMC[i], markerSizeSpectrum, colorMC[i] , colorMC[i]);
            histoEtaEfficiency[i]->Draw("p,same,e1");
            legendEfficiencyEta->AddEntry(histoEtaEfficiency[i],generatorName[i].Data(),"p");
        }    
        legendEfficiencyEta->Draw();
    
    canvasEffi->Update();
    canvasEffi->Print(Form("%s/Eta_Efficiency_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));

    canvasEffi->SetLogy(0);
    canvasEffi->SetTopMargin(0.035);
    
        histo2DPi0Effi->GetYaxis()->SetRangeUser(1e-5,3E-3);
        histo2DPi0Effi->DrawCopy(); 

        for (Int_t i = 0; i < nGenerators; i++){
            histoPi0Efficiency[i]->Draw("p,same,e1");
        }    
        legendEfficiencyPi0->Draw();
    
    canvasEffi->Update();
    canvasEffi->Print(Form("%s/Pi0_Efficiency_%s_LinY.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));

        histo2DEtaEffi->GetYaxis()->SetRangeUser(1e-5,3E-3);
        histo2DEtaEffi->DrawCopy(); 

        for (Int_t i = 0; i < nGenerators; i++){
            histoEtaEfficiency[i]->Draw("p,same,e1");
        }    
        legendEfficiencyEta->Draw();
    
    canvasEffi->Update();
    canvasEffi->Print(Form("%s/Eta_Efficiency_%s_LinY.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));
    
    //	**********************************************************************************************************************
    //	****************************************Write fits & input MC spectra to file ****************************************
    //	**********************************************************************************************************************	
    TFile fMCSpectraInput("MCSpectraInputpp.root","UPDATE");
        if (fitPi0Yield){
            fitPi0Yield->SetRange(0,100);
            fitPi0Yield->Write(Form("Pi0_Fit_Data_%s",collisionSystemForWriting.Data()),TObject::kOverwrite);
        }
        if (fitEtaYield){
            fitEtaYield->SetRange(0,100);
            fitEtaYield->Write(Form("Eta_Fit_Data_%s",collisionSystemForWriting.Data()),TObject::kOverwrite);     
        }
        
        for (Int_t i = 0; i < nGenerators; i++){
            cout << "writing everything for " << generatorName[i].Data() << endl;
            histoPi0InputMCWOWeights[i]->SetTitle(Form("Pi0_%s",outputString[i].Data()));
            histoPi0InputMCWOWeights[i]->Write(Form("Pi0_%s",outputString[i].Data()),TObject::kOverwrite);
            histoEtaInputMCWOWeights[i]->SetTitle(Form("Eta_%s",outputString[i].Data()));
            histoEtaInputMCWOWeights[i]->Write(Form("Eta_%s",outputString[i].Data()),TObject::kOverwrite);
        }
    fMCSpectraInput.Close();
            
}
