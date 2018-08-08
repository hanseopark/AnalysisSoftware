/****************************************************************************************************************************
******      Friederike Bock, friederike.bock@cern.ch                                      *****
******      Hikari Murakami, Pedro Gonzalez
******      Efficiency study 
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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/CombinationFunctions.h"

extern TSystem* gSystem;


// run afterburner with data and from MC only added particles (specify path in baseDir)
// read CorrectedYieldTrueEff from Data afterburner output
// read MCYield_Meson_oldBin from MC afterurner output


//**********************************************************************************************************************
//***************************** Main function **************************************************************************
//**********************************************************************************************************************
void ExtractInputForWeightsPbPb(    TString suffix                  = "pdf", 
				    TString fOptEnergy              = "PbPb_5.02TeV",
				    TString periodData              = "LHC15o",
				    TString periodMC                = "LHC16h4",  
				    TString folderName              = "",
				    Bool_t runDrawReweighted        = kTRUE
                            ){

    //**********************************************************************************************************************
    //*************************************************** general style settings *******************************************
    //**********************************************************************************************************************
    gROOT->Reset();   
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();  
    SetPlotStyle();

    
    const Int_t nCentClasses    = 9;
    const Int_t nCentClassesUsed = 9;
    TString baseDir             = "/home/meike/analysis/results/photonconvResults/PbPb/pTWeights/";
    TString cent[nCentClasses]  = {"0-5%", "0-10%", "5-10%", "10-20%", "20-30%", "30-40%", "40-60%", "60-80%", "80-90%"};

    TString cutStringsAdded[nCentClasses] = {"30110a23_00200009247602008250404000_0652501500000000",        
					     "10110a23_00200009247602008250404000_0652501500000000",      
					     "31210a23_00200009247602008250404000_0652501500000000",      
					     "11210a23_00200009247602008250404000_0652501500000000",      
					     "12310a23_00200009247602008250404000_0652501500000000",      
					     "13410a23_00200009247602008250404000_0652501500000000",    
					     "14610a23_00200009247602008250404000_0652501500000000",     
 					     "16810a23_00200009247602008250404000_0652501500000000",     
					     "18910a23_00200009247602008250404000_0652501500000000"};     
    
    TString fitFunctionsPi0[nCentClasses] = {"doHag", "oHag", "oHag", "oHag", "oHag", "oHag", "modkfunc", "modkfunc", "doHag"};
    TString fitFunctionsEta[nCentClasses] = {"modkfunc", "modkfunc", "modkfunc", "oHag","oHag", "doHag", "oHag", "doHag", "modkfunc"};
    
    TString fileNameDataPi0[nCentClasses];
    TString fileNameDataEta[nCentClasses];
    TString fileNameAddedPi0[nCentClasses];
    TString fileNameAddedEta[nCentClasses];

    Color_t colorData[nCentClasses];
    Color_t colorMC[nCentClasses];

    Double_t minPtPlot    = 0.3;
    Double_t maxPtPlot    = 30.0;
    Double_t minPtFitPi0  = 0.4;
    Double_t maxPtFitPi0  = 20.0;
    Double_t minPtFitEta  = 1.0;
    Double_t maxPtFitEta  = 10.0;

    
    Style_t markerStyleSpectrum = 20;
    Size_t  markerSizeSpectrum  = 1.5;

    for(Int_t i=0; i<nCentClassesUsed; i++){
      fileNameDataPi0[i]    = Form("%s/%s/%s/%s/Pi0_data_GammaConvV1Correction_%s.root",baseDir.Data(),cutStringsAdded[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStringsAdded[i].Data());
      fileNameDataEta[i]    = Form("%s/%s/%s/%s/Eta_data_GammaConvV1Correction_%s.root",baseDir.Data(),cutStringsAdded[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStringsAdded[i].Data());
      fileNameAddedPi0[i]    = Form("%s/%s/%s/%s/Pi0_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStringsAdded[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStringsAdded[i].Data());
      fileNameAddedEta[i]    = Form("%s/%s/%s/%s/Eta_MC_GammaConvV1Correction_%s.root",baseDir.Data(),cutStringsAdded[i].Data(),fOptEnergy.Data(),folderName.Data(),cutStringsAdded[i].Data());

      colorData[i]           = GetColorDefaultColor( fOptEnergy, "", cent[i], kFALSE);
      colorMC[i]             = GetColorDefaultColor( fOptEnergy, "", cent[i], kFALSE);        
    }

     
    TString collisionSystem                         = ReturnFullCollisionsSystem(fOptEnergy);
    TString collisionSystemForWriting               = ReturnCollisionEnergyStringForTheory(fOptEnergy);
    Style_t lineStyleMC                             = 1;
    Style_t lineStyleData                           = 1;
    Style_t markerStyleMC                           = 20;
        
    TString dateForOutput           = ReturnDateStringForOutput();
    TString outputDir               = Form("../%s/%s/ExtractInputForWeights",suffix.Data(),dateForOutput.Data());
    cout << "output dir: " << outputDir.Data() << endl;
    gSystem->Exec("mkdir -p "+outputDir);

                
    //**********************************************************************************************************************
    //************************************** read in data and MC files Pi0 and Eta *****************************************
    //**********************************************************************************************************************

    
    TFile* filePi0DataInput[nCentClasses];
    TFile* fileEtaDataInput[nCentClasses];
    TFile* filePi0AddedInput[nCentClasses];
    TFile* fileEtaAddedInput[nCentClasses];

    TH1D* histoYieldDataPi0[nCentClasses];                // Data histo
    TH1D* histoYieldDataEta[nCentClasses];
    TGraphAsymmErrors* graphYieldDataPi0[nCentClasses];   // Data graph
    TGraphAsymmErrors* graphYieldDataEta[nCentClasses];
    TF1* fitPi0DataYield[nCentClasses];                       // Data fit
    TF1* fitEtaDataYield[nCentClasses]; 
    
    TH1D* histoPi0InputMCWOWeights[nCentClasses];     // MC histo without weights
    TH1D* histoPi0InputMCWWeights[nCentClasses];      // MC histo with    weights
    TH1D* histoPi0Efficiency[nCentClasses];
    TH1D* histoEtaInputMCWOWeights[nCentClasses];
    TH1D* histoEtaInputMCWWeights[nCentClasses];
    TH1D* histoEtaEfficiency[nCentClasses];

    TH1D* histoPi0RatioDataToFit[nCentClasses];   // Data histo to Data fit
    TH1D* histoEtaRatioDataToFit[nCentClasses];
	
    TH1D* histoPi0Ratio[nCentClasses];            // ratio MC(Added) without weight histo to Data fit
    TH1D* histoEtaRatio[nCentClasses];            

    
    for (Int_t i = 0; i < nCentClassesUsed; i++){

      filePi0DataInput[i]            = new TFile(fileNameDataPi0[i]);
      fileEtaDataInput[i]            = new TFile(fileNameDataEta[i]);
      
      filePi0AddedInput[i]         = new TFile(fileNameAddedPi0[i]);
      fileEtaAddedInput[i]         = new TFile(fileNameAddedEta[i]);

      histoYieldDataPi0[i]            = (TH1D*)filePi0DataInput[i]->Get("CorrectedYieldTrueEff");
      histoYieldDataEta[i]            = (TH1D*)fileEtaDataInput[i]->Get("CorrectedYieldTrueEff");
      
      histoPi0InputMCWOWeights[i]    = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson_oldBinWOWeights");  
      histoPi0InputMCWWeights[i]     = (TH1D*)filePi0AddedInput[i]->Get("MCYield_Meson_oldBin");  
      histoPi0Efficiency[i]          = (TH1D*)filePi0AddedInput[i]->Get("TrueMesonEffiPt");
      
      histoEtaInputMCWOWeights[i]    = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson_oldBinWOWeights");
      histoEtaInputMCWWeights[i]     = (TH1D*)fileEtaAddedInput[i]->Get("MCYield_Meson_oldBin");
      histoEtaEfficiency[i]          = (TH1D*)fileEtaAddedInput[i]->Get("TrueMesonEffiPt");

      // convert histo to graph
      graphYieldDataPi0[i]           = new TGraphAsymmErrors(histoYieldDataPi0[i]);
      graphYieldDataEta[i]           = new TGraphAsymmErrors(histoYieldDataEta[i]);
	
      // fit data
      fitPi0DataYield[i] = FitObject(fitFunctionsPi0[i].Data(),"fitPi0DataYield","Pi0", NULL, minPtFitPi0, maxPtFitPi0);
      graphYieldDataPi0[i]->Fit(fitPi0DataYield[i],"QNRME+", "QNRME+", minPtFitPi0, maxPtFitPi0);
      //fitPi0DataYield[i]->SetRange(minPtFitPi0, maxPtFitPi0);
      fitPi0DataYield[i]->SetRange(minPtPlot, maxPtPlot);
      fitPi0DataYield[i]->SetLineColor(colorMC[i]);
      fitPi0DataYield[i]->SetLineStyle(2);

      fitEtaDataYield[i]     = FitObject(fitFunctionsEta[i].Data(),"fitEtaDataYield","Eta", NULL, minPtFitEta, maxPtFitEta);
      graphYieldDataEta[i]->Fit(fitEtaDataYield[i],"QNRME+", "QNRME+", minPtFitEta, maxPtFitEta);
      //fitEtaDataYield[i]->SetRange(minPtFitEta, maxPtFitEta);
      fitEtaDataYield[i]->SetRange(minPtPlot, maxPtPlot);
      fitEtaDataYield[i]->SetLineColor(colorMC[i]);
      fitEtaDataYield[i]->SetLineStyle(2);

      // calculate ratio MC histo to Data fit
      histoPi0Ratio[i]  = CalculateHistoRatioToFit(histoPi0InputMCWOWeights[i], fitPi0DataYield[i],kTRUE);
      histoEtaRatio[i]  = CalculateHistoRatioToFit(histoEtaInputMCWOWeights[i], fitEtaDataYield[i],kTRUE);

      // calculate ratio Data histo to fit
      histoPi0RatioDataToFit[i] = CalculateHistoRatioToFit(histoYieldDataPi0[i], fitPi0DataYield[i], kTRUE);
      histoEtaRatioDataToFit[i] = CalculateHistoRatioToFit(histoYieldDataEta[i], fitEtaDataYield[i], kTRUE);
    
    }

    delete graphYieldDataPi0;
    delete graphYieldDataEta;
    delete filePi0DataInput;
    delete fileEtaDataInput;
    delete filePi0AddedInput;
    delete fileEtaAddedInput;

    
    //**********************************************************************************************************************
    //********************************************* PLOTTING SPECTRA *******************************************************
    //**********************************************************************************************************************
    TCanvas* canvasSpectra  = new TCanvas("canvasSpectra", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasSpectra,  0.13, 0.01, 0.015, 0.08);
    canvasSpectra->SetLogy();
    canvasSpectra->SetLogx();
    TLegend* legendSpectra            = GetAndSetLegend2(0.16, 0.12, 0.73, 0.12+0.035*(nCentClassesUsed), 0.035, 1, "", 42, 0.15);
    TLatex *labelCollisionSystem        = new TLatex(0.93,0.92,collisionSystem);
    SetStyleTLatex( labelCollisionSystem, 0.035, 4, 1, 42, kTRUE, 31);
    TLatex *labelSpectraPi0Label        = new TLatex(0.925,0.87,"#pi^{0} #rightarrow #gamma #gamma");
    SetStyleTLatex( labelSpectraPi0Label, 0.035, 4, 1, 42, kTRUE, 31);
    TLatex *labelSpectraEtaLabel        = new TLatex(0.925,0.87,"#eta #rightarrow #gamma #gamma");
    SetStyleTLatex( labelSpectraEtaLabel, 0.035, 4, 1, 42, kTRUE, 31);

    TH1F * histo1DSpectra;
    histo1DSpectra          = new TH1F("histo1DSpectra", "histo1DSpectra",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( histo1DSpectra, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
                            0.03, 0.035, 0.03, 0.035, 0.83, 1.7);
    histo1DSpectra->GetYaxis()->SetRangeUser(0.0000001, 5000.);
    histo1DSpectra->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DSpectra->GetXaxis()->SetLabelOffset(-0.01);
    histo1DSpectra->GetYaxis()->SetLabelOffset(0.01);
    histo1DSpectra->DrawCopy(); 

       
    // +++ Pi0 +++ Data  histo & MC Added histo UN-weighted +++
    for (Int_t i = 0; i < nCentClassesUsed; i++){
      histoYieldDataPi0[i]->Draw("p,same,e1");
      DrawGammaSetMarker(histoYieldDataPi0[i], markerStyleSpectrum, markerSizeSpectrum, colorData[i] , colorData[i]);
      if(i==0) legendSpectra->AddEntry(histoYieldDataPi0[i], Form(" Data %s", cent[i].Data()),"pe");
      else     legendSpectra->AddEntry(histoYieldDataPi0[i], Form(" %s", cent[i].Data()),"lpe");
      SetStyleHisto(histoPi0InputMCWOWeights[i], 3., lineStyleMC[i], colorMC[i]);  
      histoPi0InputMCWOWeights[i]->Draw("same,hist,c");
      if(i==0) legendSpectra->AddEntry(histoPi0InputMCWOWeights[i], Form(" MC %s", cent[i].Data()),"l");
      fitPi0DataYield[i]->Draw("same");
      if(i==0) legendSpectra->AddEntry(fitPi0DataYield[i], Form(" fit %s", cent[i].Data()),"l");
    }    
    histoYieldDataPi0[0]->Draw("p,same,e1");
    labelCollisionSystem->Draw();
    labelSpectraPi0Label->Draw();
    legendSpectra->Draw();
    canvasSpectra->Update();
    canvasSpectra->Print(Form("%s/Pi0_Spectra_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));

    // +++ Eta +++ Data  & MC Added UN-weighted +++
    legendSpectra->Clear();
    canvasSpectra->Clear();
    histo1DSpectra->DrawCopy();
    for (Int_t i = 0; i < nCentClassesUsed; i++){
      histoYieldDataEta[i]->Draw("p,same,e1");
      DrawGammaSetMarker(histoYieldDataEta[i], markerStyleSpectrum, markerSizeSpectrum, colorData[i] , colorData[i]);
      if(i==0) legendSpectra->AddEntry(histoYieldDataEta[i], Form(" Data %s", cent[i].Data()),"pe");
      else     legendSpectra->AddEntry(histoYieldDataEta[i], Form(" %s", cent[i].Data()),"lpe");
      SetStyleHisto(histoEtaInputMCWOWeights[i], 3., lineStyleMC[i], colorMC[i]);  
      histoEtaInputMCWOWeights[i]->Draw("same,hist,c");
      if(i==0) legendSpectra->AddEntry(histoEtaInputMCWOWeights[i], Form(" MC %s", cent[i].Data()),"l");
      fitEtaDataYield[i]->Draw("same");
      if(i==0) legendSpectra->AddEntry(fitEtaDataYield[i], Form(" fit %s", cent[i].Data()),"l");
    }    
    histoYieldDataEta[0]->Draw("p,same,e1");
    labelCollisionSystem->Draw();
    labelSpectraEtaLabel->Draw();
    legendSpectra->Draw();
    canvasSpectra->Update();
    canvasSpectra->Print(Form("%s/Eta_Spectra_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));

    
    //  +++ Pi0 +++ Data  & MC Added weighted +++ 
    //legendSpectra->Clear();
    //canvasSpectra->Clear();
    //histo1DSpectra->DrawCopy();

    //for (Int_t i = 0; i < nCentClassesUsed; i++){
    //  histoYieldDataPi0[i]->Draw("p,same,e1");
    //  SetStyleHisto(histoPi0InputMCWWeights[i], 3., lineStyleMC[i], colorMC[i]);  
    //  histoPi0InputMCWWeights[i]->Draw("same,hist,c");
   // }
    //histoYieldDataPi0[0]->Draw("p,same,e1");
    //labelCollisionSystem->Draw();
    //labelSpectraPi0Label->Draw();
    //legendSpectra->Draw();
    //canvasSpectra->Update();
    //canvasSpectra->Print(Form("%s/Pi0_Spectra_MCWeighted_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));
    


    // cleanup
    delete histo1DSpectra;
    delete labelCollisionSystem;
    delete labelSpectraPi0Label;
    delete labelSpectraEtaLabel;
    delete canvasSpectra;

    
    // **********************************************************************************************************************
    // ***************************************** PLOTTING RATIOS  ***********************************************************
    // **********************************************************************************************************************
    // definition of canvas for ratios
    TCanvas* canvasRatioToFit = new TCanvas("canvasRatioToFit","",1550,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToFit,  0.08, 0.015, 0.015, 0.07);
    canvasRatioToFit->SetLogx(1);
    canvasRatioToFit->SetLogy(0);
    TLegend* legendRatio =  GetAndSetLegend2(0.12, 0.52, 0.4, 0.52+0.035*(nCentClassesUsed), 0.035, 1, "", 42, 0.2);

    TLatex *labelCollisionSystemRatio     = new TLatex(0.12,0.92,collisionSystem.Data());
    SetStyleTLatex( labelCollisionSystemRatio, 0.035 , 4, 1, 42, kTRUE, 11);  
    TLatex *labelSpectraPi0LabelRatio   = new TLatex(0.12,0.88,"#pi^{0} #rightarrow #gamma #gamma");
    SetStyleTLatex( labelSpectraPi0LabelRatio, 0.035, 4, 1, 42, kTRUE, 11);
    TLatex *labelSpectraEtaLabelRatio   = new TLatex(0.12,0.88,"#eta #rightarrow #gamma #gamma");
    SetStyleTLatex( labelSpectraEtaLabelRatio, 0.035, 4, 1, 42, kTRUE, 11);
    
    // default histo
    TH1F * histo1DRatio;
    histo1DRatio          = new TH1F("histo1DRatio", "histo1DRatio",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( histo1DRatio, "#it{p}_{T} (GeV/#it{c})", "Ratio Data histo to fit",
                            0.03, 0.035, 0.03, 0.035, 0.83, 1.2);
    histo1DRatio->GetYaxis()->SetRangeUser(0.01, 4.9.);
    histo1DRatio->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DRatio->GetXaxis()->SetLabelOffset(-0.01);
    histo1DRatio->GetYaxis()->SetLabelOffset(0.01);
    histo1DRatio->DrawCopy(); 


    // +++ Pi0 +++ Data histo / Data Fit +++
    for (Int_t i = 0; i < nCentClassesUsed; i++){  
      SetStyleHisto(histoPi0RatioDataToFit[i], 2, lineStyleMC[i], colorMC[i]);
      histoPi0RatioDataToFit[i]->SetMarkerSize(0);    
      histoPi0RatioDataToFit[i]->SetMarkerColor(colorMC[i]);
      histoPi0RatioDataToFit[i]->Draw("same");
      legendRatio->AddEntry(histoPi0RatioDataToFit[i], Form("%s",cent[i].Data()),"l");
    }
    histo1DRatio->DrawCopy("SAME"); 
    legendRatio->Draw();
    labelCollisionSystemRatio->Draw();  
    labelSpectraPi0LabelRatio->Draw();
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
    canvasRatioToFit->Update();
    canvasRatioToFit->SaveAs(Form("%s/Pi0_RatioDataHistoToFit_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));
    
    //cleanup
    legendRatio->Clear();
    canvasRatioToFit->Clear();
    //histo1DRatio->GetYaxis()->SetRangeUser(0.01, 2.2);  // for Eta
    histo1DRatio->DrawCopy(); 

    // +++ Eta +++ Data histo / Data Fit +++
    for (Int_t i = 0; i < nCentClassesUsed; i++){  
      SetStyleHisto(histoEtaRatioDataToFit[i], 2, lineStyleMC[i], colorMC[i]);
      histoEtaRatioDataToFit[i]->SetMarkerSize(0);    
      histoEtaRatioDataToFit[i]->SetMarkerColor(colorMC[i]);
      histoEtaRatioDataToFit[i]->Draw("same");
      legendRatio->AddEntry(histoEtaRatioDataToFit[i], Form("%s",cent[i].Data()),"l");
    }
    histo1DRatio->DrawCopy("SAME"); 
    legendRatio->Draw();
    labelCollisionSystemRatio->Draw();  
    labelSpectraEtaLabelRatio->Draw();
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
    canvasRatioToFit->Update();
    canvasRatioToFit->SaveAs(Form("%s/Eta_RatioDataHistoToFit_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));
    
    //cleanup
    legendRatio->Clear();
    canvasRatioToFit->Clear();
    
    canvasRatioToFit->SetLogy(1);
    
    // default histo2
    TH1F * histo1DRatio2;
    histo1DRatio2          = new TH1F("histo1DRatio2", "histo1DRatio2",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( histo1DRatio2, "#it{p}_{T} (GeV/#it{c})", "Ratio MC to Data fit",
                            0.03, 0.035, 0.03, 0.035, 0.83, 1.2);
    //histo1DRatio2->GetYaxis()->SetRangeUser(0.0001, 500000.);  // for Pi0
    histo1DRatio2->GetYaxis()->SetRangeUser(0.00000005, 900000000.);
    histo1DRatio2->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DRatio2->GetXaxis()->SetLabelOffset(-0.01);
    histo1DRatio2->GetYaxis()->SetLabelOffset(0.01);
    histo1DRatio2->DrawCopy(); 
    
    // +++ Pi0 +++ MC Added / Data fit +++
    for (Int_t i = 0; i < nCentClassesUsed; i++){  
      SetStyleHisto(histoPi0Ratio[i], 2, lineStyleMC[i], colorMC[i]);
      histoPi0Ratio[i]->SetMarkerColor(colorMC[i]);
      histoPi0Ratio[i]->Draw("same");
      legendRatio->AddEntry(histoPi0Ratio[i], Form("%s",cent[i].Data()),"l");
    }
    histo1DRatio2->DrawCopy("SAME"); 
    legendRatio->Draw();
    labelCollisionSystemRatio->Draw();  
    labelSpectraPi0LabelRatio->Draw();
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
    canvasRatioToFit->Update();
    canvasRatioToFit->SaveAs(Form("%s/Pi0_RatioMCToData_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));

    //cleanup
    legendRatio->Clear();
    canvasRatioToFit->Clear();
    histo1DRatio2->DrawCopy();

    // +++ Eta +++ MC Added / Data fit +++
    for (Int_t i = 0; i < nCentClassesUsed; i++){  
      SetStyleHisto(histoEtaRatio[i], 2, lineStyleMC[i], colorMC[i]);
      histoEtaRatio[i]->SetMarkerColor(colorMC[i]);
      histoEtaRatio[i]->Draw("same");
      legendRatio->AddEntry(histoEtaRatio[i], Form("%s",cent[i].Data()),"l");
    }
    histo1DRatio2->DrawCopy("SAME"); 
    legendRatio->Draw();
    labelCollisionSystemRatio->Draw();  
    labelSpectraEtaLabelRatio->Draw();
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
    canvasRatioToFit->Update();
    canvasRatioToFit->SaveAs(Form("%s/Eta_RatioMCToData_%s.%s",outputDir.Data(),collisionSystemForWriting.Data(),suffix.Data()));
    
    delete histo1DRatio;
    delete histo1DRatio2;
    delete labelSpectraPi0LabelRatio;
    delete labelSpectraEtaLabelRatio;
    delete labelCollisionSystemRatio;
    delete canvasRatioToFit;
        
    //	**********************************************************************************************************************
    //	****************************************Write fits & input MC spectra to file ****************************************
    //	**********************************************************************************************************************	
    TFile fMCSpectraInput("MCSpectraInputPbPb.root","RECREATE");

    for (Int_t i = 0; i < nCentClassesUsed; i++){

      cout << "writing everything for " << cent[i].Data() << endl;
      TString cutNumber = cutStringsAdded[i];
      TString centCut   = cutNumber(0,3);  // first three digits of event cut
           
      histoPi0InputMCWOWeights[i]->SetTitle(Form("Pi0_LHC16h4_AP_%s_%s",collisionSystemForWriting.Data(),centCut.Data()));
      histoPi0InputMCWOWeights[i]->Write(Form("Pi0_LHC16h4_AP_%s_%s",collisionSystemForWriting.Data(),centCut.Data()),TObject::kOverwrite);
      
      histoEtaInputMCWOWeights[i]->SetTitle(Form("Eta_LHC16h4_AP_%s_%s",collisionSystemForWriting.Data(),centCut.Data()));
      histoEtaInputMCWOWeights[i]->Write(Form("Eta_LHC16h4_AP_%s_%s",collisionSystemForWriting.Data(),centCut.Data()),TObject::kOverwrite);

      fitPi0DataYield[i]->SetTitle(Form("Pi0_Data_%s_%s",collisionSystemForWriting.Data(),centCut.Data()));
      fitPi0DataYield[i]->Write(Form("Pi0_Data_%s_%s",collisionSystemForWriting.Data(),centCut.Data()),TObject::kOverwrite);

      fitEtaDataYield[i]->SetTitle(Form("Eta_Data_%s_%s",collisionSystemForWriting.Data(),centCut.Data()));
      fitEtaDataYield[i]->Write(Form("Eta_Data_%s_%s",collisionSystemForWriting.Data(),centCut.Data()),TObject::kOverwrite);

    }
      
    fMCSpectraInput.Close();
    
}
