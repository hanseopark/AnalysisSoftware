//**************************************************************************************************
//**************************** CompareDifferentDirectoriesMerged ***********************************
//**************************************************************************************************

/***************************************************************************************************
 ******     provided by PCM Group, PWGGA,                                                       *****
 ******     Friederike Bock, friederike.bock@cern.ch                                            *****
 ****************************************************************************************************/

#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
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
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/CombinationFunctions.h"

struct SysErrorConversion {
   Double_t value;
   Double_t error;
   //   TString name;
};

void CompareDifferentDirectoriesMerged( TString FolderList          = "",
                                        TString suffix              = "gif",
                                        TString meson               = "",
                                        Bool_t kIsMC                = kFALSE,
                                        TString optionEnergy        = "",
                                        Int_t NumberOfCuts          = 1,
                                        TString optionPeriod        = "No",
                                        Int_t mode                  = 10,
                                        TString cutVariationName    = "NonLinearity",
                                        Bool_t setFullPathInInputFile   = kFALSE
                                      ){

    if (!(mode == 10 || mode == 11 )){
        cout << "incorrect mode: " << mode << endl;
        return ;
    }

    // Initialize arrays
    TString     fileDirectory[50];
    TString     cutNumber[50];
    TString     cutStringsName[50];

    // prepare nice plotting algorithms
    StyleSettingsThesis();
    SetPlotStyle();

    // Set output folder
    TString outputDir = Form("CutStudies/%s/%s",optionEnergy.Data(), cutVariationName.Data());
    TString outputDirRootFile = Form("CutStudies/%s",optionEnergy.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // Set meson name for plotting
    TString textMeson;
    if (meson.CompareTo("Pi0")==0 || meson.CompareTo("Pi0EtaBinning")==0){
        textMeson = "#pi^{0}";
    } else {
        textMeson = "#eta";
    }

    // Set MC/data string for output
    TString     prefix2;
    if (kIsMC){
        prefix2 =           "MC";
    } else {
        prefix2 =           "data";
    }

    // Determine collsision system string
    TString collisionSystem= ReturnFullCollisionsSystem(optionEnergy);
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }

    TString detProcess                                      = ReturnFullTextReconstructionProcess(mode, 0, textMeson);
    TString fNLMString                                      = "";
    Int_t fNLMmin                                           = 0;

    // Define colors for comparisons
    Color_t color[20] = {kBlack, kAzure, kGreen+2,kOrange+2,kRed, kViolet,  kBlue-9, kSpring+10,
                        kCyan+3, kCyan-10, kCyan, kGreen+4, kGreen-9,
                        kGreen,  kYellow+4, kYellow+3, kMagenta+4,
                        kMagenta-8, kGray, kGray+3};

    // read folder and name from file
    ifstream in(FolderList.Data());
    cout<<"Available Cuts:"<<endl;
    TString folderName;
    TString cutName;
    TString cutNr;
    Int_t Number = 0;
    while(!in.eof() ){
        in >> folderName >> cutNr >> cutName;
        cutName.ReplaceAll("_"," ");
        fileDirectory[Number] = folderName;
        cutStringsName[Number] = cutName;
        cutNumber[Number] = cutNr;
        cout<< fileDirectory[Number]<< "\t" << cutNumber[Number]<< "\t"<< cutStringsName[Number] <<endl;
        Number++;
    }
    cout<<"=========================="<<endl;

    // Definition of necessary histogram arrays
    const Int_t ConstNumberOfCuts = 20;
    TString FileNameCorrected[ConstNumberOfCuts];

    TFile *Cutcorrfile[ConstNumberOfCuts];

    TH1D *histoCorrectedYieldCut[ConstNumberOfCuts];
    TH1D *histoTrueEffiCut[ConstNumberOfCuts];
    TH1D *histoAcceptanceCut[ConstNumberOfCuts];
    TH1D *histoPurityCut[ConstNumberOfCuts];
    TH1D *histoRawYieldCut[ConstNumberOfCuts];

    TH1D *histoRatioCorrectedYieldCut[ConstNumberOfCuts];
    TH1D *histoRatioTrueEffiCut[ConstNumberOfCuts];
    TH1D *histoRatioAcceptanceCut[ConstNumberOfCuts];
    TH1D *histoRatioPurityCut[ConstNumberOfCuts];
    TH1D *histoRatioRawYieldCut[ConstNumberOfCuts];

    Double_t maxPt  = 0;
    for (Int_t i=0; i< NumberOfCuts; i++){

        // Decode individual cutnumber
        TString fEventCutSelection                  = "";
        TString fClusterCutSelection                = "";
        TString fClusterMergedCutSelection          = "";
        TString dummyString                         = "";
        TString fMesonCutSelection                  = "";
        ReturnSeparatedCutNumberAdvanced( cutNumber[i].Data(), fEventCutSelection, fClusterCutSelection, fClusterMergedCutSelection, dummyString, fMesonCutSelection, mode);

        // read file with corrections
        if(setFullPathInInputFile)
            FileNameCorrected[i] = Form("%s/%s_%s_GammaMergedCorrection_%s.root", fileDirectory[i].Data(), meson.Data(), prefix2.Data(), cutNumber[i].Data());
        else
            FileNameCorrected[i] = Form("%s%s/%s/%s_%s_GammaMergedCorrection_%s.root", fileDirectory[i].Data(), cutNumber[i].Data(), optionEnergy.Data(), meson.Data(), prefix2.Data(), cutNumber[i].Data());
        cout<< FileNameCorrected[i] << endl;
        Cutcorrfile[i] = new TFile(FileNameCorrected[i]);
        if (Cutcorrfile[i]->IsZombie()) return;

        if (i == 0){
            fNLMmin = ReturnClusterNLM(fClusterMergedCutSelection);
            if (fNLMmin)
                fNLMString                          = Form("%i local maximum", fNLMmin);
            else
                fNLMString                          = Form("%i local maxima", fNLMmin);
        } else {
            if ( fNLMmin != ReturnClusterNLM(fClusterMergedCutSelection)){
                fNLMString                          = "";
                fNLMmin                             = 0;
            }
        }


        // Set correct histogram name for corrected yield and efficiency
        TString nameCorrectedYield                          = "CorrectedYieldTrueEff";
        TString nameEfficiency                              = "PrimaryMesonEfficiency";
        TString namePurity                                  = "MesonPurity";
        TString nameAcceptance                              = "fHistoMCAcceptancePt";

        // Read histograms and rename them from the original files for each cut
        histoCorrectedYieldCut[i]   = (TH1D*)Cutcorrfile[i]->Get(nameCorrectedYield.Data());
        histoCorrectedYieldCut[i]->SetName(Form("%s_%s",nameCorrectedYield.Data(),cutStringsName[i].Data()));
        histoTrueEffiCut[i]         = (TH1D*)Cutcorrfile[i]->Get(nameEfficiency.Data());
        histoTrueEffiCut[i]->SetName(Form("%s_%s",nameEfficiency.Data(), cutStringsName[i].Data()));
        histoAcceptanceCut[i]       =(TH1D*)Cutcorrfile[i]->Get(nameAcceptance.Data());
        histoAcceptanceCut[i]->SetName(Form("AcceptPt_%s", cutStringsName[i].Data()));
        histoPurityCut[i]       =(TH1D*)Cutcorrfile[i]->Get(namePurity.Data());
        histoPurityCut[i]->SetName(Form("PurityPt_%s", cutStringsName[i].Data()));

        histoRawYieldCut[i]         = (TH1D*)Cutcorrfile[i]->Get("histoYieldMesonPerEvent");
        histoRawYieldCut[i]->SetName(Form("histoYieldMesonPerEvent_%s",cutStringsName[i].Data()));

        // Calculate ratios for comparisons
        histoRatioCorrectedYieldCut[i] = (TH1D*) histoCorrectedYieldCut[i]->Clone(Form("histoRatioCorrectedYieldCut_%s",cutStringsName[i].Data()));
        histoRatioCorrectedYieldCut[i]->Divide(histoRatioCorrectedYieldCut[i],histoCorrectedYieldCut[0],1.,1.,"B");
        if (i > 0){
            maxPt= histoCorrectedYieldCut[i]->GetBinCenter(histoCorrectedYieldCut[i]->GetNbinsX()) + 0.5* histoCorrectedYieldCut[i]->GetBinWidth(histoCorrectedYieldCut[i]->GetNbinsX());
        }
        histoRatioTrueEffiCut[i]    = (TH1D*) histoTrueEffiCut[i]->Clone(Form("histoRatioTrueEffiCut_%s",cutStringsName[i].Data()));
        histoRatioTrueEffiCut[i]->Divide(histoRatioTrueEffiCut[i],histoTrueEffiCut[0],1.,1.,"B");
        histoRatioAcceptanceCut[i]  = (TH1D*) histoAcceptanceCut[i]->Clone(Form("histoRatioAcceptanceCut_%s",cutStringsName[i].Data()));
        histoRatioAcceptanceCut[i]->Divide(histoRatioAcceptanceCut[i],histoAcceptanceCut[0],1.,1.,"B");
        histoRatioPurityCut[i]  = (TH1D*) histoPurityCut[i]->Clone(Form("histoRatioPurityCut_%s",cutStringsName[i].Data()));
        histoRatioPurityCut[i]->Divide(histoRatioPurityCut[i],histoPurityCut[0],1.,1.,"B");

        histoRatioRawYieldCut[i]    = (TH1D*) histoRawYieldCut[i]->Clone(Form("histoRatioRawYieldCut_%s",cutStringsName[i].Data()));
        histoRatioRawYieldCut[i]->Divide(histoRatioRawYieldCut[i],histoRawYieldCut[0],1.,1.,"B");
    }
    cout<<"=========================="<<endl;


    //**************************************************************************************
    //********************* Plotting RAW-Yield *********************************************
    //**************************************************************************************

        // Canvas Definition
        TCanvas* canvasRawYieldMeson = new TCanvas("canvasRawYieldMeson","",1350,1500);
        DrawGammaCanvasSettings( canvasRawYieldMeson,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padRawYield = new TPad("padRawYield", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padRawYield, 0.12, 0.02, 0.02, 0.);
        padRawYield->SetLogy();
        padRawYield->Draw();
        // lower pad definition
        TPad* padRawYieldRatios = new TPad("padRawYieldRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padRawYieldRatios, 0.12, 0.02, 0.0, 0.2);
        padRawYieldRatios->Draw();

        // Plot raw yield in uppper panel
        padRawYield->cd();
        TLegend* legendRawMeson = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
        if (cutVariationName.Contains("dEdxPi")){
            legendRawMeson->SetTextSize(0.02);
        }

        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
                Double_t scaleFactorRaw = 5.;
                DrawGammaSetMarker(histoRawYieldCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoRawYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s RAW Yield/(#it{N}_{ev} #it{N}_{coll})",textMeson.Data()),
                                        kTRUE, scaleFactorRaw, 10e-10, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                legendRawMeson->AddEntry(histoRawYieldCut[i],Form("standard: %s",cutStringsName[i].Data()));
            }
            else {
                if(i<20){
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRawYieldCut[i]->DrawCopy("same,e1,p");
                legendRawMeson->AddEntry(histoRawYieldCut[i],cutStringsName[i].Data());
            }
        }
        legendRawMeson->Draw();
        // Labeling of plot
        TLatex *labelCollisionSystem = new TLatex(0.55,0.91,collisionSystem.Data());
        SetStyleTLatex( labelCollisionSystem, 0.038,4);
        labelCollisionSystem->Draw();
        TLatex *labelDetProcess = NULL;
        if (optionPeriod.CompareTo("No") != 0){
            TLatex *labelPeriod = new TLatex(0.55,0.86,optionPeriod.Data());
            SetStyleTLatex( labelPeriod, 0.038,4);
            labelPeriod->Draw();
            if (detProcess.CompareTo("")!= 0){
                labelDetProcess = new TLatex(0.55,0.81,detProcess.Data());
                SetStyleTLatex( labelDetProcess, 0.038,4);
                labelDetProcess->Draw();
            }
        } else {
            if (detProcess.CompareTo("")!= 0){
                labelDetProcess = new TLatex(0.55,0.86,detProcess.Data());
                SetStyleTLatex( labelDetProcess, 0.038,4);
                labelDetProcess->Draw();
            }
        }
        // Plot ratio of raw yields in lower panel
        padRawYieldRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55;
                SetStyleHistoTH1ForGraphs(histoRatioRawYieldCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioRawYieldCut[i], 20, 1.,color[0],color[0]);
                histoRatioRawYieldCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoRatioRawYieldCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioRawYieldCut[i]->DrawCopy("same,e1,p");
            }
        }
        DrawGammaLines(0., maxPt,1., 1.,0.1);

        canvasRawYieldMeson->Update();
        canvasRawYieldMeson->SaveAs(Form("%s/%s_%s_RAWYield.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasRawYieldMeson;


    //*****************************************************************************************
    //******************* Compare Corrected Yields ********************************************
    //*****************************************************************************************
    // Define canvas
    TCanvas* canvasCorrectedYieldMeson = new TCanvas("canvasCorrectedYieldMeson","",1350,1500);
    DrawGammaCanvasSettings( canvasCorrectedYieldMeson,  0.13, 0.02, 0.02, 0.09);
    // Define upper panel
    TPad* padCorrectedYield = new TPad("padCorrectedYield", "", 0., 0.33, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYield, 0.12, 0.02, 0.02, 0.);
    padCorrectedYield->SetLogy(1);
    padCorrectedYield->Draw();
    // Define lower panel
    TPad* padCorrectedYieldRatios = new TPad("padCorrectedYieldRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYieldRatios, 0.12, 0.02, 0.0, 0.2);
    padCorrectedYieldRatios->SetLogy(0);
    padCorrectedYieldRatios->Draw();

    // Plot corrected yield in upper panel
    padCorrectedYield->cd();
        TLegend* legendCorrectedYieldMeson = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
                DrawAutoGammaMesonHistos( histoCorrectedYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)",
                                        kTRUE, 5., 5e-10,kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                DrawGammaSetMarker(histoCorrectedYieldCut[i], 20, 1., color[0], color[0]);
                histoCorrectedYieldCut[i]->DrawCopy("e1,p");
                legendCorrectedYieldMeson->AddEntry(histoCorrectedYieldCut[i], Form("standard: %s",cutStringsName[i].Data()));
            }
            else{
                if(i<20){
                    DrawGammaSetMarker(histoCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoCorrectedYieldCut[i]->DrawCopy("same,e1,p");
                legendCorrectedYieldMeson->AddEntry(histoCorrectedYieldCut[i], cutStringsName[i].Data());
            }

        }
        legendCorrectedYieldMeson->Draw();
        // Labeling of plot
        TLatex *labelCollisionSystem3 = new TLatex(0.55,0.91,collisionSystem.Data());
        SetStyleTLatex( labelCollisionSystem3, 0.038,4);
        labelCollisionSystem3->Draw();
        if (optionPeriod.CompareTo("No") != 0){
            TLatex *labelPeriod = new TLatex(0.55,0.86,optionPeriod.Data());
            SetStyleTLatex( labelPeriod, 0.038,4);
            labelPeriod->Draw();
        }
        if (labelDetProcess) labelDetProcess->Draw();
        // plot ratio of corrected yields in lower panel
        padCorrectedYieldRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.75;
                Double_t maxYRatio = 1.25;
                if(!cutVariationName.CompareTo("pThcompareMerged")){
                    minYRatio = 0.8;
                    maxYRatio = 1.05;
                }
                SetStyleHistoTH1ForGraphs(histoRatioCorrectedYieldCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20, 1.,color[0],color[0]);
                histoRatioCorrectedYieldCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoRatioCorrectedYieldCut[i]->DrawCopy("p,e1");
            }
            else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioCorrectedYieldCut[i]->DrawCopy("same,e1,p");
            }
        }
        DrawGammaLines(0., maxPt,1., 1.,0.1);

    canvasCorrectedYieldMeson->Update();
    canvasCorrectedYieldMeson->SaveAs(Form("%s/%s_%s_CorrectedYield.%s",outputDir.Data(), meson.Data(),prefix2.Data(),suffix.Data()));
    delete canvasCorrectedYieldMeson;


    //**************************************************************************************
    //********************* Plotting Efficiency *********************************************
    //**************************************************************************************
    // Define canvas
    TCanvas* canvasTrueEffiMeson = new TCanvas("canvasTrueEffiMeson","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasTrueEffiMeson,  0.13, 0.02, 0.02, 0.09);
    // Define upper panel
    TPad* padTrueEffi = new TPad("padTrueEffi", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padTrueEffi, 0.12, 0.02, 0.04, 0.);
    padTrueEffi->Draw();
    // Define lower panel
    TPad* padTrueEffiRatios = new TPad("padTrueEffiRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padTrueEffiRatios, 0.12, 0.02, 0.0, 0.2);
    padTrueEffiRatios->Draw();

      // draw efficiency in upper panel
      padTrueEffi->cd();
      TLegend* legendEffiMeson = GetAndSetLegend2(0.15,0.93-0.03*NumberOfCuts,0.3,0.93, 1500*0.75*0.032);
      for(Int_t i = 0; i< NumberOfCuts; i++){
          if(i == 0){
              DrawGammaSetMarker(histoTrueEffiCut[i], 20, 1., color[0], color[0]);
              DrawAutoGammaMesonHistos( histoTrueEffiCut[i],
                                      "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{eff,%s}",textMeson.Data()),
                                      kTRUE, 5., 10e-10,kFALSE,
                                      kTRUE, -0.1, 0.00030,
                                      kFALSE, 0., 10.);
              if (mode == 10 && fNLMmin == 1) histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(-0.01,0.2);
              if (!optionEnergy.CompareTo("pPb_8TeV")) histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(-0.01,2.65);
              if (mode == 10 && fNLMmin == 2) histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(-0.01,0.35);
              histoTrueEffiCut[i]->DrawCopy("e1,p");
              legendEffiMeson->AddEntry(histoTrueEffiCut[i],Form("standard: %s",cutStringsName[i].Data()));
          } else {
              if(i<20){
                  DrawGammaSetMarker(histoTrueEffiCut[i], 20+i, 1.,color[i],color[i]);
              } else {
                  DrawGammaSetMarker(histoTrueEffiCut[i], 20+i, 1.,color[i-20],color[i-20]);
              }
              histoTrueEffiCut[i]->DrawCopy("same,e1,p");
              legendEffiMeson->AddEntry(histoTrueEffiCut[i],cutStringsName[i].Data());
          }
      }
      legendEffiMeson->Draw();

      // Efficiency plot labeling
      TLatex *labelCollisionSystem4 = new TLatex(0.55,0.10,collisionSystem.Data());
      SetStyleTLatex( labelCollisionSystem4, 0.038,4);
      labelCollisionSystem4->Draw();
      TLatex* labelDetProcess2 = NULL;
      if (optionPeriod.CompareTo("No") != 0){
          TLatex *labelPeriod = new TLatex(0.55,0.91,optionPeriod.Data());
          SetStyleTLatex( labelPeriod, 0.038,4);
          labelPeriod->Draw();
          if (detProcess.CompareTo("")!= 0){
              labelDetProcess2 = new TLatex(0.55,0.86,detProcess.Data());
              SetStyleTLatex( labelDetProcess2, 0.038,4);
              labelDetProcess2->Draw();
          }
      } else {
          if (detProcess.CompareTo("")!= 0){
              labelDetProcess2 = new TLatex(0.55,0.05,detProcess.Data());
              SetStyleTLatex( labelDetProcess2, 0.038,4);
              labelDetProcess2->Draw();
          }
      }

      // Draw ratio of efficiencies in lower panel
      padTrueEffiRatios->cd();
      if( optionEnergy.Contains("Pb") ) padTrueEffiRatios->SetLogy(0);
      else padTrueEffiRatios->SetLogy(0);

      for(Int_t i = 0; i< NumberOfCuts; i++){
          if(i==0){
              Double_t minYRatio = 0.65;
              Double_t maxYRatio = 1.45;
              SetStyleHistoTH1ForGraphs(histoRatioTrueEffiCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
              DrawGammaSetMarker(histoRatioTrueEffiCut[i], 20, 1.,color[0],color[0]);
              histoRatioTrueEffiCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
              histoRatioTrueEffiCut[i]->DrawCopy("p,e1");
          } else{
              if(i<20){
                  DrawGammaSetMarker(histoRatioTrueEffiCut[i], 20+i, 1.,color[i],color[i]);
              } else {
                  DrawGammaSetMarker(histoRatioTrueEffiCut[i], 20+i, 1.,color[i-20],color[i-20]);
              }
              histoRatioTrueEffiCut[i]->DrawCopy("same,e1,p");
          }
          DrawGammaLines(0., maxPt,1., 1.,0.1);
          DrawGammaLines(0., maxPt, 1.1, 1.1, 0.1, kGray+1, 7);
          DrawGammaLines(0., maxPt, 0.9, 0.9, 0.1, kGray+1, 7);
      }

    canvasTrueEffiMeson->Update();
    canvasTrueEffiMeson->SaveAs(Form("%s/%s_%s_Efficiencies.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
    delete canvasTrueEffiMeson;

    //**************************************************************************************
    //********************* Plotting Acceptance *********************************************
    //**************************************************************************************
    // Define canvas
    TCanvas* canvasAcceptanceMeson = new TCanvas("canvasAcceptanceMeson","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceMeson,  0.13, 0.02, 0.02, 0.09);
    // Define upper panel
    TPad* padAcceptance = new TPad("padAcceptance", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padAcceptance, 0.12, 0.02, 0.04, 0.);
    padAcceptance->Draw();
    // Define lower panel
    TPad* padAcceptanceRatios = new TPad("padAcceptanceRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padAcceptanceRatios, 0.12, 0.02, 0.0, 0.2);
    padAcceptanceRatios->Draw();

      // draw efficiency in upper panel
      padAcceptance->cd();
      Double_t minAccept = 1;
      for (Int_t j = 1; j < histoAcceptanceCut[0]->GetNbinsX()+1; j++){
          if (histoAcceptanceCut[0]->GetBinContent(j) > 0){
            if (histoAcceptanceCut[0]->GetBinContent(j) < minAccept) minAccept = histoAcceptanceCut[0]->GetBinContent(j);
          }
      }

      TLegend* legendAcceptance = GetAndSetLegend2(0.55, 0.93-(1.1*NumberOfCuts*0.038), 0.95, 0.93,38);
      for(Int_t i = 0; i< NumberOfCuts; i++){
          if(i == 0){
              DrawGammaSetMarker(histoAcceptanceCut[i], 20, 1., color[0], color[0]);
              DrawAutoGammaMesonHistos( histoAcceptanceCut[i],
                                      "", "#it{p}_{T} (GeV/#it{c})", Form("A_{%s}",textMeson.Data()),
                                      kTRUE, 1.2, minAccept/1.5, kFALSE,
                                      kFALSE, 0., 1,
                                      kFALSE, 0., 10.);
              histoAcceptanceCut[i]->DrawCopy("e1,p");
              legendAcceptance->AddEntry(histoAcceptanceCut[i],Form("standard: %s",cutStringsName[i].Data()));
          } else {
              if(i<20){
                  DrawGammaSetMarker(histoAcceptanceCut[i], 20+i, 1.,color[i],color[i]);
              } else {
                  DrawGammaSetMarker(histoAcceptanceCut[i], 20+i, 1.,color[i-20],color[i-20]);
              }
              histoAcceptanceCut[i]->DrawCopy("same,e1,p");
              legendAcceptance->AddEntry(histoAcceptanceCut[i],cutStringsName[i].Data());
          }
      }
      legendAcceptance->Draw();

      // Efficiency plot labeling
      TLatex *labelCollisionSystem5 = new TLatex(0.55,0.10,collisionSystem.Data());
      SetStyleTLatex( labelCollisionSystem5, 0.038,4);
      labelCollisionSystem5->Draw();
      TLatex* labelDetProcess3 = NULL;
      if (optionPeriod.CompareTo("No") != 0){
          TLatex *labelPeriod = new TLatex(0.55,0.91,optionPeriod.Data());
          SetStyleTLatex( labelPeriod, 0.038,4);
          labelPeriod->Draw();
          if (detProcess.CompareTo("")!= 0){
              labelDetProcess3 = new TLatex(0.55,0.86,detProcess.Data());
              SetStyleTLatex( labelDetProcess3, 0.038,4);
              labelDetProcess3->Draw();
          }
      } else {
          if (detProcess.CompareTo("")!= 0){
              labelDetProcess3 = new TLatex(0.55,0.05,detProcess.Data());
              SetStyleTLatex( labelDetProcess3, 0.038,4);
              labelDetProcess3->Draw();
          }
      }

      // Draw ratio of efficiencies in lower panel
      padAcceptanceRatios->cd();
      if( optionEnergy.Contains("Pb") ) padAcceptanceRatios->SetLogy(0);
      else padAcceptanceRatios->SetLogy(0);

      for(Int_t i = 0; i< NumberOfCuts; i++){
          if(i==0){
              Double_t minYRatio = 0.45;
              Double_t maxYRatio = 1.55;
              SetStyleHistoTH1ForGraphs(histoRatioAcceptanceCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
              DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20, 1.,color[0],color[0]);
              histoRatioAcceptanceCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
              histoRatioAcceptanceCut[i]->DrawCopy("p,e1");
          } else{
              if(i<20){
                  DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20+i, 1.,color[i],color[i]);
              } else {
                  DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20+i, 1.,color[i-20],color[i-20]);
              }
              histoRatioAcceptanceCut[i]->DrawCopy("same,e1,p");
          }
          DrawGammaLines(0., maxPt,1., 1.,0.1);
      }

    canvasAcceptanceMeson->Update();
    canvasAcceptanceMeson->SaveAs(Form("%s/%s_%s_Acceptance.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
    delete canvasAcceptanceMeson;

    //**************************************************************************************
    //********************* Plotting Purity *********************************************
    //**************************************************************************************
    // Define canvas
    TCanvas* canvasPurityMeson = new TCanvas("canvasPurityMeson","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasPurityMeson,  0.13, 0.02, 0.02, 0.09);
    // Define upper panel
    TPad* padPurity = new TPad("padPurity", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padPurity, 0.12, 0.02, 0.04, 0.);
    padPurity->Draw();
    // Define lower panel
    TPad* padPurityRatios = new TPad("padPurityRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padPurityRatios, 0.12, 0.02, 0.0, 0.2);
    padPurityRatios->Draw();

      // draw efficiency in upper panel
      padPurity->cd();
      Double_t minPurity = 1;
      for (Int_t j = 1; j < histoPurityCut[0]->GetNbinsX()+1; j++){
          if (histoPurityCut[0]->GetBinContent(j) > 0){
            if (histoPurityCut[0]->GetBinContent(j) < minPurity) minPurity = histoPurityCut[0]->GetBinContent(j);
          }
      }

      TLegend* legendPurity = GetAndSetLegend2(0.55, 0.93-(1.1*NumberOfCuts*0.038), 0.95, 0.93,38);
      for(Int_t i = 0; i< NumberOfCuts; i++){
          if(i == 0){
              DrawGammaSetMarker(histoPurityCut[i], 20, 1., color[0], color[0]);
              DrawAutoGammaMesonHistos( histoPurityCut[i],
                                      "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{pur,%s}",textMeson.Data()),
                                      kTRUE, 1.2, minPurity/1.5, kFALSE,
                                      kFALSE, -0.1, 0.00030,
                                      kFALSE, 0., 10.);
              histoPurityCut[i]->DrawCopy("e1,p");
              legendPurity->AddEntry(histoPurityCut[i],Form("standard: %s",cutStringsName[i].Data()));
          } else {
              if(i<20){
                  DrawGammaSetMarker(histoPurityCut[i], 20+i, 1.,color[i],color[i]);
              } else {
                  DrawGammaSetMarker(histoPurityCut[i], 20+i, 1.,color[i-20],color[i-20]);
              }
              histoPurityCut[i]->DrawCopy("same,e1,p");
              legendPurity->AddEntry(histoPurityCut[i],cutStringsName[i].Data());
          }
      }
      legendPurity->Draw();

      // Efficiency plot labeling
      TLatex *labelCollisionSystem6 = new TLatex(0.55,0.10,collisionSystem.Data());
      SetStyleTLatex( labelCollisionSystem6, 0.038,4);
      labelCollisionSystem6->Draw();
      TLatex* labelDetProcess4 = NULL;
      if (optionPeriod.CompareTo("No") != 0){
          TLatex *labelPeriod = new TLatex(0.55,0.91,optionPeriod.Data());
          SetStyleTLatex( labelPeriod, 0.038,4);
          labelPeriod->Draw();
          if (detProcess.CompareTo("")!= 0){
              labelDetProcess4 = new TLatex(0.55,0.86,detProcess.Data());
              SetStyleTLatex( labelDetProcess4, 0.038,4);
              labelDetProcess4->Draw();
          }
      } else {
          if (detProcess.CompareTo("")!= 0){
              labelDetProcess4 = new TLatex(0.55,0.05,detProcess.Data());
              SetStyleTLatex( labelDetProcess4, 0.038,4);
              labelDetProcess4->Draw();
          }
      }

      // Draw ratio of efficiencies in lower panel
      padPurityRatios->cd();
      if( optionEnergy.Contains("Pb") ) padPurityRatios->SetLogy(0);
      else padPurityRatios->SetLogy(0);

      for(Int_t i = 0; i< NumberOfCuts; i++){
          if(i==0){
              Double_t minYRatio = 0.85;
              Double_t maxYRatio = 1.15;
              SetStyleHistoTH1ForGraphs(histoRatioPurityCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
              DrawGammaSetMarker(histoRatioPurityCut[i], 20, 1.,color[0],color[0]);
              histoRatioPurityCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
              histoRatioPurityCut[i]->DrawCopy("p,e1");
          } else{
              if(i<20){
                  DrawGammaSetMarker(histoRatioPurityCut[i], 20+i, 1.,color[i],color[i]);
              } else {
                  DrawGammaSetMarker(histoRatioPurityCut[i], 20+i, 1.,color[i-20],color[i-20]);
              }
              histoRatioPurityCut[i]->DrawCopy("same,e1,p");
          }
          DrawGammaLines(0., maxPt,1., 1.,0.1);
      }

    canvasPurityMeson->Update();
    canvasPurityMeson->SaveAs(Form("%s/%s_%s_Purity.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
    delete canvasPurityMeson;


  //*************************************************************************************************
  //******************** Output of the systematic Error due to Signal extraction for Pi0 ************
  //*************************************************************************************************
    // Determine number of bins
    Int_t NBinsPt = histoCorrectedYieldCut[0]->GetNbinsX();
    const Int_t NBinstPtConst = NBinsPt+1;

    // Create array of bin boundaries
    Double_t  BinsXCenter[NBinstPtConst];
    Double_t  BinsXWidth[NBinstPtConst];
    BinsXCenter[0] = 0;
    BinsXWidth[0]=0.;
    for (Int_t i = 1; i < NBinsPt +1; i++){
        BinsXCenter[i] = histoCorrectedYieldCut[0]->GetBinCenter(i);
        BinsXWidth[i]= histoCorrectedYieldCut[0]->GetBinWidth(i)/2.;
    }

    // Create array of Sys Err Objects and fill them
    SysErrorConversion SysErrCut[ConstNumberOfCuts][NBinstPtConst];
    SysErrorConversion SysErrCutRaw[ConstNumberOfCuts][NBinstPtConst];
    for (Int_t j = 0; j < NumberOfCuts; j++){
        for (Int_t i = 1; i < NBinsPt +1; i++){
            SysErrCut[j][i].value = histoCorrectedYieldCut[j]->GetBinContent(i);
            SysErrCut[j][i].error = histoCorrectedYieldCut[j]->GetBinError(i);
            SysErrCutRaw[j][i].value = histoRawYieldCut[j]->GetBinContent(i);
            SysErrCutRaw[j][i].error = histoRawYieldCut[j]->GetBinError(i);
        }
    }

    // Create Difference arrays
    Double_t DifferenceCut[ConstNumberOfCuts][NBinstPtConst];
    Double_t DifferenceErrorCut[ConstNumberOfCuts][NBinstPtConst];
    Double_t RelDifferenceCut[ConstNumberOfCuts][NBinstPtConst];
    Double_t RelDifferenceErrorCut[ConstNumberOfCuts][NBinstPtConst];
    Double_t RelDifferenceRawCut[ConstNumberOfCuts][NBinstPtConst];

    // Create largest difference array
    Double_t LargestDiffNeg[NBinstPtConst];
    Double_t LargestDiffPos[NBinstPtConst];
    Double_t LargestDiffErrorNeg[NBinstPtConst];
    Double_t LargestDiffErrorPos[NBinstPtConst];
    Double_t LargestDiffRelNeg[NBinstPtConst];
    Double_t LargestDiffRelPos[NBinstPtConst];
    Double_t LargestDiffRelErrorNeg[NBinstPtConst];
    Double_t LargestDiffRelErrorPos[NBinstPtConst];

    // Initialize all differences with 0
    for (Int_t j = 1; j < NumberOfCuts; j++){
        for ( Int_t i = 1; i < NBinstPtConst; i++) {
            DifferenceCut[j][i]=0.;
            DifferenceErrorCut[j][i]=0.;
            LargestDiffNeg[i]=0.;
            LargestDiffPos[i]=0.;
            LargestDiffErrorNeg[i]=0.;
            LargestDiffErrorPos[i]=0.;
            RelDifferenceCut[j][i]=0.;
            RelDifferenceRawCut[j][i]=0.;
            RelDifferenceErrorCut[j][i]=0.;
        }
    }

    // Calculate largest difference among cut variation
    for(Int_t j = 1; j < NumberOfCuts; j++){
        for (Int_t i = 1; i < NBinsPt +1; i++){
            // Calculate difference (rel/abs) and error for corrected yield
            DifferenceCut[j][i] = SysErrCut[j][i].value - SysErrCut[0][i].value;
            DifferenceErrorCut[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCut[j][i].error,2)-TMath::Power(SysErrCut[0][i].error,2)));
            if(SysErrCut[0][i].value != 0){
                RelDifferenceCut[j][i] = DifferenceCut[j][i]/SysErrCut[0][i].value*100. ;
                RelDifferenceErrorCut[j][i] = DifferenceErrorCut[j][i]/SysErrCut[0][i].value*100. ;
            } else {
                RelDifferenceCut[j][i] = -10000.;
                RelDifferenceErrorCut[j][i] = 100. ;
            }
            // Calculate relativ difference for raw yield
            if(SysErrCutRaw[0][i].value != 0){
                RelDifferenceRawCut[j][i] = (SysErrCutRaw[j][i].value - SysErrCutRaw[0][i].value)/SysErrCutRaw[0][i].value*100. ;
            } else {
                RelDifferenceRawCut[j][i] = -10000.;
            }
            // Calculate largest differences in positiv and negative direction
            if(DifferenceCut[j][i] < 0){ // largest negativ deviation
                // Take deviation if larger than previous largest deviation
                // and relative raw yield loss less than 75%
                if (TMath::Abs(LargestDiffNeg[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                    LargestDiffNeg[i] = DifferenceCut[j][i];
                    LargestDiffErrorNeg[i] = DifferenceErrorCut[j][i];
                }
            } else { // largest positive deviation
                // Take deviation if larger than previous largest deviation
                // and relative raw yield loss less than 75%
                if (TMath::Abs(LargestDiffPos[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                    LargestDiffPos[i] = DifferenceCut[j][i];
                    LargestDiffErrorPos[i] = DifferenceErrorCut[j][i];
                }
            }
        }
    }

    // Write systematic error input to log file
    TString SysErrDatname = Form("%s/%s_%s_SystematicErrorCutStudies.dat",outputDir.Data(),meson.Data(),prefix2.Data());
    fstream SysErrDat;
    SysErrDat.open(SysErrDatname.Data(), ios::out);
    SysErrDat << "Calculation of the systematic error due to the yield cuts" << endl;

    for (Int_t l=0; l< NumberOfCuts; l++){
        if (l == 0) {
            SysErrDat << endl <<"Bin" << "\t" << cutNumber[l] << "\t" <<endl;
            for(Int_t i = 1; i < (NBinsPt +1); i++){
                SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << endl;
            }
        } else{
            SysErrDat << endl <<"Bin" << "\t" << cutNumber[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
            for(Int_t i = 1; i < (NBinsPt +1); i++){
                if (RelDifferenceRawCut[l][i] > -75.){
                    SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << "\t" <<  DifferenceCut[l][i] << "\t"<< DifferenceErrorCut[l][i] << "\t"<< RelDifferenceCut[l][i] <<  "\t" << RelDifferenceErrorCut[l][i] <<"\t" << RelDifferenceRawCut[l][i]<< endl;
                } else {
                    SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << "\t" <<  DifferenceCut[l][i] << "\t"<< DifferenceErrorCut[l][i] << "\t"<< RelDifferenceCut[l][i] <<  "\t" << RelDifferenceErrorCut[l][i] <<"\t" << RelDifferenceRawCut[l][i]  <<"\t not considered in largest dev" <<endl;
                }
            }
        }
    }
    SysErrDat << endl;
    SysErrDat << endl;
    SysErrDat << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
    for(Int_t i = 1; i < (NBinsPt +1); i++){
        SysErrDat << BinsXCenter[i]  << "\t" << LargestDiffNeg[i] << "\t" <<LargestDiffErrorNeg[i]<< "\t" << LargestDiffPos[i] << "\t" << LargestDiffErrorPos[i]<<endl;
    }
    SysErrDat << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
    // Calculate largest relative deviations
    for(Int_t i = 0; i < (NBinsPt +1); i++){
        if ( SysErrCut[0][i].value != 0.){
            LargestDiffRelNeg[i] = - LargestDiffNeg[i]/SysErrCut[0][i].value*100.;
            LargestDiffRelPos[i] = LargestDiffPos[i]/SysErrCut[0][i].value*100.;
            LargestDiffRelErrorNeg[i] = - LargestDiffErrorNeg[i]/SysErrCut[0][i].value*100.;
            LargestDiffRelErrorPos[i] = LargestDiffErrorPos[i]/SysErrCut[0][i].value*100.;
            if (i > 0){
              SysErrDat << BinsXCenter[i] << "\t" << LargestDiffNeg[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffErrorNeg[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffPos[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffErrorPos[i]/SysErrCut[0][i].value*100.<<endl;
            } else {
              LargestDiffRelNeg[i] = 0.;
              LargestDiffRelPos[i] = 0.;
              LargestDiffRelErrorNeg[i] = 0.;
              LargestDiffRelErrorPos[i] = 0.;
            }
        } else {
            LargestDiffRelNeg[i] = 0.;
            LargestDiffRelPos[i] = 0.;
            LargestDiffRelErrorNeg[i] = 0.;
            LargestDiffRelErrorPos[i] = 0.;
        }
    }
    SysErrDat.close();

    // Create sys-err graphs
    TGraphAsymmErrors* SystErrGraphNeg = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffRelNeg, BinsXWidth, BinsXWidth, LargestDiffRelErrorNeg, LargestDiffRelErrorNeg);
    SystErrGraphNeg->SetName(Form("%s_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()));
    TGraphAsymmErrors* SystErrGraphPos = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffRelPos, BinsXWidth, BinsXWidth, LargestDiffRelErrorPos, LargestDiffRelErrorPos);
    SystErrGraphPos->SetName(Form("%s_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()));

    // Write sys-err graph to root output file
    TString Outputname = Form("%s/%s_%s_SystematicErrorCuts.root",outputDirRootFile.Data(),meson.Data(),prefix2.Data());
    TFile* SystematicErrorFile = new TFile(Outputname.Data(),"UPDATE");
        SystErrGraphPos->Write(Form("%s_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
        SystErrGraphNeg->Write(Form("%s_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
    SystematicErrorFile->Write();
    SystematicErrorFile->Close();

}
