//**************************************************************************************************
//**************************** CompareDifferentDirectories *****************************************
//**************************************************************************************************

/***************************************************************************************************
 ******                         provided by PCM Group, PWGGA,                                  *****
 ******                         Friederike Bock, friederike.bock@cern.ch                       *****
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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

struct SysErrorConversion {
   Double_t value;
   Double_t error;
   //   TString name;
};

void CompareDifferentDirectories(   TString FolderList              = "", 
                                    TString suffix                  = "gif", 
                                    TString meson                   = "", 
                                    Bool_t kIsMC                    = 0, 
                                    TString optionEnergy            = "", 
                                    Int_t NumberOfCuts              = 1, 
                                    TString optionPeriod            = "No",
                                    Int_t mode                      = 9,
                                    TString cutVariationName        = "NonLinearity"
                                ){

    // Initialize arrays
    TString fileDirectory[50];
    TString cutNumber[50];
    TString cutStringsName[50];
    
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
    TString 	prefix2;
    if (kIsMC){
        prefix2 = 			"MC";
    } else {
        prefix2 = 			"data";
    }
    
    // Determine collsision system string
    TString collisionSystem= ReturnFullCollisionsSystem(optionEnergy);   
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;     
    }

    TString detProcess = "";
    if (mode == 0 || mode == 9){
        detProcess = "#gamma with PCM";
    } else if (mode == 2 ){
        detProcess = "#gamma with PCM, EMCal";
    } else if (mode == 4 ){
        detProcess = "#gamma with EMCal";
    } else if (mode == 3 ){
        detProcess = "#gamma with PCM, PHOS";
    } else if (mode == 5 ){
        detProcess = "#gamma with PHOS";		
    }	
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
    const Int_t ConstNumberOfCuts = 10;
    TString FileNameCorrected[ConstNumberOfCuts];
    TString FileNameUnCorrected[ConstNumberOfCuts];

    TFile *Cutcorrfile[ConstNumberOfCuts];
    TFile *Cutuncorrfile[ConstNumberOfCuts];

    TH1D *histoCorrectedYieldCut[ConstNumberOfCuts];
    TH1D *histoTrueEffiCut[ConstNumberOfCuts];
    TH1D *histoAcceptanceCut[ConstNumberOfCuts];
    TH1D *histoRawYieldCut[ConstNumberOfCuts];
    TH1D *histoMassCut[ConstNumberOfCuts];
    TH1D *histoWidthCut[ConstNumberOfCuts];
    TH1D *histoSBCut[ConstNumberOfCuts];
    TH1D *histoClusterE[ConstNumberOfCuts];
    Bool_t isClusterE = kTRUE;

    TH1D *histoRatioCorrectedYieldCut[ConstNumberOfCuts];
    TH1D *histoRatioTrueEffiCut[ConstNumberOfCuts];
    TH1D *histoRatioAcceptanceCut[ConstNumberOfCuts];
    TH1D *histoRatioRawYieldCut[ConstNumberOfCuts];
    TH1D *histoRatioMassCut[ConstNumberOfCuts];
    TH1D *histoRatioWidthCut[ConstNumberOfCuts];
    TH1D *histoRatioSBCut[ConstNumberOfCuts];
    TH1D *histoRatioClusterE[ConstNumberOfCuts];
    
    Double_t maxPt	= 0;
    for (Int_t i=0; i< NumberOfCuts; i++){
        
        // Decode individual cutnumber
        TString fEventCutSelection="";
        TString fGammaCutSelection="";
        TString fClusterCutSelection="";
        TString fElectronCutSelection="";
        TString fMesonCutSelection="";
        cout << cutNumber[i].Data() << endl;
        ReturnSeparatedCutNumberAdvanced(cutNumber[i].Data(),fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
        
        // read file with corrections
        FileNameCorrected[i] = Form("%s%s/%s/%s_%s_GammaConvV1Correction_%s.root", fileDirectory[i].Data(), cutNumber[i].Data(), optionEnergy.Data(), meson.Data(), prefix2.Data(), cutNumber[i].Data());
        cout<< FileNameCorrected[i] << endl;
        Cutcorrfile[i] = new TFile(FileNameCorrected[i]);	
        if (Cutcorrfile[i]->IsZombie()) return;
        // read file without corrections
        FileNameUnCorrected[i] = Form("%s%s/%s/%s_%s_GammaConvV1WithoutCorrection_%s.root", fileDirectory[i].Data(), cutNumber[i].Data(), optionEnergy.Data(), meson.Data(), prefix2.Data(), cutNumber[i].Data());
        cout<< FileNameUnCorrected[i] << endl;
        Cutuncorrfile[i] = new TFile(FileNameUnCorrected[i]);
        if (Cutuncorrfile[i]->IsZombie()) return;
        
        // Set correct histogram name for corrected yield and efficiency
        TString nameCorrectedYield;
        TString nameEfficiency;
        nameCorrectedYield = "CorrectedYieldTrueEff";
        nameEfficiency = "TrueMesonEffiPt";
        if ( mode == 2 || mode == 3 || mode == 4 || mode == 5 ){
            nameCorrectedYield = "CorrectedYieldNormEff";
            nameEfficiency = "MesonEffiPt";
        }	
        
        // Read histograms and rename them from the original files for each cut
        histoCorrectedYieldCut[i]   = (TH1D*)Cutcorrfile[i]->Get(nameCorrectedYield.Data());
        histoCorrectedYieldCut[i]->SetName(Form("%s_%s",nameCorrectedYield.Data(),cutStringsName[i].Data()));
        histoTrueEffiCut[i]         = (TH1D*)Cutcorrfile[i]->Get(nameEfficiency.Data());
        histoTrueEffiCut[i]->SetName(Form("%s_%s",nameEfficiency.Data(), cutStringsName[i].Data()));
        histoAcceptanceCut[i]       =(TH1D*)Cutcorrfile[i]->Get("fMCMesonAccepPt");
        histoAcceptanceCut[i]->SetName(Form("AcceptPt_%s", cutStringsName[i].Data()));
        histoRawYieldCut[i]         = (TH1D*)Cutuncorrfile[i]->Get("histoYieldMesonPerEvent");
        histoRawYieldCut[i]->SetName(Form("histoYieldMesonPerEvent_%s",cutStringsName[i].Data()));
        histoMassCut[i]             = (TH1D*)Cutuncorrfile[i]->Get("histoMassGaussianMeson");
        histoMassCut[i]->SetName(Form("histoMassGaussianMeson_%s",cutStringsName[i].Data()));
        histoWidthCut[i]             = (TH1D*)Cutuncorrfile[i]->Get("histoFWHMMeson");
        histoWidthCut[i]->SetName(Form("histoWidthGaussianMeson_%s",cutStringsName[i].Data()));
        histoSBCut[i]               = (TH1D*)Cutuncorrfile[i]->Get("histoSBdefaultMeson");
        histoSBCut[i]->SetName(Form("histoSBdefaultMeson_%s",cutNumber[i].Data()));
        histoClusterE[i]       =(TH1D*)Cutuncorrfile[i]->Get("ClusterEPerEvent");
        if(histoClusterE[i]) histoClusterE[i]->SetName(Form("ClusterEPerEvent_%s", cutStringsName[i].Data()));
        else isClusterE = kFALSE;
        
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

        histoRatioRawYieldCut[i]    = (TH1D*) histoRawYieldCut[i]->Clone(Form("histoRatioRawYieldCut_%s",cutStringsName[i].Data()));
        histoRatioRawYieldCut[i]->Divide(histoRatioRawYieldCut[i],histoRawYieldCut[0],1.,1.,"B");
        histoRatioMassCut[i]        = (TH1D*) histoMassCut[i]->Clone(Form("histoRatioMassCut_%s",cutStringsName[i].Data()));
        histoRatioMassCut[i]->Divide(histoRatioMassCut[i],histoMassCut[0],1.,1.,"B");
        histoRatioWidthCut[i]        = (TH1D*) histoWidthCut[i]->Clone(Form("histoRatioWidthCut_%s",cutStringsName[i].Data()));
        histoRatioWidthCut[i]->Divide(histoRatioWidthCut[i],histoWidthCut[0],1.,1.,"B");

        histoRatioSBCut[i]          = (TH1D*) histoSBCut[i]->Clone(Form("histoRatioSBCut_%s", cutStringsName[i].Data()));
        histoRatioSBCut[i]->Divide(histoRatioSBCut[i],histoSBCut[0],1.,1.,"B");

        if(isClusterE){
          histoRatioClusterE[i]       =(TH1D*)histoClusterE[i]->Clone(Form("histoRatioClusterEPerEvent_%s", cutStringsName[i].Data()));
          histoRatioClusterE[i]->Divide(histoRatioClusterE[i],histoClusterE[0],1.,1.,"B");
        }
        
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
                if (cutVariationName.Contains("PhotonQuality")){
                    minYRatio = 0.001;
                    maxYRatio = 2;
                    padRawYieldRatios->SetLogy(1);
                } else if (cutVariationName.Contains("LHC12NL") || cutVariationName.Contains("LHC12-MultWeight")){
                    minYRatio = 0.9;
                    maxYRatio = 1.1;
                }
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
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55;
                if (mode != 0 && mode!= 1 ){
                    minYRatio = 0.75;
                    maxYRatio = 1.55;
                }
                if (cutVariationName.Contains("LHC12NL") || cutVariationName.Contains("LHC12-MultWeight")){
                    minYRatio = 0.81;
                    maxYRatio = 1.19;
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
//       if (mode == 2 || mode == 3 ) padTrueEffi->SetLogy(1);
//       else padTrueEffi->SetLogy(0);

      TLegend* legendEffiMeson = GetAndSetLegend2(0.15,0.93-0.03*NumberOfCuts,0.3,0.93, 1500*0.75*0.032); 
      for(Int_t i = 0; i< NumberOfCuts; i++){
          if(i == 0){
              DrawGammaSetMarker(histoTrueEffiCut[i], 20, 1., color[0], color[0]);
              DrawAutoGammaMesonHistos( histoTrueEffiCut[i],
                                      "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{%s}",textMeson.Data()),
                                      kTRUE, 5., 10e-10,kFALSE,
                                      kTRUE, -0.1, 0.00030,
                                      kFALSE, 0., 10.);
              if (mode == 9 || mode == 0 )histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(0.0,0.003);
//               if (mode == 2 || mode == 3 )histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(1.e-7,2);
              if (mode == 2 || mode == 3 )histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(0,0.07); 
              if (mode == 4 || mode == 5 )histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(0,0.7); 
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
              Double_t minYRatio = 0.45;
              Double_t maxYRatio = 1.55;
              if (cutVariationName.Contains("Weighting")){
                    minYRatio = 0.75;
                    maxYRatio = 1.25;
              }
              if (cutVariationName.Contains("LHC12NL") || cutVariationName.Contains("LHC12-MultWeight")){
                  minYRatio = 0.81;
                  maxYRatio = 1.19;
              }
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
          if (cutVariationName.Contains("Weighting")){              
            DrawGammaLines(0., maxPt, 1.05, 1.05, 0.1, kGray+1, 9);
            DrawGammaLines(0., maxPt, 0.95, 0.95, 0.1, kGray+1, 9);
          }    
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
//       if (mode == 2 || mode == 3 ) padAcceptance->SetLogy(1);
//       else padAcceptance->SetLogy(0);

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
                                      kFALSE, -0.1, 0.00030,
                                      kFALSE, 0., 10.);
//               if (mode == 9 || mode == 0 )histoAcceptanceCut[i]->GetYaxis()->SetRangeUser(0.0,0.003);
//               if (mode == 2 || mode == 3 )histoAcceptanceCut[i]->GetYaxis()->SetRangeUser(1.e-7,2);
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
              Double_t minYRatio = 0.85;
              Double_t maxYRatio = 1.15;
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
    //********************* Plotting Mass **************************************************
    //**************************************************************************************

    TCanvas* canvasMassMeson = new TCanvas("canvasMassMeson","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasMassMeson,  0.13, 0.02, 0.02, 0.09);
    // Define upper panel
    TPad* padMass = new TPad("padMass", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padMass, 0.12, 0.02, 0.02, 0.);
    padMass->Draw();
    // Define lower panel
    TPad* padMassRatios = new TPad("padMassRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padMassRatios, 0.12, 0.02, 0.0, 0.2);
    padMassRatios->Draw();

    // draw efficiency in upper panel
    padMass->cd();
    if ( (mode == 2 || mode == 3) && optionEnergy.CompareTo("8TeV")!=0 ) padMass->SetLogy(1);
    else padMass->SetLogy(0);

    TLegend* legendMass = GetAndSetLegend2(0.15,0.93-0.03*NumberOfCuts,0.3,0.93, 1500*0.75*0.032);
    for(Int_t i = 0; i< NumberOfCuts; i++){
        if(i == 0){
            DrawGammaSetMarker(histoMassCut[i], 20, 1., color[0], color[0]);
            Double_t massMin = 0.1275;
            if(mode == 2 || mode == 3 || mode == 4) massMin = 0.12;
            DrawAutoGammaMesonHistos( histoMassCut[i],
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("#it{m}_{%s} (GeV/c^{2})",textMeson.Data()),
                                    kFALSE, 5., 10e-10,kFALSE,
                                    kTRUE, massMin, 0.1425,
                                    kFALSE, 0., 10.);
            if (optionEnergy.CompareTo("8TeV") == 0) histoMassCut[i]->GetYaxis()->SetRangeUser(massMin,0.15);
            if (optionEnergy.CompareTo("8TeV") == 0 && cutVariationName.Contains("Calo") && cutVariationName.Contains("EGA")) histoMassCut[i]->GetYaxis()->SetRangeUser(0.12,0.20);
            if ( !meson.Contains("Pi0") )histoMassCut[i]->GetYaxis()->SetRangeUser(0.4,0.6);
            histoMassCut[i]->GetYaxis()->SetTitleOffset(1.4);
            histoMassCut[i]->DrawCopy("e1,p");
            legendMass->AddEntry(histoMassCut[i],Form("standard: %s",cutStringsName[i].Data()));
        } else {			
            if(i<20){
                DrawGammaSetMarker(histoMassCut[i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarker(histoMassCut[i], 20+i, 1.,color[i-20],color[i-20]);
            }
            histoMassCut[i]->DrawCopy("same,e1,p");
            legendMass->AddEntry(histoMassCut[i],cutStringsName[i].Data());
        }
    }
    legendMass->Draw();
    
    labelCollisionSystem4->Draw();		
    if (labelDetProcess2)labelDetProcess2->Draw();

    // Draw ratio of efficiencies in lower panel
    padMassRatios->cd();
    if( optionEnergy.Contains("Pb") ) padMassRatios->SetLogy(0);
    else padMassRatios->SetLogy(0);

    for(Int_t i = 0; i< NumberOfCuts; i++){
        if(i==0){
            Double_t minYRatio = 0.95;
            Double_t maxYRatio = 1.05;
            SetStyleHistoTH1ForGraphs(histoRatioMassCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.08, 0.1, 0.75, 0.5, 510,505);
            DrawGammaSetMarker(histoRatioMassCut[i], 20, 1.,color[0],color[0]);
            histoRatioMassCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
            histoRatioMassCut[i]->DrawCopy("p,e1");
        } else{
            if(i<20){
                DrawGammaSetMarker(histoRatioMassCut[i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarker(histoRatioMassCut[i], 20+i, 1.,color[i-20],color[i-20]);
            }
            histoRatioMassCut[i]->DrawCopy("same,e1,p");
        }
        DrawGammaLines(0., maxPt,1., 1.,0.1);
    }
    
    canvasMassMeson->Update();
    canvasMassMeson->SaveAs(Form("%s/%s_%s_Mass.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
    delete canvasMassMeson;

    //**************************************************************************************
    //********************* Plotting Width **************************************************
    //**************************************************************************************

    TCanvas* canvasWidthMeson = new TCanvas("canvasWidthMeson","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasWidthMeson,  0.13, 0.02, 0.02, 0.09);
    // Define upper panel
    TPad* padWidth = new TPad("padWidth", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padWidth, 0.12, 0.02, 0.02, 0.);
    padWidth->Draw();
    // Define lower panel
    TPad* padWidthRatios = new TPad("padWidthRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padWidthRatios, 0.12, 0.02, 0.0, 0.2);
    padWidthRatios->Draw();

    // draw efficiency in upper panel
    padWidth->cd();
    if ( (mode == 2 || mode == 3) && optionEnergy.CompareTo("8TeV")!=0 ) padWidth->SetLogy(1);
    else padWidth->SetLogy(0);

    TLegend* legendWidth = GetAndSetLegend2(0.15,0.93-0.03*NumberOfCuts,0.3,0.93, 1500*0.75*0.032);
    for(Int_t i = 0; i< NumberOfCuts; i++){
        if(i == 0){
            DrawGammaSetMarker(histoWidthCut[i], 20, 1., color[0], color[0]);
            DrawAutoGammaMesonHistos( histoWidthCut[i],
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("#it{#sigma}_{%s} (GeV/c^{2})",textMeson.Data()),
                                    kFALSE, 5., 10e-10,kFALSE,
                                    kTRUE, 0., 0.05,
                                    kFALSE, 0., 10.);
            if (optionEnergy.CompareTo("8TeV") == 0 && cutVariationName.Contains("Calo") && cutVariationName.Contains("EGA")) histoWidthCut[i]->GetYaxis()->SetRangeUser(0.,0.10);
            histoWidthCut[i]->GetYaxis()->SetTitleOffset(1.4);
            histoWidthCut[i]->DrawCopy("e1,p");
            legendWidth->AddEntry(histoWidthCut[i],Form("standard: %s",cutStringsName[i].Data()));
        } else {
            if(i<20){
                DrawGammaSetMarker(histoWidthCut[i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarker(histoWidthCut[i], 20+i, 1.,color[i-20],color[i-20]);
            }
            histoWidthCut[i]->DrawCopy("same,e1,p");
            legendWidth->AddEntry(histoWidthCut[i],cutStringsName[i].Data());
        }
    }
    legendWidth->Draw();

    labelCollisionSystem4->Draw();
    if (labelDetProcess2)labelDetProcess2->Draw();

    // Draw ratio of efficiencies in lower panel
    padWidthRatios->cd();
    if( optionEnergy.Contains("Pb") ) padWidthRatios->SetLogy(0);
    else padWidthRatios->SetLogy(0);

    for(Int_t i = 0; i< NumberOfCuts; i++){
        if(i==0){
            Double_t minYRatio = 0.8;
            Double_t maxYRatio = 1.2;
            SetStyleHistoTH1ForGraphs(histoRatioWidthCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.08, 0.1, 0.75, 0.5, 510,505);
            DrawGammaSetMarker(histoRatioWidthCut[i], 20, 1.,color[0],color[0]);
            histoRatioWidthCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
            histoRatioWidthCut[i]->DrawCopy("p,e1");
        } else{
            if(i<20){
                DrawGammaSetMarker(histoRatioWidthCut[i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarker(histoRatioWidthCut[i], 20+i, 1.,color[i-20],color[i-20]);
            }
            histoRatioWidthCut[i]->DrawCopy("same,e1,p");
        }
        DrawGammaLines(0., maxPt,1., 1.,0.1);
    }

    canvasWidthMeson->Update();
    canvasWidthMeson->SaveAs(Form("%s/%s_%s_Width.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
    delete canvasWidthMeson;

    //**************************************************************************************
    //************************ Plotting SB  ************************************************
    //**************************************************************************************

        TCanvas* canvasSBMeson = new TCanvas("canvasSBMeson","",1350,1500);
        DrawGammaCanvasSettings( canvasSBMeson,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padSB = new TPad("padSB", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padSB, 0.12, 0.02, 0.02, 0.);
        padSB->SetLogy();
        padSB->Draw();
        // lower pad definition
        TPad* padSBRatios = new TPad("padSBRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padSBRatios, 0.12, 0.02, 0.0, 0.2);
        padSBRatios->Draw();

        padSB->cd();
        padSB->SetTickx();
        padSB->SetTicky();

        // Plot SB in uppper panel
        padSB->cd();

        TLegend* legendSB = GetAndSetLegend2(0.15,0.93-0.03*NumberOfCuts,0.3,0.93, 1500*0.75*0.032);
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
                DrawGammaSetMarker(histoSBCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoSBCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s S/B",textMeson.Data()),
                                        kTRUE, 10., 1e-4, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                legendSB->AddEntry(histoSBCut[i],Form("standard: %s",cutStringsName[i].Data()));
            }
            else {
                if(i<20){
                    DrawGammaSetMarker(histoSBCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoSBCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoSBCut[i]->DrawCopy("same,e1,p");
                legendSB->AddEntry(histoSBCut[i],cutStringsName[i].Data());
            }

        }
        legendSB->Draw();
        padSBRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55; //qui
                if (cutVariationName.Contains("LHC12DistBC")){
                    minYRatio = 0.81;
                    maxYRatio = 1.29;
                }
                SetStyleHistoTH1ForGraphs(histoRatioSBCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioSBCut[i], 20, 1.,color[0],color[0]);
                histoRatioSBCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoRatioSBCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioSBCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioSBCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioSBCut[i]->DrawCopy("same,e1,p");
            }
            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasSBMeson->Update();
        canvasSBMeson->SaveAs(Form("%s/%s_%s_SB.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasSBMeson;


        //**************************************************************************************
        //************************ Plotting ClusterE *******************************************
        //**************************************************************************************

        if(isClusterE){
          TCanvas* canvasClusEMeson = new TCanvas("canvasClusEMeson","",1350,1500);
          DrawGammaCanvasSettings( canvasClusEMeson,  0.13, 0.02, 0.02, 0.09);
          // Upper pad definition
          TPad* padClusE = new TPad("padClusE", "", 0., 0.33, 1., 1.,-1, -1, -2);
          DrawGammaPadSettings( padClusE, 0.12, 0.02, 0.02, 0.);
          padClusE->SetLogy();
          padClusE->Draw();
          // lower pad definition
          TPad* padClusERatios = new TPad("padClusERatios", "", 0., 0., 1., 0.33,-1, -1, -2);
          DrawGammaPadSettings( padClusERatios, 0.12, 0.02, 0.0, 0.2);
          padClusERatios->Draw();

          padClusE->cd();
          padClusE->SetTickx();
          padClusE->SetTicky();

          // Plot ClusE in uppper panel
          padClusE->cd();

          TLegend* legendClusE = GetAndSetLegend2(0.35,0.93-0.03*NumberOfCuts,0.5,0.93, 1500*0.75*0.032);
          for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
              DrawGammaSetMarker(histoClusterE[i], 20, 1., color[0], color[0]);
              DrawAutoGammaMesonHistos( histoClusterE[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", "cluster E (GeV)",
                                        kTRUE, 10., 1e-10, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
              legendClusE->AddEntry(histoClusterE[i],Form("standard: %s",cutStringsName[i].Data()));
            }
            else {
              if(i<20){
                DrawGammaSetMarker(histoClusterE[i], 20+i, 1.,color[i],color[i]);
              } else {
                DrawGammaSetMarker(histoClusterE[i], 20+i, 1.,color[i-20],color[i-20]);
              }
              histoClusterE[i]->DrawCopy("same,e1,p");
              legendClusE->AddEntry(histoClusterE[i],cutStringsName[i].Data());
            }

          }
          legendClusE->Draw();
          padClusERatios->cd();
          for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
              // Set ratio min and max
              Double_t minYRatio = 0.45;
              Double_t maxYRatio = 1.55; //qui
              SetStyleHistoTH1ForGraphs(histoRatioClusterE[i], "#it{p}_{T} (GeV/#it{c})", "cluster E (GeV)", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
              DrawGammaSetMarker(histoRatioClusterE[i], 20, 1.,color[0],color[0]);
              histoRatioClusterE[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
              histoRatioClusterE[i]->DrawCopy("p,e1");
            } else{
              if(i<20){
                DrawGammaSetMarker(histoRatioClusterE[i], 20+i, 1.,color[i],color[i]);
              } else {
                DrawGammaSetMarker(histoRatioClusterE[i], 20+i, 1.,color[i-20],color[i-20]);
              }
              histoRatioClusterE[i]->DrawCopy("same,e1,p");
            }
            DrawGammaLines(0., maxPt,1., 1.,0.1);
          }

          canvasClusEMeson->Update();
          canvasClusEMeson->SaveAs(Form("%s/%s_%s_ClusE.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
          delete canvasClusEMeson;
        }


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
