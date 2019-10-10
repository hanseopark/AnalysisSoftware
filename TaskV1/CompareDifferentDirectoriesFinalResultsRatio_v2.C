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

TGraphAsymmErrors* ScaleGraphAsym (TGraphAsymmErrors* graph, Double_t scaleFac){
	TGraphAsymmErrors* dummyGraph = (TGraphAsymmErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));
	Double_t * xValue         = dummyGraph->GetX();
	Double_t * yValue         = dummyGraph->GetY();
	Double_t* xErrorLow       = dummyGraph->GetEXlow();
	Double_t* xErrorHigh      = dummyGraph->GetEXhigh();
	Double_t* yErrorLow       = dummyGraph->GetEYlow();
	Double_t* yErrorHigh      = dummyGraph->GetEYhigh();
	Int_t nPoints             = dummyGraph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i]               = yValue[i]*scaleFac;
		yErrorLow[i]            = yErrorLow[i]*scaleFac;
		yErrorHigh[i]           = yErrorHigh[i]*scaleFac;
	}
	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
	return returnGraph;
}

TGraphAsymmErrors * DivideTGraphAsymErrorByTGraphAsymError(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2, TString nameGraph){
	TGraphAsymmErrors* dummyGraph = (TGraphAsymmErrors*)graph1->Clone(Form("%s_R1",graph1->GetName()));
	Int_t nPoints              = dummyGraph->GetN();
	Double_t * xValue1         = dummyGraph->GetX();
	Double_t * yValue1         = dummyGraph->GetY();
	Double_t* xErrorLow1       = dummyGraph->GetEXlow();
	Double_t* xErrorHigh1      = dummyGraph->GetEXhigh();
	Double_t* yErrorLow1       = dummyGraph->GetEYlow();
	Double_t* yErrorHigh1      = dummyGraph->GetEYhigh();

	Double_t * yValue2         = graph2->GetY();
	Double_t* yErrorLow2       = graph2->GetEYlow();
	Double_t* yErrorHigh2      = graph2->GetEYhigh();

	for (Int_t i = 0; i < nPoints; i++){
        if(yValue2[i]!=0)
		    yValue1[i]               = yValue1[i]/yValue2[i];
        else
		    yValue1[i]               = 0;
		// yErrorLow1[i]            = yValue1[i]*TMath::Sqrt( TMath::Power(yErrorLow1[i]/yValue1[i],2) + TMath::Power(yErrorLow2[i]/yValue1[i],2));
		yErrorHigh1[i]           = yValue1[i]*TMath::Sqrt( TMath::Power(yErrorHigh1[i]/yValue1[i],2) + TMath::Power(yErrorHigh2[i]/yValue2[i],2));
	}

	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValue1,yValue1,xErrorLow1,xErrorHigh1,yErrorHigh1,yErrorHigh1);
    returnGraph->SetName(nameGraph.Data());
	return returnGraph;
}

void CompareDifferentDirectoriesFinalResultsRatio_v2(
    TString FolderList          = "",
    TString suffix              = "gif",
    TString meson               = "",
    Bool_t kIsMC                = kFALSE,
    TString optionEnergy1       = "",
    TString optionEnergy2       = "",
    Int_t NumberOfCuts          = 1,
    TString optionPeriod        = "No",
    Int_t mode                  = 10,
    TString cutVariationName    = "NonLinearity",
    TString triggerName1         = "INT7",
    TString triggerName2         = "INT7",
    Bool_t setFullPathInInputFile   = kFALSE
){

    // if (!(mode == 10 || mode == 11 )){
    //     cout << "incorrect mode: " << mode << endl;
    //     return ;
    // }
    TString     optionEnergy[2];
    optionEnergy[0]=optionEnergy1;
    optionEnergy[1]=optionEnergy2;
    // Initialize arrays
    TString     filePathInput[2][50];
    TString     directoryInFile[2][50];
    TString     cutNumber[2][50];
    TString     cutStringsName[50];
    Double_t    scaleFacDenom[50];

    // prepare nice plotting algorithms
    StyleSettingsThesis();
    SetPlotStyle();

    if(triggerName1.CompareTo("") && triggerName2.CompareTo("")){
        cutVariationName = Form("%s%s", cutVariationName.Data(),triggerName1.Data());
    }

    // Set output folder
    TString outputDir = Form("CutStudies/%sOver%s/%s",optionEnergy[0].Data(),optionEnergy[1].Data(), cutVariationName.Data());

    TString outputDirRootFile = Form("CutStudies/%sOver%s",optionEnergy[0].Data(),optionEnergy[1].Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // Set meson name for plotting
    TString textMeson;
    if (meson.Contains("Pi0")){
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
    TString collisionSystem1= ReturnFullCollisionsSystem(optionEnergy[0]);
    if (collisionSystem1.CompareTo("") == 0){
        cout << "No correct collision system (1) specification, has been given" << endl;
        return;
    }
    // Determine collsision system string
    TString collisionSystem2= ReturnFullCollisionsSystem(optionEnergy[1]);
    if (collisionSystem2.CompareTo("") == 0){
        cout << "No correct collision system (2) specification, has been given" << endl;
        return;
    }

    TString detProcess                                      = ReturnFullTextReconstructionProcess(mode, 0, textMeson);

    // Define colors for comparisons
    Color_t color[20] = {kBlack, kAzure, kGreen+2,kOrange+2,kRed, kViolet,  kBlue-9, kSpring+10,
                        kCyan+3, kCyan-10, kCyan, kGreen+4, kGreen-9,
                        kGreen,  kYellow+4, kYellow+3, kMagenta+4,
                        kMagenta-8, kGray, kGray+3};

    // read folder and name from file
    ifstream in(FolderList.Data());
    cout<<"Available Cuts:"<<endl;
    TString folderName[2];
    TString dirInFileName[2];
    TString cutName;
    Double_t scalingfactor;
    TString cutNr[2];
    Int_t Number = 0;
    while(!in.eof() ){
        in >> folderName[0] >> cutNr[0] >> dirInFileName[0] >> folderName[1] >> cutNr[1] >> dirInFileName[1] >> scalingfactor >> cutName;
        cutName.ReplaceAll("_"," ");
        filePathInput[0][Number] = folderName[0];
        filePathInput[1][Number] = folderName[1];
        cutNumber[0][Number] = cutNr[0];
        cutNumber[1][Number] = cutNr[1];
        directoryInFile[0][Number] = dirInFileName[0];
        directoryInFile[1][Number] = dirInFileName[1];
        cutStringsName[Number] = cutName;
        scaleFacDenom[Number] = scalingfactor;
        cout << "Set " << Number << ":" << endl;
        cout<< "\t" <<filePathInput[0][Number]<< "\t" << cutNumber[0][Number]<< "\t" << directoryInFile[0][Number]<< "\t"<< cutStringsName[Number] <<endl;
        cout<< "\t" <<filePathInput[1][Number]<< "\t" << cutNumber[1][Number]<< "\t" << directoryInFile[1][Number]<< "\t"<< scaleFacDenom[Number]<< "\t"<< cutStringsName[Number] <<endl;
        Number++;
    }
    cout<<"=========================="<<endl;

    // Definition of necessary graphgram arrays
    const Int_t ConstNumberOfCuts = 20;

    // TFile *fileResoCorr[2];
    // if(cutVariationName.Contains("Resolution")){
    //     if(optionEnergy1.CompareTo("pPb_8TeV"))
    //         fileResoCorr[0] = new TFile("/media/nschmidt/local/ANALYSIS/pPb_8TeV_mEDC_Resolution/pPb8TeV/pdf/FractionsCheck/ToyMCOutput_Pi0_pPb8TeV.root");
    //     if(optionEnergy2.CompareTo("8TeV"))
    //         fileResoCorr[1] = new TFile("/media/nschmidt/local/ANALYSIS/pp_8TeV_mEMC_Resolution/pp8TeV/pdf/FractionsCheck/ToyMCOutput_Pi0_pp8TeV.root");
    //     if (fileResoCorr[0]->IsZombie() || fileResoCorr[1]->IsZombie()){
    //         cout << "could not find resolution systematics inputs!... returning!!" << endl;
    //         return;
    //     }
    // }



    TFile *CutFinalResFile[2][ConstNumberOfCuts];

    TH1D *histoCorrectedYieldCutDummy[2][ConstNumberOfCuts] = {{NULL}};
    TGraphAsymmErrors *graphCorrectedYieldCutDummy[2][ConstNumberOfCuts] = {{NULL}};
    TH1D *histoRatioCorrectedYieldCutDummy[2][ConstNumberOfCuts] = {{NULL}};
    TGraphAsymmErrors *graphRatioCorrectedYieldCutDummy[2][ConstNumberOfCuts] = {{NULL}};

    TH1D *histoNuclModFacCut[ConstNumberOfCuts] = {NULL};
    TGraphAsymmErrors *graphNuclModFacCut[ConstNumberOfCuts] = {NULL};

    TH1D *histoRatioCorrectedYieldCut[ConstNumberOfCuts] = {NULL};
    TGraphAsymmErrors *graphRatioCorrectedYieldCut[ConstNumberOfCuts] = {NULL};

    Double_t minPt  = 0;
    Double_t maxPt  = 0;
    TString nameCorrectedYieldHisto[2];
    TString nameCorrectedYield[2];
    for (Int_t i=0; i< NumberOfCuts; i++){
        // Set correct graphgram name for corrected yield and efficiency
        nameCorrectedYieldHisto[0]                     = Form("CorrectedYield%s",meson.Data());
        nameCorrectedYieldHisto[1]                     = Form("CorrectedYield%s",meson.Data());
        if(triggerName1.CompareTo("") && triggerName2.CompareTo("")){
            nameCorrectedYieldHisto[0]                     = Form("RAWYieldPerEvents%s_%s",meson.Data(),triggerName1.Data());
            nameCorrectedYieldHisto[1]                     = Form("RAWYieldPerEvents%s_%s",meson.Data(),triggerName2.Data());
        }
        nameCorrectedYield[0]                          = Form("graphCorrectedYield%s",meson.Data());
        nameCorrectedYield[1]                          = Form("graphCorrectedYield%s",meson.Data());
        if(triggerName1.CompareTo("") && triggerName2.CompareTo("")){
            nameCorrectedYield[0]                     = Form("CorrectedYield%s_%s",meson.Data(),triggerName1.Data());
            nameCorrectedYield[1]                     = Form("CorrectedYield%s_%s",meson.Data(),triggerName2.Data());
        }
        for (Int_t j=0; j< 2; j++){
            // read file with corrections
            // cout<< filePathInput[j][i] << endl;
            CutFinalResFile[j][i] = new TFile(filePathInput[j][i]);
            if (CutFinalResFile[j][i]->IsZombie()) return;
             gSystem->Exec(Form("cp %s %s/Input_%d_%d.root", filePathInput[j][i].Data(), outputDir.Data(), j , i));

            // Read graphgrams and rename them from the original files for each cut
            histoCorrectedYieldCutDummy[j][i]   = (TH1D*)CutFinalResFile[j][i]->Get(Form("%s%s/%s",meson.Data(),directoryInFile[j][i].Data(),nameCorrectedYieldHisto[j].Data()));
            graphCorrectedYieldCutDummy[j][i]   = (TGraphAsymmErrors*)CutFinalResFile[j][i]->Get(Form("%s%s/%s",meson.Data(),directoryInFile[j][i].Data(),nameCorrectedYield[j].Data()));
            histoCorrectedYieldCutDummy[j][i]->Sumw2();
            if(j==1)histoCorrectedYieldCutDummy[j][i]->Scale(scaleFacDenom[i]);
            if(j==1)graphCorrectedYieldCutDummy[j][i] = ScaleGraphAsym(graphCorrectedYieldCutDummy[j][i],scaleFacDenom[i]);
            // graphCorrectedYieldCutDummy[j][i]->Print();
            // histoCorrectedYieldCutDummy[j][i]->Print("all");
            for(Int_t bin=1;bin<histoCorrectedYieldCutDummy[j][i]->GetNbinsX();bin++){
                Double_t currentXValue = histoCorrectedYieldCutDummy[j][i]->GetBinCenter(bin);
                // cout << "currentXValue: " << currentXValue << endl;
                for (Int_t ibin = 0; ibin < graphCorrectedYieldCutDummy[j][i]->GetN(); ibin++){
                    if(graphCorrectedYieldCutDummy[j][i]->GetX()[ibin] == currentXValue){
                        // cout << "found bin in graph, setting bin " << bin << " to graphvalues: " << ibin << "\t" << graphCorrectedYieldCutDummy[j][i]->GetY()[ibin] << "\t" << graphCorrectedYieldCutDummy[j][i]->GetEYlow()[ibin] << endl;
                        histoCorrectedYieldCutDummy[j][i]->SetBinContent(bin, graphCorrectedYieldCutDummy[j][i]->GetY()[ibin]);
                        histoCorrectedYieldCutDummy[j][i]->SetBinError(bin, graphCorrectedYieldCutDummy[j][i]->GetEYlow()[ibin]);
                    }
                }
            }
            // cout << "REPRINT!!:"<<endl;
            // graphCorrectedYieldCutDummy[j][i]->Print();
            // histoCorrectedYieldCutDummy[j][i]->Print("all");

        }
        // Read graphgrams and rename them from the original files for each cut
        if(histoCorrectedYieldCutDummy[0][i] && histoCorrectedYieldCutDummy[1][i]){
            histoNuclModFacCut[i] = (TH1D*) histoCorrectedYieldCutDummy[0][i]->Clone(Form("histoCorrectedYieldCut%d",i));
            histoNuclModFacCut[i]->Divide(histoCorrectedYieldCutDummy[0][i], histoCorrectedYieldCutDummy[1][i]);

            histoRatioCorrectedYieldCutDummy[0][i] = (TH1D*) histoCorrectedYieldCutDummy[0][i]->Clone(Form("histoRatioCorrectedYieldCutDummy0%d",i));
            histoRatioCorrectedYieldCutDummy[0][i]->Divide(histoCorrectedYieldCutDummy[0][i], histoCorrectedYieldCutDummy[0][0],1.,1.,"B");
            histoRatioCorrectedYieldCutDummy[1][i] = (TH1D*) histoCorrectedYieldCutDummy[0][i]->Clone(Form("histoRatioCorrectedYieldCutDummy1%d",i));
            histoRatioCorrectedYieldCutDummy[1][i]->Divide(histoCorrectedYieldCutDummy[1][i], histoCorrectedYieldCutDummy[1][0],1.,1.,"B");
        }
        if(graphCorrectedYieldCutDummy[0][i] && graphCorrectedYieldCutDummy[1][i])
            graphNuclModFacCut[i]   = DivideTGraphAsymErrorByTGraphAsymError(graphCorrectedYieldCutDummy[0][i], graphCorrectedYieldCutDummy[1][i], Form("%s_%s",nameCorrectedYield[1].Data(),cutStringsName[i].Data()));

        // Calculate ratios for comparisons
        if(histoNuclModFacCut[i] && histoNuclModFacCut[0]){
            histoRatioCorrectedYieldCut[i] = (TH1D*) histoNuclModFacCut[i]->Clone(Form("histoRatioCorrectedYieldCut%d",i));
            histoRatioCorrectedYieldCut[i]->Divide(histoNuclModFacCut[i], histoNuclModFacCut[0],1.,1.,"B");
        }
        if(graphNuclModFacCut[i] && graphNuclModFacCut[0]){
            graphRatioCorrectedYieldCut[i] = DivideTGraphAsymErrorByTGraphAsymError(graphNuclModFacCut[i], graphNuclModFacCut[0], Form("graphRatioCorrectedYieldCut_%s",cutStringsName[i].Data()));
            if (i > 0){
                minPt= graphNuclModFacCut[i]->GetX()[0] - 3*graphNuclModFacCut[i]->GetEXlow()[0];
                maxPt= graphNuclModFacCut[i]->GetX()[graphNuclModFacCut[i]->GetN()-1] + 3*graphNuclModFacCut[i]->GetEXhigh()[graphNuclModFacCut[i]->GetN()-1];
            }
        }
    }
    cout<<"=========================="<<endl;

    //*****************************************************************************************
    //******************* Compare Corrected Yields ********************************************
    //*****************************************************************************************
    // Define canvas
    TCanvas* canvasCorrectedYieldMeson = new TCanvas("canvasCorrectedYieldMeson","",1350,2000);
    DrawGammaCanvasSettings( canvasCorrectedYieldMeson,  0.13, 0.02, 0.02, 0.09);
    // Define upper panel
    TPad* padCorrectedYield = new TPad("padCorrectedYield", "", 0., 0.5, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYield, 0.12, 0.02, 0.02, 0.0);
    padCorrectedYield->SetLogx(1);
    padCorrectedYield->SetLogy(1);
    padCorrectedYield->Draw();
    // Define middle panel
    TPad* padCorrectedYieldRatios = new TPad("padCorrectedYieldRatios", "", 0., 0.27, 1., 0.5,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYieldRatios, 0.12, 0.02, 0.0, 0.00);
    padCorrectedYieldRatios->SetLogy(0);
    padCorrectedYieldRatios->SetLogx(1);
    padCorrectedYieldRatios->Draw();
    // Define lower panel
    TPad* padCorrectedYieldRatiosRatios = new TPad("padCorrectedYieldRatiosRatios", "", 0., 0., 1., 0.27,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYieldRatiosRatios, 0.12, 0.02, 0.0, 0.2);
    padCorrectedYieldRatiosRatios->SetLogy(0);
    padCorrectedYieldRatiosRatios->SetLogx(1);
    padCorrectedYieldRatiosRatios->Draw();

    TH2F * histo2DDummyCorrYield;
    histo2DDummyCorrYield                            = new TH2F("histo2DDummyCorrYield","histo2DDummyCorrYield",11000,minPt,maxPt,1000,graphCorrectedYieldCutDummy[0][0]->GetY()[graphCorrectedYieldCutDummy[0][0]->GetN()-1]*0.1, graphCorrectedYieldCutDummy[0][0]->GetY()[0]*10);
    SetStyleHistoTH2ForGraphs(histo2DDummyCorrYield, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)",0.035,0.04, 0.035,0.04, 1.,1.);
    // Plot corrected yield in upper panel
    padCorrectedYield->cd();
    TLegend* legendCorrectedYieldMeson = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
    histo2DDummyCorrYield->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DDummyCorrYield->Draw("copy");
    for(Int_t i = 0; i< NumberOfCuts; i++){
        for (Int_t j=0; j< 2; j++){
            if(i<20){
                DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldCutDummy[j][i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldCutDummy[j][i], 20+i, 1.,color[i-20],color[i-20]);
            }
            graphCorrectedYieldCutDummy[j][i]->Draw("e1,p,same");
            if(j==0)
                legendCorrectedYieldMeson->AddEntry(graphCorrectedYieldCutDummy[j][i], Form("%s",cutStringsName[i].Data()));
        }
    }
    legendCorrectedYieldMeson->Draw();
    // Labeling of plot
    TLatex *labelCollisionSystem3 = new TLatex(0.55,0.91,collisionSystem1.Data());
    SetStyleTLatex( labelCollisionSystem3, 0.038,4);
    labelCollisionSystem3->Draw();
    if (optionPeriod.CompareTo("No") != 0){
        TLatex *labelPeriod = new TLatex(0.55,0.86,optionPeriod.Data());
        SetStyleTLatex( labelPeriod, 0.038,4);
        labelPeriod->Draw();
    }
    TH2F * histo2DDummyCorrYieldRatio;
    histo2DDummyCorrYieldRatio                            = new TH2F("histo2DDummyCorrYieldRatio","histo2DDummyCorrYieldRatio",11000,minPt,maxPt,1000,0.86, 1.14);
    SetStyleHistoTH2ForGraphs(histo2DDummyCorrYieldRatio, "#it{p}_{T} (GeV/#it{c})", "pPb #frac{modified}{standard}",0.035,0.04, 0.08,0.09, 1.,0.5);
    // Plot corrected yield in upper panel
    padCorrectedYieldRatios->cd();
    for(Int_t i = 0; i< NumberOfCuts; i++){
        if(i == 0){
            histo2DDummyCorrYieldRatio->Draw("copy");
            DrawGammaSetMarker(histoRatioCorrectedYieldCutDummy[0][i], 20, 1., color[0], color[0]);
            histoRatioCorrectedYieldCutDummy[0][i]->Draw("same,e1,p");
        }
        else{
            if(i<20){
                DrawGammaSetMarker(histoRatioCorrectedYieldCutDummy[0][i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarker(histoRatioCorrectedYieldCutDummy[0][i], 20+i, 1.,color[i-20],color[i-20]);
            }
            histoRatioCorrectedYieldCutDummy[0][i]->Draw("e1,p,same");
        }

    }
    DrawGammaLines(0., maxPt,1., 1.,1,  kGray+2, 2);
    DrawGammaLines(0., maxPt, 1.1, 1.1, 1, kGray+1, 7);
    DrawGammaLines(0., maxPt, 0.9, 0.9, 1, kGray+1, 7);

    // if (labelDetProcess) labelDetProcess->Draw();
    // plot ratio of corrected yields in lower panel
    padCorrectedYieldRatiosRatios->cd();
    TH2F * histo2DDummyCorrYieldRatioRatio;
    histo2DDummyCorrYieldRatioRatio                            = new TH2F("histo2DDummyCorrYieldRatioRatio","histo2DDummyCorrYieldRatioRatio",11000,minPt,maxPt,1000,0.86, 1.14);
    SetStyleHistoTH2ForGraphs(histo2DDummyCorrYieldRatioRatio, "#it{p}_{T} (GeV/#it{c})", "pp #frac{modified}{standard}",0.07,0.08, 0.08,0.09, 1.,0.5);
    for(Int_t i = 0; i< NumberOfCuts; i++){
        if(i==0){
            // Set ratio min and max
            Double_t minYRatio = 0.81;
            Double_t maxYRatio = 1.19;
            histo2DDummyCorrYieldRatioRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
            histo2DDummyCorrYieldRatioRatio->GetXaxis()->SetNoExponent(kTRUE);
            histo2DDummyCorrYieldRatioRatio->Draw("copy");
            DrawGammaSetMarkerTGraphAsym(graphRatioCorrectedYieldCut[i], 20, 1.,color[0],color[0]);
            // graphRatioCorrectedYieldCut[i]->Draw("e1,p,same");
            DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20, 1.,color[0],color[0]);
            histoRatioCorrectedYieldCut[i]->Draw("same,e1,p");
        }
        else{
            if(i<20){
                DrawGammaSetMarker(histoRatioCorrectedYieldCutDummy[1][i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarker(histoRatioCorrectedYieldCutDummy[1][i], 20+i, 1.,color[i-20],color[i-20]);
            }
            histoRatioCorrectedYieldCutDummy[1][i]->Draw("same,e1,p");
        }
    }
    DrawGammaLines(0., maxPt,1., 1.,1,  kGray+2, 2);
    DrawGammaLines(0., maxPt, 1.1, 1.1, 1, kGray+1, 7);
    DrawGammaLines(0., maxPt, 0.9, 0.9, 1, kGray+1, 7);
    canvasCorrectedYieldMeson->Update();
    canvasCorrectedYieldMeson->SaveAs(Form("%s/%s_%s_CorrectedYield.%s",outputDir.Data(), meson.Data(),prefix2.Data(),suffix.Data()));
    delete canvasCorrectedYieldMeson;
   

    //*****************************************************************************************
    //******************* Compare Corrected Yields ********************************************
    //*****************************************************************************************
    // Define canvas
    TCanvas* canvasNuclModFacMeson = new TCanvas("canvasNuclModFacMeson","",1350,1500);
    DrawGammaCanvasSettings( canvasNuclModFacMeson,  0.13, 0.02, 0.02, 0.09);
    // Define upper panel
    TPad* padNuclModFac = new TPad("padNuclModFac", "", 0., 0.5, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padNuclModFac, 0.12, 0.02, 0.02, 0.0);
    padNuclModFac->SetLogx(1);
    padNuclModFac->SetLogy(1);
    padNuclModFac->Draw();
    // Define middle panel
    TPad* padNuclModFacRatios = new TPad("padNuclModFacRatios", "", 0., 0.25, 1., 0.5,-1, -1, -2);
    DrawGammaPadSettings( padNuclModFacRatios, 0.12, 0.02, 0.0, 0.00);
    padNuclModFacRatios->SetLogy(0);
    padNuclModFacRatios->SetLogx(1);
    padNuclModFacRatios->Draw();
    // Define lower panel
    TPad* padNuclModFacRatiosRatios = new TPad("padNuclModFacRatiosRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padNuclModFacRatiosRatios, 0.12, 0.02, 0.0, 0.2);
    padNuclModFacRatiosRatios->SetLogy(0);
    padNuclModFacRatiosRatios->SetLogx(1);
    padNuclModFacRatiosRatios->Draw();

    TH2F * histo2DDummyNuclModFac;
    histo2DDummyNuclModFac                            = new TH2F("histo2DDummyNuclModFac","histo2DDummyNuclModFac",11000,minPt,maxPt,1000,graphCorrectedYieldCutDummy[0][0]->GetY()[graphCorrectedYieldCutDummy[0][0]->GetN()-1]*0.1, graphCorrectedYieldCutDummy[0][0]->GetY()[0]*10);
    SetStyleHistoTH2ForGraphs(histo2DDummyNuclModFac, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)",0.035,0.04, 0.035,0.04, 1.,1.);
    // Plot corrected yield in upper panel
    padNuclModFac->cd();
    TLegend* legendNuclModFacMeson = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
    histo2DDummyNuclModFac->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DDummyNuclModFac->Draw("copy");
    for(Int_t i = 0; i< NumberOfCuts; i++){
        for (Int_t j=0; j< 2; j++){
            if(i<20){
                DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldCutDummy[j][i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldCutDummy[j][i], 20+i, 1.,color[i-20],color[i-20]);
            }
            graphCorrectedYieldCutDummy[j][i]->Draw("e1,p,same");
            if(j==0)
                legendNuclModFacMeson->AddEntry(graphCorrectedYieldCutDummy[j][i], Form("%s",cutStringsName[i].Data()));
        }
    }
    legendNuclModFacMeson->Draw();
    // Labeling of plot
    labelCollisionSystem3 = new TLatex(0.55,0.91,collisionSystem1.Data());
    SetStyleTLatex( labelCollisionSystem3, 0.038,4);
    labelCollisionSystem3->Draw();
    if (optionPeriod.CompareTo("No") != 0){
        TLatex *labelPeriod = new TLatex(0.55,0.86,optionPeriod.Data());
        SetStyleTLatex( labelPeriod, 0.038,4);
        labelPeriod->Draw();
    }
    TH2F * histo2DDummyNuclModFacRatio;
    histo2DDummyNuclModFacRatio                            = new TH2F("histo2DDummyNuclModFacRatio","histo2DDummyNuclModFacRatio",11000,minPt,maxPt,1000,0.21, 1.49);
    SetStyleHistoTH2ForGraphs(histo2DDummyNuclModFacRatio, "#it{p}_{T} (GeV/#it{c})", "#it{R}_{pA}",0.035,0.04, 0.08,0.09, 1.,0.5);
    // Plot corrected yield in upper panel
    padNuclModFacRatios->cd();
    for(Int_t i = 0; i< NumberOfCuts; i++){
        if(i == 0){
            histo2DDummyNuclModFacRatio->Draw("copy");
            // DrawGammaSetMarkerTGraphAsym(graphNuclModFacCut[i], 20, 1., color[0], color[0]);
            // graphNuclModFacCut[i]->Draw("e1,p,same");
            DrawGammaSetMarker(histoNuclModFacCut[i], 20, 1., color[0], color[0]);
            histoNuclModFacCut[i]->Draw("e1,p,same");
        }
        else{
            if(i<20){
                // DrawGammaSetMarkerTGraphAsym(graphNuclModFacCut[i], 20+i, 1.,color[i],color[i]);
                DrawGammaSetMarker(histoNuclModFacCut[i], 20+i, 1.,color[i],color[i]);
            } else {
                // DrawGammaSetMarkerTGraphAsym(graphNuclModFacCut[i], 20+i, 1.,color[i-20],color[i-20]);
                DrawGammaSetMarker(histoNuclModFacCut[i], 20+i, 1.,color[i-20],color[i-20]);
            }
            // graphNuclModFacCut[i]->Draw("e1,p,same");
            histoNuclModFacCut[i]->Draw("e1,p,same");
        }

    }

    // if (labelDetProcess) labelDetProcess->Draw();
    // plot ratio of corrected yields in lower panel
    padNuclModFacRatiosRatios->cd();
    TH2F * histo2DDummyNuclModFacRatioRatio;
    histo2DDummyNuclModFacRatioRatio                            = new TH2F("histo2DDummyNuclModFacRatioRatio","histo2DDummyNuclModFacRatioRatio",11000,minPt,maxPt,1000,0.4, 1.6);
    SetStyleHistoTH2ForGraphs(histo2DDummyNuclModFacRatioRatio, "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}",0.07,0.08, 0.08,0.09, 1.,0.5);
    for(Int_t i = 0; i< NumberOfCuts; i++){
        if(i==0){
            // Set ratio min and max
            Double_t minYRatio = 0.81;
            Double_t maxYRatio = 1.19;
            histo2DDummyNuclModFacRatioRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
            histo2DDummyNuclModFacRatioRatio->GetXaxis()->SetNoExponent(kTRUE);
            histo2DDummyNuclModFacRatioRatio->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
            histo2DDummyNuclModFacRatioRatio->Draw("copy");
            DrawGammaSetMarkerTGraphAsym(graphRatioCorrectedYieldCut[i], 20, 1.,color[0],color[0]);
            // graphRatioCorrectedYieldCut[i]->Draw("e1,p,same");
            DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20, 1.,color[0],color[0]);
            histoRatioCorrectedYieldCut[i]->Draw("same,e1,p");
        }
        else{
            if(i<20){
                DrawGammaSetMarkerTGraphAsym(graphRatioCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
                DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarkerTGraphAsym(graphRatioCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
            }
            // graphRatioCorrectedYieldCut[i]->Draw("e1,p,same");
            histoRatioCorrectedYieldCut[i]->Draw("same,e1,p");
        }
    }
    DrawGammaLines(0., maxPt,1., 1.,1,  kGray+2, 2);
    DrawGammaLines(0., maxPt, 1.1, 1.1, 1, kGray+1, 7);
    DrawGammaLines(0., maxPt, 0.9, 0.9, 1, kGray+1, 7);
    canvasNuclModFacMeson->Update();
    canvasNuclModFacMeson->SaveAs(Form("%s/%s_%s_NuclModFactor.%s",outputDir.Data(), meson.Data(),prefix2.Data(),suffix.Data()));
    delete canvasNuclModFacMeson;
   

  //d*************************************************************************************************
  //d******************** Output of the systematic Error due to Signal extraction for Pi0 ************
  //d*************************************************************************************************
    // Determine number of bins
    Int_t NBinsPt = graphNuclModFacCut[0]->GetN();
    const Int_t NBinstPtConst = NBinsPt+1;

    // Create array of bin boundaries
    Double_t  BinsXCenter[NBinstPtConst];
    Double_t  BinsXWidth[NBinstPtConst];
    BinsXCenter[0] = 0;
    BinsXWidth[0]=0.;
    for (Int_t i = 0; i < NBinsPt; i++){
        BinsXCenter[i] = graphNuclModFacCut[0]->GetX()[i];
        BinsXWidth[i]= graphNuclModFacCut[0]->GetEXlow()[i];
    }

    // Create array of Sys Err Objects and fill them
    SysErrorConversion SysErrCut[ConstNumberOfCuts][NBinstPtConst];
    for (Int_t j = 0; j < NumberOfCuts; j++){
        for (Int_t i = 0; i < NBinsPt; i++){
            SysErrCut[j][i].value = graphNuclModFacCut[j]->GetY()[i];
            SysErrCut[j][i].error = graphNuclModFacCut[j]->GetEYlow()[i];
        }
    }

    // Create Difference arrays
    Double_t DifferenceCut[ConstNumberOfCuts][NBinstPtConst];
    Double_t DifferenceErrorCut[ConstNumberOfCuts][NBinstPtConst];
    Double_t RelDifferenceCut[ConstNumberOfCuts][NBinstPtConst];
    Double_t RelDifferenceErrorCut[ConstNumberOfCuts][NBinstPtConst];

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
        for ( Int_t i = 0; i < NBinstPtConst; i++) {
            DifferenceCut[j][i]=0.;
            DifferenceErrorCut[j][i]=0.;
            LargestDiffNeg[i]=0.;
            LargestDiffPos[i]=0.;
            LargestDiffErrorNeg[i]=0.;
            LargestDiffErrorPos[i]=0.;
            RelDifferenceCut[j][i]=0.;
            RelDifferenceErrorCut[j][i]=0.;
        }
    }

    // Calculate largest difference among cut variation
    for(Int_t j = 1; j < NumberOfCuts; j++){
        for (Int_t i = 0; i < NBinsPt; i++){
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
            // Calculate largest differences in positiv and negative direction
            if(DifferenceCut[j][i] < 0){ // largest negativ deviation
                // Take deviation if larger than previous largest deviation
                // and relative raw yield loss less than 75%
                if (TMath::Abs(LargestDiffNeg[i]) < TMath::Abs(DifferenceCut[j][i])){
                    LargestDiffNeg[i] = DifferenceCut[j][i];
                    LargestDiffErrorNeg[i] = DifferenceErrorCut[j][i];
                }
            } else { // largest positive deviation
                // Take deviation if larger than previous largest deviation
                // and relative raw yield loss less than 75%
                if (TMath::Abs(LargestDiffPos[i]) < TMath::Abs(DifferenceCut[j][i])){
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
            for(Int_t i = 0; i < (NBinsPt); i++){
                SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << endl;
            }
        } else{
            SysErrDat << endl <<"Bin" << "\t" << cutNumber[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
            for(Int_t i = 0; i < (NBinsPt); i++){
                SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << "\t" <<  DifferenceCut[l][i] << "\t"<< DifferenceErrorCut[l][i] << "\t"<< RelDifferenceCut[l][i] <<  "\t" << RelDifferenceErrorCut[l][i] << endl;
            }
        }
    }
    SysErrDat << endl;
    SysErrDat << endl;
    SysErrDat << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
    for(Int_t i = 0; i < (NBinsPt); i++){
        SysErrDat << BinsXCenter[i]  << "\t" << LargestDiffNeg[i] << "\t" <<LargestDiffErrorNeg[i]<< "\t" << LargestDiffPos[i] << "\t" << LargestDiffErrorPos[i]<<endl;
    }
    SysErrDat << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
    // Calculate largest relative deviations
    for(Int_t i = 0; i < (NBinsPt); i++){
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
    TGraphAsymmErrors* SystErrGraphNeg = new TGraphAsymmErrors(NBinsPt, BinsXCenter, LargestDiffRelNeg, BinsXWidth, BinsXWidth, LargestDiffRelErrorNeg, LargestDiffRelErrorNeg);
    SystErrGraphNeg->SetName(Form("%sRatio_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()));
    TGraphAsymmErrors* SystErrGraphPos = new TGraphAsymmErrors(NBinsPt, BinsXCenter, LargestDiffRelPos, BinsXWidth, BinsXWidth, LargestDiffRelErrorPos, LargestDiffRelErrorPos);
    SystErrGraphPos->SetName(Form("%sRatio_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()));

    // Write sys-err graph to root output file
    TString Outputname = Form("%s/%sRatio_%s_SystematicErrorCuts.root",outputDirRootFile.Data(),meson.Data(),prefix2.Data());
    TFile* SystematicErrorFile = new TFile(Outputname.Data(),"UPDATE");
        SystErrGraphPos->Write(Form("%sRatio_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
        SystErrGraphNeg->Write(Form("%sRatio_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
    SystematicErrorFile->Write();
    SystematicErrorFile->Close();
// */
}

