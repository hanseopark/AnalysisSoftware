//***********************************************************************************************
//**************************** CompDiffCutsOverview ***********************************************
//***********************************************************************************************
/************************************************************************************************
 ******     provided by Gamma Conversion Group, PWG4,                                         *****
 ******        Ana Marin, marin@physi.uni-heidelberg.de                                        *****
 ******        Friederike Bock, friederike.bock@cern.ch                                        *****
 ***********************************************************************************************/

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

struct SysErrorConversion {
   Double_t value;
   Double_t error;
   //    TString name;
};


void CompareDiffCutsTrigg(  TString CombineCutsName                 = "CombineCuts.dat",
                            TString suffix                          = "gif",
                            TString optionEnergy                    = "",
                            TString cutVariationName                = "",
                            Int_t NumberOfCuts                      = 1,
                            TString optionPeriod                    = "No",
                            TString fileName                        = ""
                        ){

    // Define global arrays
    TString     cutNumber               [50];
    TString     cutNumberAdv            [50];

    // Set common default plot style
    StyleSettingsThesis();
    SetPlotStyle();

    // Set cutvariation-name to "" for no explicit name
    if (cutVariationName.CompareTo("None")==0){
        cutVariationName                                        = "";
    }

    // Define Output Directory
    TString outputDir                                           = Form("CompDiffCuts/%s",optionEnergy.Data());
    if (cutVariationName.CompareTo("None")!=0) outputDir        = Form("CompDiffCuts/%s/%s",optionEnergy.Data(),cutVariationName.Data());
    TString outputDirRootFile                                   = Form("CompDiffCuts/%s",optionEnergy.Data());
    gSystem->Exec("mkdir -p "+outputDir);


    // Set collisions system
    TString collisionSystem     = ReturnFullCollisionsSystem(optionEnergy);
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }

    // Define colors for differnt cuts
    Color_t color[20]                                           = { kBlack, kAzure, kGreen+2,kOrange+2,kRed, kViolet,  kBlue-9, kMagenta+4,
                                                                    kCyan+3, kCyan-10, kCyan, kGreen+4, kGreen-9,
                                                                    kGreen,  kYellow+4, kYellow+3, kSpring+10,
                                                                    kMagenta-8, kGray, kGray+3};
    Double_t minCent = -1;
    Double_t maxCent = -1;
    Double_t minCentPlot = 0;
    Double_t maxCentPlot = 100;
    if (cutVariationName.Contains("Cent0010")){
        minCent = 0;
        maxCent = 10;
        minCentPlot = 0;
        maxCentPlot = 15;
    } else if (cutVariationName.Contains("Cent2050")){
        minCent = 20;
        maxCent = 50;
        minCentPlot = 15;
        maxCentPlot = 55;
    } else if (cutVariationName.Contains("Cent3050")){
        minCent = 30;
        maxCent = 50;
        minCentPlot = 25;
        maxCentPlot = 55;
    }

    // Read cuts from CutSelection file
    ifstream in(CombineCutsName.Data());
    cout<<"Available Cuts:"<<endl;
    string TempCutNumber;
    Int_t Number = 0;
    while(getline(in, TempCutNumber)){
        TString tempCutNumberAdv                                = TempCutNumber;
        cutNumberAdv[Number]                                    = tempCutNumberAdv;
        cutNumber[Number]                                       = tempCutNumberAdv;
        cout<< cutNumber[Number]<<endl;
        Number++;
    }
    cout<<"=========================="<<endl;

    cout << "analysing " << cutVariationName << " cut variations" << endl;
    cout << " " << endl;

    // Define necessary histogram/file/string arrays
    const Int_t ConstNumberOfCuts                               = NumberOfCuts;
    const Int_t MaxNumberOfCuts = 20;
    if(ConstNumberOfCuts > MaxNumberOfCuts){
        cout << "Too many cuts, beware!" << endl;
        return;
    }
    TFile*  inputFile           = new TFile(fileName);
    TDirectory*   dirCut        [ConstNumberOfCuts];
    TH1D*   histoCentCut        [ConstNumberOfCuts];
    TH1D*   histoV0MultCut      [ConstNumberOfCuts];
    TH1D*   histoV0TriggCut     [ConstNumberOfCuts];
    TH1D*   histoCentCutRatio   [ConstNumberOfCuts];
    TH1D*   histoV0MultCutRatio [ConstNumberOfCuts];
    TH1D*   histoV0TriggCutRatio[ConstNumberOfCuts];
    Bool_t  isEmulatedTrigg     [ConstNumberOfCuts];
    Double_t nEventTrig         [ConstNumberOfCuts];
    Double_t nEventTrigInCent   [ConstNumberOfCuts];
    Double_t fracNEventTrigInCent[ConstNumberOfCuts];


    for (Int_t i=0; i< NumberOfCuts; i++){
        isEmulatedTrigg[i]          = kFALSE;
        TObjArray *arr              = cutNumber[i].Tokenize("_");
        TString specialTriggCut     = ((TObjString*)arr->At(2))->GetString();
        Int_t specialTrigg          = ((TString)specialTriggCut(0,1)).Atoi();
        if (specialTrigg == 1)
            isEmulatedTrigg[i] = kTRUE;
        cutNumberAdv[i]             = Form("%i%s-%i%s", ((TString)specialTriggCut(1,2)).Atoi(), "%", ((TString)specialTriggCut(3,3)).Atoi(), "%");
        dirCut[i]                   = (TDirectory*)inputFile->Get(Form("Trigger_%s", cutNumber[i].Data()));
        histoCentCut[i]             = (TH1D*)dirCut[i]->Get("histoCent");
        histoCentCut[i]->SetName(Form("histoCent_%s", cutNumber[i].Data()));
        histoV0MultCut[i]           = (TH1D*)dirCut[i]->Get("histoV0Mult");
        histoV0MultCut[i]->SetName(Form("histoV0Mult_%s", cutNumber[i].Data()));
        histoV0TriggCut[i]           = (TH1D*)dirCut[i]->Get("histoV0Trigg");
        histoV0TriggCut[i]->SetName(Form("histoV0Trigg_%s", cutNumber[i].Data()));
        histoCentCutRatio[i]         = (TH1D*)histoCentCut[i]->Clone(Form("histoCentRatio_%s", cutNumber[i].Data()));
        histoCentCutRatio[i]->Sumw2();
        histoCentCutRatio[i]->Divide(histoCentCutRatio[i],histoCentCut[0],1.,1.,"B");
        histoV0MultCutRatio[i]         = (TH1D*)histoV0MultCut[i]->Clone(Form("histoV0MultCutRatio_%s", cutNumber[i].Data()));
        histoV0MultCutRatio[i]->Sumw2();
        histoV0MultCutRatio[i]->Divide(histoV0MultCutRatio[i],histoV0MultCut[0],1.,1.,"B");
        histoV0TriggCutRatio[i]         = (TH1D*)histoV0TriggCut[i]->Clone(Form("histoV0TriggCutRatio_%s", cutNumber[i].Data()));
        histoV0TriggCutRatio[i]->Sumw2();
        histoV0TriggCutRatio[i]->Divide(histoV0TriggCutRatio[i],histoV0TriggCut[0],1.,1.,"B");
        if (isEmulatedTrigg[i]){
            nEventTrig[i]           = histoCentCut[i]->Integral(0,100);
            nEventTrigInCent[i]     = histoCentCut[i]->Integral(minCent,maxCent);
            fracNEventTrigInCent[i] = nEventTrigInCent[i]/nEventTrig[i];
            cutNumberAdv[i]         = Form("emulated %s, %0.2f%s outside des. cent window", cutNumberAdv[i].Data(), (1-fracNEventTrigInCent[i])*100, "%");
        }
    }

    cout<<"=========================="<<endl;


    //**************************************************************************************
    //********************* Plotting Cent Dist *********************************************
    //**************************************************************************************

        TCanvas* canvasCent = new TCanvas("canvasCent","",1350,1500);
        DrawGammaCanvasSettings( canvasCent,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padCent = new TPad("padCent", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padCent, 0.12, 0.02, 0.03, 0.);
//         padCent->SetLogy();
        padCent->Draw();
        // lower pad definition
        TPad* padCentRatios = new TPad("padCentRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padCentRatios, 0.12, 0.02, 0.0, 0.2);
        padCentRatios->Draw();

        padCent->cd();
        padCent->SetTickx();
        padCent->SetTicky();

        // Plot raw yield in uppper panel
        padCent->cd();
        TLegend* legendRawMeson = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
                Double_t scaleFactorRaw = 1.1;
                DrawGammaSetMarker(histoCentCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoCentCut[i],
                                          "", "Cent (%)", "events",
                                        kTRUE, scaleFactorRaw, 1, kFALSE,
                                        kFALSE, 0.0, 0.030,
                                        kTRUE, minCentPlot, maxCentPlot);
                legendRawMeson->AddEntry(histoCentCut[i],Form("standard: %s",cutNumberAdv[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoCentCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoCentCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoCentCut[i]->DrawCopy("same,e1,p");
                legendRawMeson->AddEntry(histoCentCut[i],cutNumberAdv[i].Data());
            }

        }
        legendRawMeson->Draw();
        // Labeling of plot
//         PutProcessLabelAndEnergyOnPlot( 0.94, 0.95, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod, 1, 1.25, 31);

        padCentRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.;
                Double_t maxYRatio = 1.1; //qui
                SetStyleHistoTH1ForGraphs(histoCentCutRatio[i], "Cent (%)", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoCentCutRatio[i], 20, 1.,color[0],color[0]);
                histoCentCutRatio[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoCentCutRatio[i]->GetXaxis()->SetRangeUser(minCentPlot,maxCentPlot);
                for(Int_t b = 0; b< histoCentCutRatio[i]->GetNbinsX(); b++){
                    histoCentCutRatio[i]->SetBinError(b+1,histoCentCut[i]->GetBinError(b+1)/histoCentCut[i]->GetBinContent(b+1));
                }
                histoCentCutRatio[i]->SetFillColor(kGray+2);
                histoCentCutRatio[i]->SetFillStyle(0);
                histoCentCutRatio[i]->DrawCopy("p,e2");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoCentCutRatio[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoCentCutRatio[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoCentCutRatio[i]->DrawCopy("same,e1,p");
            }
            DrawGammaLines(minCentPlot, maxCentPlot, 1., 1.,0.1);
        }

        canvasCent->Update();
        canvasCent->SaveAs(Form("%s/Cent.%s",outputDir.Data(),suffix.Data()));
        delete canvasCent;

    return;
}
