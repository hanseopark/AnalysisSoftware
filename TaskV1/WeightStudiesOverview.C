//***********************************************************************************************
//**************************** WeightStudiesOverview ***********************************************
//***********************************************************************************************
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
#include "./CommonHeaders/PlottingGammaConversionHistos.h"
#include "./CommonHeaders/PlottingGammaConversionAdditional.h"
#include "./CommonHeaders/FittingGammaConversion.h"
#include "./CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "./CommonHeaders/ConversionFunctions.h"

struct SysErrorConversion {
   Double_t value;
   Double_t error;
   //    TString name;
};


void WeightStudiesOverview(TString CombineFilesName             = "CombineCuts.dat",
                        TString suffix                          = "gif",
                        TString optionEnergy                    = "",
                        TString cutVariationName                = "",
                        Int_t NumberOfCuts                      = 1,
                        Int_t mode                              = 9,
                        Bool_t doBarlow                         = kFALSE
                       ){

  // order of cuts should be
  // 0: Onfly Pythia   (default)
  // 1: Onfly Phojet
  // 2: Offline Pythia
  // 3: Offline Phojet

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    StyleSettingsThesis();
    SetPlotStyle();


    // Define global arrays
    TString     cutNumber               [50];
    TString     cutNumberAdv            [50];

    // Set common default plot style
    StyleSettingsThesis();
    SetPlotStyle();

    // Set cutvariation-name to "" for no explicit name
    if (cutVariationName.CompareTo("None")==0) cutVariationName = "";

    // Define Output Directory
    TString outputDir                                           = Form("WeightStudies/%s",optionEnergy.Data());
    if (cutVariationName.CompareTo("None")!=0) outputDir        = Form("WeightStudies/%s/%s",optionEnergy.Data(),cutVariationName.Data());
    TString outputDirRootFile                                   = Form("WeightStudies/%s",optionEnergy.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    Double_t maxR = 180.;
    // Define colors for differnt cuts
    Color_t color[20]                                           = { kBlack, kAzure, kGreen+2,kOrange+2,kRed, kViolet,  kBlue-9, kMagenta+4,
                                                                    kCyan+3, kCyan-10, kCyan, kGreen+4, kGreen-9,
                                                                    kGreen,  kYellow+4, kYellow+3, kSpring+10,
                                                                    kMagenta-8, kGray, kGray+3};

    // Set collisions system
    TString collisionSystem     = ReturnFullCollisionsSystem(optionEnergy);
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    TString detectionProcess    = ReturnFullTextReconstructionProcess(mode);

    // Define necessary histogram/file/string arrays
    const Int_t ConstNumberOfCuts = NumberOfCuts;
    const Int_t MaxNumberOfCuts = 20;
    if(ConstNumberOfCuts > MaxNumberOfCuts){
        cout << "Too many cuts, beware!" << endl;
        return;
    }
    TString FileNames[MaxNumberOfCuts];

    // Read cuts from CutSelection file
    ifstream in(CombineFilesName.Data());
    cout<<"Available Cuts:"<<endl;
    string TempFile;
    Int_t Number = 0;
    while(getline(in, TempFile)){
        TString tempFile                                = TempFile;
        //cout<< tempFile<<endl;
        FileNames[Number]=tempFile;
        cutNumber[Number]=tempFile;
        cutNumber[Number].Remove(0,32);
        cout<< cutNumber[Number].Length()-5<< endl;
        cutNumber[Number].Remove(cutNumber[Number].Length()-5,5);
        cout<<	FileNames[Number]<< "  " << cutNumber[Number].Data() << endl;
        Number++;
    }

    TString periodName[50];
    TString clusterName[50];
    TString generatorName[50];
    TString cutString[50];
    TString V0ReaderName[50];
    for(Int_t i=0; i<Number; i++){
        TObjArray *arr;
        arr = cutNumber[i].Tokenize("_");
        TObjString* string1;
        TObjString* string2;
        TObjString* string3;
        TObjString* string4;
        string1 = (TObjString*)arr->At(0);
        periodName[i] = string1->GetString();
        string2 = (TObjString*)arr->At(1);
        clusterName[i] = string2->GetString();
        if(clusterName[i].Contains("fast") || clusterName[i].Contains("wSDD") || clusterName[i].Contains("woSDD")){
            string3 = (TObjString*)arr->At(2);
            string4 = (TObjString*)arr->At(4);
        } else {
            clusterName[i] = "";
            string3 = (TObjString*)arr->At(1);
            string4 = (TObjString*)arr->At(3);
        }
        generatorName[i] = string3->GetString();
        cutString[i] = string4->GetString();
        TString fV0Reader                                   = cutString[i](GetPhotonV0FinderCutPosition(cutString[i]),1);
        V0ReaderName[i]                                   = AnalyseV0ReaderCut(CutNumberToInteger(fV0Reader));

//         cout << "period " << periodName[i] << " " << clusterName[i] << endl;
//         cout << "generator " << generatorName[i] << endl;
//         cout << "V0Reader " << V0ReaderName[i] << endl;
    }

    cout<<"=========================="<<endl;
    cout << "analysing " << cutVariationName << " cut variations" << endl;
    cout << " " << endl;
    TFile*  WeightFile                 [ConstNumberOfCuts];
    TProfile* profileWeight            [ConstNumberOfCuts];
    TH1F * histWeight                  [ConstNumberOfCuts];
    TH1F * histoRatioWeightCut         [ConstNumberOfCuts];
    TH1F * histoDiffWeightCut          [ConstNumberOfCuts];
    TH1F * histRelUncWeightCut         [ConstNumberOfCuts];
    TH1F * rData [ConstNumberOfCuts];

    for (Int_t i=0; i< NumberOfCuts; i++){

        WeightFile[i] = new TFile(FileNames[i]);
        profileWeight[i]=  (TProfile*)WeightFile[i]->Get("profileContainingMaterialBudgetWeights_manyRadialBins");
//         profileWeight[i]->Print();
        histWeight[i]   = (TH1F*)WeightFile[i]->Get("histoDataMCRatioRinPtBinScaledToGas03");
        rData[i]        = (TH1F*)WeightFile[i]->Get("Data");
//         histWeight[i] = profileWeight[i]->ProjectionX();

        // Calculate ratios
        //      histoRatioWeightCut[i]         = (TH1D*) histWeight[i]->Clone(Form("histoRatioWeights_%d", i));
        histoRatioWeightCut[i] = (TH1F*) histWeight[i]->Clone("histoRatioWeights");
        if (i != 3){
            histoRatioWeightCut[i]->Sumw2();
            histoRatioWeightCut[i]->Divide(histWeight[i],histWeight[0],1.,1.,"");

            //      histoDiffWeightCut[i]         = (TH1D*) histWeight[i]->Clone(Form("histoDiffWeights%d",i));
            histoDiffWeightCut[i] = (TH1F*)histWeight[i]->Clone("histoDiffWeights");
            histoDiffWeightCut[i]->Sumw2();
            histoDiffWeightCut[i]->Add(histWeight[0],-1.);
        } else {
            histoRatioWeightCut[i]->Sumw2();
            histoRatioWeightCut[i]->Divide(histWeight[i],histWeight[2],1.,1.,"");

            //      histoDiffWeightCut[i]         = (TH1D*) histWeight[i]->Clone(Form("histoDiffWeights%d",i));
            histoDiffWeightCut[i] = (TH1F*)histWeight[i]->Clone("histoDiffWeights");
            histoDiffWeightCut[i]->Sumw2();
            histoDiffWeightCut[i]->Add(histWeight[2],-1.);
        }


        histRelUncWeightCut[i]                               = (TH1F*)histWeight[i]->Clone("histRelUncWeight");
        for (Int_t k = 1; k < histRelUncWeightCut[i]->GetNbinsX()+1; k++){
            if (histRelUncWeightCut[i]->GetBinContent(k) != 0){
                histRelUncWeightCut[i]->SetBinContent(k, histRelUncWeightCut[i]->GetBinError(k)/histRelUncWeightCut[i]->GetBinContent(k)*100);
                histRelUncWeightCut[i]->SetBinError(k, 0);
            } else {
                histRelUncWeightCut[i]->SetBinContent(k, 0);
                histRelUncWeightCut[i]->SetBinError(k, 0);
            }
        }

    }
    cout << "line " << __LINE__ << endl;

    // Calculation of total systematics per bin
    Double_t totalWeightError[histoDiffWeightCut[0]->GetNbinsX()];
    Double_t totalWeightSys[histoDiffWeightCut[0]->GetNbinsX()];
    Double_t weightSta[histoDiffWeightCut[0]->GetNbinsX()];
    Double_t weightValue[histoDiffWeightCut[0]->GetNbinsX()];
    Double_t xError[histoDiffWeightCut[0]->GetNbinsX()];

    for (Int_t j=1; j<histoDiffWeightCut[0]->GetNbinsX()+1; j++){
        weightValue[j-1] = histWeight[0]->GetBinContent(j);
        weightSta[j-1] = histWeight[0]->GetBinError(j);
        if (histoDiffWeightCut[1]->GetBinCenter(j) < 5){
            totalWeightSys[j-1] = TMath::Sqrt(histoDiffWeightCut[1]->GetBinContent(j)* histoDiffWeightCut[1]->GetBinContent(j));
        } else if(histoDiffWeightCut[1]->GetBinCenter(j) > 5 &&  histoDiffWeightCut[1]->GetBinCenter(j) < 13) {
            totalWeightSys[j-1] = TMath::Sqrt(1.*histoDiffWeightCut[1]->GetBinContent(j)* histoDiffWeightCut[1]->GetBinContent(j) +
                                              1.*histoDiffWeightCut[2]->GetBinContent(j)* histoDiffWeightCut[2]->GetBinContent(j));
        } else {
        //					histoDiffWeightCut[2]->GetBinContent(j)* histoDiffWeightCut[2]->GetBinContent(j));
            totalWeightSys[j-1] = TMath::Sqrt(histoDiffWeightCut[1]->GetBinContent(j)* histoDiffWeightCut[1]->GetBinContent(j) +
                                              histoDiffWeightCut[2]->GetBinContent(j)* histoDiffWeightCut[2]->GetBinContent(j));
        }
        totalWeightError[j-1] = TMath::Sqrt( weightSta[j-1]* weightSta[j-1]+ totalWeightSys[j-1] *totalWeightSys[j-1] );

        cout << "R = " << histoDiffWeightCut[1]->GetBinCenter(j) << "cm, Weight = " << weightValue[j-1] << ", Sys = " << totalWeightSys[j-1] << " (Sys = " << (totalWeightSys[j-1]/histWeight[0]->GetBinContent(j))*100 << " %),  Stat = " <<  weightSta[j-1] << " (Stat = "<<  (weightSta[j-1]/histWeight[0]->GetBinContent(j))*100 << "%), Tot = "  <<
        totalWeightError[j-1] << " (Tot = "<< (totalWeightError[j-1]/histWeight[0]->GetBinContent(j))*100 << "%)" << endl ;
    }

    // Calculation of average sys. for 1 photon based on the R weight beyond 5 cm
    Int_t       nBins = histWeight[0]->GetNbinsX();
    Double_t*   binsR = new Double_t[nBins+1];
    for (Int_t i=0; i<nBins+1; i++) {
        if (i<nBins)    binsR[i]   = histWeight[0]->GetXaxis()->GetBinLowEdge(i+1);
        else            binsR[i]   = histWeight[0]->GetXaxis()->GetBinUpEdge(i);
    }

    TH1F * histWeightStaErr = new TH1F("weightStaErr","weightStaErr",nBins,binsR);
    TH1F * histWeightSysErr = new TH1F("weightSysErr","weightSysErr",nBins,binsR);
    TH1F * histWeightTotErr = new TH1F("weightTotErr","weightTotErr",nBins,binsR);
    TProfile* fProfileContainingMaterialBudgetWeights = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins","profileContainingMaterialBudgetWeights_manyRadialBins",nBins,binsR);
    TProfile* fProfileContainingMaterialBudgetWeightsUpTotErr = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBinsUpTotErr","profileContainingMaterialBudgetWeights_manyRadialBinsUpTotErr",nBins,binsR);
    TProfile* fProfileContainingMaterialBudgetWeightsDownTotErr = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBinsDownTotErr","profileContainingMaterialBudgetWeights_manyRadialBinsDownTotErr",nBins,binsR);

    for (Int_t j=1; j<histoDiffWeightCut[0]->GetNbinsX()+1; j++){
        histWeightStaErr->SetBinContent(j,weightValue[j-1]);
        histWeightStaErr->SetBinError(j,weightSta[j-1]);

        histWeightSysErr->SetBinContent(j,weightValue[j-1]);
        histWeightSysErr->SetBinError(j,totalWeightSys[j-1]);

        histWeightTotErr->SetBinContent(j,weightValue[j-1]);
        histWeightTotErr->SetBinError(j,totalWeightError[j-1]);
        fProfileContainingMaterialBudgetWeights->Fill(binsR[j-1], weightValue[j-1]);
        fProfileContainingMaterialBudgetWeightsUpTotErr->Fill(binsR[j-1], weightValue[j-1]+totalWeightError[j-1]);
        fProfileContainingMaterialBudgetWeightsDownTotErr->Fill(binsR[j-1], weightValue[j-1]-totalWeightError[j-1]);
    }

    // cout<< "Bins::"<< rData[0]->GetXaxis()->FindBin(5.+0.001)<< endl;
    // cout<< "Bins::"<< rData[0]->GetXaxis()->FindBin(180.-0.001)<< endl;
    // cout << "Bin width::"<<     rData[0]->GetBinWidth(1)<< " " <<rData[0]->GetBinWidth(10)<< " " << endl;
    Double_t totInRRange= rData[0]->Integral(rData[0]->GetXaxis()->FindBin(5.+0.001),rData[0]->GetXaxis()->FindBin(180.-0.001),"width");
    Double_t averageWeightSys=0.;
    Double_t averageWeightSta=0.;
    Double_t averageWeightTot=0.;
    Double_t averageWeight=0.;
    Int_t rminBin=0;
    Int_t rmaxBin=0;

    for (Int_t j=2; j<nBins; j++){

        rminBin=rData[0]->GetXaxis()->FindBin(binsR[j]+0.0001);
        rmaxBin=rData[0]->GetXaxis()->FindBin(binsR[j+1]-0.0001);
        Double_t totInRBin= rData[0]->Integral(rminBin,rmaxBin,"width");
        //     cout<<"R Limits"<< binsR[j] << "  " << binsR[j+1]<< " " << rminBin<<"  " << rmaxBin<< " " << totInRBin  << "  " << totInRRange<< "  " << histWeight[0]->GetBinContent(j+1)<< endl;
        //      cout<<"R Limits::"<<  binsR[j] << "  " << binsR[j+1]<< " " << totInRBin  << "  " << totInRBin/totInRRange<< "  "<<  histWeight[0]->GetBinContent(j+1)<< endl;
        averageWeightSta+=weightSta[j]*(totInRBin/totInRRange);
        averageWeightSys+=totalWeightSys[j]*(totInRBin/totInRRange);
        averageWeightTot+=totalWeightError[j]*(totInRBin/totInRRange);
        averageWeight+=weightValue[j] *(totInRBin/totInRRange);
    }
    cout << "Average Weight = " << averageWeight << ", Sys = " << (averageWeightSys/averageWeight)*100 << "%, Stat = " << (averageWeightSta/averageWeight)*100 << "%, Tot = " << (averageWeightTot/averageWeight)*100 << "%" << endl;


    // plotting
    Double_t minYRatio = 0.9;
    Double_t maxYRatio = 1.5; //qui
    TCanvas* canvasWeight = new TCanvas("canvasWeight","",1350,1500);
    DrawGammaCanvasSettings( canvasWeight,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padWeight = new TPad("padWeight", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padWeight, 0.12, 0.03, 0.02, 0.);
        //  padWeight->SetLogy();
        padWeight->Draw();
        // lower pad definition
        TPad* padWeightRatios = new TPad("padWeightRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padWeightRatios, 0.12, 0.03, 0.0, 0.2);
        padWeightRatios->Draw();
        padWeight->cd();
        padWeight->SetTickx();
        padWeight->SetTicky();
        // Plot raw yield in uppper panel
        padWeight->cd();       // Set legend

        TH2F * histo2DDummy = new TH2F("","",1000,0.,180.,1000,0.,1.5);
        SetStyleHistoTH2ForGraphs(histo2DDummy, "R (cm)", "Weights",0.04,0.04, 0.04,0.04, 1.,1.);
//         histo2DDummy->GetXaxis()->SetMoreLogLabels();
//         histo2DDummy->GetXaxis()->SetLabelOffset(-0.01);
        histo2DDummy->GetYaxis()->SetRangeUser(0.81,1.39);
        histo2DDummy->Draw("copy");

        TLegend* legendWeight = GetAndSetLegend2(0.3,0.93-1.15*0.04*NumberOfCuts,0.65,0.93,36);
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
//                 SetStyleHistoTH1ForGraphs(histWeight[i], "R (cm)", "Weight", 0.04, 0.06, 0.03, 0.05, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histWeight[i], 20, 1.5,color[i],color[i]);

                histWeight[i]->Draw("same,l");

                legendWeight->AddEntry(histWeight[i],Form("%s, %s (%s) - %s",periodName[i].Data(),clusterName[i].Data(),generatorName[i].Data(),V0ReaderName[i].Data()));   //" %s",cutNumber[i].Data()));
            } else {
//                 SetStyleHistoTH1ForGraphs(histWeight[i], "R (cm)", "Weight", 0.04, 0.06, 0.03, 0.05, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histWeight[i], 20, 1.5,color[i],color[i]);
                legendWeight->AddEntry(histWeight[i],Form("%s, %s (%s) - %s",periodName[i].Data(),clusterName[i].Data(),generatorName[i].Data(),V0ReaderName[i].Data()));   //" %s",cutNumber[i].Data()));
                histWeight[i]->Draw("same,l");
            }
            DrawGammaLines(0., maxR,1., 1.,1.1,kBlack,2);
            legendWeight->Draw();
        }
        padWeightRatios->cd();
        minYRatio = 0.935;
        maxYRatio = 1.05; //qui
        TH2F * histo2DDummyRatio = new TH2F("","",1000,0.,180.,1000,0.,2.);
        SetStyleHistoTH2ForGraphs(histo2DDummyRatio, "R (cm)", "#frac{modified}{standard}",0.08,0.085, 0.08,0.085,0.9,0.67);
        histo2DDummyRatio->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
        histo2DDummyRatio->Draw("copy");

        DrawGammaLines(0., maxR,1., 1.,1.1,kBlack,2);
        DrawGammaLines(0., maxR,1.02, 1.02,1.1,kBlack,2);
        DrawGammaLines(0., maxR,0.98, 0.98,1.1,kBlack,2);

        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i<20)  DrawGammaSetMarker(histoRatioWeightCut[i], 20, 1.5,color[i],color[i]);
            else      DrawGammaSetMarker(histoRatioWeightCut[i], 20, 1.5,color[i-20],color[i-20]);
            histoRatioWeightCut[i]->Draw("same");
        }

    canvasWeight->Update();
    canvasWeight->SaveAs(Form("%s/Weight_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));
    delete canvasWeight;


    TCanvas* canvasRelUncWeight = new TCanvas("canvasRelUncWeight","",1300,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelUncWeight,  0.08, 0.03, 0.02, 0.08);
    canvasRelUncWeight->SetGridy();

        Double_t maxYRel    = 0;
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if (histRelUncWeightCut[i]->GetMaximum() > maxYRel)
                maxYRel = histRelUncWeightCut[i]->GetMaximum();
        }
        if (maxYRel > 60) maxYRel = 60;
        TH2F* histo1DDummy2              = new TH2F("histo1DDummy2","histo1DDummy2",1000, 0., 180.,1000,0,100.);
        SetStyleHistoTH1ForGraphs(histo1DDummy2,  "R (cm)", Form("Weight rel unc. (%s)", "%"), 0.035 ,0.04, 0.035,0.04, 0.9, 1.,510,540);
        histo1DDummy2->GetYaxis()->SetRangeUser(0,maxYRel*1.05);
        histo1DDummy2->DrawCopy();

        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i<20) DrawGammaSetMarker(histRelUncWeightCut[i], 20, 1.5,color[i],color[i]);
            else     DrawGammaSetMarker(histRelUncWeightCut[i], 20, 1.5,color[i-20],color[i-20]);
            histRelUncWeightCut[i]->DrawCopy("same,c,p");
        }
        legendWeight->Draw();

    canvasRelUncWeight->Update();
    canvasRelUncWeight->SaveAs(Form("%s/WeightRelUncertainty_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));
    delete canvasRelUncWeight;


    TCanvas* canvasWeightDiff = new TCanvas("canvasWeightDiff","",1300,900);
    DrawGammaCanvasSettings( canvasWeightDiff, 0.08, 0.03, 0.035, 0.09);

        Double_t minYDiff = -0.15;
        Double_t maxYDiff = 0.15; //qui
        TH2F * histo2DDummyDiff = new TH2F("","",1000,0.,180.,1000,-1.5,1.);
        SetStyleHistoTH2ForGraphs(histo2DDummyDiff, "R (cm)", "modified - standard",0.04,0.04, 0.04,0.04, 1.,1.);
        histo2DDummyDiff->GetYaxis()->SetRangeUser(minYDiff,maxYDiff);
        histo2DDummyDiff->Draw("copy");

        DrawGammaLines(0., maxR,1., 1.,1.1,kBlack,2);
        DrawGammaLines(0., maxR,0.02, 0.02,1.1,kBlack,2);
        DrawGammaLines(0., maxR,-0.02, -0.02,1.1,kBlack,2);

        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i<20) DrawGammaSetMarker(histoDiffWeightCut[i], 20, 1.5,color[i],color[i]);
            else     DrawGammaSetMarker(histoDiffWeightCut[i], 20, 1.5,color[i-20],color[i-20]);
            histoDiffWeightCut[i]->Draw("same");
        }
        legendWeight->Draw();

    canvasWeightDiff->Update();
    canvasWeightDiff->SaveAs(Form("%s/WeightDiff_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));
    delete canvasWeightDiff;

    TFile outFile(Form("%s/weightsWithErrors_Test.root",outputDir.Data()) ,"RECREATE");
    histWeightStaErr->Write();
    histWeightSysErr->Write();
    histWeightTotErr->Write();
    fProfileContainingMaterialBudgetWeights->Write();
    fProfileContainingMaterialBudgetWeightsUpTotErr->Write();
    fProfileContainingMaterialBudgetWeightsDownTotErr->Write();
    outFile.Close();

 }
