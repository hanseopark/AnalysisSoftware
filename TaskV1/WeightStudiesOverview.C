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


void WeightStudiesOverview(TString CombineFilesName             = "CombineCuts.dat",
                        TString suffix                          = "gif",
                        TString optionEnergy                    = "",
                        TString cutVariationName                = "",
                        Int_t NumberOfCuts                      = 1,
                        Int_t mode                              = 9,
                        Int_t sequence                          = 0,
                        Bool_t doBarlow                         = kFALSE
                       ){

  // order of cuts should be if sequence =
  // 0: Onfly Pythia (default), Onfly Phojet, Offline Pythia, Offline Phojet
  // 1: Onfly Pythia 13TeV, Onfly Phojet 13TeV, Onfly Pythia 5TeV, Onfly Phojet 5TeV, idem with Offline
  // 2: Onfly Pythia 5TeV2017, Offline Pythia 5TeV2017, Onfly Pythia 5TeV2015, Offline Pythia 5TeV2015

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
    TString outputDir                                    = Form("WeightStudies/%s",optionEnergy.Data());
    if (cutVariationName.CompareTo("None")!=0) outputDir = Form("WeightStudies/%s/%s",optionEnergy.Data(),cutVariationName.Data());
    TString outputDirRootFile                            = Form("WeightStudies/%s",optionEnergy.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    Double_t maxR = 180.;
    // Define colors for differnt cuts
    Color_t color[12] = { kBlack, kAzure, kGreen+2, kOrange+2,
                          kRed+1, kCyan+2, kYellow+2, kViolet-3,
                          kSpring+10, kMagenta-8, kGray, kGray+3};

    TString label5cm = "#rightarrow R = 5 cm";

    // Set collisions system
    TString collisionSystem     = ReturnFullCollisionsSystem(optionEnergy);
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    TString detectionProcess    = ReturnFullTextReconstructionProcess(mode);

    // Define necessary histogram/file/string arrays
    const Int_t ConstNumberOfCuts = NumberOfCuts;
    const Int_t MaxNumberOfCuts = 12;
    const Int_t nBinsR = 12;
    if(ConstNumberOfCuts > MaxNumberOfCuts){
        cout << "Too many cuts, beware!" << endl;
        return;
    }
    TString FileNames[MaxNumberOfCuts];
    TString FileNamesAdditional[MaxNumberOfCuts];

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
        FileNamesAdditional[Number] = Form("AdditionalMBHistos%s",cutNumber[Number].Data());
        cout<< cutNumber[Number].Length()-5<< endl;
        cutNumber[Number].Remove(cutNumber[Number].Length()-5,5);
        cout<<	FileNames[Number]<< "  " << cutNumber[Number].Data() << endl;
	cout<<  FileNamesAdditional[Number] << endl;
        Number++;
    }

    TString period[50];
    TString periodName[50];
    TString clusterName[50];
    TString generatorName[50];
    TString cutString[50];
    TString V0ReaderName[50];
    Marker_t marker[50];
    Int_t counterOnfly = 0;
    Int_t counterOffline = 0;
    for(Int_t i=0; i<Number; i++){
        TObjArray *arr;
        arr = cutNumber[i].Tokenize("_");
        TObjString* string1;
        TObjString* string2;
        TObjString* string3;
        TObjString* string4;
        string1 = (TObjString*)arr->At(0);
        period[i] = string1->GetString();
        if(period[i].Contains("LHC16")) periodName[i] = "13 TeV (Ar)";
        else if(period[i].Contains("LHC17")) periodName[i] = "  5 TeV (Ne)";
        else periodName[i] = "  5 TeV (Ar)";
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
        if(V0ReaderName[i].Contains("Onfly")){
            V0ReaderName[i] = "On-the-Fly V0";
            marker[i] = 20;
            counterOnfly++;
        } else if (V0ReaderName[i].Contains("Offline")){
            V0ReaderName[i] = "Offline V0";
            marker[i] = 21;
            counterOffline++;
        }

	//	cout << "period " << periodName[i] << " " << clusterName[i] << endl;
	//        cout << "generator " << generatorName[i] << endl;
	//	cout << "V0Reader " << V0ReaderName[i] << endl;
    }

    cout<<"=========================="<<endl;
    cout << "analysing " << cutVariationName << " cut variations" << endl;
    cout << " " << endl;
    TFile*  WeightFile                 [ConstNumberOfCuts];
    TProfile* profileWeight            [ConstNumberOfCuts];
    TH1F * histWeight                  [ConstNumberOfCuts];
    TH1F * histoRatioWeightCut         [ConstNumberOfCuts];
    TH1F * histoRatioWeightEnergyCut   [ConstNumberOfCuts];
    TH1F * histoDiffWeightCut          [ConstNumberOfCuts];
    TH1F * histRelUncWeightCut         [ConstNumberOfCuts];
    TH1F * rData [ConstNumberOfCuts];

    TFile*  WeightFileAdditional       [ConstNumberOfCuts];
    TH1F * histoWeightsEachRPtMin[nBinsR][ConstNumberOfCuts];
    TH1F * histoWeightsEachRPtMinSecSub[nBinsR][ConstNumberOfCuts];   // Still a test
    TH1F * histoWeightsEachRPtMinSecSubUsingCocktail[nBinsR][ConstNumberOfCuts];   // Still a test
    TH1F * histoPurityPtEachRBin[nBinsR][ConstNumberOfCuts];

    for (Int_t i=0; i< NumberOfCuts; i++){

        WeightFile[i] = TFile::Open(FileNames[i]);
        profileWeight[i]=  (TProfile*)WeightFile[i]->Get("profileContainingMaterialBudgetWeights_manyRadialBins");
	//histWeight[i]   = (TH1F*)WeightFile[i]->Get("histoDataMCRatioRinPtBinScaledToGasPtBin3");
	//AM: Take the scondary subtracted histogram
        histWeight[i]   = (TH1F*)WeightFile[i]->Get("histoDataMCRatioRScaledToGasSecSubPtBin3");
        rData[i]        = (TH1F*)WeightFile[i]->Get("Data");

	WeightFileAdditional[i] = TFile::Open(FileNamesAdditional[i]);

        // Calculate ratios
        histoRatioWeightCut[i] = (TH1F*) histWeight[i]->Clone("histoRatioWeights");
        histoRatioWeightCut[i]->Sumw2();
        histoDiffWeightCut[i] = (TH1F*)histWeight[i]->Clone("histoDiffWeights");
        histoDiffWeightCut[i]->Sumw2();
        histoRatioWeightEnergyCut[i] = (TH1F*) histWeight[i]->Clone("histoRatioWeightsByEnergy");
        histoRatioWeightEnergyCut[i]->Sumw2();
	for (Int_t j=0; j< nBinsR; j++){
	  //	  cout<< "name:"<< Form("histoWeightsEachRPtMin%i",j)<< endl;
	  histoWeightsEachRPtMin[j][i] = (TH1F*)WeightFileAdditional[i]->Get(Form("histoWeightsEachRPtMin%i",j));
	  histoWeightsEachRPtMinSecSub[j][i] = (TH1F*)WeightFileAdditional[i]->Get(Form("histoWeightsEachRPtMinSecSub%i",j));
	  histoWeightsEachRPtMinSecSubUsingCocktail[j][i] = (TH1F*)WeightFileAdditional[i]->Get(Form("histoWeightsEachRPtMinSecSubUsingCocktail%i",j));
	  histoPurityPtEachRBin[j][i] = (TH1F*)WeightFileAdditional[i]->Get(Form("histoPurityPtEachRBin_%i",j)); ;
	  //	  cout<< "mean" <<histoWeightsEachRPtMin[i][j]->GetMean() <<endl;
	}
        if(sequence==0){
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

        } else if(sequence==1){
            if ( i <= 4 || V0ReaderName[i].Contains("On-the-Fly") ){
                histoRatioWeightCut[i]->Divide(histWeight[i],histWeight[0],1.,1.,"");
                histoDiffWeightCut[i]->Add(histWeight[0],-1.);
            } else if( i > 4 && V0ReaderName[i].Contains("Offline") ){
                histoRatioWeightCut[i]->Divide(histWeight[i],histWeight[4],1.,1.,"");
                histoDiffWeightCut[i]->Add(histWeight[4],-1.);
            }
            if(period[i].Contains("LHC16")){
                cout << i <<  " 13 TeV" << endl;
                if ( i != 5 || V0ReaderName[i].Contains("On-the-Fly") )
                    histoRatioWeightEnergyCut[i]->Divide(histWeight[i],histWeight[0],1.,1.,"");
                else if( i == 5 && V0ReaderName[i].Contains("Offline") )
                    histoRatioWeightEnergyCut[i]->Divide(histWeight[i],histWeight[4],1.,1.,"");
            } else if(period[i].Contains("LHC17")){
                cout << i <<  " 5 TeV" << endl;
                if ( i != 7 || V0ReaderName[i].Contains("On-the-Fly") )
                    histoRatioWeightEnergyCut[i]->Divide(histWeight[i],histWeight[2],1.,1.,"");
                else if( i == 7 && V0ReaderName[i].Contains("Offline") )
                    histoRatioWeightEnergyCut[i]->Divide(histWeight[i],histWeight[6],1.,1.,"");
            }

        } else if(sequence==2){

            if (V0ReaderName[i].Contains("On-the-Fly") ){
                histoRatioWeightCut[i]->Divide(histWeight[i],histWeight[0],1.,1.,"");
                histoDiffWeightCut[i]->Add(histWeight[0],-1.);
            } else if(V0ReaderName[i].Contains("Offline") ){
                histoRatioWeightCut[i]->Divide(histWeight[i],histWeight[4],1.,1.,"");
                histoDiffWeightCut[i]->Add(histWeight[2],-1.);
            }
            if(period[i].Contains("LHC17")){
                cout << i <<  " 5 TeV" << endl;
                if ( i != 7 || V0ReaderName[i].Contains("On-the-Fly") )
                    histoRatioWeightEnergyCut[i]->Divide(histWeight[i],histWeight[2],1.,1.,"");
                else if( i == 7 && V0ReaderName[i].Contains("Offline") )
                    histoRatioWeightEnergyCut[i]->Divide(histWeight[i],histWeight[6],1.,1.,"");
            } else  if(period[i].Contains("LHC15")){
                cout << i <<  " 5 TeV" << endl;
                if ( i != 4 || V0ReaderName[i].Contains("On-the-Fly") )
                    histoRatioWeightEnergyCut[i]->Divide(histWeight[i],histWeight[0],1.,1.,"");
                else if( i == 4 && V0ReaderName[i].Contains("Offline") )
                    histoRatioWeightEnergyCut[i]->Divide(histWeight[i],histWeight[1],1.,1.,"");
            }

        }

        histRelUncWeightCut[i] = (TH1F*)histWeight[i]->Clone("histRelUncWeight");
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

        cout << "R = " << histoDiffWeightCut[1]->GetBinCenter(j) << "cm, Weight = " << weightValue[j-1] << ", Stat = " <<  weightSta[j-1] << " (Stat = "<<  (weightSta[j-1]/histWeight[0]->GetBinContent(j))*100 << "%), Sys = " << totalWeightSys[j-1] << " (Sys = " << (totalWeightSys[j-1]/histWeight[0]->GetBinContent(j))*100 << "%), Tot = "  << totalWeightError[j-1] << " (Tot = "<< (totalWeightError[j-1]/histWeight[0]->GetBinContent(j))*100 << "%)" << endl ;
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
    Double_t totInRRange = rData[0]->Integral(rData[0]->GetXaxis()->FindBin(5.+0.001),rData[0]->GetXaxis()->FindBin(180.-0.001),"width");
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
    cout << "Average Weight = " << averageWeight << ", Stat = " << (averageWeightSta/averageWeight)*100 << "%, Sys = " << (averageWeightSys/averageWeight)*100 << "%, Tot = " << (averageWeightTot/averageWeight)*100 << "%" << endl;


    // plotting
    Double_t minYRatio = 0.9;
    Double_t maxYRatio = 1.5;
    TCanvas* canvasWeight = new TCanvas("canvasWeight","",1500,1300);
    DrawGammaCanvasSettings( canvasWeight,  0.13, 0.02, 0.02, 0.09);
        TPad* padWeight = new TPad("padWeight", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padWeight, 0.1, 0.03, 0.02, 0.);
        padWeight->Draw();

        TPad* padWeightRatios = new TPad("padWeightRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padWeightRatios, 0.1, 0.03, 0.0, 0.2);
        padWeightRatios->Draw();
        padWeight->cd();
        padWeight->SetTickx();
        padWeight->SetTicky();
        padWeight->cd();

        TH2F * histo2DDummy = new TH2F("","",1000,0.,180.,1000,0.,1.5);
        SetStyleHistoTH2ForGraphs(histo2DDummy, "R (cm)", "Weights",0.04,0.04, 0.04,0.04, 1.,1.);
//         histo2DDummy->GetXaxis()->SetMoreLogLabels();
//         histo2DDummy->GetXaxis()->SetLabelOffset(-0.01);
        histo2DDummy->GetYaxis()->SetRangeUser(0.875,1.35);
        histo2DDummy->Draw("copy");

        DrawGammaLines(5., 5.,0.875,1.35,1.1,kGray+2,1);
        TLatex *labelR5cm = new TLatex(0.13,0.10,label5cm.Data());
        SetStyleTLatex( labelR5cm,0.04,4);
        labelR5cm->Draw();

        TLegend* legendWeightOnfly= new TLegend(0.45,0.93-(1+counterOnfly)*0.04,0.7,0.93); //0.17,0.13,0.5,0.24);
        legendWeightOnfly->SetFillColor(0);
        legendWeightOnfly->SetMargin(0.17);
        legendWeightOnfly->SetLineColor(0);
        legendWeightOnfly->SetTextFont(42);
        legendWeightOnfly->SetTextSize(0.035);
        legendWeightOnfly->SetHeader("On-the-Fly V0 finder");

        TLegend* legendWeightOffline= new TLegend(0.7,0.93-(1+counterOffline)*0.04,0.93,0.93); //0.17,0.13,0.5,0.24);
        legendWeightOffline->SetFillColor(0);
        legendWeightOffline->SetMargin(0.17);
        legendWeightOffline->SetLineColor(0);
        legendWeightOffline->SetTextFont(42);
        legendWeightOffline->SetTextSize(0.035);
        legendWeightOffline->SetHeader("Offline V0 finder");

        TLegend* legendWeight = GetAndSetLegend2(0.3,0.93-1.15*0.04*NumberOfCuts,0.65,0.93,40);
        DrawGammaLines(0., maxR,1., 1.,1.1,kGray+1,2);
        for(Int_t i = 0; i< NumberOfCuts; i++){

            DrawGammaSetMarker(histWeight[i], marker[i], 1.5,color[i],color[i]);
            histWeight[i]->Draw("same,c,p");

            if(V0ReaderName[i].Contains("Offline"))
                legendWeightOffline->AddEntry(histWeight[i],Form("%s, %s",periodName[i].Data(),generatorName[i].Data()));
            else
                legendWeightOnfly->AddEntry(histWeight[i],Form("%s, %s",periodName[i].Data(),generatorName[i].Data()));

            legendWeightOnfly->Draw();
            legendWeightOffline->Draw();


                legendWeight->AddEntry(histWeight[i],Form("%s, (%s) - %s",periodName[i].Data(),generatorName[i].Data(),V0ReaderName[i].Data()));   //" %s",cutNumber[i].Data()));
//                 legendWeight->AddEntry(histWeight[i],Form("%s, %s (%s) - %s",periodName[i].Data(),clusterName[i].Data(),generatorName[i].Data(),V0ReaderName[i].Data()));   //" %s",cutNumber[i].Data()));
//             legendWeight->Draw();
        }

        histo2DDummy->DrawCopy("axis,same");
        padWeightRatios->cd();
        minYRatio = 0.92;
        maxYRatio = 1.065; //qui
        TH2F * histo2DDummyRatio = new TH2F("","",1000,0.,180.,1000,0.,2.);
        SetStyleHistoTH2ForGraphs(histo2DDummyRatio, "R (cm)", "#frac{modified}{standard}",0.08,0.085, 0.08,0.085,0.9,0.55);
        histo2DDummyRatio->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
        histo2DDummyRatio->Draw("copy");

        DrawGammaLines(0., maxR,1., 1.,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,1.02, 1.02,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,0.98, 0.98,1.1,kGray+1,2);

        DrawGammaLines(5., 5.,minYRatio,maxYRatio,1.1,kGray+2,1);

        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i<12)  DrawGammaSetMarker(histoRatioWeightCut[i], marker[i], 1.5,color[i],color[i]);
            else      DrawGammaSetMarker(histoRatioWeightCut[i], marker[i], 1.5,color[i-12],color[i-12]);
            histoRatioWeightCut[i]->Draw("same,c,p");
        }

        histo2DDummyRatio->DrawCopy("axis,same");
    canvasWeight->Update();
    canvasWeight->SaveAs(Form("%s/Weight_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));

    canvasWeight->cd();
    padWeight->cd();

        histo2DDummy->Draw("copy");
        DrawGammaLines(0., maxR,1., 1.,1.1,kGray+1,2);
        DrawGammaLines(5., 5.,0.875,1.35,1.1,kGray+2,1);
        labelR5cm->Draw();

        TLegend* legendStd= new TLegend(0.6,0.93-3.5*0.04,0.93,0.93); //0.17,0.13,0.5,0.24);
        legendStd->SetFillColor(0);
        legendStd->SetMargin(0.17);
        legendStd->SetLineColor(0);
        legendStd->SetTextFont(42);
        legendStd->SetTextSize(0.04);
        legendStd->SetHeader("Pythia, On-the-Fly V0 finder");

        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(generatorName[i].Contains("Pythia") && V0ReaderName[i].Contains("On-the-Fly")){

                DrawGammaSetMarker(histWeight[i], marker[i], 1.5,color[i],color[i]);
                histWeight[i]->Draw("same,c,p");

                legendStd->AddEntry(histWeight[i],Form("%s",periodName[i].Data()));
            }
        }
        legendStd->Draw();

        histo2DDummy->DrawCopy("axis,same");
        padWeightRatios->cd();
        histo2DDummyRatio->Draw("copy");

        DrawGammaLines(0., maxR,1., 1.,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,1.02, 1.02,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,0.98, 0.98,1.1,kGray+1,2);
        DrawGammaLines(5., 5.,minYRatio,maxYRatio,1.1,kGray+2,1);
//         TLatex *labelR5cm = new TLatex(0.13,0.9,label5cm.Data());
        SetStyleTLatex( labelR5cm,0.04,4);

        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(generatorName[i].Contains("Pythia") && V0ReaderName[i].Contains("On-the-Fly")){
                DrawGammaSetMarker(histoRatioWeightCut[i], marker[i], 1.5,color[i],color[i]);
                histoRatioWeightCut[i]->Draw("same,c,p");
            }
        }

        histo2DDummyRatio->DrawCopy("axis,same");
    canvasWeight->Update();
    canvasWeight->SaveAs(Form("%s/WeightStandard_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));
    delete canvasWeight;


    TCanvas* canvasWeightDoublePad = new TCanvas("canvasWeightDoublePad","",1500,1300);
    DrawGammaCanvasSettings( canvasWeightDoublePad,  0.13, 0.02, 0.02, 0.03);

    TPad* padWeight2PadRatio1 = new TPad("padWeight2PadRatio1", "", 0., 0.5, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padWeight2PadRatio1, 0.12, 0.03, 0.02, 0.);
    padWeight2PadRatio1->Draw();

    TPad* padWeight2PadRatio2 = new TPad("padWeight2PadRatio2", "", 0., 0., 1., 0.5,-1, -1, -2);
    DrawGammaPadSettings( padWeight2PadRatio2, 0.12, 0.03, 0.0, 0.15);
    padWeight2PadRatio2->Draw();


    padWeight2PadRatio1->cd();
    padWeight2PadRatio1->SetTickx();
    padWeight2PadRatio1->SetTicky();
    padWeight2PadRatio1->cd();

        TH2F * histo2DDummyUpper = new TH2F("","",1000,0.,180.,1000,0.,2);
        SetStyleHistoTH2ForGraphs(histo2DDummyUpper, "R (cm)", "#frac{modified}{standard}",0.06,0.07, 0.06,0.07,0.9,0.8);
        histo2DDummyUpper->GetYaxis()->SetRangeUser(0.93,1.13);
        histo2DDummyUpper->Draw("copy");

        DrawGammaLines(0., maxR,1., 1.,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,1.02, 1.02,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,0.98, 0.98,1.1,kGray+1,2);

        DrawGammaLines(5., 5.,0.93,1.13,1.1,kGray+2,1);
        labelR5cm = new TLatex(0.15,0.88,label5cm.Data());
        SetStyleTLatex( labelR5cm,0.06,4);
        labelR5cm->Draw();

        TLegend* legend13TeV = new TLegend(0.6,0.93-(1+counterOnfly)*0.06,0.93,0.93); //0.17,0.13,0.5,0.24);
        legend13TeV->SetFillColor(0);
        legend13TeV->SetMargin(0.17);
        legend13TeV->SetLineColor(0);
//         legend13TeV->SetNColumns(2);
        legend13TeV->SetTextFont(42);
        legend13TeV->SetTextSize(0.06);
        legend13TeV->SetHeader("pp, #sqrt{s} = 13 TeV");
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(period[i].Contains("LHC16")){
                DrawGammaSetMarker(histoRatioWeightEnergyCut[i], marker[i], 1.5,color[i],color[i]);
                histoRatioWeightEnergyCut[i]->Draw("same,c,p");

                legend13TeV->AddEntry(histoRatioWeightEnergyCut[i],Form("%s, %s",generatorName[i].Data(), V0ReaderName[i].Data()));
            }
        }
        legend13TeV->Draw();

        histo2DDummyUpper->DrawCopy("axis,same");
    padWeight2PadRatio2->cd();

        TH2F * histo2DDummyLower = new TH2F("","",1000,0.,180.,1000,0.,2.);
        SetStyleHistoTH2ForGraphs(histo2DDummyLower, "R (cm)", "#frac{modified}{standard}",0.06,0.07, 0.06,0.07,0.9,0.8);
        histo2DDummyLower->GetYaxis()->SetRangeUser(0.93,1.13);
        histo2DDummyLower->Draw("copy");

        DrawGammaLines(0., maxR,1., 1.,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,1.02, 1.02,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,0.98, 0.98,1.1,kGray+1,2);

        DrawGammaLines(5., 5.,0.93,1.13,1.1,kGray+2,1);
        labelR5cm = new TLatex(0.15,0.88,label5cm.Data());
        SetStyleTLatex( labelR5cm,0.06,4);
        labelR5cm->Draw();

        TLegend* legend5TeV = new TLegend(0.6,0.93-(1+counterOnfly)*0.06,0.93,0.93); //0.17,0.13,0.5,0.24);
        legend5TeV->SetFillColor(0);
        legend5TeV->SetMargin(0.17);
        legend5TeV->SetLineColor(0);
//         legend5TeV->SetNColumns(2);
        legend5TeV->SetTextFont(42);
        legend5TeV->SetTextSize(0.06);
        legend5TeV->SetHeader("pp, #sqrt{s} = 5 TeV");
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(period[i].Contains("LHC17")){
                DrawGammaSetMarker(histoRatioWeightEnergyCut[i], marker[i], 1.5,color[i-2],color[i-2]);
                histoRatioWeightEnergyCut[i]->Draw("same,c,p");

                legend5TeV->AddEntry(histoRatioWeightEnergyCut[i],Form("%s, %s",generatorName[i].Data(), V0ReaderName[i].Data()));
            }
        }
        legend5TeV->Draw();

        histo2DDummyRatio->DrawCopy("axis,same");
    canvasWeightDoublePad->Update();
    canvasWeightDoublePad->SaveAs(Form("%s/Weight2PadRatio_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));
    delete canvasWeightDoublePad;

    TCanvas* canvasRelUncWeight = new TCanvas("canvasRelUncWeight","",1300,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelUncWeight,  0.08, 0.03, 0.02, 0.08);
//     canvasRelUncWeight->SetGridy();

        Double_t maxYRel    = 0;
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if (histRelUncWeightCut[i]->GetMaximum() > maxYRel)
                maxYRel = histRelUncWeightCut[i]->GetMaximum();
        }
        if (maxYRel > 60) maxYRel = 60;
        TH2F* histo1DDummy2              = new TH2F("histo1DDummy2","histo1DDummy2",1000, 0., 180.,1000,0,10.);
        SetStyleHistoTH1ForGraphs(histo1DDummy2,  "R (cm)", Form("Weight rel unc. (%s)", "%"), 0.035 ,0.04, 0.035,0.04, 0.9, .9);//,510,540);
        histo1DDummy2->GetYaxis()->SetRangeUser(0.,4.5/*maxYRel*1.05*/);
        histo1DDummy2->DrawCopy();

        DrawGammaLines(0., maxR,1., 1.,1.1,kGray+2,2);
        DrawGammaLines(0., maxR,1.5, 1.5,1.1,kGray+2,2);
        DrawGammaLines(0., maxR,2., 2.,1.1,kGray+2,2);
        DrawGammaLines(0., maxR,2.5, 2.5,1.1,kGray+2,2);

        DrawGammaLines(5., 5.,0.,4.5,1.1,kGray+2,1);
        labelR5cm = new TLatex(0.11,0.88,label5cm.Data());
        SetStyleTLatex( labelR5cm,0.035,4);
        labelR5cm->Draw();

        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i<12) DrawGammaSetMarker(histRelUncWeightCut[i], marker[i], 1.5,color[i],color[i]);
            else     DrawGammaSetMarker(histRelUncWeightCut[i], marker[i], 1.5,color[i-12],color[i-12]);
            histRelUncWeightCut[i]->DrawCopy("same,c,p");
        }
//         legendWeight->Draw();
            legendWeightOnfly->Draw();
            legendWeightOffline->Draw();

        histo1DDummy2->DrawCopy("axis,same");
    canvasRelUncWeight->Update();
    canvasRelUncWeight->SaveAs(Form("%s/WeightRelUncertainty_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));
    delete canvasRelUncWeight;
//
    TCanvas* canvasWeightDiff = new TCanvas("canvasWeightDiff","",1300,900);
    DrawGammaCanvasSettings( canvasWeightDiff, 0.08, 0.03, 0.035, 0.09);

        Double_t minYDiff = -0.15;
        Double_t maxYDiff = 0.15; //qui
        TH2F * histo2DDummyDiff = new TH2F("","",1000,0.,180.,1000,-1.5,1.);
        SetStyleHistoTH2ForGraphs(histo2DDummyDiff, "R (cm)", "modified - standard",0.04,0.04, 0.04,0.04, 1.,1.);
        histo2DDummyDiff->GetYaxis()->SetRangeUser(minYDiff,maxYDiff);
        histo2DDummyDiff->Draw("copy");

        DrawGammaLines(0., maxR,1., 1.,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,0.02, 0.02,1.1,kGray+1,2);
        DrawGammaLines(0., maxR,-0.02, -0.02,1.1,kGray+1,2);

        DrawGammaLines(5., 5.,minYDiff,maxYDiff,1.1,kGray+2,1);
        labelR5cm = new TLatex(0.11,0.88,label5cm.Data());
        SetStyleTLatex( labelR5cm,0.04,4);
        labelR5cm->Draw();

        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i<12) DrawGammaSetMarker(histoDiffWeightCut[i], marker[i], 1.5,color[i],color[i]);
            else     DrawGammaSetMarker(histoDiffWeightCut[i], marker[i], 1.5,color[i-12],color[i-12]);
            histoDiffWeightCut[i]->Draw("same");
        }
//         legendWeight->Draw();
        legendWeightOnfly->Draw();
        legendWeightOffline->Draw();

        histo2DDummyDiff->DrawCopy("axis,same");
    canvasWeightDiff->Update();
    canvasWeightDiff->SaveAs(Form("%s/WeightDiff_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));
    delete canvasWeightDiff;


   //_______________________ Ploting weights vs pT in each R bin On-the-fly finder__________________________

    TCanvas *canvasMBWeightEachROnfly          = new TCanvas("canvasMBWeighEachROnfly","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasMBWeightEachROnfly, 0, 0, 0, 0);
    canvasMBWeightEachROnfly->cd();
    TPad * padMBWeightEachROnfly               = new TPad("padMBWeightEachROnfly","",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    DrawGammaPadSettings( padMBWeightEachROnfly, 0, 0, 0, 0);
    padMBWeightEachROnfly->Divide(4,3,0.0,0.0);
    padMBWeightEachROnfly->Draw();
    Int_t place  = 0;
    TH2F *histoDummyWeightEachROnfly =  new TH2F("histoDummyWeightEachROnfly","histoDummyWeightEachROnfly",1000,-0.05,1.,1000,0.5,2.);
    SetStyleHistoTH2ForGraphs(histoDummyWeightEachROnfly, "#it{p}_{T}^{Min} (GeV/c)","w_i", 0.05,0.05, 0.05,0.05);

    TLegend* legendWeightPtOnfly= new TLegend(0.5,0.91-(1+counterOnfly)*0.05,0.8,0.91); //0.17,0.13,0.5,0.24);
    legendWeightPtOnfly->SetFillColor(0);
    legendWeightPtOnfly->SetMargin(0.17);
    legendWeightPtOnfly->SetLineColor(0);
    legendWeightPtOnfly->SetTextFont(42);
    legendWeightPtOnfly->SetTextSize(0.05);
    legendWeightPtOnfly->SetHeader("On-the-Fly V0 finder");
    for(Int_t i = 0; i< 2; i++){
        if(V0ReaderName[i].Contains("On-the-Fly"))
            legendWeightPtOnfly->AddEntry(histWeight[i],Form("%s, %s",periodName[i].Data(),generatorName[i].Data()));
    }

    for(Int_t j=0; j < nBinsR; j++){
        place  = place + 1;
        padMBWeightEachROnfly->cd(place);
        padMBWeightEachROnfly->cd(place)->SetTopMargin(0.04);
        padMBWeightEachROnfly->cd(place)->SetBottomMargin(0.15);
        padMBWeightEachROnfly->cd(place)->SetLeftMargin(0.15);
        padMBWeightEachROnfly->cd(place)->SetRightMargin(0.05);
        histoDummyWeightEachROnfly->GetYaxis()->SetNdivisions(504);

//         if(j>1 && j<=3)      histoDummyWeightEachROnfly->GetYaxis()->SetRangeUser(0.9,1.3);
//         else if(j>3 && j<=6) histoDummyWeightEachROnfly->GetYaxis()->SetRangeUser(0.9,1.2);
//         else if(j>7 )        histoDummyWeightEachROnfly->GetYaxis()->SetRangeUser(0.9,1.1);
        histoDummyWeightEachROnfly->GetYaxis()->SetRangeUser(0.82,1.5);

        histoDummyWeightEachROnfly->DrawCopy();
        DrawGammaLines(0.,1.,1., 1.,1.,kGray,1);
        DrawGammaLines(0.4,0.4,0.82,1.3,1.,kGray,2);
        for(Int_t i = 0; i< 2; i++){
            if(i<12) DrawGammaSetMarker(histoWeightsEachRPtMin[j][i], marker[i], 0.5,color[i],color[i]);
            else     DrawGammaSetMarker(histoWeightsEachRPtMin[j][i], marker[i], 0.5,color[i-12],color[i-12]);
	    if(j<2) histoWeightsEachRPtMin[j][i]->SetLineStyle(3);
	    histoWeightsEachRPtMin[j][i]->Draw("same");
	    // if(i<12) DrawGammaSetMarker(histoWeightsEachRPtMinSecSub[j][i], marker[i], 0.5,6,color[i]);
            // else     DrawGammaSetMarker(histoWeightsEachRPtMinSecSub[j][i], marker[i], 0.5,6,color[i-12]);
            // histoWeightsEachRPtMinSecSub[j][i]->Draw("same");
        }
        if (j==0)  legendWeightPtOnfly->Draw();
    }

    canvasMBWeightEachROnfly->Print(Form("%s/MBWeightVSPtMinOnflyEachR_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));



 //_______________________ Ploting weights vs pT in each R bin Offline finder__________________________

    TCanvas *canvasMBWeightEachROffline          = new TCanvas("canvasMBWeighEachROffline","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasMBWeightEachROffline, 0, 0, 0, 0);
    canvasMBWeightEachROffline->cd();
    TPad * padMBWeightEachROffline               = new TPad("padMBWeightEachROffline","",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    DrawGammaPadSettings( padMBWeightEachROffline, 0, 0, 0, 0);
    padMBWeightEachROffline->Divide(4,3,0.0,0.0);
    padMBWeightEachROffline->Draw();
     place  = 0;
    TH2F *histoDummyWeightEachROffline =  new TH2F("histoDummyWeightEachROffline","histoDummyWeightEachROffline",1000,-0.05,1.,1000,0.5,2.);
    SetStyleHistoTH2ForGraphs(histoDummyWeightEachROffline, "#it{p}_{T}^{Min} (GeV/c)","w_i", 0.05,0.05, 0.05,0.05);

    TLegend* legendWeightPtOffline= new TLegend(0.5,0.31-(1+counterOnfly)*0.05,0.8,0.31); //0.17,0.13,0.5,0.24);
    legendWeightPtOffline->SetFillColor(0);
    legendWeightPtOffline->SetMargin(0.17);
    legendWeightPtOffline->SetLineColor(0);
    legendWeightPtOffline->SetTextFont(42);
    legendWeightPtOffline->SetTextSize(0.05);
    legendWeightPtOffline->SetHeader("Offline V0 finder");
    for(Int_t i = 2; i< NumberOfCuts; i++){
        if(V0ReaderName[i].Contains("Offline"))
            legendWeightPtOffline->AddEntry(histWeight[i],Form("%s, %s",periodName[i].Data(),generatorName[i].Data()));
    }

    for(Int_t j=0; j < nBinsR; j++){
        place  = place + 1;
        padMBWeightEachROffline->cd(place);
        padMBWeightEachROffline->cd(place)->SetTopMargin(0.04);
        padMBWeightEachROffline->cd(place)->SetBottomMargin(0.15);
        padMBWeightEachROffline->cd(place)->SetLeftMargin(0.15);
        padMBWeightEachROffline->cd(place)->SetRightMargin(0.05);
        histoDummyWeightEachROffline->GetYaxis()->SetNdivisions(504);

//         if(j>1 && j<=3)      histoDummyWeightEachROffline->GetYaxis()->SetRangeUser(0.9,1.3);
//         else if(j>3 && j<=6) histoDummyWeightEachROffline->GetYaxis()->SetRangeUser(0.9,1.2);
//         else if(j>7 )       histoDummyWeightEachROffline->GetYaxis()->SetRangeUser(0.9,1.1);
        histoDummyWeightEachROffline->GetYaxis()->SetRangeUser(0.82,1.5);

        histoDummyWeightEachROffline->DrawCopy();
        DrawGammaLines(0.,1.,1., 1.,1.,kGray,1);
        DrawGammaLines(0.4,0.4,0.85, 1.5,1.,kGray,2);
        for(Int_t i = 2; i< NumberOfCuts; i++){
	  if(i<12) DrawGammaSetMarker(histoWeightsEachRPtMin[j][i], marker[i], 0.5,color[i],color[i]);
	  else     DrawGammaSetMarker(histoWeightsEachRPtMin[j][i], marker[i], 0.5,color[i-12],color[i-12]);
	  if(j<2) histoWeightsEachRPtMin[j][i]->SetLineStyle(3);
	  histoWeightsEachRPtMin[j][i]->Draw("same");
	  // if(i<12) DrawGammaSetMarker(histoWeightsEachRPtMinSecSub[j][i], marker[i], 0.5,6,color[i]);
	  // else     DrawGammaSetMarker(histoWeightsEachRPtMinSecSub[j][i], marker[i], 0.5,6,color[i-12]);
	  // histoWeightsEachRPtMinSecSub[j][i]->Draw("same");
        }

        if (j==0)  legendWeightPtOffline->Draw();
    }
    canvasMBWeightEachROffline->Print(Form("%s/MBWeightVSPtMinOfflineEachR_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));


   //____________________________________________________________________________________

   //_______________________ Ploting weights vs pT in each R bin On-the-fly finder Secondary Subtracted usingCocktail__________________________

    TCanvas *canvasMBWeightEachROnflySecSubCock          = new TCanvas("canvasMBWeighEachROnflySecSubCock","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasMBWeightEachROnflySecSubCock, 0, 0, 0, 0);
    canvasMBWeightEachROnflySecSubCock->cd();
    TPad * padMBWeightEachROnflySecSubCock               = new TPad("padMBWeightEachROnflySecSubCock","",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    DrawGammaPadSettings( padMBWeightEachROnflySecSubCock, 0, 0, 0, 0);
    padMBWeightEachROnflySecSubCock->Divide(4,3,0.0,0.0);
    padMBWeightEachROnflySecSubCock->Draw();
    place  = 0;
    TH2F *histoDummyWeightEachROnflySecSubCock =  new TH2F("histoDummyWeightEachROnflySecSubCock","histoDummyWeightEachROnflySecSubCock",1000,-0.05,1.,1000,0.5,2.);
    SetStyleHistoTH2ForGraphs(histoDummyWeightEachROnflySecSubCock, "#it{p}_{T}^{Min} (GeV/c)","purity,w_i", 0.05,0.05, 0.05,0.05);

    TLegend* legendWeightPtOnflySecSubCock= new TLegend(0.4,0.91-(1+counterOnfly)*0.05,0.8,0.91); //0.17,0.13,0.5,0.24);
    legendWeightPtOnflySecSubCock->SetFillColor(0);
    legendWeightPtOnflySecSubCock->SetMargin(0.17);
    legendWeightPtOnflySecSubCock->SetLineColor(0);
    legendWeightPtOnflySecSubCock->SetTextFont(42);
    legendWeightPtOnflySecSubCock->SetTextSize(0.05);
    legendWeightPtOnflySecSubCock->SetHeader("On-the-Fly V0 finder SecSubCock");
    for(Int_t i = 0; i< 2; i++){
        if(V0ReaderName[i].Contains("On-the-Fly"))
            legendWeightPtOnflySecSubCock->AddEntry(histWeight[i],Form("%s, %s",periodName[i].Data(),generatorName[i].Data()));
    }

    for(Int_t j=0; j < nBinsR; j++){
        place  = place + 1;
        padMBWeightEachROnflySecSubCock->cd(place);
        padMBWeightEachROnflySecSubCock->cd(place)->SetTopMargin(0.04);
        padMBWeightEachROnflySecSubCock->cd(place)->SetBottomMargin(0.15);
        padMBWeightEachROnflySecSubCock->cd(place)->SetLeftMargin(0.15);
        padMBWeightEachROnflySecSubCock->cd(place)->SetRightMargin(0.05);
        histoDummyWeightEachROnflySecSubCock->GetYaxis()->SetNdivisions(504);

//         if(j>1 && j<=3)      histoDummyWeightEachROnfly->GetYaxis()->SetRangeUser(0.9,1.3);
//         else if(j>3 && j<=6) histoDummyWeightEachROnfly->GetYaxis()->SetRangeUser(0.9,1.2);
//         else if(j>7 )        histoDummyWeightEachROnfly->GetYaxis()->SetRangeUser(0.9,1.1);
        histoDummyWeightEachROnflySecSubCock->GetYaxis()->SetRangeUser(0.82,1.35);

        histoDummyWeightEachROnflySecSubCock->DrawCopy();
        DrawGammaLines(0.,1.,1., 1.,1.,kGray,1);
        DrawGammaLines(0.4,0.4,0.82,1.3,1.,kGray,2);
        for(Int_t i = 0; i< 2; i++){
            if(i<12) DrawGammaSetMarker(histoWeightsEachRPtMinSecSubUsingCocktail[j][i], marker[i], 0.5,color[i],color[i]);
            else     DrawGammaSetMarker(histoWeightsEachRPtMinSecSubUsingCocktail[j][i], marker[i], 0.5,color[i-12],color[i-12]);
	    if(j<2) histoWeightsEachRPtMinSecSubUsingCocktail[j][i]->SetLineStyle(3);
	    histoWeightsEachRPtMinSecSubUsingCocktail[j][i]->Draw("same");

	    if(i<12) DrawGammaSetMarker(histoPurityPtEachRBin[j][i], marker[i], 0.5,color[i],color[i]);
            else     DrawGammaSetMarker(histoPurityPtEachRBin[j][i], marker[i], 0.5,color[i-12],color[i-12]);
	    histoPurityPtEachRBin[j][i]->SetLineStyle(3);
	    histoPurityPtEachRBin[j][i]->Draw("same");
	    // if(i<12) DrawGammaSetMarker(histoWeightsEachRPtMinSecSub[j][i], marker[i], 0.5,6,color[i]);
            // else     DrawGammaSetMarker(histoWeightsEachRPtMinSecSub[j][i], marker[i], 0.5,6,color[i-12]);
            // histoWeightsEachRPtMinSecSub[j][i]->Draw("same");
        }
        if (j==0)  legendWeightPtOnflySecSubCock->Draw();
    }

    canvasMBWeightEachROnflySecSubCock->Print(Form("%s/MBWeightVSPtMinOnflyEachRSecSubCock_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));



 //_______________________ Ploting weights vs pT in each R bin Offline findersecondary subtracted using cocktail__________________________

    TCanvas *canvasMBWeightEachROfflineSecSubCock          = new TCanvas("canvasMBWeighEachROfflineSecSubCock","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasMBWeightEachROfflineSecSubCock, 0, 0, 0, 0);
    canvasMBWeightEachROfflineSecSubCock->cd();
    TPad * padMBWeightEachROfflineSecSubCock               = new TPad("padMBWeightEachROfflineSecSubCock","",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    DrawGammaPadSettings( padMBWeightEachROfflineSecSubCock, 0, 0, 0, 0);
    padMBWeightEachROfflineSecSubCock->Divide(4,3,0.0,0.0);
    padMBWeightEachROfflineSecSubCock->Draw();
    place  = 0;
    TH2F *histoDummyWeightEachROfflineSecSubCock =  new TH2F("histoDummyWeightEachROfflineSecSubCock","histoDummyWeightEachROfflineSecSubCock",1000,-0.05,1.,1000,0.5,2.);
    SetStyleHistoTH2ForGraphs(histoDummyWeightEachROfflineSecSubCock, "#it{p}_{T}^{Min} (GeV/c)","purity,w_i", 0.05,0.05, 0.05,0.05);

    TLegend* legendWeightPtOfflineSecSubCock = new TLegend(0.4,0.31-(1+counterOnfly)*0.05,0.8,0.31); //0.17,0.13,0.5,0.24);
    legendWeightPtOfflineSecSubCock->SetFillColor(0);
    legendWeightPtOfflineSecSubCock->SetMargin(0.17);
    legendWeightPtOfflineSecSubCock->SetLineColor(0);
    legendWeightPtOfflineSecSubCock->SetTextFont(42);
    legendWeightPtOfflineSecSubCock->SetTextSize(0.05);
    legendWeightPtOfflineSecSubCock->SetHeader("Offline V0 finder SecSubCock");
    for(Int_t i = 2; i< NumberOfCuts; i++){
        if(V0ReaderName[i].Contains("Offline"))
            legendWeightPtOfflineSecSubCock->AddEntry(histWeight[i],Form("%s, %s",periodName[i].Data(),generatorName[i].Data()));
    }

    for(Int_t j=0; j < nBinsR; j++){
        place  = place + 1;
        padMBWeightEachROfflineSecSubCock->cd(place);
        padMBWeightEachROfflineSecSubCock->cd(place)->SetTopMargin(0.04);
        padMBWeightEachROfflineSecSubCock->cd(place)->SetBottomMargin(0.15);
        padMBWeightEachROfflineSecSubCock->cd(place)->SetLeftMargin(0.15);
        padMBWeightEachROfflineSecSubCock->cd(place)->SetRightMargin(0.05);
        histoDummyWeightEachROfflineSecSubCock->GetYaxis()->SetNdivisions(504);

//         if(j>1 && j<=3)      histoDummyWeightEachROffline->GetYaxis()->SetRangeUser(0.9,1.3);
//         else if(j>3 && j<=6) histoDummyWeightEachROffline->GetYaxis()->SetRangeUser(0.9,1.2);
//         else if(j>7 )       histoDummyWeightEachROffline->GetYaxis()->SetRangeUser(0.9,1.1);
        histoDummyWeightEachROfflineSecSubCock->GetYaxis()->SetRangeUser(0.82,1.35);

        histoDummyWeightEachROfflineSecSubCock->DrawCopy();
        DrawGammaLines(0.,1.,1., 1.,1.,kGray,1);
        DrawGammaLines(0.4,0.4,0.85, 1.5,1.,kGray,2);
        for(Int_t i = 2; i< NumberOfCuts; i++){
	  if(i<12) DrawGammaSetMarker(histoWeightsEachRPtMinSecSubUsingCocktail[j][i], marker[i], 0.5,color[i],color[i]);
	  else     DrawGammaSetMarker(histoWeightsEachRPtMinSecSubUsingCocktail[j][i], marker[i], 0.5,color[i-12],color[i-12]);
	  if(j<2)histoWeightsEachRPtMinSecSubUsingCocktail[j][i]->SetLineStyle(3);
	  histoWeightsEachRPtMinSecSubUsingCocktail[j][i]->Draw("same");
	  if(i<12) DrawGammaSetMarker(histoPurityPtEachRBin[j][i], marker[i], 0.5,color[i],color[i]);
	  else     DrawGammaSetMarker(histoPurityPtEachRBin[j][i], marker[i], 0.5,color[i-12],color[i-12]);
	  histoPurityPtEachRBin[j][i]->SetLineStyle(3);
	  histoPurityPtEachRBin[j][i]->Draw("same");

	  // if(i<12) DrawGammaSetMarker(histoWeightsEachRPtMinSecSub[j][i], marker[i], 0.5,6,color[i]);
	  // else     DrawGammaSetMarker(histoWeightsEachRPtMinSecSub[j][i], marker[i], 0.5,6,color[i-12]);
	  // histoWeightsEachRPtMinSecSub[j][i]->Draw("same");
        }

        if (j==0)  legendWeightPtOfflineSecSubCock->Draw();
    }
    canvasMBWeightEachROfflineSecSubCock->Print(Form("%s/MBWeightVSPtMinOfflineEachRSecSubCock_%s.%s",outputDir.Data(),cutVariationName.Data(),suffix.Data()));


   //____________________________________________________________________________________









    TFile outFile(Form("%s/weightsWithErrors_Test.root",outputDir.Data()) ,"RECREATE");
    histWeightStaErr->Write();
    histWeightSysErr->Write();
    histWeightTotErr->Write();
    fProfileContainingMaterialBudgetWeights->Write();
    fProfileContainingMaterialBudgetWeightsUpTotErr->Write();
    fProfileContainingMaterialBudgetWeightsDownTotErr->Write();
    outFile.Close();

 }
