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
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

void FinaliseSystematicErrorsConv_pPb(  const char* nameDataFileErrors      = "", 
                                        TString energy                      = "", 
                                        TString meson                       = "", 
                                        Int_t numberOfPtBins                = 1 ,
                                        Int_t numberCutStudies              = 1,
                                        Int_t offSetBeginning               = 0, 
                                        TString additionalName              = "",
                                        TString additionalNameOutput        = "",
                                        TString suffix = "eps"
                                     ){
    
    StyleSettingsThesis();	
    SetPlotStyle();
    
    TString date            = ReturnDateString();
    TString dateForOutput   = ReturnDateStringForOutput();
    TString energyForOutput = energy;
    energyForOutput.ReplaceAll(".","_");    
    
    Int_t color[20] = {860,894,807,880,418,403,802,923,634,432,404,435,420,407,416,830,404,608,920,1};
    
    Int_t 	numberOfEntriesPos = 0;
    Int_t 	numberOfEntriesNeg = 0;
    
    TFile* fileErrorInput= new TFile(nameDataFileErrors);
    const Int_t nPtBins = numberOfPtBins;
    const Int_t nCuts = numberCutStudies;
    Double_t* ptBins;
    Double_t* ptBinsErr;
    
    TString nameCutVariationpPb5023GeV[14] = {"Yield extraction","Pile-up","dE/dx e-line","dE/dx #pi-line","TPC cluster","Single e^{#pm} p_{T}",  "#chi^{2} #gamma","q_{T}","#alpha meson", "#psi_{pair}","cos(#Theta_{point})","BG method","MC smearing","#eta_{#gamma, e^{#pm}}"};  
    TString nameCutVariation[14];
    TString nameCutVariationSC[14];
    
    TString nameCutVariationSCpPb5023GeV[14] = {"YieldExtraction","BGEstimate","dEdxE","dEdxPi", "Cluster", "SinglePt", "Chi2", "Qt", "Alpha", "PsiPair","CosPoint","BG","Smearing","Eta"};
        
    if (energy.CompareTo("pPb_5.023TeV") == 0) {
        for (Int_t i = 0; i < numberCutStudies; i++){
            nameCutVariation[i] = nameCutVariationpPb5023GeV[i];
            nameCutVariationSC[i] = nameCutVariationSCpPb5023GeV[i];
        }
    }
    
    gSystem->Exec("mkdir -p SystematicErrorsCalculatedConv");
    
    Double_t* errorsNeg[nCuts];
    Double_t errorsNegCorr[nCuts][nPtBins];
    Double_t errorsNegSummed[nPtBins];
    Double_t errorsNegCorrSummed[nPtBins];
    Double_t errorsNegCorrMatSummed[nPtBins];
    
    Double_t* errorsNegErr[nCuts];
    Double_t errorsNegErrCorr[nCuts][nPtBins];
    Double_t errorsNegErrSummed[nPtBins];
    Double_t errorsNegErrCorrSummed[nPtBins];
    
    Double_t* errorsPos[nCuts];
    Double_t errorsPosCorr[nCuts][nPtBins];
    Double_t errorsPosSummed[nPtBins];
    Double_t errorsPosCorrSummed[nPtBins];
    Double_t errorsPosCorrMatSummed[nPtBins];
    
    Double_t* errorsPosErr[nCuts];
    Double_t errorsPosErrSummed[nPtBins];
    Double_t errorsPosErrCorr[nCuts][nPtBins];
    Double_t errorsPosErrCorrSummed[nPtBins];
    
    Double_t errorsMean[nCuts][nPtBins];
    Double_t errorsMeanCorr[nCuts][nPtBins];
    Double_t errorsMeanSummed[nPtBins];
    Double_t errorsMeanCorrSummed[nPtBins];
    Double_t errorsMeanCorrMatSummed[nPtBins];
    
    Double_t errorsMeanErr[nCuts][nPtBins];
    Double_t errorsMeanErrCorr[nCuts][nPtBins];
    Double_t errorsMeanErrSummed[nPtBins];
    Double_t errorsMeanErrCorrSummed[nPtBins];
    Double_t errorsMeanErrCorrMatSummed[nPtBins];
    
    TGraphErrors* negativeErrors[nCuts];
    TGraphErrors* negativeErrorsSummed;
    TGraphErrors* positiveErrors[nCuts];
    TGraphErrors* positiveErrorsSummed;
    TGraphErrors* negativeErrorsCorr[nCuts];
    TGraphErrors* negativeErrorsCorrSummed;
    TGraphErrors* positiveErrorsCorr[nCuts];
    TGraphErrors* positiveErrorsCorrSummed;
    TGraphErrors* meanErrors[nCuts];
    TGraphErrors* meanErrorsSummed;
    TGraphErrors* meanErrorsCorr[nCuts];
    TGraphErrors* meanErrorsCorrSummed;
    TGraphErrors* meanErrorsCorrSummedIncMat;
    
    for (Int_t l = 0; l < nPtBins; l++){
            errorsPosSummed[l] = 0.;
            errorsNegSummed[l] = 0.;
            errorsMeanSummed[l] = 0.;
            errorsPosCorrSummed[l] = 0.;
            errorsNegCorrSummed[l] = 0.;
            errorsMeanCorrSummed[l] = 0.;
        } 
    
    for (Int_t i = 0; i < nCuts; i++){
        TGraphAsymmErrors* graphPosErrors;
        TGraphAsymmErrors* graphNegErrors;
        if (i == 0){
            TString nameGraphPos = Form("%s_SystErrorRelPos_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
            TString nameGraphNeg = Form("%s_SystErrorRelNeg_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("BGEstimate")==0  ){
            cout << "including pileup for PbPb" << endl;
            TString nameGraphPos;
            TString nameGraphNeg;
            nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),"dEdxE",additionalNameOutput.Data()  );
            nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),"dEdxE",additionalNameOutput.Data()  );
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else {
            TString nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
            TString nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        }
        
        for (Int_t j = 0; j < offSetBeginning; j++){
            graphPosErrors->RemovePoint(0);
            graphNegErrors->RemovePoint(0);
        }
        if (i == 0) {
            ptBins = graphNegErrors->GetX();
            ptBinsErr = graphNegErrors->GetEXhigh();
        }
        errorsNeg[i] = graphNegErrors->GetY();
        errorsNegErr[i] = graphNegErrors->GetEYhigh();
        errorsPos[i] = graphPosErrors->GetY();
        errorsPosErr[i] = graphPosErrors->GetEYhigh();
            
        cout << nameCutVariationSC[i].Data() << endl;
        CalculateMeanSysErr(errorsMean[i], errorsMeanErr[i], errorsPos[i], errorsNeg[i], nPtBins);	
        CorrectSystematicErrorsWithMean(errorsPos[i],errorsPosErr[i], errorsPosCorr[i], errorsPosErrCorr[i], nPtBins);
        CorrectSystematicErrorsWithMean(errorsNeg[i],errorsNegErr[i], errorsNegCorr[i], errorsNegErrCorr[i], nPtBins);
        CorrectSystematicErrorsWithMean(errorsMean[i], errorsMeanErr[i], errorsMeanCorr[i], errorsMeanErrCorr[i], nPtBins);
        
        if (meson.CompareTo("Eta")==0 ){
            for (Int_t k = 0; k < nPtBins; k++){
                if (ptBins[k] == 0.8){
                    errorsMean[i][k] = errorsMean[i][k+1];
                    errorsMeanErr[i][k] = errorsMeanErr[i][k+1];
                    errorsMeanCorr[i][k] = errorsMeanCorr[i][k+1];
                    errorsMeanErrCorr[i][k] = errorsMeanErrCorr[i][k+1];
                }
            }	
        }
        
        
        if (nameCutVariationSC[i].CompareTo("PsiPair")==0 && meson.CompareTo("Pi0")==0){
            for (Int_t k = 0; k < nPtBins; k++){
                if (ptBins[k] > 3.0){
                    errorsMean[i][k] = 2.5;
                    errorsMeanErr[i][k] = 0.025;
                    errorsMeanCorr[i][k] = 2.5;
                    errorsMeanErrCorr[i][k] = 0.025;
                }
            }	
        }
        if (nameCutVariationSC[i].CompareTo("PsiPair")==0 && meson.CompareTo("Eta")==0 ){
            for (Int_t k = 0; k < nPtBins; k++){
// 				if (ptBins[k] > 3.0){
                    errorsMean[i][k] = 3;
                    errorsMeanErr[i][k] = 0.03;
                    errorsMeanCorr[i][k] = 3;
                    errorsMeanErrCorr[i][k] = 0.03;
// 				}
            }	
        }
        if (nameCutVariationSC[i].CompareTo("PsiPair")==0 && meson.CompareTo("Pi0EtaBinning")==0 ){
            for (Int_t k = 0; k < nPtBins; k++){
// 				if (ptBins[k] > 3.0){
                    errorsMean[i][k] = 1.5;
                    errorsMeanErr[i][k] = 0.015;
                    errorsMeanCorr[i][k] = 1.5;
                    errorsMeanErrCorr[i][k] = 0.015;
// 				}
            }	
        }

        if (nameCutVariationSC[i].CompareTo("CosPoint")==0  ){
            for (Int_t k = 0; k < nPtBins; k++){
// 				if (ptBins[k] > 3.0){
                    errorsMean[i][k] = 0.5;
                    errorsMeanErr[i][k] = 0.005;
                    errorsMeanCorr[i][k] = 0.5;
                    errorsMeanErrCorr[i][k] = 0.005;
// 				}
            }	
        }
        if (nameCutVariationSC[i].CompareTo("Alpha")==0){
            for (Int_t k = 0; k < nPtBins; k++){
// 				if (ptBins[k] > 3.0){
                    errorsMean[i][k] = 0.3;
                    errorsMeanErr[i][k] = 0.003;
                    errorsMeanCorr[i][k] = 0.3;
                    errorsMeanErrCorr[i][k] = 0.003;
// 				}
            }	
        }
        if (nameCutVariationSC[i].CompareTo("Eta")==0 && meson.Contains("Pi0")  ){
            for (Int_t k = 0; k < nPtBins; k++){
// 				if (ptBins[k] > 3.0){
                    errorsMean[i][k] = 1;
                    errorsMeanErr[i][k] = 0.01;
                    errorsMeanCorr[i][k] = 1;
                    errorsMeanErrCorr[i][k] = 0.01;
// 				}
            }	
        }
        if (nameCutVariationSC[i].CompareTo("Eta")==0 && meson.CompareTo("Eta")==0  ){
            for (Int_t k = 0; k < nPtBins; k++){
// 				if (ptBins[k] > 3.0){
                    errorsMean[i][k] = 3;
                    errorsMeanErr[i][k] = 0.03;
                    errorsMeanCorr[i][k] = 3;
                    errorsMeanErrCorr[i][k] = 0.03;
// 				}
            }	
        }

        if (nameCutVariationSC[i].CompareTo("Qt")==0 && meson.CompareTo("Eta")==0  ){
            for (Int_t k = 0; k < nPtBins; k++){
                if (ptBins[k] == 7.0){
                    errorsMean[i][k] = 3.5;
                    errorsMeanErr[i][k] = 0.035;
                    errorsMeanCorr[i][k] = 3.5;
                    errorsMeanErrCorr[i][k] = 0.035;
                }
            }	
        }
        
        if (nameCutVariationSC[i].CompareTo("BGEstimate")==0  ){
            for (Int_t k = 0; k < nPtBins; k++){
// 				if (ptBins[k] > 3.0){
                    errorsMean[i][k] = 0.3;
                    errorsMeanErr[i][k] = 0.003;
                    errorsMeanCorr[i][k] = 0.3;
                    errorsMeanErrCorr[i][k] = 0.003;
// 				}
            }	
        }
        
        for (Int_t l = 0; l < nPtBins; l++){
            errorsPosSummed[l] = errorsPosSummed[l]+pow(errorsPos[i][l],2);
            errorsNegSummed[l] = errorsNegSummed[l]+ pow(errorsNeg[i][l],2);
            errorsMeanSummed[l] = errorsMeanSummed[l]+ pow(errorsMean[i][l],2);
            errorsPosCorrSummed[l] = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
            errorsNegCorrSummed[l] = errorsNegCorrSummed[l] +pow(errorsNegCorr[i][l],2);
            errorsMeanCorrSummed[l] =errorsMeanCorrSummed[l]+ pow(errorsMeanCorr[i][l],2);
        }
        negativeErrors[i] = new TGraphErrors(nPtBins,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
        meanErrors[i] = new TGraphErrors(nPtBins,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
        positiveErrors[i] = new TGraphErrors(nPtBins,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
        negativeErrorsCorr[i] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
        meanErrorsCorr[i] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
        positiveErrorsCorr[i] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );
        
    }
    
    Double_t errorMaterial = 4.50;
    
    for (Int_t l = 0; l < nPtBins; l++){
        errorsPosSummed[l] = pow(errorsPosSummed[l],0.5);
        errorsMeanSummed[l] = pow(errorsMeanSummed[l],0.5);
        errorsPosErrSummed[l] = errorsPosSummed[l]*0.001;
        errorsMeanErrSummed[l] = errorsMeanSummed[l]*0.001;
        errorsNegSummed[l] = -pow(errorsNegSummed[l],0.5);
        errorsNegErrSummed[l] = errorsNegSummed[l]*0.001;
        errorsPosCorrMatSummed[l] = pow(errorsPosCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
        errorsMeanCorrMatSummed[l] = pow(errorsMeanCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
        errorsNegCorrMatSummed[l] = -pow(errorsNegCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
        
        errorsPosCorrSummed[l] = pow(errorsPosCorrSummed[l],0.5);
        errorsMeanCorrSummed[l] = pow(errorsMeanCorrSummed[l],0.5);
        errorsPosErrCorrSummed[l] = errorsPosCorrSummed[l]*0.001;
        errorsMeanErrCorrSummed[l] = errorsMeanCorrSummed[l]*0.001;
        errorsMeanErrCorrMatSummed[l] = errorsMeanCorrMatSummed[l]*0.001;
        errorsNegCorrSummed[l] = -pow(errorsNegCorrSummed[l],0.5);
        errorsNegErrCorrSummed[l] = errorsNegCorrSummed[l]*0.001;
// 		cout << errorsMeanSummed[l] << "\t" << errorsMeanCorrMatSummed[l]<< endl;
    }
    
    Double_t errorsMat[nPtBins];
    for (Int_t l = 0; l < nPtBins; l++){
        errorsMat[l] = 2*errorMaterial;
        
    }
    TGraphErrors* graphMaterialError = new TGraphErrors(nPtBins,ptBins ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );
    
    cout << "here" << endl;
    negativeErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
    negativeErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
    positiveErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
    positiveErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
    meanErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
    meanErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
    meanErrorsCorrSummedIncMat = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );

    
    
    cout << "here" << endl;
    
    TCanvas* canvasSysErrMean = new TCanvas("canvasSysErrMean","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasSysErrMean, 0.08, 0.01, 0.015, 0.09);
    TH2D *histo2DSysErrMean ;
// 	if (meson.CompareTo("Pi0")==0 ){
        histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,30.);
// 	} else {
// 		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,50.);
// 	}	
    histo2DSysErrMean->SetYTitle("mean systematic Err %");
    histo2DSysErrMean->SetXTitle("p_{T} (GeV/c)");
    histo2DSysErrMean->GetYaxis()->SetNdivisions(510,kTRUE);
    histo2DSysErrMean->GetYaxis()->SetLabelSize(0.03);
    histo2DSysErrMean->GetYaxis()->SetTitleSize(0.04);
    histo2DSysErrMean->GetYaxis()->SetDecimals();
    histo2DSysErrMean->GetYaxis()->SetTitleOffset(0.9);
    histo2DSysErrMean->GetXaxis()->SetTitleOffset(1.);
    histo2DSysErrMean->GetXaxis()->SetLabelSize(0.03);
    histo2DSysErrMean->GetXaxis()->SetTitleSize(0.04);
    histo2DSysErrMean->SetTitle("");
    histo2DSysErrMean->Draw();
    
    TLegend* legendMean; 
    if (meson.Contains("Pi0")){
        legendMean= new TLegend(0.15,0.6,0.57,0.95);
    } else {
        legendMean= new TLegend(0.20,0.6,0.62,0.95);
    }
    legendMean->SetTextSize(0.035);
    legendMean->SetFillColor(0);
    legendMean->SetBorderSize(0);
    if (numberCutStudies> 9) legendMean->SetNColumns(2);
    for(Int_t i = 0; i< numberCutStudies ; i++){
        DrawGammaSetMarkerTGraphErr(meanErrors[i], 20+i, 1.,color[i],color[i]);
        meanErrors[i]->Draw("p,csame");
        legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
    }
    DrawGammaSetMarkerTGraphErr(meanErrorsSummed, 20, 1.,1,1);
    meanErrorsSummed->Draw("p,csame");
    legendMean->AddEntry(meanErrorsSummed,"quad. sum.","p");
    legendMean->Draw();
    canvasSysErrMean->Update();
    canvasSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysMean_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    delete canvasSysErrMean;
    
    TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);
    
    TH2D *histo2DNewSysErrMean ;
// 	if (meson.CompareTo("Pi0")==0 ){
        histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
// 	} else { 
// 		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
// 	}
    histo2DNewSysErrMean->SetYTitle("mean smoothed systematic Err %");
    histo2DNewSysErrMean->SetXTitle("p_{T} (GeV/c)");
    histo2DNewSysErrMean->GetYaxis()->SetNdivisions(510,kTRUE);
    histo2DNewSysErrMean->GetYaxis()->SetLabelSize(0.03);
    histo2DNewSysErrMean->GetYaxis()->SetTitleSize(0.04);
    histo2DNewSysErrMean->GetYaxis()->SetDecimals();
    histo2DNewSysErrMean->GetYaxis()->SetTitleOffset(0.9);
    histo2DNewSysErrMean->GetXaxis()->SetTitleOffset(1.);
    histo2DNewSysErrMean->GetXaxis()->SetLabelSize(0.03);
    histo2DNewSysErrMean->GetXaxis()->SetTitleSize(0.04);
    histo2DNewSysErrMean->SetTitle("");
    histo2DNewSysErrMean->Draw();
    
    TLegend* legendMeanNew;
    if (meson.Contains("Pi0")){
        legendMeanNew= new TLegend(0.15,0.6,0.57,0.95);
    } else {
        legendMeanNew= new TLegend(0.20,0.6,0.62,0.95);
    }	
    legendMeanNew->SetTextSize(0.035);
    legendMeanNew->SetFillColor(0);
    legendMeanNew->SetBorderSize(0);
    if (numberCutStudies> 9) legendMeanNew->SetNColumns(2);
    for(Int_t i = 0; i< numberCutStudies ; i++){
        DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], 20+i, 1.,color[i],color[i]);
        meanErrorsCorr[i]->Draw("p,csame");
        legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");	
    }
    
    DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummed, 20, 1.,1,1);
    meanErrorsCorrSummed->Draw("p,csame");
    legendMeanNew->AddEntry(meanErrorsCorrSummed,"quad. sum.","p");
    
    DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kRed,kRed);
    meanErrorsCorrSummedIncMat->Draw("p,csame");
    legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum., inc. mat.","p");
    legendMeanNew->Draw();
    
    TLatex *labelMeson;
    if (meson.Contains("Pi0")){
        labelMeson= new TLatex(0.65,0.89,Form("#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}"));
    } else {
        labelMeson= new TLatex(0.65,0.89,Form("#eta #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}"));
    }
    SetStyleTLatex( labelMeson, 0.038,4);
    labelMeson->Draw();

    TLatex *labelCentrality = new TLatex(0.65,0.93,Form("%s p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV"	,additionalName.Data()	));
    SetStyleTLatex( labelCentrality, 0.038,4);
    labelCentrality->Draw();

    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysMeanNewWithMean_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
// 	delete canvasNewSysErrMean;

    const char *SysErrDatname = Form("SystematicErrorsCalculatedConv/SystematicErrorPCM_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDat <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }

    SysErrDat.close();

    const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedConv/SystematicErrorAveragedPCM_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDatAver << ptBins[l] << "\t"<< "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }

    SysErrDatAver.close();
    
    const char *SysErrDatnameMeanSingleErr = Form("SystematicErrorsCalculatedConv/SystematicErrorAveragedSinglePCM_%s_%s%s_%s.dat", meson.Data(), energyForOutput.Data(), additionalNameOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAverSingle;
    cout << SysErrDatnameMeanSingleErr << endl;
    SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
    SysErrDatAverSingle << "Pt bin\t" ; 
    for (Int_t i= 0; i< numberCutStudies; i++){
        SysErrDatAverSingle << nameCutVariationSC[i] << "\t";
    }
    if(meson.CompareTo("EtaToPi0"))
        SysErrDatAverSingle << "Material";
    SysErrDatAverSingle << endl; 
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDatAverSingle << ptBins[l] << "\t";
        for (Int_t i= 0; i< numberCutStudies; i++){
            SysErrDatAverSingle << errorsMeanCorr[i][l] << "\t";
        }  
        if(meson.CompareTo("EtaToPi0"))
            SysErrDatAverSingle << 9 << "\t";
        SysErrDatAverSingle << errorsMeanCorrMatSummed[l] << endl;
    }
    SysErrDatAverSingle.close();
    
    Double_t errorsMeanCorrPID[nPtBins];	
    Double_t errorsMeanCorrSignalExtraction[nPtBins];	
    Double_t errorsMeanCorrTrackReco[nPtBins];	
    Double_t errorsMeanCorrPhotonReco[nPtBins];	
    
    for (Int_t l=0; l< nPtBins; l++){
// 			"YieldExtraction"-0,"BGEstimate"-1,"dEdxE"-2,"dEdxPi"-3, "Cluster"-4, "SinglePt"-5, "Chi2"-6, "Qt"-7, "Alpha"-8, "PsiPair"-9, "CosPoint"-10, "Eta"-11, "BG"-12, "MCSmearing"-13, "RCut" -14
        // grouping:
        // Pileup : BGEstimate 1
        // Signal extraction: Yield extraction 0, Alpha 8,  BG 12,MC Smearing 13
        if (numberCutStudies>12){
            errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+ errorsMeanCorr[12][l]*errorsMeanCorr[12][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
        } else if (numberCutStudies>9){
            errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[9][l]*errorsMeanCorr[9][l]);
            
        } else {
            errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]);
        }
        // PID: dEdxE 2, dEdxPi 3
        errorsMeanCorrPID[l] =TMath::Sqrt(errorsMeanCorr[2][l]*errorsMeanCorr[2][l]+ errorsMeanCorr[3][l]*errorsMeanCorr[3][l]);
        // photon reco: Chi2 6, Qt 7, PsiPair 9, CosPoint 10, MinR 14, 
        if (numberCutStudies>14){
            errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[6][l]* errorsMeanCorr[6][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l]+errorsMeanCorr[9][l]*errorsMeanCorr[9][l]+errorsMeanCorr[10][l]*errorsMeanCorr[10][l]+errorsMeanCorr[14][l]*errorsMeanCorr[14][l]);
        } else if (numberCutStudies>10){
            errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[6][l]* errorsMeanCorr[6][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l] +errorsMeanCorr[9][l]*errorsMeanCorr[9][l]+errorsMeanCorr[10][l]*errorsMeanCorr[10][l]);
        } else if (numberCutStudies>9){
            errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[6][l]* errorsMeanCorr[6][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l] +errorsMeanCorr[9][l]*errorsMeanCorr[9][l]);	
        } else {
            errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[6][l]* errorsMeanCorr[6][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l]);
        }
        // track reconstruction: Cluster 4, Single pt 5, Eta 11
        errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[5][l]*errorsMeanCorr[5][l]);				
        if (numberCutStudies>11) errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[5][l]*errorsMeanCorr[5][l]+errorsMeanCorr[11][l]*errorsMeanCorr[11][l]);				
            
    }
    TGraphErrors* meanErrorsPID = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPhotonReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsSignalExtraction = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsTrackReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    
    cout << "here" << endl;


    TCanvas* canvasSummedErrMean = new TCanvas("canvasSummedErrMean","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasSummedErrMean, 0.08, 0.01, 0.015, 0.09);
    TH2D *histo2DSummedErrMean ;
    
// 	if (meson.Contains("Pi0") ){	
        histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
// 	} else {
// 		histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
// 	}
    histo2DSummedErrMean->SetYTitle("mean smoothed systematic Err %");
    histo2DSummedErrMean->SetXTitle("p_{T} (GeV/c)");
    histo2DSummedErrMean->GetYaxis()->SetNdivisions(510,kTRUE);
    histo2DSummedErrMean->GetYaxis()->SetLabelSize(0.03);
    histo2DSummedErrMean->GetYaxis()->SetTitleSize(0.04);
    histo2DSummedErrMean->GetYaxis()->SetDecimals();
    histo2DSummedErrMean->GetYaxis()->SetTitleOffset(0.9);
    histo2DSummedErrMean->GetXaxis()->SetTitleOffset(1.);
    histo2DSummedErrMean->GetXaxis()->SetLabelSize(0.03);
    histo2DSummedErrMean->GetXaxis()->SetTitleSize(0.04);
    histo2DSummedErrMean->SetTitle("");
    histo2DSummedErrMean->Draw();
    
    TLegend* legendSummedMeanNew;
    if (meson.Contains("Pi0")){
        legendSummedMeanNew= new TLegend(0.15,0.7,0.45,0.95);
    } else {
        legendSummedMeanNew= new TLegend(0.20,0.7,0.5,0.95);
    }
    legendSummedMeanNew->SetTextSize(0.035);
    legendSummedMeanNew->SetFillColor(0);
    legendSummedMeanNew->SetBorderSize(0);
    
    
    DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, 1.,color[0],color[0]);
    meanErrorsSignalExtraction->Draw("p,csame");
    cout << "here" << endl;
    DrawGammaSetMarkerTGraphErr(meanErrorsPID, 21, 1.,color[1],color[1]);
    meanErrorsPID->Draw("p,csame");
    cout << "here" << endl;
    DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, 22, 1.,color[2],color[2]);
    meanErrorsTrackReco->Draw("p,csame");
    cout << "here" << endl;
    DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, 23, 1.,color[3],color[3]);
    meanErrorsPhotonReco->Draw("p,csame");
    cout << "here" << endl;
    DrawGammaSetMarkerTGraphErr(meanErrorsCorr[1], 25, 1.,color[5],color[5]);
    meanErrorsCorr[1]->Draw("p,csame");
    cout << "here" << endl;
    DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[4],color[4]);
    graphMaterialError->Draw("p,csame");
    cout << "here" << endl;
    
    legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
    legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
    legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reconstruction","p");
    legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reconstruction","p");
    legendSummedMeanNew->AddEntry(meanErrorsCorr[1],"Pileup Estimate","p");
    legendSummedMeanNew->AddEntry(graphMaterialError,"Material","p");
    DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
    meanErrorsCorrSummedIncMat->Draw("p,csame");
    legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
    legendSummedMeanNew->Draw();
    
    labelMeson->Draw();
    labelCentrality->Draw();
    
    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysErrorSummedVisu_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    delete canvasSummedErrMean;

    
    
    const char *SysErrDatnameMeanPaper = Form("SystematicErrorsCalculatedConv/SystematicErrorAveragedPaper_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAverPaper;
    cout << SysErrDatnameMeanPaper << endl;
    SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
    SysErrDatAverPaper  << "p_{T}" << "\t Material \t Yield Extraction \t PID \t photon reco \t track recon \t summed" <<  endl;
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial*2 << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
    }

    SysErrDatAverPaper.close();
        
}