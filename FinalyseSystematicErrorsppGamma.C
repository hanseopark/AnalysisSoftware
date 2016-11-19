// ------------------------------------------------
// This code finalizes the systematic errors for --
// the direct photon analysis. Gamma spectrum as --
// well as double ratio systematics are produced --
// ------------------------------------------------
// Code is maintained by:                        --
//          - Nicolas Schmidt                    --
//          - Lucas Altenkämper                  --
// for the Photon Conversion Group Heidelberg    --
// ------------------------------------------------
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

void FinalyseSystematicErrorsppGamma(const char* nameDataFileErrors ="", TString energy="", TString spectrumName = "", Int_t numberOfPtBins =1 ,Int_t numberCutStudies=1, Int_t offSetBeginning = 0, TString suffix = "pdf", Int_t mode = 0){

    // Set plotting style and labels
    StyleSettingsThesis();
    SetPlotStyle();

    Double_t textSizeSpectra                    = 0.04;
    TString date = ReturnDateString();
    TString collisionSystem                     = ReturnFullCollisionsSystem(energy);
    TString detectionProcess                    = ReturnFullTextReconstructionProcess(mode);
    TString dateForOutput                       = ReturnDateStringForOutput();
    TString energyForOutput                     = energy;

    Int_t color[20]                             = {860,894,807,880,418,403,802,923,634,432,404,435,420,407,416,830,404,608,920,1};
    TLatex *labelGamma                          = new TLatex(0.75,0.88,detectionProcess);
    SetStyleTLatex( labelGamma, 0.038,4);

    TLatex *labelEnergy                         = new TLatex(0.75,0.92,collisionSystem);
    SetStyleTLatex( labelEnergy, 0.038,4);

    Int_t numberOfEntriesPos                    = 0;
    Int_t numberOfEntriesNeg                    = 0;

    TFile* fileErrorInput                       = new TFile(nameDataFileErrors);

    const Int_t nPtBins                         = numberOfPtBins;
    const Int_t nCuts                           = numberCutStudies;
    Double_t* ptBins                            = 0;
    Double_t* ptBinsErr                         = 0;

    // Set names of cut variations for legends
    TString nameCutVariation[16]           = {"R_{min}","dE/dx e-line","TPC cluster","Single e^{#pm} p_{T}",  "#chi^{2} #gamma & #psi_{pair}","q_{T}","Double Count","BG method","Periods","SPD pileup","Out of bunch pileup","cos(#Theta_{point})","dE/dx #pi-line","#alpha meson","Cocktail"};

    // Set names of cut variations for file input
    TString nameCutVariationSC[16]         = {"Rcut","dEdxE", "TPCCluster", "SinglePt", "Chi2" , "Qt" , "DoubleCount" ,"BG" ,"Periods" ,"SPD" ,"Pileup","CosPoint","dEdxPi", "Alpha" ,"Cocktail"};

    // Create output folder
    gSystem->Exec("mkdir -p GammaSystematicErrorsCalculated");

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
        errorsPosSummed[l]                      = 0.;
        errorsNegSummed[l]                      = 0.;
        errorsMeanSummed[l]                     = 0.;
        errorsPosCorrSummed[l]                  = 0.;
        errorsNegCorrSummed[l]                  = 0.;
        errorsMeanCorrSummed[l]                 = 0.;
    }

    for (Int_t i = 0; i < nCuts; i++){
        TGraphAsymmErrors* graphPosErrors       = NULL;
        TGraphAsymmErrors* graphNegErrors       = NULL;

        // Set currently undetermined uncertainties
        if ( nameCutVariationSC[i].CompareTo("Periods")==0  ){
            TString nameGraphPos;
            TString nameGraphNeg;
            nameGraphPos                        = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[i-1].Data()  );
            nameGraphNeg                        = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[i-1].Data()  );
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("Pileup")==0  ){
            TString nameGraphPos;
            TString nameGraphNeg;
            nameGraphPos                        = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[i-1].Data()  );
            nameGraphNeg                        = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[i-1].Data()  );
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("Cocktail")==0  ){
            TString nameGraphPos;
            TString nameGraphNeg;
            nameGraphPos                        = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[i-1].Data()  );
            nameGraphNeg                        = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[i-1].Data()  );
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else {
            // Load input graphs from systematics file
            TString nameGraphPos                = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[i].Data() );
            TString nameGraphNeg                = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[i].Data() );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        }
        // Remove first points depending on chosen offset
        for (Int_t j = 0; j < offSetBeginning; j++){
            graphPosErrors->RemovePoint(0);
            graphNegErrors->RemovePoint(0);
        }
        if (i == 0) {
            ptBins                              = graphNegErrors->GetX();
            ptBinsErr                           = graphNegErrors->GetEXhigh();
        }
        errorsNeg[i]                            = graphNegErrors->GetY();
        errorsNegErr[i]                         = graphNegErrors->GetEYhigh();
        errorsPos[i]                            = graphPosErrors->GetY();
        errorsPosErr[i]                         = graphPosErrors->GetEYhigh();

        // Calculate systematic error from input spectrum
        cout << nameCutVariationSC[i].Data() << endl;
        CalculateMeanSysErr(            errorsMean[i],  errorsMeanErr[i],   errorsPos[i],       errorsNeg[i],           nPtBins);	
        CorrectSystematicErrorsWithMean(errorsPos[i],   errorsPosErr[i],    errorsPosCorr[i],   errorsPosErrCorr[i],    nPtBins);
        CorrectSystematicErrorsWithMean(errorsNeg[i],   errorsNegErr[i],    errorsNegCorr[i],   errorsNegErrCorr[i],    nPtBins);
        CorrectSystematicErrorsWithMean(errorsMean[i],  errorsMeanErr[i],   errorsMeanCorr[i],  errorsMeanErrCorr[i],   nPtBins);

        // Add systematic error contribution from current cutvariation to total summed error
        for (Int_t l = 0; l < nPtBins; l++){
            errorsPosSummed[l]                  = errorsPosSummed[l]+pow(errorsPos[i][l],2);
            errorsNegSummed[l]                  = errorsNegSummed[l]+ pow(errorsNeg[i][l],2);
            errorsMeanSummed[l]                 = errorsMeanSummed[l]+ pow(errorsMean[i][l],2);
            errorsPosCorrSummed[l]              = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
            errorsNegCorrSummed[l]              = errorsNegCorrSummed[l] +pow(errorsNegCorr[i][l],2);
            errorsMeanCorrSummed[l]             = errorsMeanCorrSummed[l]+ pow(errorsMeanCorr[i][l],2);
        }
        negativeErrors[i]                       = new TGraphErrors(nPtBins,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
        meanErrors[i]                           = new TGraphErrors(nPtBins,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
        positiveErrors[i]                       = new TGraphErrors(nPtBins,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
        negativeErrorsCorr[i]                   = new TGraphErrors(nPtBins,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
        meanErrorsCorr[i]                       = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
        positiveErrorsCorr[i]                   = new TGraphErrors(nPtBins,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );

    }

    // Set the material budget error
    Double_t errorMaterial                      = 4.50;
    if (spectrumName.CompareTo("NoMatRatio") == 0)
        errorMaterial                           = 0.;

    for (Int_t l = 0; l < nPtBins; l++){
        errorsPosSummed[l]                      = pow(errorsPosSummed[l],0.5);
        errorsMeanSummed[l]                     = pow(errorsMeanSummed[l],0.5);
        errorsPosErrSummed[l]                   = errorsPosSummed[l]*0.001;
        errorsMeanErrSummed[l]                  = errorsMeanSummed[l]*0.001;
        errorsNegSummed[l]                      = -pow(errorsNegSummed[l],0.5);
        errorsNegErrSummed[l]                   = errorsNegSummed[l]*0.001;
        errorsPosCorrMatSummed[l]               = pow(errorsPosCorrSummed[l]+ pow(errorMaterial ,2.),0.5);
        errorsMeanCorrMatSummed[l]              = pow(errorsMeanCorrSummed[l]+ pow(errorMaterial ,2.),0.5);
        errorsNegCorrMatSummed[l]               = -pow(errorsNegCorrSummed[l]+ pow(errorMaterial ,2.),0.5);

        errorsPosCorrSummed[l]                  = pow(errorsPosCorrSummed[l],0.5);
        errorsMeanCorrSummed[l]                 = pow(errorsMeanCorrSummed[l],0.5);
        errorsPosErrCorrSummed[l]               = errorsPosCorrSummed[l]*0.001;
        errorsMeanErrCorrSummed[l]              = errorsMeanCorrSummed[l]*0.001;
        errorsMeanErrCorrMatSummed[l]           = errorsMeanCorrMatSummed[l]*0.001;
        errorsNegCorrSummed[l]                  = -pow(errorsNegCorrSummed[l],0.5);
        errorsNegErrCorrSummed[l]               = errorsNegCorrSummed[l]*0.001;
    }

    Double_t errorsMat[nPtBins];
    for (Int_t l = 0; l < nPtBins; l++){
        errorsMat[l]                            = errorMaterial;
    }
    TGraphErrors* graphMaterialError            = new TGraphErrors(nPtBins,ptBins ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );

    negativeErrorsSummed                        = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
    negativeErrorsCorrSummed                    = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
    positiveErrorsSummed                        = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
    positiveErrorsCorrSummed                    = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
    meanErrorsSummed                            = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
    meanErrorsCorrSummed                        = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
    meanErrorsCorrSummedIncMat                  = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );

//++++++++++++++++++++++++++++++ PLOTTING OF SYSMEAN +++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TCanvas* canvasSysErrMean                   = new TCanvas("canvasSysErrMean","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasSysErrMean, 0.08, 0.01, 0.015, 0.09);
    TH2D *histo2DSysErrMean                     = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+2,1000.,0.,25.);
    SetStyleHistoTH2ForGraphs(histo2DSysErrMean, "#it{p}_{T} (GeV/#it{c})","mean systematic Err %",
                              0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.0,0.9);
    histo2DSysErrMean->Draw();

    TLegend* legendMean                         = new TLegend(0.20,0.6,0.62,0.95);
    legendMean->SetTextSize(0.035);
    legendMean->SetFillColor(0);
    legendMean->SetBorderSize(0);
    legendMean->SetNColumns(2);
    for(Int_t i = 0; i< numberCutStudies ; i++){
        DrawGammaSetMarkerTGraphErr(meanErrors[i], 20+i, 1.,color[i],color[i]);
        meanErrors[i]->Draw("p,csame");
        legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
    }
    DrawGammaSetMarkerTGraphErr(meanErrorsSummed, 20, 1.,1,1);
    meanErrorsSummed->Draw("p,csame");
    legendMean->AddEntry(meanErrorsSummed,"quad. sum.","p");
    legendMean->Draw();
    labelEnergy->Draw();
    labelGamma->Draw();
    canvasSysErrMean->Update();
    canvasSysErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysMean_%s_%s_%s.%s",spectrumName.Data(), energy.Data(),dateForOutput.Data(),suffix.Data()));

    delete canvasSysErrMean;

//+++++++++++++++++++++++++ PLOTTING OF SYSMEANNEWITHMEAN ++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TCanvas* canvasNewSysErrMean                = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);

    TH2D *histo2DNewSysErrMean ;
    histo2DNewSysErrMean                        = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+2,1000.,0,25.);
    SetStyleHistoTH2ForGraphs(histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})","mean smoothed systematic Err %",
                              0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.0,0.9);
    histo2DNewSysErrMean->Draw();

    TLegend* legendMeanNew;
    legendMeanNew                               = new TLegend(0.15,0.7,0.67,0.95);
    legendMeanNew->SetTextSize(0.035);
    legendMeanNew->SetFillColor(0);
    legendMeanNew->SetBorderSize(0);
    legendMeanNew->SetNColumns(3);
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

    labelGamma->Draw();
    labelEnergy->Draw();

    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysMeanNewWithMean_%s_%s_%s.%s",spectrumName.Data(), energy.Data(),dateForOutput.Data(),suffix.Data()));

//+++++++++++++++++++++++++ SAVING SYSTEMATICS TO DAT FILE +++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatname                   = Form("GammaSystematicErrorsCalculated/SystematicError_%s_%s_%s.dat",spectrumName.Data(),energy.Data(),dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDat <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }
    SysErrDat.close();

//+++++++++++++++++++++++++ SAVING AVERAGE SYSTEMATICS TO DAT ++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMean               = Form("GammaSystematicErrorsCalculated/SystematicErrorAveraged_%s_%s_%s.dat",spectrumName.Data(),energy.Data(),dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDatAver  << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
        // SysErrDatAver << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }
    SysErrDatAver.close();

//+++++++++++++++++++++++++ SAVING DETAILED SYSTEMATICS ++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanSingleErr      = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedSinglePCM_%s_%s_%s.dat",spectrumName.Data(),energyForOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAverSingle;
    cout << SysErrDatnameMeanSingleErr << endl;
    SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
    SysErrDatAverSingle << "Pt bin\t" ; 
    for (Int_t i= 0; i< numberCutStudies; i++){
        SysErrDatAverSingle << nameCutVariationSC[i] << "\t";
    }
    SysErrDatAverSingle << "Material";
    SysErrDatAverSingle << endl; 
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDatAverSingle << ptBins[l] << "\t";
        for (Int_t i= 0; i< numberCutStudies; i++){
            SysErrDatAverSingle << errorsMeanCorr[i][l] << "\t";
        }
        SysErrDatAverSingle << 4.5 << "\t";
        SysErrDatAverSingle << errorsMeanCorrMatSummed[l] << endl;
    }
    SysErrDatAverSingle.close();

    Double_t errorsMeanCorrPID[nPtBins];	
    Double_t errorsMeanCorrSignalExtraction[nPtBins];	
    Double_t errorsMeanCorrTrackReco[nPtBins];	
    Double_t errorsMeanCorrPhotonReco[nPtBins];	

    for (Int_t l=0; l< nPtBins; l++){
        //"0 Rcut"
        //"1 dEdxE"
        //"2 TPCCluster"
        //"3 SinglePt"
        //"4 Chi2"
        //"5 Qt"
        //"6 DoubleCount"
        //"7 BG"
        //"8 Periods"
        //"9 SPD"
        //"10 Pileup"
        //"11 CosPoint"
        //"12 dEdxPi"
        //"13 Alpha"
        //"14 Cocktail"

        // Signal extraction: Yield extraction, Alpha ,BG, MC Smearing
        errorsMeanCorrSignalExtraction[l]       =   TMath::Sqrt(pow(errorsMeanCorr[13][l],2)+   // Alpha
                                                                pow(errorsMeanCorr[7][l],2)+    // BG
                                                                pow(errorsMeanCorr[9][l],2)+    // SPD
                                                                pow(errorsMeanCorr[10][l],2));  // Pileup

        // PID: dEdxE, dEdxPi
        errorsMeanCorrPID[l]                    =   TMath::Sqrt(pow(errorsMeanCorr[1][l],2)+    // dEdxE
                                                                pow(errorsMeanCorr[12][l],2));  // dEdxPi

        // photon reco: Chi2+PsiPair, Qt, CosPoint, MinR
        errorsMeanCorrPhotonReco[l]             =   TMath::Sqrt(pow(errorsMeanCorr[0][l],2)+    // Rcut
                                                                pow(errorsMeanCorr[4][l],2)+    // Chi2 PsiPair
                                                                pow(errorsMeanCorr[5][l],2)+    // Qt
                                                                pow(errorsMeanCorr[6][l],2)+    // DoubleCount
                                                                pow(errorsMeanCorr[11][l],2));  // CosPoint

        // track reconstruction: TPCCluster, Single pT, Eta
        errorsMeanCorrTrackReco[l]              =   TMath::Sqrt(pow(errorsMeanCorr[2][l],2)+    // TPCCluster
                                                                pow(errorsMeanCorr[3][l],2));   // Single pT
    }
    TGraphErrors* meanErrorsPID                 = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPhotonReco          = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsSignalExtraction    = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsTrackReco           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );

//++++++++++++++++++++++++ PLOTTING OF SYSERRORSUMMEDVISU ++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TCanvas* canvasSummedErrMean                = new TCanvas("canvasSummedErrMean","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasSummedErrMean, 0.08, 0.01, 0.015, 0.09);
    TH2D *histo2DSummedErrMean                  = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+2,1000.,0,25.);
    SetStyleHistoTH2ForGraphs(histo2DSummedErrMean, "#it{p}_{T} (GeV/#it{c})","mean smoothed systematic Err %",
                              0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.0,0.9);
    histo2DSummedErrMean->Draw();

    TLegend* legendSummedMeanNew                = new TLegend(0.20,0.7,0.5,0.95);
    legendSummedMeanNew->SetTextSize(0.035);
    legendSummedMeanNew->SetFillColor(0);
    legendSummedMeanNew->SetBorderSize(0);

    DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, 1.,color[0],color[0]);
    meanErrorsSignalExtraction->Draw("p,csame");
    DrawGammaSetMarkerTGraphErr(meanErrorsPID, 21, 1.,color[1],color[1]);
    meanErrorsPID->Draw("p,csame");
    DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, 22, 1.,color[2],color[2]);
    meanErrorsTrackReco->Draw("p,csame");
    DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, 23, 1.,color[3],color[3]);
    meanErrorsPhotonReco->Draw("p,csame");
    DrawGammaSetMarkerTGraphErr(meanErrorsCorr[1], 25, 1.,color[5],color[5]);
    meanErrorsCorr[1]->Draw("p,csame");
    DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[4],color[4]);
    graphMaterialError->Draw("p,csame");

    legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
    legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
    legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reconstruction","p");
    legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reconstruction","p");
    legendSummedMeanNew->AddEntry(meanErrorsCorr[10],"Pileup Estimate","p");
    if(numberCutStudies>14)
        legendSummedMeanNew->AddEntry(meanErrorsCorr[14],"Cocktail Uncertainty","p");
    legendSummedMeanNew->AddEntry(graphMaterialError,"Material","p");
    DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
    meanErrorsCorrSummedIncMat->Draw("p,csame");
    legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
    legendSummedMeanNew->Draw();

    labelGamma->Draw();
    labelEnergy->Draw();

    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysErrorSummedVisu_%s_%s_%s.%s",spectrumName.Data(), energy.Data(),dateForOutput.Data(),suffix.Data()));

    delete canvasSummedErrMean;

//+++++++++++++++++++++++++ SAVING SYSTEMATICS PAPER STYLE +++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanPaper          = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedPaper_%s_%s_%s.dat",spectrumName.Data(),energy.Data(),dateForOutput.Data());
    fstream SysErrDatAverPaper;
    cout << SysErrDatnameMeanPaper << endl;
    SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
    SysErrDatAverPaper  << "p_{T}" << "\t Material \t Yield Extraction \t PID \t photon reco \t track recon \t summed" <<  endl;
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial*2 << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
    }
    SysErrDatAverPaper.close();
}
