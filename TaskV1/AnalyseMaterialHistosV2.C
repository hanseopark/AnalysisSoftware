// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
//This file is not supposed to be run on outputfiles of the GammaConv-Software before the 30th Sept 2010.

#include <stdlib.h>
#include <iostream>
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
#include "TProfile2D.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "../CommonHeaders/PlottingMeson.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "AnalyseMaterialHistos.h"


TH1F* ScaleByIntegralWithinLimits(const TH1F* hist, Double_t rmin, Double_t rmax, Double_t mcGasCorrectionFactor )
{
    TH1F * histToBeScaled = (TH1F*)hist->Clone();
    Double_t nconvInGas =  histToBeScaled->Integral(hist->GetXaxis()->FindBin(rmin),hist->GetXaxis()->FindBin(rmax),"width");
    GammaScalingHistogramm(histToBeScaled,1./(mcGasCorrectionFactor*nconvInGas));
    cout << "nconvInGas = " << nconvInGas << "  sqrt(nconvInGas) = " << TMath::Sqrt(nconvInGas) << " sqrt(nconvInGas)/nconvInGas = "<<  TMath::Sqrt(nconvInGas)/nconvInGas<<endl;
    return histToBeScaled;
}


void AnalyseMaterialHistosV2( TString fileName         = "",
                            TString fileNameMC       = "",
                            TString cutSelection     = "",
                            TString optionEnergy     = "",
                            TString optionPeriod     = "",
                            TString suffix           = "pdf",
                            TString outputFolderName = "MaterialBudgetHistos",
                            Int_t mode               = 0)
{

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    StyleSettingsThesis();
    SetPlotStyle();

    fEnergyFlag         = optionEnergy;
    fCollisionSystem    = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    fDetectionProcess   = ReturnFullTextReconstructionProcess(mode);

    fCutSelection = cutSelection;
    TString fCutSelectionRead = cutSelection;
    cout << "Cut selection:: " << fCutSelectionRead.Data()<< endl;

    TString outputDirectory= Form("%s/%s/%s",outputFolderName.Data(),fEnergyFlag.Data(),fCutSelectionRead.Data());
    gSystem->Exec("mkdir -p "+outputDirectory);
    if(optionPeriod.CompareTo(""))
        outputDirectory= Form("%s/%s/%s/%s",outputFolderName.Data(),fEnergyFlag.Data(),optionPeriod.Data(),fCutSelectionRead.Data());


    // factor needed for Run1 simulations. LHC10b
    //static const Float_t mcGasCorrectionFactor = 0.960035693454;
    Double_t mcGasCorrectionFactor = 1.;
    if (optionPeriod.CompareTo("LHC10b") == 0 || optionPeriod.CompareTo("LHC10bc") == 0 || optionEnergy.Contains("PbPb") ){
      mcGasCorrectionFactor = 0.960035693454;
    }

    TFile fData(fileName.Data());
    TString nameMainDir = "GammaConvMaterial";
    TList *TopDir =(TList*)fData.Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    TList *HistosGammaConversion       = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversion == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
        return;
    }
    TList *ESDContainer                = (TList*)HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));


    TFile fMC(fileNameMC.Data());
    TList *TopDirMC =(TList*)fMC.Get(nameMainDir.Data());
    if(TopDirMC == NULL){
        cout<<"ERROR: TopDirMC not Found"<<endl;
        return;
    }
    TList *HistosGammaConversionMC       = (TList*)TopDirMC->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversionMC == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in FileMC"<<endl;
        return;
    }
    TList *ESDContainerMC             = (TList*)HistosGammaConversionMC->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    TList *MCContainer                = (TList*)HistosGammaConversionMC->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
    TList *TrueConversionContainerMC  = (TList*)HistosGammaConversionMC->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));

    TString nameNTrackHisto = "GoodESDTracks09";
    // TString nameNTrackHisto = "GoodESDTracks14";
    // TString nameNTrackHisto = "GoodESDTracks09_14";

    TH1F* histoEventQualityData          = (TH1F*)ESDContainer->FindObject("NEvents");
    TH2F *histoMappingRPhiData           = (TH2F*)ESDContainer->FindObject("ESD_ConversionMapping_RPhi");
    TH2F *histoMappingRZData             = (TH2F*)ESDContainer->FindObject("ESD_ConversionMapping_RZ");
    TH1F *histoMappingRData              = (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_R");
    TH1F *histoMappingMidPtRData         = (TH1F*)ESDContainer->FindObject("ESD_ConversionMappingMidPt_R");
    TH1F *histoMappingPtData             = (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_Pt");
    TH1F *histoMappingPt5cmData          = (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_Pt5cm");
    TH2F *histoMappingPtRData            = (TH2F*)ESDContainer->FindObject("ESD_ConversionMapping_Pt_R");
    TH1F *histoMappingDCAData            = (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_DCA");
    TH1F *histoMappingPsiPairData        = (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_PsiPair");
    TH1F *histoMappingChi2Data           = (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_Chi2");
    TH1F *histoMappingMassData           = (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_Mass");
    TH1F *histoNumberOfGoodESDTracksData = (TH1F*)ESDContainer->FindObject(nameNTrackHisto.Data());
    TH1F *histoMappingHighPtRData        = (TH1F*)ESDContainer->FindObject("ESD_ConversionMappingHighPt_R");

    TH1F* histoEventQualityMC            = (TH1F*)ESDContainerMC->FindObject("NEvents");
    TH2F *histoMappingRPhiMC             = (TH2F*)ESDContainerMC->FindObject("ESD_ConversionMapping_RPhi");
    TH2F *histoMappingRZMC               = (TH2F*)ESDContainerMC->FindObject("ESD_ConversionMapping_RZ");
    TH1F *histoMappingRMC                = (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_R");
    TH1F *histoMappingMidPtRMC           = (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMappingMidPt_R");
    TH1F *histoMappingPtMC               = (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Pt");
    TH1F *histoMappingPt5cmMC            = (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Pt5cm");
    TH2F *histoMappingPtRMC              = (TH2F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Pt_R");
    TH1F *histoMappingDCAMC              = (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_DCA");
    TH1F *histoMappingPsiPairMC          = (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_PsiPair");
    TH1F *histoMappingChi2MC             = (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Chi2");
    TH1F *histoMappingMassMC             = (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Mass");

    TH2F *histoMappingPtTrueMC           = (TH2F*)TrueConversionContainerMC->FindObject("ESD_TrueConversionMapping_Pt");
    TH2F *histoMappingPtTrue5cmMC        = (TH2F*)TrueConversionContainerMC->FindObject("ESD_TrueConversionMapping_Pt5cm");
    TH1F *histoMappingRTrueMC            = (TH1F*)TrueConversionContainerMC->FindObject("ESD_TrueConversionMapping_R");
    TH1F *histoMappingMidPtRTrueMC       = (TH1F*)TrueConversionContainerMC->FindObject("ESD_TrueConversionMappingMidPt_R");
    TH1F *histoMappingRTruePi0DalMC      = (TH1F*)TrueConversionContainerMC->FindObject("ESD_TruePi0DalConversionMapping_R");
    TH1F *histoMappingRTrueEtaDalMC      = (TH1F*)TrueConversionContainerMC->FindObject("ESD_TrueEtaDalConversionMapping_R");
    TH1F *histoMappingRTrueCombMC        = (TH1F*)TrueConversionContainerMC->FindObject("ESD_TrueCombConversionMapping_R");
    TH1F *histoNumberOfGoodESDTracksMC   = (TH1F*)ESDContainerMC->FindObject(nameNTrackHisto.Data());
    TH1F *histoMappingHighPtRMC          = (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMappingHighPt_R");

    Int_t startBin  = histoMappingPtRData->GetXaxis()->FindBin(fMinPt+0.001); //0.
    Int_t endBin    = histoMappingPtRData->GetXaxis()->FindBin(fMaxPt-0.001);  //20.
    cout << "bin range:: " << startBin << " " << endBin << endl;

    TH1F *histoMappingRFullPtData        = (TH1F*)histoMappingPtRData->ProjectionY("histoMappingRFullPtData",startBin,endBin,"e");
    TH1F * histoMappingRFullPtMC         = (TH1F*)histoMappingPtRMC->ProjectionY("histoMappingRFullPtMC",startBin,endBin,"e");

    Int_t startBin1  = histoMappingPtRData->GetXaxis()->FindBin(fMinPt1+0.001); //0.8
    TH1F * histoMappingRMinPt1Data           = (TH1F*)histoMappingPtRData->ProjectionY("histoMappingRMinPt1Data",startBin1,endBin,"e");

    //_____________________________________normalize and rebin histograms___________________________________
    Float_t numberGoodEventsData         =  histoEventQualityData->GetBinContent(1);
    Double_t meanMultiplitcityData       =  histoNumberOfGoodESDTracksData->GetMean();
    Float_t normFactorReconstData        = 1./numberGoodEventsData*1./meanMultiplitcityData;
    // normFactorReconstData=          1./numberGoodEventsData;
    //normFactorReconstData=1;

    Float_t numberGoodEventsMC           = histoEventQualityMC->GetBinContent(1);
    Double_t meanMultiplitcityMC         = histoNumberOfGoodESDTracksMC->GetMean();
    Float_t normFactorReconstMC          = 1./numberGoodEventsMC*1./meanMultiplitcityMC;
    // normFactorReconstMC=          1./numberGoodEventsMC;
    //    normFactorReconstMC=1;
    cout << "Normalization factor Data/MC:: data " << meanMultiplitcityData<< " MC " << meanMultiplitcityMC<< " -> " << meanMultiplitcityData/meanMultiplitcityMC<< endl;

    GammaScalingHistogramm(histoNumberOfGoodESDTracksData,normFactorReconstData);
    GammaScalingHistogramm(histoNumberOfGoodESDTracksMC,normFactorReconstMC);

    GammaScalingHistogramm(histoMappingRFullPtData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingRMinPt1Data,normFactorReconstData);
    GammaScalingHistogramm(histoMappingRFullPtMC,normFactorReconstMC);

    GammaScalingHistogramm(histoMappingRData,normFactorReconstData);
    ConvGammaRebinWithBinCorrection(histoMappingRData,rebinRPlots);
    GammaScalingHistogramm(histoMappingRMC,normFactorReconstMC);
    ConvGammaRebinWithBinCorrection(histoMappingRMC,rebinRPlots);
    GammaScalingHistogramm(histoMappingRTrueMC,normFactorReconstMC);
    ConvGammaRebinWithBinCorrection(histoMappingRTrueMC,rebinRPlots);
    GammaScalingHistogramm(histoMappingRTruePi0DalMC,normFactorReconstMC);
    ConvGammaRebinWithBinCorrection(histoMappingRTruePi0DalMC,rebinRPlots);
    GammaScalingHistogramm(histoMappingRTrueEtaDalMC,normFactorReconstMC);
    ConvGammaRebinWithBinCorrection(histoMappingRTrueEtaDalMC,rebinRPlots);
    GammaScalingHistogramm(histoMappingRTrueCombMC,normFactorReconstMC);
    ConvGammaRebinWithBinCorrection(histoMappingRTrueCombMC,rebinRPlots);

    GammaScalingHistogramm(histoMappingMidPtRData,normFactorReconstData);
    ConvGammaRebinWithBinCorrection(histoMappingMidPtRData,rebinRPlots);
    GammaScalingHistogramm(histoMappingMidPtRMC,normFactorReconstMC);
    ConvGammaRebinWithBinCorrection(histoMappingMidPtRMC,rebinRPlots);
    GammaScalingHistogramm(histoMappingMidPtRTrueMC,normFactorReconstMC);
    ConvGammaRebinWithBinCorrection(histoMappingMidPtRTrueMC,rebinRPlots);

    GammaScalingHistogramm(histoMappingPtData,normFactorReconstData);
    ConvGammaRebinWithBinCorrection(histoMappingPtData,rebinPtPlots);
    GammaScalingHistogramm(histoMappingPtMC,normFactorReconstMC);
    ConvGammaRebinWithBinCorrection(histoMappingPtMC,rebinPtPlots);
    GammaScalingHistogramm(histoMappingPtTrueMC,normFactorReconstMC);
    ConvGammaRebinWithBinCorrection(histoMappingPtTrueMC,rebinPtPlots);
    GammaScalingHistogramm(histoMappingPt5cmData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingPt5cmMC,normFactorReconstMC);

    GammaScalingHistogramm(histoMappingPtTrue5cmMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingDCAData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingDCAMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingPsiPairData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingPsiPairMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingChi2Data,normFactorReconstData);
    GammaScalingHistogramm(histoMappingChi2MC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingMassData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingMassMC,normFactorReconstMC);


    //normalize with Nconversions in the gas
    cout << "scaling histoMappingRData:" << endl;
    histoMappingRDataScaledToGas       = ScaleByIntegralWithinLimits(histoMappingRData, rMinGas, rMaxGas, 1. );
    cout << "scaling histoMappingRFullPtData:" << endl;
    histoMappingRFullPtDataScaledToGas = ScaleByIntegralWithinLimits(histoMappingRFullPtData, rMinGas, rMaxGas, 1. );
    cout << "scaling histoMappingRMC:" << endl;
    histoMappingRMCScaledToGas         = ScaleByIntegralWithinLimits(histoMappingRMC, rMinGas, rMaxGas, mcGasCorrectionFactor );
    cout << "scaling histoMappingRFullPtMC:" << endl;
    histoMappingRFullPtMCScaledToGas   = ScaleByIntegralWithinLimits(histoMappingRFullPtMC, rMinGas, rMaxGas, mcGasCorrectionFactor );

    histoMappingRDataScaledToGasRebin  = (TH1F*)histoMappingRDataScaledToGas->Rebin(nBinsR,"histoMappingRDataScaledToGasRebin", arrayRBins);
    histoMappingRMCScaledToGasRebin    = (TH1F*)histoMappingRMCScaledToGas->Rebin(nBinsR,"histoMappingRMCScaledToGasRebin", arrayRBins);
    histoMappingDataMCRatioRScaledToGas           = (TH1F*)histoMappingRDataScaledToGasRebin->Clone("histoMappingDataMCRatioRScaledToGas");
    histoMappingDataMCRatioRScaledToGas->Divide(histoMappingRDataScaledToGasRebin,histoMappingRMCScaledToGasRebin,1.,1.,"B");

    histoMappingRFullPtDataScaledToGasRebin  = (TH1F*)histoMappingRFullPtDataScaledToGas->Rebin(nBinsR,"histoMappingRFullPtDataScaledToGasRebin", arrayRBins);
    histoMappingRFullPtMCScaledToGasRebin    = (TH1F*)histoMappingRFullPtMCScaledToGas->Rebin(nBinsR,"histoMappingRFullPtMCScaledToGasRebin", arrayRBins);
    histoMappingDataMCRatioRFullPtScaledToGas           = (TH1F*)histoMappingRFullPtDataScaledToGasRebin->Clone("histoMappingDataMCRatioRFullPtScaledToGas");
    histoMappingDataMCRatioRFullPtScaledToGas->Divide(histoMappingRFullPtDataScaledToGasRebin,histoMappingRFullPtMCScaledToGasRebin,1.,1.,"B");

    histoMappingDataMCRatioR = (TH1F*)histoMappingRData->Clone("histoMappingDataMCRatioR");
    histoMappingDataMCRatioR->Divide(histoMappingDataMCRatioR,histoMappingRMC,1.,1.,"B");
    histoMappingMidPtDataMCRatioR = (TH1F*)histoMappingMidPtRData->Clone("histoMappingDataMCRatioR");
    histoMappingMidPtDataMCRatioR->Divide(histoMappingMidPtRData,histoMappingMidPtRMC,1.,1.,"B");

    histoMappingPurityR = (TH1F*)histoMappingRTrueMC->Clone("histoMappingPurityR");
    histoMappingPurityR->Sumw2();
    histoMappingPurityR->Divide(histoMappingRTrueMC,histoMappingRMC,1.,1.,"B");

    histoMappingPurityPt = (TH1F*)histoMappingPtTrueMC->Clone("histoMappingPurityPt");
    histoMappingPurityPt->Sumw2();
    histoMappingPurityPt->Divide(histoMappingPtTrueMC,histoMappingPtMC,1.,1.,"B");

    histoMappingPurityPt5cm = (TH1F*)histoMappingPtTrue5cmMC->Clone("histoMappingPurityPt5cm");
    histoMappingPurityPt5cm->Sumw2();
    histoMappingPurityPt5cm->Divide(histoMappingPtTrue5cmMC,histoMappingPt5cmMC,1.,1.,"B");


    for(Int_t iR = 0; iR < nBinsR; iR++){

      histoMappingPhiInRData[iR]        = (TH1D*)histoMappingRPhiData->ProjectionX(Form("histoMappingPhiInRData_%i",iR), histoMappingRPhiData->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingRPhiData->GetYaxis()->FindBin(arrayRBins[iR+1]));
      GammaScalingHistogramm(histoMappingPhiInRData[iR],normFactorReconstData);
      ConvGammaRebinWithBinCorrection(histoMappingPhiInRData[iR],rebinPhiPlots);

      histoMappingPhiInRMC[iR]          = (TH1D*)histoMappingRPhiMC->ProjectionX(Form("histoMappingPhiInRMC_%i",iR), histoMappingRPhiMC->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingRPhiMC->GetYaxis()->FindBin(arrayRBins[iR+1]));
      GammaScalingHistogramm(histoMappingPhiInRMC[iR],normFactorReconstMC);
      ConvGammaRebinWithBinCorrection(histoMappingPhiInRMC[iR],rebinPhiPlots);

      histoMappingDataMCRatioPhiInR[iR] = (TH1D*)histoMappingPhiInRData[iR]->Clone(Form("histoMappingDataMCRatioPhiInR_%02d",iR));
      histoMappingDataMCRatioPhiInR[iR]->Divide(histoMappingPhiInRData[iR],histoMappingPhiInRMC[iR]);

      if (iR > 9) rebinZPlots = 4;
      histoMappingZInRData[iR]          =  (TH1D*)histoMappingRZData->ProjectionX( Form("histoMappingZInRData_%i",iR), histoMappingRZData->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingRZData->GetYaxis()->FindBin(arrayRBins[iR+1]));
      GammaScalingHistogramm(histoMappingZInRData[iR],normFactorReconstData);
      ConvGammaRebinWithBinCorrection(histoMappingZInRData[iR],rebinZPlots);

      histoMappingZInRMC[iR]            =  (TH1D*)histoMappingRZMC->ProjectionX( Form("histoMappingZInRMC_%i",iR), histoMappingRZMC->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingRZMC->GetYaxis()->FindBin(arrayRBins[iR+1]));
      GammaScalingHistogramm(histoMappingZInRMC[iR],normFactorReconstMC);
      ConvGammaRebinWithBinCorrection(histoMappingZInRMC[iR],rebinZPlots);

      histoMappingDataMCRatioZInR[iR]         = (TH1D*)histoMappingZInRData[iR]->Clone(Form("histoMappingDataMCRatioZInR_%02d",iR));
      histoMappingDataMCRatioZInR[iR]->Divide(histoMappingDataMCRatioZInR[iR],histoMappingZInRMC[iR]);

    }

    //______________________________________ Multiplicity _________________________________________
    TCanvas * canvasNTracks = new TCanvas("canvasNTracks","",1200,1000);
    DrawGammaCanvasSettings( canvasNTracks, 0.1, 0.03, 0.05, 0.09);
    //canvasNTracks->Divide(2,2);
    canvasNTracks->SetLogy(1);
    TH2F * histoDummyNTracks = new TH2F("histoDummyNTracks","histoDummyNTracks",1000,0.,2000.,1000,1.e-10,10);
		SetStyleHistoTH2ForGraphs(histoDummyNTracks, "Good TPC tracks","Counts", 0.035,0.04,0.035,0.04,1.,1.);
// 		histoDummyNTracks->GetYaxis()->SetRangeUser(2.e-9, 1.);
		histoDummyNTracks->DrawCopy();

		DrawGammaSetMarker(histoNumberOfGoodESDTracksData, 20, 1.5, colorData, colorData);
        histoNumberOfGoodESDTracksData->Draw("same,hist");
		DrawGammaSetMarker(histoNumberOfGoodESDTracksMC, 20, 1.5, colorMC, colorMC);
        histoNumberOfGoodESDTracksMC->Draw("same,hist");

        TLegend* legend = GetAndSetLegend(0.75,0.75,2);
        legend->AddEntry(histoNumberOfGoodESDTracksData,"Data","l");
        legend->AddEntry(histoNumberOfGoodESDTracksMC,"MC","l");
        legend->Draw();

    canvasNTracks->Print(Form("%s/NGoodESDTracks%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));



    //______________________________________ R distrib scaled to gas __________________________________
    Double_t arrayX[2];
    Double_t arrayY[3];
    Double_t relX[3];
    Double_t relY[3];
    Int_t textSizeLabels = 30;
    Double_t textsizeLabelsDown = 0;
    Double_t textsizeFacDown = 0;
    Double_t textsizeLabelsUp = 0;
    Double_t textsizeFacUp = 0;
    ReturnCorrectValuesForCanvasScaling(1200, 1000, 1, 2, 0.1, 0.025, 0.025, 0.08, arrayX, arrayY, relX, relY);
    TCanvas * canvasTwoPanels = new TCanvas("canvasTwoPanels","",1200,1000);

    TPad* padUpper = new TPad("","",arrayX[0], arrayY[1], arrayX[1], arrayY[0],-1, -1, -2);
    DrawGammaPadSettings( padUpper, relX[0], relX[2], relY[0], relY[1]);
    TPad* padLower = new TPad("","",arrayX[0], arrayY[2], arrayX[1], arrayY[1],-1, -1, -2);
    DrawGammaPadSettings( padLower, relX[0], relX[2], relY[1], relY[2]);
    Double_t marginXRatio = relX[0]*1200;
        padUpper->Draw();
        if (padUpper->XtoPixel(padUpper->GetX2()) < padUpper->YtoPixel(padUpper->GetY1())){
            textsizeLabelsUp = (Double_t)textSizeLabels/padUpper->XtoPixel(padUpper->GetX2()) ;
            textsizeFacUp = (Double_t)1./padUpper->XtoPixel(padUpper->GetX2()) ;
        } else {
            textsizeLabelsUp = (Double_t)textSizeLabels/padUpper->YtoPixel(padUpper->GetY1());
            textsizeFacUp = (Double_t)1./padUpper->YtoPixel(padUpper->GetY1());
        }
        padLower->Draw();
        if (padLower->XtoPixel(padLower->GetX2()) < padLower->YtoPixel(padLower->GetY1())){
            textsizeLabelsDown = (Double_t)textSizeLabels/padLower->XtoPixel(padLower->GetX2()) ;
            textsizeFacDown = (Double_t)1./padLower->XtoPixel(padLower->GetX2()) ;
        } else {
            textsizeLabelsDown = (Double_t)textSizeLabels/padLower->YtoPixel(padLower->GetY1());
            textsizeFacDown = (Double_t)1./padLower->YtoPixel(padLower->GetY1());
        }

        TH2F *histoDummyTwoPanelsUp = new TH2F("histoDummyTwoPanelsUp","histoDummyTwoPanelsUp",1000,0.,180.,1000,-0.1,5.);
        SetStyleHistoTH2ForGraphs(histoDummyTwoPanelsUp, "R (cm)","Counts", 0.9*textsizeLabelsUp, textsizeLabelsUp,0.9*textsizeLabelsUp,textsizeLabelsUp, 1,0.15/(textsizeFacUp*marginXRatio));
//         histoDummyTwoPanelsUp->GetYaxis()->SetLabelOffset(+0.01);
        TH2F *histoDummyTwoPanelsDown =  new TH2F("histoDummyTwoPanelsDown","histoDummyTwoPanelsDown",1000,0.,180.,1000,0.,2.);
        SetStyleHistoTH2ForGraphs(histoDummyTwoPanelsDown, "R (cm)","#frac{Data}{MC} ", 0.9*textsizeLabelsDown, textsizeLabelsDown,0.9*textsizeLabelsDown,textsizeLabelsDown, 1,0.15/(textsizeFacDown*marginXRatio));
//         histoDummyTwoPanelsDown->GetYaxis()->SetLabelOffset(+0.01);
        histoDummyTwoPanelsDown->GetYaxis()->SetRangeUser(0.3,1.55);

    padUpper->cd();
    histoDummyTwoPanelsUp->DrawCopy();

        DrawGammaSetMarker(histoMappingRDataScaledToGas, 20, 1.5, colorData, colorData);
        histoMappingRDataScaledToGas->Draw("same,hist");
		DrawGammaSetMarker(histoMappingRMCScaledToGas, 20, 1.5, colorMC, colorMC);
        histoMappingRMCScaledToGas->Draw("same,hist");

        TLegend* legenUpPanel = GetAndSetLegend2(0.7, 0.8-(2*0.9*textsizeLabelsUp), 0.9, 0.8, textSizeLabels);
        legenUpPanel->AddEntry(histoMappingRDataScaledToGas,"Data","l");
        legenUpPanel->AddEntry(histoMappingRMCScaledToGas,"MC","l");
        legenUpPanel->Draw();

    padLower->cd();
    histoDummyTwoPanelsDown->DrawCopy();

        DrawGammaLines(0.,180,1., 1.,0.1,kGray);

        DrawGammaSetMarker(histoMappingDataMCRatioRScaledToGas, 20, 1.5, colorData, colorData);
        histoMappingDataMCRatioRScaledToGas->Draw("same");
		DrawGammaSetMarker(histoMappingDataMCRatioRFullPtScaledToGas, 24, 1.5, colorData, colorData);
        histoMappingDataMCRatioRFullPtScaledToGas->Draw("same");

        TLegend* legenLowPanel = GetAndSetLegend2(0.7, 0.85-(2*0.9*textsizeLabelsDown), 0.9, 0.85, textSizeLabels);
        legenLowPanel->AddEntry(histoMappingDataMCRatioRScaledToGas,"From proj","lp");
        legenLowPanel->AddEntry(histoMappingDataMCRatioRFullPtScaledToGas,"Full pT","lp");
        legenLowPanel->Draw();


    histoDummyTwoPanelsDown->Draw("same,axis");
    canvasTwoPanels->Print(Form("%s/PhotonConvRScaledToGas%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    // TCanvas * canvasRScaledToGas= new TCanvas("canvasRPlot","",1200,1000);
    // canvasRScaledToGas->SetLogy(1);
    // histoMappingRData->SetXTitle("R(cm)");
    // histoMappingRData->SetMaximum(2e-3);
    // histoMappingRData->SetMinimum(1e-6);
    // histoMappingRData->SetLineColor(4);
    // histoMappingRData->DrawCopy("hist");
    // histoMappingRFullPtData->SetLineColor(2);
    // histoMappingRFullPtData->DrawCopy("hist,same");
    // histoMappingRMinPt1Data->SetLineColor(6);
    // histoMappingRMinPt1Data->DrawCopy("hist,same");
    // canvasRScaledToGas->Print(Form("%s/RScaledtoGasData_%s.pdf",outputDirectory.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //----------------------------------------------

    TLatex* latexCutName= new TLatex(doubleLatexNamingCutX, doubleLatexNamingCutY,fCutSelectionRead.Data());
    SetStyleTLatex( latexCutName, sizeTextNameBins,2);

    //______________________________________ R distribution ______________________________
    Double_t arrayXRdistrib[3];
    Double_t arrayYRdistrib[3];
    Double_t relXRdistrib[3];
    Double_t relYRdistrib[3];
    textSizeLabels = 50;
    textsizeLabelsUp = 0;
    textsizeFacUp = 0;
    textsizeLabelsDown = 0;
    textsizeFacDown = 0;
    ReturnCorrectValuesForCanvasScaling(2500,2000, 2, 2,0.06, 0.025, 0.025,0.06,arrayXRdistrib,arrayYRdistrib,relXRdistrib,relYRdistrib);
    TCanvas* canvasPhotonR = new TCanvas("canvasPhotonR","",0,0,2500,2000);

    TPad* padRDistribLowerLeft = new TPad("padRDistribLowerLeft", "", arrayXRdistrib[0], arrayYRdistrib[2], arrayXRdistrib[1], arrayYRdistrib[1],-1, -1, -2);
    DrawGammaPadSettings( padRDistribLowerLeft, relXRdistrib[0], relXRdistrib[1], relYRdistrib[1], relYRdistrib[2]);
    padRDistribLowerLeft->Draw();
    TPad* padRDistribLowerRight = new TPad("padRDistribLowerRight", "", arrayXRdistrib[1], arrayYRdistrib[2], arrayXRdistrib[2], arrayYRdistrib[1],-1, -1, -2);
    DrawGammaPadSettings( padRDistribLowerRight, relXRdistrib[1], relXRdistrib[2], relYRdistrib[1], relYRdistrib[2]);
    padRDistribLowerRight->Draw();
    TPad* padRDistribUpperLeft = new TPad("padRDistribUpperLeft", "", arrayXRdistrib[0], arrayYRdistrib[1], arrayXRdistrib[1], arrayYRdistrib[0],-1, -1, -2);
    DrawGammaPadSettings( padRDistribUpperLeft, relXRdistrib[0], relXRdistrib[1], relYRdistrib[0], relYRdistrib[1]);
    padRDistribUpperLeft->Draw();
    TPad* padRDistribUpperRight = new TPad("padRDistribUpperRight", "", arrayXRdistrib[1], arrayYRdistrib[1], arrayXRdistrib[2], arrayYRdistrib[0],-1, -1, -2);
    DrawGammaPadSettings( padRDistribUpperRight, relXRdistrib[1], relXRdistrib[2], relYRdistrib[0], relYRdistrib[1]);
    padRDistribUpperRight->Draw();
    if (padRDistribUpperLeft->XtoPixel(padRDistribUpperLeft->GetX2()) < padRDistribUpperLeft->YtoPixel(padRDistribUpperLeft->GetY1())){
        textsizeLabelsUp = (Double_t)textSizeLabels/padRDistribUpperLeft->XtoPixel(padRDistribUpperLeft->GetX2()) ;
        textsizeFacUp = (Double_t)1./padRDistribUpperLeft->XtoPixel(padRDistribUpperLeft->GetX2()) ;
    } else {
        textsizeLabelsUp = (Double_t)textSizeLabels/padRDistribUpperLeft->YtoPixel(padRDistribUpperLeft->GetY1());
        textsizeFacUp = (Double_t)1./padRDistribUpperLeft->YtoPixel(padRDistribUpperLeft->GetY1());
    }
    if (padRDistribLowerLeft->XtoPixel(padRDistribLowerLeft->GetX2()) < padRDistribLowerLeft->YtoPixel(padRDistribLowerLeft->GetY1())){
        textsizeLabelsDown = (Double_t)textSizeLabels/padRDistribLowerLeft->XtoPixel(padRDistribLowerLeft->GetX2()) ;
        textsizeFacDown = (Double_t)1./padRDistribLowerLeft->XtoPixel(padRDistribLowerLeft->GetX2()) ;
    } else {
        textsizeLabelsDown = (Double_t)textSizeLabels/padRDistribLowerLeft->YtoPixel(padRDistribLowerLeft->GetY1());
        textsizeFacDown = (Double_t)1./padRDistribLowerLeft->YtoPixel(padRDistribLowerLeft->GetY1());
    }

    TH2F *histoDummy4PanelsUp = new TH2F("histoDummy4PanelsUp","histoDummy4PanelsUp",1000,0.,184.,1000,1.5e-9,1.);
    SetStyleHistoTH2ForGraphs(histoDummy4PanelsUp, "R (cm)","Counts", 0.85*textsizeLabelsUp, textsizeLabelsUp,0.85*textsizeLabelsUp,textsizeLabelsUp, 1,0.9);
    histoDummy4PanelsUp->GetYaxis()->SetRangeUser(1.5e-9,2e-3);

    TH2F * histoDummy4PanelsDown = new TH2F("histoDummy4PanelsDown","histoDummy4PanelsDown",1000,0.,184.,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histoDummy4PanelsDown, "R (cm)","#frac{Data}{MC} ", 0.85*textsizeLabelsDown, textsizeLabelsDown,0.85*textsizeLabelsDown,textsizeLabelsDown, 1,0.9);
    histoDummy4PanelsDown->GetYaxis()->SetRangeUser(0.35,1.65);

    padRDistribUpperLeft->cd();
    padRDistribUpperLeft->SetLogy();
    histoDummy4PanelsUp->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++)
            DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,0.1,kGray);

        DrawGammaSetMarker(histoMappingRData, 20, 1.5, colorData, colorData);
        histoMappingRData->Draw("same,hist");
		DrawGammaSetMarker(histoMappingRMC, 20, 1.5, colorMC, colorMC);
        histoMappingRMC->Draw("same,hist");
		DrawGammaSetMarker(histoMappingRTrueMC, 20, 1.5, colorTrueMC, colorTrueMC);
        histoMappingRTrueMC->Draw("same,hist");
		DrawGammaSetMarker(histoMappingRTrueCombMC, 20, 1.5, colorTrueCombMC, colorTrueCombMC);
        histoMappingRTrueCombMC->DrawCopy("same,hist");
		DrawGammaSetMarker(histoMappingRTruePi0DalMC, 20, 1.5, kBlue-9, kBlue-9);
        histoMappingRTruePi0DalMC->SetFillColor(kBlue-9);
        histoMappingRTruePi0DalMC->SetFillStyle(3244);
        histoMappingRTruePi0DalMC->DrawCopy("same,hist");
		DrawGammaSetMarker(histoMappingRTrueEtaDalMC, 20, 1.5, kBlue-3, kBlue-3);
        histoMappingRTrueEtaDalMC->SetFillColor(kBlue-3);
        histoMappingRTrueEtaDalMC->SetFillStyle(3002);
        histoMappingRTrueEtaDalMC->DrawCopy("same,hist");

        TLegend* legenRdistrib = GetAndSetLegend2(0.63, 0.9-(6*0.9*textsizeLabelsUp), 0.9, 0.9, textSizeLabels);
        legenRdistrib->AddEntry(histoMappingRData,"Data","l");
        legenRdistrib->AddEntry(histoMappingRMC,"MC","l");
        legenRdistrib->AddEntry(histoMappingRTrueMC,"True MC","l");
        legenRdistrib->AddEntry(histoMappingRTrueCombMC,"True comb MC","l");
        legenRdistrib->AddEntry(histoMappingRTruePi0DalMC,"True MC #pi^{0} Dal.","lf");
        legenRdistrib->AddEntry(histoMappingRTrueEtaDalMC,"True MC #eta Dal.","fl");
        legenRdistrib->Draw();

    histoDummy4PanelsUp->Draw("axis,same");
    padRDistribLowerLeft->cd();
    histoDummy4PanelsDown->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++)
            DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,0.1,kGray);
        DrawGammaLines(0.,180,1., 1.,0.1,kGray);

        DrawGammaSetMarker(histoMappingDataMCRatioR, 20, 1.5, colorData, colorData);
        histoMappingDataMCRatioR->Draw("same,histo");

    histoDummy4PanelsDown->Draw("axis,same");
    padRDistribUpperRight->cd();
    padRDistribUpperRight->SetLogy();
    histoDummy4PanelsUp->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++)
            DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,0.1,kGray);

        histoMappingRData->Draw("same,histo");
        histoMappingRMC->Draw("same,histo");
        histoMappingRTrueMC->Draw("same,histo");

		DrawGammaSetMarker(histoMappingMidPtRData, 20, 1.5, colorData, colorData);
        histoMappingMidPtRData->SetLineStyle(2);
        histoMappingMidPtRData->DrawCopy("same,histo");

		DrawGammaSetMarker(histoMappingMidPtRMC, 20, 1.5, colorMC, colorMC);
        histoMappingMidPtRMC->SetLineStyle(2);
        histoMappingMidPtRMC->DrawCopy("same,histo");

		DrawGammaSetMarker(histoMappingMidPtRTrueMC, 20, 1.5, colorTrueMC, colorTrueMC);
        histoMappingMidPtRTrueMC->SetLineStyle(2);
        histoMappingMidPtRTrueMC->DrawCopy("same,histo");

        TLegend* legenRdistribPtcuts = GetAndSetLegend2(0.05, 0.35-(6.2*0.9*textsizeLabelsUp), 0.35, 0.35, textSizeLabels);
        legenRdistribPtcuts->AddEntry(histoMappingRData,"Data","l");
        legenRdistribPtcuts->AddEntry(histoMappingRMC,"MC","l");
        legenRdistribPtcuts->AddEntry(histoMappingRTrueMC,"True MC","l");
        legenRdistribPtcuts->AddEntry(histoMappingMidPtRData,"Data, 0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","l");
        legenRdistribPtcuts->AddEntry(histoMappingMidPtRMC,"MC, 0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","lf");
        legenRdistribPtcuts->AddEntry(histoMappingMidPtRTrueMC,"True MC, 0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","fl");
        legenRdistribPtcuts->Draw();

    histoDummy4PanelsUp->Draw("axis,same");
    padRDistribLowerRight->cd();
    histoDummy4PanelsDown->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++)
            DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,0.1,kGray);
        DrawGammaLines(0.,180,1., 1.,0.1,kGray);

        histoMappingDataMCRatioR->Draw("same,histo");
		DrawGammaSetMarker(histoMappingMidPtDataMCRatioR, 20, 1.5, colorComparisonMC, colorComparisonMC);
        histoMappingMidPtDataMCRatioR->SetLineStyle(2);
        histoMappingMidPtDataMCRatioR->Draw("same,histo");

        TLegend* legenRdistribPtcutsRatio = GetAndSetLegend2(0.05, 0.27-(2*0.9*textsizeLabelsDown), 0.35, 0.27, textSizeLabels);
        legenRdistribPtcutsRatio->AddEntry(histoMappingDataMCRatioR,"No #it{p}_{T} cut","l");
        legenRdistribPtcutsRatio->AddEntry(histoMappingMidPtDataMCRatioR,"0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","l");
        legenRdistribPtcutsRatio->Draw();

    histoDummy4PanelsDown->Draw("axis,same");
    canvasPhotonR->Print(Form("%s/PhotonConvR%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //______________________________________ Purity _________________________________
    Double_t arrayXpurity[2];
    Double_t arrayYpurity[3];
    Double_t relXpurity[3];
    Double_t relYpurity[3];
    textSizeLabels = 30;
    Double_t textsizeLabelsDownpurity = 0;
    Double_t textsizeFacDownpurity = 0;
    Double_t textsizeLabelsUppurity = 0;
    Double_t textsizeFacUppurity = 0;
    ReturnCorrectValuesForCanvasScaling(1200, 1200, 1, 2, 0.08, 0.025, 0.025, 0.08, arrayXpurity, arrayYpurity, relXpurity, relYpurity);
    TCanvas * canvasPurity = new TCanvas("canvasPurity","",1200,1200);

    TPad* padUpperPurity = new TPad("","",arrayXpurity[0], arrayYpurity[1], arrayXpurity[1], arrayYpurity[0],-1, -1, -2);
    DrawGammaPadSettings( padUpperPurity, relXpurity[0], relXpurity[2], relYpurity[0], relYpurity[1]-0.2);
    TPad* padLowerPurity = new TPad("","",arrayXpurity[0], arrayYpurity[2], arrayXpurity[1], arrayYpurity[1],-1, -1, -2);
    DrawGammaPadSettings( padLowerPurity, relXpurity[0], relXpurity[2], relYpurity[1]-0.2, relYpurity[2]);
    marginXRatio = relXpurity[0]*1200;
        if (padUpperPurity->XtoPixel(padUpperPurity->GetX2()) < padUpperPurity->YtoPixel(padUpperPurity->GetY1())){
            textsizeLabelsUppurity = (Double_t)textSizeLabels/padUpperPurity->XtoPixel(padUpperPurity->GetX2()) ;
            textsizeFacUppurity = (Double_t)1./padUpperPurity->XtoPixel(padUpperPurity->GetX2()) ;
        } else {
            textsizeLabelsUppurity = (Double_t)textSizeLabels/padUpperPurity->YtoPixel(padUpperPurity->GetY1());
            textsizeFacUppurity = (Double_t)1./padUpperPurity->YtoPixel(padUpperPurity->GetY1());
        }
        padUpperPurity->Draw();
        if (padLowerPurity->XtoPixel(padLowerPurity->GetX2()) < padLowerPurity->YtoPixel(padLowerPurity->GetY1())){
            textsizeLabelsDownpurity = (Double_t)textSizeLabels/padLowerPurity->XtoPixel(padLowerPurity->GetX2()) ;
            textsizeFacDownpurity = (Double_t)1./padLowerPurity->XtoPixel(padLowerPurity->GetX2()) ;
        } else {
            textsizeLabelsDownpurity = (Double_t)textSizeLabels/padLowerPurity->YtoPixel(padLowerPurity->GetY1());
            textsizeFacDownpurity = (Double_t)1./padLowerPurity->YtoPixel(padLowerPurity->GetY1());
        }
        padLowerPurity->Draw();

        TH2F *histoDummyPurityR = new TH2F("histoDummyPurityR","histoDummyPurityR",1000,0.,180.,1000,0.,1.2);
        SetStyleHistoTH2ForGraphs(histoDummyPurityR, "R (cm)","Purity", 0.9*textsizeLabelsUppurity, textsizeLabelsUppurity,0.9*textsizeLabelsUppurity,textsizeLabelsUppurity, 0.9,0.1/(textsizeFacUppurity*marginXRatio));
//         histoDummyPurityR->GetYaxis()->SetLabelOffset(+0.01);
        TH2F *histoDummyPurityPt =  new TH2F("histoDummyPurityPt","histoDummyPurityPt",1000,0.,8.,1000,0.,1.2);
        SetStyleHistoTH2ForGraphs(histoDummyPurityPt, "p_{T} (GeV/c)","Purity", 0.9*textsizeLabelsDownpurity, textsizeLabelsDownpurity,0.9*textsizeLabelsDownpurity,textsizeLabelsDownpurity, 0.9,0.1/(textsizeFacDownpurity*marginXRatio));
//         histoDummyPurityPt->GetYaxis()->SetLabelOffset(+0.01);
//         histoDummyPurityPt->GetYaxis()->SetRangeUser(0.3,1.55);

    padUpperPurity->cd();
    histoDummyPurityR->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++){
            DrawGammaLines(arrayRBins[iR],arrayRBins[iR],0.,1.2,0.1,kGray,2);
        }
        DrawGammaLines(0.,180,1., 1.,0.1,kGray+1);
        DrawGammaSetMarker(histoMappingPurityR, 20, 1.5, colorData, colorData);

        histoMappingPurityR->Draw("same,hist");

        TLegend* legenPurity = GetAndSetLegend2(0.15, 0.9-(1*0.9*textsizeLabelsUppurity), 0.7, 0.9, textSizeLabels);
        legenPurity->AddEntry((TObject*)0,Form("Cut: %s",fCutSelectionRead.Data()),"");
        legenPurity->Draw();

    histoDummyPurityR->Draw("same,axis");

    padLowerPurity->cd();
    histoDummyPurityPt->DrawCopy();

        DrawGammaLines(0.,8.,1., 1.,0.1,kGray+1);
        DrawGammaSetMarker(histoMappingPurityPt, 20, 1.5, colorData, colorData);
        histoMappingPurityPt->Draw("same,hist");
        DrawGammaSetMarker(histoMappingPurityPt5cm, 20, 1.5, kGreen+2, kGreen+2);
        histoMappingPurityPt5cm->Draw("same,hist");

        TLegend* legenPurityPt = GetAndSetLegend2(0.75, 0.55-(2*0.9*textsizeLabelsUppurity), 0.9, 0.55, textSizeLabels);
        legenPurityPt->AddEntry(histoMappingPurityPt,"no R cut","l");
        legenPurityPt->AddEntry(histoMappingPurityPt5cm,"R > 5 cm","l");
        legenPurityPt->Draw();
        legenPurityPt->Draw();

    histoDummyPurityPt->Draw("same,axis");
    canvasPurity->Print(Form("%s/PhotonPurity%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //______________________________________ Photon characteristics _________________________________
    Double_t arrayXPhotonChar[3];
    Double_t arrayYPhotonChar[3];
    Double_t relXPhotonChar[3];
    Double_t relYPhotonChar[3];
    textSizeLabels = 50;
    textsizeLabelsUp = 0;
    textsizeFacUp = 0;
    textsizeLabelsDown = 0;
    textsizeFacDown = 0;
    ReturnCorrectValuesForCanvasScaling(2500,2000, 2, 2,0.06, 0.025, 0.025,0.06,arrayXPhotonChar,arrayYPhotonChar,relXPhotonChar,relYPhotonChar);
    TCanvas* canvasPhotonChar = new TCanvas("canvasPhotonChar","",0,0,2500,2000);

    TPad* padPhotonCharLowerLeft = new TPad("padPhotonCharLowerLeft", "", arrayXPhotonChar[0], arrayYPhotonChar[2], arrayXPhotonChar[1], arrayYPhotonChar[1],-1, -1, -2);
    DrawGammaPadSettings( padPhotonCharLowerLeft, relXPhotonChar[0], relXPhotonChar[1]-0.001, relYPhotonChar[1]-0.005, relYPhotonChar[2]);
    padPhotonCharLowerLeft->Draw();

    TPad* padPhotonCharUpperLeft = new TPad("padPhotonCharUpperLeft", "", arrayXPhotonChar[0], arrayYPhotonChar[1], arrayXPhotonChar[1], arrayYPhotonChar[0],-1, -1, -2);
    DrawGammaPadSettings( padPhotonCharUpperLeft, relXPhotonChar[0], relXPhotonChar[1]-0.001, relYPhotonChar[0], relYPhotonChar[1]-0.005);
    padPhotonCharUpperLeft->Draw();

    TPad* padPhotonCharLowerRight = new TPad("padPhotonCharLowerRight", "", arrayXPhotonChar[1], arrayYPhotonChar[2], arrayXPhotonChar[2], arrayYPhotonChar[1],-1, -1, -2);
    DrawGammaPadSettings( padPhotonCharLowerRight, relXPhotonChar[1]-0.001, relXPhotonChar[2], relYPhotonChar[1]-0.005, relYPhotonChar[2]);
    padPhotonCharLowerRight->Draw();

    TPad* padPhotonCharUpperRight = new TPad("padPhotonCharUpperRight", "", arrayXPhotonChar[1], arrayYPhotonChar[1], arrayXPhotonChar[2], arrayYPhotonChar[0],-1, -1, -2);
    DrawGammaPadSettings( padPhotonCharUpperRight, relXPhotonChar[1]-0.001, relXPhotonChar[2], relYPhotonChar[0], relYPhotonChar[1]-0.005);
    padPhotonCharUpperRight->Draw();

    if (padPhotonCharUpperLeft->XtoPixel(padPhotonCharUpperLeft->GetX2()) < padPhotonCharUpperLeft->YtoPixel(padPhotonCharUpperLeft->GetY1())){
        textsizeLabelsUp = (Double_t)textSizeLabels/padPhotonCharUpperLeft->XtoPixel(padPhotonCharUpperLeft->GetX2()) ;
        textsizeFacUp = (Double_t)1./padPhotonCharUpperLeft->XtoPixel(padPhotonCharUpperLeft->GetX2()) ;
    } else {
        textsizeLabelsUp = (Double_t)textSizeLabels/padPhotonCharUpperLeft->YtoPixel(padPhotonCharUpperLeft->GetY1());
        textsizeFacUp = (Double_t)1./padPhotonCharUpperLeft->YtoPixel(padPhotonCharUpperLeft->GetY1());
    }
    if (padPhotonCharLowerLeft->XtoPixel(padPhotonCharLowerLeft->GetX2()) < padPhotonCharLowerLeft->YtoPixel(padPhotonCharLowerLeft->GetY1())){
        textsizeLabelsDown = (Double_t)textSizeLabels/padPhotonCharLowerLeft->XtoPixel(padPhotonCharLowerLeft->GetX2()) ;
        textsizeFacDown = (Double_t)1./padPhotonCharLowerLeft->XtoPixel(padPhotonCharLowerLeft->GetX2()) ;
    } else {
        textsizeLabelsDown = (Double_t)textSizeLabels/padPhotonCharLowerLeft->YtoPixel(padPhotonCharLowerLeft->GetY1());
        textsizeFacDown = (Double_t)1./padPhotonCharLowerLeft->YtoPixel(padPhotonCharLowerLeft->GetY1());
    }

    padPhotonCharUpperLeft->cd();
    padPhotonCharUpperLeft->SetLogy();
    TH2F *histoDummyDCA = new TH2F("histoDummyDCA","histoDummyDCA",1000,0.,2.,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyDCA, "DCA (cm)","Counts", 0.85*textsizeLabelsUp, textsizeLabelsUp,0.85*textsizeLabelsUp,textsizeLabelsUp, 0.9,0.9);
    histoDummyDCA->GetYaxis()->SetRangeUser(1.e-6,2e-3);
    histoDummyDCA->DrawCopy();

        DrawGammaSetMarker(histoMappingDCAData, 20, 1.5, colorData, colorData);
        histoMappingDCAData->Draw("same,hist");
		DrawGammaSetMarker(histoMappingDCAMC, 20, 1.5, colorMC, colorMC);
        histoMappingDCAMC->Draw("same,hist");

    histoDummyDCA->Draw("axis,same");
    padPhotonCharLowerLeft->cd();
    padPhotonCharLowerLeft->SetLogy();
    TH2F * histoDummyMass = new TH2F("histoDummyMass","histoDummyMass",1000,0.,0.2,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyMass, "Mass (GeV)","Counts ", 0.85*textsizeLabelsDown, textsizeLabelsDown,0.85*textsizeLabelsDown,textsizeLabelsDown, 0.9,0.92);
    histoDummyMass->GetYaxis()->SetRangeUser(1.e-6,2e-3);
    histoDummyMass->DrawCopy();

        DrawGammaSetMarker(histoMappingMassData, 20, 1.5, colorData, colorData);
        histoMappingMassData->Draw("same,hist");
		DrawGammaSetMarker(histoMappingMassMC, 20, 1.5, colorMC, colorMC);
        histoMappingMassMC->Draw("same,hist");

    histoDummyMass->Draw("axis,same");
    padPhotonCharUpperRight->cd();
    padPhotonCharUpperRight->SetLogy();
    TH2F *histoDummyPsiPair = new TH2F("histoDummyPsiPair","histoDummyPsiPair",1000,0.,0.3,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyPsiPair, "#Psi_{pair}","Counts", 0.85*textsizeLabelsUp, textsizeLabelsUp,0.85*textsizeLabelsUp,textsizeLabelsUp, 0.9,0.91);
    histoDummyPsiPair->GetYaxis()->SetRangeUser(1.e-6,2.e-2);
    histoDummyPsiPair->DrawCopy();

        DrawGammaSetMarker(histoMappingPsiPairData, 20, 1.5, colorData, colorData);
        histoMappingPsiPairData->Draw("same,hist");
		DrawGammaSetMarker(histoMappingPsiPairMC, 20, 1.5, colorMC, colorMC);
        histoMappingPsiPairMC->Draw("same,hist");

    histoDummyPsiPair->Draw("axis,same");
    padPhotonCharLowerRight->cd();
    padPhotonCharLowerRight->SetLogy();
    TH2F * histoDummyChi2 = new TH2F("histoDummyChi2","histoDummyChi2",1000,0.,40.,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyChi2, "#chi^{2}","Counts ", 0.85*textsizeLabelsDown, textsizeLabelsDown,0.85*textsizeLabelsDown,textsizeLabelsDown, 0.9,0.92);
    histoDummyChi2->GetYaxis()->SetRangeUser(1.e-6,2.e-2);
    histoDummyChi2->DrawCopy();

        DrawGammaSetMarker(histoMappingChi2Data, 20, 1.5, colorData, colorData);
        histoMappingChi2Data->Draw("same,hist");
		DrawGammaSetMarker(histoMappingChi2MC, 20, 1.5, colorMC, colorMC);
        histoMappingChi2MC->Draw("same,hist");

    histoDummyChi2->Draw("axis,same");
    canvasPhotonChar->Print(Form("%s/PhotonChar%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    Double_t arrayXPhiInRbins[2];
    Double_t arrayYPhiInRbins[3];
    Double_t relXPhiInRbins[3];
    Double_t relYPhiInRbins[3];
    textSizeLabels = 30;
    textsizeLabelsDown = 0;
    textsizeFacDown = 0;
    textsizeLabelsUp = 0;
    textsizeFacUp = 0;
    ReturnCorrectValuesForCanvasScaling(1200, 1000, 1, 2, 0.1, 0.025, 0.025, 0.08, arrayX, arrayY, relX, relY);
    TCanvas * canvasPhiInRbins = new TCanvas("canvasPhiInRbins","",1200,1000);

    TPad* padPhiInRbinsUpper = new TPad("","",arrayX[0], arrayY[1], arrayX[1], arrayY[0],-1, -1, -2);
    DrawGammaPadSettings( padPhiInRbinsUpper, relX[0], relX[2], relY[0], relY[1]);
    TPad* padPhiInRbinsLower = new TPad("","",arrayX[0], arrayY[2], arrayX[1], arrayY[1],-1, -1, -2);
    DrawGammaPadSettings( padPhiInRbinsLower, relX[0], relX[2], relY[1], relY[2]);
    marginXRatio = relX[0]*1200;
        padPhiInRbinsUpper->Draw();
        if (padPhiInRbinsUpper->XtoPixel(padPhiInRbinsUpper->GetX2()) < padPhiInRbinsUpper->YtoPixel(padPhiInRbinsUpper->GetY1())){
            textsizeLabelsUp = (Double_t)textSizeLabels/padPhiInRbinsUpper->XtoPixel(padPhiInRbinsUpper->GetX2()) ;
            textsizeFacUp = (Double_t)1./padPhiInRbinsUpper->XtoPixel(padPhiInRbinsUpper->GetX2()) ;
        } else {
            textsizeLabelsUp = (Double_t)textSizeLabels/padPhiInRbinsUpper->YtoPixel(padPhiInRbinsUpper->GetY1());
            textsizeFacUp = (Double_t)1./padPhiInRbinsUpper->YtoPixel(padPhiInRbinsUpper->GetY1());
        }
        padPhiInRbinsLower->Draw();
        if (padPhiInRbinsLower->XtoPixel(padPhiInRbinsLower->GetX2()) < padPhiInRbinsLower->YtoPixel(padPhiInRbinsLower->GetY1())){
            textsizeLabelsDown = (Double_t)textSizeLabels/padPhiInRbinsLower->XtoPixel(padPhiInRbinsLower->GetX2()) ;
            textsizeFacDown = (Double_t)1./padPhiInRbinsLower->XtoPixel(padPhiInRbinsLower->GetX2()) ;
        } else {
            textsizeLabelsDown = (Double_t)textSizeLabels/padPhiInRbinsLower->YtoPixel(padPhiInRbinsLower->GetY1());
            textsizeFacDown = (Double_t)1./padPhiInRbinsLower->YtoPixel(padPhiInRbinsLower->GetY1());
        }

    TH2F *histoDummyPhiInRbinsUp = new TH2F("histoDummyPhiInRbinsUp","histoDummyPhiInRbinsUp",1000,0.,2*TMath::Pi(),10000,2e-7,.1);
    SetStyleHistoTH2ForGraphs(histoDummyPhiInRbinsUp, "#varphi (rad)","Counts", 0.9*textsizeLabelsUp, textsizeLabelsUp,0.9*textsizeLabelsUp,textsizeLabelsUp, 1,0.15/(textsizeFacUp*marginXRatio));

    TH2F *histoDummyPhiInRbinsDown =  new TH2F("histoDummyPhiInRbinsDown","histoDummyPhiInRbinsDown",1000,0.,2*TMath::Pi(),1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histoDummyPhiInRbinsDown, "#varphi (rad)","#frac{Data}{MC} ", 0.9*textsizeLabelsDown, textsizeLabelsDown,0.9*textsizeLabelsDown,textsizeLabelsDown, 1,0.15/(textsizeFacDown*marginXRatio));
    histoDummyPhiInRbinsDown->GetYaxis()->SetRangeUser(0.3,1.55);

    for(Int_t iR = 0; iR < nBinsR; iR++){

        canvasPhiInRbins->cd();
        padPhiInRbinsUpper->cd();
        padPhiInRbinsUpper->SetLogy();
        histoDummyPhiInRbinsUp->GetYaxis()->SetRangeUser(histoMappingPhiInRData[iR]->GetMinimum(),histoMappingPhiInRData[iR]->GetMaximum()*10);
        histoDummyPhiInRbinsUp->DrawCopy();

            DrawGammaSetMarker(histoMappingPhiInRData[iR], 20, 1.5, colorData, colorData);
            histoMappingPhiInRData[iR]->Draw("same,hist");
            DrawGammaSetMarker(histoMappingPhiInRMC[iR], 20, 1.5, colorMC, colorMC);
            histoMappingPhiInRMC[iR]->Draw("same,hist");

            TLatex *latexBinning = new TLatex(0.15,0.85,Form("#varphi in range %s (%s)",arrayRangesRBins[iR].Data(), arrayNamesRBins[iR].Data()));
            SetStyleTLatex( latexBinning, sizeTextNameBins,2);
            latexBinning->Draw();

            TLegend* legenPhiInRbins = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
            legenPhiInRbins->AddEntry(histoMappingPhiInRData[iR],"Data","l");
            legenPhiInRbins->AddEntry(histoMappingPhiInRMC[iR],"MC","l");
            legenPhiInRbins->Draw();

        histoDummyPhiInRbinsUp->Draw("same,axis");
        padPhiInRbinsLower->cd();
        histoDummyPhiInRbinsDown->DrawCopy();

            DrawGammaLines(0.,2*TMath::Pi(),1., 1.,0.1,kGray+1);

            DrawGammaSetMarker(histoMappingDataMCRatioPhiInR[iR], 20, 1.5, colorData, colorData);
            histoMappingDataMCRatioPhiInR[iR]->Draw("same,hist");

        histoDummyPhiInRbinsDown->Draw("same,axis");
        canvasPhiInRbins->Update();
        canvasPhiInRbins->SaveAs(Form("%s/PhiInRbins%s_%i_%s.%s",outputDirectory.Data(),optionPeriod.Data(),iR,fCutSelectionRead.Data(),suffix.Data()));

    }
    delete canvasPhiInRbins;

    Double_t arrayXZInRbins[2];
    Double_t arrayYZInRbins[3];
    Double_t relXZInRbins[3];
    Double_t relYZInRbins[3];
    textSizeLabels = 30;
    textsizeLabelsDown = 0;
    textsizeFacDown = 0;
    textsizeLabelsUp = 0;
    textsizeFacUp = 0;
    ReturnCorrectValuesForCanvasScaling(1200, 1000, 1, 2, 0.1, 0.025, 0.03, 0.08, arrayX, arrayY, relX, relY);
    TCanvas * canvasZInRbins = new TCanvas("canvasZInRbins","",1200,1000);

    TPad* padZInRbinsUpper = new TPad("","",arrayX[0], arrayY[1], arrayX[1], arrayY[0],-1, -1, -2);
    DrawGammaPadSettings( padZInRbinsUpper, relX[0], relX[2], relY[0], relY[1]);
    TPad* padZInRbinsLower = new TPad("","",arrayX[0], arrayY[2], arrayX[1], arrayY[1],-1, -1, -2);
    DrawGammaPadSettings( padZInRbinsLower, relX[0], relX[2], relY[1], relY[2]);
    marginXRatio = relX[0]*1200;
        padZInRbinsUpper->Draw();
        if (padZInRbinsUpper->XtoPixel(padZInRbinsUpper->GetX2()) < padZInRbinsUpper->YtoPixel(padZInRbinsUpper->GetY1())){
            textsizeLabelsUp = (Double_t)textSizeLabels/padZInRbinsUpper->XtoPixel(padZInRbinsUpper->GetX2()) ;
            textsizeFacUp = (Double_t)1./padZInRbinsUpper->XtoPixel(padZInRbinsUpper->GetX2()) ;
        } else {
            textsizeLabelsUp = (Double_t)textSizeLabels/padZInRbinsUpper->YtoPixel(padZInRbinsUpper->GetY1());
            textsizeFacUp = (Double_t)1./padZInRbinsUpper->YtoPixel(padZInRbinsUpper->GetY1());
        }
        padZInRbinsLower->Draw();
        if (padZInRbinsLower->XtoPixel(padZInRbinsLower->GetX2()) < padZInRbinsLower->YtoPixel(padZInRbinsLower->GetY1())){
            textsizeLabelsDown = (Double_t)textSizeLabels/padZInRbinsLower->XtoPixel(padZInRbinsLower->GetX2()) ;
            textsizeFacDown = (Double_t)1./padZInRbinsLower->XtoPixel(padZInRbinsLower->GetX2()) ;
        } else {
            textsizeLabelsDown = (Double_t)textSizeLabels/padZInRbinsLower->YtoPixel(padZInRbinsLower->GetY1());
            textsizeFacDown = (Double_t)1./padZInRbinsLower->YtoPixel(padZInRbinsLower->GetY1());
        }

    TH2F *histoDummyZInRbinsUp = new TH2F("histoDummyZInRbinsUp","histoDummyZInRbinsUp",1000,-160,160,10000,-0.001,.1);
    SetStyleHistoTH2ForGraphs(histoDummyZInRbinsUp, "Z (cm)","Counts", 0.9*textsizeLabelsUp, textsizeLabelsUp,0.9*textsizeLabelsUp,textsizeLabelsUp, 1,0.15/(textsizeFacUp*marginXRatio));

    TH2F *histoDummyZInRbinsDown =  new TH2F("histoDummyZInRbinsDown","histoDummyZInRbinsDown",1000,-160,160,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histoDummyZInRbinsDown, "Z (cm)","#frac{Data}{MC} ", 0.9*textsizeLabelsDown, textsizeLabelsDown,0.9*textsizeLabelsDown,textsizeLabelsDown, 1,0.15/(textsizeFacDown*marginXRatio));
    histoDummyZInRbinsDown->GetYaxis()->SetRangeUser(0.3,1.55);

    for(Int_t iR = 0; iR < nBinsR; iR++){

        canvasZInRbins->cd();
        padZInRbinsUpper->cd();
//         histoDummyZInRbinsUp->GetYaxis()->SetRangeUser(histoMappingZInRData[iR]->GetMinimum()*0.1,histoMappingZInRData[iR]->GetMaximum()*10);
        histoDummyZInRbinsUp->GetYaxis()->SetRangeUser(0.,histoMappingZInRData[iR]->GetMaximum()*1.5);
        histoDummyZInRbinsUp->DrawCopy();

            DrawGammaLines(-10,-10, histoMappingZInRData[iR]->GetMinimum()*0.1,histoMappingZInRData[iR]->GetMaximum()*1.5,0.1,kGray+1,2);
            DrawGammaLines(10,10, histoMappingZInRData[iR]->GetMinimum()*0.1,histoMappingZInRData[iR]->GetMaximum()*1.5,0.1,kGray+1,2);

            DrawGammaSetMarker(histoMappingZInRData[iR], 20, 1.5, colorData, colorData);
            histoMappingZInRData[iR]->Draw("same,hist");
            DrawGammaSetMarker(histoMappingZInRMC[iR], 20, 1.5, colorMC, colorMC);
            histoMappingZInRMC[iR]->Draw("same,hist");

            TLatex *latexBinning = new TLatex(0.15,0.85,Form("Z in range %s (%s)",arrayRangesRBins[iR].Data(), arrayNamesRBins[iR].Data()));
            SetStyleTLatex( latexBinning, sizeTextNameBins,2);
            latexBinning->Draw();

            TLegend* legenZInRbins = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
            legenZInRbins->AddEntry(histoMappingZInRData[iR],"Data","l");
            legenZInRbins->AddEntry(histoMappingZInRMC[iR],"MC","l");
            legenZInRbins->Draw();

        histoDummyZInRbinsUp->Draw("same,axis");
        padZInRbinsLower->cd();
        histoDummyZInRbinsDown->DrawCopy();

            DrawGammaLines(-10,-10, 0.3,1.55,0.1,kGray+1,2);
            DrawGammaLines(10,10,0.3,1.55,0.1,kGray+1,2);
            DrawGammaLines(-20,20,1., 1.,0.1,kGray+1);

            DrawGammaSetMarker(histoMappingDataMCRatioZInR[iR], 20, 1.5, colorData, colorData);
            histoMappingDataMCRatioZInR[iR]->Draw("same,hist");

        histoDummyZInRbinsDown->Draw("same,axis");
        canvasZInRbins->Update();
        canvasZInRbins->SaveAs(Form("%s/ZInRbins%s_%i_%s.%s",outputDirectory.Data(),optionPeriod.Data(),iR,fCutSelectionRead.Data(),suffix.Data()));

    }
    delete canvasZInRbins;


    TProfile* fProfileContainingMaterialBudgetWeights = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins","profileContainingMaterialBudgetWeights_manyRadialBins",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeights->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeights->GetYaxis()->SetTitle("weight_per_gamma");

    for(Int_t i=0; i<nBinsR; i++){
        cout<< arrayRBins[i] << " " << histoMappingDataMCRatioR->GetBinContent(i+1) <<endl;
        fProfileContainingMaterialBudgetWeights->Fill(arrayRBins[i],histoMappingDataMCRatioR->GetBinContent(i+1));
    }

    cout<< " pT 300"<< endl;
    TProfile* fProfileContainingMaterialBudgetWeightsFullPt = new TProfile("profileContainingMaterialBudgetWeightsFullPt_manyRadialBins","profileContainingMaterialBudgetWeightsFullPt_manyRadialBins",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsFullPt->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsFullPt->GetYaxis()->SetTitle("weight_per_gamma");

    for(Int_t i=0; i<nBinsR; i++){
        cout << "R bin: " << arrayRBins[i] << " ratio Data/MC scaled to gas " << histoMappingDataMCRatioRFullPtScaledToGas->GetBinContent(i+1) <<endl;
        fProfileContainingMaterialBudgetWeightsFullPt->Fill(arrayRBins[i],histoMappingDataMCRatioRFullPtScaledToGas->GetBinContent(i+1));
    }
    //Form("%s/PhotonChar_%s.%s",optionPeriod.Data(),fCutSelectionRead.Data() ))
    cout<< Form("MCInputFileMaterialBudgetWeights%s_%s.root",optionPeriod.Data(),fCutSelectionRead.Data())<< endl;
    TFile outFile(Form("MCInputFileMaterialBudgetWeights%s_%s.root",optionPeriod.Data(),fCutSelectionRead.Data()) ,"RECREATE");

        fProfileContainingMaterialBudgetWeights->Write();
        fProfileContainingMaterialBudgetWeightsFullPt->Write();

        histoMappingRData->Write("Data");
        histoMappingRMC->Write("MC");
        histoMappingMidPtRData->Write("DataMidtPt");
        histoMappingMidPtRMC->Write("MCMidPt");
        histoMappingRDataScaledToGas->Write("DataScaledToGas");
        histoMappingRMCScaledToGas->Write("MCScaledToGas");
    outFile.Close();


}
