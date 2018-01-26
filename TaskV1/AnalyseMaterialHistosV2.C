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
#include "AnalyseMaterialHistosV2.h"


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

    //_______________________ Data _______________________
    TFile fData(fileName.Data());
    TString nameMainDir = "GammaConvMaterial";
    TList *TopDir =(TList*)fData.Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversion == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
        return;
    }
    TList *ESDContainer           = (TList*)HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));

    TH1F* histoEventQualityData   = (TH1F*)ESDContainer->FindObject("NEvents");
    TH1F *histoGoodESDTracksData  = (TH1F*)ESDContainer->FindObject("GoodESDTracks09");
    TH2F *histoRPhiData           = (TH2F*)ESDContainer->FindObject("ESD_Conversion_RPhi");
    TH2F *histoREtaData           = (TH2F*)ESDContainer->FindObject("ESD_Conversion_REta");
    TH2F *histoRZData             = (TH2F*)ESDContainer->FindObject("ESD_Conversion_RZ");
    TH2F *histoRPtData            = (TH2F*)ESDContainer->FindObject("ESD_Conversion_RPt");
    TH1F *histoDCAData            = (TH1F*)ESDContainer->FindObject("ESD_Conversion_DCA");
    TH1F *histoPsiPairData        = (TH1F*)ESDContainer->FindObject("ESD_Conversion_PsiPair");
    TH1F *histoChi2Data           = (TH1F*)ESDContainer->FindObject("ESD_Conversion_Chi2");
    TH1F *histoMassData           = (TH1F*)ESDContainer->FindObject("ESD_Conversion_Mass");

    //________________________ MC ________________________
    TFile fMC(fileNameMC.Data());
    TList *TopDirMC =(TList*)fMC.Get(nameMainDir.Data());
    if(TopDirMC == NULL){
        cout<<"ERROR: TopDirMC not Found"<<endl;
        return;
    }
    TList *HistosGammaConversionMC = (TList*)TopDirMC->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversionMC == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in FileMC"<<endl;
        return;
    }
    TList *ESDContainerMC         = (TList*)HistosGammaConversionMC->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));

    TH1F* histoEventQualityMC     = (TH1F*)ESDContainerMC->FindObject("NEvents");
    TH1F *histoGoodESDTracksMC    = (TH1F*)ESDContainerMC->FindObject("GoodESDTracks09");
    TH2F *histoRPhiMC             = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_RPhi");
    TH2F *histoREtaMC             = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_REta");
    TH2F *histoRZMC               = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_RZ");
    TH2F *histoRPtMC              = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_RPt");
    TH1F *histoDCAMC              = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_DCA");
    TH1F *histoPsiPairMC          = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_PsiPair");
    TH1F *histoChi2MC             = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_Chi2");
    TH1F *histoMassMC             = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_Mass");

    TList *MCContainer            = (TList*)HistosGammaConversionMC->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
    TList *TrueMCContainer        = (TList*)HistosGammaConversionMC->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));

    TH2F *histoRPhiTrueMC         = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RPhi");
    TH2F *histoREtaTrueMC         = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_REta");
    TH2F *histoRZTrueMC           = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RZ");
    TH2F *histoRPtTrueMC          = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RPt");
    TH1F *histoDCATrueMC          = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_DCA");
    TH1F *histoPsiPairTrueMC      = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_PsiPair");
    TH1F *histoChi2TrueMC         = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_Chi2");
    TH1F *histoMassTrueMC         = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_Mass");

    TH1F *histoPi0DalRPtTrueMC    = (TH1F*)TrueMCContainer->FindObject("ESD_TruePi0DalConversion_RPt");
    TH1F *histoPi0DalREtaTrueMC   = (TH1F*)TrueMCContainer->FindObject("ESD_TruePi0DalConversion_REta");
    TH1F *histoEtaDalRPtTrueMC    = (TH1F*)TrueMCContainer->FindObject("ESD_TrueEtaDalConversion_RPt");
    TH1F *histoEtaDalREtaTrueMC   = (TH1F*)TrueMCContainer->FindObject("ESD_TrueEtaDalConversion_REta");
    TH1F *histoCombRPtTrueMC      = (TH1F*)TrueMCContainer->FindObject("ESD_TrueCombinatorialConversion_RPt");
    TH1F *histoCombREtaTrueMC     = (TH1F*)TrueMCContainer->FindObject("ESD_TrueCombinatorialConversion_REta");


    //________________________ normalisation _______________________
    //using the charged particles (~pions) to normalise instead of the neutral pions (~tot gammas)

    Float_t numberGoodEventsData  = histoEventQualityData->GetBinContent(1);
    Double_t meanMultData         = histoGoodESDTracksData->GetMean();
    Float_t normFactorReconstData = 1./(numberGoodEventsData*meanMultData);

    Float_t numberGoodEventsMC    = histoEventQualityMC->GetBinContent(1);
    Double_t meanMultMC           = histoGoodESDTracksMC->GetMean();
    Float_t normFactorReconstMC   = 1./numberGoodEventsMC*1./meanMultMC;

    cout << "Normalization factor data " << meanMultData << " and MC " << meanMultMC << " -> Data/MC: " << meanMultData/meanMultMC << endl;

    histoGoodESDTracksData->Scale(1./numberGoodEventsData);
//     GammaScalingHistogramm(histoGoodESDTracksData,normFactorReconstData);
    histoGoodESDTracksMC->Scale(1./numberGoodEventsMC);
//     GammaScalingHistogramm(histoGoodESDTracksMC,normFactorReconstMC);

    GammaScalingHistogramm(histoRPtData,normFactorReconstData);
    GammaScalingHistogramm(histoRPtMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRPtTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPi0DalRPtTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoEtaDalRPtTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoCombRPtTrueMC,normFactorReconstMC);

    GammaScalingHistogramm(histoRPhiData,normFactorReconstData);
    GammaScalingHistogramm(histoRPhiMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRPhiTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRZData,normFactorReconstData);
    GammaScalingHistogramm(histoRZMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRZTrueMC,normFactorReconstMC);

    GammaScalingHistogramm(histoDCAData,normFactorReconstData);
    GammaScalingHistogramm(histoDCAMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPsiPairData,normFactorReconstData);
    GammaScalingHistogramm(histoPsiPairMC,normFactorReconstMC);
    GammaScalingHistogramm(histoChi2Data,normFactorReconstData);
    GammaScalingHistogramm(histoChi2MC,normFactorReconstMC);
    GammaScalingHistogramm(histoMassData,normFactorReconstData);
    GammaScalingHistogramm(histoMassMC,normFactorReconstMC);


    //___________________ Projecting histos ___________________

    cout << "******************** projection data histograms... ******************** " << endl;
    TH1F *histoRData        = (TH1F*)histoRPtData->ProjectionX("histoRData");
    Double_t *projRBins[4]  = {0., 5., 35., 180.};
    TH1F *histoPtinRBinData[4];
    for(Int_t i=0; i<4; i++){
        histoPtinRBinData[i] = (TH1F*)histoRPtData->ProjectionX(Form("histoPtinRBinData_%i",i),
                                                                histoRPtData->GetXaxis()->FindBin(projRBins[i]+0.001),
                                                                histoRPtData->GetXaxis()->FindBin(projRBins[i+1]+0.001),"e");
        cout << "Projecting Pt in R bin " << projRBins[i] << " - " << projRBins[i+1] << endl;
    }
    TH1F *histoRDataRebin = (TH1F*)histoRData->Clone("histoRDataRebin");
    ConvGammaRebinWithBinCorrection(histoRDataRebin,rebinRPlots);

    TH1F *histoPtData       = (TH1F*)histoRPtData->ProjectionY("histoPtData");
    Double_t *projPtBins[5]  = {0.4, 4., 8., 12., 20.};
    TH1F *histoRinPtBinData[5];
    for(Int_t i=0; i<5; i++){
        histoRinPtBinData[i] = (TH1F*)histoRPtData->ProjectionY(Form("histoRinPtBinData_%i",i),
                                                                histoRPtData->GetYaxis()->FindBin(projPtBins[i]+0.001),
                                                                histoRPtData->GetYaxis()->FindBin(projPtBins[i+1]+0.001),"e");
        cout << "Projecting R in Pt bin " << projPtBins[i] << " - " << projPtBins[i+1] << endl;
    }
    TH1F *histoPtDataRebin = (TH1F*)histoPtData->Clone("histoPtDataRebin");
    ConvGammaRebinWithBinCorrection(histoPtDataRebin,rebinPtPlots);


    cout << "******************** projection MC histograms... ******************** " << endl;
    TH1F *histoRMC        = (TH1F*)histoRPtMC->ProjectionX("histoRMC");
    Double_t *projRBins[4]  = {0., 5., 35., 180.};
    TH1F *histoPtinRBinMC[4];
    for(Int_t i=0; i<4; i++){
        histoPtinRBinMC[i] = (TH1F*)histoRPtMC->ProjectionX(Form("histoPtinRBinMC_%i",i),
                                                                histoRPtMC->GetXaxis()->FindBin(projRBins[i]+0.001),
                                                                histoRPtMC->GetXaxis()->FindBin(projRBins[i+1]+0.001),"e");
        cout << "Projecting Pt in R bin " << projRBins[i] << " - " << projRBins[i+1] << endl;
    }
    TH1F *histoRMCRebin = (TH1F*)histoRMC->Clone("histoRMCRebin");
    ConvGammaRebinWithBinCorrection(histoRMCRebin,rebinRPlots);

    TH1F *histoPtMC       = (TH1F*)histoRPtMC->ProjectionY("histoPtMC");
    Double_t *projPtBins[5]  = {0.4, 4., 8., 12., 20.};
    TH1F *histoRinPtBinMC[5];
    for(Int_t i=0; i<5; i++){
        histoRinPtBinMC[i] = (TH1F*)histoRPtMC->ProjectionY(Form("histoRinPtBinMC_%i",i),
                                                                histoRPtMC->GetYaxis()->FindBin(projPtBins[i]+0.001),
                                                                histoRPtMC->GetYaxis()->FindBin(projPtBins[i+1]+0.001),"e");
        cout << "Projecting R in Pt bin " << projPtBins[i] << " - " << projPtBins[i+1] << endl;
    }
    TH1F *histoPtMCRebin = (TH1F*)histoPtMC->Clone("histoPtMCRebin");
    ConvGammaRebinWithBinCorrection(histoPtMCRebin,rebinPtPlots);


    cout << "******************** projection true MC histograms... ******************** " << endl;
    TH1F *histoRTrueMC       = (TH1F*)histoRPtTrueMC->ProjectionX("histoRTrueMC");
    TH1F *histoPi0DalRTrueMC = (TH1F*)histoPi0DalRPtTrueMC->ProjectionX("histoPi0DalRTrueMC");
    TH1F *histoEtaDalRTrueMC = (TH1F*)histoEtaDalRPtTrueMC->ProjectionX("histoEtaDalRTrueMC");
    TH1F *histoCombRTrueMC   = (TH1F*)histoCombRPtTrueMC->ProjectionX("histoCombRTrueMC");
    Double_t *projRBins[4]  = {0., 5., 35., 180.};
    TH1F *histoPtinRBinTrueMC[4];
    for(Int_t i=0; i<4; i++){
        histoPtinRBinTrueMC[i] = (TH1F*)histoRPtTrueMC->ProjectionX(Form("histoPtinRBinTrueMC_%i",i),
                                                                histoRPtTrueMC->GetXaxis()->FindBin(projRBins[i]+0.001),
                                                                histoRPtTrueMC->GetXaxis()->FindBin(projRBins[i+1]+0.001),"e");
        cout << "Projecting Pt in R bin " << projRBins[i] << " - " << projRBins[i+1] << endl;
    }
    TH1F *histoRTrueMCRebin = (TH1F*)histoRTrueMC->Clone("histoRTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoRTrueMCRebin,rebinRPlots);
    TH1F *histoPi0DalRTrueMCRebin = (TH1F*)histoPi0DalRTrueMC->Clone("histoPi0DalRTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoPi0DalRTrueMCRebin,rebinRPlots);
    TH1F *histoEtaDalRTrueMCRebin = (TH1F*)histoEtaDalRTrueMC->Clone("histoEtaDalRTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoEtaDalRTrueMCRebin,rebinRPlots);
    TH1F *histoCombRTrueMCRebin = (TH1F*)histoCombRTrueMC->Clone("histoCombRTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoCombRTrueMCRebin,rebinRPlots);

    TH1F *histoPtTrueMC       = (TH1F*)histoRPtTrueMC->ProjectionY("histoPtTrueMC");
    TH1F *histoPi0DalPtTrueMC = (TH1F*)histoPi0DalRPtTrueMC->ProjectionY("histoPi0DalPtTrueMC");
    TH1F *histoEtaDalPtTrueMC = (TH1F*)histoEtaDalRPtTrueMC->ProjectionY("histoEtaDalPtTrueMC");
    TH1F *histoCombPtTrueMC   = (TH1F*)histoCombRPtTrueMC->ProjectionY("histoCombPtTrueMC");
    Double_t *projPtBins[5]  = {0.4, 4., 8., 12., 20.};
    TH1F *histoRinPtBinTrueMC[5];
    for(Int_t i=0; i<5; i++){
        histoRinPtBinTrueMC[i] = (TH1F*)histoRPtTrueMC->ProjectionY(Form("histoRinPtBinTrueMC_%i",i),
                                                                histoRPtTrueMC->GetYaxis()->FindBin(projPtBins[i]+0.001),
                                                                histoRPtTrueMC->GetYaxis()->FindBin(projPtBins[i+1]+0.001),"e");
        cout << "Projecting R in Pt bin " << projPtBins[i] << " - " << projPtBins[i+1] << endl;
    }
    TH1F *histoPtTrueMCRebin        = (TH1F*)histoPtTrueMC->Clone("histoPtTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoPtTrueMCRebin,rebinPtPlots);
    TH1F *histoPi0DalRTrueMCRebin   = (TH1F*)histoPi0DalPtTrueMC->Clone("histoPi0DalRTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoPi0DalRTrueMCRebin,rebinRPlots);
    TH1F *histoEtaDalRTrueMCRebin   = (TH1F*)histoEtaDalPtTrueMC->Clone("histoEtaDalRTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoEtaDalRTrueMCRebin,rebinRPlots);
    TH1F *histoCombRTrueMCRebin   = (TH1F*)histoCombPtTrueMC->Clone("histoCombRTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoCombRTrueMCRebin,rebinRPlots);


    //normalize with Nconversions in the gas
    cout << "scaling to Nconv in gas:" << endl;
    histoRDataScaledToGas        = ScaleByIntegralWithinLimits(histoRData, rMinGas, rMaxGas, 1. );
    histoRDataScaledToGasRebin   = (TH1F*)histoRDataScaledToGas->Rebin(nBinsR,"histoRDataScaledToGasRebin", arrayRBins);

    histoRMCScaledToGas          = ScaleByIntegralWithinLimits(histoRMC, rMinGas, rMaxGas, mcGasCorrectionFactor );
    histoRMCScaledToGasRebin     = (TH1F*)histoRMCScaledToGas->Rebin(nBinsR,"histoRMCScaledToGasRebin", arrayRBins);

    histoDataMCRatioRScaledToGas = (TH1F*)histoRDataScaledToGasRebin->Clone("histoDataMCRatioRScaledToGas");
    histoDataMCRatioRScaledToGas->Divide(histoRDataScaledToGasRebin,histoRMCScaledToGasRebin,1.,1.,"B");

    histoDataMCRatioR            = (TH1F*)histoRData->Clone("histoDataMCRatioR");
    histoDataMCRatioR->Divide(histoRData,histoRMC,1.,1.,"B");

    histoPurityR                 = (TH1F*)histoRTrueMC->Clone("histoPurityR");
    histoPurityR->Sumw2();
    histoPurityR->Divide(histoRTrueMC,histoRMC,1.,1.,"B");

    histoPurityPt                = (TH1F*)histoPtTrueMC->Clone("histoPurityPt");
    histoPurityPt->Sumw2();
    histoPurityPt->Divide(histoPtTrueMC,histoPtMC,1.,1.,"B");

    TH1F* histoPt5cmMC = (TH1F*)histoPtinRBinMC[1]->Clone("histoPt5cmMC");
      histoPt5cmMC->Sumw2();
      histoPt5cmMC->Add(histoPtinRBinMC[2]);
      histoPt5cmMC->Add(histoPtinRBinMC[3]);
    TH1F* histoPt5cmTrueMC = (TH1F*)histoPtinRBinTrueMC[1]->Clone("histoPt5cmTrueMC");
      histoPt5cmTrueMC->Sumw2();
      histoPt5cmTrueMC->Add(histoPtinRBinTrueMC[2]);
      histoPt5cmTrueMC->Add(histoPtinRBinTrueMC[3]);

    histoPurityPt5cm = (TH1F*)histoPt5cmTrueMC->Clone("histoPurityPt5cm");
    histoPurityPt5cm->Sumw2();
    histoPurityPt5cm->Divide(histoPt5cmTrueMC,histoPt5cmMC,1.,1.,"B");


    for(Int_t iR = 0; iR < nBinsR; iR++){

      histoPhiInRData[iR]        = (TH1D*)histoRPhiData->ProjectionX(Form("histoPhiInRData_%i",iR),
                                                                     histoRPhiData->GetYaxis()->FindBin(arrayRBins[iR]),
                                                                     histoRPhiData->GetYaxis()->FindBin(arrayRBins[iR+1]),"e");
      ConvGammaRebinWithBinCorrection(histoPhiInRData[iR],rebinPhiPlots);

      histoPhiInRMC[iR]          = (TH1D*)histoRPhiMC->ProjectionX(Form("histoPhiInRMC_%i",iR),
                                                                   histoRPhiMC->GetYaxis()->FindBin(arrayRBins[iR]),
                                                                   histoRPhiMC->GetYaxis()->FindBin(arrayRBins[iR+1]),"e");
      ConvGammaRebinWithBinCorrection(histoPhiInRMC[iR],rebinPhiPlots);

      histoPhiInRTrueMC[iR]      = (TH1D*)histoRPhiTrueMC->ProjectionX(Form("histoPhiInRTrueMC_%i",iR),
                                                                       histoRPhiTrueMC->GetYaxis()->FindBin(arrayRBins[iR]),
                                                                       histoRPhiTrueMC->GetYaxis()->FindBin(arrayRBins[iR+1]),"e");
      ConvGammaRebinWithBinCorrection(histoPhiInRTrueMC[iR],rebinPhiPlots);

      histoDataMCRatioPhiInR[iR] = (TH1D*)histoPhiInRData[iR]->Clone(Form("histoDataMCRatioPhiInR_%02d",iR));
      histoDataMCRatioPhiInR[iR]->Divide(histoPhiInRData[iR],histoPhiInRMC[iR]);


      if (iR > 9) rebinZPlots = 4;
      histoZInRData[iR]          =  (TH1D*)histoRZData->ProjectionX( Form("histoZInRData_%i",iR),
                                                                     histoRZData->GetYaxis()->FindBin(arrayRBins[iR]),
                                                                     histoRZData->GetYaxis()->FindBin(arrayRBins[iR+1]),"e");
      ConvGammaRebinWithBinCorrection(histoZInRData[iR],rebinZPlots);

      histoZInRMC[iR]            =  (TH1D*)histoRZMC->ProjectionX( Form("histoZInRMC_%i",iR),
                                                                   histoRZMC->GetYaxis()->FindBin(arrayRBins[iR]),
                                                                   histoRZMC->GetYaxis()->FindBin(arrayRBins[iR+1]),"e");
      ConvGammaRebinWithBinCorrection(histoZInRMC[iR],rebinZPlots);

      histoZInRTrueMC[iR]        =  (TH1D*)histoRZTrueMC->ProjectionX( Form("histoZInRTrueMC_%i",iR),
                                                                       histoRZTrueMC->GetYaxis()->FindBin(arrayRBins[iR]),
                                                                       histoRZTrueMC->GetYaxis()->FindBin(arrayRBins[iR+1]),"e");
      ConvGammaRebinWithBinCorrection(histoZInRTrueMC[iR],rebinZPlots);

      histoDataMCRatioZInR[iR]   = (TH1D*)histoZInRData[iR]->Clone(Form("histoDataMCRatioZInR_%02d",iR));
      histoDataMCRatioZInR[iR]->Divide(histoDataMCRatioZInR[iR],histoZInRMC[iR]);

    }

    //______________________________________ Multiplicity _________________________________________
    TCanvas * canvasNTracks = new TCanvas("canvasNTracks","",1200,1000);
    DrawGammaCanvasSettings( canvasNTracks, 0.1, 0.03, 0.05, 0.09);
    //canvasNTracks->Divide(2,2);
    canvasNTracks->SetLogy(1);
    TH2F * histoDummyNTracks = new TH2F("histoDummyNTracks","histoDummyNTracks",1000,0.,2000.,1000,1.e-10,10);
		SetStyleHistoTH2ForGraphs(histoDummyNTracks, "Good TPC tracks","Counts", 0.035,0.04,0.035,0.04,1.,1.);
        histoDummyNTracks->DrawCopy();

		DrawGammaSetMarker(histoGoodESDTracksData, 20, 1.5, colorData, colorData);
        histoGoodESDTracksData->Draw("same,hist");
		DrawGammaSetMarker(histoGoodESDTracksMC, 20, 1.5, colorMC, colorMC);
        histoGoodESDTracksMC->Draw("same,hist");

        TLegend* legend = GetAndSetLegend(0.75,0.75,2);
        legend->AddEntry(histoGoodESDTracksData,"Data","l");
        legend->AddEntry(histoGoodESDTracksMC,"MC","l");
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

        DrawGammaSetMarker(histoRDataScaledToGas, 20, 1.5, colorData, colorData);
        histoRDataScaledToGas->Draw("same,hist");
		DrawGammaSetMarker(histoRMCScaledToGas, 20, 1.5, colorMC, colorMC);
        histoRMCScaledToGas->Draw("same,hist");

        TLegend* legenUpPanel = GetAndSetLegend2(0.7, 0.8-(2*0.9*textsizeLabelsUp), 0.9, 0.8, textSizeLabels);
        legenUpPanel->AddEntry(histoRDataScaledToGas,"Data","l");
        legenUpPanel->AddEntry(histoRMCScaledToGas,"MC","l");
        legenUpPanel->Draw();

    padLower->cd();
    histoDummyTwoPanelsDown->DrawCopy();

        DrawGammaLines(0.,180,1., 1.,0.1,kGray);

        DrawGammaSetMarker(histoDataMCRatioRScaledToGas, 20, 1.5, colorData, colorData);
        histoDataMCRatioRScaledToGas->Draw("same");
		DrawGammaSetMarker(histoDataMCRatioRFullPtScaledToGas, 24, 1.5, colorData, colorData);
        histoDataMCRatioRFullPtScaledToGas->Draw("same");

        TLegend* legenLowPanel = GetAndSetLegend2(0.7, 0.85-(2*0.9*textsizeLabelsDown), 0.9, 0.85, textSizeLabels);
        legenLowPanel->AddEntry(histoDataMCRatioRScaledToGas,"From proj","lp");
        legenLowPanel->AddEntry(histoDataMCRatioRFullPtScaledToGas,"Full pT","lp");
        legenLowPanel->Draw();


    histoDummyTwoPanelsDown->Draw("same,axis");
    canvasTwoPanels->Print(Form("%s/PhotonConvRScaledToGas%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));



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

        DrawGammaSetMarker(histoRData, 20, 1.5, colorData, colorData);
        histoRData->Draw("same,hist");
		DrawGammaSetMarker(histoRMC, 20, 1.5, colorMC, colorMC);
        histoRMC->Draw("same,hist");
		DrawGammaSetMarker(histoRTrueMC, 20, 1.5, colorTrueMC, colorTrueMC);
        histoRTrueMC->Draw("same,hist");
		DrawGammaSetMarker(histoRTrueCombMC, 20, 1.5, colorTrueCombMC, colorTrueCombMC);
        histoRTrueCombMC->DrawCopy("same,hist");
		DrawGammaSetMarker(histoRTruePi0DalMC, 20, 1.5, kBlue-9, kBlue-9);
        histoRTruePi0DalMC->SetFillColor(kBlue-9);
        histoRTruePi0DalMC->SetFillStyle(3244);
        histoRTruePi0DalMC->DrawCopy("same,hist");
		DrawGammaSetMarker(histoRTrueEtaDalMC, 20, 1.5, kBlue-3, kBlue-3);
        histoRTrueEtaDalMC->SetFillColor(kBlue-3);
        histoRTrueEtaDalMC->SetFillStyle(3002);
        histoRTrueEtaDalMC->DrawCopy("same,hist");

        TLegend* legenRdistrib = GetAndSetLegend2(0.63, 0.9-(6*0.9*textsizeLabelsUp), 0.9, 0.9, textSizeLabels);
        legenRdistrib->AddEntry(histoRData,"Data","l");
        legenRdistrib->AddEntry(histoRMC,"MC","l");
        legenRdistrib->AddEntry(histoRTrueMC,"True MC","l");
        legenRdistrib->AddEntry(histoRTrueCombMC,"True comb MC","l");
        legenRdistrib->AddEntry(histoRTruePi0DalMC,"True MC #pi^{0} Dal.","lf");
        legenRdistrib->AddEntry(histoRTrueEtaDalMC,"True MC #eta Dal.","fl");
        legenRdistrib->Draw();

    histoDummy4PanelsUp->Draw("axis,same");
    padRDistribLowerLeft->cd();
    histoDummy4PanelsDown->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++)
            DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,0.1,kGray);
        DrawGammaLines(0.,180,1., 1.,0.1,kGray);

        DrawGammaSetMarker(histoDataMCRatioR, 20, 1.5, colorData, colorData);
        histoDataMCRatioR->Draw("same,histo");

    histoDummy4PanelsDown->Draw("axis,same");
    padRDistribUpperRight->cd();
    padRDistribUpperRight->SetLogy();
    histoDummy4PanelsUp->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++)
            DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,0.1,kGray);

        histoRData->Draw("same,histo");
        histoRMC->Draw("same,histo");
        histoRTrueMC->Draw("same,histo");

        for(Int_t i=0; i<5; i++){
            DrawGammaSetMarker(histoRinPtBinData[i], 20, 1.5, colorData, colorData);
            histoRinPtBinData[i]->SetLineStyle(2);
            histoRinPtBinData[i]->DrawCopy("same,histo");

            DrawGammaSetMarker(histoRinPtBinMC[i], 20, 1.5, colorMC, colorMC);
            histoRinPtBinMC[i]->SetLineStyle(2);
            histoRinPtBinMC[i]->DrawCopy("same,histo");

            DrawGammaSetMarker(histoRinPtBinTrueMC[i], 20, 1.5, colorTrueMC, colorTrueMC);
            histoRinPtBinTrueMC[i]->SetLineStyle(2);
            histoRinPtBinTrueMC[i]->DrawCopy("same,histo");
        }

        TLegend* legenRdistribPtcuts = GetAndSetLegend2(0.05, 0.35-(6.2*0.9*textsizeLabelsUp), 0.35, 0.35, textSizeLabels);
        legenRdistribPtcuts->AddEntry(histoRData,"Data","l");
        legenRdistribPtcuts->AddEntry(histoRMC,"MC","l");
        legenRdistribPtcuts->AddEntry(histoRTrueMC,"True MC","l");
        legenRdistribPtcuts->AddEntry(histoRinPtBinData[1],"Data, 0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","l");
        legenRdistribPtcuts->AddEntry(histoRinPtBinMC[1],"MC, 0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","lf");
        legenRdistribPtcuts->AddEntry(histoRinPtBinTrueMC[1],"True MC, 0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","fl");
        legenRdistribPtcuts->Draw();

    histoDummy4PanelsUp->Draw("axis,same");
    padRDistribLowerRight->cd();
    histoDummy4PanelsDown->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++)
            DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,0.1,kGray);
        DrawGammaLines(0.,180,1., 1.,0.1,kGray);

        histoDataMCRatioR->Draw("same,histo");

        TLegend* legenRdistribPtcutsRatio = GetAndSetLegend2(0.05, 0.27-(2*0.9*textsizeLabelsDown), 0.35, 0.27, textSizeLabels);
        legenRdistribPtcutsRatio->AddEntry(histoDataMCRatioR,"No #it{p}_{T} cut","l");
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
        DrawGammaSetMarker(histoPurityR, 20, 1.5, colorData, colorData);

        histoPurityR->Draw("same,hist");

        TLegend* legenPurity = GetAndSetLegend2(0.15, 0.9-(1*0.9*textsizeLabelsUppurity), 0.7, 0.9, textSizeLabels);
        legenPurity->AddEntry((TObject*)0,Form("Cut: %s",fCutSelectionRead.Data()),"");
        legenPurity->Draw();

    histoDummyPurityR->Draw("same,axis");

    padLowerPurity->cd();
    histoDummyPurityPt->DrawCopy();

        DrawGammaLines(0.,8.,1., 1.,0.1,kGray+1);
        DrawGammaSetMarker(histoPurityPt, 20, 1.5, colorData, colorData);
        histoPurityPt->Draw("same,hist");
        DrawGammaSetMarker(histoPurityPt5cm, 20, 1.5, kGreen+2, kGreen+2);
        histoPurityPt5cm->Draw("same,hist");

        TLegend* legenPurityPt = GetAndSetLegend2(0.75, 0.55-(2*0.9*textsizeLabelsUppurity), 0.9, 0.55, textSizeLabels);
        legenPurityPt->AddEntry(histoPurityPt,"no R cut","l");
        legenPurityPt->AddEntry(histoPurityPt5cm,"R > 5 cm","l");
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

        DrawGammaSetMarker(histoDCAData, 20, 1.5, colorData, colorData);
        histoDCAData->Draw("same,hist");
		DrawGammaSetMarker(histoDCAMC, 20, 1.5, colorMC, colorMC);
        histoDCAMC->Draw("same,hist");

    histoDummyDCA->Draw("axis,same");
    padPhotonCharLowerLeft->cd();
    padPhotonCharLowerLeft->SetLogy();
    TH2F * histoDummyMass = new TH2F("histoDummyMass","histoDummyMass",1000,0.,0.2,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyMass, "Mass (GeV)","Counts ", 0.85*textsizeLabelsDown, textsizeLabelsDown,0.85*textsizeLabelsDown,textsizeLabelsDown, 0.9,0.92);
    histoDummyMass->GetYaxis()->SetRangeUser(1.e-6,2e-3);
    histoDummyMass->DrawCopy();

        DrawGammaSetMarker(histoMassData, 20, 1.5, colorData, colorData);
        histoMassData->Draw("same,hist");
		DrawGammaSetMarker(histoMassMC, 20, 1.5, colorMC, colorMC);
        histoMassMC->Draw("same,hist");

    histoDummyMass->Draw("axis,same");
    padPhotonCharUpperRight->cd();
    padPhotonCharUpperRight->SetLogy();
    TH2F *histoDummyPsiPair = new TH2F("histoDummyPsiPair","histoDummyPsiPair",1000,0.,0.3,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyPsiPair, "#Psi_{pair}","Counts", 0.85*textsizeLabelsUp, textsizeLabelsUp,0.85*textsizeLabelsUp,textsizeLabelsUp, 0.9,0.91);
    histoDummyPsiPair->GetYaxis()->SetRangeUser(1.e-6,2.e-2);
    histoDummyPsiPair->DrawCopy();

        DrawGammaSetMarker(histoPsiPairData, 20, 1.5, colorData, colorData);
        histoPsiPairData->Draw("same,hist");
		DrawGammaSetMarker(histoPsiPairMC, 20, 1.5, colorMC, colorMC);
        histoPsiPairMC->Draw("same,hist");

    histoDummyPsiPair->Draw("axis,same");
    padPhotonCharLowerRight->cd();
    padPhotonCharLowerRight->SetLogy();
    TH2F * histoDummyChi2 = new TH2F("histoDummyChi2","histoDummyChi2",1000,0.,40.,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyChi2, "#chi^{2}","Counts ", 0.85*textsizeLabelsDown, textsizeLabelsDown,0.85*textsizeLabelsDown,textsizeLabelsDown, 0.9,0.92);
    histoDummyChi2->GetYaxis()->SetRangeUser(1.e-6,2.e-2);
    histoDummyChi2->DrawCopy();

        DrawGammaSetMarker(histoChi2Data, 20, 1.5, colorData, colorData);
        histoChi2Data->Draw("same,hist");
		DrawGammaSetMarker(histoChi2MC, 20, 1.5, colorMC, colorMC);
        histoChi2MC->Draw("same,hist");

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
        histoDummyPhiInRbinsUp->GetYaxis()->SetRangeUser(histoPhiInRData[iR]->GetMinimum(),histoPhiInRData[iR]->GetMaximum()*10);
        histoDummyPhiInRbinsUp->DrawCopy();

            DrawGammaSetMarker(histoPhiInRData[iR], 20, 1.5, colorData, colorData);
            histoPhiInRData[iR]->Draw("same,hist");
            DrawGammaSetMarker(histoPhiInRMC[iR], 20, 1.5, colorMC, colorMC);
            histoPhiInRMC[iR]->Draw("same,hist");

            TLatex *latexBinning = new TLatex(0.15,0.85,Form("#varphi in range %s (%s)",arrayRangesRBins[iR].Data(), arrayNamesRBins[iR].Data()));
            SetStyleTLatex( latexBinning, sizeTextNameBins,2);
            latexBinning->Draw();

            TLegend* legenPhiInRbins = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
            legenPhiInRbins->AddEntry(histoPhiInRData[iR],"Data","l");
            legenPhiInRbins->AddEntry(histoPhiInRMC[iR],"MC","l");
            legenPhiInRbins->Draw();

        histoDummyPhiInRbinsUp->Draw("same,axis");
        padPhiInRbinsLower->cd();
        histoDummyPhiInRbinsDown->DrawCopy();

            DrawGammaLines(0.,2*TMath::Pi(),1., 1.,0.1,kGray+1);

            DrawGammaSetMarker(histoDataMCRatioPhiInR[iR], 20, 1.5, colorData, colorData);
            histoDataMCRatioPhiInR[iR]->Draw("same,hist");

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
//         histoDummyZInRbinsUp->GetYaxis()->SetRangeUser(histoZInRData[iR]->GetMinimum()*0.1,histoZInRData[iR]->GetMaximum()*10);
        histoDummyZInRbinsUp->GetYaxis()->SetRangeUser(0.,histoZInRData[iR]->GetMaximum()*1.5);
        histoDummyZInRbinsUp->DrawCopy();

            DrawGammaLines(-10,-10, histoZInRData[iR]->GetMinimum()*0.1,histoZInRData[iR]->GetMaximum()*1.5,0.1,kGray+1,2);
            DrawGammaLines(10,10, histoZInRData[iR]->GetMinimum()*0.1,histoZInRData[iR]->GetMaximum()*1.5,0.1,kGray+1,2);

            DrawGammaSetMarker(histoZInRData[iR], 20, 1.5, colorData, colorData);
            histoZInRData[iR]->Draw("same,hist");
            DrawGammaSetMarker(histoZInRMC[iR], 20, 1.5, colorMC, colorMC);
            histoZInRMC[iR]->Draw("same,hist");

            TLatex *latexBinning = new TLatex(0.15,0.85,Form("Z in range %s (%s)",arrayRangesRBins[iR].Data(), arrayNamesRBins[iR].Data()));
            SetStyleTLatex( latexBinning, sizeTextNameBins,2);
            latexBinning->Draw();

            TLegend* legenZInRbins = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
            legenZInRbins->AddEntry(histoZInRData[iR],"Data","l");
            legenZInRbins->AddEntry(histoZInRMC[iR],"MC","l");
            legenZInRbins->Draw();

        histoDummyZInRbinsUp->Draw("same,axis");
        padZInRbinsLower->cd();
        histoDummyZInRbinsDown->DrawCopy();

            DrawGammaLines(-10,-10, 0.3,1.55,0.1,kGray+1,2);
            DrawGammaLines(10,10,0.3,1.55,0.1,kGray+1,2);
            DrawGammaLines(-20,20,1., 1.,0.1,kGray+1);

            DrawGammaSetMarker(histoDataMCRatioZInR[iR], 20, 1.5, colorData, colorData);
            histoDataMCRatioZInR[iR]->Draw("same,hist");

        histoDummyZInRbinsDown->Draw("same,axis");
        canvasZInRbins->Update();
        canvasZInRbins->SaveAs(Form("%s/ZInRbins%s_%i_%s.%s",outputDirectory.Data(),optionPeriod.Data(),iR,fCutSelectionRead.Data(),suffix.Data()));

    }
    delete canvasZInRbins;


    TProfile* fProfileContainingMaterialBudgetWeights = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins","profileContainingMaterialBudgetWeights_manyRadialBins",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeights->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeights->GetYaxis()->SetTitle("weight_per_gamma");

    for(Int_t i=0; i<nBinsR; i++){
        cout<< arrayRBins[i] << " " << histoDataMCRatioR->GetBinContent(i+1) <<endl;
        fProfileContainingMaterialBudgetWeights->Fill(arrayRBins[i],histoDataMCRatioR->GetBinContent(i+1));
    }

    cout<< " pT 300"<< endl;
    TProfile* fProfileContainingMaterialBudgetWeightsFullPt = new TProfile("profileContainingMaterialBudgetWeightsFullPt_manyRadialBins","profileContainingMaterialBudgetWeightsFullPt_manyRadialBins",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsFullPt->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsFullPt->GetYaxis()->SetTitle("weight_per_gamma");

    for(Int_t i=0; i<nBinsR; i++){
        cout << "R bin: " << arrayRBins[i] << " ratio Data/MC scaled to gas " << histoDataMCRatioRFullPtScaledToGas->GetBinContent(i+1) <<endl;
        fProfileContainingMaterialBudgetWeightsFullPt->Fill(arrayRBins[i],histoDataMCRatioRFullPtScaledToGas->GetBinContent(i+1));
    }
    //Form("%s/PhotonChar_%s.%s",optionPeriod.Data(),fCutSelectionRead.Data() ))
    cout<< Form("MCInputFileMaterialBudgetWeights%s_%s.root",optionPeriod.Data(),fCutSelectionRead.Data())<< endl;
    TFile outFile(Form("MCInputFileMaterialBudgetWeights%s_%s.root",optionPeriod.Data(),fCutSelectionRead.Data()) ,"RECREATE");

        fProfileContainingMaterialBudgetWeights->Write();
        fProfileContainingMaterialBudgetWeightsFullPt->Write();

        histoRData->Write("Data");
        histoRMC->Write("MC");
        histoRDataScaledToGas->Write("DataScaledToGas");
        histoRMCScaledToGas->Write("MCScaledToGas");
    outFile.Close();


}
