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

Double_t CalculateIntegralWithinLimits(const TH1F* hist, Double_t rmin, Double_t rmax, Double_t mcGasCorrectionFactor )
{
    TH1F * histToBeIntegrated = (TH1F*)hist->Clone();
    Double_t nconvInRange =  histToBeIntegrated->Integral(hist->GetXaxis()->FindBin(rmin),hist->GetXaxis()->FindBin(rmax),"width");
    nconvInRange*=mcGasCorrectionFactor;
    cout << "nconvInRange = " << nconvInRange << "  "  << hist->GetXaxis()->FindBin(rmin) << " " << hist->GetXaxis()->FindBin(rmax)<< endl;
    return nconvInRange;
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



    TString outputDirectory= Form("%s/%s/%s/%s",outputFolderName.Data(),fEnergyFlag.Data(),optionPeriod.Data(),fCutSelectionRead.Data());
    gSystem->Exec("mkdir -p "+outputDirectory);
    if(optionPeriod.CompareTo(""))
        outputDirectory= Form("%s/%s/%s/%s",outputFolderName.Data(),fEnergyFlag.Data(),optionPeriod.Data(),fCutSelectionRead.Data());

    cout<< "Output directory is::"<<  outputDirectory.Data()<<endl;

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
    TH1F *histoGoodESDTracksData  = (TH1F*)ESDContainer->FindObject("GoodESDTracksEta09");
    TH2F *histoRPhiData           = (TH2F*)ESDContainer->FindObject("ESD_Conversion_RPhi");
    TH2F *histoREtaData           = (TH2F*)ESDContainer->FindObject("ESD_Conversion_REta");
    TH2F *histoRZData             = (TH2F*)ESDContainer->FindObject("ESD_Conversion_RZ");
    TH2F *histoRPtData            = (TH2F*)ESDContainer->FindObject("ESD_Conversion_RPt");
    TH2F *histoRElecdEdxData      = (TH2F*)ESDContainer->FindObject("Electron_RdEdx");
    TH2F *histoRElecNSdEdxData    = (TH2F*)ESDContainer->FindObject("Electron_RNSigmadEdx");
    TH1F *histoDCAData            = (TH1F*)ESDContainer->FindObject("ESD_Conversion_DCA");
    TH1F *histoPsiPairData        = (TH1F*)ESDContainer->FindObject("ESD_Conversion_PsiPair");
    TH1F *histoChi2Data           = (TH1F*)ESDContainer->FindObject("ESD_Conversion_Chi2");
    TH1F *histoMassData           = (TH1F*)ESDContainer->FindObject("ESD_Conversion_Mass");
    TH1F *histoRRejLargeData      = (TH1F*)ESDContainer->FindObject("ESD_Conversion_RLarge");
    TH1F *histoRRejSmallData      = (TH1F*)ESDContainer->FindObject("ESD_Conversion_RSmall");
    TH2F *histoAsymPData          = (TH2F*)ESDContainer->FindObject("ESD_ConversionMapping_AsymP");

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
    TH1F *histoGoodESDTracksMC    = (TH1F*)ESDContainerMC->FindObject("GoodESDTracksEta09");
    TH2F *histoRPhiMC             = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_RPhi");
    TH2F *histoREtaMC             = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_REta");
    TH2F *histoRZMC               = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_RZ");
    TH2F *histoRPtMC              = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_RPt");
    TH2F *histoRElecdEdxMC        = (TH2F*)ESDContainerMC->FindObject("Electron_RdEdx");
    TH2F *histoRElecNSdEdxMC      = (TH2F*)ESDContainerMC->FindObject("Electron_RNSigmadEdx");
    TH1F *histoDCAMC              = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_DCA");
    TH1F *histoPsiPairMC          = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_PsiPair");
    TH1F *histoChi2MC             = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_Chi2");
    TH1F *histoMassMC             = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_Mass");
    TH1F *histoRRejLargeMC        = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_RLarge");
    TH1F *histoRRejSmallMC        = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_RSmall");


    TList *MCContainer            = (TList*)HistosGammaConversionMC->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
    TH2F * histoConvRPtMC          = (TH2F*)MCContainer->FindObject("MC_Conversion_RPt");      // All Converted
    TH2F * histoConvRPhiMC          = (TH2F*)MCContainer->FindObject("MC_Conversion_RPhi");      // All Converted
    TH2F * histoConvREtaMC          = (TH2F*)MCContainer->FindObject("MC_Conversion_REta");      // All Converted


    TList *TrueMCContainer        = (TList*)HistosGammaConversionMC->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
    TH2F *histoRPhiTrueMC         = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RPhi");
    TH2F *histoREtaTrueMC         = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_REta");
    TH2F *histoRZTrueMC           = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RZ");
    TH2F *histoRPtTrueMC          = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RPt");
    TH2F *histoRPtMCRPtTrueMC     = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RPtMCRPt");
    TH1F *histoDCATrueMC          = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_DCA");
    TH1F *histoPsiPairTrueMC      = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_PsiPair");
    TH1F *histoChi2TrueMC         = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_Chi2");
    TH1F *histoMassTrueMC         = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_Mass");
    TH1F *histoRRejLargeTrueMC    = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_RLarge");
    TH1F *histoRRejSmallTrueMC    = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_RSmall");
    TH2F *histoAsymPTrueMC        = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversionMapping_AsymP");

    TH2F *histoPi0DalRPtTrueMC    = (TH2F*)TrueMCContainer->FindObject("ESD_TruePi0DalConversion_RPt");
    TH2F *histoPi0DalREtaTrueMC   = (TH2F*)TrueMCContainer->FindObject("ESD_TruePi0DalConversion_REta");
    TH2F *histoEtaDalRPtTrueMC    = (TH2F*)TrueMCContainer->FindObject("ESD_TrueEtaDalConversion_RPt");
    TH2F *histoEtaDalREtaTrueMC   = (TH2F*)TrueMCContainer->FindObject("ESD_TrueEtaDalConversion_REta");
    TH2F *histoCombRPtTrueMC      = (TH2F*)TrueMCContainer->FindObject("ESD_TrueCombinatorialConversion_RPt");
    TH2F *histoCombREtaTrueMC     = (TH2F*)TrueMCContainer->FindObject("ESD_TrueCombinatorialConversion_REta");


    //________________________ normalisation _______________________
    //using the charged particles (~pions) to normalise instead of the neutral pions (~tot gammas)

    Float_t numberGoodEventsData  = histoEventQualityData->GetBinContent(1);

    Double_t meanMultData         = histoGoodESDTracksData->GetMean();


    Float_t normFactorReconstData = 1./(numberGoodEventsData*meanMultData);
    //    normFactorReconstData = 1;
    Float_t numberGoodEventsMC    = histoEventQualityMC->GetBinContent(1);
    Double_t meanMultMC           = histoGoodESDTracksMC->GetMean();
    Float_t normFactorReconstMC   = 1./numberGoodEventsMC*1./meanMultMC;
    //    normFactorReconstMC = 1;
    cout<< "Number of Good Events Data, histoEventQualityData->GetBinContent(1):: "<< numberGoodEventsData<< endl;
    cout<< "Number of Good Events MC, histoEventQualityMC->GetBinContent(1):: "<< numberGoodEventsMC<< endl;

    cout << "Normalization factor data " << meanMultData << " and MC " << meanMultMC << " -> Data/MC: " << meanMultData/meanMultMC << endl;

    histoGoodESDTracksData->Scale(1./numberGoodEventsData);
    histoGoodESDTracksMC->Scale(1./numberGoodEventsMC);


    //-AM   Don't do scaling prior to calculation of errors
    // GammaScalingHistogramm(histoRPtData,normFactorReconstData);
    // GammaScalingHistogramm(histoRPtMC,normFactorReconstMC);
    // GammaScalingHistogramm(histoRPtTrueMC,normFactorReconstMC);
    // GammaScalingHistogramm(histoPi0DalRPtTrueMC,normFactorReconstMC);
    // GammaScalingHistogramm(histoEtaDalRPtTrueMC,normFactorReconstMC);
    // GammaScalingHistogramm(histoCombRPtTrueMC,normFactorReconstMC);

    GammaScalingHistogramm(histoRElecdEdxData,normFactorReconstData);
    GammaScalingHistogramm(histoRElecdEdxMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRElecNSdEdxData,normFactorReconstData);
    GammaScalingHistogramm(histoRElecNSdEdxMC,normFactorReconstMC);

    GammaScalingHistogramm(histoDCAData,normFactorReconstData);
    GammaScalingHistogramm(histoDCAMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPsiPairData,normFactorReconstData);
    GammaScalingHistogramm(histoPsiPairMC,normFactorReconstMC);
    GammaScalingHistogramm(histoChi2Data,normFactorReconstData);
    GammaScalingHistogramm(histoChi2MC,normFactorReconstMC);
    GammaScalingHistogramm(histoMassData,normFactorReconstData);
    GammaScalingHistogramm(histoMassMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRRejLargeData,normFactorReconstData);
    GammaScalingHistogramm(histoRRejLargeMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRRejSmallData,normFactorReconstData);
    GammaScalingHistogramm(histoRRejSmallMC,normFactorReconstMC);


    //___________________ Projecting histos ___________________

    cout << "******************** projection data histograms... ******************** " << endl;
    TH1F *histoRData        = (TH1F*)histoRPtData->ProjectionY("histoRData");
    Double_t projRBins[4]  = {0., 5., 35., 180.};
    TH1F *histoPtinRBinData[3];
    Double_t epsilon=0.001;
    Double_t eps=0.001;
    for(Int_t i=0; i<3; i++){
        histoPtinRBinData[i] = (TH1F*)histoRPtData->ProjectionX(Form("histoPtinRBinData_%i",i),
                                                                histoRPtData->GetYaxis()->FindBin(projRBins[i]+0.001),
                                                                histoRPtData->GetYaxis()->FindBin(projRBins[i+1]-0.001),"e");
        cout << "Projecting Pt in R bin " << projRBins[i] << " - " << projRBins[i+1] << endl;
	cout << "Bins::"<<  histoRPtData->GetYaxis()->FindBin(projRBins[i]+0.001) << "  "<<  histoRPtData->GetYaxis()->FindBin(projRBins[i+1]-0.001)<< endl;

    }
    // TH1F *histoRDataRebin = (TH1F*)histoRData->Clone("histoRDataRebin");
    // ConvGammaRebinWithBinCorrection(histoRDataRebin,rebinRPlots);

    TH1F *histoPtData       = (TH1F*)histoRPtData->ProjectionX("histoPtData");
    Double_t projPtBins[6]  = {0.15,0.3, 0.4, 0.5, 0.6, 0.7};
    TH1F *histoRinPtBinData[6];
    for(Int_t i=0; i<6; i++){
        histoRinPtBinData[i] = (TH1F*)histoRPtData->ProjectionY(Form("histoRinPtBinData_%i",i),
                                                                histoRPtData->GetXaxis()->FindBin(projPtBins[i]+0.001),
                                                                histoRPtData->GetNbinsX(),"e");
	//                                                                histoRPtData->GetXaxis()->FindBin(projPtBins[i+1]-0.001),"e");

        cout << "Projecting R in Pt bin, lower limit to maximum " << projPtBins[i] << " - " << histoRPtData->GetXaxis()->GetBinCenter( histoRPtData->GetNbinsX()) << " "<< histoRinPtBinData[i]->GetMean()<< endl;
    }
    TH1F *histoPtDataRebin = (TH1F*)histoPtData->Clone("histoPtDataRebin");
    ConvGammaRebinWithBinCorrection(histoPtDataRebin,rebinPtPlots);

    TH1F *histoEtaData        = (TH1F*)histoREtaData->ProjectionX("histoEtaData");

    TH1F *histoAsymPDataLowP;
    TH1F *histoAsymPDataHighP;
    histoAsymPDataLowP = (TH1F*)histoAsymPData->ProjectionY("histoAsymPDataLowP",
							    1,
							    histoAsymPData->GetXaxis()->FindBin(0.4-0.001),"e");

    histoAsymPDataHighP = (TH1F*)histoAsymPData->ProjectionY("histoAsymPDataHighP",
							     histoAsymPData->GetXaxis()->FindBin(0.4+0.001),
							     histoAsymPData->GetNbinsX(),"e");
    ConvGammaRebinWithBinCorrection(histoAsymPDataLowP,4);
    ConvGammaRebinWithBinCorrection(histoAsymPDataHighP,4);


    cout<< " Projection of MC Histograms , all converted"<< endl;

    TH1F *histoConvRMC       = (TH1F*)histoConvRPtMC->ProjectionY("histoConvRMC");
    TH1F *histoConvPtMC      = (TH1F*)histoConvRPtMC->ProjectionX("histoConvPtMC");
    TH1F *histoConvPhiMC      = (TH1F*)histoConvRPhiMC->ProjectionX("histoConvPhiMC");
    TH1F *histoConvEtaMC      = (TH1F*)histoConvREtaMC->ProjectionX("histoConvEtaMC");

    cout << "******************** projection MC histograms... ******************** " << endl;
    TH1F *histoRMC        = (TH1F*)histoRPtMC->ProjectionY("histoRMC");
    //Double_t *   projRBins  = {0., 5., 35., 180.};
    TH1F *histoPtinRBinMC[3];
    for(Int_t i=0; i<3; i++){
        histoPtinRBinMC[i] = (TH1F*)histoRPtMC->ProjectionX(Form("histoPtinRBinMC_%i",i),
                                                                histoRPtMC->GetYaxis()->FindBin(projRBins[i]+0.001),
                                                                histoRPtMC->GetYaxis()->FindBin(projRBins[i+1]-0.001),"e");
        cout << "Projecting Pt in R bin " << projRBins[i] << " - " << projRBins[i+1] << endl;
    }
    TH1F *histoRMCRebin = (TH1F*)histoRMC->Clone("histoRMCRebin");
    ConvGammaRebinWithBinCorrection(histoRMCRebin,rebinRPlots);

    TH1F *histoPtMC       = (TH1F*)histoRPtMC->ProjectionX("histoPtMC");
    //    Double_t * projPtBins[5]  = {0.4, 4., 8., 12., 20.};
    TH1F *histoRinPtBinMC[6];
    for(Int_t i=0; i<6; i++){
        histoRinPtBinMC[i] = (TH1F*)histoRPtMC->ProjectionY(Form("histoRinPtBinMC_%i",i),
							    histoRPtMC->GetXaxis()->FindBin(projPtBins[i]+0.001),
							    histoRPtMC->GetNbinsX(),"e");
	//                                                                histoRPtMC->GetXaxis()->FindBin(projPtBins[i+1]-0.001),"e");
        cout << "Projecting R in Pt bin , lower limit to maximum" << projPtBins[i] << " - " << histoRPtMC->GetXaxis()->GetBinCenter( histoRPtMC->GetNbinsX()) << endl;
    }
    TH1F *histoPtMCRebin = (TH1F*)histoPtMC->Clone("histoPtMCRebin");
    ConvGammaRebinWithBinCorrection(histoPtMCRebin,rebinPtPlots);

    TH1F *histoEtaMC       = (TH1F*)histoREtaMC->ProjectionX("histoEtaMC");

    cout << "******************** projection true MC histograms... ******************** " << endl;
    TH1F *histoRTrueMC       = (TH1F*)histoRPtTrueMC->ProjectionY("histoRTrueMC");
    TH1F *histoRMCRTrueMC       = (TH1F*)histoRPtMCRPtTrueMC->ProjectionY("histoRMCRTrueMC");
    TH1F *histoPi0DalRTrueMC = (TH1F*)histoPi0DalRPtTrueMC->ProjectionY("histoPi0DalRTrueMC");
    TH1F *histoEtaDalRTrueMC = (TH1F*)histoEtaDalRPtTrueMC->ProjectionY("histoEtaDalRTrueMC");
    TH1F *histoCombRTrueMC   = (TH1F*)histoCombRPtTrueMC->ProjectionY("histoCombRTrueMC");
    //    Double_t *projRBins[4]  = {0., 5., 35., 180.};
    TH1F *histoPtinRBinTrueMC[3];
    for(Int_t i=0; i<3; i++){
        histoPtinRBinTrueMC[i] = (TH1F*)histoRPtTrueMC->ProjectionX(Form("histoPtinRBinTrueMC_%i",i),
                                                                histoRPtTrueMC->GetYaxis()->FindBin(projRBins[i]+0.001),
                                                                histoRPtTrueMC->GetYaxis()->FindBin(projRBins[i+1]-0.001),"e");
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

    TH1F *histoPtTrueMC       = (TH1F*)histoRPtTrueMC->ProjectionX("histoPtTrueMC");
    TH1F *histoPtMCPtTrueMC       = (TH1F*)histoRPtMCRPtTrueMC->ProjectionX("histoPtMCPtTrueMC");
    TH1F *histoPi0DalPtTrueMC = (TH1F*)histoPi0DalRPtTrueMC->ProjectionX("histoPi0DalPtTrueMC");
    TH1F *histoEtaDalPtTrueMC = (TH1F*)histoEtaDalRPtTrueMC->ProjectionX("histoEtaDalPtTrueMC");
    TH1F *histoCombPtTrueMC   = (TH1F*)histoCombRPtTrueMC->ProjectionX("histoCombPtTrueMC");
    // Double_t *projPtBins[5]  = {0.4, 4., 8., 12., 20.};
    TH1F *histoRinPtBinTrueMC[6];
    for(Int_t i=0; i<6; i++){
        histoRinPtBinTrueMC[i] = (TH1F*)histoRPtTrueMC->ProjectionY(Form("histoRinPtBinTrueMC_%i",i),
                                                                histoRPtTrueMC->GetXaxis()->FindBin(projPtBins[i]+0.001),
                                                                histoRPtTrueMC->GetNbinsX(),"e");
	//                                                                histoRPtTrueMC->GetXaxis()->FindBin(projPtBins[i+1]-0.001),"e");
        cout << "Projecting R in Pt bin " << projPtBins[i] << " - " << projPtBins[i+1] << endl;
    }
    TH1F *histoPtTrueMCRebin        = (TH1F*)histoPtTrueMC->Clone("histoPtTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoPtTrueMCRebin,rebinPtPlots);
    // histograms are re-defined. To be cheked/understood! Should it be Pt AM. 21.03.18 ???
    TH1F *histoPi0DalPtTrueMCRebin   = (TH1F*)histoPi0DalPtTrueMC->Clone("histoPi0DalPtTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoPi0DalPtTrueMCRebin,rebinPtPlots);
    TH1F *histoEtaDalPtTrueMCRebin   = (TH1F*)histoEtaDalPtTrueMC->Clone("histoEtaDalPtTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoEtaDalPtTrueMCRebin,rebinPtPlots);
    TH1F *histoCombPtTrueMCRebin   = (TH1F*)histoCombPtTrueMC->Clone("histoCombPtTrueMCRebin");
    ConvGammaRebinWithBinCorrection(histoCombPtTrueMCRebin,rebinPtPlots);

    TH1F *histoAsymPTrueMCLowP;
    TH1F *histoAsymPTrueMCHighP;
    histoAsymPTrueMCLowP = (TH1F*)histoAsymPTrueMC->ProjectionY("histoAsymPTrueMCLowP",
							    1,
							    histoAsymPTrueMC->GetXaxis()->FindBin(0.4-0.001),"e");

    histoAsymPTrueMCHighP = (TH1F*)histoAsymPTrueMC->ProjectionY("histoAsymPTrueMCHighP",
							     histoAsymPTrueMC->GetXaxis()->FindBin(0.4+0.001),
							     histoAsymPTrueMC->GetNbinsX(),"e");


    ConvGammaRebinWithBinCorrection(histoAsymPTrueMCLowP,4);
    ConvGammaRebinWithBinCorrection(histoAsymPTrueMCHighP,4);
    //normalize with Nconversions in the gas
    cout << "Data scaling to Nconv in gas:" << endl;
    Double_t nconvInRangeData = CalculateIntegralWithinLimits(histoRData, rMinGas+eps, rMaxGas-eps,1.);
    Double_t dataStatErrorGas = TMath::Sqrt(nconvInRangeData);
    Double_t dataStatRelErrorGas = TMath::Sqrt(nconvInRangeData)/nconvInRangeData;

    Double_t nconvInRangeData01 = CalculateIntegralWithinLimits(histoRinPtBinData[0], rMinGas+eps, rMaxGas-eps,1.);
    Double_t dataStatErrorGas01 = TMath::Sqrt(nconvInRangeData01);
    Double_t dataStatRelErrorGas01 = TMath::Sqrt(nconvInRangeData01)/nconvInRangeData01;

    Double_t nconvInRangeData02 = CalculateIntegralWithinLimits(histoRinPtBinData[1], rMinGas+eps, rMaxGas-eps,1.);
    Double_t dataStatErrorGas02 = TMath::Sqrt(nconvInRangeData02);
    Double_t dataStatRelErrorGas02 = TMath::Sqrt(nconvInRangeData02)/nconvInRangeData02;

    Double_t nconvInRangeData03 = CalculateIntegralWithinLimits(histoRinPtBinData[2], rMinGas+eps, rMaxGas-eps,1.);
    Double_t dataStatErrorGas03 = TMath::Sqrt(nconvInRangeData03);
    Double_t dataStatRelErrorGas03 = TMath::Sqrt(nconvInRangeData03)/nconvInRangeData03;


    histoIntegralGasData      = new TH1D("histoIntegralGasData", "histoIntegralGasData", nBinsR, arrayRBins);
    histoIntegralGasData01      = new TH1D("histoIntegralGasData01", "histoIntegralGasData01", nBinsR, arrayRBins);
    histoIntegralGasData02      = new TH1D("histoIntegralGasData02", "histoIntegralGasData02", nBinsR, arrayRBins);
    histoIntegralGasData03      = new TH1D("histoIntegralGasData03", "histoIntegralGasData03", nBinsR, arrayRBins);



    for(Int_t i=1; i<nBinsR+1; i++){
      histoIntegralGasData->SetBinContent(i,nconvInRangeData);
      histoIntegralGasData->SetBinError(i,dataStatErrorGas);

      histoIntegralGasData01->SetBinContent(i,nconvInRangeData01);
      histoIntegralGasData01->SetBinError(i,dataStatErrorGas01);

      histoIntegralGasData02->SetBinContent(i,nconvInRangeData02);
      histoIntegralGasData02->SetBinError(i,dataStatErrorGas02);

      histoIntegralGasData03->SetBinContent(i,nconvInRangeData03);
      histoIntegralGasData03->SetBinError(i,dataStatErrorGas03);
    }

    TH1F* histoIntegralGasDataWide = (TH1F*) histoRData->Clone("histoIntegralGasDataWide");
    for(Int_t i=1; i<histoIntegralGasDataWide->GetNbinsX()+1; i++){
      histoIntegralGasDataWide->SetBinContent(i,nconvInRangeData);
      histoIntegralGasDataWide->SetBinError(i,dataStatErrorGas);
    }

    //new TH1D("histoIntegralGasData", "histoIntegralGasData", nBinsR, arrayRBins);

    histoRDataScaledToGas = (TH1F*) histoRData->Clone("histoRDataScaledToGas");
    histoRDataScaledToGas->Sumw2();
    histoRDataScaledToGas->Divide(histoRDataScaledToGas , histoIntegralGasDataWide,1.0,1.0,"E");

    histoRDataRebin   = (TH1F*)histoRData->Rebin(nBinsR,"histoRDataRebin", arrayRBins);
    histoRDataRebin->Sumw2();
    histoRDataScaledToGasRebin =(TH1F*)histoRDataRebin->Clone("histoRDataScaledToGasRebin");
    histoRDataScaledToGasRebin->Divide(histoRDataRebin,histoIntegralGasData,1.0,1.0,"E");

    histoRinPtBinDataRebin01 = (TH1F*)histoRinPtBinData[0]->Rebin(nBinsR,"histoRinPtBinDataRebin01", arrayRBins);
    histoRinPtBinDataRebin01->Sumw2();
    histoRinPtBinDataScaledToGasRebin01= (TH1F*)histoRinPtBinDataRebin01->Clone("histoRinPtBinDataScaledToGasRebin01");
    histoRinPtBinDataScaledToGasRebin01->Divide(histoRinPtBinDataRebin01,histoIntegralGasData01,1.0,1.0,"E");

    histoRinPtBinDataRebin02 = (TH1F*)histoRinPtBinData[1]->Rebin(nBinsR,"histoRinPtBinDataRebin02", arrayRBins);
    histoRinPtBinDataRebin02->Sumw2();
    histoRinPtBinDataScaledToGasRebin02= (TH1F*)histoRinPtBinDataRebin02->Clone("histoRinPtBinDataScaledToGasRebin02");
    histoRinPtBinDataScaledToGasRebin02->Divide(histoRinPtBinDataRebin02,histoIntegralGasData02,1.0,1.0,"E");


    histoRinPtBinDataRebin03 = (TH1F*)histoRinPtBinData[2]->Rebin(nBinsR,"histoRinPtBinDataRebin03", arrayRBins);
    histoRinPtBinDataRebin03->Sumw2();
    histoRinPtBinDataScaledToGasRebin03= (TH1F*)histoRinPtBinDataRebin03->Clone("histoRinPtBinDataScaledToGasRebin03");
    histoRinPtBinDataScaledToGasRebin03->Divide(histoRinPtBinDataRebin03,histoIntegralGasData03,1.0,1.0,"E");



    cout << "MC scaling to Nconv in gas:" << endl;
    Double_t nconvInRangeMC = CalculateIntegralWithinLimits(histoRMC, rMinGas+eps, rMaxGas-eps, mcGasCorrectionFactor);
    Double_t mcStatErrorGas = TMath::Sqrt(nconvInRangeMC);
    Double_t mcStatRelErrorGas = TMath::Sqrt(nconvInRangeMC)/nconvInRangeMC;

    Double_t nconvInRangeMC01 = CalculateIntegralWithinLimits(histoRinPtBinMC[0], rMinGas+eps, rMaxGas-eps, mcGasCorrectionFactor);
    Double_t mcStatErrorGas01 = TMath::Sqrt(nconvInRangeMC01);
    Double_t mcStatRelErrorGas01 = TMath::Sqrt(nconvInRangeMC01)/nconvInRangeMC01;

    Double_t nconvInRangeMC02 = CalculateIntegralWithinLimits(histoRinPtBinMC[1], rMinGas+eps, rMaxGas-eps, mcGasCorrectionFactor);
    Double_t mcStatErrorGas02 = TMath::Sqrt(nconvInRangeMC02);
    Double_t mcStatRelErrorGas02 = TMath::Sqrt(nconvInRangeMC02)/nconvInRangeMC02;

    Double_t nconvInRangeMC03 = CalculateIntegralWithinLimits(histoRinPtBinMC[2], rMinGas+eps, rMaxGas-eps, mcGasCorrectionFactor);
    Double_t mcStatErrorGas03 = TMath::Sqrt(nconvInRangeMC03);
    Double_t mcStatRelErrorGas03 = TMath::Sqrt(nconvInRangeMC03)/nconvInRangeMC03;

    histoIntegralGasMC      = new TH1D("histoIntegralGasMC", "histoIntegralGasMC", nBinsR, arrayRBins);
    histoIntegralGasMC01      = new TH1D("histoIntegralGasMC01", "histoIntegralGasMC01", nBinsR, arrayRBins);
    histoIntegralGasMC02      = new TH1D("histoIntegralGasMC02", "histoIntegralGasMC02", nBinsR, arrayRBins);
    histoIntegralGasMC03      = new TH1D("histoIntegralGasMC03", "histoIntegralGasMC03", nBinsR, arrayRBins);

    for(Int_t i=1; i<nBinsR+1; i++){
      histoIntegralGasMC->SetBinContent(i,nconvInRangeMC);
      histoIntegralGasMC->SetBinError(i,mcStatErrorGas);

      histoIntegralGasMC01->SetBinContent(i,nconvInRangeMC01);
      histoIntegralGasMC01->SetBinError(i,mcStatErrorGas01);

      histoIntegralGasMC02->SetBinContent(i,nconvInRangeMC02);
      histoIntegralGasMC02->SetBinError(i,mcStatErrorGas02);

      histoIntegralGasMC03->SetBinContent(i,nconvInRangeMC03);
      histoIntegralGasMC03->SetBinError(i,mcStatErrorGas03);
    }

    TH1F* histoIntegralGasMCWide = (TH1F*) histoRData->Clone("histoIntegralGasMCWide");
    for(Int_t i=1; i<histoIntegralGasMCWide->GetNbinsX()+1; i++){
      histoIntegralGasMCWide->SetBinContent(i,nconvInRangeMC);
      histoIntegralGasMCWide->SetBinError(i,mcStatErrorGas);
    }

    histoRMCScaledToGas = (TH1F*) histoRMC->Clone("histoRMCScaledToGas");
    histoRMCScaledToGas->Sumw2();
    histoRMCScaledToGas->Divide(histoRMCScaledToGas , histoIntegralGasMCWide,1.0,1.0,"E");


    histoRMCRebin   = (TH1F*)histoRMC->Rebin(nBinsR,"histoRMCRebin", arrayRBins);
    histoRMCRebin->Sumw2();

    histoRMCScaledToGasRebin =(TH1F*)histoRMCRebin->Clone("histoRMCScaledToGasRebin");
    histoRMCScaledToGasRebin->Divide(histoRMCRebin,histoIntegralGasMC,1.0,1.0,"E");

    histoRinPtBinMCRebin01 = (TH1F*)histoRinPtBinMC[0]->Rebin(nBinsR,"histoRinPtBinMCRebin01", arrayRBins);
    histoRinPtBinMCRebin01->Sumw2();
    histoRinPtBinMCScaledToGasRebin01= (TH1F*)histoRinPtBinMCRebin01->Clone("histoRinPtBinMCScaledToGasRebin01");
    histoRinPtBinMCScaledToGasRebin01->Divide(histoRinPtBinMCRebin01,histoIntegralGasMC01,1.0,1.0,"E");

    histoRinPtBinMCRebin02 = (TH1F*)histoRinPtBinMC[1]->Rebin(nBinsR,"histoRinPtBinMCRebin02", arrayRBins);
    histoRinPtBinMCRebin02->Sumw2();
    histoRinPtBinMCScaledToGasRebin02= (TH1F*)histoRinPtBinMCRebin02->Clone("histoRinPtBinMCScaledToGasRebin02");
    histoRinPtBinMCScaledToGasRebin02->Divide(histoRinPtBinMCRebin02,histoIntegralGasMC02,1.0,1.0,"E");


    histoRinPtBinMCRebin03 = (TH1F*)histoRinPtBinMC[2]->Rebin(nBinsR,"histoRinPtBinMCRebin03", arrayRBins);
    histoRinPtBinMCRebin03->Sumw2();
    histoRinPtBinMCScaledToGasRebin03= (TH1F*)histoRinPtBinMCRebin03->Clone("histoRinPtBinMCScaledToGasRebin03");
    histoRinPtBinMCScaledToGasRebin03->Divide(histoRinPtBinMCRebin03,histoIntegralGasMC03,1.0,1.0,"E");





    histoDataMCRatioRScaledToGas = (TH1F*)histoRDataScaledToGasRebin->Clone("histoDataMCRatioRScaledToGas");
    histoDataMCRatioRScaledToGas->Divide(histoRDataScaledToGasRebin,histoRMCScaledToGasRebin,1.,1.,"E");

    cout << "scaling to Nconv in gas:pT>0.15" << endl;
    histoDataMCRatioRinPtBinScaledToGas01  = (TH1F*)histoRinPtBinDataScaledToGasRebin01->Clone("histoDataMCRatioRinPtBinScaledToGas01");
    histoDataMCRatioRinPtBinScaledToGas01->Divide(histoRinPtBinDataScaledToGasRebin01 ,histoRinPtBinMCScaledToGasRebin01 ,1.,1.,"E");

    cout << "scaling to Nconv in gas:pT>0.3" << endl;
    histoDataMCRatioRinPtBinScaledToGas02  = (TH1F*)histoRinPtBinDataScaledToGasRebin02->Clone("histoDataMCRatioRinPtBinScaledToGas02");
    histoDataMCRatioRinPtBinScaledToGas02->Divide(histoRinPtBinDataScaledToGasRebin02 ,histoRinPtBinMCScaledToGasRebin02 ,1.,1.,"E");

    cout << "scaling to Nconv in gas:pT>0.4" << endl;
    histoDataMCRatioRinPtBinScaledToGas03  = (TH1F*)histoRinPtBinDataScaledToGasRebin03->Clone("histoDataMCRatioRinPtBinScaledToGas03");
    histoDataMCRatioRinPtBinScaledToGas03->Divide(histoRinPtBinDataScaledToGasRebin03 ,histoRinPtBinMCScaledToGasRebin03 ,1.,1.,"E");


   //  for(Int_t i=1; i<nBinsR+1; i++){
   //    cout<<"Test-0::"<<  histoRDataRebin->GetBinCenter(i) << " " << histoRDataScaledToGasRebin->GetBinCenter(i)<< endl;
   //    cout<<"Test-00::"<< histoIntegralGasData->GetBinContent(i) << " "<< histoIntegralGasData->GetBinError(i)<< endl;
   //    cout<<"Test-1::"<< histoRDataRebin->GetBinContent(i)  << "  "<< histoRDataRebin->GetBinError(i)<< "   "<<  histoRDataRebin->GetBinError(i)/ histoRDataRebin ->GetBinContent(i)  << endl;
   //    cout<<"Test-2::"<< histoRDataScaledToGasRebin->GetBinContent(i) << "  "<< histoRDataScaledToGasRebin->GetBinError(i) << "  "<<  histoRDataScaledToGasRebin->GetBinError(i)/ histoRDataScaledToGasRebin->GetBinContent(i) <<endl;
   //    cout<<"Test-3::"<< histoRMCScaledToGasRebin->GetBinContent(i) << "  "<< histoRMCScaledToGasRebin->GetBinError(i) << "  "<<  histoRMCScaledToGasRebin->GetBinError(i)/ histoRMCScaledToGasRebin->GetBinContent(i) <<endl;
   //    cout<<"Test-4::"<< histoDataMCRatioRScaledToGas-> GetBinContent(i) << "  "<<histoDataMCRatioRScaledToGas->GetBinError(i) << "  "<< histoDataMCRatioRScaledToGas-> GetBinError(i)/histoDataMCRatioRScaledToGas->GetBinContent(i) << endl;
   // }


    //-AM    Normalization included here
    GammaScalingHistogramm(histoRTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPtTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPtMC,normFactorReconstMC);
    GammaScalingHistogramm(histoEtaMC,normFactorReconstMC);


    GammaScalingHistogramm(histoRMCRTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPtMCPtTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoConvRMC,normFactorReconstMC);
    GammaScalingHistogramm(histoConvPtMC,normFactorReconstMC);


    GammaScalingHistogramm(histoPi0DalRTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoEtaDalRTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoCombRTrueMC,normFactorReconstMC);



    GammaScalingHistogramm(histoEtaData,normFactorReconstData);
    GammaScalingHistogramm(histoPtData,normFactorReconstData);
    GammaScalingHistogramm(histoRData,normFactorReconstData);
    GammaScalingHistogramm(histoRMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRPhiData,normFactorReconstData);
    GammaScalingHistogramm(histoRPhiMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRPhiTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRZData,normFactorReconstData);
    GammaScalingHistogramm(histoRZMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRZTrueMC,normFactorReconstMC);


    GammaScalingHistogramm(histoAsymPDataLowP,normFactorReconstData);
    GammaScalingHistogramm(histoAsymPDataHighP,normFactorReconstData);

    GammaScalingHistogramm(histoAsymPTrueMCLowP,normFactorReconstMC);
    GammaScalingHistogramm(histoAsymPTrueMCHighP,normFactorReconstMC);


    for(Int_t i=0; i<3; i++){
      GammaScalingHistogramm(histoPtinRBinMC[i],normFactorReconstMC);
      GammaScalingHistogramm(histoPtinRBinTrueMC[i],normFactorReconstMC);
    }


    // - AM calculation of efficiencies
    histoEffiR = (TH1F*)histoRMCRTrueMC->Clone("histoEffiR");
    histoEffiR->Sumw2();
    histoEffiR->Divide(histoRMCRTrueMC,histoConvRMC,1.,1.,"B");

    Int_t nBinsPtNew=63;
    Double_t  arrPtBins[nBinsPtNew];
    for(Int_t i=0;i<nBinsPtNew+1;i++){
        if(i<20){       // 0.05    , 0.-1.
            arrPtBins[i]=0.05*i;
        }else if(i<40){     // 0.1  1-3.
            arrPtBins[i]=1.+0.1*(i-20);
        }else if(i<45){     // 0.2  3.-4.
            arrPtBins[i]=3.+0.2*(i-40);
        }else if(i<49){     // 0.5  4.-6.
            arrPtBins[i]=4.+0.5*(i-45);
        }else if(i<nBinsPtNew+1){     // 0.1  1-3.
            arrPtBins[i]=6.+1.*(i-49);
        }
        // cout << "i, pT(i)::"<< i  << "   "<<  arrPtBins[i]<< endl;
    }
    TH1F * histoPtMCPtTrueMCRebin = (TH1F *) histoPtMCPtTrueMC->Rebin(nBinsPtNew,"histoPtMCPtTrueMCRebin",arrPtBins);
    TH1F * histoConvPtMCRebin = (TH1F *) histoConvPtMC->Rebin(nBinsPtNew,"histoConvPtMCRebin",arrPtBins);
    histoEffiPt                 = (TH1F*)histoPtMCPtTrueMCRebin->Clone("histoEffiPt");
    histoEffiPt->Sumw2();
    histoEffiPt->Divide(histoPtMCPtTrueMCRebin,histoConvPtMCRebin,1.,1.,"B");

    TH1F *histoPhiTrueMC       = (TH1F*)histoRPhiTrueMC->ProjectionX("histoPhiTrueMC");
    histoEffiPhi = (TH1F*)histoPhiTrueMC->Clone("histoEffiR");
    histoEffiPhi->Sumw2();
    histoEffiPhi->Divide(histoPhiTrueMC,histoConvPhiMC,1.,1.,"B");

    TH1F *histoEtaTrueMC       = (TH1F*)histoREtaTrueMC->ProjectionX("histoEtaTrueMC");
    histoEffiEta = (TH1F*)histoEtaTrueMC->Clone("histoEffiR");
    histoEffiEta->Sumw2();
    histoEffiEta->Divide(histoEtaTrueMC,histoConvEtaMC,1.,1.,"B");


    //- AM  calculation of quantities only related to MC.
    histoPurityR                 = (TH1F*)histoRTrueMC->Clone("histoPurityR");
    histoPurityR->Sumw2();
    histoPurityR->Divide(histoRTrueMC,histoRMC,1.,1.,"B");

    histoPurityPt                = (TH1F*)histoPtTrueMC->Clone("histoPurityPt");
    histoPurityPt->Sumw2();
    histoPurityPt->Divide(histoPtTrueMC,histoPtMC,1.,1.,"B");

    TH1F* histoPt5cmMC = (TH1F*)histoPtinRBinMC[1]->Clone("histoPt5cmMC");
    histoPt5cmMC->Sumw2();
    histoPt5cmMC->Add(histoPtinRBinMC[2]);

    TH1F* histoPt5cmTrueMC = (TH1F*)histoPtinRBinTrueMC[1]->Clone("histoPt5cmTrueMC");
    histoPt5cmTrueMC->Sumw2();
    histoPt5cmTrueMC->Add(histoPtinRBinTrueMC[2]);

    histoPurityPt5cm = (TH1F*)histoPt5cmTrueMC->Clone("histoPurityPt5cm");
    histoPurityPt5cm->Sumw2();
    histoPurityPt5cm->Divide(histoPt5cmTrueMC,histoPt5cmMC,1.,1.,"B");

    //- AM  Start comparison Data-MC without scaling to gas. Normalization to number of events should be inserted here

    histoDataMCRatioR            = (TH1F*)histoRData->Clone("histoDataMCRatioR");
    histoDataMCRatioR->Divide(histoRData,histoRMC,1.,1.,"E");


    for(Int_t iR = 0; iR < nBinsR; iR++){
      histoPhiInRData[iR]        = (TH1D*)histoRPhiData->ProjectionX(Form("histoPhiInRData_%i",iR),
                                                                     histoRPhiData->GetYaxis()->FindBin(arrayRBins[iR]+eps),
                                                                     histoRPhiData->GetYaxis()->FindBin(arrayRBins[iR+1]-eps),"e");
      ConvGammaRebinWithBinCorrection(histoPhiInRData[iR],rebinPhiPlots);

      histoPhiInRMC[iR]          = (TH1D*)histoRPhiMC->ProjectionX(Form("histoPhiInRMC_%i",iR),
                                                                   histoRPhiMC->GetYaxis()->FindBin(arrayRBins[iR]+eps),
                                                                   histoRPhiMC->GetYaxis()->FindBin(arrayRBins[iR+1]-eps),"e");
      ConvGammaRebinWithBinCorrection(histoPhiInRMC[iR],rebinPhiPlots);

      histoPhiInRTrueMC[iR]      = (TH1D*)histoRPhiTrueMC->ProjectionX(Form("histoPhiInRTrueMC_%i",iR),
                                                                       histoRPhiTrueMC->GetYaxis()->FindBin(arrayRBins[iR]+eps),
                                                                       histoRPhiTrueMC->GetYaxis()->FindBin(arrayRBins[iR+1]-eps),"e");
      ConvGammaRebinWithBinCorrection(histoPhiInRTrueMC[iR],rebinPhiPlots);

      histoDataMCRatioPhiInR[iR] = (TH1D*)histoPhiInRData[iR]->Clone(Form("histoDataMCRatioPhiInR_%02d",iR));
      histoDataMCRatioPhiInR[iR]->Divide(histoPhiInRData[iR],histoPhiInRMC[iR]);


      if (iR > 9) rebinZPlots = 4;
      histoZInRData[iR]          =  (TH1D*)histoRZData->ProjectionX( Form("histoZInRData_%i",iR),
                                                                     histoRZData->GetYaxis()->FindBin(arrayRBins[iR]+eps),
                                                                     histoRZData->GetYaxis()->FindBin(arrayRBins[iR+1]-eps),"e");
      ConvGammaRebinWithBinCorrection(histoZInRData[iR],rebinZPlots);

      histoZInRMC[iR]            =  (TH1D*)histoRZMC->ProjectionX( Form("histoZInRMC_%i",iR),
                                                                   histoRZMC->GetYaxis()->FindBin(arrayRBins[iR]+eps),
                                                                   histoRZMC->GetYaxis()->FindBin(arrayRBins[iR+1]-eps),"e");
      ConvGammaRebinWithBinCorrection(histoZInRMC[iR],rebinZPlots);

      histoZInRTrueMC[iR]        =  (TH1D*)histoRZTrueMC->ProjectionX( Form("histoZInRTrueMC_%i",iR),
                                                                       histoRZTrueMC->GetYaxis()->FindBin(arrayRBins[iR]+eps),
                                                                       histoRZTrueMC->GetYaxis()->FindBin(arrayRBins[iR+1]-eps),"e");
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
    if ( (optionEnergy.CompareTo("13TeV") == 0) || (optionEnergy.CompareTo("5TeV") == 0)  ) {
       histoDummyNTracks->GetXaxis()->SetRangeUser(0.,100);
    }
    histoDummyNTracks->DrawCopy();

    DrawGammaSetMarker(histoGoodESDTracksData, 20, markerSize, colorData, colorData);
    histoGoodESDTracksData->Draw("same,hist");
    DrawGammaSetMarker(histoGoodESDTracksMC, 20, markerSize, colorMC, colorMC);
    histoGoodESDTracksMC->Draw("same,hist");

    TLegend* legend = GetAndSetLegend(0.75,0.75,2);
    legend->AddEntry(histoGoodESDTracksData,"Data","l");
    legend->AddEntry(histoGoodESDTracksMC,"MC","l");
    legend->Draw();

    canvasNTracks->Print(Form("%s/NGoodESDTracks%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    //______________________________________ Pt _________________________________________
    TCanvas * canvasPt = new TCanvas("canvasPt","",1200,1000);
    DrawGammaCanvasSettings( canvasPt, 0.1, 0.03, 0.05, 0.09);
    //canvasNTracks->Divide(2,2);
    canvasPt->SetLogy(1);
    TH2F * histoDummyPt = new TH2F("histoDummyPt","histoDummyPt",1000,0.,15.,1000,1.e-10,1.e-2);
    SetStyleHistoTH2ForGraphs(histoDummyPt, "pT (GeV)","Counts/(Nev. MeanMult)", 0.035,0.04,0.035,0.04,1.,1.);
    histoDummyPt->DrawCopy();

    DrawGammaSetMarker(histoPtData, 20, markerSize, colorData, colorData);
    histoPtData->Draw("same,hist");
    DrawGammaSetMarker(histoPtMC, 20, markerSize, colorMC, colorMC);
    histoPtMC->Draw("same,hist");

    TLegend* legendPt = GetAndSetLegend(0.75,0.75,2);
    legendPt->AddEntry(histoPtData,"Data","l");
    legendPt->AddEntry(histoPtMC,"MC","l");
    legendPt->Draw();

    canvasPt->Print(Form("%s/pT%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

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
    if ((optionEnergy.CompareTo("13TeV") == 0) || (optionEnergy.CompareTo("5TeV") == 0)) {
      histoDummyTwoPanelsUp->GetYaxis()->SetRangeUser(-0.1,2.5);
    }
    SetStyleHistoTH2ForGraphs(histoDummyTwoPanelsUp, "R (cm)","Counts", 0.9*textsizeLabelsUp, textsizeLabelsUp,0.9*textsizeLabelsUp,textsizeLabelsUp, 1,0.15/(textsizeFacUp*marginXRatio));
    //         histoDummyTwoPanelsUp->GetYaxis()->SetLabelOffset(+0.01);
    TH2F *histoDummyTwoPanelsDown =  new TH2F("histoDummyTwoPanelsDown","histoDummyTwoPanelsDown",1000,0.,180.,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histoDummyTwoPanelsDown, "R (cm)","#frac{Data}{MC} ", 0.9*textsizeLabelsDown, textsizeLabelsDown,0.9*textsizeLabelsDown,textsizeLabelsDown, 1,0.15/(textsizeFacDown*marginXRatio));
    //         histoDummyTwoPanelsDown->GetYaxis()->SetLabelOffset(+0.01);
    histoDummyTwoPanelsDown->GetYaxis()->SetRangeUser(0.3,1.55);
    if ((optionEnergy.CompareTo("13TeV") == 0 ) || (optionEnergy.CompareTo("5TeV") == 0)){
      histoDummyTwoPanelsDown->GetYaxis()->SetRangeUser(0.8,1.55);
    }

    padUpper->cd();
    histoDummyTwoPanelsUp->DrawCopy();


    DrawGammaSetMarker(histoRDataScaledToGas, 20, markerSize, colorData, colorData);
    histoRDataScaledToGas->Draw("same,hist");

    DrawGammaSetMarker(histoRMCScaledToGas, 20, markerSize, colorMC, colorMC);
    histoRMCScaledToGas->Draw("same,hist");

    TLegend* legenUpPanel = GetAndSetLegend2(0.7, 0.8-(2*0.9*textsizeLabelsUp), 0.9, 0.8, textSizeLabels);
    legenUpPanel->AddEntry(histoRDataScaledToGas,"Data","l");
    legenUpPanel->AddEntry(histoRMCScaledToGas,"MC","l");
    legenUpPanel->Draw();

    padLower->cd();
    histoDummyTwoPanelsDown->DrawCopy();

    DrawGammaLines(0.,180,1., 1.,1.,kGray,1);
    DrawGammaLines(0.,180,1.15, 1.15,1.,kGray,4);
    DrawGammaLines(0.,180,1.1, 1.1,1.,kGray,2);
    DrawGammaLines(0.,180,1.05, 1.05,1.,kGray,3);
    DrawGammaLines(0.,180,0.95, 0.95,1.,kGray,3);

    DrawGammaSetMarker(histoDataMCRatioRScaledToGas, 20, markerSize, colorData, colorData);
    DrawGammaSetMarker(histoDataMCRatioRinPtBinScaledToGas01, 25, markerSize, kViolet, kViolet);
    DrawGammaSetMarker(histoDataMCRatioRinPtBinScaledToGas02, 25, markerSize, kBlue, kBlue);
    DrawGammaSetMarker(histoDataMCRatioRinPtBinScaledToGas03, 25, markerSize, kGreen, kGreen);

    histoDataMCRatioRScaledToGas->Draw("same,pE");
    histoDataMCRatioRinPtBinScaledToGas01->Draw("same,pE");
    histoDataMCRatioRinPtBinScaledToGas02->Draw("same,pE");
    histoDataMCRatioRinPtBinScaledToGas03->Draw("same,pE");

    TLegend* legenLowPanel = GetAndSetLegend2(0.5, 0.85-(4*0.9*textsizeLabelsDown), 0.7, 0.85, textSizeLabels);
    TString legendWithPeriod = Form("From proj, Period %s",optionPeriod.Data());

    legenLowPanel->AddEntry(histoDataMCRatioRScaledToGas,legendWithPeriod.Data(),"lp");
    legenLowPanel->AddEntry(histoDataMCRatioRinPtBinScaledToGas01,"pT>0.15 GeV","lp");
    legenLowPanel->AddEntry(histoDataMCRatioRinPtBinScaledToGas02,"pT>0.3 GeV","lp");
    legenLowPanel->AddEntry(histoDataMCRatioRinPtBinScaledToGas03,"pT>0.4 GeV","lp");
    legenLowPanel->Draw();



    histoDummyTwoPanelsDown->Draw("same,axis");
    canvasTwoPanels->Print(Form("%s/PhotonConvRScaledToGas%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    TCanvas * canvasPhiR = new TCanvas("canvasPhiR","",480,440);  // gives the page size
    DrawGammaCanvasSettings( canvasPhiR, 0.09, 0.105, 0.02, 0.085);
    canvasPhiR->SetLogz(1);
    canvasPhiR->cd();

    SetStyleHistoTH2ForGraphs(  histoRPhiData, "#varphi (rad)","R (cm)", 0.035, 0.04,
				0.035, 0.04, 0.9, 1.0, 510, 510);
    histoRPhiData->GetYaxis()->SetRangeUser(0,200.);
    histoRPhiData->GetXaxis()->SetRangeUser(0.,2*TMath::Pi());

    histoRPhiData->Draw("colz");

    canvasPhiR->Print(Form("%s/PhotonConvRPhi%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasEtaR = new TCanvas("canvasEtaR","",480,440);  // gives the page size
    DrawGammaCanvasSettings( canvasEtaR, 0.09, 0.105, 0.02, 0.085);
    canvasEtaR->SetLogz(1);
    canvasEtaR->cd();

    SetStyleHistoTH2ForGraphs(  histoREtaData, "#eta","R (cm)", 0.035, 0.04,
				0.035, 0.04, 0.9, 1.0, 510, 510);
    histoREtaData->GetYaxis()->SetRangeUser(0,200.);
    histoREtaData->GetXaxis()->SetRangeUser(-2., 2.);

    histoREtaData->Draw("colz");

    canvasEtaR->Print(Form("%s/PhotonConvREta%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasZR = new TCanvas("canvasZR","",480,440);  // gives the page size
    DrawGammaCanvasSettings( canvasZR, 0.09, 0.105, 0.02, 0.085);
    canvasZR->SetLogz(1);
    canvasZR->cd();

    SetStyleHistoTH2ForGraphs(  histoRZData, "Z (cm)","R (cm)", 0.035, 0.04,
				0.035, 0.04, 0.9, 1.0, 510, 510);
    histoRZData->GetYaxis()->SetRangeUser(0,200.);
    histoRZData->GetXaxis()->SetRangeUser(-180., 180.);

    histoRZData->Draw("colz");

    canvasZR->Print(Form("%s/PhotonConvRZ%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasPtR = new TCanvas("canvasPtR","",480,440);  // gives the page size
    DrawGammaCanvasSettings( canvasPtR, 0.09, 0.105, 0.02, 0.085);
    canvasPtR->SetLogz(1);
    canvasPtR->cd();

    SetStyleHistoTH2ForGraphs(  histoRPtData, "#it{p}_{T} (GeV/#it{c})","R (cm)", 0.035, 0.04,
				0.035, 0.04, 0.9, 1.0, 510, 510);
    histoRPtData->GetYaxis()->SetRangeUser(0,200.);
    histoRPtData->GetXaxis()->SetRangeUser(0., 20.);

    histoRPtData->Draw("colz");

    canvasPtR->Print(Form("%s/PhotonConvRPt%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasElecdEdxR = new TCanvas("canvasElecdEdxR","",480,440);  // gives the page size
    DrawGammaCanvasSettings( canvasElecdEdxR, 0.09, 0.105, 0.02, 0.085);
    canvasElecdEdxR->SetLogz(1);
    canvasElecdEdxR->cd();

    SetStyleHistoTH2ForGraphs(  histoRElecdEdxData, "electron dEdx","R (cm)", 0.035, 0.04,
				0.035, 0.04, 0.9, 1.0, 510, 510);
    histoRElecdEdxData->GetYaxis()->SetRangeUser(0,200.);
    histoRElecdEdxData->GetXaxis()->SetRangeUser(0.,200.);

    histoRElecdEdxData->Draw("colz");

    canvasElecdEdxR->Print(Form("%s/PhotonConvRElecdEdx%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasElecNSdEdxR = new TCanvas("canvasElecNSdEdxR","",480,440);  // gives the page size
    DrawGammaCanvasSettings( canvasElecNSdEdxR, 0.09, 0.105, 0.02, 0.085);
    canvasElecNSdEdxR->SetLogz(1);
    canvasElecNSdEdxR->cd();

    SetStyleHistoTH2ForGraphs(  histoRElecNSdEdxData, "electron n#sigma dEdx","R (cm)", 0.035, 0.04,
				0.035, 0.04, 0.9, 1.0, 510, 510);
    histoRElecNSdEdxData->GetYaxis()->SetRangeUser(0,200.);
    histoRElecNSdEdxData->GetXaxis()->SetRangeUser(-10., 10.);

    histoRElecNSdEdxData->Draw("colz");

    canvasElecNSdEdxR->Print(Form("%s/PhotonConvRElecNSdEdx%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


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
      DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.,kGray);

    DrawGammaSetMarker(histoRData, 20, markerSize, colorData, colorData);
    histoRData->Draw("same,hist");
    DrawGammaSetMarker(histoRMC, 20, markerSize, colorMC, colorMC);
    histoRMC->Draw("same,hist");
    DrawGammaSetMarker(histoRTrueMC, 20, markerSize, colorTrueMC, colorTrueMC);
    histoRTrueMC->Draw("same,hist");
    //histo not defined . To be checked
    DrawGammaSetMarker(histoCombRTrueMC, 20, markerSize, colorTrueCombMC, colorTrueCombMC);
    histoCombRTrueMC->DrawCopy("same,hist");
    DrawGammaSetMarker(histoPi0DalRTrueMC, 20, markerSize, kBlue-9, kBlue-9);
    histoPi0DalRTrueMC->SetFillColor(kBlue-9);
    histoPi0DalRTrueMC->SetFillStyle(3244);
    histoPi0DalRTrueMC->DrawCopy("same,hist");
    DrawGammaSetMarker(histoEtaDalRTrueMC, 20, markerSize, kBlue-3, kBlue-3);
    histoEtaDalRTrueMC->SetFillColor(kBlue-3);
    histoEtaDalRTrueMC->SetFillStyle(3002);
    histoEtaDalRTrueMC->DrawCopy("same,hist");

    TLegend* legenRdistrib = GetAndSetLegend2(0.6, 0.93-(6*0.85*textsizeLabelsUp), 0.9, 0.93, textSizeLabels);
    legenRdistrib->AddEntry(histoRData,"Data","l");
    legenRdistrib->AddEntry(histoRMC,"MC","l");
    legenRdistrib->AddEntry(histoRTrueMC,"True MC","l");
    legenRdistrib->AddEntry(histoCombRTrueMC,"True MC comb.","l");
    legenRdistrib->AddEntry(histoPi0DalRTrueMC,"True MC #pi^{0} Dal.","lf");
    legenRdistrib->AddEntry(histoEtaDalRTrueMC,"True MC #eta Dal.","fl");
    legenRdistrib->Draw();

    histoDummy4PanelsUp->Draw("axis,same");
    padRDistribLowerLeft->cd();
    histoDummy4PanelsDown->DrawCopy();

    for(Int_t iR = 1; iR < nBinsR; iR++)
      DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.,kGray);
    DrawGammaLines(0.,180,1., 1.,1.,kGray);

    DrawGammaSetMarker(histoDataMCRatioR, 20, markerSize, colorData, colorData);
    histoDataMCRatioR->Draw("same,histo");

    histoDummy4PanelsDown->Draw("axis,same");
    padRDistribUpperRight->cd();
    padRDistribUpperRight->SetLogy();
    histoDummy4PanelsUp->DrawCopy();

    for(Int_t iR = 1; iR < nBinsR; iR++)
      DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.,kGray);

    histoRData->Draw("same,histo");
    histoRMC->Draw("same,histo");
    histoRTrueMC->Draw("same,histo");

    for(Int_t i=0; i<6; i++){
      Int_t colorDataPt=i+1;
      DrawGammaSetMarker(histoRinPtBinData[i], 20, markerSize, colorDataPt, colorDataPt);
      histoRinPtBinData[i]->SetLineStyle(2);
      histoRinPtBinData[i]->DrawCopy("same,histo");

      DrawGammaSetMarker(histoRinPtBinMC[i], 20, markerSize, colorMC, colorMC);
      histoRinPtBinMC[i]->SetLineStyle(2);
      histoRinPtBinMC[i]->DrawCopy("same,histo");

      DrawGammaSetMarker(histoRinPtBinTrueMC[i], 20, markerSize, colorTrueMC, colorTrueMC);
      histoRinPtBinTrueMC[i]->SetLineStyle(2);
      histoRinPtBinTrueMC[i]->DrawCopy("same,histo");
    }
    // AM- pT limits wrong !!!!
    TLegend* legenRdistribPtcuts = GetAndSetLegend2(0.05, 0.35-(6.2*0.9*textsizeLabelsUp), 0.35, 0.35, textSizeLabels);
    legenRdistribPtcuts->AddEntry(histoRData,"Data","l");
    legenRdistribPtcuts->AddEntry(histoRMC,"MC","l");
    legenRdistribPtcuts->AddEntry(histoRTrueMC,"True MC","l");
    legenRdistribPtcuts->AddEntry(histoRinPtBinData[1],"Data, 0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","l");
    legenRdistribPtcuts->AddEntry(histoRinPtBinMC[1],"MC, 0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","lf");
    legenRdistribPtcuts->AddEntry(histoRinPtBinTrueMC[1],"True MC, 0.4 < #it{p}_{T} < 1.5 GeV/#it{c}","fl");
    // legenRdistribPtcuts->Draw();

    histoDummy4PanelsUp->Draw("axis,same");
    padRDistribLowerRight->cd();

    histoDummy4PanelsDown->DrawCopy();
     legenRdistribPtcuts->Draw();
    for(Int_t iR = 1; iR < nBinsR; iR++)
      DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.1,kGray);
    DrawGammaLines(0.,180,1., 1.,1.,kGray);

    histoDataMCRatioR->Draw("same,histo");

    TLegend* legenRdistribPtcutsRatio = GetAndSetLegend2(0.05, 0.27-(2*0.9*textsizeLabelsDown), 0.35, 0.27, textSizeLabels);
    legenRdistribPtcutsRatio->AddEntry(histoDataMCRatioR,"No #it{p}_{T} cut","l");
    legenRdistribPtcutsRatio->Draw();

    TLatex *latexPeriod = new TLatex(0.05,0.85,Form("Period %s (%s)",optionPeriod.Data(), fCutSelection.Data()));
    SetStyleTLatex( latexPeriod, sizeTextNameBins,2);
    latexPeriod->Draw();


    histoDummy4PanelsDown->Draw("axis,same");
    canvasPhotonR->Print(Form("%s/PhotonConvR%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    //_____________________________________________________________________________________________________________
    TCanvas* canvasPhotonR_siglepad = new TCanvas("canvasPhotonR_siglepad","",0,0,1300,1000);
    DrawGammaCanvasSettings( canvasPhotonR_siglepad,  0.08, 0.02, 0.02, 0.09);
    canvasPhotonR_siglepad->SetLogy();

    TH2F *histoDummySinglePad = new TH2F("histoDummySinglePad","histoDummySinglePad",1000,0.,184.,1000,1.5e-9,1.);
    SetStyleHistoTH2ForGraphs(histoDummySinglePad, "R (cm)","Counts", 0.035 ,0.04, 0.035,0.04, 0.9, 1.);
    histoDummySinglePad->GetYaxis()->SetRangeUser(1.5e-9,2e-3);
    histoDummySinglePad->Draw("copy");

        for(Int_t iR = 1; iR < nBinsR; iR++) DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.,kGray);

        histoRData->SetLineWidth(2);
        histoRData->Draw("same,hist");
        histoRMC->SetLineWidth(2);
        histoRMC->Draw("same,hist");
        histoRTrueMC->SetLineWidth(2);
        histoRTrueMC->Draw("same,hist");
        histoCombRTrueMC->SetLineWidth(2);
        histoCombRTrueMC->DrawCopy("same,hist");
        histoPi0DalRTrueMC->SetLineWidth(2);
        histoPi0DalRTrueMC->DrawCopy("same,hist");
        histoEtaDalRTrueMC->SetLineWidth(2);
        histoEtaDalRTrueMC->DrawCopy("same,hist");

        legenRdistrib->Draw();

    histoDummySinglePad->Draw("axis,same");
    canvasPhotonR_siglepad->Print(Form("%s/PhotonConvRsingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

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
        histoDummyPurityPt->GetYaxis()->SetRangeUser(0.6,1.1);

    padUpperPurity->cd();
    histoDummyPurityR->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++){
            DrawGammaLines(arrayRBins[iR],arrayRBins[iR],0.,1.2,1.,kGray,2);
        }
        DrawGammaLines(0.,180,1., 1.,1.,kGray+1);
        DrawGammaSetMarker(histoPurityR, 20, markerSize, colorData, colorData);

        histoPurityR->Draw("same,hist");

        TLegend* legenPurity = GetAndSetLegend2(0.15, 0.9-(1*0.9*textsizeLabelsUppurity), 0.7, 0.9, textSizeLabels);
        legenPurity->AddEntry((TObject*)0,Form("Cut: %s",fCutSelectionRead.Data()),"");
//         legenPurity->Draw();

    histoDummyPurityR->Draw("same,axis");

    padLowerPurity->cd();
    histoDummyPurityPt->DrawCopy();

        DrawGammaLines(0.,8.,1., 1.,1.1,kGray+1);
        DrawGammaSetMarker(histoPurityPt, 20, markerSize, colorData, colorData);
        histoPurityPt->Draw("same,hist");
        DrawGammaSetMarker(histoPurityPt5cm, 20, markerSize, kGreen+2, kGreen+2);
        histoPurityPt5cm->Draw("same,hist");

        TLegend* legenPurityPt = GetAndSetLegend2(0.13, 0.2+(2*0.9*textsizeLabelsUppurity), 0.5, 0.2, textSizeLabels);
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

        DrawGammaSetMarker(histoDCAData, 20, markerSize, colorData, colorData);
        histoDCAData->Draw("same,hist");
		DrawGammaSetMarker(histoDCAMC, 20, markerSize, colorMC, colorMC);
        histoDCAMC->Draw("same,hist");

    histoDummyDCA->Draw("axis,same");
    padPhotonCharLowerLeft->cd();
    padPhotonCharLowerLeft->SetLogy();
    TH2F * histoDummyMass = new TH2F("histoDummyMass","histoDummyMass",1000,0.,0.2,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyMass, "Mass (GeV)","Counts ", 0.85*textsizeLabelsDown, textsizeLabelsDown,0.85*textsizeLabelsDown,textsizeLabelsDown, 0.9,0.92);
    histoDummyMass->GetYaxis()->SetRangeUser(1.e-6,2e-3);
    if ( (optionEnergy.CompareTo("13TeV") == 0) || (optionEnergy.CompareTo("5TeV") == 0)  ) {
       histoDummyMass->GetYaxis()->SetRangeUser(1.e-6,2e-2);
    }


    histoDummyMass->DrawCopy();

        DrawGammaSetMarker(histoMassData, 20, markerSize, colorData, colorData);
        histoMassData->Draw("same,hist");
		DrawGammaSetMarker(histoMassMC, 20, markerSize, colorMC, colorMC);
        histoMassMC->Draw("same,hist");

    histoDummyMass->Draw("axis,same");
    padPhotonCharUpperRight->cd();
    padPhotonCharUpperRight->SetLogy();
    TH2F *histoDummyPsiPair = new TH2F("histoDummyPsiPair","histoDummyPsiPair",1000,0.,0.3,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyPsiPair, "#Psi_{pair}","Counts", 0.85*textsizeLabelsUp, textsizeLabelsUp,0.85*textsizeLabelsUp,textsizeLabelsUp, 0.9,0.91);
    histoDummyPsiPair->GetYaxis()->SetRangeUser(1.e-6,2.e-2);
    histoDummyPsiPair->DrawCopy();

        DrawGammaSetMarker(histoPsiPairData, 20, markerSize, colorData, colorData);
        histoPsiPairData->Draw("same,hist");
		DrawGammaSetMarker(histoPsiPairMC, 20, markerSize, colorMC, colorMC);
        histoPsiPairMC->Draw("same,hist");

    histoDummyPsiPair->Draw("axis,same");
    padPhotonCharLowerRight->cd();
    padPhotonCharLowerRight->SetLogy();
    TH2F * histoDummyChi2 = new TH2F("histoDummyChi2","histoDummyChi2",1000,0.,40.,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyChi2, "#chi^{2}","Counts ", 0.85*textsizeLabelsDown, textsizeLabelsDown,0.85*textsizeLabelsDown,textsizeLabelsDown, 0.9,0.92);
    histoDummyChi2->GetYaxis()->SetRangeUser(1.e-6,2.e-2);
    histoDummyChi2->DrawCopy();

        DrawGammaSetMarker(histoChi2Data, 20, markerSize, colorData, colorData);
        histoChi2Data->Draw("same,hist");
		DrawGammaSetMarker(histoChi2MC, 20, markerSize, colorMC, colorMC);
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

            DrawGammaSetMarker(histoPhiInRData[iR], 20, markerSize, colorData, colorData);
            histoPhiInRData[iR]->Draw("same,hist");
            DrawGammaSetMarker(histoPhiInRMC[iR], 20, markerSize, colorMC, colorMC);
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

            DrawGammaLines(0.,2*TMath::Pi(),1., 1.,1.1,kGray+1);

            DrawGammaSetMarker(histoDataMCRatioPhiInR[iR], 20, markerSize, colorData, colorData);
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

            DrawGammaLines(-10,-10, histoZInRData[iR]->GetMinimum()*0.1,histoZInRData[iR]->GetMaximum()*1.5,1.1,kGray+1,2);
            DrawGammaLines(10,10, histoZInRData[iR]->GetMinimum()*0.1,histoZInRData[iR]->GetMaximum()*1.5,1.1,kGray+1,2);

            DrawGammaSetMarker(histoZInRData[iR], 20, markerSize, colorData, colorData);
            histoZInRData[iR]->Draw("same,hist");
            DrawGammaSetMarker(histoZInRMC[iR], 20, markerSize, colorMC, colorMC);
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

            DrawGammaLines(-10,-10, 0.3,1.55,1.1,kGray+1,2);
            DrawGammaLines(10,10,0.3,1.55,1.1,kGray+1,2);
            DrawGammaLines(-20,20,1., 1.,1.1,kGray+1);

            DrawGammaSetMarker(histoDataMCRatioZInR[iR], 20, markerSize, colorData, colorData);
            histoDataMCRatioZInR[iR]->Draw("same,hist");

        histoDummyZInRbinsDown->Draw("same,axis");
        canvasZInRbins->Update();
        canvasZInRbins->SaveAs(Form("%s/ZInRbins%s_%i_%s.%s",outputDirectory.Data(),optionPeriod.Data(),iR,fCutSelectionRead.Data(),suffix.Data()));

    }
    delete canvasZInRbins;


    // Ploting single Canvas
    cout << "Plotting..." << endl;
    Double_t EtaRange[2] = {-1.5,1.5};
    Double_t PhiRange[2] = {0.,6.28};
    Double_t PtRange[2]  = {0.,15.};
    Double_t RRange[2]   = {0.,180.};

    Double_t minYRatio = 0.9;
    Double_t maxYRatio = 1.5;

    Color_t colorOnFly = kRed+1;
    Color_t colorOffline = kBlue+1;

    TCanvas* canvasV0Finder = new TCanvas("canvasV0Finder","",1300,1000);
    DrawGammaCanvasSettings( canvasV0Finder,  0.13, 0.02, 0.02, 0.09);
    TH2F * histo2DDummy = new TH2F("","",1000,0.,200.,1000,0.,2.);

    // draw efficiency vs Pt, smaller eta range
    canvasV0Finder->SetLogy(1);
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#it{p}_{T} (GeV/#it{c})","Efficiency",0.04,0.04, 0.04,0.04, 1.,1.);
    histo2DDummy->GetXaxis()->SetRangeUser(PtRange[0],PtRange[1]);
    histo2DDummy->GetYaxis()->SetRangeUser(1e-2,5.);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoEffiPt, 20, 1.5,colorOnFly,colorOnFly);
        histoEffiPt->Draw("same,l");

        TLegend* legendV0Finder = GetAndSetLegend2(0.2,0.14,0.5,0.14+0.04*3.5,36);
        legendV0Finder->SetHeader(fCollisionSystem.Data());
        legendV0Finder->AddEntry(histoEffiPt,"","pl");
//         legendV0Finder->Draw();


    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/PhotonEffiPtSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    // draw efficiency vs R
    canvasV0Finder->SetLogy(0);
    SetStyleHistoTH2ForGraphs(histo2DDummy,"R (cm)","Efficiency",0.04,0.04, 0.04,0.04, 1.,1.);
    histo2DDummy->GetXaxis()->SetRangeUser(RRange[0],RRange[1]+5.);
    histo2DDummy->GetYaxis()->SetRangeUser(0.,.55);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoEffiR, 20, 1.5,colorOnFly,colorOnFly);
        histoEffiR->Draw("same,l");

        TLegend* legendV0Finder_R = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legendV0Finder_R->SetHeader(fCollisionSystem.Data());
        legendV0Finder_R->AddEntry(histoEffiR,"","pl");
//         legendV0Finder_R->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/PhotonEffiRSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    // draw efficiency vs phi
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#phi (rad)","Efficiency",0.04,0.04, 0.04,0.04, 1.,1.3);
    histo2DDummy->GetXaxis()->SetRangeUser(PhiRange[0],PhiRange[1]);
    histo2DDummy->GetYaxis()->SetRangeUser(0.,.4);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoEffiPhi, 20, 1.5,colorOnFly,colorOnFly);
        histoEffiPhi->Draw("same,l");

        TLegend* legendV0Finder_Phi = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legendV0Finder_Phi->SetHeader(fCollisionSystem.Data());
        legendV0Finder_Phi->AddEntry(histoEffiPhi,"","pl");
//         legendV0Finder_Phi->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/PhotonEffiEtaSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    // draw efficiency vs eta
    histo2DDummy = new TH2F("","",1000,-2.,2.,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#eta","Efficiency",0.04,0.04, 0.04,0.04, 1.,1.3);
    histo2DDummy->GetXaxis()->SetRangeUser(EtaRange[0],EtaRange[1]);
    histo2DDummy->GetYaxis()->SetRangeUser(0.,.4);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoEffiEta, 20, 1.5,colorOnFly,colorOnFly);
        histoEffiEta->Draw("same,l");

        TLegend* legendV0Finder_Eta = GetAndSetLegend2(0.17,0.93-0.04*2.5,0.5,0.93,36);
        legendV0Finder_Eta->SetHeader(fCollisionSystem.Data());
        legendV0Finder_Eta->AddEntry(histoEffiEta,"","pl");
//         legendV0Finder_Eta->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/PhotonEffiEtaSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));
    delete canvasV0Finder;


    TCanvas * canvasSinglePtPurity = new TCanvas("canvasSinglePtPurity","",1200,1000);
    DrawGammaCanvasSettings( canvasSinglePtPurity,  0.12, 0.02, 0.02, 0.12);
    DrawGammaSetMarker(histoPurityPt5cm, 20, markerSize, colorMC, colorMC);
    TH2F *histoDummyPuritySinglePt =  new TH2F("histoDummyPuritySinglePt","histoDummyPuritySinglePt",1000,0.,10.,1000,0.5,1.1);
    SetStyleHistoTH2ForGraphs(histoDummyPuritySinglePt, "p_{T} (GeV/c)","Purity",0.05,0.05, 0.05,0.05);
    // 0.9*textsizeLabelsDownpurity, textsizeLabelsDownpurity,0.9*textsizeLabelsDownpurity,textsizeLabelsDownpurity, 0.9,0.1/(textsizeFacDownpurity*marginXRatio));
    histoDummyPuritySinglePt->DrawCopy();
    histoPurityPt5cm->DrawCopy("same");
    canvasSinglePtPurity->Print(Form("%s/PhotonPuritySingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasSingleEtaDataMC = new TCanvas("canvasSingleEtaDataMC","",1200,1000);
    DrawGammaCanvasSettings( canvasSingleEtaDataMC,  0.15, 0.02, 0.07, 0.12);
    DrawGammaSetMarker(histoEtaData, 20, markerSize, colorData, colorData);
    DrawGammaSetMarker(histoEtaMC, 20, markerSize, colorMC, colorMC);
    TH2F *histoDummyEtaDataMC =  new TH2F("histoDummyEtaDataMC","histoDummyEtaDataMC",1000,-1.1,1.1,1000,0.0,1.2*histoEtaData->GetMaximum());
    SetStyleHistoTH2ForGraphs(histoDummyEtaDataMC, "#eta","#frac{dN_{#gamma}}{dEta}",0.05,0.05, 0.05,0.05);
    histoDummyEtaDataMC->DrawCopy();
    histoEtaData->Draw("same,histo");
    histoEtaMC->Draw("same,histo");
    TLegend* legenEta = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
    legenEta->AddEntry(histoEtaData,"Data","p");
    legenEta->AddEntry(histoEtaMC,"MC","p");
    legenEta->Draw();
    canvasSingleEtaDataMC->Print(Form("%s/PhotonEtaDataMCSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasSingleAsymDataMCLowP = new TCanvas("canvasSingleAsymDataMCLowP","",1200,1000);
    DrawGammaCanvasSettings( canvasSingleAsymDataMCLowP,  0.15, 0.02, 0.07, 0.12);
    DrawGammaSetMarker(histoAsymPDataLowP, 20, markerSize, colorData, colorData);
    DrawGammaSetMarker(histoAsymPTrueMCLowP, 20, markerSize, colorMC, colorMC);
    TH2F *histoDummyAsymDataMCLowP =  new TH2F("histoDummyAsymDataMCLowP","histoDummyAsymDataMCLowP",1000,0.,1.0,1000,0.0,1.2*histoAsymPDataLowP->GetMaximum());
    SetStyleHistoTH2ForGraphs(histoDummyAsymDataMCLowP, "#frac{p_{e}}{p_{#gamma}}","#frac{dN_{#gamma}}{dAsym}",0.05,0.05, 0.05,0.05);
    histoDummyAsymDataMCLowP->DrawCopy();
    histoAsymPDataLowP->Draw("same");
    histoAsymPTrueMCLowP->Draw("same");
    //   histoAsymPTrueMCLowP->Draw("same,hist");
    TLegend* legenAsymLowP = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
    legenAsymLowP->AddEntry(histoAsymPDataLowP,"Data","p");
    legenAsymLowP->AddEntry(histoAsymPTrueMCLowP,"MC","p");
    legenAsymLowP->Draw();
    canvasSingleAsymDataMCLowP->Print(Form("%s/PhotonAsymDataMCLowPSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasSingleAsymDataMCHighP = new TCanvas("canvasSingleAsymDataMCHighP","",1200,1000);
    DrawGammaCanvasSettings( canvasSingleAsymDataMCHighP,  0.15, 0.02, 0.07, 0.12);
    DrawGammaSetMarker(histoAsymPDataHighP, 20, markerSize, colorData, colorData);
    DrawGammaSetMarker(histoAsymPTrueMCHighP, 20, markerSize, colorMC, colorMC);
    TH2F *histoDummyAsymDataMCHighP =  new TH2F("histoDummyAsymDataMCHighP","histoDummyAsymDataMCHighP",1000,0.,1.0,1000,0.0,1.2*histoAsymPDataHighP->GetMaximum());
    SetStyleHistoTH2ForGraphs(histoDummyAsymDataMCHighP, "#frac{p_{e}}{p_{#gamma}}","#frac{dN_{#gamma}}{dAsym}",0.05,0.05, 0.05,0.05);
      histoDummyAsymDataMCHighP->DrawCopy();
    histoAsymPDataHighP->Draw("same");
    histoAsymPTrueMCHighP->Draw("same");
    histoAsymPTrueMCHighP->Draw("same,hist");
    TLegend* legenAsymHighP = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
    legenAsymHighP->AddEntry(histoAsymPDataHighP,"Data","p");
    legenAsymHighP->AddEntry(histoAsymPTrueMCHighP,"MC","p");
    legenAsymHighP->Draw();
    canvasSingleAsymDataMCHighP->Print(Form("%s/PhotonAsymDataMCHighPSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    TProfile* fProfileContainingMaterialBudgetWeightsFullPt = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBinsFullPt","profileContainingMaterialBudgetWeights_manyRadialBinsFullPt",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsFullPt->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsFullPt->GetYaxis()->SetTitle("weight_per_gamma");

    for(Int_t i=0; i<nBinsR; i++){
      cout<< arrayRBins[i] << " " << histoDataMCRatioRScaledToGas->GetBinContent(i+1) <<"  "<< histoDataMCRatioRScaledToGas->GetBinError(i+1)<< endl;
      fProfileContainingMaterialBudgetWeightsFullPt->Fill(arrayRBins[i], histoDataMCRatioRScaledToGas->GetBinContent(i+1));
    }
    cout<< endl;

    TProfile* fProfileContainingMaterialBudgetWeights01 = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins01","profileContainingMaterialBudgetWeights_manyRadialBins01",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeights01->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeights01->GetYaxis()->SetTitle("weight_per_gamma, p_{T} > 0.15 GeV/c");

    for(Int_t i=0; i<nBinsR; i++){
        cout<< arrayRBins[i] << " " << histoDataMCRatioRinPtBinScaledToGas01->GetBinContent(i+1) << "  " << histoDataMCRatioRinPtBinScaledToGas01->GetBinError(i+1) << endl;
        fProfileContainingMaterialBudgetWeights01->Fill(arrayRBins[i], histoDataMCRatioRinPtBinScaledToGas01->GetBinContent(i+1));
    }

    cout<< endl;
    TProfile* fProfileContainingMaterialBudgetWeights02 = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins02","profileContainingMaterialBudgetWeights_manyRadialBins02",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeights02->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeights02->GetYaxis()->SetTitle("weight_per_gamma, p_{T} > 0.3 GeV/c");

    for(Int_t i=0; i<nBinsR; i++){
        cout<< arrayRBins[i] << " " << histoDataMCRatioRinPtBinScaledToGas02->GetBinContent(i+1) << "  " << histoDataMCRatioRinPtBinScaledToGas02->GetBinError(i+1) << endl;
        fProfileContainingMaterialBudgetWeights02->Fill(arrayRBins[i], histoDataMCRatioRinPtBinScaledToGas02->GetBinContent(i+1));
    }
    cout<< endl;

    TProfile* fProfileContainingMaterialBudgetWeights03 = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins03","profileContainingMaterialBudgetWeights_manyRadialBins03",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeights03->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeights03->GetYaxis()->SetTitle("weight_per_gamma, p_{T} > 0.4 GeV/c");

    for(Int_t i=0; i<nBinsR; i++){
        cout<< arrayRBins[i] << " " << histoDataMCRatioRinPtBinScaledToGas03->GetBinContent(i+1) << "  " << histoDataMCRatioRinPtBinScaledToGas03->GetBinError(i+1) << endl;
        fProfileContainingMaterialBudgetWeights03->Fill(arrayRBins[i], histoDataMCRatioRinPtBinScaledToGas03->GetBinContent(i+1));
    }
    cout<< endl;
    // pT>0.4 declared as default

    TProfile* fProfileContainingMaterialBudgetWeights = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins","profileContainingMaterialBudgetWeights_manyRadialBins",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeights->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeights->GetYaxis()->SetTitle("weight_per_gamma, p_{T} > 0.4 GeV/c");

    for(Int_t i=0; i<nBinsR; i++){
      cout<< arrayRBins[i] << " " << histoDataMCRatioRinPtBinScaledToGas03->GetBinContent(i+1) <<"  " << histoDataMCRatioRinPtBinScaledToGas03->GetBinError(i+1) << endl;
        fProfileContainingMaterialBudgetWeights->Fill(arrayRBins[i], histoDataMCRatioRinPtBinScaledToGas03->GetBinContent(i+1));

	//        cout<< arrayRBins[i] << " " << histoDataMCRatioR->GetBinContent(i+1) <<endl;
	//        fProfileContainingMaterialBudgetWeights->Fill(arrayRBins[i], histoDataMCRatioR->GetBinContent(i+1));
    }


    //Form("%s/PhotonChar_%s.%s",optionPeriod.Data(),fCutSelectionRead.Data() ))
    cout<< Form("MCInputFileMaterialBudgetWeights%s_%s.root",optionPeriod.Data(),fCutSelectionRead.Data())<< endl;
    TFile outFile(Form("MCInputFileMaterialBudgetWeights%s_%s.root",optionPeriod.Data(),fCutSelectionRead.Data()) ,"RECREATE");

    fProfileContainingMaterialBudgetWeightsFullPt->Write();
    fProfileContainingMaterialBudgetWeights->Write();
    fProfileContainingMaterialBudgetWeights01->Write();
    fProfileContainingMaterialBudgetWeights02->Write();
    fProfileContainingMaterialBudgetWeights03->Write();

    histoRData->Write("Data");
    histoRMC->Write("MC");
    histoRDataScaledToGas->Write("DataScaledToGas");
    histoRMCScaledToGas->Write("MCScaledToGas");

    histoDataMCRatioRinPtBinScaledToGas01->Write();
    histoDataMCRatioRinPtBinScaledToGas02->Write();
    histoDataMCRatioRinPtBinScaledToGas03->Write();
    histoPurityPt5cm->Write();
    histoEffiR->Write();
    histoEffiPt->Write();
    histoAsymPDataLowP->Write();
    histoAsymPDataHighP->Write();
    histoAsymPTrueMCLowP->Write();
    histoAsymPTrueMCHighP->Write();

    outFile.Close();


}
