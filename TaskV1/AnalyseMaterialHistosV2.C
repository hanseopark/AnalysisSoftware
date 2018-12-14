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
#include "TH3F.h"
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

// nBinsR: 12 Rbins to calculate the weights. Optimized based on the material borders
// nBinsPtFine: To store minimum value of pT to calculate the weights pTMin-infinity
// nBinsPtMin: as previous but in the middle of the pt bin.
// Cocktail to be used in wideBinning 0.1 GeV/Bin
// MaterialHistos need a rebin of 2 . Original: 0-20 GeV/c and 400 Bins

TH1F* ScaleByIntegralWithinLimits(const TH1F* hist, Double_t rmin, Double_t rmax, Double_t mcGasCorrectionFactor )
{
    TH1F * histToBeScaled = (TH1F*)hist->Clone();
    Double_t nconvInGas =  histToBeScaled->Integral(hist->GetXaxis()->FindBin(rmin),hist->GetXaxis()->FindBin(rmax));
    GammaScalingHistogramm(histToBeScaled,1./(mcGasCorrectionFactor*nconvInGas));
    cout << "NconvInGas = " << nconvInGas << "  sqrt(nconvInGas) = " << TMath::Sqrt(nconvInGas) << " sqrt(nconvInGas)/nconvInGas = "<<  TMath::Sqrt(nconvInGas)/nconvInGas<<endl;
    return histToBeScaled;
}

Double_t CalculateIntegralWithinLimits(const TH1F* hist, Double_t rmin, Double_t rmax, Double_t mcGasCorrectionFactor )
{
    TH1F * histToBeIntegrated = (TH1F*)hist->Clone();
    // Option "width" removed becase it does not give proper number of entries
    Double_t nconvInRange =  histToBeIntegrated->Integral(hist->GetXaxis()->FindBin(rmin),hist->GetXaxis()->FindBin(rmax));
    nconvInRange*=mcGasCorrectionFactor;
    cout << "nconvInRange = " << nconvInRange << endl;
    //cout << "nconvInRange = " << nconvInRange << " for bins "  << hist->GetXaxis()->FindBin(rmin) << " - " << hist->GetXaxis()->FindBin(rmax)<< endl;
    return nconvInRange;
}

Double_t CalculateIntegralFromPtMin(const TH1F* hist, Double_t ptMin)
{
    TH1F * histToBeIntegrated = (TH1F*)hist->Clone();
    Double_t nconvInRangeFromPtMin =  histToBeIntegrated->Integral(hist->GetXaxis()->FindBin(ptMin),hist->GetXaxis()->GetNbins());
    //    cout << "Inside CalculateIntegralFromPtMin nconvInRangeFromPtMin = " << ptMin<< " " << hist->GetXaxis()->FindBin(ptMin) << " " << nconvInRangeFromPtMin << endl;
    //cout << "nconvInRange = " << nconvInRange << " for bins "  << hist->GetXaxis()->FindBin(rmin) << " - " << hist->GetXaxis()->FindBin(rmax)<< endl;
    return nconvInRangeFromPtMin;
}



void AnalyseMaterialHistosV2( TString fileName         = "",
                            TString fileNameMC       = "",
                            TString cutSelection     = "",
                            TString optionEnergy     = "",
                            TString optionPeriod     = "",
                            TString suffix           = "pdf",
                            TString outputFolderName = "MaterialBudgetHistos",
                            TString createDummy      = "",
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
    cout << "Cut selection: " << fCutSelectionRead.Data()<< endl;

    TString outputDirectory = Form("%s/%s/%s/%s",outputFolderName.Data(),fEnergyFlag.Data(),optionPeriod.Data(),fCutSelectionRead.Data());
    if(optionPeriod.CompareTo("")==0) outputDirectory= Form("%s/%s/%s",outputFolderName.Data(),fEnergyFlag.Data(),fCutSelectionRead.Data());
    gSystem->Exec("mkdir -p "+outputDirectory);
    cout<< "Output directory is: "<<  outputDirectory.Data()<<endl;

    // factor needed for Run1 simulations
    Double_t mcGasCorrectionFactor = 1.;
    if (optionPeriod.CompareTo("LHC10b") == 0 || optionPeriod.CompareTo("LHC10bc") == 0 || optionEnergy.Contains("PbPb") || createDummy.Contains("GasCorr")){
        mcGasCorrectionFactor = 0.960035693454;
    }

    // Cocktail file
    TString cutSelectionWithMesonCut;
    cutSelectionWithMesonCut=cutSelection;
    cutSelectionWithMesonCut.Append("_0152103500000000");
    cout<< " cutSelectionWithMesonCut::"<< cutSelectionWithMesonCut.Data()<< endl;

    if (!LoadSecondariesFromCocktailFile(cutSelectionWithMesonCut.Data(), fEnergyFlag.Data())){
      cout<< "Problems with the cocktail"<< endl;
    } else{
      cout<< "cocktail successfully loaded"<< endl;
    } 

    //The histograms  fHistoSecGammaCocktailFromXPt[k]  contain the secondaries from the 4 sources;

    TString periodData, periodMC, clusterName, generatorName, V0ReaderName;
    if(optionPeriod.Contains("Pythia")) generatorName = "Pythia";
    else if(optionPeriod.Contains("Phojet")) generatorName = "Phojet";
    else generatorName = "MC";

    if(outputFolderName.Contains("Onfly") || outputFolderName.Contains("onfly") ) V0ReaderName = "On-the-Fly V0";
    else if(outputFolderName.Contains("Offline") || outputFolderName.Contains("offline") ) generatorName = "Offline V0";

    // pt dependent binning for pT histograms 69 bins . 0-2. (40 bins a 0.05), 2-3 (10 bins a0.1) 3-10 (14 bins a 0.5) 10-20 (5 a 2GeV)
 

    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 40) arrayPtBins[i]              = 0.05*i;
      else if(i<50) arrayPtBins[i]           = 2.+ 0.1*(i-40);
      else if(i<64) arrayPtBins[i]          = 3. + 0.5*(i-50);
      else if(i<69) arrayPtBins[i]          = 10.+ 2.0*(i-64);
      else arrayPtBins[i]                    = fMaxPt;
      //      cout<< i<< " "<<  arrayPtBins[i]<< endl;
    }

    for(Int_t i=0; i<nBinsPtTwo+1;i++){
      if (i < 20) arrayPtBinsTwo[i]            = 0.1*i;
      else if(i<28) arrayPtBinsTwo[i]          = 2. + 0.5*(i-20);
      else if(i<30) arrayPtBinsTwo[i]          = 6. + 2.0*(i-28);
      else if(i<31) arrayPtBinsTwo[i]          = 10. + 10*(i-30);
      else arrayPtBinsTwo[i]                   = fMaxPt;
      //        cout<< i<< " "<<  arrayPtBinsTwo[i]<< endl;
    }


    //_______________________ Data _______________________
    TFile fData(fileName.Data());
    TList *TopDir                 = (TList*)fData.Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    TList *HistosGammaConversion  = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversion == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
        return;
    }

    TList *ESDContainer           = (TList*)HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    TH1F *histoEventQualityData   = (TH1F*)ESDContainer->FindObject("NEvents");
    TH1F *histoGoodESDTracksData  = (TH1F*)ESDContainer->FindObject("GoodESDTracksEta08");
    TH2F *histoRPhiData           = (TH2F*)ESDContainer->FindObject("ESD_Conversion_RPhi");
    TH2F *histoREtaData           = (TH2F*)ESDContainer->FindObject("ESD_Conversion_REta");
    TH2F *histoRZData             = (TH2F*)ESDContainer->FindObject("ESD_Conversion_RZ");
    TH2F *histoRPtData            = (TH2F*)ESDContainer->FindObject("ESD_Conversion_RPt");
    TH2F *histoRElecdEdxData      = (TH2F*)ESDContainer->FindObject("Electron_RdEdx");
    TH2F *histoRElecNSdEdxData    = (TH2F*)ESDContainer->FindObject("Electron_RNSigmadEdx");
    TH2F *histoAsymPData          = (TH2F*)ESDContainer->FindObject("ESD_ConversionMapping_AsymP");
    TH1F *histoDCAData            = (TH1F*)ESDContainer->FindObject("ESD_Conversion_DCA");
    TH1F *histoPsiPairData        = (TH1F*)ESDContainer->FindObject("ESD_Conversion_PsiPair");
    TH1F *histoChi2Data           = (TH1F*)ESDContainer->FindObject("ESD_Conversion_Chi2");
    TH1F *histoMassData           = (TH1F*)ESDContainer->FindObject("ESD_Conversion_Mass");
    TH1F *histoRRejLargeData      = (TH1F*)ESDContainer->FindObject("ESD_Conversion_RLarge");
    TH1F *histoRRejSmallData      = (TH1F*)ESDContainer->FindObject("ESD_Conversion_RSmall");

    //________________________ MC ________________________
    TFile fMC(fileNameMC.Data());
    TList *TopDirMC               = (TList*)fMC.Get(nameMainDir.Data());
    if(TopDirMC == NULL){
        cout<<"ERROR: TopDirMC not Found"<<endl;
        return;
    }
    TList *HistosGammaConversionMC= (TList*)TopDirMC->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversionMC == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in FileMC"<<endl;
        return;
    }

    TList *ESDContainerMC         = (TList*)HistosGammaConversionMC->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    TH1F *histoEventQualityMC     = (TH1F*)ESDContainerMC->FindObject("NEvents");
    TH1F *histoGoodESDTracksMC    = (TH1F*)ESDContainerMC->FindObject("GoodESDTracksEta08");

    TH1F *histoGoodESDTracksWeightedMC = NULL;
    histoGoodESDTracksWeightedMC  = (TH1F*)ESDContainerMC->FindObject("GoodESDTracksWeightedEta08");
    TH2F *histoRPhiMC             = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_RPhi");
    TH2F *histoREtaMC             = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_REta");
    TH2F *histoRZMC               = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_RZ");
    TH2F *histoRPtMC              = (TH2F*)ESDContainerMC->FindObject("ESD_Conversion_RPt");
    TH2F *histoRElecdEdxMC        = (TH2F*)ESDContainerMC->FindObject("Electron_RdEdx");
    TH2F *histoRElecNSdEdxMC      = (TH2F*)ESDContainerMC->FindObject("Electron_RNSigmadEdx");
    TH2F *histoAsymPMC            = (TH2F*)ESDContainerMC->FindObject("ESD_ConversionMapping_AsymP");
    TH1F *histoDCAMC              = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_DCA");
    TH1F *histoPsiPairMC          = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_PsiPair");
    TH1F *histoChi2MC             = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_Chi2");
    TH1F *histoMassMC             = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_Mass");
    TH1F *histoRRejLargeMC        = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_RLarge");
    TH1F *histoRRejSmallMC        = (TH1F*)ESDContainerMC->FindObject("ESD_Conversion_RSmall");

    TList *MCContainer            = (TList*)HistosGammaConversionMC->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
    TH1F * histoMCAllGammaPt      = (TH1F*)MCContainer->FindObject("MC_AllGamma_Pt");  // All Primary photons
 
    TH2F * histoConvRPhiMC        = (TH2F*)MCContainer->FindObject("MC_Conversion_RPhi");  // All Converted
    TH2F * histoConvREtaMC        = (TH2F*)MCContainer->FindObject("MC_Conversion_REta");  // All Converted
    TH2F * histoConvRPtMC         = (TH2F*)MCContainer->FindObject("MC_Conversion_RPt");   // All Converted

    TH2F *histoSecGammaPtMC       = (TH2F*) MCContainer->FindObject("MC_AllSecondaryGamma_Pt");
    TH3F *histoSecConvGammaPtRMC  = (TH3F*) MCContainer->FindObject("MC_SecondaryConvGamma_PtR");

    TList *TrueMCContainer        = (TList*)HistosGammaConversionMC->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
    TH2F *histoRPhiTrueMC         = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RPhi");
    TH2F *histoREtaTrueMC         = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_REta");
    TH2F *histoRZTrueMC           = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RZ");
    TH2F *histoRPtTrueMC          = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RPt");
    TH2F *histoRPtTruePrimMC      = (TH2F*)TrueMCContainer->FindObject("ESD_TruePrimConversion_RPt");
    TH2F *histoRPtTrueSecMC       = (TH2F*)TrueMCContainer->FindObject("ESD_TrueSecConversion_RPt");

    TH3F *histoTrueSecConvGammaRPtMC  = (TH3F*)TrueMCContainer->FindObject("ESD_TrueSecondaryConvGamma_Pt");
    TH3F *histoTrueSecConvGammaMCRPtMC = (TH3F*)TrueMCContainer->FindObject("ESD_TrueSecondaryConvGamma_MCPt");

    TH2F *histoRPtMCRPtTrueMC     = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversion_RPtMCRPt");
    TH2F *histoAsymPTrueMC        = (TH2F*)TrueMCContainer->FindObject("ESD_TrueConversionMapping_AsymP");
    TH1F *histoDCATrueMC          = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_DCA");
    TH1F *histoPsiPairTrueMC      = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_PsiPair");
    TH1F *histoChi2TrueMC         = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_Chi2");
    TH1F *histoMassTrueMC         = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_Mass");
    TH1F *histoRRejLargeTrueMC    = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_RLarge");
    TH1F *histoRRejSmallTrueMC    = (TH1F*)TrueMCContainer->FindObject("ESD_TrueConversion_RSmall");

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

    Float_t numberGoodEventsMC    = histoEventQualityMC->GetBinContent(1);
    Double_t meanMultMC;
    Float_t normFactorReconstMC;
    if( histoGoodESDTracksWeightedMC != NULL) {
        meanMultMC = histoGoodESDTracksWeightedMC->GetMean();
    } else {
        meanMultMC = histoGoodESDTracksMC->GetMean();
    }
    normFactorReconstMC = 1./(numberGoodEventsMC*meanMultMC);

    cout << "*********************************************" << endl;
    cout << "Number of good events in data -> " << numberGoodEventsData << endl;
    cout << "Mean multiplicity in data -> " << meanMultData << endl;
    cout << "Normalization factor in data -> " << normFactorReconstData << endl;
    cout << "-----" << endl;
    cout << "Number of good events in MC -> " << numberGoodEventsMC << endl;
    cout << "Mean multiplicity in MC -> " << meanMultMC << endl;
    cout << "Normalization factor in MC -> " << normFactorReconstMC << endl;
    cout << "-----" << endl;
    cout << "Ratio multiplicity data/MC = " << meanMultData/meanMultMC << endl;
    cout << "*********************************************" << endl;

    histoGoodESDTracksData->Scale(1./numberGoodEventsData);
    histoGoodESDTracksMC->Scale(1./numberGoodEventsMC);
    if( histoGoodESDTracksWeightedMC != 0x0) {
        histoGoodESDTracksWeightedMC->Scale(1./numberGoodEventsMC);
    }

    cout << __LINE__ << endl;

    //-AM   Don't do scaling to some histograms prior to calculation of errors
 
    GammaScalingHistogramm(histoRElecdEdxData, normFactorReconstData);
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
    cout << "\nProjection data histograms..." << endl;
    TH1F *histoRData      = (TH1F*)histoRPtData->ProjectionY("histoRData");
    TH1F *histoRDataRebin = (TH1F*)histoRData->Clone("histoRDataRebin");
    ConvGammaRebinWithBinCorrection(histoRDataRebin,rebinRPlots);

    for(Int_t i=0; i<3; i++){
        histoPtinRBinData[i]   = (TH1F*)histoRPtData->ProjectionX(Form("histoPtinRBinData_%i",i),
                                                                  histoRPtData->GetYaxis()->FindBin(projRBins[i]+eps),
                                                                  histoRPtData->GetYaxis()->FindBin(projRBins[i+1]-eps),"e");
        cout << "Projecting pT in " << projRBins[i] << " < R < " << projRBins[i+1] << " cm" << endl;
    }
 
    TH1F *histoPtData          = (TH1F*)histoRPtData->ProjectionX("histoPtData");

    for(Int_t i=0; i<6; i++){
      histoRinPtBinData[i]   = (TH1F*)histoRPtData->ProjectionY(Form("histoRinPtBinData_%i",i),
                                                                  histoRPtData->GetXaxis()->FindBin(projPtBins[i]+eps),
                                                                  histoRPtData->GetNbinsX(),"e");
                                                                //histoRPtData->GetXaxis()->FindBin(projPtBins[i+1]-eps),"e");

      // cout << "Projecting R from " << projPtBins[i] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c "<< endl;
      // cout << "Mean value of the projected distribution: " << histoRinPtBinData[i]->GetMean() << endl;
    }

    TH1F* histoPtDataRebin = (TH1F*) histoPtData->Rebin(nBinsPt,"histoPtDataRebin",arrayPtBins);
    DivideTH1ByBinWidth(histoPtDataRebin);


    //projecting in finer pT bins
    for(Int_t i=0; i < nBinsPtFine; i++){
        histoRinPtBinDataFine[i] = (TH1F*)histoRPtData->ProjectionY(Form("histoRinPtBinDataFine_%i",i),
                                                                    histoRPtData->GetXaxis()->FindBin(projPtBinsFine[i]+0.001),
                                                                    histoRPtData->GetNbinsX(),"e");
        cout << "Projecting R from " << projPtBinsFine[i] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c "<< endl;
    }

    TH1F *histoEtaData        = (TH1F*)histoREtaData->ProjectionX("histoEtaData");
    TH1F *histoAsymPDataLowP  = (TH1F*)histoAsymPData->ProjectionY("histoAsymPDataLowP", 1, histoAsymPData->GetXaxis()->FindBin(0.4-eps),"e");
    ConvGammaRebinWithBinCorrection(histoAsymPDataLowP,4);
    TH1F *histoAsymPDataHighP = (TH1F*)histoAsymPData->ProjectionY("histoAsymPDataHighP", histoAsymPData->GetXaxis()->FindBin(0.4+eps), histoAsymPData->GetNbinsX(),"e");
    ConvGammaRebinWithBinCorrection(histoAsymPDataHighP,4);


    cout << "Projection MC reconstructed histograms..." << endl;
    TH1F *histoRMC            = (TH1F*)histoRPtMC->ProjectionY("histoRMC");
    for(Int_t i=0; i<3; i++){
        histoPtinRBinMC[i]    = (TH1F*)histoRPtMC->ProjectionX(Form("histoPtinRBinMC_%i",i),
                                                                histoRPtMC->GetYaxis()->FindBin(projRBins[i]+eps),
                                                                histoRPtMC->GetYaxis()->FindBin(projRBins[i+1]-eps),"e");
        //cout << "Projecting pT in " << projRBins[i] << " < R < " << projRBins[i+1] << " cm" << endl;
    }
    TH1F *histoRMCRebin = (TH1F*)histoRMC->Clone("histoRMCRebin");
    ConvGammaRebinWithBinCorrection(histoRMCRebin,rebinRPlots);

    TH1F* histoPt5cmMC = (TH1F*)histoPtinRBinMC[1]->Clone("histoPt5cmMC");
    histoPt5cmMC->Sumw2();
    histoPt5cmMC->Add(histoPtinRBinMC[2]);


    TH1F *histoPtMC          = (TH1F*)histoRPtMC->ProjectionX("histoPtMC");
    for(Int_t i=0; i<6; i++){
        histoRinPtBinMC[i]   = (TH1F*)histoRPtMC->ProjectionY(Form("histoRinPtBinMC_%i",i),
                                                              histoRPtMC->GetXaxis()->FindBin(projPtBins[i]+eps),
                                                              histoRPtMC->GetNbinsX(),"e");
        //cout << "Projecting R from " << projPtBins[i] << " GeV/c to " << histoRPtMC->GetXaxis()->GetBinUpEdge(histoRPtMC->GetNbinsX()) << " GeV/c "<< endl;
    }

    TH1F* histoPtMCRebin = (TH1F*) histoPtMC->Rebin(nBinsPt,"histoPtMCRebin",arrayPtBins);
    DivideTH1ByBinWidth(histoPtMCRebin);


    //Projecting in finer pT bins
    for(Int_t i=0; i<nBinsPtFine; i++){
        histoRinPtBinMCFine[i] = (TH1F*)histoRPtMC->ProjectionY(Form("histoRinPtBinMCFine_%i",i),
                                                                histoRPtMC->GetXaxis()->FindBin(projPtBinsFine[i]+0.001),
                                                                histoRPtMC->GetNbinsX(),"e");
        cout << "Projecting R from " << projPtBinsFine[i] << " GeV/c to " << histoRPtMC->GetXaxis()->GetBinUpEdge(histoRPtMC->GetNbinsX()) << " GeV/c "<< endl;
    }

    TH1F *histoEtaMC        = (TH1F*)histoREtaMC->ProjectionX("histoEtaMC");
    TH1F *histoAsymPMCLowP  = (TH1F*)histoAsymPMC->ProjectionY("histoAsymPMCLowP", 1, histoAsymPMC->GetXaxis()->FindBin(0.4-eps),"e");
    ConvGammaRebinWithBinCorrection(histoAsymPMCLowP,4);
    TH1F *histoAsymPMCHighP = (TH1F*)histoAsymPMC->ProjectionY("histoAsymPMCHighP", histoAsymPMC->GetXaxis()->FindBin(0.4+eps), histoAsymPMC->GetNbinsX(),"e");
    ConvGammaRebinWithBinCorrection(histoAsymPMCHighP,4);


    cout << "Projection MC generated histograms..." << endl;
    TH1F *histoConvRMC      = (TH1F*)histoConvRPtMC->ProjectionY("histoConvRMC");
    TH1F *histoConvPtMC     = (TH1F*)histoConvRPtMC->ProjectionX("histoConvPtMC");
    TH1F *histoConvPhiMC    = (TH1F*)histoConvRPhiMC->ProjectionX("histoConvPhiMC");
    TH1F *histoConvEtaMC    = (TH1F*)histoConvREtaMC->ProjectionX("histoConvEtaMC");

    // all secondary Photons within the eta acceptance
    TH1F* histoSecGammaPtFromK0S = (TH1F*)histoSecGammaPtMC->ProjectionX("histoSecGammaPtFromK0S",1,1,"e");
    TH1F* histoSecGammaPtAllSources = (TH1F*)histoSecGammaPtMC->ProjectionX("histoSecGammaPtAllSources",1,4,"e");
    histoSecGammaPtAllSources->Rebin(rebinPtPlots);

    // all secondary Photons that converted within the eta acceptance
    TH1F * histoSecConvGammaPtFromK0S = (TH1F*)histoSecConvGammaPtRMC->ProjectionX("histoSecConvGammaPtFromK0S",1,histoSecConvGammaPtRMC->GetNbinsY(),1,1,"e");
    TH1F * histoSecConvGammaPtAllSources = (TH1F*)histoSecConvGammaPtRMC->ProjectionX("histoSecConvGammaPtAllSources",1,histoSecConvGammaPtRMC->GetNbinsY(),1,4,"e");

 
    //AM: I need a rebin of 2 in pT to match the pT from the cocktail
    Bool_t  hasCocktailInput                                       = kTRUE;

    cout << "Projection MC true histograms..." << endl;
    TH1F *histoRTrueMC      = (TH1F*)histoRPtTrueMC->ProjectionY("histoRTrueMC");
    TH1F *histoRTruePrimMC  = (TH1F*)histoRPtTruePrimMC->ProjectionY("histoRTruePrimMC");
    TH1F *histoRTrueSecMC   = (TH1F*)histoRPtTrueSecMC->ProjectionY("histoRTrueSecMC");

    TH1F * histoPtTrueSecFromK0S    = (TH1F*) histoTrueSecConvGammaRPtMC->ProjectionX("histoPtTrueSecFromK0S",0,-1,1,1,"e");
    TH1F * histoPtTrueSecRest       = (TH1F*) histoTrueSecConvGammaRPtMC->ProjectionX("histoPtTrueSecRest",0,-1,4,4,"e");
    TH1F * histoPtTrueSecAllsources = (TH1F*) histoTrueSecConvGammaRPtMC->ProjectionX("histoPtTrueSecAllSources",0,-1,1,4,"e");

    for(Int_t j=0; j < nBinsPtMin+1; j++){
        if (j==(nBinsPtMin)) arrayBinsPtMin[j] = projPtBinsFine[j-1]+0.05;
        else                 arrayBinsPtMin[j] = projPtBinsFine[j]-0.05;
    }

    cout << __LINE__ << endl;

    for(Int_t k=0; k < 4; k++){
      histoMCSecGammaPtBySource[k]     = (TH1F*) histoSecGammaPtMC->ProjectionX(Form("histoSecGammaPtFrom%s",fSecondaries[k].Data()),k+1,k+1,"e");
      histoMCSecConvGammaPtBySource[k] = (TH1F*) histoSecConvGammaPtRMC->ProjectionX(Form("histoSecConvGammaPtFrom%s",fSecondaries[k].Data()),1,histoSecConvGammaPtRMC->GetNbinsY(),k+1,k+1,"e");
      histoMCSecGammaPtBySource[k]->Rebin(rebinPtPlots);
      histoMCSecConvGammaPtBySource[k]->Rebin(rebinPtPlots);
      if( k == 0 ){
	nBins = histoMCSecGammaPtBySource[k]->GetNbinsX();
	xMin  = histoMCSecGammaPtBySource[k]->GetXaxis()->GetXmin();
	xMax  = histoMCSecGammaPtBySource[k]->GetXaxis()->GetXmax();
      }
      histoConvProbSecBySource[k] = (TH1F*) histoMCSecConvGammaPtBySource[k]->Clone(Form("histoSecConvProbGammaPtFrom%s",fSecondaries[k].Data()));
      histoConvProbSecBySource[k]->Sumw2();
      histoConvProbSecBySource[k] ->Divide(histoMCSecConvGammaPtBySource[k],histoMCSecGammaPtBySource[k] ,1,1,"B");
    }

 
    for(Int_t j=0; j<nBinsR; j++){
      //      cout<< "testing secondary weighting with MBW:: "<<fMaterialWeightsForSecEffCor[j]<< endl;
      // Data reconstructed
      histoPtEachRBinData[j] = (TH1F*)histoRPtData->ProjectionX(Form("histoPtEachRBinData_%i",j),
								histoRPtData->GetYaxis()->FindBin(arrayRBins[j]+eps),
								histoRPtData->GetYaxis()->FindBin(arrayRBins[j+1]-eps),"e");
      histoPtEachRBinData[j]->Rebin(rebinPtPlots);

      histoPtEachRBinMC[j] = (TH1F*)histoRPtMC->ProjectionX(Form("histoPtEachRBinMC_%i",j),
							    histoRPtMC->GetYaxis()->FindBin(arrayRBins[j]+eps),
							    histoRPtMC->GetYaxis()->FindBin(arrayRBins[j+1]-eps),"e");
      histoPtEachRBinMC[j]->Rebin(rebinPtPlots);
 
      // cout<<"Entries in each r bin MC-rec::"<<j<<endl;
      // CalculateIntegralFromPtMin( histoPtEachRBinMC[j],0.05);
      // MC truth
      histoPtTrueMCEachRBin[j] = (TH1F*)histoRPtTrueMC->ProjectionX(Form("histoPtTrueMCEachRBin_%i",j),
								    histoRPtTrueMC->GetYaxis()->FindBin(arrayRBins[j]+eps),
								    histoRPtTrueMC->GetYaxis()->FindBin(arrayRBins[j+1]-eps),"e");

      histoPtTrueMCEachRBin[j]->Rebin(rebinPtPlots); 
      // MC truth Primary
      histoPtTruePrimMCEachRBin[j] = (TH1F*)histoRPtTruePrimMC->ProjectionX(Form("histoPtTruePrimMCEachRBin_%i",j),
									    histoRPtTruePrimMC->GetYaxis()->FindBin(arrayRBins[j]+eps),
									    histoRPtTruePrimMC->GetYaxis()->FindBin(arrayRBins[j+1]-eps),"e");
      histoPtTruePrimMCEachRBin[j]->Rebin(rebinPtPlots); 


      // MC truth Secondary from K0S
      histoPtTrueSecEachRBinFromK0S[j] = 
	(TH1F*) histoTrueSecConvGammaRPtMC->ProjectionX(Form("histoPtTrueSecEachRBinFromK0S_%i",j), 
							      histoTrueSecConvGammaRPtMC->GetYaxis()->FindBin(arrayRBins[j]+eps),
							      histoTrueSecConvGammaRPtMC->GetYaxis()->FindBin(arrayRBins[j+1]-eps),1,1,"e");
      histoPtTrueSecEachRBinFromK0S[j]->Rebin(rebinPtPlots); 
      // MC truth Secondary from All Sources
      histoPtTrueSecEachRBinAllSources[j] = 
	(TH1F*) histoTrueSecConvGammaRPtMC->ProjectionX(Form("histoPtTrueSecEachRBinAllSources_%i",j), 
							      histoTrueSecConvGammaRPtMC->GetYaxis()->FindBin(arrayRBins[j]+eps),
							      histoTrueSecConvGammaRPtMC->GetYaxis()->FindBin(arrayRBins[j+1]-eps),1,4,"e");

      histoPtTrueSecEachRBinAllSources[j]->Rebin(rebinPtPlots); 

      histoPtEachRBinDataSecSubtractedUsingCocktail[j]  =  (TH1F*)histoPtEachRBinData[j]->Clone(Form("histoPtEachRBinDataSecSubtractedUsingCocktail_%i",j));
      histoPtEachRBinDataSecSubtractedUsingCocktail[j]->Sumw2();


      for(Int_t k=0; k < 4; k++){
	histoMCTrueRecSecGammaPtEachRBinBySource[k][j] = (TH1F*) histoTrueSecConvGammaRPtMC->ProjectionX(Form("histoMCTrueRecSecGammaPtEachRBinBySource%s_%i",fSecondaries[k].Data(),j), 
								 histoTrueSecConvGammaRPtMC->GetYaxis()->FindBin(arrayRBins[j]+eps),
								 histoTrueSecConvGammaRPtMC->GetYaxis()->FindBin(arrayRBins[j+1]-eps),k+1,k+1,"e");

	histoMCTrueRecSecGammaPtEachRBinBySource[k][j]->Rebin(rebinPtPlots);

	// AM: I need to divide by converted photons in EachRBin ; but it is not yet available . This efficiency is only reconstruction. Conversion probability is taken out.
	histoTrueRecEffSecEachRBinBySource[k][j] = (TH1F*)histoMCTrueRecSecGammaPtEachRBinBySource[k][j]->Clone(Form("histoTrueRecEffSecPtEachRBinBySource%s_%i",fSecondaries[k].Data(),j));
	histoTrueRecEffSecEachRBinBySource[k][j] ->Sumw2();
	histoTrueRecEffSecEachRBinBySource[k][j] ->Divide( histoMCTrueRecSecGammaPtEachRBinBySource[k][j], histoMCSecConvGammaPtBySource[k] ,1,1,"B");

	// AM: This efficiency includes recostruction and conversion probability .
	histoTrueEffSecEachRBinBySource[k][j] = (TH1F*)histoMCTrueRecSecGammaPtEachRBinBySource[k][j]->Clone(Form("histoTrueEffSecPtEachRBinBySource%s_%i",fSecondaries[k].Data(),j));
	histoTrueEffSecEachRBinBySource[k][j] ->Sumw2();
	histoTrueEffSecEachRBinBySource[k][j] ->Divide( histoMCTrueRecSecGammaPtEachRBinBySource[k][j], histoMCSecGammaPtBySource[k] ,1,1,"B");

	// AM: Secondary frac from MC for each source   !!!!!!! Check!!!!!!!!!!!!
	histoFracSecPtEachRBinBySource[k][j] = (TH1F*)histoMCTrueRecSecGammaPtEachRBinBySource[k][j] ->Clone(Form("histoFracSecPtEachRBinBySource%s_%i",fSecondaries[k].Data(),j));
	histoFracSecPtEachRBinBySource[k][j]->Sumw2();
	histoFracSecPtEachRBinBySource[k][j]->Divide(histoMCTrueRecSecGammaPtEachRBinBySource[k][j],histoPtEachRBinMC[j],1,1,"B"); 


	// To get the expected secondaries from cocktail: cocktail*TrueEff*Nevt_data  or cocktail*TrueRecEff*convProb *Nevt_data
	// fHistoSecGammaCocktailFromXPt[k]
	// ConvertCocktailSecondaryToRaw from CorrectGammaV2
	// ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[0], histoGammaSecondaryFromXConvProb_MCPt[0],
	//		                         histoGammaSecFromXRecoEff_RecPt[0], histoGammaTrueSecondaryFromX_MCPt_recPt[0], nEvt,
	//		                         kFALSE);
	
	if(k<3)fHistoSecGammaCocktailFromXPt[k]->SetBins(nBins,xMin,xMax);
	if( k<3 ){
	  histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[k][j] = (TH1F*)fHistoSecGammaCocktailFromXPt[k]->Clone(Form("histoDataFromCocktailSecConvGammaPtEachRBinBySource%s_Raw%i",fSecondaries[k].Data(),j));
	  histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[k][j]->Sumw2();
	  hasCocktailInput = ConvertCocktailSecondaryToRaw(histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[k][j], 
							   histoConvProbSecBySource[k],histoTrueEffSecEachRBinBySource[k][j],histoRPtMCRPtTrueMC,numberGoodEventsData,kFALSE);
	  //							 histoConvProbSecBySource[k],histoTrueRecEffSecEachRBinBySource[k][j],histoRPtMCRPtTrueMC,numberGoodEventsData,kFALSE);
	  //AM: I need to think if numberGoodEventsData is ok, or I need some scaling ; to be compared to how Nevt for V0And is calculated for pp

	}else{
	  // For Source Rest the sec contribution is the same using cocktail or secondary fraction
	  histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[k][j] =  (TH1F*)histoPtEachRBinData[j]->Clone(Form("histoDataFromCocktailSecConvGammaPtEachRBinBySource%s_Raw%i",fSecondaries[k].Data(),j));
	  histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[k][j] ->Sumw2();
	  histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[k][j]->Multiply(histoFracSecPtEachRBinBySource[k][j]);	   
	}

	histoPtEachRBinDataSecYieldFromSecFracBySourceRaw[k][j]  =  (TH1F*)histoPtEachRBinData[j]->Clone(Form("histoPtEachRBinDataSecYieldFromSecFracBySource%s_Raw%i",fSecondaries[k].Data(),j));
	histoPtEachRBinDataSecYieldFromSecFracBySourceRaw[k][j] ->Sumw2();
	histoPtEachRBinDataSecYieldFromSecFracBySourceRaw[k][j] ->Multiply(histoFracSecPtEachRBinBySource[k][j]);
	
	histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[k][j]->Scale(fMaterialWeightsForSecEffCor[j]);

	//AM: subtracting each of the sources one by one
	histoPtEachRBinDataSecSubtractedUsingCocktail[j]->Add(histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[k][j],-1.);

      }
           

     // Secondary Efficiency: each R Bin compared to all secondary photons (all RBins together)

      histoEfficiencySecEachRBinAllSources[j] = (TH1F*)histoPtTrueSecEachRBinAllSources[j]->Clone(Form("histoEfficiencySecPtEachRBinAllSources_%i",j));
      histoEfficiencySecEachRBinAllSources[j] ->Sumw2();
      histoEfficiencySecEachRBinAllSources[j] ->Divide(histoPtTrueSecEachRBinAllSources[j],histoSecGammaPtAllSources,1,1,"B");


      // Purity in each R Bin

      histoPurityPtEachRBin[j] = (TH1F*)histoPtTrueMCEachRBin[j] ->Clone(Form("histoPurityPtEachRBin_%i",j));
      histoPurityPtEachRBin[j]->Sumw2();
      histoPurityPtEachRBin[j]->Divide(histoPtTrueMCEachRBin[j],histoPtEachRBinMC[j],1,1,"B"); 

      histoPurityPrimPtEachRBin[j] = (TH1F*)histoPtTruePrimMCEachRBin[j] ->Clone(Form("histoPurityPrimPtEachRBin_%i",j));
      histoPurityPrimPtEachRBin[j]->Sumw2();
      histoPurityPrimPtEachRBin[j]->Divide(histoPtTruePrimMCEachRBin[j],histoPtEachRBinMC[j],1,1,"B"); 

      // Fraction from true secondary photons
      histoFracSecPtEachRBin[j] = (TH1F*)histoPtTrueSecEachRBinAllSources[j] ->Clone(Form("histoFracSecPtEachRBinAllSources_%i",j));
      histoFracSecPtEachRBin[j]->Sumw2();
      histoFracSecPtEachRBin[j]->Divide(histoPtTrueSecEachRBinAllSources[j],histoPtEachRBinMC[j],1,1,"B"); 
 
      histoPtEachRBinMCSecYieldFromSecFrac[j] =  (TH1F*)histoPtEachRBinMC[j]->Clone(Form("histoPtEachRBinMCSecYieldFromSecFrac_%i",j));
      histoPtEachRBinMCSecYieldFromSecFrac[j]->Sumw2();
      histoPtEachRBinMCSecYieldFromSecFrac[j]->Multiply(histoFracSecPtEachRBin[j]);

      histoPtEachRBinDataSecYieldFromSecFrac[j] =  (TH1F*)histoPtEachRBinData[j]->Clone(Form("histoPtEachRBinDataSecYieldFromSecFrac_%i",j));
      histoPtEachRBinDataSecYieldFromSecFrac[j]->Sumw2();
      histoPtEachRBinDataSecYieldFromSecFrac[j]->Multiply(histoFracSecPtEachRBin[j]);

      histoPtEachRBinDataSecSubtracted[j]  =  (TH1F*)histoPtEachRBinData[j]->Clone(Form("histoPtEachRBinDataSecSubtracted_%i",j));
      histoPtEachRBinDataSecSubtracted[j]->Sumw2();
      histoPtEachRBinDataSecSubtracted[j]->Add(histoPtEachRBinDataSecYieldFromSecFrac[j],-1.);

      histoPtEachRBinMCSecSubtracted[j]  =  (TH1F*)histoPtEachRBinMC[j]->Clone(Form("histoPtEachRBinMCSecSubtracted_%i",j));
      histoPtEachRBinMCSecSubtracted[j]->Sumw2();
      histoPtEachRBinMCSecSubtracted[j]->Add(histoPtEachRBinMCSecYieldFromSecFrac[j],-1.);

      for(Int_t k=0; k < nBinsPtFine; k++){
	nConvInRangeFromPtMinSecSubtractedData[j][k] =  CalculateIntegralFromPtMin( histoPtEachRBinDataSecSubtracted[j],projPtBinsFine[k]);
	nConvInRangeFromPtMinSecSubtractedDataRelErr[j][k] = TMath::Power(nConvInRangeFromPtMinSecSubtractedData[j][k],0.5)/nConvInRangeFromPtMinSecSubtractedData[j][k] ;
	nConvInRangeFromPtMinSecSubtractedDataUsingCocktail[j][k] =  CalculateIntegralFromPtMin( histoPtEachRBinDataSecSubtractedUsingCocktail[j],projPtBinsFine[k]);
	nConvInRangeFromPtMinSecSubtractedDataUsingCocktailRelErr[j][k] = TMath::Power(nConvInRangeFromPtMinSecSubtractedDataUsingCocktail[j][k],0.5)/nConvInRangeFromPtMinSecSubtractedDataUsingCocktail[j][k] ;
	nConvInRangeFromPtMinSecSubtractedMC[j][k] = CalculateIntegralFromPtMin( histoPtEachRBinMCSecSubtracted[j],projPtBinsFine[k]);
	nConvInRangeFromPtMinSecSubtractedMCRelErr[j][k] = TMath::Power(nConvInRangeFromPtMinSecSubtractedMC[j][k],0.5)/nConvInRangeFromPtMinSecSubtractedMC[j][k];
	//	cout << " Relative errors data and MC,j,k::"<< j<< " "<< k<< " " << nConvInRangeFromPtMinSecSubtractedDataRelErr[j][k] << "  " << nConvInRangeFromPtMinSecSubtractedMCRelErr[j][k]<< endl;
      }

    }

    // Create a pT spectrum after secondary subtraction
    TH1F* histoPtSumRMCSecSubtracted = (TH1F*)histoPtEachRBinMCSecSubtracted[0]->Clone("histoPtSumRMCSecSubtracted");
    histoPtSumRMCSecSubtracted->Sumw2();
    TH1F* histoPtSumRBinDataSecSubtractedUsingCocktail = (TH1F*)histoPtEachRBinDataSecSubtractedUsingCocktail[0]->Clone("histoPtSumRBinDataSecSubtractedUsingCocktail");
    histoPtSumRBinDataSecSubtractedUsingCocktail->Sumw2();
    for(Int_t j=1; j < nBinsR; j++){
      histoPtSumRMCSecSubtracted->Add(histoPtEachRBinMCSecSubtracted[j]);
      histoPtSumRBinDataSecSubtractedUsingCocktail->Add(histoPtEachRBinDataSecSubtractedUsingCocktail[j]);
    }

    TH1F* histoPtSumRMCSecSubtractedRebin = (TH1F*) histoPtSumRMCSecSubtracted->Rebin(nBinsPtTwo,"histoPtSumRMCSecSubtractedRebin",arrayPtBinsTwo);
    DivideTH1ByBinWidth(histoPtSumRMCSecSubtractedRebin);

    TH1F* histoPtSumRBinDataSecSubtractedUsingCocktailRebin = (TH1F*)histoPtSumRBinDataSecSubtractedUsingCocktail->Rebin(nBinsPtTwo,"histoPtSumRBinDataSecSubtractedUsingCocktailRebin",arrayPtBinsTwo);
    DivideTH1ByBinWidth(histoPtSumRBinDataSecSubtractedUsingCocktailRebin);



    for(Int_t j=0; j < nBinsR; j++){
      for(Int_t k=0; k < nBinsPtFine; k++){
	nConvInRangeFromPtMinSecSubtractedDataToGas[j][k] =  nConvInRangeFromPtMinSecSubtractedData[j][k]/nConvInRangeFromPtMinSecSubtractedData[10][k];
	nConvInRangeFromPtMinSecSubtractedDataUsingCocktailToGas[j][k] =  nConvInRangeFromPtMinSecSubtractedDataUsingCocktail[j][k]/nConvInRangeFromPtMinSecSubtractedDataUsingCocktail[10][k];


	//  mcGasCorrectionFactor applied to the gas here when integrating pT histograms. it is set to 1 except for LHC10b,c and PbPb at 2.76
	nConvInRangeFromPtMinSecSubtractedMCToGas[j][k] = nConvInRangeFromPtMinSecSubtractedMC[j][k]/(mcGasCorrectionFactor*nConvInRangeFromPtMinSecSubtractedMC[10][k]);
	weightInRangeFromPtMinSecSubtracted[j][k] = nConvInRangeFromPtMinSecSubtractedDataToGas[j][k] /	nConvInRangeFromPtMinSecSubtractedMCToGas[j][k];
	weightInRangeFromPtMinSecSubtractedUsingCocktail[j][k] = nConvInRangeFromPtMinSecSubtractedDataUsingCocktailToGas[j][k] /nConvInRangeFromPtMinSecSubtractedMCToGas[j][k];


	// AM: For the gas we count the stat error only once. Otherwise error sqrt(2) too large
	if(j!=10){

	  nConvInRangeFromPtMinSecSubtractedDataToGasRelErr[j][k] =TMath::Power( (TMath::Power(nConvInRangeFromPtMinSecSubtractedDataRelErr[j][k],2)+
										TMath::Power(nConvInRangeFromPtMinSecSubtractedDataRelErr[10][k],2) ),0.5);
	  nConvInRangeFromPtMinSecSubtractedDataUsingCocktailToGasRelErr[j][k] =TMath::Power( (TMath::Power(nConvInRangeFromPtMinSecSubtractedDataUsingCocktailRelErr[j][k],2)+
										TMath::Power(nConvInRangeFromPtMinSecSubtractedDataUsingCocktailRelErr[10][k],2) ),0.5);

	  nConvInRangeFromPtMinSecSubtractedMCToGasRelErr[j][k] =TMath::Power( (TMath::Power(nConvInRangeFromPtMinSecSubtractedMCRelErr[j][k],2)+
										TMath::Power(nConvInRangeFromPtMinSecSubtractedMCRelErr[10][k],2) ),0.5);
	}else{

	  nConvInRangeFromPtMinSecSubtractedDataToGasRelErr[j][k] =TMath::Power( (TMath::Power(nConvInRangeFromPtMinSecSubtractedDataRelErr[j][k],2) ),0.5);
	  nConvInRangeFromPtMinSecSubtractedDataUsingCocktailToGasRelErr[j][k] =TMath::Power( (TMath::Power(nConvInRangeFromPtMinSecSubtractedDataUsingCocktailRelErr[j][k],2)),0.5);
	  nConvInRangeFromPtMinSecSubtractedMCToGasRelErr[j][k] =TMath::Power( (TMath::Power(nConvInRangeFromPtMinSecSubtractedMCRelErr[j][k],2)),0.5);
	}

	  weightInRangeFromPtMinSecSubtractedRelErr[j][k] =  TMath::Power( (TMath::Power(nConvInRangeFromPtMinSecSubtractedDataToGasRelErr[j][k],2)+
									  TMath::Power(nConvInRangeFromPtMinSecSubtractedMCToGasRelErr[j][k],2)),0.5);

	  weightInRangeFromPtMinSecSubtractedUsingCocktailRelErr[j][k] =  TMath::Power( (TMath::Power(nConvInRangeFromPtMinSecSubtractedDataUsingCocktailToGasRelErr[j][k],2)+
									  TMath::Power(nConvInRangeFromPtMinSecSubtractedMCToGasRelErr[j][k],2)),0.5);

      }
    }

    // for(Int_t j=0; j < nBinsR; j++){
    //   for(Int_t k=0; k < nBinsPtFine; k++){
    // 	cout<<"SUM ptBin:: "<< k<< " Rbin:: "<<j<<" " <<  weightInRangeFromPtMinSecSubtracted[j][k]<< " " << weightInRangeFromPtMinSecSubtractedRelErr[j][k]  << endl;
    // 	cout<<"SUC ptBin:: "<< k<< " Rbin:: "<<j<<" " <<  weightInRangeFromPtMinSecSubtractedUsingCocktail[j][k]<< " " << weightInRangeFromPtMinSecSubtractedUsingCocktailRelErr[j][k]  << endl;
    //   }
    //   cout<< endl;
    // }

    for(Int_t j=0; j < nBinsR; j++){
      histoWeightsEachRPtMinSecSub[j] = new TH1F(Form("histoWeightsEachRPtMinSecSub%i",j), Form("histoWeightsEachRPtMinSecSub%i",j), nBinsPtMin,  arrayBinsPtMin);
      histoWeightsEachRPtMinSecSubUsingCocktail[j] = new TH1F(Form("histoWeightsEachRPtMinSecSubUsingCocktail%i",j), Form("histoWeightsEachRPtMinSecSubUsingCocktail%i",j), nBinsPtMin,  arrayBinsPtMin);
      for(Int_t k=0; k < nBinsPtMin; k++){
	histoWeightsEachRPtMinSecSub[j]->SetBinContent(k+1,weightInRangeFromPtMinSecSubtracted[j][k]);
	histoWeightsEachRPtMinSecSub[j]->SetBinError(k+1,weightInRangeFromPtMinSecSubtractedRelErr[j][k]*weightInRangeFromPtMinSecSubtracted[j][k]);
	histoWeightsEachRPtMinSecSubUsingCocktail[j]->SetBinContent(k+1,weightInRangeFromPtMinSecSubtractedUsingCocktail[j][k]);
	histoWeightsEachRPtMinSecSubUsingCocktail[j]->SetBinError(k+1,weightInRangeFromPtMinSecSubtractedUsingCocktailRelErr[j][k]*weightInRangeFromPtMinSecSubtractedUsingCocktail[j][k]);

      }
    }


    for(Int_t i=0; i<3; i++){
        histoPtinRBinTrueMC[i]  = (TH1F*)histoRPtTrueMC->ProjectionX(Form("histoPtinRBinTrueMC_%i",i),
                                                                     histoRPtTrueMC->GetYaxis()->FindBin(projRBins[i]+eps),
                                                                     histoRPtTrueMC->GetYaxis()->FindBin(projRBins[i+1]-eps),"e");

        histoPtinRBinTruePrimMC[i]  = (TH1F*)histoRPtTruePrimMC->ProjectionX(Form("histoPtinRBinTruePrimMC_%i",i),
                                                                     histoRPtTruePrimMC->GetYaxis()->FindBin(projRBins[i]+eps),
                                                                     histoRPtTruePrimMC->GetYaxis()->FindBin(projRBins[i+1]-eps),"e");

        histoPtinRBinTrueSecMC[i]  = (TH1F*)histoRPtTrueSecMC->ProjectionX(Form("histoPtinRBinTrueSecMC_%i",i),
                                                                     histoRPtTrueSecMC->GetYaxis()->FindBin(projRBins[i]+eps),
                                                                     histoRPtTrueSecMC->GetYaxis()->FindBin(projRBins[i+1]-eps),"e");
    }
    TH1F *histoRMCRTrueMC       = (TH1F*)histoRPtMCRPtTrueMC->ProjectionY("histoRMCRTrueMC");

    // Summing (5,35) cm + (35,180)cm
    TH1F* histoPt5cmTrueMC = (TH1F*)histoPtinRBinTrueMC[1]->Clone("histoPt5cmTrueMC");
    histoPt5cmTrueMC->Sumw2();
    histoPt5cmTrueMC->Add(histoPtinRBinTrueMC[2]);

    TH1F* histoPt5cmTruePrimMC = (TH1F*)histoPtinRBinTruePrimMC[1]->Clone("histoPt5cmTruePrimMC");
    histoPt5cmTruePrimMC->Sumw2();
    histoPt5cmTruePrimMC->Add(histoPtinRBinTruePrimMC[2]);

    TH1F* histoPt5cmTrueSecMC = (TH1F*)histoPtinRBinTrueSecMC[1]->Clone("histoPt5cmTrueSecMC");
    histoPt5cmTrueSecMC->Sumw2();
    histoPt5cmTrueSecMC->Add(histoPtinRBinTrueSecMC[2]);


    TH1F *histoPi0DalRTrueMC    = (TH1F*)histoPi0DalRPtTrueMC->ProjectionY("histoPi0DalRTrueMC");
    TH1F *histoEtaDalRTrueMC    = (TH1F*)histoEtaDalRPtTrueMC->ProjectionY("histoEtaDalRTrueMC");
    TH1F *histoCombRTrueMC      = (TH1F*)histoCombRPtTrueMC->ProjectionY("histoCombRTrueMC");
    TH1F *histoPtTrueMC         = (TH1F*)histoRPtTrueMC->ProjectionX("histoPtTrueMC");

    for(Int_t i=0; i<6; i++){
        histoRinPtBinTrueMC[i]   = (TH1F*)histoRPtTrueMC->ProjectionY(Form("histoRinPtBinTrueMC_%i",i),
                                                                      histoRPtTrueMC->GetXaxis()->FindBin(projPtBins[i]+eps),
                                                                      histoRPtTrueMC->GetNbinsX(),"e");


       histoRinPtBinTruePrimMC[i] = (TH1F*)histoRPtTruePrimMC->ProjectionY(Form("histoRinPtBinTruePrimMC_%i",i),
                                                                      histoRPtTruePrimMC->GetXaxis()->FindBin(projPtBins[i]+eps),
                                                                      histoRPtTruePrimMC->GetNbinsX(),"e");


       histoRinPtBinTrueSecMC[i]  = (TH1F*)histoRPtTrueSecMC->ProjectionY(Form("histoRinPtBinTrueSecMC_%i",i),
                                                                      histoRPtTrueSecMC->GetXaxis()->FindBin(projPtBins[i]+eps),
                                                                      histoRPtTrueSecMC->GetNbinsX(),"e");
    }
    TH1F *histoPtMCPtTrueMC   = (TH1F*)histoRPtMCRPtTrueMC->ProjectionX("histoPtMCPtTrueMC");

    // histograms are re-defined. To be cheked/understood! Should it be Pt AM. 21.03.18 ???
    TH1F *histoPi0DalPtTrueMC = (TH1F*)histoPi0DalRPtTrueMC->ProjectionX("histoPi0DalPtTrueMC");
    TH1F *histoEtaDalPtTrueMC = (TH1F*)histoEtaDalRPtTrueMC->ProjectionX("histoEtaDalPtTrueMC");
    TH1F *histoCombPtTrueMC   = (TH1F*)histoCombRPtTrueMC->ProjectionX("histoCombPtTrueMC");


    TH1F *histoPhiTrueMC        = (TH1F*)histoRPhiTrueMC->ProjectionX("histoPhiTrueMC");
    TH1F *histoEtaTrueMC        = (TH1F*)histoREtaTrueMC->ProjectionX("histoEtaTrueMC");
    TH1F *histoAsymPTrueMCLowP  = (TH1F*)histoAsymPTrueMC->ProjectionY("histoAsymPTrueMCLowP", 1, histoAsymPTrueMC->GetXaxis()->FindBin(0.4-eps),"e");
    ConvGammaRebinWithBinCorrection(histoAsymPTrueMCLowP,4);
    TH1F *histoAsymPTrueMCHighP = (TH1F*)histoAsymPTrueMC->ProjectionY("histoAsymPTrueMCHighP", histoAsymPTrueMC->GetXaxis()->FindBin(0.4+eps), histoAsymPTrueMC->GetNbinsX(),"e");
    ConvGammaRebinWithBinCorrection(histoAsymPTrueMCHighP,4);


    //____________ normalize with Nconversions in the gas __________________
    cout << __LINE__ << endl;
    cout << "*******************************\nData: " << endl;
    cout << "Nconv in gas for full pT range" << endl;
    Double_t nconvInRangeData          = CalculateIntegralWithinLimits(histoRData, rMinGas+eps, rMaxGas-eps,1.);
    Double_t dataStatErrorGas          = TMath::Sqrt(nconvInRangeData);
    Double_t dataStatRelErrorGas       = TMath::Sqrt(nconvInRangeData)/nconvInRangeData;
    histoIntegralGasDataFullRange      = new TH1D("histoIntegralGasDataFullRange", "histoIntegralGasDataFullRange", nBinsR, arrayRBins);

    cout << "Nconv in gas from " << projPtBins[0] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c " << endl;
    Double_t nconvInRangeDataPtBin1    = CalculateIntegralWithinLimits(histoRinPtBinData[0], rMinGas+eps, rMaxGas-eps,1.);
    Double_t dataStatErrorGasPtBin1    = TMath::Sqrt(nconvInRangeDataPtBin1);
    Double_t dataStatRelErrorGasPtBin1 = TMath::Sqrt(nconvInRangeDataPtBin1)/nconvInRangeDataPtBin1;
    histoIntegralGasDataPtBin1         = new TH1D("histoIntegralGasDataPtBin1", "histoIntegralGasDataPtBin1", nBinsR, arrayRBins);

    cout << "Nconv in gas from " << projPtBins[1] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c " << endl;
    Double_t nconvInRangeDataPtBin2    = CalculateIntegralWithinLimits(histoRinPtBinData[1], rMinGas+eps, rMaxGas-eps,1.);
    Double_t dataStatErrorGasPtBin2    = TMath::Sqrt(nconvInRangeDataPtBin2);
    Double_t dataStatRelErrorGasPtBin2 = TMath::Sqrt(nconvInRangeDataPtBin2)/nconvInRangeDataPtBin2;
    histoIntegralGasDataPtBin2         = new TH1D("histoIntegralGasDataPtBin2", "histoIntegralGasDataPtBin2", nBinsR, arrayRBins);

    cout << "Nconv in gas from " << projPtBins[2] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c " << endl;
    Double_t nconvInRangeDataPtBin3    = CalculateIntegralWithinLimits(histoRinPtBinData[2], rMinGas+eps, rMaxGas-eps,1.);
    Double_t dataStatErrorGasPtBin3    = TMath::Sqrt(nconvInRangeDataPtBin3);
    Double_t dataStatRelErrorGasPtBin3 = TMath::Sqrt(nconvInRangeDataPtBin3)/nconvInRangeDataPtBin3;
    histoIntegralGasDataPtBin3         = new TH1D("histoIntegralGasDataPtBin3", "histoIntegralGasDataPtBin3", nBinsR, arrayRBins);


    for(Int_t i=0; i < nBinsPtFine; i++){
        cout << "Nconv in gas from " << projPtBinsFine[i] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c " << endl;
        nconvInRangeDataFine[i] = CalculateIntegralWithinLimits(histoRinPtBinDataFine[i], rMinGas+eps, rMaxGas-eps,1.);
        dataStatErrorGasFine[i] = TMath::Sqrt(nconvInRangeDataFine[i]);
        dataStatRelErrorGasFine[i] = TMath::Sqrt(nconvInRangeDataFine[i])/nconvInRangeDataFine[i];
        histoIntegralGasDataFine[i] = new TH1D(Form("histoIntegralGasDataFine%d",i), Form("histoIntegralGasDataFine%d",i), nBinsR, arrayRBins);
        for(Int_t j=1; j<nBinsR+1; j++){
            histoIntegralGasDataFine[i]->SetBinContent(j,nconvInRangeDataFine[i]);
            if(j!=11){
	      histoIntegralGasDataFine[i]->SetBinError(j,dataStatErrorGasFine[i]);
	    }else{
	      histoIntegralGasDataFine[i]->SetBinError(j,0.);
	    } 
        }

        histoRinPtBinDataRebinFine[i]  = (TH1F*)histoRinPtBinDataFine[i]->Rebin(nBinsR,Form("histoRinPtBinDataRebinFine%d",i), arrayRBins);
        histoRinPtBinDataRebinFine[i]->Sumw2();
        histoRinPtBinDataScaledToGasRebinFine[i] = (TH1F*)histoRinPtBinDataRebinFine[i]->Clone(Form("histoRinPtBinDataScaledToGasRebinFine%d",i));
        histoRinPtBinDataScaledToGasRebinFine[i]->Divide(histoRinPtBinDataRebinFine[i],histoIntegralGasDataFine[i],1.,1.,"");
    }


    for(Int_t i=1; i<nBinsR+1; i++){
      histoIntegralGasDataFullRange->SetBinContent(i,nconvInRangeData);
      histoIntegralGasDataPtBin1->SetBinContent(i,nconvInRangeDataPtBin1);
      histoIntegralGasDataPtBin2->SetBinContent(i,nconvInRangeDataPtBin2);
      histoIntegralGasDataPtBin3->SetBinContent(i,nconvInRangeDataPtBin3);
      if(i!=11){
	histoIntegralGasDataFullRange->SetBinError(i,dataStatErrorGas);
	histoIntegralGasDataPtBin1->SetBinError(i,dataStatErrorGasPtBin1);
	histoIntegralGasDataPtBin2->SetBinError(i,dataStatErrorGasPtBin2);
	histoIntegralGasDataPtBin3->SetBinError(i,dataStatErrorGasPtBin3);
      }else{
	histoIntegralGasDataFullRange->SetBinError(i,0.);
	histoIntegralGasDataPtBin1->SetBinError(i,0.);
	histoIntegralGasDataPtBin2->SetBinError(i,0.);
	histoIntegralGasDataPtBin3->SetBinError(i,0.);
      } 
    }

    TH1F* histoIntegralGasDataWide = (TH1F*) histoRData->Clone("histoIntegralGasDataWide");
    for(Int_t i=1; i<histoIntegralGasDataWide->GetNbinsX()+1; i++){
      histoIntegralGasDataWide->SetBinContent(i,nconvInRangeData);
      if(i!=11){
	histoIntegralGasDataWide->SetBinError(i,dataStatErrorGas);
      }else{
	histoIntegralGasDataWide->SetBinError(i,0.);
      } 
    }

    histoRDataScaledToGas = (TH1F*) histoRData->Clone("histoRDataScaledToGas");
    histoRDataScaledToGas->Sumw2();
    histoRDataScaledToGas->Divide(histoRDataScaledToGas , histoIntegralGasDataWide,1.,1.,"");

    histoRDataRebin   = (TH1F*)histoRData->Rebin(nBinsR,"histoRDataRebin", arrayRBins);
    histoRDataRebin->Sumw2();
    histoRDataScaledToGasRebin =(TH1F*)histoRDataRebin->Clone("histoRDataScaledToGasRebin");
    histoRDataScaledToGasRebin->Divide(histoRDataRebin,histoIntegralGasDataFullRange,1.,1.,"");

    cout<< "testing Data::"<< histoRDataRebin->GetBinContent(11)<< " " <<histoRDataRebin->GetBinError(11)<<endl; 
    cout<< "testing Data-Gas::"<< histoIntegralGasDataFullRange->GetBinContent(11)<< " " << histoIntegralGasDataFullRange->GetBinError(11)<< "  "<<histoIntegralGasDataFullRange->GetBinError(11)/histoIntegralGasDataFullRange->GetBinContent(11) << endl;
    cout<< "Testing::"<< histoRDataScaledToGasRebin->GetBinContent(11)<< " "<< histoRDataScaledToGasRebin->GetBinError(11)<<endl;

    histoRinPtBinDataPtBin1Rebin = (TH1F*)histoRinPtBinData[0]->Rebin(nBinsR,"histoRinPtBinDataPtBin1Rebin", arrayRBins);
    histoRinPtBinDataPtBin1Rebin->Sumw2();
    histoRinPtBinDataScaledToGasPtBin1Rebin= (TH1F*)histoRinPtBinDataPtBin1Rebin->Clone("histoRinPtBinDataScaledToGasPtBin1Rebin");
    histoRinPtBinDataScaledToGasPtBin1Rebin->Divide(histoRinPtBinDataPtBin1Rebin,histoIntegralGasDataPtBin1,1.,1.,"");

    histoRinPtBinDataPtBin2Rebin = (TH1F*)histoRinPtBinData[1]->Rebin(nBinsR,"histoRinPtBinDataPtBin2Rebin", arrayRBins);
    histoRinPtBinDataPtBin2Rebin->Sumw2();
    histoRinPtBinDataScaledToGasPtBin2Rebin= (TH1F*)histoRinPtBinDataPtBin2Rebin->Clone("histoRinPtBinDataScaledToGasPtBin2Rebin");
    histoRinPtBinDataScaledToGasPtBin2Rebin->Divide(histoRinPtBinDataPtBin2Rebin,histoIntegralGasDataPtBin2,1.,1.,"");

    histoRinPtBinDataPtBin3Rebin = (TH1F*)histoRinPtBinData[2]->Rebin(nBinsR,"histoRinPtBinDataPtBin3Rebin", arrayRBins);
    histoRinPtBinDataPtBin3Rebin->Sumw2();
    histoRinPtBinDataScaledToGasPtBin3Rebin= (TH1F*)histoRinPtBinDataPtBin3Rebin->Clone("histoRinPtBinDataScaledToGasPtBin3Rebin");
    histoRinPtBinDataScaledToGasPtBin3Rebin->Divide(histoRinPtBinDataPtBin3Rebin,histoIntegralGasDataPtBin3,1.,1.,"");


    cout << __LINE__ << endl;
    cout << "*******************************\nMC: " << endl;
    cout << "Nconv in gas for full pT range" << endl;
    Double_t nconvInRangeMC          = CalculateIntegralWithinLimits(histoRMC, rMinGas+eps, rMaxGas-eps, mcGasCorrectionFactor);
    Double_t mcStatErrorGas          = TMath::Sqrt(nconvInRangeMC);
    Double_t mcStatRelErrorGas       = TMath::Sqrt(nconvInRangeMC)/nconvInRangeMC;
    histoIntegralGasMCFullRange      = new TH1D("histoIntegralGasMCFullRange", "histoIntegralGasMCFullRange", nBinsR, arrayRBins);

    cout << "Nconv in gas from " << projPtBins[0] << " GeV/c to " << histoRPtMC->GetXaxis()->GetBinUpEdge(histoRPtMC->GetNbinsX()) << " GeV/c " << endl;
    Double_t nconvInRangeMCPtBin1        = CalculateIntegralWithinLimits(histoRinPtBinMC[0], rMinGas+eps, rMaxGas-eps, mcGasCorrectionFactor);
    Double_t mcStatErrorGasPtBin1    = TMath::Sqrt(nconvInRangeMCPtBin1);
    Double_t mcStatRelErrorGasPtBin1 = TMath::Sqrt(nconvInRangeMCPtBin1)/nconvInRangeMCPtBin1;
    histoIntegralGasMCPtBin1         = new TH1D("histoIntegralGasMCPtBin1", "histoIntegralGasMCPtBin1", nBinsR, arrayRBins);

    cout << "Nconv in gas from " << projPtBins[1] << " GeV/c to " << histoRPtMC->GetXaxis()->GetBinUpEdge(histoRPtMC->GetNbinsX()) << " GeV/c " << endl;
    Double_t nconvInRangeMCPtBin2        = CalculateIntegralWithinLimits(histoRinPtBinMC[1], rMinGas+eps, rMaxGas-eps, mcGasCorrectionFactor);
    Double_t mcStatErrorGasPtBin2    = TMath::Sqrt(nconvInRangeMCPtBin2);
    Double_t mcStatRelErrorGasPtBin2 = TMath::Sqrt(nconvInRangeMCPtBin2)/nconvInRangeMCPtBin2;
    histoIntegralGasMCPtBin2         = new TH1D("histoIntegralGasMCPtBin2", "histoIntegralGasMCPtBin2", nBinsR, arrayRBins);

    cout << "Nconv in gas from " << projPtBins[2] << " GeV/c to " << histoRPtMC->GetXaxis()->GetBinUpEdge(histoRPtMC->GetNbinsX()) << " GeV/c " << endl;
    Double_t nconvInRangeMCPtBin3        = CalculateIntegralWithinLimits(histoRinPtBinMC[2], rMinGas+eps, rMaxGas-eps, mcGasCorrectionFactor);
    Double_t mcStatErrorGasPtBin3    = TMath::Sqrt(nconvInRangeMCPtBin3);
    Double_t mcStatRelErrorGasPtBin3 = TMath::Sqrt(nconvInRangeMCPtBin3)/nconvInRangeMCPtBin3;
    histoIntegralGasMCPtBin3         = new TH1D("histoIntegralGasMCPtBin3", "histoIntegralGasMCPtBin3", nBinsR, arrayRBins);

    for(Int_t i=0; i < nBinsPtFine; i++){
        cout << "Nconv in gas from " << projPtBinsFine[i] << " GeV/c to " << histoRPtMC->GetXaxis()->GetBinUpEdge(histoRPtMC->GetNbinsX()) << " GeV/c " << endl;
        nconvInRangeMCFine[i] = CalculateIntegralWithinLimits(histoRinPtBinMCFine[i], rMinGas+eps, rMaxGas-eps,1.);
        dataStatErrorGasMCFine[i] = TMath::Sqrt(nconvInRangeMCFine[i]);
        dataStatRelErrorGasMCFine[i] = TMath::Sqrt(nconvInRangeMCFine[i])/nconvInRangeMCFine[i];
        histoIntegralGasMCFine[i] = new TH1D(Form("histoIntegralGasMCFine%d",i), Form("histoIntegralGasMCFine%d",i), nBinsR, arrayRBins);
        for(Int_t j=1; j<nBinsR+1; j++){
	  histoIntegralGasMCFine[i]->SetBinContent(j,nconvInRangeMCFine[i]);
	  if(j!=11){
	    histoIntegralGasMCFine[i]->SetBinError(j,dataStatErrorGasMCFine[i]);
	  }else{
	    //	    cout<< "testing setting error to 0::"<< histoIntegralGasMCFine[i]->GetBinCenter(j) << endl;
	    histoIntegralGasMCFine[i]->SetBinError(j,0.);
	  } 
        }

        histoRinPtBinMCRebinFine[i]  = (TH1F*)histoRinPtBinMCFine[i]->Rebin(nBinsR,Form("histoRinPtBinMCRebinFine%d",i), arrayRBins);
        histoRinPtBinMCRebinFine[i]->Sumw2();
        histoRinPtBinMCScaledToGasRebinFine[i] = (TH1F*)histoRinPtBinMCRebinFine[i]->Clone(Form("histoRinPtBinMCScaledToGasRebinFine%d",i));
        histoRinPtBinMCScaledToGasRebinFine[i]->Divide(histoRinPtBinMCRebinFine[i],histoIntegralGasMCFine[i],1.0,1.0,"");
    }
    cout << "*******************************" << endl;

    for(Int_t i=1; i<nBinsR+1; i++){
        histoIntegralGasMCFullRange->SetBinContent(i,nconvInRangeMC);
        histoIntegralGasMCPtBin1->SetBinContent(i,nconvInRangeMCPtBin1);
        histoIntegralGasMCPtBin2->SetBinContent(i,nconvInRangeMCPtBin2);
        histoIntegralGasMCPtBin3->SetBinContent(i,nconvInRangeMCPtBin3);
	if(i!=11){
	  histoIntegralGasMCFullRange->SetBinError(i,mcStatErrorGas);
	  histoIntegralGasMCPtBin1->SetBinError(i,mcStatErrorGasPtBin1);
	  histoIntegralGasMCPtBin2->SetBinError(i,mcStatErrorGasPtBin2);
	  histoIntegralGasMCPtBin3->SetBinError(i,mcStatErrorGasPtBin3);
	}else{
	  histoIntegralGasMCFullRange->SetBinError(i,0.);
	  histoIntegralGasMCPtBin1->SetBinError(i,0.);
	  histoIntegralGasMCPtBin2->SetBinError(i,0.);
	  histoIntegralGasMCPtBin3->SetBinError(i,0.);
	}
    }

    // -AM The histogram is cloned to get the binning. Taking histoRData is ok
    TH1F* histoIntegralGasMCWide = (TH1F*) histoRData->Clone("histoIntegralGasMCWide");
    for(Int_t i=1; i<histoIntegralGasMCWide->GetNbinsX()+1; i++){
        histoIntegralGasMCWide->SetBinContent(i,nconvInRangeMC);
        histoIntegralGasMCWide->SetBinError(i,mcStatErrorGas);
    }

    histoRMCScaledToGas = (TH1F*) histoRMC->Clone("histoRMCScaledToGas");
    histoRMCScaledToGas->Sumw2();
    histoRMCScaledToGas->Divide(histoRMCScaledToGas , histoIntegralGasMCWide,1.,1.,"");

    histoRMCRebin   = (TH1F*)histoRMC->Rebin(nBinsR,"histoRMCRebin", arrayRBins);
    histoRMCRebin->Sumw2();
    histoRMCScaledToGasRebin =(TH1F*)histoRMCRebin->Clone("histoRMCScaledToGasRebin");
    histoRMCScaledToGasRebin->Divide(histoRMCRebin,histoIntegralGasMCFullRange,1.,1.,"");

    histoRinPtBinMCPtBin1Rebin = (TH1F*)histoRinPtBinMC[0]->Rebin(nBinsR,"histoRinPtBinMCPtBin1Rebin", arrayRBins);
    histoRinPtBinMCPtBin1Rebin->Sumw2();
    histoRinPtBinMCScaledToGasPtBin1Rebin= (TH1F*)histoRinPtBinMCPtBin1Rebin->Clone("histoRinPtBinMCScaledToGasPtBin1Rebin");
    histoRinPtBinMCScaledToGasPtBin1Rebin->Divide(histoRinPtBinMCPtBin1Rebin,histoIntegralGasMCPtBin1,1.,1.,"");

    histoRinPtBinMCPtBin2Rebin = (TH1F*)histoRinPtBinMC[1]->Rebin(nBinsR,"histoRinPtBinMCPtBin2Rebin", arrayRBins);
    histoRinPtBinMCPtBin2Rebin->Sumw2();
    histoRinPtBinMCScaledToGasPtBin2Rebin= (TH1F*)histoRinPtBinMCPtBin2Rebin->Clone("histoRinPtBinMCScaledToGasPtBin2Rebin");
    histoRinPtBinMCScaledToGasPtBin2Rebin->Divide(histoRinPtBinMCPtBin2Rebin,histoIntegralGasMCPtBin2,1.,1.,"");

    histoRinPtBinMCPtBin3Rebin = (TH1F*)histoRinPtBinMC[2]->Rebin(nBinsR,"histoRinPtBinMCPtBin3Rebin", arrayRBins);
    histoRinPtBinMCPtBin3Rebin->Sumw2();
    histoRinPtBinMCScaledToGasPtBin3Rebin= (TH1F*)histoRinPtBinMCPtBin3Rebin->Clone("histoRinPtBinMCScaledToGasPtBin3Rebin");
    histoRinPtBinMCScaledToGasPtBin3Rebin->Divide(histoRinPtBinMCPtBin3Rebin,histoIntegralGasMCPtBin3,1.,1.,"");


    //_______________ ratio data to MC after scaling __________________
    cout << __LINE__ << endl;
    histoDataMCRatioRScaledToGas = (TH1F*)histoRDataScaledToGasRebin->Clone("histoDataMCRatioRScaledToGas");
    histoDataMCRatioRScaledToGas->Divide(histoRDataScaledToGasRebin,histoRMCScaledToGasRebin,1.,1.,"");
    //   cout<< "checking stat errors Ratio::"<<  histoDataMCRatioRScaledToGas->GetBinCenter(11)<< " "<< histoDataMCRatioRScaledToGas->GetBinContent(11)<< " " << histoDataMCRatioRScaledToGas->GetBinError(11)<< endl;
    //    cout<< "checking stat errors Data::"<<  histoRDataScaledToGasRebin->GetBinContent(11)<< " " << histoRDataScaledToGasRebin->GetBinError(11)<< endl;
    //    cout<< "checking stat errors MC::"<<  histoRMCScaledToGasRebin->GetBinContent(11)<< " " << histoRMCScaledToGasRebin->GetBinError(11)<< endl;


    histoDataMCRatioRScaledToGasSecSub = (TH1F*)histoRDataScaledToGasRebin->Clone("histoDataMCRatioRScaledToGasSecSub");
    for(Int_t i=0;i<histoDataMCRatioRScaledToGasSecSub->GetNbinsX()+1;i++){
      //      cout<< "testing new histograms" << i<< "   "<< histoDataMCRatioRScaledToGasSecSub->GetBinCenter(i+1)<< endl;
      histoDataMCRatioRScaledToGasSecSub->SetBinContent(i+1,weightInRangeFromPtMinSecSubtractedUsingCocktail[i][0]);
      histoDataMCRatioRScaledToGasSecSub->SetBinError(i+1,weightInRangeFromPtMinSecSubtractedUsingCocktail[i][0]*weightInRangeFromPtMinSecSubtractedUsingCocktailRelErr[i][0]);
    }

    cout << "scaling to Nconv in gas for pT>0.1" << endl;
    histoDataMCRatioRinPtBinScaledToGasPtBin1  = (TH1F*)histoRinPtBinDataScaledToGasPtBin1Rebin->Clone("histoDataMCRatioRinPtBinScaledToGasPtBin1");
    histoDataMCRatioRinPtBinScaledToGasPtBin1->Divide(histoRinPtBinDataScaledToGasPtBin1Rebin ,histoRinPtBinMCScaledToGasPtBin1Rebin ,1.,1.,"");
    histoDataMCRatioRScaledToGasSecSubPtBin1 = (TH1F*)histoRinPtBinDataScaledToGasPtBin1Rebin->Clone("histoDataMCRatioRScaledToGasSecSubPtBin1");
    for(Int_t i=0;i<histoDataMCRatioRScaledToGasSecSubPtBin1->GetNbinsX()+1;i++){
      histoDataMCRatioRScaledToGasSecSubPtBin1->SetBinContent(i+1,weightInRangeFromPtMinSecSubtractedUsingCocktail[i][1]);
      histoDataMCRatioRScaledToGasSecSubPtBin1->SetBinError(i+1,weightInRangeFromPtMinSecSubtractedUsingCocktail[i][1]*weightInRangeFromPtMinSecSubtractedUsingCocktailRelErr[i][1]);
    }

    cout << "scaling to Nconv in gas for pT>0.3" << endl;
    histoDataMCRatioRinPtBinScaledToGasPtBin2  = (TH1F*)histoRinPtBinDataScaledToGasPtBin2Rebin->Clone("histoDataMCRatioRinPtBinScaledToGasPtBin2");
    histoDataMCRatioRinPtBinScaledToGasPtBin2->Divide(histoRinPtBinDataScaledToGasPtBin2Rebin ,histoRinPtBinMCScaledToGasPtBin2Rebin ,1.,1.,"");
    histoDataMCRatioRScaledToGasSecSubPtBin2 = (TH1F*)histoRinPtBinDataScaledToGasPtBin2Rebin->Clone("histoDataMCRatioRScaledToGasSecSubPtBin2");
    for(Int_t i=0;i<histoDataMCRatioRScaledToGasSecSubPtBin2->GetNbinsX()+1;i++){
      histoDataMCRatioRScaledToGasSecSubPtBin2->SetBinContent(i+1,weightInRangeFromPtMinSecSubtractedUsingCocktail[i][3]);
      histoDataMCRatioRScaledToGasSecSubPtBin2->SetBinError(i+1,weightInRangeFromPtMinSecSubtractedUsingCocktail[i][3]*weightInRangeFromPtMinSecSubtractedUsingCocktailRelErr[i][3]);
    }

    cout << "scaling to Nconv in gas for pT>0.4" << endl;
    histoDataMCRatioRinPtBinScaledToGasPtBin3  = (TH1F*)histoRinPtBinDataScaledToGasPtBin3Rebin->Clone("histoDataMCRatioRinPtBinScaledToGasPtBin3");
    histoDataMCRatioRinPtBinScaledToGasPtBin3->Divide(histoRinPtBinDataScaledToGasPtBin3Rebin ,histoRinPtBinMCScaledToGasPtBin3Rebin ,1.,1.,"");
    histoDataMCRatioRScaledToGasSecSubPtBin3 = (TH1F*)histoRinPtBinDataScaledToGasPtBin3Rebin->Clone("histoDataMCRatioRScaledToGasSecSubPtBin3");
    for(Int_t i=0;i<histoDataMCRatioRScaledToGasSecSubPtBin3->GetNbinsX()+1;i++){
      histoDataMCRatioRScaledToGasSecSubPtBin3->SetBinContent(i+1,weightInRangeFromPtMinSecSubtractedUsingCocktail[i][4]);
      histoDataMCRatioRScaledToGasSecSubPtBin3->SetBinError(i+1,weightInRangeFromPtMinSecSubtractedUsingCocktail[i][4]*weightInRangeFromPtMinSecSubtractedUsingCocktailRelErr[i][4]);
    }


    for(Int_t i=0; i < nBinsPtFine; i++){
        histoDataMCRatioRinPtBinScaledToGasFine[i] = (TH1F*)histoRinPtBinDataScaledToGasRebinFine[i]->Clone(Form("histoDataMCRatioRinPtBinScaledToGasFine%d",i));
        histoDataMCRatioRinPtBinScaledToGasFine[i]->Divide(histoRinPtBinDataScaledToGasRebinFine[i] ,histoRinPtBinMCScaledToGasRebinFine[i] ,1.,1.,"");
    }

    for(Int_t i=0; i < nBinsR; i++){
        histoWeightsEachRPtMin[i] = new TH1F(Form("histoWeightsEachRPtMin%i",i), Form("histoWeightsEachRPtMin%i",i), nBinsPtMin,  arrayBinsPtMin);
        for(Int_t j=0; j < nBinsPtMin; j++){
            histoWeightsEachRPtMin[i]->SetBinContent(j+1,histoDataMCRatioRinPtBinScaledToGasFine[j]->GetBinContent(i+1));
            histoWeightsEachRPtMin[i]->SetBinError(j+1,histoDataMCRatioRinPtBinScaledToGasFine[j]->GetBinError(i+1));
        }
    }

    //-AM    Normalization included here
    GammaScalingHistogramm(histoRTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRTruePrimMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRTrueSecMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPtTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPtMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPtMCRebin,normFactorReconstMC);
    GammaScalingHistogramm(histoEtaMC,normFactorReconstMC);

    GammaScalingHistogramm(histoRMCRTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoPtMCPtTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoConvRMC,normFactorReconstMC);
    GammaScalingHistogramm(histoConvPtMC,normFactorReconstMC);

    GammaScalingHistogramm(histoPtTrueSecFromK0S,normFactorReconstMC);
    GammaScalingHistogramm(histoPtTrueSecRest,normFactorReconstMC);
    GammaScalingHistogramm(histoPtTrueSecAllsources,normFactorReconstMC);

    GammaScalingHistogramm(histoPi0DalRTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoEtaDalRTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoCombRTrueMC,normFactorReconstMC);

    GammaScalingHistogramm(histoEtaData,normFactorReconstData);
    GammaScalingHistogramm(histoPtData,normFactorReconstData);
    GammaScalingHistogramm(histoPtDataRebin,normFactorReconstData);
    GammaScalingHistogramm(histoRData,normFactorReconstData);
    GammaScalingHistogramm(histoRDataRebin,normFactorReconstData);
    GammaScalingHistogramm(histoRMC,normFactorReconstMC);
    GammaScalingHistogramm(histoRMCRebin,normFactorReconstMC);
    GammaScalingHistogramm(histoPtSumRMCSecSubtracted,normFactorReconstMC);
    GammaScalingHistogramm(histoPtSumRBinDataSecSubtractedUsingCocktail,normFactorReconstData);
    GammaScalingHistogramm(histoPtSumRMCSecSubtractedRebin,normFactorReconstMC);
    GammaScalingHistogramm(histoPtSumRBinDataSecSubtractedUsingCocktailRebin,normFactorReconstData);

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
        GammaScalingHistogramm(histoPtinRBinTruePrimMC[i],normFactorReconstMC);
        GammaScalingHistogramm(histoPtinRBinTrueSecMC[i],normFactorReconstMC);
    }

    TH1F * histoPtDataToMCSecSubtractedRebin = (TH1F*)histoPtSumRBinDataSecSubtractedUsingCocktailRebin->Clone("histoPtDataToMCSecSubtractedRebin");
    histoPtDataToMCSecSubtractedRebin->Sumw2();
    histoPtDataToMCSecSubtractedRebin->Divide(histoPtSumRBinDataSecSubtractedUsingCocktailRebin,histoPtSumRMCSecSubtractedRebin);



    //__________________ calculation of efficiencies __________________________
    histoEffiR = (TH1F*)histoRMCRTrueMC->Clone("histoEffiR");
    histoEffiR->Sumw2();
    histoEffiR->Divide(histoRMCRTrueMC,histoConvRMC,1.,1.,"B");

    histoEffiPhi = (TH1F*)histoPhiTrueMC->Clone("histoEffiPhi");
    histoEffiPhi->Sumw2();
    histoEffiPhi->Divide(histoPhiTrueMC,histoConvPhiMC,1.,1.,"B");

    histoEffiEta = (TH1F*)histoEtaTrueMC->Clone("histoEffiEta");
    histoEffiEta->Sumw2();
    histoEffiEta->Divide(histoEtaTrueMC,histoConvEtaMC,1.,1.,"B");

    Int_t nBinsPtNew = 63;
    Double_t arrPtBins[nBinsPtNew];
    for(Int_t i=0; i<nBinsPtNew+1; i++){
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
    TH1F * histoConvPtMCRebin     = (TH1F *) histoConvPtMC->Rebin(nBinsPtNew,"histoConvPtMCRebin",arrPtBins);
    histoEffiPt                   = (TH1F*)histoPtMCPtTrueMCRebin->Clone("histoEffiPt");
    histoEffiPt->Sumw2();
    histoEffiPt->Divide(histoPtMCPtTrueMCRebin,histoConvPtMCRebin,1.,1.,"B");


    //_______________________ Photon purity __________________________
    //versus R, all photons
    histoPurityR = (TH1F*)histoRTrueMC->Clone("histoPurityR");
    histoPurityR->Sumw2();
    histoPurityR->Divide(histoRTrueMC,histoRMC,1.,1.,"B");

    //versus R, primary photons
    histoPurityPrimR = (TH1F*)histoRTruePrimMC->Clone("histoPurityPrimR");
    histoPurityPrimR->Sumw2();
    histoPurityPrimR->Divide(histoRTruePrimMC,histoRMC,1.,1.,"B");

    //versus R, primary photons above 0.4 GeV/c
    histoPurityPrimRPtBin3 = (TH1F*)histoRinPtBinTruePrimMC[2]->Clone("histoPurityPrimRPtBin3");
    histoPurityPrimRPtBin3->Sumw2();
    histoPurityPrimRPtBin3->Divide(histoRinPtBinTruePrimMC[2],histoRinPtBinMC[2],1.,1.,"B");

    //secondary photons
    histoPuritySecR = (TH1F*)histoRTrueSecMC->Clone("histoPuritySecR");
    histoPuritySecR->Sumw2();
    histoPuritySecR->Divide(histoRTrueSecMC,histoRMC,1.,1.,"B");

    //secondary photons above 0.4 GeV/c
    histoPuritySecRPtBin3 = (TH1F*)histoRinPtBinTrueSecMC[2]->Clone("histoPuritySecRPtBin3");
    histoPuritySecRPtBin3->Sumw2();
    histoPuritySecRPtBin3->Divide(histoRinPtBinTrueSecMC[2],histoRinPtBinMC[2],1.,1.,"B");

    //verus pT
    histoPurityPt = (TH1F*)histoPtTrueMC->Clone("histoPurityPt");
    histoPurityPt->Sumw2();
    histoPurityPt->Divide(histoPtTrueMC,histoPtMC,1.,1.,"B");

    //versus pT, R > 5cm
    histoPurityPt5cm = (TH1F*)histoPt5cmTrueMC->Clone("histoPurityPt5cm");
    histoPurityPt5cm->Sumw2();
    histoPurityPt5cm->Divide(histoPt5cmTrueMC,histoPt5cmMC,1.,1.,"B");

    histoPurityPrimPt5cm = (TH1F*)histoPt5cmTruePrimMC->Clone("histoPurityPrimPt5cm");
    histoPurityPrimPt5cm->Sumw2();
    histoPurityPrimPt5cm->Divide(histoPt5cmTruePrimMC,histoPt5cmMC,1.,1.,"B");

    histoPuritySecPt5cm = (TH1F*)histoPt5cmTrueSecMC->Clone("histoPuritySecPt5cm");
    histoPuritySecPt5cm->Sumw2();
    histoPuritySecPt5cm->Divide(histoPt5cmTrueSecMC,histoPt5cmMC,1.,1.,"B");




    //__________________ Comparison Data-MC WITHOUT scaling to gas _______________________
    //- AM  Normalization to number of events should be inserted here

    histoDataMCRatioR            = (TH1F*)histoRData->Clone("histoDataMCRatioR");
    histoDataMCRatioR->Divide(histoRData,histoRMC,1.,1.,"");

    histoDataMCRatioRRebin       = (TH1F*)histoRDataRebin->Clone("histoDataMCRatioRRebin");
    histoDataMCRatioRRebin->Divide(histoRDataRebin,histoRMCRebin,1.,1.,"");

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

    cout << "\nPlotting histograms..." << endl;

    //______________________________________ Multiplicity _________________________________________
    TCanvas * canvasNTracks = new TCanvas("canvasNTracks","",1200,1000);
    DrawGammaCanvasSettings( canvasNTracks, 0.1, 0.03, 0.05, 0.09);
    canvasNTracks->SetLogy(1);

    TH2F * histoDummyNTracks = new TH2F("histoDummyNTracks","histoDummyNTracks",1000,0.,2000.,1000,1.e-10,10);
    SetStyleHistoTH2ForGraphs(histoDummyNTracks, "Good TPC tracks","Counts", 0.035,0.04,0.035,0.04,1.,1.);
    if ( optionEnergy.CompareTo("13TeV") == 0 || optionEnergy.Contains("5TeV") ) histoDummyNTracks->GetXaxis()->SetRangeUser(0.,100);
    histoDummyNTracks->DrawCopy();

    DrawGammaSetMarker(histoGoodESDTracksData, 20, markerSize, colorData, colorData);
    histoGoodESDTracksData->Draw("same,hist");
    DrawGammaSetMarker(histoGoodESDTracksMC, 20, markerSize, colorMC, colorMC);
    histoGoodESDTracksMC->Draw("same,hist");
    
    if( histoGoodESDTracksWeightedMC != 0x0) histoGoodESDTracksWeightedMC->Draw("same,hist");

    TLegend* legend = GetAndSetLegend(0.5,0.75,2.5);
    legend->AddEntry(histoGoodESDTracksData,"Data","l");
    legend->AddEntry(histoGoodESDTracksMC,generatorName.Data(),"l");
    if( histoGoodESDTracksWeightedMC != 0x0) legend->AddEntry(histoGoodESDTracksWeightedMC,Form("%s, mult. weighted",generatorName.Data()),"l");
    legend->Draw();
    
    histoDummyNTracks->Draw("same,axis");
    canvasNTracks->Print(Form("%s/NGoodESDTracks%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //______________________________________ Pt _________________________________________
    TCanvas * canvasPt = new TCanvas("canvasPt","",1200,1000);
    DrawGammaCanvasSettings( canvasPt, 0.1, 0.03, 0.05, 0.09);
    canvasPt->SetLogy(1);
    canvasPt->SetLogx(1);

    TH2F * histoDummyPt = new TH2F("histoDummyPt","histoDummyPt",1000,0.,15.,1000,1.e-10,1.e-1);
    SetStyleHistoTH2ForGraphs(histoDummyPt, "#it{p}_{T} (GeV/#it{c})","#frac{Counts}{(N_{ev.}*<N_{ch.}>)}", 0.035,0.04,0.035,0.04,1.,1.);
    histoDummyPt->DrawCopy();

    DrawGammaSetMarker(histoPtDataRebin, 20, markerSize, colorData, colorData);
    histoPtDataRebin->Draw("same,hist");
    DrawGammaSetMarker(histoPtMCRebin, 20, markerSize, colorMC, colorMC);
    histoPtMCRebin->Draw("same,hist");
    DrawGammaLines(0.4, 0.4,1.e-10, 1.e-1,1.1,kGray+1,2);
 
        TLegend* legendPt = GetAndSetLegend(0.75,0.75,2);
        legendPt->AddEntry(histoPtData,"Data","l");
        legendPt->AddEntry(histoPtMC,generatorName.Data(),"l");
        legendPt->Draw();

    histoDummyPt->Draw("same,axis");
    canvasPt->Print(Form("%s/pT%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    //______________________________________ Pt Ratio_________________________________________
    TCanvas * canvasPtRatio = new TCanvas("canvasPtRatio","",1200,1000);
    DrawGammaCanvasSettings( canvasPtRatio, 0.1, 0.03, 0.05, 0.09);
    //    canvasPtRatio->SetLogy(1);
    canvasPtRatio->SetLogx(1);

    TH2F * histoDummyPtRatio = new TH2F("histoDummyPtRatio","histoDummyPtRatio",1000,0.,15.,1000,0.,5.);
    SetStyleHistoTH2ForGraphs(histoDummyPtRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{MC}", 0.035,0.04,0.035,0.04,1.,1.);

    histoDummyPtRatio->GetYaxis()->SetRangeUser(0.5,4.); 
    histoDummyPtRatio->DrawCopy();
    TH1F * histoPtDataToMC = (TH1F*)histoPtDataRebin->Clone("histoPtDataToMC");
    histoPtDataToMC->Sumw2();
    histoPtDataToMC->Divide(histoPtDataRebin,histoPtMCRebin);

    DrawGammaSetMarker(histoPtDataToMC, 20, markerSize, colorData, colorData);
 
    histoPtDataToMC->Draw("same,hist");
    DrawGammaLines(0., 10.,1., 1.,1.1,kGray+1,2);
    DrawGammaLines(0.4, 0.4,0.5, 4.,1.1,kGray+1,2);
 

        TLegend* legendPtRatio = GetAndSetLegend(0.75,0.75,2);
        legendPtRatio->AddEntry(histoPtDataToMC,"Data/MC","l");
	//       legendPtRatio->AddEntry(histoPtMC,generatorName.Data(),"l");
        legendPtRatio->Draw();

    histoDummyPtRatio->Draw("same,axis");
    canvasPtRatio->Print(Form("%s/pTRatio%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //___________________________________ dEdx vs R  ____________________________
    TCanvas * canvasElecdEdxR = new TCanvas("canvasElecdEdxR","",480,440);
    DrawGammaCanvasSettings( canvasElecdEdxR, 0.09, 0.105, 0.02, 0.085);
    canvasElecdEdxR->SetLogz(1);
    canvasElecdEdxR->cd();
        SetStyleHistoTH2ForGraphs(  histoRElecdEdxData, "electron dEdx","R (cm)", 0.035, 0.04, 0.035, 0.04, 0.9, 1.0, 510, 510);
        histoRElecdEdxData->GetYaxis()->SetRangeUser(0,200.);
        histoRElecdEdxData->GetXaxis()->SetRangeUser(0.,200.);
        histoRElecdEdxData->Draw("colz");
    canvasElecdEdxR->Print(Form("%s/PhotonConvRElecdEdx%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //___________________________________ nsigma dEdx vs R  ____________________________
    TCanvas * canvasElecNSdEdxR = new TCanvas("canvasElecNSdEdxR","",480,440);
    DrawGammaCanvasSettings( canvasElecNSdEdxR, 0.09, 0.105, 0.02, 0.085);
    canvasElecNSdEdxR->SetLogz(1);
    canvasElecNSdEdxR->cd();
        SetStyleHistoTH2ForGraphs(  histoRElecNSdEdxData, "electron n#sigma dEdx","R (cm)", 0.035, 0.04, 0.035, 0.04, 0.9, 1.0, 510, 510);
        histoRElecNSdEdxData->GetYaxis()->SetRangeUser(0,200.);
        histoRElecNSdEdxData->GetXaxis()->SetRangeUser(-10., 10.);
        histoRElecNSdEdxData->Draw("colz");
    canvasElecNSdEdxR->Print(Form("%s/PhotonConvRElecNSdEdx%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //___________________________________ R distrib vs Phi  ____________________________
    TCanvas * canvasPhiR = new TCanvas("canvasPhiR","",480,440);
    DrawGammaCanvasSettings( canvasPhiR, 0.09, 0.105, 0.02, 0.085);
    canvasPhiR->SetLogz(1);
    canvasPhiR->cd();
        SetStyleHistoTH2ForGraphs(  histoRPhiData, "#varphi (rad)","R (cm)", 0.035, 0.04, 0.035, 0.04, 0.9, 1.0, 510, 510);
        histoRPhiData->GetYaxis()->SetRangeUser(0,200.);
        histoRPhiData->GetXaxis()->SetRangeUser(0.,2*TMath::Pi());
        histoRPhiData->Draw("colz");
    canvasPhiR->Print(Form("%s/PhotonConvRPhi%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //___________________________________ R distrib vs Eta  ____________________________
    TCanvas * canvasEtaR = new TCanvas("canvasEtaR","",480,440);
    DrawGammaCanvasSettings( canvasEtaR, 0.09, 0.105, 0.02, 0.085);
    canvasEtaR->SetLogz(1);
    canvasEtaR->cd();
        SetStyleHistoTH2ForGraphs(  histoREtaData, "#eta","R (cm)", 0.035, 0.04, 0.035, 0.04, 0.9, 1.0, 510, 510);
        histoREtaData->GetYaxis()->SetRangeUser(0,200.);
        histoREtaData->GetXaxis()->SetRangeUser(-2., 2.);
        histoREtaData->Draw("colz");
    canvasEtaR->Print(Form("%s/PhotonConvREta%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //___________________________________ R distrib vs Z  ____________________________
    TCanvas * canvasZR = new TCanvas("canvasZR","",480,440);
    DrawGammaCanvasSettings( canvasZR, 0.09, 0.105, 0.02, 0.085);
    canvasZR->SetLogz(1);
    canvasZR->cd();
        SetStyleHistoTH2ForGraphs(  histoRZData, "Z (cm)","R (cm)", 0.035, 0.04, 0.035, 0.04, 0.9, 1.0, 510, 510);
        histoRZData->GetYaxis()->SetRangeUser(0,200.);
        histoRZData->GetXaxis()->SetRangeUser(-180., 180.);
        histoRZData->Draw("colz");
    canvasZR->Print(Form("%s/PhotonConvRZ%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //___________________________________ R distrib vs pT  ____________________________
    TCanvas * canvasPtR = new TCanvas("canvasPtR","",480,440);
    DrawGammaCanvasSettings( canvasPtR, 0.09, 0.105, 0.02, 0.085);
    canvasPtR->SetLogz(1);
    canvasPtR->cd();
        SetStyleHistoTH2ForGraphs(  histoRPtData, "#it{p}_{T} (GeV/#it{c})","R (cm)", 0.035, 0.04, 0.035, 0.04, 0.9, 1.0, 510, 510);
        histoRPtData->GetYaxis()->SetRangeUser(0,200.);
        histoRPtData->GetXaxis()->SetRangeUser(0., 20.);
        histoRPtData->Draw("colz");
    canvasPtR->Print(Form("%s/PhotonConvRPt%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //___________________________________ R distrib scaled to gas ____________________________
    ReturnCorrectValuesForCanvasScaling(1200, 1000, 1, 2, 0.1, 0.025, 0.025, 0.08, arrayX, arrayY, relX, relY,kFALSE);
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
    if (optionEnergy.CompareTo("13TeV") == 0 || optionEnergy.Contains("5TeV")) histoDummyTwoPanelsUp->GetYaxis()->SetRangeUser(-0.1,2.5);
    SetStyleHistoTH2ForGraphs(histoDummyTwoPanelsUp, "R (cm)","Counts", 0.9*textsizeLabelsUp, textsizeLabelsUp,0.9*textsizeLabelsUp,textsizeLabelsUp, 1,0.15/(textsizeFacUp*marginXRatio));

    TH2F *histoDummyTwoPanelsDown =  new TH2F("histoDummyTwoPanelsDown","histoDummyTwoPanelsDown",1000,0.,180.,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histoDummyTwoPanelsDown, "R (cm)","#frac{Data}{MC} ", 0.9*textsizeLabelsDown, textsizeLabelsDown,0.9*textsizeLabelsDown,textsizeLabelsDown, 1,0.15/(textsizeFacDown*marginXRatio));
    histoDummyTwoPanelsDown->GetYaxis()->SetRangeUser(0.3,1.55);
    if (optionEnergy.CompareTo("13TeV") == 0 || optionEnergy.Contains("5TeV")) histoDummyTwoPanelsDown->GetYaxis()->SetRangeUser(0.8,1.55);

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
        DrawGammaSetMarker(histoDataMCRatioRinPtBinScaledToGasPtBin1, 25, markerSize, color[2], color[2]);
        DrawGammaSetMarker(histoDataMCRatioRinPtBinScaledToGasPtBin2, 25, markerSize, color[4], color[4]);
        DrawGammaSetMarker(histoDataMCRatioRinPtBinScaledToGasPtBin3, 25, markerSize, color[6], color[6]);

        histoDataMCRatioRScaledToGas->Draw("same,pE");
        histoDataMCRatioRinPtBinScaledToGasPtBin1->Draw("same,pE");
        histoDataMCRatioRinPtBinScaledToGasPtBin2->Draw("same,pE");
        histoDataMCRatioRinPtBinScaledToGasPtBin3->Draw("same,pE");

        TLegend* legenLowPanel = GetAndSetLegend2(0.5, 0.85-(4*0.9*textsizeLabelsDown), 0.7, 0.85, textSizeLabels);
        TString legendWithPeriod = Form("From proj, Period %s",optionPeriod.Data());

        legenLowPanel->AddEntry(histoDataMCRatioRScaledToGas,"full #it{p}_{T} range");//legendWithPeriod.Data(),"lp");
        legenLowPanel->AddEntry(histoDataMCRatioRinPtBinScaledToGasPtBin1,"#it{p}_{T} > 0.1 GeV/#it{c}","lp");
        legenLowPanel->AddEntry(histoDataMCRatioRinPtBinScaledToGasPtBin2,"#it{p}_{T} > 0.3 GeV/#it{c}","lp");
        legenLowPanel->AddEntry(histoDataMCRatioRinPtBinScaledToGasPtBin3,"#it{p}_{T} > 0.4 GeV/#it{c}","lp");
        legenLowPanel->Draw();

    histoDummyTwoPanelsDown->Draw("same,axis");
    canvasTwoPanels->Print(Form("%s/PhotonConvRScaledToGas%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    //---------------------------
    // Weights vs R . Different pT Min sec. subtracted
    //---------------------------

    TCanvas * canvasWeightsSecSub_singlepad = new TCanvas("canvasWeightsSecSub_singlepad","",1200,1000);
    DrawGammaCanvasSettings( canvasWeightsSecSub_singlepad, 0.095, 0.02, 0.02, 0.08);

    TH2F *histoDummyWeightsSecSub =  new TH2F("histoDummyWeightsSecSub","histoDummyWeightsSecSub",1000,0.,180.,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histoDummyWeightsSecSub, "R (cm)","#frac{Data}{MC} (Sec. sub.)", 0.9*textsizeLabelsDown, textsizeLabelsDown,0.9*textsizeLabelsDown,textsizeLabelsDown, 1,0.15/(textsizeFacDown*marginXRatio));
    histoDummyWeightsSecSub->GetYaxis()->SetRangeUser(0.3,1.55);
    if (optionEnergy.CompareTo("13TeV") == 0 || optionEnergy.Contains("5TeV")) histoDummyWeightsSecSub->GetYaxis()->SetRangeUser(0.8,1.55);

    histoDummyWeightsSecSub->DrawCopy();

    DrawGammaLines(0.,180,1., 1.,1.,kGray,1);
    DrawGammaLines(0.,180,1.15, 1.15,1.,kGray,4);
    DrawGammaLines(0.,180,1.1, 1.1,1.,kGray,2);
    DrawGammaLines(0.,180,1.05, 1.05,1.,kGray,3);
    DrawGammaLines(0.,180,0.95, 0.95,1.,kGray,3);
    
    DrawGammaSetMarker(histoDataMCRatioRScaledToGasSecSub, 20, markerSize, colorData, colorData);
    DrawGammaSetMarker(histoDataMCRatioRScaledToGasSecSubPtBin1, 25, markerSize, color[2], color[2]);
    DrawGammaSetMarker(histoDataMCRatioRScaledToGasSecSubPtBin2, 25, markerSize, color[4], color[4]);
    DrawGammaSetMarker(histoDataMCRatioRScaledToGasSecSubPtBin3, 25, markerSize, color[6], color[6]);
    
    histoDataMCRatioRScaledToGasSecSub->Draw("same,pE");
    histoDataMCRatioRScaledToGasSecSubPtBin1->Draw("same,pE");
    histoDataMCRatioRScaledToGasSecSubPtBin2->Draw("same,pE");
    histoDataMCRatioRScaledToGasSecSubPtBin3->Draw("same,pE");
    
    TLegend* legenWeightSecSub = GetAndSetLegend2(0.5, 0.85-(4*0.9*textsizeLabelsDown), 0.7, 0.85, textSizeLabels);
    TString legenWeightSecSubWithPeriod = Form("From proj, Period %s",optionPeriod.Data());
	
    legenWeightSecSub->AddEntry(histoDataMCRatioRScaledToGasSecSub,"full #it{p}_{T} range");//legendWithPeriod.Data(),"lp");
    legenWeightSecSub->AddEntry(histoDataMCRatioRScaledToGasSecSubPtBin1,"#it{p}_{T} > 0.1 GeV/#it{c}","lp");
    legenWeightSecSub->AddEntry(histoDataMCRatioRScaledToGasSecSubPtBin2,"#it{p}_{T} > 0.3 GeV/#it{c}","lp");
    legenWeightSecSub->AddEntry(histoDataMCRatioRScaledToGasSecSubPtBin3,"#it{p}_{T} > 0.4 GeV/#it{c}","lp");
    legenWeightSecSub->Draw();
    
    histoDummyWeightsSecSub->Draw("same,axis");
    canvasWeightsSecSub_singlepad->Print(Form("%s/PhotonConvRScaledToGasSecSub%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //---------------------------

    TCanvas * canvasComparisonOldMethod = new TCanvas("canvasComparisonOldMethod","",1200,1000);
    DrawGammaCanvasSettings( canvasComparisonOldMethod, 0.095, 0.02, 0.02, 0.08);

    TH2F * histoDummyComparisonOld = new TH2F("histoDummyComparisonOld","histoDummyComparisonOld",1000,0.,180.,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histoDummyComparisonOld, "R (cm)","#frac{Data}{MC} ",0.03,0.035,0.03,0.035,1.,1.);
    if (optionEnergy.CompareTo("13TeV") == 0 || optionEnergy.Contains("5TeV")) histoDummyComparisonOld->GetYaxis()->SetRangeUser(0.8,1.55);
    histoDummyComparisonOld->DrawCopy();

        DrawGammaLines(0.,180,1., 1.,1.,kGray,1);
        DrawGammaLines(0.,180,1.15, 1.15,1.,kGray,4);
        DrawGammaLines(0.,180,1.1, 1.1,1.,kGray,2);
        DrawGammaLines(0.,180,1.05, 1.05,1.,kGray,3);
        DrawGammaLines(0.,180,0.95, 0.95,1.,kGray,3);

        DrawGammaSetMarker(histoDataMCRatioRRebin, 20, markerSize, kGray+1, kGray+1);
        histoDataMCRatioRScaledToGas->Draw("same,pE");
        histoDataMCRatioRinPtBinScaledToGasPtBin3->Draw("same,pE");
        histoDataMCRatioRRebin->Draw("same,pE");

        TLegend* legenOldMethod = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legenOldMethod->AddEntry(histoDataMCRatioRScaledToGas,legendWithPeriod.Data(),"lp");
        legenOldMethod->AddEntry(histoDataMCRatioRinPtBinScaledToGasPtBin3,"#it{p}_{T} > 0.4 GeV/#it{c}","lp");
        legenOldMethod->AddEntry(histoDataMCRatioRRebin,"Data and MC norm. to  their (N_{ev.}*<N_{ch.}>)","lp");
        legenOldMethod->Draw();

    histoDummyComparisonOld->Draw("same,axis");
    canvasComparisonOldMethod->Print(Form("%s/PhotonConvRScaledToGasVSOldMethod%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //______________________________________ R distribution ______________________________
    textSizeLabels = 50;
    ReturnCorrectValuesForCanvasScaling(2500,2000, 2, 2,0.06, 0.025, 0.025,0.06,arrayXRdistrib,arrayYRdistrib,relXRdistrib,relYRdistrib,kFALSE);
    TCanvas* canvasPhotonR = new TCanvas("canvasPhotonR","",0,0,2500,2000);
    canvasPhotonR->cd();

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

        for(Int_t iR = 1; iR < nBinsR; iR++) DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.,kGray);

        DrawGammaSetMarker(histoRData, 20, markerSize, colorData, colorData);
        histoRData->Draw("same,hist");
        DrawGammaSetMarker(histoRMC, 20, markerSize, colorMC, colorMC);
        histoRMC->Draw("same,hist");
        DrawGammaSetMarker(histoRTrueMC, 20, markerSize, colorTrueMC, colorTrueMC);
        histoRTrueMC->Draw("same,hist");
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
        DrawGammaSetMarker(histoRTrueSecMC, 20, markerSize, color[5], color[5]);
	histoRTrueSecMC->SetLineWidth(2);
        histoRTrueSecMC->Draw("same,hist");

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

        for(Int_t iR = 1; iR < nBinsR; iR++) DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.,kGray);
        DrawGammaLines(0.,180,1., 1.,1.,kGray);

        DrawGammaSetMarker(histoDataMCRatioR, 20, markerSize, colorData, colorData);
        histoDataMCRatioR->Draw("same,histo");

    histoDummy4PanelsDown->Draw("axis,same");
    padRDistribUpperRight->cd();
    padRDistribUpperRight->SetLogy();
    histoDummy4PanelsUp->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++) DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.,kGray);

        histoRData->Draw("same,histo");
        histoRMC->Draw("same,histo");
        histoRTrueMC->Draw("same,histo");

        for(Int_t i=0; i<6; i++){
            DrawGammaSetMarker(histoRinPtBinData[i], 20, markerSize, color[i], color[i]);
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
	// AM- pT limits corrected
        TLegend* legenRdistribPtcuts = GetAndSetLegend2(0.05, 0.35-(6.2*0.9*textsizeLabelsUp), 0.35, 0.35, textSizeLabels);
        legenRdistribPtcuts->AddEntry(histoRData,"Data","l");
        legenRdistribPtcuts->AddEntry(histoRMC,"MC","l");
        legenRdistribPtcuts->AddEntry(histoRTrueMC,"True MC","l");
        legenRdistribPtcuts->AddEntry(histoRinPtBinData[1],"Data,  #it{p}_{T} > 0.3 GeV/#it{c}","l");
        legenRdistribPtcuts->AddEntry(histoRinPtBinMC[1],"MC,  #it{p}_{T} >0.3 GeV/#it{c}","lf");
        legenRdistribPtcuts->AddEntry(histoRinPtBinTrueMC[1],"True MC,  #it{p}_{T} > 0.3 GeV/#it{c}","fl");
        legenRdistribPtcuts->Draw();

    histoDummy4PanelsUp->Draw("axis,same");
    padRDistribLowerRight->cd();

        for(Int_t iR = 1; iR < nBinsR; iR++) DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.1,kGray);
        DrawGammaLines(0.,180,1., 1.,1.,kGray);

        histoDummy4PanelsDown->DrawCopy();
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
        DrawGammaSetMarker(histoRTruePrimMC, 20, markerSize, color[3], color[3]);
        histoRTruePrimMC->SetLineWidth(2);
        histoRTruePrimMC->Draw("same,hist");
        DrawGammaSetMarker(histoRTrueSecMC, 20, markerSize, color[5], color[5]);
        histoRTrueSecMC->SetLineWidth(2);
        histoRTrueSecMC->Draw("same,hist");
        histoCombRTrueMC->SetLineWidth(2);
        histoCombRTrueMC->DrawCopy("same,hist");
        histoPi0DalRTrueMC->SetLineWidth(2);
        histoPi0DalRTrueMC->DrawCopy("same,hist");
        histoEtaDalRTrueMC->SetLineWidth(2);
        histoEtaDalRTrueMC->DrawCopy("same,hist");

        TLegend* legenRdistribSingle = GetAndSetLegend2(0.6, 0.93-(8*0.85*textsizeLabelsUp), 0.9, 0.93, textSizeLabels);
        legenRdistribSingle->AddEntry(histoRData,"Data","l");
        legenRdistribSingle->AddEntry(histoRMC,"MC","l");
        legenRdistribSingle->AddEntry(histoRTrueMC,"True MC","l");
        legenRdistribSingle->AddEntry(histoRTruePrimMC,"True prim. MC","l");
        legenRdistribSingle->AddEntry(histoRTrueSecMC,"True sec. MC","l");
        legenRdistribSingle->AddEntry(histoCombRTrueMC,"True MC comb.","l");
        legenRdistribSingle->AddEntry(histoPi0DalRTrueMC,"True MC #pi^{0} Dal.","lf");
        legenRdistribSingle->AddEntry(histoEtaDalRTrueMC,"True MC #eta Dal.","fl");
        legenRdistribSingle->Draw();

    histoDummySinglePad->Draw("axis,same");
    canvasPhotonR_siglepad->Print(Form("%s/PhotonConvRsingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    textSizeLabels = 30;
    ReturnCorrectValuesForCanvasScaling(1300, 1000, 1, 2, 0.08, 0.02, 0.02, 0.05, arrayXRConv2Pad, arrayYRConv2Pad, relXRConv2Pad, relYRConv2Pad,kFALSE);
    TCanvas * canvasRConv2Pad = new TCanvas("canvasRConv2Pad","",1200,1200);

    TPad* padUpperRConv2Pad = new TPad("","",arrayXRConv2Pad[0], arrayYRConv2Pad[1], arrayXRConv2Pad[1], arrayYRConv2Pad[0],-1, -1, -2);
    DrawGammaPadSettings( padUpperRConv2Pad, relXRConv2Pad[0], relXRConv2Pad[2], relYRConv2Pad[0], relYRConv2Pad[1]);
    TPad* padLowerRConv2Pad = new TPad("","",arrayXRConv2Pad[0], arrayYRConv2Pad[2], arrayXRConv2Pad[1], arrayYRConv2Pad[1],-1, -1, -2);
    DrawGammaPadSettings( padLowerRConv2Pad, relXRConv2Pad[0], relXRConv2Pad[2], relYRConv2Pad[1], relYRConv2Pad[2]);
    marginXRatio = relXRConv2Pad[0]*1200;
    if (padUpperRConv2Pad->XtoPixel(padUpperRConv2Pad->GetX2()) < padUpperRConv2Pad->YtoPixel(padUpperRConv2Pad->GetY1())){
        textsizeLabelsUpRConv2Pad = (Double_t)textSizeLabels/padUpperRConv2Pad->XtoPixel(padUpperRConv2Pad->GetX2()) ;
        textsizeFacUpRConv2Pad = (Double_t)1./padUpperRConv2Pad->XtoPixel(padUpperRConv2Pad->GetX2()) ;
    } else {
        textsizeLabelsUpRConv2Pad = (Double_t)textSizeLabels/padUpperRConv2Pad->YtoPixel(padUpperRConv2Pad->GetY1());
        textsizeFacUpRConv2Pad = (Double_t)1./padUpperRConv2Pad->YtoPixel(padUpperRConv2Pad->GetY1());
    }
    padUpperRConv2Pad->Draw();
    if (padLowerRConv2Pad->XtoPixel(padLowerRConv2Pad->GetX2()) < padLowerRConv2Pad->YtoPixel(padLowerRConv2Pad->GetY1())){
        textsizeLabelsDownRConv2Pad = (Double_t)textSizeLabels/padLowerRConv2Pad->XtoPixel(padLowerRConv2Pad->GetX2()) ;
        textsizeFacDownRConv2Pad = (Double_t)1./padLowerRConv2Pad->XtoPixel(padLowerRConv2Pad->GetX2()) ;
    } else {
        textsizeLabelsDownRConv2Pad = (Double_t)textSizeLabels/padLowerRConv2Pad->YtoPixel(padLowerRConv2Pad->GetY1());
        textsizeFacDownRConv2Pad = (Double_t)1./padLowerRConv2Pad->YtoPixel(padLowerRConv2Pad->GetY1());
    }
    padLowerRConv2Pad->Draw();

    TH2F *histoDummyRConv2PadR = new TH2F("histoDummyRConv2PadR","histoDummyRConv2PadR",1000,0.,184.,1000,1.5e-9,1.);
    SetStyleHistoTH2ForGraphs(histoDummyRConv2PadR, "R (cm)","Counts",0.85*textsizeLabelsUp, textsizeLabelsUp,0.85*textsizeLabelsUp,textsizeLabelsUp,.9, .7);
    histoDummyRConv2PadR->GetYaxis()->SetRangeUser(1.5e-9,2e-3);

    TH2F *histoDummyRConv2PadPt =  new TH2F("histoDummyRConv2PadPt","histoDummyRConv2PadPt",1000,0.,184.,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histoDummyRConv2PadPt,"R (cm)","#frac{Data}{MC}",0.85*textsizeLabelsDown, textsizeLabelsDown,0.85*textsizeLabelsDown,textsizeLabelsDown, 0.9,0.7);
    histoDummyRConv2PadPt->GetYaxis()->SetRangeUser(0.35,1.65);

    padUpperRConv2Pad->cd();
    padUpperRConv2Pad->SetLogy();
    histoDummyRConv2PadR->Draw("copy");

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
	histoRTrueSecMC->DrawCopy("same,hist");

        TLegend* legenRdistrib2 = GetAndSetLegend2(0.6, 0.93-(8*0.85*textsizeLabelsUp), 0.9, 0.93, textSizeLabels);
        legenRdistrib2->AddEntry(histoRData,"Data","l");
        legenRdistrib2->AddEntry(histoRMC,"MC","l");
        legenRdistrib2->AddEntry(histoRTrueMC,"True MC","l");
        legenRdistrib2->AddEntry(histoRTrueSecMC,"True Sec MC","l");
        legenRdistrib2->AddEntry(histoCombRTrueMC,"True MC comb.","l");
        legenRdistrib2->AddEntry(histoPi0DalRTrueMC,"True MC #pi^{0} Dal.","lf");
        legenRdistrib2->AddEntry(histoEtaDalRTrueMC,"True MC #eta Dal.","fl");
        legenRdistrib2->Draw();

    histoDummyRConv2PadR->Draw("axis,same");
    padLowerRConv2Pad->cd();
    histoDummyRConv2PadPt->Draw("copy");

        for(Int_t iR = 1; iR < nBinsR; iR++) DrawGammaLines(arrayRBins[iR],arrayRBins[iR],1.5e-9,2e-3,1.,kGray);
        DrawGammaLines(0.,180,1., 1.,1.,kGray);

        DrawGammaSetMarker(histoDataMCRatioR, 20, markerSize, colorData, colorData);
        histoDataMCRatioR->Draw("same,histo");

    histoDummyRConv2PadPt->Draw("same,axis");
    canvasRConv2Pad->Print(Form("%s/PhotonRConv2Pad%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //______________________________________ Photon characteristics _________________________________
    textSizeLabels = 50;
    ReturnCorrectValuesForCanvasScaling(2500,2000, 2, 2,0.06, 0.025, 0.025,0.06,arrayXPhotonChar,arrayYPhotonChar,relXPhotonChar,relYPhotonChar,kFALSE);
    TCanvas* canvasPhotonChar = new TCanvas("canvasPhotonChar","",0,0,2500,2000);

    TPad* padPhotonCharLowerLeft = new TPad("", "", arrayXPhotonChar[0], arrayYPhotonChar[2], arrayXPhotonChar[1], arrayYPhotonChar[1],-1, -1, -2);
    DrawGammaPadSettings( padPhotonCharLowerLeft, relXPhotonChar[0], relXPhotonChar[1]-0.001, relYPhotonChar[1]-0.005, relYPhotonChar[2]);
    padPhotonCharLowerLeft->Draw();
    TPad* padPhotonCharUpperLeft = new TPad("", "", arrayXPhotonChar[0], arrayYPhotonChar[1], arrayXPhotonChar[1], arrayYPhotonChar[0],-1, -1, -2);
    DrawGammaPadSettings( padPhotonCharUpperLeft, relXPhotonChar[0], relXPhotonChar[1]-0.001, relYPhotonChar[0], relYPhotonChar[1]-0.005);
    padPhotonCharUpperLeft->Draw();
    TPad* padPhotonCharLowerRight = new TPad("", "", arrayXPhotonChar[1], arrayYPhotonChar[2], arrayXPhotonChar[2], arrayYPhotonChar[1],-1, -1, -2);
    DrawGammaPadSettings( padPhotonCharLowerRight, relXPhotonChar[1]-0.001, relXPhotonChar[2], relYPhotonChar[1]-0.005, relYPhotonChar[2]);
    padPhotonCharLowerRight->Draw();
    TPad* padPhotonCharUpperRight = new TPad("", "", arrayXPhotonChar[1], arrayYPhotonChar[1], arrayXPhotonChar[2], arrayYPhotonChar[0],-1, -1, -2);
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
    if ( optionEnergy.CompareTo("13TeV") == 0 || optionEnergy.Contains("5TeV")) histoDummyMass->GetYaxis()->SetRangeUser(1.e-6,2e-2);
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


    textSizeLabels = 30;
    ReturnCorrectValuesForCanvasScaling(1200, 1000, 1, 2, 0.1, 0.025, 0.025, 0.08, arrayX, arrayY, relX, relY,kFALSE);
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

    textSizeLabels = 30;
    ReturnCorrectValuesForCanvasScaling(1200, 1000, 1, 2, 0.1, 0.025, 0.03, 0.08, arrayX, arrayY, relX, relY,kFALSE);
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


    //_____________________________ single canvas plotting ________________________
    TCanvas* canvasV0Finder = new TCanvas("canvasV0Finder","",1300,1000);
    DrawGammaCanvasSettings( canvasV0Finder,  0.13, 0.02, 0.02, 0.09);
    TH2F * histo2DDummy = new TH2F("","",1000,0.,200.,1000,0.,2.);

    // draw efficiency vs Pt, smaller eta range
    canvasV0Finder->SetLogy(1);
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#it{p}_{T} (GeV/#it{c})","Efficiency ()",0.04,0.04, 0.04,0.04, 1.,1.);
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
    SetStyleHistoTH2ForGraphs(histo2DDummy,"R (cm)","Efficiency | #eta | < 0.8",0.04,0.04, 0.04,0.04, 1.,1.);
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
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#phi (rad)","Efficiency | #eta | < 0.8",0.04,0.04, 0.04,0.04, 1.,1.3);
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
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#eta","Efficiency | #eta | < 0.8",0.04,0.04, 0.04,0.04, 1.,1.3);
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

    TCanvas* canvasChi2 = new TCanvas("canvasChi2","",1300,1000);
    DrawGammaCanvasSettings( canvasChi2,  0.12, 0.02, 0.02, 0.12);
    canvasChi2->SetLogy(1);

    TH2F * histoDummyChi2Single = new TH2F("histoDummyChi2Single","histoDummyChi2Single",1000,0.,30.,10000,1.e-6,1.);
    SetStyleHistoTH2ForGraphs(histoDummyChi2Single, "#chi^{2}","Counts ", 0.85*textsizeLabelsDown, textsizeLabelsDown,0.85*textsizeLabelsDown,textsizeLabelsDown, 0.9,0.92);
    histoDummyChi2Single->GetYaxis()->SetRangeUser(1.e-6,2.e-2);
    histoDummyChi2Single->DrawCopy();

        DrawGammaSetMarker(histoChi2Data, 20, markerSize, colorData, colorData);
        histoChi2Data->Draw("same,hist");
		DrawGammaSetMarker(histoChi2MC, 20, markerSize, colorMC, colorMC);
        histoChi2MC->Draw("same,hist");

    histoDummyChi2Single->Draw("axis,same");
    canvasChi2->Update();
    canvasChi2->SaveAs(Form("%s/PhotonChi2Single%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));



    //______________________________________ Purity _________________________________
    textSizeLabels = 30;
    ReturnCorrectValuesForCanvasScaling(1200, 1200, 1, 2, 0.08, 0.025, 0.025, 0.08, arrayXpurity, arrayYpurity, relXpurity, relYpurity,kFALSE);
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

    TH2F *histoDummyPurityPt =  new TH2F("histoDummyPurityPt","histoDummyPurityPt",1000,0.,10.,1000,0.,1.2);
    SetStyleHistoTH2ForGraphs(histoDummyPurityPt, "#it{p}_{T} (GeV/c)","Purity", 0.9*textsizeLabelsDownpurity, textsizeLabelsDownpurity,0.9*textsizeLabelsDownpurity,textsizeLabelsDownpurity, 0.9,0.1/(textsizeFacDownpurity*marginXRatio));
    histoDummyPurityPt->GetYaxis()->SetRangeUser(0.6,1.1);

    padUpperPurity->cd();
    histoDummyPurityR->DrawCopy();

        for(Int_t iR = 1; iR < nBinsR; iR++) DrawGammaLines(arrayRBins[iR],arrayRBins[iR],0.,1.2,1.,kGray,2);
        DrawGammaLines(0.,180,1., 1.,1.,kGray+1);

        //all photons
        DrawGammaSetMarker(histoPurityR, 20, markerSize, colorData, colorData);
        histoPurityR->Draw("same,hist");
        //primary photons
        DrawGammaSetMarker(histoPurityPrimR, 20, markerSize, color[3], color[3]);
        histoPurityPrimR->SetLineStyle(2);
        histoPurityPrimR->Draw("same,hist");
        //primary photons from 0.4 GeV
        DrawGammaSetMarker(histoPurityPrimR, 20, markerSize, color[3], color[3]);
        histoPurityPrimRPtBin3->Draw("same,hist");
        //secondaries photons
        DrawGammaSetMarker(histoPuritySecR, 20, markerSize, colorData, colorData);
        histoPuritySecR->SetLineStyle(3);
//         histoPuritySecR->Draw("same,hist");
        //secondaries photons from 0.4 GeV
//         histoPuritySecRPtBin3->Draw("same,hist");

        TLegend* legenPurity = GetAndSetLegend2(0.5, 0.6-(3*0.98*textsizeLabelsUppurity), 0.8, 0.6, textSizeLabels);
        legenPurity->AddEntry(histoPurityR,"All #gamma","l");
        legenPurity->AddEntry(histoPurityPrimR,"Primary #gamma","l");
        legenPurity->AddEntry(histoPurityPrimRPtBin3,"Primary #gamma, #it{p}_{T} > 0.4 GeV/#it{c}","l");
//         legenPurity->AddEntry(histoPuritySecR,"Secondary #gamma","l");
//         legenPurity->AddEntry(histoPuritySecRPtBin3,"Secondary #gamma, #it{p}_{T} > 0.4 GeV/#it{c}","l");
        legenPurity->Draw();

    histoDummyPurityR->Draw("same,axis");
    padLowerPurity->cd();
    histoDummyPurityPt->GetXaxis()->SetRangeUser(0.04,10.);
    histoDummyPurityPt->DrawCopy();

        DrawGammaLines(0.,8.,1., 1.,1.1,kGray+1);
        DrawGammaSetMarker(histoPurityPt, 20, markerSize, colorData, colorData);
        histoPurityPt->Draw("same,hist");
        DrawGammaSetMarker(histoPurityPt5cm, 20, markerSize, kGreen+2, kGreen+2);
        histoPurityPt5cm->Draw("same,hist");
        DrawGammaSetMarker(histoPurityPrimPt5cm, 20, markerSize, kAzure+2, kAzure+2);
        histoPurityPrimPt5cm->Draw("same,hist");

        TLegend* legenPurityPt = GetAndSetLegend2(0.13, 0.2+(3*0.9*textsizeLabelsUppurity), 0.5, 0.2, textSizeLabels);
        legenPurityPt->AddEntry(histoPurityPt,"no R cut","l");
        legenPurityPt->AddEntry(histoPurityPt5cm,"R > 5 cm","l");
        legenPurityPt->AddEntry(histoPurityPrimPt5cm,"R > 5 cm, primary #gamma","l");
        legenPurityPt->Draw();

    histoDummyPurityPt->Draw("same,axis");
    canvasPurity->Print(Form("%s/PhotonPurity%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    padLowerPurity->SetLogx();
    DrawGammaLines(0.4,0.4,0.6,1.1,1.,kGray,2);
 
    canvasPurity->Print(Form("%s/PhotonPurity_logX_%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));



    TCanvas * canvasSinglePtPurity = new TCanvas("canvasSinglePtPurity","",1200,1000);
    DrawGammaCanvasSettings( canvasSinglePtPurity,  0.12, 0.02, 0.02, 0.12);

    TH2F *histoDummyPuritySinglePt =  new TH2F("histoDummyPuritySinglePt","histoDummyPuritySinglePt",1000,0.,10.,1000,0.5,1.1);
    SetStyleHistoTH2ForGraphs(histoDummyPuritySinglePt, "#it{p}_{T} (GeV/c)","Purity",0.05,0.05, 0.05,0.05);
    histoDummyPuritySinglePt->GetXaxis()->SetRangeUser(0.04,5.);
    histoDummyPuritySinglePt->DrawCopy();

        DrawGammaSetMarker(histoPurityPt5cm, 20, markerSize, colorMC, colorMC);
        histoPurityPt5cm->DrawCopy("same");

        DrawGammaSetMarker(histoPurityPrimPt5cm, 24, markerSize, kAzure+2, kAzure+2);
        histoPurityPrimPt5cm->Draw("same");

        TLegend* legenPurityPtSingle = GetAndSetLegend2(0.22, 0.2+(2.5*0.9*textsizeLabelsUppurity), 0.5, 0.2, textSizeLabels);
        legenPurityPtSingle->SetHeader("R > 5 cm");
        legenPurityPtSingle->AddEntry(histoPurityPt5cm,"All #gamma","pl");
        legenPurityPtSingle->AddEntry(histoPurityPrimPt5cm,"Primary #gamma Frac","pl");
        legenPurityPtSingle->Draw();

    histoDummyPuritySinglePt->Draw("same,axis");
    canvasSinglePtPurity->Print(Form("%s/PhotonPuritySingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    canvasSinglePtPurity->SetLogx();
    DrawGammaLines(0.4,0.4,0.6,1.1,1.,kGray,2);
    canvasSinglePtPurity->Print(Form("%s/PhotonPuritySingle_LogX_%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));



    TCanvas * canvasSingleEtaDataMC = new TCanvas("canvasSingleEtaDataMC","",1200,1000);
    DrawGammaCanvasSettings( canvasSingleEtaDataMC,  0.15, 0.02, 0.07, 0.12);
    TH2F *histoDummyEtaDataMC =  new TH2F("histoDummyEtaDataMC","histoDummyEtaDataMC",1000,-1.1,1.1,1000,0.0,1.2*histoEtaData->GetMaximum());
    SetStyleHistoTH2ForGraphs(histoDummyEtaDataMC, "#eta","#frac{dN_{#gamma}}{dEta}",0.05,0.05, 0.05,0.05);
    histoDummyEtaDataMC->DrawCopy();
        DrawGammaSetMarker(histoEtaData, 20, markerSize, colorData, colorData);
        DrawGammaSetMarker(histoEtaMC, 20, markerSize, colorMC, colorMC);
        histoEtaData->Draw("same,histo");
        histoEtaMC->Draw("same,histo");
        TLegend* legenEta = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
        legenEta->AddEntry(histoEtaData,"Data","p");
        legenEta->AddEntry(histoEtaMC,"MC","p");
        legenEta->Draw();
    histoDummyEtaDataMC->Draw("same,axis");
    canvasSingleEtaDataMC->Print(Form("%s/PhotonEtaDataMCSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasSingleAsymDataMCLowP = new TCanvas("canvasSingleAsymDataMCLowP","",1200,1000);
    DrawGammaCanvasSettings( canvasSingleAsymDataMCLowP,  0.15, 0.02, 0.07, 0.12);
    TH2F *histoDummyAsymDataMCLowP =  new TH2F("histoDummyAsymDataMCLowP","histoDummyAsymDataMCLowP",1000,0.,1.0,1000,0.0,1.2*histoAsymPDataLowP->GetMaximum());
    SetStyleHistoTH2ForGraphs(histoDummyAsymDataMCLowP, "#frac{p_{e}}{p_{#gamma}}","#frac{dN_{#gamma}}{dAsym}",0.05,0.05, 0.05,0.05);
    histoDummyAsymDataMCLowP->DrawCopy();
        DrawGammaSetMarker(histoAsymPDataLowP, 20, markerSize, colorData, colorData);
        DrawGammaSetMarker(histoAsymPTrueMCLowP, 20, markerSize, colorMC, colorMC);
        histoAsymPDataLowP->Draw("same");
        histoAsymPMCLowP->Draw("same");
        histoAsymPTrueMCLowP->Draw("same");
        TLegend* legenAsymLowP = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
        legenAsymLowP->AddEntry(histoAsymPDataLowP,"Data","p");
        legenAsymLowP->AddEntry(histoAsymPTrueMCLowP,"MC","p");
        legenAsymLowP->Draw();
    histoDummyAsymDataMCLowP->Draw("same,axis");
    canvasSingleAsymDataMCLowP->Print(Form("%s/PhotonAsymDataMCLowPSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));

    TCanvas * canvasSingleAsymDataMCHighP = new TCanvas("canvasSingleAsymDataMCHighP","",1200,1000);
    DrawGammaCanvasSettings( canvasSingleAsymDataMCHighP,  0.15, 0.02, 0.07, 0.12);
    TH2F *histoDummyAsymDataMCHighP =  new TH2F("histoDummyAsymDataMCHighP","histoDummyAsymDataMCHighP",1000,0.,1.0,1000,0.0,1.2*histoAsymPDataHighP->GetMaximum());
    SetStyleHistoTH2ForGraphs(histoDummyAsymDataMCHighP, "#frac{p_{e}}{p_{#gamma}}","#frac{dN_{#gamma}}{dAsym}",0.05,0.05, 0.05,0.05);
    histoDummyAsymDataMCHighP->DrawCopy();
        DrawGammaSetMarker(histoAsymPDataHighP, 20, markerSize, colorData, colorData);
        DrawGammaSetMarker(histoAsymPTrueMCHighP, 20, markerSize, colorMC, colorMC);
        histoAsymPDataHighP->Draw("same");
        histoAsymPMCHighP->Draw("same");
        histoAsymPTrueMCHighP->Draw("same,hist");
        TLegend* legenAsymHighP = GetAndSetLegend2(0.15, 0.8-(2*0.9*textsizeLabelsUp), 0.4, 0.8, textSizeLabels);
        legenAsymHighP->AddEntry(histoAsymPDataHighP,"Data","p");
        legenAsymHighP->AddEntry(histoAsymPTrueMCHighP,"MC","p");
        legenAsymHighP->Draw();
    histoDummyAsymDataMCHighP->Draw("same,axis");
    canvasSingleAsymDataMCHighP->Print(Form("%s/PhotonAsymDataMCHighPSingle%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));



    //_______________________ Ploting weights vs pT in each R bin __________________________

    TCanvas *canvasMBWeightEachR          = new TCanvas("canvasMBWeighEachR","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasMBWeightEachR, 0, 0, 0, 0);
    canvasMBWeightEachR->cd();

    TPad * padMBWeightEachR               = new TPad("padMBWeightEachR","",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    DrawGammaPadSettings( padMBWeightEachR, 0, 0, 0, 0);
    padMBWeightEachR->Divide(4,3,0.0,0.0);
    padMBWeightEachR->Draw();
    Int_t place  = 0;
    TH2F *histoDummyWeightEachR =  new TH2F("histoDummyWeightEachR","histoDummyWeightEachR",1000,0.,1.,1000,0.9,1.3);
    SetStyleHistoTH2ForGraphs(histoDummyWeightEachR, "#it{p}_{T}^{Min} (GeV/c)","w_i", 0.05,0.05, 0.05,0.05);
    for(Int_t i=0; i < nBinsR; i++){
      place  = place + 1;
      padMBWeightEachR->cd(place);
      padMBWeightEachR->cd(place)->SetTopMargin(0.12);
      padMBWeightEachR->cd(place)->SetBottomMargin(0.15);
      padMBWeightEachR->cd(place)->SetLeftMargin(0.1);
      padMBWeightEachR->cd(place)->SetRightMargin(0.1);
      histoDummyWeightEachR->DrawCopy();
      DrawGammaLines(0.,1.5,1., 1.,1.,kGray,1);

      DrawGammaSetMarker(histoWeightsEachRPtMin[i], 20, markerSize, colorData, colorData);
      histoWeightsEachRPtMin[i]->Draw("same");
      DrawGammaSetMarker(histoWeightsEachRPtMinSecSub[i], 25, markerSize, colorData, colorData);
      histoWeightsEachRPtMinSecSub[i]->Draw("same");
      DrawGammaSetMarker(histoWeightsEachRPtMinSecSubUsingCocktail[i], 26, markerSize, colorData, colorData);
      histoWeightsEachRPtMinSecSubUsingCocktail[i]->Draw("same");

    }
    canvasMBWeightEachR->Print(Form("%s/MBWeightVSPtMinEachR%s_%s.%s",outputDirectory.Data(),optionPeriod.Data(),fCutSelectionRead.Data(),suffix.Data()));


    //_______________________ creation TProfile for MB weights Secondary subtracted___________________________________
    TProfile* fProfileContainingMaterialBudgetWeightsFullPt = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBinsFullPt","profileContainingMaterialBudgetWeights_manyRadialBinsFullPt",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsFullPt->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsFullPt->GetYaxis()->SetTitle("weight_per_gamma");

    cout << "MB weights for full pT range:" << endl;
    for(Int_t i=0; i<nBinsR; i++){
      cout << arrayRBins[i] << " < R < " << arrayRBins[i+1] << " cm, weight = " << histoDataMCRatioRScaledToGasSecSub->GetBinContent(i+1) << " +- " << histoDataMCRatioRScaledToGasSecSub->GetBinError(i+1)<< endl;
      fProfileContainingMaterialBudgetWeightsFullPt->Fill(arrayRBins[i],histoDataMCRatioRScaledToGasSecSub->GetBinContent(i+1));
    }
    cout<< endl;

    TProfile* fProfileContainingMaterialBudgetWeightsPtBin1 = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins01","profileContainingMaterialBudgetWeights_manyRadialBins01",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsPtBin1->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsPtBin1->GetYaxis()->SetTitle("weight_per_gamma, #it{p}_{T} > 0.1 GeV/c");

    cout << "MB weights for " << projPtBins[0] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c" << endl;
    for(Int_t i=0; i<nBinsR; i++){
        cout << arrayRBins[i] << " < R < " << arrayRBins[i+1] << " cm, weight = " << histoDataMCRatioRScaledToGasSecSubPtBin1->GetBinContent(i+1) << " +- " << histoDataMCRatioRScaledToGasSecSubPtBin1->GetBinError(i+1)<< endl;
        fProfileContainingMaterialBudgetWeightsPtBin1->Fill(arrayRBins[i], histoDataMCRatioRScaledToGasSecSubPtBin1->GetBinContent(i+1));
    }
    cout<< endl;

    TProfile* fProfileContainingMaterialBudgetWeightsPtBin2 = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins02","profileContainingMaterialBudgetWeights_manyRadialBins02",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsPtBin2->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsPtBin2->GetYaxis()->SetTitle("weight_per_gamma, #it{p}_{T} > 0.3 GeV/c");

    cout << "MB weights for " << projPtBins[1] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c" << endl;
    for(Int_t i=0; i<nBinsR; i++){
        cout << arrayRBins[i] << " < R < " << arrayRBins[i+1] << " cm, weight = " << histoDataMCRatioRScaledToGasSecSubPtBin2->GetBinContent(i+1) << " +- " << histoDataMCRatioRScaledToGasSecSubPtBin2->GetBinError(i+1)<< endl;
        fProfileContainingMaterialBudgetWeightsPtBin2->Fill(arrayRBins[i], histoDataMCRatioRScaledToGasSecSubPtBin2->GetBinContent(i+1));
    }
    cout<< endl;

    TProfile* fProfileContainingMaterialBudgetWeightsPtBin3 = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins03","profileContainingMaterialBudgetWeights_manyRadialBins03",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsPtBin3->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsPtBin3->GetYaxis()->SetTitle("weight_per_gamma, #it{p}_{T} > 0.4 GeV/c");

    cout << "MB weights for " << projPtBins[2] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c (-> DEFAULT)" << endl;
    for(Int_t i=0; i<nBinsR; i++){
        cout << arrayRBins[i] << " < R < " << arrayRBins[i+1] << " cm, weight = " << histoDataMCRatioRScaledToGasSecSubPtBin3->GetBinContent(i+1) << " +- " << histoDataMCRatioRScaledToGasSecSubPtBin3->GetBinError(i+1)<< endl;
        fProfileContainingMaterialBudgetWeightsPtBin3->Fill(arrayRBins[i], histoDataMCRatioRScaledToGasSecSubPtBin3->GetBinContent(i+1));
    }
    cout<< endl;

    //_______________________ creation TProfile for MB weights with secondaries ___________________________________
    TProfile* fProfileContainingMaterialBudgetWeightsWithSecFullPt = new TProfile("profileContainingMaterialBudgetWeightsWithSec_manyRadialBinsFullPt","profileContainingMaterialBudgetWeightsWithSec_manyRadialBinsFullPt",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsWithSecFullPt->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsWithSecFullPt->GetYaxis()->SetTitle("weight_per_gamma");

    cout << "MB weights for full pT range:" << endl;
    for(Int_t i=0; i<nBinsR; i++){
      cout << arrayRBins[i] << " < R < " << arrayRBins[i+1] << " cm, weight = " << histoDataMCRatioRScaledToGas->GetBinContent(i+1) << " +- " << histoDataMCRatioRScaledToGas->GetBinError(i+1)<< endl;
      fProfileContainingMaterialBudgetWeightsWithSecFullPt->Fill(arrayRBins[i],histoDataMCRatioRScaledToGas->GetBinContent(i+1));
    }
    cout<< endl;

    TProfile* fProfileContainingMaterialBudgetWeightsWithSecPtBin1 = new TProfile("profileContainingMaterialBudgetWeightsWithSec_manyRadialBins01","profileContainingMaterialBudgetWeightsWithSec_manyRadialBins01",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsWithSecPtBin1->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsWithSecPtBin1->GetYaxis()->SetTitle("weight_per_gamma, #it{p}_{T} > 0.1 GeV/c");

    cout << "MB weights for " << projPtBins[0] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c" << endl;
    for(Int_t i=0; i<nBinsR; i++){
        cout << arrayRBins[i] << " < R < " << arrayRBins[i+1] << " cm, weight = " << histoDataMCRatioRinPtBinScaledToGasPtBin1->GetBinContent(i+1) << " +- " << histoDataMCRatioRinPtBinScaledToGasPtBin1->GetBinError(i+1)<< endl;
        fProfileContainingMaterialBudgetWeightsWithSecPtBin1->Fill(arrayRBins[i], histoDataMCRatioRinPtBinScaledToGasPtBin1->GetBinContent(i+1));
    }
    cout<< endl;

    TProfile* fProfileContainingMaterialBudgetWeightsWithSecPtBin2 = new TProfile("profileContainingMaterialBudgetWeightsWithSec_manyRadialBins02","profileContainingMaterialBudgetWeightsWithSec_manyRadialBins02",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsWithSecPtBin2->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsWithSecPtBin2->GetYaxis()->SetTitle("weight_per_gamma, #it{p}_{T} > 0.3 GeV/c");

    cout << "MB weights for " << projPtBins[1] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c" << endl;
    for(Int_t i=0; i<nBinsR; i++){
        cout << arrayRBins[i] << " < R < " << arrayRBins[i+1] << " cm, weight = " << histoDataMCRatioRinPtBinScaledToGasPtBin2->GetBinContent(i+1) << " +- " << histoDataMCRatioRinPtBinScaledToGasPtBin2->GetBinError(i+1)<< endl;
        fProfileContainingMaterialBudgetWeightsWithSecPtBin2->Fill(arrayRBins[i], histoDataMCRatioRinPtBinScaledToGasPtBin2->GetBinContent(i+1));
    }
    cout<< endl;

    TProfile* fProfileContainingMaterialBudgetWeightsWithSecPtBin3 = new TProfile("profileContainingMaterialBudgetWeightsWithSec_manyRadialBins03","profileContainingMaterialBudgetWeightsWithSec_manyRadialBins03",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsWithSecPtBin3->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsWithSecPtBin3->GetYaxis()->SetTitle("weight_per_gamma, #it{p}_{T} > 0.4 GeV/c");

    cout << "MB weights for " << projPtBins[2] << " GeV/c to " << histoRPtData->GetXaxis()->GetBinUpEdge(histoRPtData->GetNbinsX()) << " GeV/c (-> DEFAULT)" << endl;
    for(Int_t i=0; i<nBinsR; i++){
        cout << arrayRBins[i] << " < R < " << arrayRBins[i+1] << " cm, weight = " << histoDataMCRatioRinPtBinScaledToGasPtBin3->GetBinContent(i+1) << " +- " << histoDataMCRatioRinPtBinScaledToGasPtBin3->GetBinError(i+1)<< endl;
        fProfileContainingMaterialBudgetWeightsWithSecPtBin3->Fill(arrayRBins[i], histoDataMCRatioRinPtBinScaledToGasPtBin3->GetBinContent(i+1));
    }
    cout<< endl;



    //________________________ pT > 0.4 taken as default secondary subtracted: __________________________
    TProfile* fProfileContainingMaterialBudgetWeights = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins","profileContainingMaterialBudgetWeights_manyRadialBins",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeights->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeights->GetYaxis()->SetTitle("weight_per_gamma, #it{p}_{T} > 0.4 GeV/c");
    for(Int_t i=0; i<nBinsR; i++){
      fProfileContainingMaterialBudgetWeights->Fill(arrayRBins[i], histoDataMCRatioRScaledToGasSecSubPtBin3->GetBinContent(i+1));
    }
    cout<< endl;


    if(createDummy.CompareTo("onlyGasCorr")==0){

        // file with weights all at 1. except for the gas (last 2 bins)
        // -> the gas correction is set at the beginning of the macro
        TString dummyfilename = Form("MCInputFileMaterialBudgetWeightsDummy%s_%s.root",createDummy.Data(),optionEnergy.Data());
        cout << "Creating dummy file with only gas corr: " << dummyfilename.Data();
        TFile outFile(dummyfilename,"RECREATE");
        TProfile* fProfileMaterialBudgetWeightsDummy = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins", "profileContainingMaterialBudgetWeights_manyRadialBins", nBinsR, arrayRBins);
        fProfileMaterialBudgetWeightsDummy->GetXaxis()->SetTitle("R (cm)");
        fProfileMaterialBudgetWeightsDummy->GetYaxis()->SetTitle("weight_per_gamma");

        for(Int_t i=0; i<nBinsR; i++){
            if(arrayRBins[i] < 95.) fProfileMaterialBudgetWeightsDummy->Fill(arrayRBins[i], 1.);
            else if(arrayRBins[i] >= 95.) fProfileMaterialBudgetWeightsDummy->Fill(arrayRBins[i], mcGasCorrectionFactor);
        }

        fProfileMaterialBudgetWeightsDummy->Write();
        outFile.Close();

    } else if(createDummy.CompareTo("plusGasCorr")==0){

        // file with reference weights (pp 13 TeV for example) * gas correction (only in last 2 bins)
        // -> this file needs to be produced running the reference weights with the option on
        // -> the gas correction is set at the beginning of the macro
        TString filename = Form("MCInputFileMaterialBudgetWeights%s_%s.root",createDummy.Data(),optionEnergy.Data());
        cout << "Creating file with MB weights AND gas corr: " << filename.Data();
        TFile outFile(filename,"RECREATE");
        TProfile* fProfileMaterialBudgetWeightsWithGasCorr = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins", "profileContainingMaterialBudgetWeights_manyRadialBins", nBinsR, arrayRBins);
        fProfileMaterialBudgetWeightsWithGasCorr->GetXaxis()->SetTitle("R (cm)");
        fProfileMaterialBudgetWeightsWithGasCorr->GetYaxis()->SetTitle("weight_per_gamma");
	//Am: Take secondary subtracted weights
        for(Int_t i=0; i<nBinsR; i++){
            if(arrayRBins[i] < 95.) fProfileMaterialBudgetWeightsWithGasCorr->Fill(arrayRBins[i], histoDataMCRatioRScaledToGasSecSubPtBin3->GetBinContent(i+1));
            else if(arrayRBins[i] >= 95.) fProfileMaterialBudgetWeightsWithGasCorr->Fill(arrayRBins[i], mcGasCorrectionFactor*histoDataMCRatioRScaledToGasSecSubPtBin3->GetBinContent(i+1));
        }

        fProfileMaterialBudgetWeightsWithGasCorr->Write();
        outFile.Close();

    } else {

        // file with MB weights
        TString filenameMBweights = Form("MCInputFileMaterialBudgetWeights%s_%s.root",optionPeriod.Data(),fCutSelectionRead.Data());
        cout << "Creating file with MB weights: " << filenameMBweights.Data() << endl;
        TFile outFile(filenameMBweights,"RECREATE");

        fProfileContainingMaterialBudgetWeightsFullPt->Write();
        fProfileContainingMaterialBudgetWeights->Write();
        fProfileContainingMaterialBudgetWeightsPtBin1->Write();
        fProfileContainingMaterialBudgetWeightsPtBin2->Write();
        fProfileContainingMaterialBudgetWeightsPtBin3->Write();
        histoDataMCRatioRinPtBinScaledToGasPtBin3->Write();
        histoDataMCRatioRScaledToGasSecSubPtBin3->Write();
        histoRData->Write("Data");

        outFile.Close();

    }

    cout << "Creating file with additional MB histograms: " << Form("AdditionalMBHistos%s_%s.root",optionPeriod.Data(),fCutSelectionRead.Data()) << endl;
    TFile outAddHistoFile(Form("AdditionalMBHistos%s_%s.root",optionPeriod.Data(),fCutSelectionRead.Data()),"RECREATE");

    histoRData->Write("Data");
    histoRDataRebin->Write("DataRebin");
    histoRMC->Write("MC");
    histoRMCRebin->Write("MCRebin");
    histoRDataScaledToGas->Write("DataScaledToGas");
    histoRMCScaledToGas->Write("MCScaledToGas");
    histoDataMCRatioRScaledToGas ->Write();
    histoDataMCRatioRinPtBinScaledToGasPtBin1->Write();
    histoDataMCRatioRinPtBinScaledToGasPtBin2->Write();
    histoDataMCRatioRinPtBinScaledToGasPtBin3->Write();
    histoDataMCRatioRScaledToGasSecSub-> Write();
    histoDataMCRatioRScaledToGasSecSubPtBin1->Write();
    histoDataMCRatioRScaledToGasSecSubPtBin2->Write();
    histoDataMCRatioRScaledToGasSecSubPtBin3->Write();

    histoPurityR->Write();
    histoPurityPrimR->Write();
    histoPuritySecR->Write();
    histoPurityPrimRPtBin3->Write();
    histoPuritySecRPtBin3->Write();
    histoPurityPt5cm->Write();
    histoPurityPrimPt5cm->Write();
    histoPuritySecPt5cm->Write();
    histoEffiR->Write();
    histoEffiPt->Write();
    histoAsymPDataLowP->Write();
    histoAsymPDataHighP->Write();
    histoAsymPTrueMCLowP->Write();
    histoAsymPTrueMCHighP->Write();
    for(Int_t i=0; i < nBinsPtFine; i++){
      //histoRinPtBinDataScaledToGasRebinFine[i]->Write();
      //histoRinPtBinMCScaledToGasRebinFine[i]->Write();
      histoDataMCRatioRinPtBinScaledToGasFine[i]->Write();
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoWeightsEachRPtMin[i]->Write();
    }

    for(Int_t j=0; j < nBinsR; j++){
      histoWeightsEachRPtMinSecSub[j]->Write();
    }
    for(Int_t j=0; j < nBinsR; j++){
      histoWeightsEachRPtMinSecSubUsingCocktail[j]->Write();
    }

    histoSecConvGammaPtAllSources ->Write();
    histoSecConvGammaPtFromK0S ->Write();
    histoPtTrueSecFromK0S  ->Write();
    histoPtTrueSecAllsources  ->Write();
    histoPtTrueSecRest->Write();

    for(Int_t i=0; i < nBinsR; i++){
      histoPtEachRBinData[i] ->Write();
    }

    for(Int_t i=0; i < nBinsR; i++){
      histoPtEachRBinMC[i] ->Write();
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoPtTrueMCEachRBin[i]->Write();
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoPtTruePrimMCEachRBin[i] ->Write();
    }

    for(Int_t i=0; i < nBinsR; i++){
      histoPtTrueSecEachRBinFromK0S[i]->Write();
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoPurityPtEachRBin[i]->Write();
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoPurityPrimPtEachRBin[i]->Write();
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoFracSecPtEachRBin[i]->Write();
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoEfficiencySecEachRBinAllSources[i] ->Write();
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoPtEachRBinDataSecYieldFromSecFrac[i]->Write(); 
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoPtEachRBinMCSecYieldFromSecFrac[i]->Write(); 
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoPtEachRBinDataSecSubtracted[i] ->Write(); 
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoPtEachRBinDataSecSubtractedUsingCocktail[i] ->Write(); 
    }
    for(Int_t i=0; i < nBinsR; i++){
      histoPtEachRBinMCSecSubtracted[i] ->Write(); 
    }

    for(Int_t k=0; k < 4; k++){
      for(Int_t i=0; i < nBinsR; i++){
	histoFracSecPtEachRBinBySource[k][i]->Write();
      }
    }
    for(Int_t k=0; k < 4; k++){
      for(Int_t i=0; i < nBinsR; i++){
	histoMCTrueRecSecGammaPtEachRBinBySource[k][i]->Write();
      }
    }

    for(Int_t k=0; k < 4; k++){
      for(Int_t i=0; i < nBinsR; i++){
	histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[k][i]->Write();
      }
    }
    for(Int_t k=0; k < 4; k++){
      for(Int_t i=0; i < nBinsR; i++){
	histoPtEachRBinDataSecYieldFromSecFracBySourceRaw[k][i]->Write();
      }
    }

    for(Int_t k=0; k < 3; k++){
      for(Int_t i=0; i < nBinsR; i++){
	histoTrueRecEffSecEachRBinBySource[k][i]->Write();
      }
    }
   for(Int_t k=0; k < 3; k++){
      for(Int_t i=0; i < nBinsR; i++){
	histoTrueEffSecEachRBinBySource[k][i]->Write();
      }
    }
   histoPtSumRMCSecSubtracted->Write();
   histoPtSumRBinDataSecSubtractedUsingCocktail->Write();
   histoPtSumRMCSecSubtractedRebin->Write();
   histoPtSumRBinDataSecSubtractedUsingCocktailRebin->Write();
   histoPtDataRebin->Write();
   histoPtDataToMCSecSubtractedRebin->Write();
   outAddHistoFile.Close();


}



//******************************************************************************
// Function taken from ExtractGammaSignalV2.C
// Load secondary gamma histos from cocktail file
// - put them in proper scaling
// - rebin them according to current gamma binning
//******************************************************************************
Bool_t LoadSecondariesFromCocktailFile(TString cutSelection, TString optionEnergy){

    // search for cocktail file
    // For the time being the rapidity range is fixed to 0.8, not passed as parameter
    TString nameCocktailFile                                    = Form("%s/%s/SecondaryGamma_0.80_%s.root",cutSelection.Data(),optionEnergy.Data(),cutSelection.Data());
    fFileCocktailInput                                          = new TFile(nameCocktailFile.Data());
    if (fFileCocktailInput->IsZombie()) fFileCocktailInput      = NULL;

    // break if no cocktail file was found
    if (!fFileCocktailInput) {
      cout << "Cocktail file: " << nameCocktailFile.Data() << " not found!" << endl;
      return kFALSE;
    }

    // get secondary spectra from cocktail file
    cout << "Found cocktail file: " << nameCocktailFile.Data() << " -> will add cocktail histos to output" << endl;
    for (Int_t k = 0; k < 3; k++){
      fHistoSecGammaCocktailFromXPt[k]              = (TH1F*)fFileCocktailInput->Get(Form("Gamma_From_X_From_%s_Pt_OrBin", fSecondaries[k].Data()));
      if(!fHistoSecGammaCocktailFromXPt[k]){
	cout << "Gamma_From_X_From_" << fSecondaries[k].Data() << "_Pt_OrBin not found in cocktail file! Cocktail will not be used." << endl;
	return kFALSE;
      }
      fHistoSecGammaCocktailFromXPt[k]->Sumw2();
      fHistoSecGammaCocktailFromXPtOrBin[k]         = (TH1F*)fHistoSecGammaCocktailFromXPt[k]->Clone(Form("CocktailSecGammaFromXFrom%s_PtOrBin", fSecondaries[k].Data()));
      //         RebinSpectrum(fHistoSecGammaCocktailFromXPt[k]);
      //         fHistoSecGammaCocktailFromXPtOrBin[k]->SetBins(nBins,xMin,xMax);
    }
    
    // all spectra found
    return kTRUE;
 }
