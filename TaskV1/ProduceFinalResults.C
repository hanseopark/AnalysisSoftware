/*******************************************************************************************************************
 ******     provided by Gamma Conversion Group, PWG4,                                                          *****
 ******     Ana Marin, marin@physi.uni-heidelberg.de                                                           *****
 ******     Friederike Bock, friederike.bock@cern.ch                                                           *****
 *******************************************************************************************************************/

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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "ProduceFinalResults.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //    TString name;
};

void  ProduceFinalResults( const char *fileNamePi0 = "myOutput",
                           const char *fileNameEta = "", 
                           TString cutSelection = "", 
                           const char *suffix = "gif", 
                           TString isMC= "", 
                           TString makeBinShiftWithFunction = "", 
                           TString useSameBinningPi0Eta ="", 
                           TString optionEnergy = "",
                           TString conferencePlots ="", 
                           TString multFlag= "",
                           TString thesisPlots="", 
                           TString optDalitz = "", 
                           TString optNoBinShift="kFALSE", 
                           Bool_t pileUpApplied=kTRUE, 
                           Int_t mode = 9){    
    
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();
    SetPlotStyle();
    
    Bool_t kDalitz = kFALSE;
    if (optDalitz.CompareTo("Dalitz")==0){
        kDalitz = kTRUE;
    }
    
    date = ReturnDateString();
    TString fileNameSysErrPi0 ="SystematicErrorsInput/SystematicErrorAveraged_Pi0_7TeV_24_Apr_2012.dat"; // default
    collisionSystem= ReturnFullCollisionsSystem(optionEnergy);

    // original systematic errors were from 18_May_2011 inversion of material error caused changing
    if(optionEnergy.CompareTo("7TeV") == 0){
        minPtForFits=0.3;
        minPtForFitsEta=0.4;
        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fileNameSysErrPi0 = "SystematicErrorsInput/SystematicErrorAveraged_Pi0_7TeV_24_Apr_2012.dat"; //SystematicError_Pi0_7TeV_12_Mar_2012.dat
        } else {
            fileNameSysErrPi0 = "SystematicErrorsInput/SystematicErrorAveraged_Pi0EtaBinning_7TeV_24_Apr_2012.dat"; // SystematicError_Pi0EtaBinning_7TeV_12_Mar_2012.dat
            minPtForFits=0.4;
        }    
        fileNameSysErrEta = "SystematicErrorsInput/SystematicErrorAveraged_Eta_7TeV_24_Apr_2012.dat";//SystematicError_Eta_7TeV_12_Mar_2012.dat
        cout << "You have choosen 7TeV" << endl;
    } else if( optionEnergy.CompareTo("8TeV") == 0) {
        minPtForFits=0.4;
        minPtForFitsEta=0.4;

        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0_8TeV_16_May_2015.dat";
        } else {
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_8TeV_16_May_2015.dat"; 
            minPtForFits=0.4;
        }    
        fileNameSysErrEta = "SystematicErrorsNew/SystematicErrorAveraged_Eta_8TeV_16_May_2015.dat"; 
        cout << "You have choosen 8TeV" << endl;
    } else if( optionEnergy.CompareTo("13TeV") == 0) {
        minPtForFits=0.4;
        minPtForFitsEta=0.4;

        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0_2.76TeV_5_Aug_2013.dat";
        } else {
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_2.76TeV_5_Aug_2013.dat"; 
            minPtForFits=0.4;
        }    
        fileNameSysErrEta = "SystematicErrorsNew/SystematicErrorAveraged_Eta_2.76TeV_5_Aug_2013.dat"; 
        cout << "You have choosen 13TeV. Note that you will use 2.76TeV systematic errors." << endl;
        
    } else if( optionEnergy.CompareTo("2.76TeV") == 0) {
        minPtForFits=0.4;
        minPtForFitsEta=0.6;
        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0_2.76TeV_5_Aug_2013.dat";
        } else {
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_2.76TeV_5_Aug_2013.dat"; 
            minPtForFits=0.6;
        }
        fileNameSysErrEta = "SystematicErrorsNew/SystematicErrorAveraged_Eta_2.76TeV_5_Aug_2013.dat"; 
        if (mode == 2){
            minPtForFits=0.6;
            minPtForFitsEta=1;

//             SystematicErrorsNew/SystematicErrorAveraged_Pi0_2.76TeVEG1_2015_10_19.dat
//             SystematicErrorsNew/SystematicErrorAveraged_Pi0_2.76TeVEG2_2015_10_19.dat
//             SystematicErrorsNew/SystematicErrorAveraged_Pi0_2.76TeVEMC1_2015_10_19.dat
//             SystematicErrorsNew/SystematicErrorAveraged_Pi0_2.76TeVEMC7_2015_10_19.dat
//             SystematicErrorsNew/SystematicErrorAveraged_Pi0_2.76TeVINT7_2015_10_19.dat
            if (useSameBinningPi0Eta.CompareTo("")==0){         
                fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0_2.76TeV_2015_10_19.dat";
            } else {
                fileNameSysErrPi0 = "SystematicErrorsCalculatedConvCalo/SystematicErrorAveraged_Pi0EtaBinning_2.76TeV_2015_05_18.dat"; 
                minPtForFits=1;
            }    
            fileNameSysErrEta = "SystematicErrorsCalculatedConvCalo/SystematicErrorAveraged_Eta_2.76TeV_2015_05_18.dat"; 
        }    
        
        
    } else if( optionEnergy.CompareTo("900GeV") == 0) {
        minPtForFits=0.4;
        minPtForFitsEta=0.9;
        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fileNameSysErrPi0 = "SystematicErrorsInput/SystematicErrorAveraged_Pi0_900GeV_24_Apr_2012.dat"; //SystematicError_Pi0_900GeV_12_Mar_2012.dat
        } else {
            fileNameSysErrPi0 = "SystematicErrorsInput/SystematicErrorAveraged_Pi0EtaBinning_900GeV_24_Apr_2012.dat"; //SystematicError_Pi0EtaBinning_900GeV_12_Mar_2012.dat
            minPtForFits=0.9;
        }
        fileNameSysErrEta = "SystematicErrorsInput/SystematicErrorAveraged_Eta_900GeV_24_Apr_2012.dat"; //SystematicError_Eta_900GeV_12_Mar_2012.dat
    } else if( optionEnergy.CompareTo("PbPb_2.76TeV") == 0 || optionEnergy.CompareTo("pPb_5.023TeV") == 0  ) {
        cout << "This macro will not work for pPb or PbPb" << endl;
      return;
    } else {
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    
    Int_t fExampleBinPi0     = 7;
    Int_t fExampleBinEta     = 6;
    if(optionEnergy.CompareTo("7TeV") == 0){
        fExampleBinPi0     = 7;
        fExampleBinEta     = 6;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("8TeV") == 0) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 6;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("2.76TeV") == 0) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 4;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("900GeV") == 0) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 2;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("pPb_5.023TeV") == 0 ) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 4;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("PbPb_2.76TeV") == 0 ) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 6;
    }
    
    TString outputDir = Form("%s/%s/%s/ProduceFinalResults%s",cutSelection.Data(),optionEnergy.Data(),suffix, optDalitz.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    
    if(conferencePlots.CompareTo("conference") == 0){// means we want to plot values for the pi0
        conference = kTRUE;
        //-AM
        // minPtForFits = 0.4;
    }    

    if(thesisPlots.CompareTo("thesis") == 0){// means we want to plot values for the pi0
        thesis = kTRUE;
    }    
    
    TString fEventCutSelection="";
    TString fGammaCutSelection="";
    TString fClusterCutSelection="";
    TString fElectronCutSelection="";
    TString fMesonCutSelection="";
    
    if (mode == 9){
        if (kDalitz){
            ReturnSeparatedCutNumber(cutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection,kTRUE);
        } else {
            ReturnSeparatedCutNumber(cutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
        }
        
        fEventCutSelection = fGammaCutSelection(0,7);
        fGammaCutSelection = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
        cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
    } else {
        ReturnSeparatedCutNumberAdvanced(cutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
    }
   
    deltaEta =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
            
    //declaration for printing logo 
    
    if (isMC.CompareTo("kTRUE")==0){ 
        prefix2 = "MC";
        pictDrawingOptions[1] = kTRUE;
    } else {    
        prefix2 = "data";
        pictDrawingOptions[1] = kFALSE;
    }
    
    if (makeBinShiftWithFunction.CompareTo("Hagedorn")==0){ 
        kHag = kTRUE;
        kLevy = kFALSE;
    }else { 
        if(makeBinShiftWithFunction.CompareTo("Levy")==0){
            kHag = kFALSE;
            kLevy = kTRUE;
        }else {
            cout << "No known fit for Binshift defined, Levy will be taken!" <<endl;
            kHag = kFALSE;
            kLevy = kTRUE;
        }
    }
    
    mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    
    TString nameCorrectedYield;
    TString nameCorrectedYieldEta;
    TString nameEfficiency;
    TString nameEfficiencyEta;
    if ( mode == 0 || mode == 9 || mode == 1  ){
        nameCorrectedYield = "CorrectedYieldTrueEff";
        nameEfficiency = "TrueMesonEffiPt";
        nameCorrectedYieldEta = "CorrectedYieldTrueEff";
        nameEfficiencyEta = "TrueMesonEffiPt";
    } else if (mode == 2){
        nameCorrectedYield = "CorrectedYieldTrueEff";
        nameEfficiency = "TrueMesonEffiPt";
        nameCorrectedYieldEta = "CorrectedYieldTrueEff";
        nameEfficiencyEta = "TrueMesonEffiPt";
    } else {
        nameCorrectedYield = "CorrectedYieldNormEff";
        nameEfficiency = "MesonEffiPt";
        nameCorrectedYieldEta = "CorrectedYieldNormEff";
        nameEfficiencyEta = "MesonEffiPt";
    }    
    
    cout << "reading pi0 file" << endl;
    // File definitions
    fileNamePi0ch                           = fileNamePi0;
    filePi0                                 = new TFile(fileNamePi0ch);
    histoCorrectedYieldPi0                  = (TH1D*)filePi0->Get(nameCorrectedYield.Data());
    histoUncorrectedYieldPi0                = (TH1D*)filePi0->Get("histoYieldMeson");
    histoFWHMMesonPi0                       = (TH1D*)filePi0->Get("histoFWHMMeson");
    histoMassMesonPi0                       = (TH1D*)filePi0->Get("histoMassMeson");
    histoAccPi0                             = (TH1D*)filePi0->Get("fMCMesonAccepPt");
    histoTrueEffPtPi0                       = (TH1D*)filePi0->Get(nameEfficiency.Data()); 
    histoTrueFWHMMesonPi0                   = (TH1D*)filePi0->Get("histoTrueFWHMMeson");
    histoTrueMassMesonPi0                   = (TH1D*)filePi0->Get("histoTrueMassMeson");
    histoEventQualtityPi0                   = (TH1F*)filePi0->Get("NEvents");
    histoMCInputPi0                         = (TH1D*)filePi0->Get("MCYield_Meson_oldBin");
    TH1D* histoMCInputPi0WOWeight           = (TH1D*)filePi0->Get("MCYield_Meson_oldBinWOWeights");
    TH1D* histoMCInputPi0Weights            = (TH1D*)filePi0->Get("WeightsMeson");
    TH1D* histoMCInputPi0AddedSig           = (TH1D*)filePi0->Get("MCYield_Meson_oldBin_AddedSig");
    TH1D* histoMCInputPi0WOWeightAddedSig   = (TH1D*)filePi0->Get("MCYield_Meson_oldBinWOWeights_AddedSig");
    TH1D* histoMCInputPi0WeightsAddedSig    = (TH1D*)filePi0->Get("WeightsMeson_AddedSig");
    TH1D* histoSignalPlusBGPi0              = (TH1D*)filePi0->Get(Form("InvMassSigPlusBG_PtBin%02d",fExampleBinPi0));
    TH1D* histoBGPi0                        = (TH1D*)filePi0->Get(Form("InvMassBG_PtBin%02d",fExampleBinPi0));
    TH1D* histoSignalPi0                    = (TH1D*)filePi0->Get(Form("InvMassSig_PtBin%02d",fExampleBinPi0));
    TF1* fitSignalPi0                       = (TF1*)filePi0->Get(Form("FitInvMassSig_PtBin%02d",fExampleBinPi0));
    
    histoMassMesonPi0MinusExp= CalculateMassMinusExpectedMass(histoMassMesonPi0,mesonMassExpectPi0);
    histoTrueMassMesonPi0MinusExp = CalculateMassMinusExpectedMass(histoTrueMassMesonPi0,mesonMassExpectPi0);
    histoFWHMMesonPi0MeV = (TH1D*) histoFWHMMesonPi0->Clone();
    histoFWHMMesonPi0MeV->Scale(1000.);
    histoTrueFWHMMesonPi0MeV = (TH1D*) histoTrueFWHMMesonPi0->Clone();
    histoTrueFWHMMesonPi0MeV->Scale(1000.);
    
    if (pileUpApplied && optionEnergy.CompareTo("7TeV") == 0){
        cout << "correcting pi0 for PileUp with factor: " << fPileUpCorrectionConv7TeV << endl;
        histoCorrectedYieldPi0->Scale(fPileUpCorrectionConv7TeV);
        histoUncorrectedYieldPi0->Scale(fPileUpCorrectionConv7TeV);
    }
    
//     widthMeV = histoMesonSignalFullPtInvMass->GetBinWidth(10)*1000;
    nEvt =  GetNEvents(histoEventQualtityPi0);
    pictDrawingCoordinates[8] = nEvt;
    
    cout << "reading eta file" << endl;
    fileEta                                 = new TFile(fileNameEta);
    histoCorrectedYieldEta                     = (TH1D*)fileEta->Get(nameCorrectedYieldEta.Data()); 
    histoUnCorrectedYieldEta                 = (TH1D*)fileEta->Get("histoYieldMeson");
    histoFWHMMesonEta                        = (TH1D*)fileEta->Get("histoFWHMMeson");
    histoMassMesonEta                        = (TH1D*)fileEta->Get("histoMassMeson");
    histoAccEta                                = (TH1D*)fileEta->Get("fMCMesonAccepPt");
    histoTrueEffPtEta                         = (TH1D*)fileEta->Get(nameEfficiencyEta.Data()); //not yet correct MesonEffiPt
    histoTrueFWHMMesonEta                     = (TH1D*)fileEta->Get("histoTrueFWHMMeson");
    histoTrueMassMesonEta                     = (TH1D*)fileEta->Get("histoTrueMassMeson"); 
    histoMCInputEta                         = (TH1D*)fileEta->Get("MCYield_Meson_oldBin");
    TH1D* histoMCInputEtaWOWeight             = (TH1D*)fileEta->Get("MCYield_Meson_oldBinWOWeights");
    TH1D* histoMCInputEtaWeights             = (TH1D*)fileEta->Get("WeightsMeson");
    TH1D* histoMCInputEtaAddedSig             = (TH1D*)fileEta->Get("MCYield_Meson_oldBin_AddedSig");
    TH1D* histoMCInputEtaWOWeightAddedSig     = (TH1D*)fileEta->Get("MCYield_Meson_oldBinWOWeights_AddedSig");
    TH1D* histoMCInputEtaWeightsAddedSig     = (TH1D*)fileEta->Get("WeightsMeson_AddedSig");
    TH1D* histoSignalPlusBGEta                 = (TH1D*)fileEta->Get(Form("InvMassSigPlusBG_PtBin%02d",fExampleBinEta));
    TH1D* histoBGEta                         = (TH1D*)fileEta->Get(Form("InvMassBG_PtBin%02d",fExampleBinEta));
    TH1D* histoSignalEta                    = (TH1D*)fileEta->Get(Form("InvMassSig_PtBin%02d",fExampleBinEta));
    TF1* fitSignalEta                        = (TF1*)fileEta->Get(Form("FitInvMassSig_PtBin%02d",fExampleBinEta));

    histoMassMesonEtaMinusExp= CalculateMassMinusExpectedMass(histoMassMesonEta,mesonMassExpectEta);
    histoTrueMassMesonEtaMinusExp = CalculateMassMinusExpectedMass(histoTrueMassMesonEta,mesonMassExpectEta);
    histoFWHMMesonEtaMeV = (TH1D*) histoFWHMMesonEta->Clone();
    histoFWHMMesonEtaMeV->Scale(1000.);
    histoTrueFWHMMesonEtaMeV = (TH1D*) histoTrueFWHMMesonEta->Clone();
    histoTrueFWHMMesonEtaMeV->Scale(1000.);
    
    if (pileUpApplied && optionEnergy.CompareTo("7TeV") == 0){
        cout << "correcting eta for PileUp with factor: " << fPileUpCorrectionConv7TeV << endl;
        histoCorrectedYieldEta->Scale(fPileUpCorrectionConv7TeV);
        histoUnCorrectedYieldEta->Scale(fPileUpCorrectionConv7TeV);
    }
    
    fileSysErrEta.open(fileNameSysErrEta,ios_base::in);
    cout << fileNameSysErrEta << endl;
    
    while(!fileSysErrEta.eof() && nPointsEta < 100){
        fileSysErrEta >> relSystErrorEtaDown[nPointsEta] >> relSystErrorEtaUp[nPointsEta] >> relSystErrorWOMaterialEtaDown[nPointsEta] >> relSystErrorWOMaterialEtaUp[nPointsEta];
        cout << nPointsEta << "\t"  << relSystErrorEtaDown[nPointsEta] << "\t"  <<relSystErrorEtaUp[nPointsEta] <<  "\t"  << relSystErrorWOMaterialEtaDown[nPointsEta] << "\t"  <<relSystErrorWOMaterialEtaUp[nPointsEta] << endl;
        nPointsEta++;
    }
    fileSysErrEta.close();
    nPointsEta = nPointsEta-1;
    
    fileSysErrPi0.open(fileNameSysErrPi0,ios_base::in);
    cout << fileNameSysErrPi0 << endl;

    while(!fileSysErrPi0.eof() && nPointsPi0 < 100){
        fileSysErrPi0 >> relSystErrorPi0Down[nPointsPi0] >> relSystErrorPi0Up[nPointsPi0]>>    relSystErrorWOMaterialPi0Down[nPointsPi0] >> relSystErrorWOMaterialPi0Up[nPointsPi0];
        cout << nPointsPi0 << "\t"  << relSystErrorPi0Down[nPointsPi0] << "\t"  <<relSystErrorPi0Up[nPointsPi0] << "\t" << relSystErrorWOMaterialPi0Down[nPointsPi0] << "\t"  <<relSystErrorWOMaterialPi0Up[nPointsPi0] << endl;;
        nPointsPi0++;
    }
    fileSysErrPi0.close();
    nPointsPi0 = nPointsPi0-1;
    
    maxPtPi0 = histoCorrectedYieldPi0->GetXaxis()->GetBinUpEdge(histoCorrectedYieldPi0->GetNbinsX());
    maxPtEta = histoCorrectedYieldEta->GetXaxis()->GetBinUpEdge(histoCorrectedYieldEta->GetNbinsX());
    
    TFile* fileMCGenerated =     new TFile("ExternalInput/PCM/MCGeneratedSpectra.root");
    if(optionEnergy.CompareTo("8TeV") == 0){
        histoEtaToPi0Phojet = (TH1D*)fileMCGenerated->Get(Form("EtaToPi0_generatedSpectrum_%s_Phojet","7TeV"));
        histoEtaToPi0Pythia = (TH1D*)fileMCGenerated->Get(Form("EtaToPi0_generatedSpectrum_%s_Pythia","7TeV"));
    
        histoPi0ToChargedPhojet = (TH1D*)fileMCGenerated->Get(Form("Pi0ToCharged_generatedSpectrum_%s_Phojet","7TeV"));
        histoPi0ToChargedPythia = (TH1D*)fileMCGenerated->Get(Form("Pi0ToCharged_generatedSpectrum_%s_Pythia","7TeV"));
    } else if(optionEnergy.CompareTo("13TeV")==0){
        cout << "Caution!!! use 7TeV MC generated EtaToPi0 and Pi0ToCharged spectra" << endl;
        histoEtaToPi0Phojet = (TH1D*)fileMCGenerated->Get(Form("EtaToPi0_generatedSpectrum_7TeV_Phojet"));
        histoEtaToPi0Pythia = (TH1D*)fileMCGenerated->Get(Form("EtaToPi0_generatedSpectrum_7TeV_Pythia"));
        histoPi0ToChargedPhojet = (TH1D*)fileMCGenerated->Get(Form("Pi0ToCharged_generatedSpectrum_7TeV_Phojet"));
        histoPi0ToChargedPythia = (TH1D*)fileMCGenerated->Get(Form("Pi0ToCharged_generatedSpectrum_7TeV_Pythia"));
    } else {
        histoEtaToPi0Phojet = (TH1D*)fileMCGenerated->Get(Form("EtaToPi0_generatedSpectrum_%s_Phojet",optionEnergy.Data()));
        histoEtaToPi0Pythia = (TH1D*)fileMCGenerated->Get(Form("EtaToPi0_generatedSpectrum_%s_Pythia",optionEnergy.Data()));
    
        histoPi0ToChargedPhojet = (TH1D*)fileMCGenerated->Get(Form("Pi0ToCharged_generatedSpectrum_%s_Phojet",optionEnergy.Data()));
        histoPi0ToChargedPythia = (TH1D*)fileMCGenerated->Get(Form("Pi0ToCharged_generatedSpectrum_%s_Pythia",optionEnergy.Data()));
    }
        
    histoNumberOfEvents =  new TH1D("histoNumberOfEvents","histoNumberOfEvents",2,0.,2.);
    histoNumberOfEvents->SetBinContent(1,nEvt);

    Int_t isV0AND           = 0; 
    if (histoEventQualtityPi0->GetNbinsX() > 7){
        if (histoEventQualtityPi0->GetBinContent(9) > 0){
            isV0AND         = 1; 
            histoNumberOfEvents->SetBinContent(2,2);
        }
    }
    if (optionEnergy.CompareTo("8TeV") == 0){
        isV0AND             = 1;
    }    
    xSection                = ReturnCorrectXSection( optionEnergy, isV0AND);
    
    
   //**********************************************************************************
   //******************** FWHM Plot ***************************************************
   //**********************************************************************************
    TCanvas* canvasFWHM = new TCanvas("canvasFWHM","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasFWHM, 0.08, 0.02, 0.04, 0.09);
    
    DrawAutoGammaMesonHistos( histoFWHMMesonPi0, 
                                 "", "p_{T} (GeV/c)", "FWHM/2.36 (GeV/c^{2})", 
                                 kFALSE, 3.,0., kFALSE,
                                 kTRUE, -0.004, 0.020, 
                                 kFALSE, 0., 10.);
    histoFWHMMesonPi0->GetYaxis()->SetNdivisions(510); 
    DrawGammaSetMarker(histoFWHMMesonEta, 20, 0.8, kBlue, kBlue);
    histoFWHMMesonEta->DrawCopy("same,e1,p"); 
    
    DrawGammaSetMarker(histoTrueFWHMMesonEta, 24, 0.8, kRed+2, kRed+2);
    histoTrueFWHMMesonEta->DrawCopy("same,e1,p"); 
    DrawGammaSetMarker(histoFWHMMesonPi0, 22, 0.8, kBlack, kBlack); 
    histoFWHMMesonPi0->DrawCopy("same,e1,p"); 
    
    DrawGammaSetMarker(histoTrueFWHMMesonPi0, 26, 0.8, kRed+2, kRed+2);
    histoTrueFWHMMesonPi0->DrawCopy("same,e1,p"); 
    
    TLegend* legendFWHM = new TLegend(0.3,0.1,0.7,0.2);
    legendFWHM->SetTextSize(0.02);
    legendFWHM->SetFillColor(0);
    legendFWHM->SetLineColor(0);
    if(conference){
        legendFWHM->AddEntry(histoFWHMMesonPi0,Form("#pi^{0}"));
    }else{
        legendFWHM->AddEntry(histoFWHMMesonPi0,Form("#pi^{0}  %s", cutSelection.Data()));
    }
    legendFWHM->AddEntry(histoFWHMMesonEta,"#eta");
    legendFWHM->AddEntry(histoTrueFWHMMesonPi0,"True reconstructed #pi^{0}");
    legendFWHM->AddEntry(histoTrueFWHMMesonEta,"True reconstructed #eta");
    legendFWHM->Draw();
//     if(!thesis)DrawAliceLogoCombined(pictDrawingCoordinatesFWHM[0], pictDrawingCoordinatesFWHM[1], pictDrawingCoordinatesFWHM[2], pictDrawingCoordinatesFWHM[3], pictDrawingCoordinatesFWHM[4], pictDrawingCoordinatesFWHM[5], pictDrawingCoordinatesFWHM[6], pictDrawingCoordinatesFWHM[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,kDalitz);
    
    canvasFWHM->Update();
    
    canvasFWHM->SaveAs(Form("%s/%s_FWHMcombined_%s_%s.%s",outputDir.Data(),prefix2.Data(),useSameBinningPi0Eta.Data(),cutSelection.Data(),suffix));
    delete canvasFWHM;
    
    
    //**********************************************************************************
    //******************** Mass_Pi0 Plot *********************************************
    //**********************************************************************************
    TCanvas* canvasMassPi0 = new TCanvas("canvasMassPi0","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasMassPi0, 0.1, 0.02, 0.04, 0.09);
    
    DrawAutoGammaMesonHistos( histoMassMesonPi0, 
                             "", "p_{T} (GeV/c)", Form("Mass (GeV/c^{2})"), 
                             kFALSE, 3.,0.,  kFALSE,
                             kTRUE, 0.130,0.140, 
                             kFALSE, 0., 10.);
    histoMassMesonPi0->GetYaxis()->SetNdivisions(510); 
    
    DrawGammaSetMarker(histoMassMesonPi0, 22, 0.8, kBlack, kBlack); 
    histoMassMesonPi0->DrawCopy("e1,p"); 
    
    DrawGammaSetMarker(histoTrueMassMesonPi0, 26, 0.8, kRed+2, kRed+2);
    histoTrueMassMesonPi0->DrawCopy("same,e1,p"); 
    
    
    TLegend* legendMassPi0 = new TLegend(0.15,0.1,0.5,0.2);
    legendMassPi0->SetTextSize(0.02);
    legendMassPi0->SetFillColor(0);
    legendMassPi0->AddEntry(histoMassMesonPi0,"Data");
    legendMassPi0->AddEntry(histoTrueMassMesonPi0,"MonteCarlo");
    legendMassPi0->Draw();
//     if(!thesis)DrawAliceLogoPi0Performance(pictDrawingCoordinatesMass[0], pictDrawingCoordinatesMass[1], pictDrawingCoordinatesMass[2], pictDrawingCoordinatesMass[3], pictDrawingCoordinatesMass[4], pictDrawingCoordinatesMass[5], pictDrawingCoordinatesMass[6], pictDrawingCoordinatesMass[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1350,900,date, "MinBias",kDalitz);
    
    canvasMassPi0->Update();
     if(!conference)canvasMassPi0->SaveAs(Form("%s/%s_Mass_Pi0_combined_%s_%s.%s",outputDir.Data(),prefix2.Data(),useSameBinningPi0Eta.Data(),cutSelection.Data(),suffix));
    delete canvasMassPi0;
    
    //**********************************************************************************
    //******************** Mass_Eta Plot *********************************************
    //**********************************************************************************
    TCanvas* canvasMassEta = new TCanvas("canvasMassEta","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasMassEta, 0.1, 0.02, 0.04, 0.09);
    
    DrawAutoGammaMesonHistos( histoMassMesonEta, 
                            "", "p_{T} (GeV/c)", "Mass (GeV/c^{2})", 
                            kFALSE, 3.,0., kFALSE,
                            kTRUE, 0.545,0.560, 
                            kFALSE, 0., 10.);
    histoMassMesonEta->GetYaxis()->SetNdivisions(510); 
    
    DrawGammaSetMarker(histoMassMesonEta, 22, 0.8, kBlack, kBlack); 
    histoMassMesonEta->DrawCopy("e1,p"); 
    
    DrawGammaSetMarker(histoTrueMassMesonEta, 26, 0.8, kRed+2, kRed+2);
    histoTrueMassMesonEta->DrawCopy("same,e1,p"); 
    
    TLegend* legendMassEta = new TLegend(0.15,0.1,0.5,0.2);
    legendMassEta->SetTextSize(0.02);
    legendMassEta->SetFillColor(0);
    legendMassEta->SetLineColor(0);
        legendMassEta->AddEntry(histoMassMesonEta,"Data");
        legendMassEta->AddEntry(histoTrueMassMesonEta,"MonteCarlo");
    legendMassEta->Draw();
//     if(!thesis)DrawAliceLogoPi0Performance(pictDrawingCoordinatesMass[0], pictDrawingCoordinatesMass[1], pictDrawingCoordinatesMass[2], pictDrawingCoordinatesMass[3], pictDrawingCoordinatesMass[4], pictDrawingCoordinatesMass[5], pictDrawingCoordinatesMass[6], pictDrawingCoordinatesMass[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1350,900,date, "MinBias",kDalitz);
    
    canvasMassEta->Update();
    if(!conference)canvasMassEta->SaveAs(Form("%s/%s_Mass_Eta_combined_%s_%s.%s",outputDir.Data(),prefix2.Data(),useSameBinningPi0Eta.Data(),cutSelection.Data(),suffix));
    delete canvasMassEta;

    
     if (isMC.CompareTo("kTRUE")==0){                    
        //**********************************************************************************
        //******************** Acceptance Plot *********************************************
        //**********************************************************************************
        
        TCanvas* canvasAcc = new TCanvas("canvasAcc","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasAcc, 0.1, 0.02, 0.02, 0.09);
        
        histoAccPi0->SetXTitle("p_{T} (GeV/c)");
        histoAccPi0->SetYTitle( Form("A (|y|< %s)",rapidityRange.Data()));
        histoAccPi0->GetXaxis()->SetLabelSize(0.03);
        histoAccPi0->GetYaxis()->SetLabelSize(0.03);
    if(mode == 4){
      histoAccPi0->GetYaxis()->SetRangeUser(0.0,1.5);
    }
    else{
      histoAccPi0->GetYaxis()->SetRangeUser(0.4,1.1);
    }
    histoAccPi0->GetYaxis()->SetTitleOffset(1.2);
        DrawGammaSetMarker(histoAccPi0, 22, 1.2, kBlack, kBlack);
        histoAccPi0->DrawCopy("e1"); 
        
        DrawGammaSetMarker(histoAccEta, 20, 1.2, kBlue, kBlue);
        histoAccEta->DrawCopy("e1,same"); 
        TLegend* legendAcc = new TLegend(0.3,0.15,0.5,0.25);
        legendAcc->SetTextSize(0.04);
        legendAcc->SetFillColor(0);
        legendAcc->SetLineColor(0);
        legendAcc->SetNColumns(2);
        legendAcc->AddEntry(histoAccPi0,"#pi^{0}");
        legendAcc->AddEntry(histoAccEta,"#eta");
        legendAcc->Draw();
        
//         if(!thesis)DrawAliceLogoCombined(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,kDalitz);
        canvasAcc->Update();
        
        canvasAcc->SaveAs(Form("%s/AcceptanceCombined_%s_%s.%s",outputDir.Data(),cutSelection.Data(),useSameBinningPi0Eta.Data(),suffix));
        delete canvasAcc;
        
        //**********************************************************************************
        //******************** Efficiency Simple Plot **************************************
        //**********************************************************************************
        TCanvas* canvasEffiSimple = new TCanvas("canvasEffiSimple","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasEffiSimple, 0.1, 0.02, 0.02, 0.09);
    if(mode != 4){
      canvasEffiSimple->SetLogy(1);
    }
    
    if(mode == 4){
      DrawAutoGammaMesonHistos( histoTrueEffPtPi0,
                               "", "p_{T} (GeV/c)", "#epsilon_{#pi^{0}(#eta)}",
                               kTRUE, 1., 0.001, kFALSE,
                               kFALSE, 0., 0.7,
                               kFALSE, 0., 10.);
      histoTrueEffPtPi0->GetYaxis()->SetRangeUser(0.001,0.5);
    }
    else{
      DrawAutoGammaMesonHistos( histoTrueEffPtPi0,
                               "", "p_{T} (GeV/c)", "#epsilon_{#pi^{0}(#eta)}",
                               kTRUE, 1., 3e-7, kFALSE,
                               kFALSE, 0., 0.7,
                               kFALSE, 0., 10.);
      histoTrueEffPtPi0->GetYaxis()->SetRangeUser(3e-7,1.1e-2);
    }
    
    DrawGammaSetMarker(histoTrueEffPtPi0, 22, 1.2, kBlack, kBlack);
        histoTrueEffPtPi0->DrawCopy("e1"); 
        
        DrawGammaSetMarker(histoTrueEffPtEta, 22, 1.2, kBlue, kBlue); 
        histoTrueEffPtEta->DrawCopy("e1,same"); 
        
        TLegend* legendEffi = new TLegend(0.3,0.15,0.5,0.22);
        legendEffi->SetTextSize(0.04);
        legendEffi->SetFillColor(0);
        legendEffi->SetLineColor(0);
        legendEffi->SetNColumns(2);
        legendEffi->AddEntry(histoTrueEffPtPi0,Form("#pi^{0}"));
        legendEffi->AddEntry(histoTrueEffPtEta,"#eta");
        legendEffi->Draw();
        
//         if(!thesis)DrawAliceLogoPi0MC(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,kDalitz);
        canvasEffiSimple->Update();
        
        canvasEffiSimple->SaveAs(Form("%s/EfficiencyCombined_%s_%s.%s",outputDir.Data(),cutSelection.Data(),useSameBinningPi0Eta.Data(),suffix));
        delete canvasEffiSimple;
        
    }//end if (isMC.CompareTo("kTRUE")==0)
    
    //**********************************************************************************
    //******************** RAW Yield spectrum ******************************************
    //**********************************************************************************
    TCanvas* canvasRawYield = new TCanvas("canvasRawYield","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRawYield, 0.08, 0.02, 0.02, 0.09);
    canvasRawYield->SetLogy(1);
    
    histoUncorrectedYieldPi0->Scale(1./nEvt);
    DrawAutoGammaMesonHistos( histoUncorrectedYieldPi0, 
                             "", "p_{T} (GeV/c)", "RAW Yield/ N_{Evt}", 
                         kTRUE, 4., 4e-10, kTRUE,
                             kFALSE, 0., 0.7, 
                             kFALSE, 0., 10.);
    
    DrawGammaSetMarker(histoUncorrectedYieldPi0, 20, 0.4, kBlack, kBlack); 
    histoUncorrectedYieldPi0->DrawCopy("e1"); 
    
    
    histoUnCorrectedYieldEta->Scale(1./nEvt);
    DrawGammaSetMarker(histoUnCorrectedYieldEta, 22, 1., kBlue, kBlue); 
    histoUnCorrectedYieldEta->DrawCopy("e1,same"); 
    
    TLegend* legendRawYield = new TLegend(0.25,0.11,0.4,0.2);
    legendRawYield->SetTextSize(0.02);
    legendRawYield->SetFillColor(0);
    legendRawYield->AddEntry(histoUncorrectedYieldPi0,Form("#pi^{0}"));
     legendRawYield->AddEntry(histoUnCorrectedYieldEta,"#eta");
    legendRawYield->Draw();
    
//     if(!thesis)DrawAliceLogoCombined(pictDrawingCoordinatesInv[0], pictDrawingCoordinatesInv[1], pictDrawingCoordinatesInv[2], pictDrawingCoordinatesInv[3], pictDrawingCoordinatesInv[4], pictDrawingCoordinatesInv[5], pictDrawingCoordinatesInv[6], pictDrawingCoordinatesInv[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,kDalitz);
    canvasRawYield->Update();
    
    canvasRawYield->SaveAs(Form("%s/%s_RAWYieldPtCombined_%s_%s.%s",outputDir.Data(),prefix2.Data(),useSameBinningPi0Eta.Data(),cutSelection.Data(),suffix));
    delete canvasRawYield;
    
    
    TCanvas* canvasCorrectedSame = new TCanvas("canvasCorrectedSame","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasCorrectedSame, 0.08, 0.02, 0.02, 0.09);
    canvasCorrectedSame->SetLogy();
    
    DrawAutoGammaMesonHistos( histoCorrectedYieldPi0, 
                             "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
                             kTRUE, 3., 4e-10, kTRUE,
                             kFALSE, 3e-8,10, 
                             kFALSE, 0., 10.);
    
    DrawGammaSetMarker(histoCorrectedYieldPi0, 22, 1., kBlack, kBlack);
    histoCorrectedYieldPi0->DrawCopy("e1"); 
    
//     if(!thesis)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1350,900,kDalitz);
    canvasCorrectedSame->Update();
    
    canvasCorrectedSame->SaveAs(Form("%s/%s_CorrectedYieldOnlyPi0_%s.%s",outputDir.Data(),prefix2.Data(),cutSelection.Data(),suffix));
    delete canvasCorrectedSame;
    
    
    
    nameFinalResDat = Form("%s/%s/%s_FinalExtraction_%s_%s_%s.dat",cutSelection.Data(),optionEnergy.Data(), prefix2.Data(), cutSelection.Data(), useSameBinningPi0Eta.Data(), makeBinShiftWithFunction.Data());
    
    fileFinalResults.open(nameFinalResDat, ios::out);
    
    
    //*****************************************************************************************************
    //**************************** Binshifted Spectra *****************************************************
    //*****************************************************************************************************
    
    TCanvas* canvasBinShifted = new TCanvas("canvasBinShifted","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasBinShifted, 0.13, 0.02, 0.02, 0.09);
    canvasBinShifted->SetLogy();
    
    TPad* padBinShiftedHistos = new TPad("padBinShiftedHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padBinShiftedHistos, 0.12, 0.02, 0.02, 0.);
    padBinShiftedHistos->Draw();
    
    TPad* padBinShiftedRatios = new TPad("padBinShiftedRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padBinShiftedRatios, 0.12, 0.02, 0., 0.2);
    padBinShiftedRatios->Draw();
    
    padBinShiftedHistos->cd();
    padBinShiftedHistos->SetLogy();
    
    histoRecBinShiftPi0 = (TH1D*) histoCorrectedYieldPi0->Clone();
    nb = histoRecBinShiftPi0->GetXaxis()->GetNbins();
    histoRecBinShiftPi0->SetMarkerStyle(20);
    histoRecBinShiftPi0->SetMarkerSize(0.9);
    histoRecBinShiftPi0->SetMarkerColor(kGray+1);
    //histoRecBinShiftPi0->DrawCopy("e1"); 
    
    histoPi0CorrYieldBinShifted= (TH1D*)histoRecBinShiftPi0->Clone();
    if(kHag){
        fitBinShifting = FitObject("h","fitBinShiftingPi0","Pi0");
        legendEntryFunction= "Hagedorn (FitRange: 0.4 - maximum pt)";
    }
    else if(kLevy){
        fitBinShifting = FitObject("l","fitBinShiftingPi0","Pi0");
        legendEntryFunction= "Levy (FitRange: 0.4 - maximum pt)";
    }    
    fitBinShifting->SetRange(minPtForFits,maxPtPi0);

    TGraphAsymmErrors* graphCorrectedYieldXShifted = new TGraphAsymmErrors(histoRecBinShiftPi0);
    graphCorrectedYieldXShifted->RemovePoint(0);
    Int_t nPointsPi0X = graphCorrectedYieldXShifted->GetN();
    graphCorrectedYieldXShifted->Print();

    Double_t ratio[50];
    Double_t ratioErrUp[50];
    Double_t ratioErrDown[50];
    
    TF1* fitBinShiftingUpErr = FitObject("l","fitBinShiftingPi0UpErr","Pi0");
    TF1* fitBinShiftingDownErr = FitObject("l","fitBinShiftingPi0DownErr","Pi0");
    
    while(TMath::Abs(globalRatio-testGlobalRatio) > 0.00001 ){    
        testGlobalRatio = globalRatio;
        globalRatio = 0.;
        color++;
        // fit acc+eff corrected yield with Hagedorn function
        // apply bin shift correction
        histoPi0CorrYieldBinShifted->Fit(fitBinShifting,"QRME0");
        for ( Int_t i = 0; i < 3; i++){
            fitBinShiftingUpErr->SetParameter(i, fitBinShifting->GetParameter(i)+fitBinShifting->GetParError(i));
            fitBinShiftingDownErr->SetParameter(i, fitBinShifting->GetParameter(i)-fitBinShifting->GetParError(i));
        }
    
        TFitResultPtr r = histoPi0CorrYieldBinShifted->Fit(fitBinShifting, "S E M 0", "", minPtForFits,maxPtPi0);
        for (int ib=2; ib<=nb; ib++) {    
            Double_t ptMin = histoRecBinShiftPi0->GetBinLowEdge(ib);
            Double_t ptMax = ptMin + histoRecBinShiftPi0->GetBinWidth(ib);
        
            // the bin shift affected value of the fit function in current bin
            Double_t shiftedValue = fitBinShifting->Integral(ptMin,ptMax) / (ptMax - ptMin);
            Double_t shiftedValueErr = fitBinShifting->IntegralError(ptMin, ptMax)/ (ptMax - ptMin);
            
            // the correct value at the bin center
            Double_t trueValue = fitBinShifting->Eval((ptMax + ptMin)/2.);
            Double_t trueValueErrUp = trueValue - fitBinShiftingUpErr->Eval((ptMax + ptMin)/2.);
            Double_t trueValueErrDown = trueValue - fitBinShiftingDownErr->Eval((ptMax + ptMin)/2.);
            
            cout << shiftedValue << "\t"<< shiftedValueErr << "\t"<< trueValue << "\t"<< trueValueErrUp << "\t" << trueValueErrDown << endl;
            
            // the bin shift correction factor
            ratio[ib] =  shiftedValue / trueValue;
            ratioErrUp[ib] = TMath::Sqrt( pow(shiftedValueErr/trueValue,2)  + pow( trueValueErrUp*shiftedValue/pow(trueValue,2),2) );
            ratioErrDown[ib] = TMath::Sqrt( pow(shiftedValueErr/trueValue,2)  + pow( trueValueErrDown*shiftedValue/pow(trueValue,2),2) );
            
            float pt = histoRecBinShiftPi0->GetBinCenter(ib);
            
            fileFinalResults<< "BinShiftCorrection, pt = " << pt << " GeV : " << ratio[ib] << "\t+"<< ratioErrUp[ib] << "\t-" << ratioErrDown[ib] << endl; 
            
            histoPi0CorrYieldBinShifted->SetBinContent(ib,  histoRecBinShiftPi0->GetBinContent(ib) / ratio[ib]);
            histoPi0CorrYieldBinShifted->SetBinError(ib,  histoRecBinShiftPi0->GetBinError(ib) / ratio[ib]);
            
            globalRatio = globalRatio + ratio[ib];
        }
        
        globalRatio = globalRatio/(nb-1);
        fileFinalResults << color << "\t globalRatio: \t"<< globalRatio  << "\t testGlobalRatio:\t "<< testGlobalRatio << endl;
        fileFinalResults << endl << endl;

        histoPi0CorrYieldBinShifted->SetMarkerStyle(24);
        histoPi0CorrYieldBinShifted->SetMarkerSize(0.9);
        histoPi0CorrYieldBinShifted->SetMarkerColor(kBlack);
        DrawAutoGammaMesonHistos( histoPi0CorrYieldBinShifted, 
                             "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
                             kFALSE, 3., 4e-10, kTRUE,
                             kFALSE, 3e-8,10, 
                             kFALSE, 0., 10.);
    
        DrawGammaSetMarkerTF1( fitBinShifting, 1, 0.4, kBlue-4);
        fitBinShifting->DrawCopy("same");
        
        Double_t parameter[3];
        fitBinShifting->GetParameters(parameter);
        fitBinShifting->SetParameter(0,parameter[0]);  
        fitBinShifting->SetParameter(1,parameter[1]); 
        fitBinShifting->SetParameter(2,parameter[2]); 
    }

    histoRecBinShiftPi0->DrawCopy("same,e1"); 
    
    Double_t* xValues = graphCorrectedYieldXShifted->GetX();
//     Double_t* YValues = graphCorrectedYieldXShifted->GetY();
    Double_t* xErrorsLow = graphCorrectedYieldXShifted->GetEXlow();
    Double_t* xErrorsHigh = graphCorrectedYieldXShifted->GetEXhigh();
    for (Int_t i = 0; i < nPointsPi0X; i++){
        Double_t xValuesNew = xValues[i]+ xErrorsLow[i]*2*(1-ratio[i+2]);
        Double_t xErrorsLowNew =xErrorsLow[i]*2*(ratioErrUp[i+2]);
        Double_t xErrorsHighNew = xErrorsLow[i]*2*(ratioErrDown[i+2]);
        xValues[i] = xValuesNew;
        xErrorsLow[i] = xErrorsLowNew;
        xErrorsHigh[i] = xErrorsHighNew;
        
    }
    cout << "und hier korrigiert" << endl;
    graphCorrectedYieldXShifted->Print();
    
    
    cout << "fitting xShifted with Levy" << endl;
    TF1* fitPtLevyPi0X = FitObject("l","fitPtLevyPi0X","Pi0",graphCorrectedYieldXShifted,minPtForFits,maxPtPi0,NULL,"QNRME+");
    DrawGammaSetMarkerTF1(fitPtLevyPi0X, 1, 1.5, kBlue);
    kLevySuccPi0 = kTRUE;
    fileFinalResults << "fitting xShifted with Levy" << endl;
    forOutput= WriteParameterToFile(fitPtLevyPi0X);
    fileFinalResults<< forOutput.Data()<< endl;
    
    
    forOutput= WriteParameterToFile(fitBinShifting);
    fileFinalResults<< forOutput.Data()<< endl;
    
    /*DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldXShifted, 25, 0.7, kRed, kRed);
    graphCorrectedYieldXShifted->Draw("same,pe");
    */
    
    
    TLegend* legendYieldBinShifted = new TLegend(0.2,0.05,0.5,0.2);
    legendYieldBinShifted->SetFillColor(0);
    legendYieldBinShifted->SetLineColor(0);
    legendYieldBinShifted->AddEntry(histoRecBinShiftPi0,"Yield","p");
    legendYieldBinShifted->AddEntry(histoPi0CorrYieldBinShifted,"Yield Corr (iter)","p");
    legendYieldBinShifted->AddEntry(fitBinShifting,legendEntryFunction,"l");
    legendYieldBinShifted->Draw();
    
//     if(!thesis)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1350,1125,kDalitz);
    
    padBinShiftedRatios->cd();
    padBinShiftedRatios->SetLogy(0);
    
    histoRatioBinShiftedPi0 = (TH1D*) histoRecBinShiftPi0->Clone();
    histoRatioBinShiftedPi0->Divide(histoRatioBinShiftedPi0,histoPi0CorrYieldBinShifted,1.,1.,"");
    for (Int_t i = 0; i < (histoRatioBinShiftedPi0->GetNbinsX()+1); i++){
        histoRatioBinShiftedPi0->SetBinError(i,0.);
    }
    histoRatioBinShiftedPi0->SetYTitle("#frac{value}{shifted value}");
    if(useSameBinningPi0Eta.CompareTo("same") == 0) histoRatioBinShiftedPi0->GetYaxis()->SetRangeUser(0.95,1.42);
    else histoRatioBinShiftedPi0->GetYaxis()->SetRangeUser(0.95,1.22);
    histoRatioBinShiftedPi0->GetYaxis()->SetNdivisions(505);
    histoRatioBinShiftedPi0->GetYaxis()->SetLabelSize(0.08);
    histoRatioBinShiftedPi0->GetYaxis()->SetTitleSize(0.1);
    histoRatioBinShiftedPi0->GetYaxis()->SetDecimals();
    histoRatioBinShiftedPi0->GetYaxis()->SetTitleOffset(0.42);
    histoRatioBinShiftedPi0->GetXaxis()->SetTitleSize(0.11);
    histoRatioBinShiftedPi0->GetXaxis()->SetLabelSize(0.08);
    DrawGammaSetMarker(histoRatioBinShiftedPi0, 20, 0.5, kBlack, kBlack);
    histoRatioBinShiftedPi0->DrawCopy("p"); 
    
    DrawGammaLines(0., maxPtPi0,1., 1.,0.1);
    
    canvasBinShifted->Update();
    
    canvasBinShifted->SaveAs(Form("%s/Pi0_%s_CorrectedYieldBinShifted_%s_%s_%s.%s",outputDir.Data(),prefix2.Data(),makeBinShiftWithFunction.Data(),cutSelection.Data(),useSameBinningPi0Eta.Data(),suffix));
    
    //**************************************************************************************************
    //**************************** Bin shifting with Integral optionEnergy ***********************************
    //**************************************************************************************************
    
    
    
    
    //***************************************************************************************************
    //*************************** Fitting  Pi0 Spectrum *************************************************
    //***************************************************************************************************
    
    fileFinalResults << "Final Results for Pi0" << endl << endl;
    
    //****************************************************************************************************
    //************************** Fitting of corrected normal spectrum ************************************
    //****************************************************************************************************
    
    TCanvas* canvasFittingPi0 = new TCanvas("canvasFittingPi0","",200,10,1350,900);  // gives the page size
    
    histoFittingPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
    
    //****************************** Fit histCorr with Levy *****************************************
    cout << "fitting with Levy" << endl;
    fitPtLevyPi0 = FitObject("l","fitPtLevyPi0","Pi0",histoFittingPi0,minPtForFits,maxPtPi0,NULL,"QNRME+");
    DrawGammaSetMarkerTF1(fitPtLevyPi0, 1, 1.5, kBlue);
    kLevySuccPi0 = kTRUE;
    forOutput= WriteParameterToFile(fitPtLevyPi0);
    fileFinalResults<< forOutput.Data()<< endl;
    
    //**************************** Fit histCorr with Hagedorn  **********************************
    cout << "fitting with Hagedorn" << endl;
    fitPtHagedornPi0 = FitObject("h","fitPtHagedornPi0","Pi0",histoFittingPi0,minPtForFits,maxPtPi0,NULL,"QNRME+");
    DrawGammaSetMarkerTF1( fitPtHagedornPi0, 1, 1.5, kGreen+2);
    fitPtHagedornPi0->SetLineStyle(7);
    kHagSuccPi0 = kTRUE;
    forOutput= WriteParameterToFile(fitPtHagedornPi0);
    fileFinalResults<< forOutput.Data()<< endl;
    
    //************************** Fit histCorr with Boltzmann  ***********************************
    if (useSameBinningPi0Eta.CompareTo("")==0){         
        fitPtBoltzmannPi0 = FitObject("b","fitPtBoltzmannPi0","Pi0",histoFittingPi0,0.3,2.,NULL,"QNRME+");
        DrawGammaSetMarkerTF1( fitPtBoltzmannPi0, 1, 1.5, kMagenta+4);
        kBoltzSuccPi0 = kTRUE;
        forOutput= WriteParameterToFile(fitPtBoltzmannPi0);
        fileFinalResults<< forOutput.Data()<< endl;
    }    
    
    //*************************** Fit histCorr with Exponential **********************************
    if (useSameBinningPi0Eta.CompareTo("")==0){
        fitPtExpPi0 = FitObject("e","fitPtExpPi0","Pi0",histoFittingPi0,0.3,2.,NULL,"QNRME+");
        DrawGammaSetMarkerTF1( fitPtExpPi0, 1, 1.5, kOrange+7);
        kExpSuccPi0 = kTRUE;
        forOutput= WriteParameterToFile(fitPtExpPi0);
        fileFinalResults<< forOutput.Data()<< endl;
    }
    
    //**************************** Fit histCorr with Powerlaw  *********************************
    fitPtPowerlawPi0 = FitObject("p","fitPtPowerlawPi0","Pi0",histoFittingPi0,1.5,maxPtPi0,NULL,"QNRME+");
    DrawGammaSetMarkerTF1( fitPtPowerlawPi0, 1, 1.5, kTeal);
    kPowSuccPi0 = kTRUE;
    forOutput= WriteParameterToFile(fitPtPowerlawPi0);
    fileFinalResults<< forOutput.Data()<< endl;
    
    //**************************** Fit histCorr with ModPowerlaw  *********************************
    fitPtModPowerlawPi0 = FitObject("m","fitPtModPowerlawPi0","Pi0",histoFittingPi0,0.3,maxPtPi0,NULL,"QNRME+");
    DrawGammaSetMarkerTF1( fitPtModPowerlawPi0, 1, 1.5, kMagenta+2);
    kModPowSuccPi0 = kTRUE;
    forOutput= WriteParameterToFile(fitPtModPowerlawPi0);
    fileFinalResults<< forOutput.Data()<< endl;

    //**************************** Fit histCorr with two component model (Bylinkin)  *********************************
    //fitPtTwoCompModelPi0 = FitObject("tcmp","fitPtTwoCompModelPi0","Pi0",histoFittingPi0,minPtForFits,maxPtPi0,NULL,"QNRME+");
    //DrawGammaSetMarkerTF1( fitPtTwoCompModelPi0, 1, 1.5, kRed);
    //kTwoCompModelSuccPi0 = kTRUE;
    //forOutput= WriteParameterToFile(fitPtTwoCompModelPi0);
    //fileFinalResults<< forOutput.Data()<< endl;

    //*************************** Calculating Ratios *******************************************
    histoRatioFitLevyPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
    histoRatioFitHagPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
    histoRatioFitBoltzPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
    histoRatioFitExpPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
    histoRatioFitPowPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
    histoRatioFitModPowPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
    //histoRatioFitTwoCompModelPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();

    if(kLevySuccPi0) {
        histoRatioFitLevyPi0 = CalculateHistoRatioToFit (histoRatioFitLevyPi0, fitPtLevyPi0); 
    }
    if(kHagSuccPi0) {
        histoRatioFitHagPi0 = CalculateHistoRatioToFit (histoRatioFitHagPi0, fitPtHagedornPi0);
    }
    if (useSameBinningPi0Eta.CompareTo("")==0){
        if(kBoltzSuccPi0) {
            histoRatioFitBoltzPi0 = CalculateHistoRatioToFit (histoRatioFitBoltzPi0, fitPtBoltzmannPi0);
        }
        if(kExpSuccPi0) {
            histoRatioFitExpPi0 = CalculateHistoRatioToFit (histoRatioFitExpPi0, fitPtExpPi0);
        }
    }
    if(kPowSuccPi0) {
        histoRatioFitPowPi0 = CalculateHistoRatioToFit (histoRatioFitPowPi0, fitPtPowerlawPi0);
    }    
    if(kModPowSuccPi0) {
        histoRatioFitModPowPi0 = CalculateHistoRatioToFit (histoRatioFitModPowPi0, fitPtModPowerlawPi0);
    }
    //if (kTwoCompModelSuccPi0){
    //    histoRatioFitTwoCompModelPi0 = CalculateHistoRatioToFit (histoRatioFitTwoCompModelPi0, fitPtTwoCompModelPi0);
    //}

    delete canvasFittingPi0;
    //**********************************************************************************************
    //********************** Plotting of fitted spectrum *******************************************
    //**********************************************************************************************    
    TCanvas* canvasFittingSpectra = new TCanvas("canvasFittingSpectra","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasFittingSpectra, 0.13, 0.02, 0.02, 0.09);
    canvasFittingSpectra->SetLogy();
    
    TPad* padFittedSpectraHistos = new TPad("padFittedSpectraHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padFittedSpectraHistos, 0.12, 0.02, 0.02, 0.);
    padFittedSpectraHistos->Draw();
    
    TPad* padFittedSpectraRatios = new TPad("padFittedSpectraRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padFittedSpectraRatios, 0.12, 0.02, 0., 0.26);
    padFittedSpectraRatios->Draw();
    
    padFittedSpectraHistos->cd();
    padFittedSpectraHistos->SetLogy();

    DrawAutoGammaMesonHistos( histoPi0CorrYieldBinShifted, 
                             "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
                             kTRUE, 3., 4e-10, kTRUE,
                             kFALSE, 3e-8,10, 
                             kFALSE, 0., 10.);
    
    DrawGammaSetMarker(histoPi0CorrYieldBinShifted, 22, 1., kBlack, kBlack);
    histoCorrectedYieldPi0->Draw("e1,x0"); 
    
    fitPtLevyPi0->Draw("same");
    fitPtHagedornPi0->Draw("same");
    if (useSameBinningPi0Eta.CompareTo("")==0){
        fitPtBoltzmannPi0->Draw("same");
        fitPtExpPi0->Draw("same");
    }
    fitPtPowerlawPi0->Draw("same");
    //fitPtModPowerlawPi0->Draw("same");
    //fitPtTwoCompModelPi0->Draw("same");
    
    TLegend* legendFit = new TLegend(0.15,0.02,0.65,0.15);
    legendFit->SetTextSize(0.02);
    legendFit->SetFillColor(0);
    legendFit->AddEntry(histoPi0CorrYieldBinShifted,"#pi^{0}");
    legendFit->AddEntry(fitPtLevyPi0,"Levy fit");
    legendFit->AddEntry(fitPtHagedornPi0,"Hagedorn fit");
    if (useSameBinningPi0Eta.CompareTo("")==0){
        legendFit->AddEntry(fitPtBoltzmannPi0,"Boltzman fit");
        legendFit->AddEntry(fitPtExpPi0,"Exponential fit");
    }
    legendFit->AddEntry(fitPtPowerlawPi0,"Powerlaw fit");
    //legendFit->AddEntry(fitPtModPowerlawPi0,"ModPowerlaw fit");
    //legendFit->AddEntry(fitPtTwoCompModelPi0, "Two component model (Bylinkin-Rostovtsev)");
    legendFit->Draw();
    
//     if(!thesis)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1350,1125,kDalitz);
    canvasFittingSpectra->Update();
    
    padFittedSpectraRatios->cd();
    
    histoRatioFitLevyPi0->SetYTitle("#frac{value}{fit}");
    histoRatioFitLevyPi0->GetYaxis()->SetRangeUser(0.5,1.55);
    histoRatioFitLevyPi0->GetYaxis()->SetLabelSize(0.08);
    histoRatioFitLevyPi0->GetYaxis()->SetTitleSize(0.1);
    histoRatioFitLevyPi0->GetYaxis()->SetDecimals();
    histoRatioFitLevyPi0->GetYaxis()->SetTitleOffset(0.4);
    histoRatioFitLevyPi0->GetXaxis()->SetTitleSize(0.11);
    histoRatioFitLevyPi0->GetXaxis()->SetLabelSize(0.08);
    
    DrawGammaSetMarker(histoRatioFitLevyPi0, 21, 0.5, kBlue, kBlue);
    histoRatioFitLevyPi0->Draw("e1,x0");
    
    DrawGammaSetMarker(histoRatioFitHagPi0, 21, 0.5, kGreen+2, kGreen+2);
    if(kHagSuccPi0)histoRatioFitHagPi0->Draw("e1,x0,same");
    
    if (useSameBinningPi0Eta.CompareTo("")==0){
        DrawGammaSetMarker(histoRatioFitBoltzPi0, 21, 0.5, kMagenta+4, kMagenta+4);
        if(kBoltzSuccPi0)histoRatioFitBoltzPi0->Draw("e1,x0,same");
        
        DrawGammaSetMarker(histoRatioFitExpPi0, 21, 0.5, kOrange+7, kOrange+7);
        if(kExpSuccPi0)histoRatioFitExpPi0->Draw("e1,x0,same");
    }    
    DrawGammaSetMarker(histoRatioFitPowPi0, 21, 0.5, kTeal, kTeal);
    if(kPowSuccPi0)histoRatioFitPowPi0->Draw("e1,x0,same");
    
    DrawGammaSetMarker(histoRatioFitModPowPi0, 21, 0.5, kMagenta+2, kMagenta+2);
    //if(kModPowSuccPi0)histoRatioFitModPowPi0->Draw("e1,same");
    
    //DrawGammaSetMarker(histoRatioFitTwoCompModelPi0, 21, 0.5, kRed, kRed);
    //if(kTwoCompModelSuccPi0)histoRatioFitTwoCompModelPi0->Draw("e1,x0,same");

    DrawGammaLines(0., maxPtPi0 ,1., 1.,0.1);
    
    canvasFittingSpectra->SaveAs(Form("%s/Pi0_%s_CorrectedYieldFitted_%s_%s.%s",outputDir.Data(),prefix2.Data(),cutSelection.Data(),useSameBinningPi0Eta.Data(),suffix));
    delete canvasFittingSpectra;
    

        
    //-------------------------- Fit histCorr with Levy  ----------------------------------------
    fitPtLevyEta = (TF1*)fitPtLevyPi0->Clone("fitPtLevyEta");
    if(conference){
        fitPtLevyEta->SetRange(0.4,maxPtEta);
    }else{
        fitPtLevyEta->SetRange(0.6,maxPtEta);
    }

    //-------------------------- Fit histCorr with Hagedorn  ------------------------------------
    fitPtHagedornEta = (TF1*)fitPtHagedornPi0->Clone("fitPtHagedornEta");
    if(conference){
        fitPtHagedornEta->SetRange(0.4,maxPtEta);
    }else{
        fitPtHagedornEta->SetRange(0.6,maxPtEta);
    }
    
    //*****************************************************************************************************
    //**************************** Binshifted Spectra *****************************************************
    //*****************************************************************************************************

    TCanvas* canvasBinShiftedEta = new TCanvas("canvasBinShiftedEta","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasBinShiftedEta, 0.13, 0.02, 0.02, 0.09);
    canvasBinShiftedEta->SetLogy();
    
    TPad* padBinShiftedEtaHistos = new TPad("padBinShiftedEtaHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padBinShiftedEtaHistos, 0.12, 0.02, 0.02, 0.);
    padBinShiftedEtaHistos->Draw();
    
    TPad* padBinShiftedEtaRatios = new TPad("padBinShiftedEtaRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padBinShiftedEtaRatios, 0.12, 0.02, 0., 0.2);
    padBinShiftedEtaRatios->Draw();
    
    padBinShiftedEtaHistos->cd();
    padBinShiftedEtaHistos->SetLogy();
    
    histoRecBinShiftEta = (TH1D*) histoCorrectedYieldEta->Clone();
    nbEta = histoRecBinShiftEta->GetXaxis()->GetNbins();
    
    histoRecBinShiftEta->SetMarkerStyle(20);
    histoRecBinShiftEta->SetMarkerSize(0.9);
    histoRecBinShiftEta->SetMarkerColor(kGray+1);
    
    
    histoEtaCorrYieldBinShifted=(TH1D*)histoRecBinShiftEta->Clone();
    if(kHag){
        fitBinShiftingEta = (TF1*)fitPtHagedornEta->Clone("fitBinShiftingEta");
        legendEntryFunctionEta= "Hagedorn (FitRange: 0.4 - max pt)";
    }
    if(kLevy){
        fitBinShiftingEta = (TF1*)fitPtLevyEta->Clone("fitBinShiftingEta");
        legendEntryFunctionEta= "Levy (FitRange: 0.4 - max pt)";
    }    
    if( optionEnergy.CompareTo("900GeV") != 0) {
    
        color = 1;
        while(TMath::Abs(globalRatioEta-testGlobalRatioEta) > 0.00001 ){    
            
            testGlobalRatioEta = globalRatioEta;
            globalRatioEta = 0.;
            color++;
            // fit acc+eff corrected yield with Hagedorn function
            histoEtaCorrYieldBinShifted->Fit(fitBinShiftingEta,"QRME0");
            
            // apply bin shift correction
            for (int ib=2; ib<=nbEta; ib++) {
                
                double ptMin = histoRecBinShiftEta->GetBinLowEdge(ib);
                double ptMax = ptMin + histoRecBinShiftEta->GetBinWidth(ib);
                
                // the bin shift affected value of the fit function in current bin            
                double shiftedValue = fitBinShiftingEta->Integral(ptMin,ptMax) / (ptMax - ptMin);
                
                // the correct value at the bin center
                double trueValue = fitBinShiftingEta->Eval((ptMax + ptMin)/2.);
                
                // the bin shift correction factor
                double ratioBS =  shiftedValue / trueValue;
                
                float pt = histoRecBinShiftEta->GetBinCenter(ib);
                
                //                 cout << "BinShiftCorrection, pt = " << pt << " GeV : " << ratioBS << endl; 
                fileFinalResults<< "BinShiftCorrection, pt = " << pt << " GeV : " << ratioBS << endl; 
                
                histoEtaCorrYieldBinShifted->SetBinContent(ib,  histoRecBinShiftEta->GetBinContent(ib) / ratioBS);
                histoEtaCorrYieldBinShifted->SetBinError(ib,  histoRecBinShiftEta->GetBinError(ib) / ratioBS);
                
                globalRatioEta = globalRatioEta + ratioBS;
            }
            
            globalRatioEta = globalRatioEta/(nbEta-1);
            
            //             cout << color << "\t globalRatioEta: \t"<< globalRatioEta  << "\t testGlobalRatioEta:\t "<< testGlobalRatioEta << endl;
            fileFinalResults << color << "\t globalRatioEta: \t"<< globalRatioEta  << "\t testGlobalRatioEta:\t "<< testGlobalRatioEta << endl;
            fileFinalResults << endl << endl;
            
            histoEtaCorrYieldBinShifted->SetMarkerStyle(24);
            histoEtaCorrYieldBinShifted->SetMarkerSize(0.9);
            histoEtaCorrYieldBinShifted->SetMarkerColor(kBlack);
            DrawAutoGammaMesonHistos( histoEtaCorrYieldBinShifted, 
                             "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
                             kFALSE, 3., 4e-10, kTRUE,
                             kTRUE, 3e-8,10, 
                             kFALSE, 0., 10.);
            
            DrawGammaSetMarkerTF1( fitBinShiftingEta, 1, 0.4, kBlue-4);
            fitBinShiftingEta->DrawCopy("same");
            
            Double_t parameter[3];
            fitBinShiftingEta->GetParameters(parameter);
            fitBinShiftingEta->SetParameter(0,parameter[0]);  
            fitBinShiftingEta->SetParameter(1,parameter[1]); 
            fitBinShiftingEta->SetParameter(2,parameter[2]); 
        }
        
        forOutput = WriteParameterToFile(fitBinShiftingEta);
        fileFinalResults<< forOutput.Data()<< endl;
    }
    
    histoRecBinShiftEta->DrawCopy("same,e1"); 

    TLegend* legendYieldBinShiftedEta = new TLegend(0.2,0.05,0.5,0.2);
    legendYieldBinShiftedEta->SetFillColor(0);
    legendYieldBinShiftedEta->AddEntry(histoRecBinShiftEta,"Yield","p");
    legendYieldBinShiftedEta->AddEntry(histoEtaCorrYieldBinShifted,"Yield Corr","p");
    legendYieldBinShiftedEta->AddEntry(fitBinShiftingEta,legendEntryFunctionEta,"l");
    legendYieldBinShiftedEta->Draw();
    
//     if(!thesis)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1350,1125,kDalitz);
    
    padBinShiftedEtaRatios->cd();
    padBinShiftedEtaRatios->SetLogy(0);

    histoRatioBinShiftedEta = (TH1D*) histoRecBinShiftEta->Clone();
    histoRatioBinShiftedEta->Divide(histoRatioBinShiftedEta,histoEtaCorrYieldBinShifted,1.,1.,"");
    for (Int_t i = 0; i < (histoRatioBinShiftedEta->GetNbinsX()+1); i++){
        histoRatioBinShiftedEta->SetBinError(i,0.);
    }
    histoRatioBinShiftedEta->SetYTitle("#frac{value}{shifted value}");
    histoRatioBinShiftedEta->GetYaxis()->SetRangeUser(0.95,1.42);
    histoRatioBinShiftedEta->GetYaxis()->SetNdivisions(505);
    histoRatioBinShiftedEta->GetYaxis()->SetLabelSize(0.08);
    histoRatioBinShiftedEta->GetYaxis()->SetTitleSize(0.1);
    histoRatioBinShiftedEta->GetYaxis()->SetDecimals();
    histoRatioBinShiftedEta->GetYaxis()->SetTitleOffset(0.42);
    histoRatioBinShiftedEta->GetXaxis()->SetTitleSize(0.11);
    histoRatioBinShiftedEta->GetXaxis()->SetLabelSize(0.08);
    DrawGammaSetMarker(histoRatioBinShiftedEta, 20, 0.5, kBlack, kBlack);
    histoRatioBinShiftedEta->DrawCopy("p"); 
    
    DrawGammaLines(0., maxPtEta,1., 1.,0.1);
    
    canvasBinShiftedEta->Update();
    
    canvasBinShiftedEta->SaveAs(Form("%s/Eta_%s_CorrectedYieldBinShifted_%s_%s.%s",outputDir.Data(),prefix2.Data(),makeBinShiftWithFunction.Data(),cutSelection.Data(),suffix));
    
    //***************************************************************************************************
    //*************************** Fitting  Eta Spectrum *************************************************
    //***************************************************************************************************

    fileFinalResults << "Final Results for Eta" << endl << endl;
    
    //****************************************************************************************************
    //************************** Fitting of corrected normal spectrum ************************************
    //****************************************************************************************************
    
    TCanvas* canvasFittingEta = new TCanvas("canvasFittingEta","",200,10,1350,900);  // gives the page size
    
    histoFittingEta = (TH1D*) histoEtaCorrYieldBinShifted->Clone();
    //-------------------------- Fit histCorr with Levy  ----------------------------------------
    fitPtLevyEta = (TF1*)fitPtLevyPi0->Clone("fitPtLevyEta");
    fitPtLevyEta->SetRange(minPtForFitsEta,maxPtEta);
    histoFittingEta->Fit(fitPtLevyEta,"QRME+"); 
    DrawGammaSetMarkerTF1( fitPtLevyEta, 1, 1.5, kBlue);
    kLevySuccEta = kTRUE;
    forOutput = WriteParameterToFile(fitPtLevyEta);
    fileFinalResults<< forOutput.Data()<< endl;
    
    //-------------------------- Fit histCorr with Hagedorn  ------------------------------------
    fitPtHagedornEta = (TF1*)fitPtHagedornPi0->Clone("fitPtHagedornEta");
    fitPtHagedornEta->SetRange(minPtForFitsEta,maxPtEta);
    histoFittingEta->Fit(fitPtHagedornEta,"QRME+"); 
    DrawGammaSetMarkerTF1( fitPtHagedornEta, 1, 1.5, kGreen+2);
    fitPtHagedornEta->SetLineStyle(7);
    kHagSuccEta = kTRUE;
    forOutput= WriteParameterToFile(fitPtHagedornEta);
    fileFinalResults<< forOutput.Data()<< endl;
    
    //-------------------------- Fit histCorr with Powerlaw  ------------------------------------
    fitPtPowerlawEta = (TF1*)fitPtPowerlawPi0->Clone("fitPtPowerlawEta");
    fitPtPowerlawEta->SetRange(minPtForFitsEta,maxPtEta);
    histoFittingEta->Fit(fitPtPowerlawEta,"QRME+"); 
    DrawGammaSetMarkerTF1( fitPtPowerlawEta, 1, 1.5, kTeal);
    kPowSuccEta = kTRUE;
    forOutput= WriteParameterToFile(fitPtPowerlawEta);
    fileFinalResults<< forOutput.Data()<< endl;
    
    //**************************** Fit histCorr with ModPowerlaw  *********************************
    fitPtModPowerlawEta = (TF1*)fitPtModPowerlawPi0->Clone("fitPtModPowerlawEta");
    fitPtModPowerlawEta->SetRange(minPtForFitsEta,maxPtEta);
    histoFittingEta->Fit(fitPtModPowerlawEta,"QRME+"); 
    DrawGammaSetMarkerTF1( fitPtModPowerlawEta, 1, 1.5, kMagenta+2);
    kModPowSuccEta = kTRUE;
    forOutput = WriteParameterToFile(fitPtModPowerlawEta);
    fileFinalResults<< forOutput.Data()<< endl;

    //**************************** Fit histCorr with two component model  *********************************
    //fitPtTwoCompModelEta = (TF1*)fitPtTwoCompModelPi0->Clone("fitPtTwoCompModelEta");
    //fitPtTwoCompModelEta->SetRange(minPtForFitsEta,maxPtEta);
    //histoFittingEta->Fit(fitPtTwoCompModelEta,"QRME+"); 
    //DrawGammaSetMarkerTF1(fitPtTwoCompModelEta, 1, 1.5, kRed);
    //kTwoCompSuccEta = kTRUE;
    //forOutput = WriteParameterToFile(fitPtTwoCompModelEta);
    //fileFinalResults<< forOutput.Data()<< endl;
    //----
    //fitPtTwoCompModelEta = FitObject("tcmp","fitPtTwoCompModelEta","Eta",histoFittingEta,minPtForFitsEta,maxPtEta,NULL,"QNRME+");
    //DrawGammaSetMarkerTF1( fitPtTwoCompModelEta, 1, 1.5, kRed);
    //kTwoCompSuccEta = kTRUE;
    //forOutput= WriteParameterToFile(fitPtTwoCompModelEta);
    //fileFinalResults<< forOutput.Data()<< endl;
    
    //************************** Calculating Ratios ***********************************************
    histoRatioFitLevyEta = (TH1D*) histoEtaCorrYieldBinShifted->Clone();
    histoRatioFitHagEta = (TH1D*) histoEtaCorrYieldBinShifted->Clone();
    histoRatioFitPowEta = (TH1D*) histoEtaCorrYieldBinShifted->Clone();
    histoRatioFitModPowEta = (TH1D*) histoEtaCorrYieldBinShifted->Clone();
    //histoRatioFitTwoCompEta = (TH1D*) histoEtaCorrYieldBinShifted->Clone();
    
    if(kLevySuccEta) {
        histoRatioFitLevyEta = CalculateHistoRatioToFit (histoRatioFitLevyEta, fitPtLevyEta); 
    }    
    if(kHagSuccEta) {
        histoRatioFitHagEta = CalculateHistoRatioToFit (histoRatioFitHagEta, fitPtHagedornEta);
    }
    if(kPowSuccEta) {
        histoRatioFitPowEta = CalculateHistoRatioToFit (histoRatioFitPowEta, fitPtPowerlawEta);
    }
    if(kModPowSuccEta) {
        histoRatioFitModPowEta = CalculateHistoRatioToFit (histoRatioFitModPowEta, fitPtModPowerlawEta);
    }
    //if(kTwoCompSuccEta) {
    //    histoRatioFitTwoCompEta = CalculateHistoRatioToFit (histoRatioFitTwoCompEta, fitPtTwoCompModelEta);
    //}
    delete canvasFittingEta;
    
    //**********************************************************************************************
    //********************** Plotting of fitted spectrum *******************************************
    //**********************************************************************************************    
    TCanvas* canvasFittingSpectraEta = new TCanvas("canvasFittingSpectraEta","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasFittingSpectraEta, 0.13, 0.02, 0.02, 0.09);
    canvasFittingSpectraEta->SetLogy();
    
    TPad* padFittedSpectraEtaHistos = new TPad("padFittedSpectraEtaHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padFittedSpectraEtaHistos, 0.13, 0.02, 0.02, 0.);
    padFittedSpectraEtaHistos->Draw();
    
    TPad* padFittedSpectraEtaRatios = new TPad("padFittedSpectraEtaRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padFittedSpectraEtaRatios, 0.13, 0.02, 0., 0.26);
    padFittedSpectraEtaRatios->Draw();
    
    padFittedSpectraEtaHistos->cd();
    padFittedSpectraEtaHistos->SetLogy();
    
    DrawAutoGammaMesonHistos( histoEtaCorrYieldBinShifted, 
                                "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
                                kTRUE, 3., 4e-10,kTRUE, 
                                kFALSE, 3e-8,10, 
			        kFALSE, 0., 10.);
    
    DrawGammaSetMarker(histoEtaCorrYieldBinShifted, 22, 1., kBlack, kBlack);
    histoCorrectedYieldEta->Draw("e1,x0"); 
    
    fitPtLevyEta->Draw("same");
    fitPtHagedornEta->Draw("same");
    fitPtPowerlawEta->Draw("same");
    //fitPtModPowerlawEta->Draw("same");
    //fitPtTwoCompModelEta->Draw("same");
    
    TLegend* legendFitEta = new TLegend(0.15,0.02,0.65,0.15);
    legendFitEta->SetTextSize(0.02);
    legendFitEta->SetFillColor(0);
    legendFitEta->AddEntry(histoCorrectedYieldEta,"#eta");
    legendFitEta->AddEntry(fitPtLevyEta,"Levy fit");
    legendFitEta->AddEntry(fitPtHagedornEta,"Hagedorn fit");
    legendFitEta->AddEntry(fitPtPowerlawEta,"Powerlaw fit");
    //legendFitEta->AddEntry(fitPtModPowerlawEta,"ModPowerlaw fit");
    //legendFitEta->AddEntry(fitPtTwoCompModelEta, "Two component model (Bylinkin-Rostovtsev)");
    legendFitEta->Draw();
    
//     if(!thesis)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1350,1125,kDalitz);
    canvasFittingSpectraEta->Update();
    
    padFittedSpectraEtaRatios->cd();
    
    histoRatioFitLevyEta->SetYTitle("#frac{value}{fit}");
    histoRatioFitLevyEta->GetYaxis()->SetRangeUser(0.5,1.55);
    histoRatioFitLevyEta->GetYaxis()->SetLabelSize(0.08);
    histoRatioFitLevyEta->GetYaxis()->SetTitleSize(0.1);
    histoRatioFitLevyEta->GetYaxis()->SetDecimals();
    histoRatioFitLevyEta->GetYaxis()->SetTitleOffset(0.4);
    histoRatioFitLevyEta->GetXaxis()->SetTitleSize(0.11);
    histoRatioFitLevyEta->GetXaxis()->SetLabelSize(0.08);
    
    DrawGammaSetMarker(histoRatioFitLevyEta, 21, 0.5, kBlue, kBlue);
    histoRatioFitLevyEta->Draw("e1,x0");
    
    DrawGammaSetMarker(histoRatioFitHagEta, 21, 0.5, kGreen+2, kGreen+2);
    if(kHagSuccEta)histoRatioFitHagEta->Draw("e1,x0,same");
    
    DrawGammaSetMarker(histoRatioFitPowEta, 21, 0.5, kTeal, kTeal);
    if(kPowSuccEta)histoRatioFitPowEta->Draw("e1,x0,same");
    
    DrawGammaSetMarker(histoRatioFitModPowEta, 21, 0.5, kMagenta+2, kMagenta+2);
    //if(kModPowSuccEta)histoRatioFitModPowEta->Draw("e1,same");

    //DrawGammaSetMarker(histoRatioFitTwoCompEta, 21, 0.5, kRed, kRed);
    //if(kTwoCompSuccEta)histoRatioFitTwoCompEta->Draw("e1,x0,same");
    
    DrawGammaLines(0., maxPtEta,1., 1.,0.1);
    
    canvasFittingSpectraEta->SaveAs(Form("%s/Eta_%s_CorrectedYieldFitted_%s.%s",outputDir.Data(),prefix2.Data(),cutSelection.Data(),suffix));
    delete canvasFittingSpectraEta;
    cout << "here" << endl;

    //***************************** Reading systematic error and plotting Pi0 ***********************************************
    if( useSameBinningPi0Eta.CompareTo("")==0){
                    
        if (optionEnergy.CompareTo("900GeV") == 0){
            cout << "900GeV" << endl;
            graphCorrectedYieldPi0SysErr = CalculateSysErrFromRelSysHisto( histoCorrectedYieldPi0, "Pi0SystError",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0SysErrBinShifted = CalculateSysErrFromRelSysHisto( histoPi0CorrYieldBinShifted, "Pi0SystErrorBinShifted",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0SysErrA = CalculateSysErrFromRelSysHisto( histoCorrectedYieldPi0, "Pi0SystErrorA",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0SysErrABinShifted = CalculateSysErrFromRelSysHisto( histoPi0CorrYieldBinShifted, "Pi0SystErrorABinShifted",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);

            graphCorrectedYieldPi0StatPlusSys = CalculateSysErrFromRelSysHistoComplete( histoCorrectedYieldPi0, "Pi0ComplError",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0StatPlusSysBinShifted = CalculateSysErrFromRelSysHistoComplete( histoPi0CorrYieldBinShifted, "Pi0ComplErrorBinShifted",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
        
            if (optNoBinShift.CompareTo("kFALSE")==0){
                histoInvCrossSectionPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
            } else {
                histoInvCrossSectionPi0 = (TH1D*) histoCorrectedYieldPi0->Clone();
            }
            histoInvCrossSectionPi0->Scale(xSection*recalcBarn);
            graphInvCrossSectionSysPi0 = CalculateSysErrFromRelSysHisto( histoInvCrossSectionPi0 , "Pi0InvCrossSectionSys",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphInvCrossSectionSysAPi0 = CalculateSysErrFromRelSysHisto( histoInvCrossSectionPi0 , "Pi0InvCrossSectionSysA",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
            
        } else if (optionEnergy.CompareTo("2.76TeV") == 0){
            cout << "2.76TeV" << endl;
            graphCorrectedYieldPi0SysErr = CalculateSysErrFromRelSysHisto( histoCorrectedYieldPi0, "Pi0SystError",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0SysErrBinShifted = CalculateSysErrFromRelSysHisto( histoPi0CorrYieldBinShifted, "Pi0SystErrorBinShifted",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0SysErrA = CalculateSysErrFromRelSysHisto( histoCorrectedYieldPi0, "Pi0SystErrorA",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0SysErrABinShifted = CalculateSysErrFromRelSysHisto( histoPi0CorrYieldBinShifted, "Pi0SystErrorABinShifted",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
            
            graphCorrectedYieldPi0StatPlusSys = CalculateSysErrFromRelSysHistoComplete( histoCorrectedYieldPi0, "Pi0ComplError",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0StatPlusSysBinShifted = CalculateSysErrFromRelSysHistoComplete( histoPi0CorrYieldBinShifted, "Pi0ComplErrorBinShifted",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
                        
            if (optNoBinShift.CompareTo("kFALSE")==0){
                histoInvCrossSectionPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
            } else {
                histoInvCrossSectionPi0 = (TH1D*) histoCorrectedYieldPi0->Clone();
            }
            
            histoInvCrossSectionPi0->Scale(xSection*recalcBarn);
            graphInvCrossSectionSysPi0 = CalculateSysErrFromRelSysHisto( histoInvCrossSectionPi0 , "Pi0InvCrossSectionSys",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphInvCrossSectionSysAPi0 = CalculateSysErrFromRelSysHisto( histoInvCrossSectionPi0 , "Pi0InvCrossSectionSysA",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
            
        } else {cout << "else" << endl;
            graphCorrectedYieldPi0SysErr = CalculateSysErrFromRelSysHisto( histoCorrectedYieldPi0, "Pi0SystError",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0SysErrBinShifted = CalculateSysErrFromRelSysHisto( histoPi0CorrYieldBinShifted, "Pi0SystErrorBinShifted",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0SysErrA = CalculateSysErrFromRelSysHisto( histoCorrectedYieldPi0, "Pi0SystErrorA",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0SysErrABinShifted = CalculateSysErrFromRelSysHisto( histoPi0CorrYieldBinShifted, "Pi0SystErrorABinShifted",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0StatPlusSys = CalculateSysErrFromRelSysHistoComplete( histoCorrectedYieldPi0, "Pi0ComplError",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphCorrectedYieldPi0StatPlusSysBinShifted = CalculateSysErrFromRelSysHistoComplete( histoPi0CorrYieldBinShifted, "Pi0ComplErrorBinShifted",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            
            if (optNoBinShift.CompareTo("kFALSE")==0){
                histoInvCrossSectionPi0 = (TH1D*) histoPi0CorrYieldBinShifted->Clone();
            } else {
                histoInvCrossSectionPi0 = (TH1D*) histoCorrectedYieldPi0->Clone();
            }
            histoInvCrossSectionPi0->Scale(xSection*recalcBarn);
            graphInvCrossSectionSysPi0 = CalculateSysErrFromRelSysHisto( histoInvCrossSectionPi0 , "Pi0InvCrossSectionSys",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
            graphInvCrossSectionSysAPi0 = CalculateSysErrFromRelSysHisto( histoInvCrossSectionPi0 , "Pi0InvCrossSectionSysA",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
        }

        //********************************************************************************************
        //*********************** Systematic Error ***************************************************
        //********************************************************************************************
        
        //Draw Pi0 with systematic error
        TCanvas* canvasSysErrorConversion2 = new TCanvas("canvasSysErrorConversion2","",1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasSysErrorConversion2, 0.13, 0.02, 0.02, 0.09);
        canvasSysErrorConversion2->SetLogy();
        
        TH2F * histo2DSysErrorConversion = new TH2F("histo2DSysErrorConversion","histo2DSysErrorConversion",10,0.,maxPtPi0,1000.,histoCorrectedYieldPi0->GetMinimum(),5);
        histo2DSysErrorConversion->SetXTitle("p_{T} [GeV/c]");
        histo2DSysErrorConversion->SetYTitle("#frac{1}{2#pi N_{ev}} #frac{d^{2}N^{#pi^{0}}}{p_{T}dp_{T}dy} (c/GeV)^{2}");
        histo2DSysErrorConversion->SetTitle("");
        histo2DSysErrorConversion->GetYaxis()->SetDecimals();
        histo2DSysErrorConversion->GetYaxis()->SetLabelSize(0.03);
        histo2DSysErrorConversion->GetYaxis()->SetTitleSize(0.03);
        histo2DSysErrorConversion->GetYaxis()->SetTitleOffset(1.6);
        histo2DSysErrorConversion->GetXaxis()->SetTitleSize(0.03);
        histo2DSysErrorConversion->GetXaxis()->SetTitleOffset(1.);
        histo2DSysErrorConversion->GetXaxis()->SetNdivisions(510,kTRUE);
        histo2DSysErrorConversion->DrawCopy();  
        
        DrawGammaSetMarkerTGraph(graphCorrectedYieldPi0SysErr, 20, 0.2, kBlack, kBlack);
        graphCorrectedYieldPi0SysErr->SetFillColor(kGray+1);
        graphCorrectedYieldPi0SysErr->Draw("same,2,p");
        DrawGammaSetMarker(histoCorrectedYieldPi0, 20, 0.5, kBlack, kBlack);
        histoCorrectedYieldPi0->Draw("same,e1,p");
        
        histo2DSysErrorConversion->DrawCopy("same");  
        
//         if(!thesis)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], 0.03, pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], 0.04, pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1350,900,kDalitz);
        
        canvasSysErrorConversion2->Update();
        canvasSysErrorConversion2->SaveAs(Form("%s/Pi0_%s_CorrectedYieldwithSysErrorConversion_%s_%s.%s",outputDir.Data(),prefix2.Data(),cutSelection.Data(),useSameBinningPi0Eta.Data(),suffix));
        delete canvasSysErrorConversion2;
        
        
        //****************** Pi0 corrected including systematic error with binshift ***************************************
        TCanvas* canvasSysErrorConversion3 = new TCanvas("canvasSysErrorConversion3","",1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasSysErrorConversion3, 0.13, 0.02, 0.02, 0.09);
        canvasSysErrorConversion3->SetLogy();
        
        TH2F * histo2DSysErrorConversionBinShifted = new TH2F("histo2DSysErrorConversionBinShifted","histo2DSysErrorConversionBinShifted",10,0.,maxPtPi0,1000.,histoCorrectedYieldPi0->GetMinimum(),5);
        histo2DSysErrorConversionBinShifted->SetXTitle("p_{T} [GeV/c]");
        histo2DSysErrorConversionBinShifted->SetYTitle("#frac{1}{2#pi N_{ev}} #frac{d^{2}N^{#pi^{0}}}{p_{T}dp_{T}dy} (c/GeV)^{2}");
        histo2DSysErrorConversionBinShifted->SetTitle("");
        histo2DSysErrorConversionBinShifted->GetYaxis()->SetDecimals();
        histo2DSysErrorConversionBinShifted->GetYaxis()->SetLabelSize(0.03);
        histo2DSysErrorConversionBinShifted->GetYaxis()->SetTitleSize(0.03);
        histo2DSysErrorConversionBinShifted->GetYaxis()->SetTitleOffset(1.6);
        histo2DSysErrorConversionBinShifted->GetXaxis()->SetTitleSize(0.03);
        histo2DSysErrorConversionBinShifted->GetXaxis()->SetTitleOffset(1.);
        histo2DSysErrorConversionBinShifted->GetXaxis()->SetNdivisions(510,kTRUE);
        histo2DSysErrorConversionBinShifted->DrawCopy();  
        
        DrawGammaSetMarkerTGraph(graphCorrectedYieldPi0SysErrBinShifted, 20, 0.2, kBlack, kBlack);
        graphCorrectedYieldPi0SysErrBinShifted->SetFillColor(kGray+1);
        graphCorrectedYieldPi0SysErrBinShifted->Draw("same,2,p");
        DrawGammaSetMarker(histoPi0CorrYieldBinShifted, 20, 0.5, kBlack, kBlack);
        histoPi0CorrYieldBinShifted->Draw("same,e1,x0,p");
        
        histo2DSysErrorConversionBinShifted->DrawCopy("same");  
                
//         if(!thesis)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], 0.03, pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], 0.04, pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1350,900,kDalitz);
        
        canvasSysErrorConversion3->Update();
        canvasSysErrorConversion3->SaveAs(Form("%s/Pi0_%s_CorrectedYieldwithSysErrorConversionBinShifted_%s_%s.%s", outputDir.Data(), prefix2.Data(), cutSelection.Data(), useSameBinningPi0Eta.Data(), suffix));
        delete canvasSysErrorConversion3;
        
        
        
        
        //***************************************************************************************************
        //*************************** Fitting  Pi0 Spectrum *************************************************
        //***************************************************************************************************
        
        fileFinalResults << "Final Results for Pi0 with Systematic Errors" << endl << endl;
        
        //****************************************************************************************************
        //************************** Fitting of corrected normal spectrum ************************************
        //****************************************************************************************************
        
        TCanvas* canvasFittingPi0SysErr = new TCanvas("canvasFittingPi0SysErr","",200,10,1350,900);  // gives the page size
        
        graphCorrectedYieldPi0SysErrForFit = (TGraphAsymmErrors*) graphCorrectedYieldPi0StatPlusSysBinShifted->Clone();
        
        //****************************** Fit histCorr with Levy *****************************************
        fitPtLevyPi0SysErr = FitObject("l","fitPtLevyPi0SysErr","Pi0",graphCorrectedYieldPi0SysErrForFit,minPtForFits,maxPtPi0,NULL,"QNRME+");
        DrawGammaSetMarkerTF1(fitPtLevyPi0SysErr, 1, 1.5, kBlue);
        forOutput= WriteParameterToFile(fitPtLevyPi0SysErr);
        fileFinalResults<< forOutput.Data()<< endl;
        
        //**************************** Fit histCorr with Hagedorn  **********************************
        fitPtHagedornPi0SysErr = FitObject("h","fitPtHagedornPi0SysErr","Pi0",graphCorrectedYieldPi0SysErrForFit,minPtForFits,maxPtPi0,NULL,"QNRME+");
        DrawGammaSetMarkerTF1( fitPtHagedornPi0SysErr, 1, 1.5, kGreen+2);
        fitPtHagedornPi0SysErr->SetLineStyle(7);
        forOutput= WriteParameterToFile(fitPtHagedornPi0SysErr);
        fileFinalResults<< forOutput.Data()<< endl;
        
        //************************** Fit histCorr with Boltzmann  ***********************************
        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fitPtBoltzmannPi0SysErr = FitObject("b","fitPtBoltzmannPi0SysErr","Pi0",graphCorrectedYieldPi0SysErrForFit,0.3,2.,NULL,"QNRME+");
            DrawGammaSetMarkerTF1( fitPtBoltzmannPi0SysErr, 1, 1.5, kMagenta+4);
            forOutput= WriteParameterToFile(fitPtBoltzmannPi0SysErr);
            fileFinalResults<< forOutput.Data()<< endl;
        }    
        
        //*************************** Fit histCorr with Exponential **********************************
        if (useSameBinningPi0Eta.CompareTo("")==0){
            fitPtExpPi0SysErr = FitObject("e","fitPtExpPi0SysErr","Pi0",graphCorrectedYieldPi0SysErrForFit,0.3,2.,NULL,"QNRME+");
            DrawGammaSetMarkerTF1( fitPtExpPi0SysErr, 1, 1.5, kOrange+7);
            forOutput= WriteParameterToFile(fitPtExpPi0SysErr);
            fileFinalResults<< forOutput.Data()<< endl;
        }
        
        //**************************** Fit histCorr with Powerlaw  *********************************
        fitPtPowerlawPi0SysErr = FitObject("p","fitPtPowerlawPi0SysErr","Pi0",graphCorrectedYieldPi0SysErrForFit,1.5,maxPtPi0,NULL,"QNRME+");
        DrawGammaSetMarkerTF1( fitPtPowerlawPi0SysErr, 1, 1.5, kTeal);
        forOutput= WriteParameterToFile(fitPtPowerlawPi0SysErr);
        fileFinalResults<< forOutput.Data()<< endl;
        
        //**************************** Fit histCorr with ModPowerlaw  *********************************
        fitPtModPowerlawPi0SysErr = FitObject("m","fitPtModPowerlawPi0SysErr","Pi0",graphCorrectedYieldPi0SysErrForFit,0.3,maxPtPi0,NULL,"QNRME+");
        DrawGammaSetMarkerTF1( fitPtModPowerlawPi0SysErr, 1, 1.5, kMagenta+2);
        forOutput= WriteParameterToFile(fitPtModPowerlawPi0SysErr);
        fileFinalResults<< forOutput.Data()<< endl;
        
        delete canvasFittingPi0SysErr;
        
//         *************************************** Systematic Error of the Eta *************************************************
        graphCorrectedYieldEtaSysErr = CalculateSysErrFromRelSysHisto( histoCorrectedYieldEta, "EtaSystError",relSystErrorEtaDown , relSystErrorEtaUp, 2, nPointsEta);
        graphCorrectedYieldEtaSysErrBinShifted = CalculateSysErrFromRelSysHisto( histoEtaCorrYieldBinShifted, "EtaSystErrorBinShifted",relSystErrorEtaDown , relSystErrorEtaUp, 2, nPointsEta);
        graphCorrectedYieldEtaSysErrA = CalculateSysErrFromRelSysHisto( histoCorrectedYieldEta, "EtaSystErrorA",relSystErrorWOMaterialEtaDown , relSystErrorWOMaterialEtaUp, 2, nPointsEta);
        graphCorrectedYieldEtaSysErrABinShifted = CalculateSysErrFromRelSysHisto( histoEtaCorrYieldBinShifted, "EtaSystErrorABinShifted",relSystErrorWOMaterialEtaDown , relSystErrorWOMaterialEtaUp, 2, nPointsEta);
        graphCorrectedYieldEtaStatPlusSys = CalculateSysErrFromRelSysHistoComplete( histoCorrectedYieldEta, "EtaComplError",relSystErrorEtaDown, relSystErrorEtaUp, 2, nPointsEta);
        graphCorrectedYieldEtaStatPlusSysBinShifted = CalculateSysErrFromRelSysHistoComplete( histoEtaCorrYieldBinShifted , "EtaComplErrorBinShifted",relSystErrorEtaDown , relSystErrorEtaUp, 2, nPointsEta);
        
        if (optNoBinShift.CompareTo("kFALSE")==0){
            histoInvCrossSectionEta = (TH1D*) histoEtaCorrYieldBinShifted->Clone();
        } else {
            histoInvCrossSectionEta = (TH1D*) histoCorrectedYieldEta->Clone();
        }
        histoInvCrossSectionEta->Scale(xSection*recalcBarn);
        graphInvCrossSectionSysEta = CalculateSysErrFromRelSysHisto( histoInvCrossSectionEta , "Pi0InvCrossSectionSys",relSystErrorEtaDown , relSystErrorEtaUp, 2, nPointsEta);
        graphInvCrossSectionSysAEta = CalculateSysErrFromRelSysHisto( histoInvCrossSectionEta , "Pi0InvCrossSectionSysA",relSystErrorWOMaterialEtaDown , relSystErrorWOMaterialEtaUp, 2, nPointsEta);
        
//         ********************************************************************************************
//         *********************** Systematic Error ***************************************************
//         ********************************************************************************************
        
        //Draw Eta with systematic error
        TCanvas* canvasSysErrorConversion2Eta = new TCanvas("canvasSysErrorConversion2Eta","",1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasSysErrorConversion2Eta, 0.13, 0.02, 0.02, 0.09);
        canvasSysErrorConversion2Eta->SetLogy();
        
        TH2F * histo2DSysErrorConversionEta = new TH2F("histo2DSysErrorConversionEta","histo2DSysErrorConversionEta",10,0.,maxPtEta,1000.,histoCorrectedYieldEta->GetMinimum(),5);
        histo2DSysErrorConversionEta->SetXTitle("p_{T} [GeV/c]");
        histo2DSysErrorConversionEta->SetYTitle("#frac{1}{2#pi N_{ev}} #frac{d^{2}N^{#pi^{0}}}{p_{T}dp_{T}dy} (c/GeV)^{2}");
        histo2DSysErrorConversionEta->SetTitle("");
        histo2DSysErrorConversionEta->GetYaxis()->SetDecimals();
        histo2DSysErrorConversionEta->GetYaxis()->SetLabelSize(0.03);
        histo2DSysErrorConversionEta->GetYaxis()->SetTitleSize(0.03);
        histo2DSysErrorConversionEta->GetYaxis()->SetTitleOffset(1.6);
        histo2DSysErrorConversionEta->GetXaxis()->SetTitleSize(0.03);
        histo2DSysErrorConversionEta->GetXaxis()->SetTitleOffset(1.);
        histo2DSysErrorConversionEta->GetXaxis()->SetNdivisions(510,kTRUE);
        histo2DSysErrorConversionEta->DrawCopy();  
        
        DrawGammaSetMarkerTGraph(graphCorrectedYieldEtaSysErr, 20, 0.2, kBlack, kBlack);
        graphCorrectedYieldEtaSysErr->SetFillColor(kGray+1);
        graphCorrectedYieldEtaSysErr->Draw("same,2,p");
        DrawGammaSetMarker(histoCorrectedYieldEta, 20, 0.5, kBlack, kBlack);
        histoCorrectedYieldEta->Draw("same,e1,p");
        
        histo2DSysErrorConversionEta->DrawCopy("same");  
        
        TLegend* legendYieldSysEta = new TLegend(0.15,0.1,0.5,0.19);
        legendYieldSysEta->SetTextSize(0.02);
        legendYieldSysEta->SetFillColor(0);
        legendYieldSysEta->AddEntry(histoCorrectedYieldEta,"Corrected yield");
        legendYieldSysEta->AddEntry(graphCorrectedYieldEtaSysErr,"Systematic uncertainty","f");
        legendYieldSysEta->Draw();
        
//         if(!thesis)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], 0.03, pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], 0.04, pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1350,900,kDalitz);
        
        canvasSysErrorConversion2Eta->Update();
        canvasSysErrorConversion2Eta->SaveAs(Form("%s/Eta_%s_CorrectedYieldwithSysErrorConversion_%s_%s.%s",outputDir.Data(),prefix2.Data(),cutSelection.Data(),useSameBinningPi0Eta.Data(),suffix));
        delete canvasSysErrorConversion2Eta;
        
        
        //****************** Eta corrected including systematic error with binshift ***************************************
        TCanvas* canvasSysErrorConversion3Eta = new TCanvas("canvasSysErrorConversion3Eta","",1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasSysErrorConversion3Eta, 0.13, 0.02, 0.02, 0.09);
        canvasSysErrorConversion3Eta->SetLogy();
        
        TH2F * histo2DSysErrorConversionEtaBinShifted = new TH2F("histo2DSysErrorConversionEtaBinShifted","histo2DSysErrorConversionEtaBinShifted",10,0.,maxPtEta,1000.,histoCorrectedYieldEta->GetMinimum(),5);
        histo2DSysErrorConversionEtaBinShifted->SetXTitle("p_{T} [GeV/c]");
        histo2DSysErrorConversionEtaBinShifted->SetYTitle("#frac{1}{2#pi N_{ev}} #frac{d^{2}N^{#pi^{0}}}{p_{T}dp_{T}dy} (c/GeV)^{2}");
        histo2DSysErrorConversionEtaBinShifted->SetTitle("");
        histo2DSysErrorConversionEtaBinShifted->GetYaxis()->SetDecimals();
        histo2DSysErrorConversionEtaBinShifted->GetYaxis()->SetLabelSize(0.03);
        histo2DSysErrorConversionEtaBinShifted->GetYaxis()->SetTitleSize(0.03);
        histo2DSysErrorConversionEtaBinShifted->GetYaxis()->SetTitleOffset(1.6);
        histo2DSysErrorConversionEtaBinShifted->GetXaxis()->SetTitleSize(0.03);
        histo2DSysErrorConversionEtaBinShifted->GetXaxis()->SetTitleOffset(1.);
        histo2DSysErrorConversionEtaBinShifted->GetXaxis()->SetNdivisions(510,kTRUE);
        histo2DSysErrorConversionEtaBinShifted->DrawCopy();  
        
        DrawGammaSetMarkerTGraph(graphCorrectedYieldEtaSysErrBinShifted, 20, 0.2, kBlack, kBlack);
        graphCorrectedYieldEtaSysErrBinShifted->SetFillColor(kGray+1);
        graphCorrectedYieldEtaSysErrBinShifted->Draw("same,2,p");
        DrawGammaSetMarker(histoEtaCorrYieldBinShifted, 20, 0.5, kBlack, kBlack);
        histoEtaCorrYieldBinShifted->Draw("same,e1,x0,p");
        
        histo2DSysErrorConversionEtaBinShifted->DrawCopy("same");  
        
        TLegend* legendYieldSysEtaBinShifted = new TLegend(0.15,0.1,0.5,0.19);
        legendYieldSysEtaBinShifted->SetTextSize(0.02);
        legendYieldSysEtaBinShifted->SetFillColor(0);
        legendYieldSysEtaBinShifted->AddEntry(histoEtaCorrYieldBinShifted,"Corrected yield");
        legendYieldSysEtaBinShifted->AddEntry(graphCorrectedYieldEtaSysErrBinShifted,"Systematic uncertainty","f");
        legendYieldSysEtaBinShifted->Draw();
        
//         if(!thesis)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], 0.03, pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], 0.04, pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1350,900,kDalitz);
        
        canvasSysErrorConversion3Eta->Update();
        canvasSysErrorConversion3Eta->SaveAs(Form("%s/Eta_%s_CorrectedYieldwithSysErrorConversionBinShifted_%s_%s.%s",outputDir.Data(),prefix2.Data(),cutSelection.Data(),useSameBinningPi0Eta.Data(),suffix));
        delete canvasSysErrorConversion3Eta;
        
        //********************************************************************************************
        //*********************** Systematic Error ***************************************************
        //********************************************************************************************
        
        
        fileFinalResults << "Final Results for Pi0 with Systematic Errors" << endl << endl;
        
        //****************************************************************************************************
        //************************** Fitting of corrected normal spectrum ************************************
        //****************************************************************************************************
        
        TCanvas* canvasFittingEtaSysErr = new TCanvas("canvasFittingEtaSysErr","",200,10,1350,900);  // gives the page size
        
        graphCorrectedYieldEtaSysErrForFit = (TGraphAsymmErrors*) graphCorrectedYieldEtaStatPlusSysBinShifted->Clone();
        
        //****************************** Fit histCorr with Levy *****************************************
        fitPtLevyEtaSysErr = FitObject("l","fitPtLevyEtaSysErr","Eta",graphCorrectedYieldEtaSysErrForFit,minPtForFits,maxPtEta,NULL,"QNRME+");
        DrawGammaSetMarkerTF1(fitPtLevyEtaSysErr, 1, 1.5, kBlue);
        forOutput= WriteParameterToFile(fitPtLevyEtaSysErr);
        fileFinalResults<< forOutput.Data()<< endl;
        
        //**************************** Fit histCorr with Hagedorn  **********************************
        fitPtHagedornEtaSysErr = FitObject("h","fitPtHagedornEtaSysErr","Eta",graphCorrectedYieldEtaSysErrForFit,minPtForFits,maxPtEta,NULL,"QNRME+");
        DrawGammaSetMarkerTF1( fitPtHagedornEtaSysErr, 1, 1.5, kGreen+2);
        fitPtHagedornEtaSysErr->SetLineStyle(7);
        forOutput= WriteParameterToFile(fitPtHagedornEtaSysErr);
        fileFinalResults<< forOutput.Data()<< endl;
        
        //**************************** Fit histCorr with Powerlaw  *********************************
        fitPtPowerlawEtaSysErr = FitObject("p","fitPtPowerlawEtaSysErr","Eta",graphCorrectedYieldEtaSysErrForFit,1.5,maxPtEta,NULL,"QNRME+");
        DrawGammaSetMarkerTF1( fitPtPowerlawEtaSysErr, 1, 1.5, kTeal);
        forOutput= WriteParameterToFile(fitPtPowerlawEtaSysErr);
        fileFinalResults<< forOutput.Data()<< endl;
        
        //**************************** Fit histCorr with ModPowerlaw  *********************************
        fitPtModPowerlawEtaSysErr = FitObject("m","fitPtModPowerlawEtaSysErr","Eta",graphCorrectedYieldEtaSysErrForFit,0.3,maxPtEta,NULL,"QNRME+");
        DrawGammaSetMarkerTF1( fitPtModPowerlawEtaSysErr, 1, 1.5, kMagenta+2);
        forOutput= WriteParameterToFile(fitPtModPowerlawEtaSysErr);
        fileFinalResults<< forOutput.Data()<< endl;
        
        delete canvasFittingEtaSysErr;
     }

    //************************************************************************************************
    //******************************** Combined Plot Pi0 and Eta *************************************
    //************************************************************************************************
    if( useSameBinningPi0Eta.CompareTo("")!=0){
        canvasCorrectedSame = new TCanvas("canvasCorrectedSame","",1350,1500);  // gives the page size
        DrawGammaCanvasSettings( canvasCorrectedSame, 0.13, 0.02, 0.02, 0.09);
        canvasCorrectedSame->SetLogy();
        
        TPad* padCorrectedHistos = new TPad("padCorrectedHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padCorrectedHistos, 0.12, 0.02, 0.02, 0.);
        padCorrectedHistos->Draw();
        
        TPad* padCorrectedRatios = new TPad("padCorrectedRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padCorrectedRatios, 0.12, 0.02, 0., 0.26);
        padCorrectedRatios->Draw();
        
        padCorrectedHistos->cd();
        padCorrectedHistos->SetLogy();
        
        DrawAutoGammaMesonHistos( histoCorrectedYieldPi0, 
                            "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
                            kTRUE, 3., 4e-10, kTRUE,
                            kFALSE, 3e-8,10, 
                            kFALSE, 0., 10.);
        
        DrawGammaSetMarker(histoCorrectedYieldPi0, 22, 1., kBlack, kBlack);
        histoCorrectedYieldPi0->DrawCopy("e1,x0"); 
        
        DrawGammaSetMarker(histoCorrectedYieldEta, 26, 1., kAzure, kAzure);
        histoCorrectedYieldEta->DrawCopy("e1,same"); 
        TLegend* legendYield = new TLegend(0.15,0.02,0.5,0.12);
        legendYield->SetTextSize(0.02);
        legendYield->SetFillColor(0);
        legendYield->AddEntry(histoCorrectedYieldPi0,"#pi^{0}");
        legendYield->AddEntry(histoCorrectedYieldEta,"#eta");
        legendYield->Draw();
        
//         if(!thesis)DrawAliceLogoCombinedWorkInProgress(pictDrawingCoordinatesRat[0], pictDrawingCoordinatesRat[1], pictDrawingCoordinatesRat[2], pictDrawingCoordinatesRat[3], pictDrawingCoordinatesRat[4], pictDrawingCoordinatesRat[5], pictDrawingCoordinatesRat[6], pictDrawingCoordinatesRat[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,1125,kDalitz);

        padCorrectedRatios->cd();
        padCorrectedRatios->SetLogy(0);
        
        histoRatioEtaPi0 = (TH1D*) histoCorrectedYieldEta->Clone();
        histoRatioEtaPi0->Divide(histoRatioEtaPi0,histoCorrectedYieldPi0,1.,1.,"");

        histoRatioEtaPi0BinShifted = (TH1D*) histoEtaCorrYieldBinShifted->Clone();
        histoRatioEtaPi0BinShifted->Divide(histoEtaCorrYieldBinShifted,histoPi0CorrYieldBinShifted,1.,1.,"");

        DrawAutoGammaMesonHistos( histoRatioEtaPi0, 
                             "", "p_{T} (GeV/c)", "#frac{#eta}{#pi^{0}}", 
                             kFALSE, 3., 4e-10, kTRUE,
                             kTRUE, 0.0,1.05,
                             kFALSE, 0., 10.);
        histoRatioEtaPi0->GetYaxis()->SetNdivisions(505);
        histoRatioEtaPi0->GetYaxis()->SetLabelSize(0.08);
        histoRatioEtaPi0->GetYaxis()->SetTitleSize(0.1);
        histoRatioEtaPi0->GetYaxis()->SetDecimals();
        histoRatioEtaPi0->GetYaxis()->SetTitleOffset(0.45);
        histoRatioEtaPi0->GetXaxis()->SetTitleSize(0.11);
        histoRatioEtaPi0->GetXaxis()->SetLabelSize(0.08);
        DrawGammaSetMarker(histoRatioEtaPi0, 22, 1., kBlack, kBlack);
        histoRatioEtaPi0->DrawCopy("e1,x0"); 
        
        if(optionEnergy.CompareTo("8TeV") != 0){
            histoEtaToPi0Phojet->SetLineColor(kRed+2);
            histoEtaToPi0Phojet->DrawCopy("same,x0,C");
        
            histoEtaToPi0Pythia->SetLineColor(kRed-2);
            histoEtaToPi0Pythia->DrawCopy("same,x0,C");
        }
        
//         DrawGammaLines(0., maxPtEta,0.45, 0.45,0.1);
        
        canvasCorrectedSame->Update();
        canvasCorrectedSame->SaveAs(Form("%s/%s_CombinedCorrectedYield_%s_%s.%s",outputDir.Data(),prefix2.Data(),useSameBinningPi0Eta.Data(),cutSelection.Data(),suffix));
        
        delete canvasCorrectedSame;

                             
        //***************************** Only ratio *************************************************
        TCanvas* canvasRatioEtaPi0 = new TCanvas("canvasRatioEtaPi0","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioEtaPi0, 0.08, 0.01, 0.015, 0.09);
        
        
        DrawAutoGammaMesonHistos( histoRatioEtaPi0, 
                            "", "p_{T} (GeV/c)", "#frac{#eta}{#pi^{0}}", 
                            kFALSE, 3., 4e-10, kTRUE,
                            kTRUE, 0.0,1.0,
                            kFALSE, 0., 10.);
        histoRatioEtaPi0->GetYaxis()->SetLabelSize(0.03);
        histoRatioEtaPi0->GetYaxis()->SetTitleSize(0.04);
        histoRatioEtaPi0->GetYaxis()->SetDecimals();
        histoRatioEtaPi0->GetYaxis()->SetTitleOffset(0.9);
        histoRatioEtaPi0->GetXaxis()->SetTitleOffset(1.);
        histoRatioEtaPi0->GetXaxis()->SetLabelSize(0.03);
        histoRatioEtaPi0->GetXaxis()->SetTitleSize(0.04);
        histoRatioEtaPi0->GetYaxis()->SetNdivisions(510);
        DrawGammaSetMarker(histoRatioEtaPi0, 20, 1., kBlack, kBlack);
        histoRatioEtaPi0->DrawCopy("e1,x0"); 
        
//         DrawGammaLines(0., maxPtEta ,0.45, 0.45,0.1);

        cout << " calculation of eta/pi0" << endl;
        for (Int_t i = 0; i < nPointsEta + 1; i++){
            ratioXValue[i] = histoRatioEtaPi0->GetBinCenter(i+2);
            ratioYValue[i] = histoRatioEtaPi0->GetBinContent(i+2);
            ratioXError[i] = histoRatioEtaPi0->GetBinWidth(i+2)/2.;

            ratioSysUpError[i]= TMath::Sqrt(TMath::Power(relSystErrorWOMaterialEtaUp[i]/100* histoCorrectedYieldEta->GetBinContent(i+2)/histoCorrectedYieldPi0->GetBinContent(i+2),2) + TMath::Power(histoCorrectedYieldEta->GetBinContent(i+2)*relSystErrorWOMaterialPi0Up[i]/100/histoCorrectedYieldPi0->GetBinContent(i+2),2));
            ratioSysDownError[i]= TMath::Sqrt(TMath::Power(relSystErrorWOMaterialEtaDown[i]/100* histoCorrectedYieldEta->GetBinContent(i+2)/histoCorrectedYieldPi0->GetBinContent(i+2),2) + TMath::Power(histoCorrectedYieldEta->GetBinContent(i+2)*relSystErrorWOMaterialPi0Down[i]/100/histoCorrectedYieldPi0->GetBinContent(i+2),2));

            ratioYValueBinShifted[i] = histoRatioEtaPi0BinShifted->GetBinContent(i+2);
            ratioSysUpErrorBinShifted[i]= TMath::Sqrt(TMath::Power(relSystErrorWOMaterialEtaUp[i]/100* histoEtaCorrYieldBinShifted->GetBinContent(i+2)/histoPi0CorrYieldBinShifted->GetBinContent(i+2),2) + TMath::Power(histoEtaCorrYieldBinShifted->GetBinContent(i+2)*relSystErrorWOMaterialPi0Up[i]/100/histoPi0CorrYieldBinShifted->GetBinContent(i+2),2));
            ratioSysDownErrorBinShifted[i]= TMath::Sqrt(TMath::Power(relSystErrorWOMaterialEtaDown[i]/100* histoEtaCorrYieldBinShifted->GetBinContent(i+2)/histoPi0CorrYieldBinShifted->GetBinContent(i+2),2) + TMath::Power(histoEtaCorrYieldBinShifted->GetBinContent(i+2)*relSystErrorWOMaterialPi0Down[i]/100/histoPi0CorrYieldBinShifted->GetBinContent(i+2),2));
            cout << ratioXValue[i] << "\t" << ratioYValue[i] << "\t" << ratioSysUpError[i] << "\t" << ratioSysDownError[i] << endl;
        }

        graphSystErrRatio = new TGraphAsymmErrors(nPointsEta,ratioXValue,ratioYValue,ratioXError,ratioXError,ratioSysDownError,ratioSysUpError);
        graphSystErrRatio->SetFillColor(kGray+1);

        graphSystErrRatioBinShifted = new TGraphAsymmErrors(nPointsEta,ratioXValue,ratioYValueBinShifted,ratioXError,ratioXError,ratioSysDownErrorBinShifted,ratioSysUpErrorBinShifted);
        graphSystErrRatioBinShifted->SetFillColor(kGray+1);
        
        graphSystErrRatio->Draw("p,2,same");
        histoRatioEtaPi0->DrawCopy("same,x0,e1"); 
                
        TLegend* legendRatio = new TLegend(0.6,0.12,0.97,0.22);
        legendRatio->SetTextSize(0.03);
        legendRatio->SetFillColor(0);
        legendRatio->AddEntry(histoRatioEtaPi0,Form("p+p ALICE ( #sqrt{#it{s}}= %s )",optionEnergy.Data()),"p");
        legendRatio->AddEntry(graphSystErrRatio,"systematic uncertainty","f");
        legendRatio->Draw();

//         if(!thesis)DrawAliceLogoCombinedWorkInProgress(pictDrawingCoordinatesPi0Eta[0], pictDrawingCoordinatesPi0Eta[1], pictDrawingCoordinatesPi0Eta[2], pictDrawingCoordinatesPi0Eta[3], pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6], pictDrawingCoordinatesPi0Eta[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,kDalitz);

        canvasRatioEtaPi0->Update();
        canvasRatioEtaPi0->SaveAs(Form("%s/%s_Pi0EtaRatio_%s_%s.%s",outputDir.Data(), prefix2.Data(), useSameBinningPi0Eta.Data(), cutSelection.Data(), suffix));

        delete canvasRatioEtaPi0;
        delete legendRatio;

        //******************************* Ratio + Pythia ********************************************

        TCanvas* canvasRatioEtaPi02 = new TCanvas("canvasRatioEtaPi02","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioEtaPi02, 0.08, 0.01, 0.015, 0.09);
        
        DrawAutoGammaMesonHistos( histoRatioEtaPi0, 
                            "", "p_{T} (GeV/c)", "#frac{#eta}{#pi^{0}}", 
                            kFALSE, 3., 4e-10, kTRUE,
                            kTRUE, 0.0,1.0,
                            kFALSE, 0., 10.);
        histoRatioEtaPi0->GetYaxis()->SetLabelSize(0.03);
        histoRatioEtaPi0->GetYaxis()->SetTitleSize(0.04);
        histoRatioEtaPi0->GetYaxis()->SetDecimals();
        histoRatioEtaPi0->GetYaxis()->SetTitleOffset(0.9);
        histoRatioEtaPi0->GetXaxis()->SetTitleOffset(1.);
        histoRatioEtaPi0->GetXaxis()->SetLabelSize(0.03);
        histoRatioEtaPi0->GetXaxis()->SetTitleSize(0.04);
        histoRatioEtaPi0->GetYaxis()->SetNdivisions(510);
        DrawGammaSetMarker(histoRatioEtaPi0, 20, 1., kBlack, kBlack);
        histoRatioEtaPi0->DrawCopy("e1"); 
//         DrawGammaSetMarker(histoRatioEtaPi0BinShifted, 20, 1., kBlue, kBlue);
//         histoRatioEtaPi0BinShifted->DrawCopy("e1"); 
        
        if(optionEnergy.CompareTo("8TeV") != 0){
            histoEtaToPi0Phojet->SetLineColor(kRed+2);
            histoEtaToPi0Phojet->Draw("same,e1,x0");
            histoEtaToPi0Pythia->SetLineColor(kRed-2);
            histoEtaToPi0Pythia->Draw("same,e1,x0");
        }
        
//         DrawGammaLines(0., maxPtEta ,0.45, 0.45,0.1);
        
        graphSystErrRatio->Draw("p,2,same");
        
        histoRatioEtaPi0->DrawCopy("same,axis,e1,x0"); 
        
        TLegend* legendRatio2 = new TLegend(0.5,0.1,0.95,0.2);
        legendRatio2->SetTextSize(0.02);
        legendRatio2->SetFillColor(0);
        legendRatio2->AddEntry(histoRatioEtaPi0,Form("p+p ALICE (#sqrt{#it{s}} = %s)",optionEnergy.Data()));
        if(optionEnergy.CompareTo("8TeV") != 0){
            legendRatio2->AddEntry(histoEtaToPi0Phojet,"Phojet");
            legendRatio2->AddEntry(histoEtaToPi0Pythia,"Pythia");
        }
        legendRatio2->Draw();
        
//         if(!thesis)DrawAliceLogoCombinedWorkInProgress(pictDrawingCoordinatesPi0Eta[0], pictDrawingCoordinatesPi0Eta[1], pictDrawingCoordinatesPi0Eta[2], pictDrawingCoordinatesPi0Eta[3], pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6], pictDrawingCoordinatesPi0Eta[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,kDalitz);
        
        canvasRatioEtaPi02->Update();
        canvasRatioEtaPi02->SaveAs(Form("%s/%s_Pi0EtaRatioPythia_%s_%s.%s",outputDir.Data(), prefix2.Data(), useSameBinningPi0Eta.Data(), cutSelection.Data(), suffix));
        
        delete canvasRatioEtaPi02;
        delete legendRatio2;
                
    } else {     

        canvasCorrectedSame = new TCanvas("canvasCorrectedSame","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasCorrectedSame, 0.13, 0.02, 0.02, 0.09);
        canvasCorrectedSame->SetLogy();

        DrawAutoGammaMesonHistos( histoCorrectedYieldPi0, 
                            "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
                            kTRUE, 3., 4e-10, kTRUE,
                            kFALSE, 3e-8,10, 
                            kFALSE, 0., 10.);

        DrawGammaSetMarker(histoCorrectedYieldPi0, 22, 1., kBlack, kBlack);
        histoCorrectedYieldPi0->DrawCopy("e1,x0"); 

        DrawGammaSetMarker(histoCorrectedYieldEta, 26, 1., kAzure, kAzure);
        histoCorrectedYieldEta->DrawCopy("e1,x0,same"); 

        TLegend* legendYield = new TLegend(0.15,0.13,0.5,0.22);
        legendYield->SetTextSize(0.02);
        legendYield->SetFillColor(0);
        legendYield->AddEntry(histoCorrectedYieldPi0,Form("#pi^{0}"));
        legendYield->AddEntry(histoCorrectedYieldEta,"#eta");
        legendYield->Draw();

//         if(!thesis)DrawAliceLogoCombinedWorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,kDalitz);
        canvasCorrectedSame->Update();

        canvasCorrectedSame->SaveAs(Form("%s/%s_CombinedCorrectedYield_%s.%s",outputDir.Data(),prefix2.Data(),cutSelection.Data(),suffix));

        delete canvasCorrectedSame;
    }
    cout << "skiped eta" << endl;
    fileFinalResults.close();

    
    //*********************************************************************************************************
    //********************** ComparisonFile Output ************************************************************
    //*********************************************************************************************************
    TString applyBinShift = "";
    if( optNoBinShift.CompareTo("kTRUE")==0){
        applyBinShift = "NoBinShifting";
        cout << applyBinShift << endl;
    }
    
    TString system = "PCM";
    if (mode == 2) system     = "PCM-EMCAL";
    if (mode == 3) system     = "PCM-PHOS";
    if (mode == 4) system     = "EMCAL-EMCAL";
    if (mode == 5) system     = "PHOS-PHOS";

    
    const char* fileNameOutputComp = Form("%s/%s_%sResultsFullCorrection_PP_%s.root",cutSelection.Data(),prefix2.Data(),system.Data(), applyBinShift.Data());
    TFile* fileOutputForComparisonFullyCorrected = new TFile(fileNameOutputComp,"UPDATE");
        cout << "done" << endl;
        cout << "here" << endl;
        if (useSameBinningPi0Eta.CompareTo("")==0){
            histoNumberOfEvents->Write(Form("histoNumberOfEvents%s",optionEnergy.Data()),TObject::kOverwrite);
            cout << "here" << endl;
        }
        cout << "here" << endl;
        fileOutputForComparisonFullyCorrected->mkdir(Form("Pi0%s%s",optDalitz.Data(),optionEnergy.Data()));
        cout << "here" << endl;
        TDirectoryFile* directoryPi0 = (TDirectoryFile*)fileOutputForComparisonFullyCorrected->Get(Form("Pi0%s%s",optDalitz.Data(),optionEnergy.Data())); 
        cout << "here" << endl;
        fileOutputForComparisonFullyCorrected->cd(Form("Pi0%s%s",optDalitz.Data(),optionEnergy.Data()));
        cout << "done" << endl;
        if (useSameBinningPi0Eta.CompareTo("")==0){
            histoCorrectedYieldPi0->Write("CorrectedYieldPi0",TObject::kOverwrite);
            cout << "here" << endl;
            histoUncorrectedYieldPi0->Write("RAWYieldPerEventsPi0",TObject::kOverwrite);
            cout << "here" << endl;
            histoAccPi0->Write("AcceptancePi0",TObject::kOverwrite);
            cout << "here" << endl;
            histoTrueEffPtPi0->Write("EfficiencyPi0",TObject::kOverwrite);
            cout << "here" << endl;
            histoFWHMMesonPi0->Write("FWHMPi0",TObject::kOverwrite);
            cout << "here" << endl;
            histoMassMesonPi0->Write("MassPi0",TObject::kOverwrite);
            cout << "here" << endl;
            histoMassMesonPi0MinusExp->Write("MassPi0MinusExp",TObject::kOverwrite);
            histoTrueMassMesonPi0->Write("TrueMassPi0",TObject::kOverwrite);
            histoTrueMassMesonPi0MinusExp->Write("TrueMassPi0MinusExp",TObject::kOverwrite);
            histoTrueFWHMMesonPi0MeV->Write("TrueFWHMPi0MeV",TObject::kOverwrite);
            histoFWHMMesonPi0MeV->Write("FWHMPi0MeV",TObject::kOverwrite);
            histoPi0CorrYieldBinShifted->Write("CorrectedYieldPi0BinShifted",TObject::kOverwrite);
            histoInvCrossSectionPi0->Write("InvCrossSectionPi0",TObject::kOverwrite);
            graphInvCrossSectionSysPi0->Write("InvCrossSectionPi0Sys",TObject::kOverwrite);
            graphInvCrossSectionSysAPi0->Write("InvCrossSectionPi0SysA",TObject::kOverwrite);
            graphCorrectedYieldPi0SysErr->Write("Pi0SystError",TObject::kOverwrite);
            graphCorrectedYieldPi0SysErrA->Write("Pi0SystErrorA",TObject::kOverwrite);
            graphCorrectedYieldPi0SysErrBinShifted->Write("Pi0SystErrorBinShifted",TObject::kOverwrite);
            graphCorrectedYieldPi0SysErrABinShifted->Write("Pi0SystErrorABinShifted",TObject::kOverwrite);
            graphCorrectedYieldPi0StatPlusSys->Write("Pi0ComplError",TObject::kOverwrite);
            graphCorrectedYieldPi0StatPlusSysBinShifted->Write("Pi0ComplErrorBinShifted",TObject::kOverwrite);
            histoMCInputPi0->Write("Pi0_Input_Reweighted",TObject::kOverwrite); //double saving, but I need this name for the weighting
            if (histoMCInputPi0WOWeight) histoMCInputPi0WOWeight->Write("Pi0_Input",TObject::kOverwrite);
            if (histoMCInputPi0Weights) histoMCInputPi0Weights->Write("Pi0_Weights",TObject::kOverwrite);
            if (histoMCInputPi0AddedSig)histoMCInputPi0AddedSig->Write("Pi0_Input_Reweighted_AddedSig",TObject::kOverwrite);
            if (histoMCInputPi0WOWeightAddedSig) histoMCInputPi0WOWeightAddedSig->Write("Pi0_Input_AddedSig",TObject::kOverwrite);
            if (histoMCInputPi0WeightsAddedSig) histoMCInputPi0WeightsAddedSig->Write("Pi0_Weights_AddedSig",TObject::kOverwrite);
            if (histoSignalPlusBGPi0) histoSignalPlusBGPi0->Write(Form("InvMassSigPlusBG_PtBin%02d",fExampleBinPi0),TObject::kOverwrite);
            if (histoBGPi0) histoBGPi0->Write(Form("InvMassBG_PtBin%02d",fExampleBinPi0),TObject::kOverwrite);
            if (histoSignalPi0) histoSignalPi0->Write(Form("InvMassSig_PtBin%02d",fExampleBinPi0),TObject::kOverwrite);
            if (fitSignalPi0) fitSignalPi0->Write(Form("FitInvMassSig_PtBin%02d",fExampleBinPi0),TObject::kOverwrite);
            
            
            cout<< "bis hier her hab ichs geschafft" << endl;
        }
        
        cout << "done" << endl;
        
        
        fileOutputForComparisonFullyCorrected->mkdir(Form("Eta%s%s",optDalitz.Data(),optionEnergy.Data()));
        TDirectoryFile* directoryEta = (TDirectoryFile*)fileOutputForComparisonFullyCorrected->Get(Form("Eta%s%s",optDalitz.Data(),optionEnergy.Data())); 
        fileOutputForComparisonFullyCorrected->cd(Form("Eta%s%s",optDalitz.Data(),optionEnergy.Data()));
        if (useSameBinningPi0Eta.CompareTo("")==0){
            histoCorrectedYieldEta->Write("CorrectedYieldEta",TObject::kOverwrite);
            histoUnCorrectedYieldEta->Write("RAWYieldPerEventsEta",TObject::kOverwrite);
            histoAccEta->Write("AcceptanceEta",TObject::kOverwrite);
            histoTrueEffPtEta->Write("EfficiencyEta",TObject::kOverwrite);
            histoFWHMMesonEta->Write("FWHMEta",TObject::kOverwrite);
            histoMassMesonEta->Write("MassEta",TObject::kOverwrite);
            histoMassMesonEtaMinusExp->Write("MassEtaMinusExp",TObject::kOverwrite);
            histoTrueMassMesonEtaMinusExp->Write("TrueMassEtaMinusExp",TObject::kOverwrite);
            histoTrueMassMesonEta->Write("TrueMassEta",TObject::kOverwrite);
            histoFWHMMesonEtaMeV->Write("FWHMEtaMeV",TObject::kOverwrite);
            histoTrueFWHMMesonEtaMeV->Write("TrueFWHMEtaMeV",TObject::kOverwrite);
            histoEtaCorrYieldBinShifted->Write("CorrectedYieldEtaBinShifted",TObject::kOverwrite);
            histoInvCrossSectionEta->Write("InvCrossSectionEta",TObject::kOverwrite);
            graphInvCrossSectionSysEta->Write("InvCrossSectionEtaSys",TObject::kOverwrite);
            graphInvCrossSectionSysAEta->Write("InvCrossSectionEtaSysA",TObject::kOverwrite);
            graphCorrectedYieldEtaSysErr->Write("EtaSystError",TObject::kOverwrite);
            graphCorrectedYieldEtaSysErrA->Write("EtaSystErrorA",TObject::kOverwrite);
            graphCorrectedYieldEtaSysErrBinShifted->Write("EtaSystErrorBinShifted",TObject::kOverwrite);
            graphCorrectedYieldEtaSysErrABinShifted->Write("EtaSystErrorABinShifted",TObject::kOverwrite);
            graphCorrectedYieldEtaStatPlusSys->Write("EtaComplError",TObject::kOverwrite);
            graphCorrectedYieldEtaStatPlusSysBinShifted->Write("EtaSystErrorBinShifted",TObject::kOverwrite);
            histoMCInputEta->Write("Eta_Input_Reweighted",TObject::kOverwrite);
            if (histoMCInputEtaWOWeight) histoMCInputEtaWOWeight->Write("Eta_Input",TObject::kOverwrite);
            if (histoMCInputEtaWeights) histoMCInputEtaWeights->Write("Eta_Weights",TObject::kOverwrite);
            if (histoMCInputEtaAddedSig)histoMCInputEtaAddedSig->Write("Eta_Input_Reweighted_AddedSig",TObject::kOverwrite);
            if (histoMCInputEtaWOWeightAddedSig) histoMCInputEtaWOWeightAddedSig->Write("Eta_Input_AddedSig",TObject::kOverwrite);
            if (histoMCInputEtaWeightsAddedSig) histoMCInputEtaWeightsAddedSig->Write("Eta_Weights_AddedSig",TObject::kOverwrite);
            if (histoSignalPlusBGEta) histoSignalPlusBGEta->Write(Form("InvMassSigPlusBG_PtBin%02d",fExampleBinEta),TObject::kOverwrite);
            if (histoBGEta) histoBGEta->Write(Form("InvMassBG_PtBin%02d",fExampleBinEta),TObject::kOverwrite);
            if (histoSignalEta) histoSignalEta->Write(Form("InvMassSig_PtBin%02d",fExampleBinEta),TObject::kOverwrite);
            if (fitSignalEta) fitSignalEta->Write(Form("FitInvMassSig_PtBin%02d",fExampleBinEta),TObject::kOverwrite);

            
            cout << "done" << endl;
        }
        if (useSameBinningPi0Eta.CompareTo("")!=0){
            histoRatioEtaPi0->Write("EtatoPi0RatioConversion",TObject::kOverwrite);
            graphSystErrRatio->Write("EtatoPi0RatioConversionSys",TObject::kOverwrite);
            histoRatioEtaPi0BinShifted->Write("EtatoPi0RatioConversionBinShifted",TObject::kOverwrite);
            graphSystErrRatioBinShifted->Write("EtatoPi0RatioConversionBinShiftedSys",TObject::kOverwrite);
        }
        directoryPi0->mkdir("Fits");
        directoryPi0->cd("Fits");
        if (useSameBinningPi0Eta.CompareTo("")==0){
            cout << "done" << endl;
            fitBinShifting->Write("fitBinShiftingPi0",TObject::kOverwrite);
            fitPtLevyPi0->Write("fitPtLevyPi0",TObject::kOverwrite);
            fitPtHagedornPi0->Write("fitPtHagedornPi0",TObject::kOverwrite);
            fitPtPowerlawPi0->Write("fitPtPowerlawPi0",TObject::kOverwrite);
            fitPtModPowerlawPi0->Write("fitPtModPowerlawPi0",TObject::kOverwrite);
            fitPtBoltzmannPi0->Write("fitPtBoltzmannPi0",TObject::kOverwrite);
            fitPtExpPi0->Write("fitPtExpPi0",TObject::kOverwrite);
            fitPtLevyPi0SysErr->Write("fitPtLevyPi0SysErr",TObject::kOverwrite);
            fitPtHagedornPi0SysErr->Write("fitPtHagedornPi0SysErr",TObject::kOverwrite);
            fitPtPowerlawPi0SysErr->Write("fitPtPowerlawPi0SysErr",TObject::kOverwrite);
            fitPtModPowerlawPi0SysErr->Write("fitPtModPowerlawPi0SysErr",TObject::kOverwrite);
            fitPtBoltzmannPi0SysErr->Write("fitPtBoltzmannPi0SysErr",TObject::kOverwrite);
            fitPtExpPi0SysErr->Write("fitPtExpPi0SysErr",TObject::kOverwrite);
            cout << "done" << endl;
        }
        directoryEta->mkdir("Fits");
        directoryEta->cd("Fits");
        if (useSameBinningPi0Eta.CompareTo("")==0){
            fitPtLevyEta->Write("fitPtLevyEta",TObject::kOverwrite);
            fitPtHagedornEta->Write("fitPtHagedornEta",TObject::kOverwrite);
            fitPtPowerlawEta->Write("fitPtPowerlawEta",TObject::kOverwrite);
            fitPtModPowerlawEta->Write("fitPtModPowerlawEta",TObject::kOverwrite);
            fitBinShiftingEta->Write("fitBinShiftingEta",TObject::kOverwrite);
            fitPtLevyEtaSysErr->Write("fitPtLevyEtaSysErr",TObject::kOverwrite);
            fitPtHagedornEtaSysErr->Write("fitPtHagedornEtaSysErr",TObject::kOverwrite);
            fitPtPowerlawEtaSysErr->Write("fitPtPowerlawEtaSysErr",TObject::kOverwrite);
            fitPtModPowerlawEtaSysErr->Write("fitPtModPowerlawEtaSysErr",TObject::kOverwrite);
            cout << "done" << endl;
        }
    fileOutputForComparisonFullyCorrected->Write();
    fileOutputForComparisonFullyCorrected->Close();
    
   if (multFlag){}
}
