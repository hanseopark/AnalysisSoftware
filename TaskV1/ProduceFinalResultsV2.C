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
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "ProduceFinalResultsV2.h"

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //    TString name;
};

void  ProduceFinalResultsV2( const char *fileNamePi0 = "myOutput",
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
                           Int_t mode = 9,
                           Int_t kColorScheme = 2){    //0:standard, 1:energy, 2:method
    
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();
    SetPlotStyle();
    colorScheme=kColorScheme;
    Bool_t kDalitz = kFALSE;
    if (optDalitz.CompareTo("Dalitz")==0){
        kDalitz = kTRUE;
    }
    
    detSystem = "PCM";
    if (mode == 2) detSystem     = "PCM-EMCAL";
    if (mode == 3) detSystem     = "PCM-PHOS";
    if (mode == 4) detSystem     = "EMCAL-EMCAL";
    if (mode == 5) detSystem     = "PHOS-PHOS";
    
    date = ReturnDateString();
    TString fileNameSysErrPi0 ="SystematicErrorsNew/SystematicErrorAveraged_Pi0_7TeV_24_Apr_2012.dat"; // default
    collisionSystem= ReturnFullCollisionsSystem(optionEnergy);

    // original systematic errors were from 18_May_2011 inversion of material error caused changing
    if(optionEnergy.CompareTo("7TeV") == 0){
        minPtForFits=0.3;
        minPtForFitsEta=0.4;
        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0_7TeV_24_Apr_2012.dat"; //SystematicError_Pi0_7TeV_12_Mar_2012.dat
        } else {
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_7TeV_24_Apr_2012.dat"; // SystematicError_Pi0EtaBinning_7TeV_12_Mar_2012.dat
            minPtForFits=0.4;
        }    
        fileNameSysErrEta = "SystematicErrorsNew/SystematicErrorAveraged_Eta_7TeV_24_Apr_2012.dat";//SystematicError_Eta_7TeV_12_Mar_2012.dat
        cout << "You have choosen 7TeV" << endl;
    } else if( optionEnergy.CompareTo("8TeV") == 0) {
        minPtForFits=0.4;
        minPtForFitsEta=0.4;

        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fileNameSysErrPi0 = "SystematicErrorsCalculated/SystematicErrorAveraged_Pi0_8TeV_2016_03_10.dat";
        } else {
            fileNameSysErrPi0 = "SystematicErrorsCalculated/SystematicErrorAveraged_Pi0EtaBinning_8TeV_2016_03_10.dat"; 
            minPtForFits=0.4;
        }    
        fileNameSysErrEta = "SystematicErrorsCalculated/SystematicErrorAveraged_Eta_8TeV_2016_03_10.dat"; 
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
        
    } else if( optionEnergy.CompareTo("5TeV") == 0) {
        minPtForFits=0.4;
        minPtForFitsEta=0.4;

        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0_2.76TeV_5_Aug_2013.dat";
        } else {
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_2.76TeV_5_Aug_2013.dat"; 
            minPtForFits=0.4;
        }    
        fileNameSysErrEta = "SystematicErrorsNew/SystematicErrorAveraged_Eta_2.76TeV_5_Aug_2013.dat"; 
        cout << "You have choosen 5TeV. Note that you will use 2.76TeV systematic errors." << endl;
        
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
        if (mode == 2 || mode == 4){
            cout << "wrong macro, please look at ProduceFinalResultsPatchedTriggers.C" << endl;
        }
    } else if( optionEnergy.CompareTo("900GeV") == 0) {
        minPtForFits=0.4;
        minPtForFitsEta=0.9;
        if (useSameBinningPi0Eta.CompareTo("")==0){         
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0_900GeV_24_Apr_2012.dat"; //SystematicError_Pi0_900GeV_12_Mar_2012.dat
        } else {
            fileNameSysErrPi0 = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_900GeV_24_Apr_2012.dat"; //SystematicError_Pi0EtaBinning_900GeV_12_Mar_2012.dat
            minPtForFits=0.9;
        }
        fileNameSysErrEta = "SystematicErrorsNew/SystematicErrorAveraged_Eta_900GeV_24_Apr_2012.dat"; //SystematicError_Eta_900GeV_12_Mar_2012.dat
    } else if( optionEnergy.CompareTo("PbPb_2.76TeV") == 0 || optionEnergy.CompareTo("pPb_5.023TeV") == 0  ) {
        cout << "This macro will not work for pPb or PbPb" << endl;
      return;
    } else {
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    
    // select the example invariant mass bin for the output
    SelectExampleBin(optionEnergy,useSameBinningPi0Eta);
    
    TString outputDir = Form("%s/%s/%s/ProduceFinalResults%s",cutSelection.Data(),optionEnergy.Data(),suffix, optDalitz.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    
    if(conferencePlots.CompareTo("conference") == 0){// means we want to plot values for the pi0
        conference = kTRUE;
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
    if (!isMC.CompareTo("kTRUE")){ 
        prefix2 = "MC";
    } else {    
        prefix2 = "data";
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
    
    nEvt =  GetNEvents(histoEventQualtityPi0);
    
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
    
    
    TFile* fileMCGenerated =     new TFile("ExternalInput/PCM/MCGeneratedSpectra.root");
    if(!optionEnergy.CompareTo("8TeV")||!optionEnergy.CompareTo("13TeV")||!optionEnergy.CompareTo("5TeV")){
		cout << "Caution!!! using 7TeV MC generated EtaToPi0 spectra" << endl;
		use7TeVPytPho = kTRUE;
        histoEtaToPi0Phojet = (TH1D*)fileMCGenerated->Get(Form("EtaToPi0_generatedSpectrum_%s_Phojet","7TeV"));
        histoEtaToPi0Pythia = (TH1D*)fileMCGenerated->Get(Form("EtaToPi0_generatedSpectrum_%s_Pythia","7TeV"));
    
        histoPi0ToChargedPhojet = (TH1D*)fileMCGenerated->Get(Form("Pi0ToCharged_generatedSpectrum_%s_Phojet","7TeV"));
        histoPi0ToChargedPythia = (TH1D*)fileMCGenerated->Get(Form("Pi0ToCharged_generatedSpectrum_%s_Pythia","7TeV"));
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
    
    
   
   if (!isMC.CompareTo("kFALSE")&&!useSameBinningPi0Eta.CompareTo("")){
		//**********************************************************************************
		//******************** FWHM Plot ***************************************************
		//**********************************************************************************
        Double_t maxYFWHM=2.0*histoFWHMMesonPi0->GetBinContent(histoFWHMMesonPi0->GetMaximumBin());
		PlotFinalOutput("FWHMpi0",histoFWHMMesonPi0,histoTrueFWHMMesonPi0,NULL,NULL,0.00,maxYFWHM,"#sigma_{#pi^{0}} (Gev/#it{c}^{2})","#pi^{0}",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		maxYFWHM=2.0*histoFWHMMesonEta->GetBinContent(histoFWHMMesonEta->GetMaximumBin());
		PlotFinalOutput("FWHMeta",histoFWHMMesonEta,histoTrueFWHMMesonEta,NULL,NULL,0.00,maxYFWHM,"#sigma_{#eta} (Gev/#it{c}^{2})","#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		PlotFinalOutput("FWHMComb",histoFWHMMesonPi0,histoTrueFWHMMesonPi0,histoFWHMMesonEta,histoTrueFWHMMesonEta,0.00,maxYFWHM,"#sigma_{#pi^{0}/#eta} (Gev/#it{c}^{2})","#pi^{0}/#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kTRUE);   
   
	   //**********************************************************************************
	   //******************** Mass Plot ***************************************************
	   //**********************************************************************************
	   Double_t minYMass=0.97*histoMassMesonPi0->GetBinContent(histoMassMesonPi0->GetMaximumBin());
	   Double_t maxYMass=1.03*histoMassMesonPi0->GetBinContent(histoMassMesonPi0->GetMaximumBin());
	   PlotFinalOutput("Masspi0",histoMassMesonPi0,histoTrueMassMesonPi0,NULL,NULL,minYMass,maxYMass,"#it{M}_{#pi^{0}} (GeV/c^{2})","#pi^{0}",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
	   
	   minYMass=0.97*histoMassMesonEta->GetBinContent(histoMassMesonEta->GetMaximumBin());
	   maxYMass=1.03*histoMassMesonEta->GetBinContent(histoMassMesonEta->GetMaximumBin());
	   PlotFinalOutput("Masseta",histoMassMesonEta,histoTrueMassMesonEta,NULL,NULL,minYMass,maxYMass,"#it{M}_{#eta} (GeV/c^{2})","#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
	   
	   //**********************************************************************************
	   //******************** Combined Mass and FWHM Plot *********************************
	   //**********************************************************************************
	   PlotFWHMMass("MassFWHMPi0",histoFWHMMesonPi0,histoTrueFWHMMesonPi0,histoMassMesonPi0,histoTrueMassMesonPi0,minYMass,maxYMass,"#it{M}_{#eta} (GeV/c^{2})","#eta",optionEnergy,mode,outputDir,prefix2,cutSelection,suffix,"");
	   PlotFWHMMass("MassFWHMEta",histoFWHMMesonEta,histoTrueFWHMMesonEta,histoMassMesonEta,histoTrueMassMesonEta,minYMass,maxYMass,"#it{M}_{#eta} (GeV/c^{2})","#eta",optionEnergy,mode,outputDir,prefix2,cutSelection,suffix,"");
	}

    if (!isMC.CompareTo("kTRUE")&&!useSameBinningPi0Eta.CompareTo("")){                    
        //**********************************************************************************
        //******************** Acceptance Plot *********************************************
        //**********************************************************************************
        Double_t minYacc=0.501;
        Double_t maxYacc=1.09;
        if(mode == 4){minYacc=0.0;maxYacc=1.5;}
		PlotFinalOutput("AccPi0",histoAccPi0,NULL,NULL,NULL,minYacc,maxYacc,Form("A_{#pi^{0}} (|y|< %s)",rapidityRange.Data()),"#pi^{0}",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		PlotFinalOutput("AccEta",histoAccEta,NULL,NULL,NULL,minYacc,maxYacc,Form("A_{#eta} (|y|< %s)",rapidityRange.Data()),"#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		PlotFinalOutput("AccComb",histoAccPi0,NULL,histoAccEta,NULL,minYacc,maxYacc,Form("A (|y|< %s)",rapidityRange.Data()),"#pi^{0}/#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
        
        
        //**********************************************************************************
        //********************* Efficiency Plot ********************************************
        //**********************************************************************************
		Double_t minYeff=3e-7;
        Double_t maxYeff=1.1e-2;
        if(mode == 4){minYacc=0.001;maxYacc=0.5;}
		PlotFinalOutput("EffPi0",histoTrueEffPtPi0,NULL,NULL,NULL,minYeff,maxYeff,"#epsilon_{#pi^{0}}","#pi^{0}",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		PlotFinalOutput("EffEta",histoTrueEffPtEta,NULL,NULL,NULL,minYeff,maxYeff,"#epsilon_{#eta}","#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		PlotFinalOutput("EffComb",histoTrueEffPtPi0,NULL,histoTrueEffPtEta,NULL,minYeff,maxYeff,"#epsilon_{eff}","#pi^{0}/#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
        
        //**********************************************************************************
        //********************* Acc x Eff Plot *********************************************
        //**********************************************************************************
		PlotFinalOutput("AccEffComb",histoTrueEffPtPi0,histoAccPi0,histoTrueEffPtEta,histoAccEta,minYeff,maxYeff,"Acceptance x Efficiency","#pi^{0}/#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
	}
    //**********************************************************************************
    //******************** RAW Yield spectrum ******************************************
    //**********************************************************************************
    
	histoUncorrectedYieldPi0->Scale(1./nEvt);
    histoUnCorrectedYieldEta->Scale(1./nEvt);
	Double_t minYieldRawPi0    = 0.2*histoUncorrectedYieldPi0->GetBinContent(histoUncorrectedYieldPi0->GetNbinsX());
    Double_t maxYieldRawPi0    = 5.0*histoUncorrectedYieldPi0->GetBinContent(histoUncorrectedYieldPi0->GetMaximumBin());
	Double_t minYieldRawEta    = 0.2*histoUnCorrectedYieldEta->GetBinContent(histoUnCorrectedYieldEta->GetNbinsX());
    Double_t maxYieldRawEta    = 5.0*histoUnCorrectedYieldEta->GetBinContent(histoUnCorrectedYieldEta->GetMaximumBin());
    if (!isMC.CompareTo("kFALSE")&&!useSameBinningPi0Eta.CompareTo("")){
		PlotFinalOutput("RawPi0",histoUncorrectedYieldPi0,NULL,NULL,NULL,minYieldRawPi0,maxYieldRawPi0,"#frac{d#it{N}_{#pi_{0}, raw}}{#it{N}_{evt}d#it{p}_{T}} (#it{c}/GeV)^{2}","#pi^{0}",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		PlotFinalOutput("RawEta",histoUnCorrectedYieldEta,NULL,NULL,NULL,minYieldRawEta,maxYieldRawEta,"#frac{d#it{N}_{#eta, raw}}{#it{N}_{evt}d#it{p}_{T}} (#it{c}/GeV)^{2}","#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		PlotFinalOutput("RawComb",histoUncorrectedYieldPi0,NULL,histoUnCorrectedYieldEta,NULL,minYieldRawPi0,maxYieldRawPi0,"#frac{d#it{N}_{#pi_{0}/#eta, raw}}{#it{N}_{evt}d#it{p}_{T}} (#it{c}/GeV)^{2}","#pi^{0}/#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
	}
    
    Double_t minCorrYieldPi0 = 0.2*histoCorrectedYieldPi0->GetBinContent(histoCorrectedYieldPi0->GetNbinsX());
    Double_t maxCorrYieldPi0 = 5.0*histoCorrectedYieldPi0->GetBinContent(histoCorrectedYieldPi0->GetMaximumBin());
    Double_t minCorrYieldEta = 0.2*histoCorrectedYieldEta->GetBinContent(histoCorrectedYieldEta->GetNbinsX());
    Double_t maxCorrYieldEta = 5.0*histoCorrectedYieldEta->GetBinContent(histoCorrectedYieldEta->GetMaximumBin());

    if (!isMC.CompareTo("kFALSE")&&!useSameBinningPi0Eta.CompareTo("")){
		PlotFinalOutput("CorrPi0",histoCorrectedYieldPi0,NULL,NULL,NULL,minCorrYieldPi0,maxCorrYieldPi0,"#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}","#pi^{0}",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		PlotFinalOutput("CorrEta",histoCorrectedYieldEta,NULL,NULL,NULL,minCorrYieldEta,maxCorrYieldEta,"#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}","#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
		PlotFinalOutput("CorrComb",histoCorrectedYieldPi0,NULL,histoCorrectedYieldEta,NULL,minCorrYieldPi0,maxCorrYieldPi0,"#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}","#pi^{0}/#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"");
	}

    //********************************************************************************************
    //*********************** Systematic Error and Fitting ***************************************
    //********************************************************************************************
    if( useSameBinningPi0Eta.CompareTo("")==0){
        //*************************************** Systematic Error of the Pi0 *************************************************
		graphCorrectedYieldPi0SysErr = CalculateSysErrFromRelSysHisto( histoCorrectedYieldPi0, "Pi0SystError",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
		graphCorrectedYieldPi0SysErrA = CalculateSysErrFromRelSysHisto( histoCorrectedYieldPi0, "Pi0SystErrorA",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
		graphCorrectedYieldPi0StatPlusSys = CalculateSysErrFromRelSysHistoComplete( histoCorrectedYieldPi0, "Pi0ComplError",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
		histoInvCrossSectionPi0 = (TH1D*) histoCorrectedYieldPi0->Clone();
		histoInvCrossSectionPi0->Scale(xSection*recalcBarn);
		graphInvCrossSectionSysPi0 = CalculateSysErrFromRelSysHisto( histoInvCrossSectionPi0 , "Pi0InvCrossSectionSys",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);
		graphInvCrossSectionSysAPi0 = CalculateSysErrFromRelSysHisto( histoInvCrossSectionPi0 , "Pi0InvCrossSectionSysA",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
        
		//***************************************************************************************************
		//*************************** Fitting  Pi0 Spectrum *************************************************
		//***************************************************************************************************
        if (!isMC.CompareTo("kFALSE")){
			PlotFinalOutput("CorrPi0SysErr",histoCorrectedYieldPi0,NULL,NULL,NULL,minCorrYieldPi0,maxCorrYieldPi0,"#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}","#pi^{0}",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kFALSE,graphCorrectedYieldPi0SysErr);
			PlotFinalOutput("CorrPi0SysErrFitted",histoCorrectedYieldPi0,NULL,NULL,NULL,minCorrYieldPi0,maxCorrYieldPi0,"#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}","#pi^{0}",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kFALSE,graphCorrectedYieldPi0SysErr);
			PlotFinalOutput("CorrPi0SysErrFittedRatio",histoCorrectedYieldPi0,NULL,NULL,NULL,0.55,1.85,"Data/Fit","#pi^{0}",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kFALSE,graphCorrectedYieldPi0SysErr);
		}
		
		//*************************************** Systematic Error of the Eta *************************************************
        graphCorrectedYieldEtaSysErr = CalculateSysErrFromRelSysHisto( histoCorrectedYieldEta, "EtaSystError",relSystErrorEtaDown , relSystErrorEtaUp, 2, nPointsEta);
        graphCorrectedYieldEtaSysErrA = CalculateSysErrFromRelSysHisto( histoCorrectedYieldEta, "EtaSystErrorA",relSystErrorWOMaterialEtaDown , relSystErrorWOMaterialEtaUp, 2, nPointsEta);
        graphCorrectedYieldEtaStatPlusSys = CalculateSysErrFromRelSysHistoComplete( histoCorrectedYieldEta, "EtaComplError",relSystErrorEtaDown, relSystErrorEtaUp, 2, nPointsEta); histoInvCrossSectionEta = (TH1D*) histoCorrectedYieldEta->Clone();
        histoInvCrossSectionEta->Scale(xSection*recalcBarn);
        graphInvCrossSectionSysEta = CalculateSysErrFromRelSysHisto( histoInvCrossSectionEta , "Pi0InvCrossSectionSys",relSystErrorEtaDown , relSystErrorEtaUp, 2, nPointsEta);
        graphInvCrossSectionSysAEta = CalculateSysErrFromRelSysHisto( histoInvCrossSectionEta , "Pi0InvCrossSectionSysA",relSystErrorWOMaterialEtaDown , relSystErrorWOMaterialEtaUp, 2, nPointsEta);
		
		//***************************************************************************************************
		//*************************** Fitting  Eta Spectrum *************************************************
		//***************************************************************************************************
        if (!isMC.CompareTo("kFALSE")){
			PlotFinalOutput("CorrEtaSysErr",histoCorrectedYieldEta,NULL,NULL,NULL,minCorrYieldEta,maxCorrYieldEta,"#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}","#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kFALSE,graphCorrectedYieldEtaSysErr);
			PlotFinalOutput("CorrEtaSysErrFitted",histoCorrectedYieldEta,NULL,NULL,NULL,minCorrYieldEta,maxCorrYieldEta,"#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}","#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kFALSE,graphCorrectedYieldEtaSysErr);
			PlotFinalOutput("CorrEtaSysErrFittedRatio",histoCorrectedYieldEta,NULL,NULL,NULL,0.55,1.85,"Data/Fit","#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kFALSE,graphCorrectedYieldEtaSysErr);
        }
    }

    //************************************************************************************************
    //******************************** Combined Plot Pi0 and Eta *************************************
    //************************************************************************************************
    if (!isMC.CompareTo("kFALSE")&&!useSameBinningPi0Eta.CompareTo("")){
		PlotFinalOutput("CorrCombSysErr",histoCorrectedYieldPi0,NULL,histoCorrectedYieldEta,NULL,minCorrYieldPi0,maxCorrYieldPi0,"#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}","#pi^{0}/#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kFALSE,graphCorrectedYieldPi0SysErr,graphCorrectedYieldEtaSysErr);
	}
    TCanvas * canvasCorrectedSame;
    if( useSameBinningPi0Eta.CompareTo("")!=0){
        histoRatioEtaPi0 = (TH1D*) histoCorrectedYieldEta->Clone();
        histoRatioEtaPi0->Divide(histoRatioEtaPi0,histoCorrectedYieldPi0,1.,1.,"");

        //***************************** Only ratio *************************************************
        cout << " calculation of eta/pi0" << endl;
        for (Int_t i = 0; i < nPointsEta + 1; i++){
            ratioXValue[i] = histoRatioEtaPi0->GetBinCenter(i+2);
            ratioYValue[i] = histoRatioEtaPi0->GetBinContent(i+2);
            ratioXError[i] = histoRatioEtaPi0->GetBinWidth(i+2)/2.;

            ratioSysUpError[i]= TMath::Sqrt(TMath::Power(relSystErrorWOMaterialEtaUp[i]/100* histoCorrectedYieldEta->GetBinContent(i+2)/histoCorrectedYieldPi0->GetBinContent(i+2),2) + TMath::Power(histoCorrectedYieldEta->GetBinContent(i+2)*relSystErrorWOMaterialPi0Up[i]/100/histoCorrectedYieldPi0->GetBinContent(i+2),2));
            ratioSysDownError[i]= TMath::Sqrt(TMath::Power(relSystErrorWOMaterialEtaDown[i]/100* histoCorrectedYieldEta->GetBinContent(i+2)/histoCorrectedYieldPi0->GetBinContent(i+2),2) + TMath::Power(histoCorrectedYieldEta->GetBinContent(i+2)*relSystErrorWOMaterialPi0Down[i]/100/histoCorrectedYieldPi0->GetBinContent(i+2),2));
            
            cout << ratioXValue[i] << "\t" << ratioYValue[i] << "\t" << ratioSysUpError[i] << "\t" << ratioSysDownError[i] << endl;
        }
        graphSystErrRatio = new TGraphAsymmErrors(nPointsEta,ratioXValue,ratioYValue,ratioXError,ratioXError,ratioSysDownError,ratioSysUpError);
        
		if (!isMC.CompareTo("kFALSE")){
			PlotFinalOutput("EtaToPi0Ratio",histoCorrectedYieldPi0,NULL,histoCorrectedYieldEta,NULL,0.01,1.05,"#eta/#pi^{0}","#pi^{0}/#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kFALSE,graphSystErrRatio);
			//******************************* Ratio + Pythia ********************************************
			plotPythiaPhojetInRatio=kTRUE;
			PlotFinalOutput("EtaToPi0RatioPythiaPhojet",histoCorrectedYieldPi0,NULL,histoCorrectedYieldEta,NULL,0.01,1.05,"#eta/#pi^{0}","#pi^{0}/#eta",optionEnergy,mode,outputDir,prefix2,useSameBinningPi0Eta,cutSelection,suffix,"",kFALSE,graphSystErrRatio);
		}
	}

    
    //*********************************************************************************************************
    //********************** ComparisonFile Output ************************************************************
    //*********************************************************************************************************
    TString applyBinShift = "";
    if( optNoBinShift.CompareTo("kTRUE")==0){
        applyBinShift = "NoBinShifting";
        cout << applyBinShift << endl;
    }

    const char* fileNameOutputComp = Form("%s/%s_%sResultsFullCorrection_PP_%s.root",cutSelection.Data(),prefix2.Data(),detSystem.Data(), applyBinShift.Data());
    TFile* fileOutputForComparisonFullyCorrected = new TFile(fileNameOutputComp,"UPDATE");
        if (useSameBinningPi0Eta.CompareTo("")==0){
            histoNumberOfEvents->Write(Form("histoNumberOfEvents%s",optionEnergy.Data()),TObject::kOverwrite);
        }
        fileOutputForComparisonFullyCorrected->mkdir(Form("Pi0%s%s",optDalitz.Data(),optionEnergy.Data()));
        TDirectoryFile* directoryPi0 = (TDirectoryFile*)fileOutputForComparisonFullyCorrected->Get(Form("Pi0%s%s",optDalitz.Data(),optionEnergy.Data())); 
        fileOutputForComparisonFullyCorrected->cd(Form("Pi0%s%s",optDalitz.Data(),optionEnergy.Data()));
        
        if (useSameBinningPi0Eta.CompareTo("")==0){
            histoCorrectedYieldPi0->Write("CorrectedYieldPi0",TObject::kOverwrite);
            histoUncorrectedYieldPi0->Write("RAWYieldPerEventsPi0",TObject::kOverwrite);
            histoAccPi0->Write("AcceptancePi0",TObject::kOverwrite);
            histoTrueEffPtPi0->Write("EfficiencyPi0",TObject::kOverwrite);
            histoFWHMMesonPi0->Write("FWHMPi0",TObject::kOverwrite);
            histoMassMesonPi0->Write("MassPi0",TObject::kOverwrite);
            histoMassMesonPi0MinusExp->Write("MassPi0MinusExp",TObject::kOverwrite);
            histoTrueMassMesonPi0->Write("TrueMassPi0",TObject::kOverwrite);
            histoTrueMassMesonPi0MinusExp->Write("TrueMassPi0MinusExp",TObject::kOverwrite);
            histoTrueFWHMMesonPi0MeV->Write("TrueFWHMPi0MeV",TObject::kOverwrite);
            histoFWHMMesonPi0MeV->Write("FWHMPi0MeV",TObject::kOverwrite);
            histoInvCrossSectionPi0->Write("InvCrossSectionPi0",TObject::kOverwrite);
            graphInvCrossSectionSysPi0->Write("InvCrossSectionPi0Sys",TObject::kOverwrite);
            graphInvCrossSectionSysAPi0->Write("InvCrossSectionPi0SysA",TObject::kOverwrite);
            graphCorrectedYieldPi0SysErr->Write("Pi0SystError",TObject::kOverwrite);
            graphCorrectedYieldPi0SysErrA->Write("Pi0SystErrorA",TObject::kOverwrite);
            graphCorrectedYieldPi0StatPlusSys->Write("Pi0ComplError",TObject::kOverwrite);
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
        }

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
            histoInvCrossSectionEta->Write("InvCrossSectionEta",TObject::kOverwrite);
            graphInvCrossSectionSysEta->Write("InvCrossSectionEtaSys",TObject::kOverwrite);
            graphInvCrossSectionSysAEta->Write("InvCrossSectionEtaSysA",TObject::kOverwrite);
            graphCorrectedYieldEtaSysErr->Write("EtaSystError",TObject::kOverwrite);
            graphCorrectedYieldEtaSysErrA->Write("EtaSystErrorA",TObject::kOverwrite);
            graphCorrectedYieldEtaStatPlusSys->Write("EtaComplError",TObject::kOverwrite);
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
        }
        if (useSameBinningPi0Eta.CompareTo("")!=0){
            histoRatioEtaPi0->Write("EtatoPi0RatioConversion",TObject::kOverwrite);
            graphSystErrRatio->Write("EtatoPi0RatioConversionSys",TObject::kOverwrite);
        }
        
    fileOutputForComparisonFullyCorrected->Write();
    fileOutputForComparisonFullyCorrected->Close();
    
   if (multFlag){}
}
