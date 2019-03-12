/****************************************************************************************************************************
 ******     provided by Gamma Conversion Group, PWGGA,                                                                  *****
 ******        Friederike Bock, friederike.bock@cern.ch                                                                 *****
 ******        Daniel Muehlheim, d.muehlheim@cern.ch                                                                    *****
 *****************************************************************************************************************************/

#include <Riostream.h>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <string>
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
#include "../CommonHeaders/ExtractSignalBinning.h"
// #include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //    TString name;
};

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos) return false;
    str.replace(start_pos, from.length(), to);
    return true;
}
Int_t GetBinPosInVec( std::vector<TString> *vec, Int_t size, Double_t lookup){
  for(Int_t i=0; i<size; i++){
    if(((TString)vec[i].at(0)).Atof() == lookup) return i;
  }
  return -1;
}

//***************************************************************************************************************
//************************* Get correct order of triggers for combinations **************************************
//***************************************************************************************************************
Int_t GetOrderedTrigger(TString triggerNameDummy){
    if ((triggerNameDummy.CompareTo("MB") == 0 || triggerNameDummy.CompareTo("INT1") == 0  || triggerNameDummy.CompareTo("MB_NLM2") == 0 || triggerNameDummy.CompareTo("INT1_NLM2") == 0 ) ){
        return 0;
    } else if ((triggerNameDummy.CompareTo("INT7") == 0 || triggerNameDummy.CompareTo("INT7_NLM2") == 0) ){
        return 1;
    } else if ((triggerNameDummy.CompareTo("EMC1") == 0 || triggerNameDummy.CompareTo("EMC1_NLM2") == 0 ) ){
        return 2;
    } else if ((triggerNameDummy.CompareTo("EMC7") == 0 || triggerNameDummy.CompareTo("EMC7_NLM2") == 0 ) ){
        return 3;
    } else if ((triggerNameDummy.CompareTo("EG2") == 0 || triggerNameDummy.CompareTo("EG2_NLM2") == 0 ||  triggerNameDummy.CompareTo("EGA") == 0) ){
        return 4;
    } else if ((triggerNameDummy.CompareTo("EG1") == 0 || triggerNameDummy.CompareTo("EG1_NLM2") == 0 ) ){
        return 5;
    } else if ((triggerNameDummy.CompareTo("MB_NLM1") == 0 || triggerNameDummy.CompareTo("INT1_NLM1") == 0  ) ){
        return 6;
    } else if (triggerNameDummy.CompareTo("INT7_NLM1") == 0 ){
        return 7;
    } else if (triggerNameDummy.CompareTo("EMC1_NLM1") == 0 ){
        return 8;
    } else if (triggerNameDummy.CompareTo("EMC7_NLM1") == 0 ){
        return 9;
    } else if (triggerNameDummy.CompareTo("EG2_NLM1") == 0 ){
        return 10;
    } else if (triggerNameDummy.CompareTo("EG1_NLM1") == 0 ){
        return 11;
    } else
     return -1;
}

//***************************************************************************************************************
//***************************** Main function *******************************************************************
//***************************************************************************************************************
void  ProduceFinalResultsPatchedTriggersEtaPrime(
    TString fileListNameEtaPrime    = "triggerFileListEtaPrime.txt",
    Int_t   mode                    = 4,
    Int_t   numberOfTrigg           = 6,
    TString suffix                  = "eps",
    TString isMC                    = "",
    TString optionEnergy            = "",
    TString period                  = "",
    Bool_t  pileUpApplied           = kTRUE,
    Float_t maxPtGlobalEtaPrime     = 16.,
    Bool_t  averagedEtaPrime        = kFALSE,
    TString nameFileFitsShift       = "",
    Bool_t  hasClusterOutput        = kTRUE,
    TString fileInputCorrFactors    = ""
){

    //***************************************************************************************************************
    //************************************ Layouting preparations & general setup ***********************************
    //***************************************************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();



    if (mode == 10) {
        cout << "This macro can't deal with this mode. Aborting..." << endl;
        std::terminate();
    }

    Int_t modeNormal    = mode;
    if (mode >= 100)
        modeNormal      = mode -100;

    TString dateForOutput           = ReturnDateStringForOutput();
    TString collisionSystem         = ReturnFullCollisionsSystem(optionEnergy);
    TString detectionProcess        = ReturnFullTextReconstructionProcess(mode);
    TString detectionProcessClus    = ReturnFullTextReconstructionProcess(mode,2);

    if (isMC.CompareTo("MC") == 0) collisionSystem = "MC, "+collisionSystem;

    Double_t maxPtGlobalCluster     = 25;
    if (optionEnergy.CompareTo("2.76TeV")==0){
        if (modeNormal==2){
            maxPtGlobalCluster          = 26;
        } else if (modeNormal == 4){
            maxPtGlobalCluster          = 31;
        }
    } else if (optionEnergy.CompareTo("8TeV")==0){
      if(modeNormal==2 || modeNormal==4){
        maxPtGlobalCluster          = 50;
      }
    } else if (optionEnergy.CompareTo("13TeV")==0){
        if(modeNormal==2 || modeNormal==4 || modeNormal==0){
            maxPtGlobalCluster          = 100;
        }
    } else if (optionEnergy.CompareTo("pPb_5.023TeV")==0){
        if(modeNormal==2 || modeNormal==4 || modeNormal==0){
            maxPtGlobalCluster          = 100;
        }
    } else if (optionEnergy.CompareTo("pPb_8TeV")==0){
      if(modeNormal==2 || modeNormal==4){
        maxPtGlobalCluster          = 50;
      }
    }

    Size_t textSizeSpectra          = 0.04;
    Int_t textSizePixelSpectra      = textSizeSpectra*1000;
    const Int_t MaxNumberOfFiles    = 12;

    Color_t colorTrigg      [12]    = { kBlack, kBlack, kBlack, kBlack, kBlack,
                                        kBlack, kBlack, kBlack, kBlack, kBlack,
                                        kBlack, kBlack  };
    Color_t colorTriggShade [12]    = { kGray+1, kGray+1, kGray+1, kGray+1, kGray+1,
                                        kGray+1, kGray+1, kGray+1, kGray+1, kGray+1,
                                        kGray+1, kGray+1 };
    Marker_t markerTrigg    [12]    = { 20, 20, 20, 20, 20,
                                        20, 20, 20, 20, 20,
                                        20, 20 };
    Marker_t markerTriggMC  [12]    = { 24, 24, 24, 24, 24,
                                        24, 24, 24, 24, 24,
                                        24, 24 };
    Size_t sizeTrigg        [12]    = { 1.5, 1.5, 1.5, 1.5, 1.5,
                                        1.5, 1.5, 1.5, 1.5, 1.5,
                                        1.5, 1.5 };

    TString strEG2_A = "EG2";
    if(optionEnergy.CompareTo("8TeV")==0) strEG2_A = "EGA";

    TString system              = "PCM";
    if (modeNormal == 2) system       = "PCM-EMCAL";
    if (modeNormal == 3) system       = "PCM-PHOS";
    if (modeNormal == 4) system       = "EMCAL-EMCAL";
    if (modeNormal == 5) system       = "PHOS-PHOS";


    //***************************************************************************************************************
    //********************************** setting correct histo names ************************************************
    //***************************************************************************************************************
    TString nameCorrectedYield                          = "CorrectedYieldTrueEff";
    TString nameEfficiency                              = "TrueMesonEffiPt";
    TString nameAcceptance                              = "fMCMesonAccepPt";
    TString nameAcceptanceWOEvtWeights                  = "fMCMesonAccepPtWOEvtWeights";
    TString nameMassMC                                  = "histoTrueMassMeson";
    TString nameWidthMC                                 = "histoTrueFWHMMeson";
    TString nameMCYield                                 = "MCYield_Meson_oldBinWOWeights";
    if ( modeNormal == 4 || modeNormal == 5 ){
        nameCorrectedYield                              = "CorrectedYieldNormEff";
        nameEfficiency                                  = "MesonEffiPt";
        nameMassMC                                      = "histoMassMesonRecMC";
        nameWidthMC                                     = "histoFWHMMesonRecMC";
        if (optionEnergy.CompareTo("PbPb_2.76TeV")==0 ){
            nameMassMC                                  = "histoTrueMassMeson";
            nameWidthMC                                 = "histoTrueFWHMMeson";
            nameCorrectedYield                          = "CorrectedYieldTrueEff";
            nameEfficiency                              = "TrueMesonEffiPt";
        }
    } else if ( (modeNormal == 2 || modeNormal == 3) && !(optionEnergy.CompareTo("8TeV")==0 || optionEnergy.CompareTo("pPb_5.023TeV")==0)){
        cout << "using rec quantities for PCM-EMC/PCM-PHOS" << endl;
        nameMassMC                                      = "histoMassMesonRecMC";
        nameWidthMC                                     = "histoFWHMMesonRecMC";
    }

    //***************************************************************************************************************


    TString cutNumber           [MaxNumberOfFiles];
    TString triggerName         [MaxNumberOfFiles];
    TString triggerNameLabel    [MaxNumberOfFiles];
    Float_t minPt               [MaxNumberOfFiles];
    Float_t maxPt               [MaxNumberOfFiles];
    Int_t trigSteps             [MaxNumberOfFiles][3];
    Float_t ptFromSpecEtaPrime  [MaxNumberOfFiles][2];
    Bool_t maskedFullyEtaPrime  [MaxNumberOfFiles];
    TString sysFileEtaPrime     [MaxNumberOfFiles];
    TString cutNumberBaseEff    [MaxNumberOfFiles];
    //***************************************************************************************************************
    //*************************** read setting from configuration file **********************************************
    //***************************************************************************************************************
    ifstream in(fileListNameEtaPrime.Data());
    cout<<"Available Triggers:"<<endl;
    // general number of triggers set
    Int_t nrOfTrigToBeComb      = 0;
    // number of triggers which are really used for the respective analysis
    Int_t nrOfTrigToBeCombEtaPrimeRed        = 0;
    while(!in.eof() && nrOfTrigToBeComb<numberOfTrigg ){
        in
            >> cutNumber[nrOfTrigToBeComb]
            >> minPt[nrOfTrigToBeComb]
            >> maxPt[nrOfTrigToBeComb]
            >> triggerName[nrOfTrigToBeComb]
            >> trigSteps[nrOfTrigToBeComb][0]
            >> trigSteps[nrOfTrigToBeComb][1]
            >> trigSteps[nrOfTrigToBeComb][2]
            >> ptFromSpecEtaPrime[nrOfTrigToBeComb][0]
            >> ptFromSpecEtaPrime[nrOfTrigToBeComb][1]
            >> sysFileEtaPrime[nrOfTrigToBeComb]
            >> cutNumberBaseEff[nrOfTrigToBeComb];
        std::cout
            << cutNumber[nrOfTrigToBeComb]
            << "\t"
            << triggerName[nrOfTrigToBeComb]
            << "\t transverse momentum range: \t"
            << minPt[nrOfTrigToBeComb]
            << " to "
            << maxPt[nrOfTrigToBeComb]
            << std::endl;
        std::cout
            << trigSteps[nrOfTrigToBeComb][0]
            << "\t"
            << trigSteps[nrOfTrigToBeComb][1]
            << "\t"<< trigSteps[nrOfTrigToBeComb][2]
            << std::endl;
        nrOfTrigToBeComb++;
        std::cout << cutNumberBaseEff[nrOfTrigToBeComb] << std::endl;
    }

    for (Int_t i = 0; i < nrOfTrigToBeComb; i++){
        // figure out which triggers are fully masked for the etaprime
        if ((ptFromSpecEtaPrime[i][1] == -1 && ptFromSpecEtaPrime[i][0] == -1 )|| (ptFromSpecEtaPrime[i][0] == 0 && ptFromSpecEtaPrime[i][1] == 0)){
            maskedFullyEtaPrime[i]       = kTRUE;
        } else {
            maskedFullyEtaPrime[i]       = kFALSE;
        }
    }

    Double_t minPtGlobalEtaPrime         = ptFromSpecEtaPrime[0][0];
    if (maskedFullyEtaPrime[0])
        minPtGlobalEtaPrime              = 12;


    for (Int_t j = 1; j < nrOfTrigToBeComb; j++){
        if (minPtGlobalEtaPrime > ptFromSpecEtaPrime[j][0] && !maskedFullyEtaPrime[j] )
            minPtGlobalEtaPrime = ptFromSpecEtaPrime[j][0];
    }
    cout << "global minimum pT for etaprime: " << minPtGlobalEtaPrime << endl;

    // set individual triggers to total number of availabel triggers, will be reduced later according to usage
    nrOfTrigToBeCombEtaPrimeRed      = nrOfTrigToBeComb;

    // variables to keep track of NLM
    TString fNLMString      = "";
    TString fNLMStringOutput= "";
    TString fMergedClusterCutNrExampl   = "";
    Int_t fNLMmin           = 0;
    TString fCent           = "";
    TString fCentOutput     = "";

    // put correct color setting for different triggers
    for (Int_t i = 0; i < numberOfTrigg; i++){
        colorTrigg[i]       = GetDefaultTriggerColorName(triggerName[i], 0, optionEnergy);
        colorTriggShade[i]  = GetDefaultTriggerColorName(triggerName[i], 1, optionEnergy);
        markerTrigg[i]      = GetDefaultTriggerMarkerStyleName(triggerName[i], 0);
        markerTriggMC[i]    = GetDefaultTriggerMarkerStyleName(triggerName[i], 1);
        sizeTrigg[i]        = GetDefaultTriggerMarkerSizeName(triggerName[i], 0);
        triggerNameLabel[i] = triggerName[i];
        triggerNameLabel[i].ReplaceAll("_"," ");
        triggerNameLabel[i].ReplaceAll("NLM","LM=");

        if ((optionEnergy.Contains("Pb") || optionEnergy.Contains("Xe")) && i == 0 ){
            fCent                                   = GetCentralityString(cutNumber[i]);
            fCentOutput                             = GetCentralityStringOutput(cutNumber[i]);
            collisionSystem                         = fCent+ " "+collisionSystem;
        }

    }

    //***************************************************************************************************************
    //*************************** set common binning ****************************************************************
    //***************************************************************************************************************
    Double_t binningEtaPrime[400];
    Int_t maxNBinsEtaPrimeAbs            = 0;
    Int_t maxNBinsEtaPrime               = GetBinning( binningEtaPrime, maxNBinsEtaPrimeAbs, "EtaPrime", optionEnergy, modeNormal, -1, kFALSE, fCent );
    Int_t maxNAllowedEtaPrime            = 0;
    Int_t nRealTriggers             = 0;
    cout << "binning etaprime" << endl;
    while (binningEtaPrime[maxNAllowedEtaPrime] < maxPtGlobalEtaPrime ) maxNAllowedEtaPrime++;
    for (Int_t i= 0; i< maxNAllowedEtaPrime+1; i++){
        cout << binningEtaPrime[i] << ", ";
    }
    cout << endl;

    //***************************************************************************************************************
    // defining output directory
    //***************************************************************************************************************
    TString outputDir =    Form("%s/%s/%s/FinalResultsTriggersPatched%s%s", suffix.Data(),optionEnergy.Data(),dateForOutput.Data(),fNLMStringOutput.Data(),system.Data());
    TString outputDirDay =    Form("%s/%s/%s", suffix.Data(),optionEnergy.Data(),dateForOutput.Data());
    if (optionEnergy.Contains("Pb") || optionEnergy.Contains("Xe")){
        outputDir       = outputDir+"_"+fCentOutput;
    }
    gSystem->Exec("mkdir -p "+outputDir);

    gSystem->Exec(Form("cp %s %s/configurationFile.txt", fileListNameEtaPrime.Data(), outputDir.Data()));

    TString nameFinalResDat                     = Form("%s/FitResults.dat",outputDir.Data());
    fstream  fileFitsOutput;
    fileFitsOutput.open(nameFinalResDat.Data(), ios::out);
    TString sysStringComb = "PCM";
    if(modeNormal == 2) sysStringComb = "PCMEMC";
    else if(modeNormal == 3) sysStringComb = "PCMPHOS";
    else if(modeNormal == 4) sysStringComb = "EMCEMC";
    else if(modeNormal == 5) sysStringComb = "PHOS";

    vector<TString>** ptSysDetail     = new vector<TString>*[MaxNumberOfFiles];
    for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++) ptSysDetail[iR] = new vector<TString>[400];
    //***************************************************************************************************************
    //******************************** Load EtaPrime histograms **********************************************************
    //***************************************************************************************************************
    TString FileNameCorrectedEtaPrime       [MaxNumberOfFiles];
    TFile*  fileCorrectedEtaPrime           [MaxNumberOfFiles];
    TString FileNameUnCorrectedEtaPrime     [MaxNumberOfFiles];
    TFile*  fileUnCorrectedEtaPrime         [MaxNumberOfFiles];
    TString FileNameUnCorrectedMCEtaPrime   [MaxNumberOfFiles];
    TFile*  fileUnCorrectedMCEtaPrime       [MaxNumberOfFiles];

    TH1D*   histoCorrectedYieldEtaPrime     [MaxNumberOfFiles];
    TH1D*   histoRawClusterPt               [MaxNumberOfFiles];
    TH1D*   histoRawClusterE                [MaxNumberOfFiles];
    TH1D*   histoEfficiencyEtaPrime         [MaxNumberOfFiles];
    TH1D*   histoAcceptanceEtaPrime         [MaxNumberOfFiles];

    TH1D*   histoAcceptanceEtaPrimeWOEvtWeights [MaxNumberOfFiles];
    TH1D*   histoEffTimesAccEtaPrime        [MaxNumberOfFiles];
    TH1D*   histoRawYieldEtaPrime           [MaxNumberOfFiles];
    TH1D*   histoRatioRawClusterPt          [MaxNumberOfFiles];
    TH1D*   histoRatioRawClusterE           [MaxNumberOfFiles];
    TH1D*   histoTriggerRejection           [MaxNumberOfFiles];
    TH1D*   histoMassEtaPrimeData           [MaxNumberOfFiles];
    TH1D*   histoMassEtaPrimeMC             [MaxNumberOfFiles];
    TH1D*   histoWidthEtaPrimeData          [MaxNumberOfFiles];
    TH1D*   histoWidthEtaPrimeMC            [MaxNumberOfFiles];
    TH1F*   histoEventQualtity              [MaxNumberOfFiles];
    TH1D*   histoMCInputEtaPrime            [MaxNumberOfFiles];

    Double_t triggRejecFac                  [MaxNumberOfFiles][MaxNumberOfFiles];
    Double_t triggRejecFacErr               [MaxNumberOfFiles][MaxNumberOfFiles];

    TString FileNameEffBaseEtaPrime         [MaxNumberOfFiles];
    TFile*  fileEffBaseEtaPrime             [MaxNumberOfFiles];
    TH1D*   histoEffBaseEtaPrime            [MaxNumberOfFiles];
    TH1D*   histoTriggerEffEtaPrime         [MaxNumberOfFiles];
    Bool_t  enableTriggerEffEtaPrime        [MaxNumberOfFiles];
    Bool_t  enableTriggerEffEtaPrimeAll                     = kFALSE;
    Bool_t  enableTriggerRejecCompMC                        = kFALSE;

    TH1D*   histoMCRawClusterPt             [MaxNumberOfFiles];
    TH1D*   histoMCRatioRawClusterPt        [MaxNumberOfFiles];
    TString rapidityRange                                   = "";
    Double_t deltaRapid                     [MaxNumberOfFiles];

    TH1D*   histoInvMassSigPlusBG           [MaxNumberOfFiles];
    TH1D*   histoInvMassSig                 [MaxNumberOfFiles];
    TH1D*   histoInvMassBG                  [MaxNumberOfFiles];
    TF1*    fitInvMassSig                   [MaxNumberOfFiles];


    // check if fit file for binshifting has to be adjusted for every energy
    TF1* fitBinShiftEtaPrimeTCM                          = 0x0;
    TF1* fitBinShiftEtaPrime                             = 0x0;

    for (Int_t i=0; i< nrOfTrigToBeComb; i++){
        // Define CutSelections
        TString fEventCutSelection                          = "";
        TString fGammaCutSelection                          = "";
        TString fClusterCutSelection                        = "";
        TString fElectronCutSelection                       = "";
        TString fMesonCutSelection                          = "";
        // disentangle cut selection
        ReturnSeparatedCutNumberAdvanced(cutNumber[i].Data(),fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);

        TString trigger                                     = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
        Double_t scaleFacSinBin                             = 1.0;
        Int_t exampleBin                                    = ReturnSingleInvariantMassBinPlotting ("EtaPrime", optionEnergy, modeNormal, trigger.Atoi(), scaleFacSinBin);

        FileNameCorrectedEtaPrime[i]                             = Form("%s/%s/EtaPrime_%s_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(),
                                                                    cutNumber[i].Data());
        cout<< FileNameCorrectedEtaPrime[i] << endl;
        fileCorrectedEtaPrime[i]                                 = new TFile(FileNameCorrectedEtaPrime[i]);
        if (fileCorrectedEtaPrime[i]->IsZombie()) return;
        // read uncorrected file
        FileNameUnCorrectedEtaPrime[i]                           = Form("%s/%s/EtaPrime_%s_GammaConvV1WithoutCorrection_%s.root",cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(),
                                                                    cutNumber[i].Data());

        cout<< FileNameUnCorrectedEtaPrime[i] << endl;
        fileUnCorrectedEtaPrime[i]                               = new TFile(FileNameUnCorrectedEtaPrime[i]);
        if (fileUnCorrectedEtaPrime[i]->IsZombie()) return;

        if (isMC.CompareTo("data") == 0){
            FileNameUnCorrectedMCEtaPrime[i]                     = Form("%s/%s/EtaPrime_MC_GammaConvV1WithoutCorrection_%s.root",cutNumber[i].Data(), optionEnergy.Data(), cutNumber[i].Data());

            cout<< FileNameUnCorrectedMCEtaPrime[i] << endl;
            fileUnCorrectedMCEtaPrime[i]                         = new TFile(FileNameUnCorrectedMCEtaPrime[i]);
            if (fileUnCorrectedMCEtaPrime[i]->IsZombie())
                enableTriggerRejecCompMC                    = kFALSE;
            else
                enableTriggerRejecCompMC                    = kTRUE;
        }
        deltaRapid[i]                                       =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
        cout << "For Trigger: " << triggerName[i] << " using rapidity: " <<  rapidityRange.Data() << endl;

        histoEventQualtity[i]                               = (TH1F*)fileCorrectedEtaPrime[i]->Get("NEvents");
        histoCorrectedYieldEtaPrime[i]                           = (TH1D*)fileCorrectedEtaPrime[i]->Get(nameCorrectedYield.Data());
        histoCorrectedYieldEtaPrime[i]->SetName(Form("CorrectedYield_%s",cutNumber[i].Data()));
        histoEfficiencyEtaPrime[i]                               = (TH1D*)fileCorrectedEtaPrime[i]->Get(nameEfficiency.Data());
        histoEfficiencyEtaPrime[i]->SetName(Form("Efficiency_%s",  cutNumber[i].Data()));
        histoAcceptanceEtaPrime[i]                               = (TH1D*)fileCorrectedEtaPrime[i]->Get(nameAcceptance.Data());
        histoAcceptanceEtaPrime[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
        histoAcceptanceEtaPrimeWOEvtWeights[i]                   = (TH1D*)fileCorrectedEtaPrime[i]->Get(nameAcceptanceWOEvtWeights.Data());


        if(histoAcceptanceEtaPrimeWOEvtWeights[i]) histoAcceptanceEtaPrimeWOEvtWeights[i]->SetName(Form("AcceptanceWOEvtWeights_%s",  cutNumber[i].Data()));

        histoRawYieldEtaPrime[i]                                 = (TH1D*)fileCorrectedEtaPrime[i]->Get("histoYieldMesonPerEvent");
        histoRawYieldEtaPrime[i]->SetName(Form("RAWYieldPerEvent_%s",cutNumber[i].Data()));

        histoMassEtaPrimeData[i]                                 = (TH1D*)fileCorrectedEtaPrime[i]->Get("histoMassMeson");
        histoMassEtaPrimeData[i]->SetName(Form("EtaPrime_Mass_data_%s",cutNumber[i].Data()));
        histoMassEtaPrimeMC[i]                                   = (TH1D*)fileCorrectedEtaPrime[i]->Get(nameMassMC.Data());
        histoMassEtaPrimeMC[i]->SetName(Form("EtaPrime_Mass_MC_%s",cutNumber[i].Data()));
        histoWidthEtaPrimeData[i]                                = (TH1D*)fileCorrectedEtaPrime[i]->Get("histoFWHMMeson");
        histoWidthEtaPrimeData[i]->SetName(Form("EtaPrime_Width_data_%s",cutNumber[i].Data()));
        histoWidthEtaPrimeMC[i]                                  = (TH1D*)fileCorrectedEtaPrime[i]->Get(nameWidthMC.Data());
        histoWidthEtaPrimeMC[i]->SetName(Form("EtaPrime_Width_MC_%s",cutNumber[i].Data()));
        histoInvMassSig[i]                                  = (TH1D*)fileCorrectedEtaPrime[i]->Get(Form("InvMassSig_PtBin%02d",exampleBin));
        if (histoInvMassSig[i]) histoInvMassSig[i]->SetName(Form("EtaPrime_InvMassSig_Example_%s",triggerName[i].Data()));
        histoInvMassSigPlusBG[i]                            = (TH1D*)fileCorrectedEtaPrime[i]->Get(Form("InvMassSigPlusBG_PtBin%02d",exampleBin));
        if (histoInvMassSigPlusBG[i]) histoInvMassSigPlusBG[i]->SetName(Form("EtaPrime_InvMassSigPlusBG_Example_%s",triggerName[i].Data()));
        histoInvMassBG[i]                                   = (TH1D*)fileCorrectedEtaPrime[i]->Get(Form("InvMassBG_PtBin%02d",exampleBin));
        if (histoInvMassBG[i]) histoInvMassBG[i]->SetName(Form("EtaPrime_InvMassBG_Example_%s",triggerName[i].Data()));
        fitInvMassSig[i]                                    = (TF1*)fileCorrectedEtaPrime[i]->Get(Form("FitInvMassSig_PtBin%02d",exampleBin));
        if (fitInvMassSig[i]) fitInvMassSig[i]->SetName(Form("EtaPrime_InvMassSigFit_Example_%s",triggerName[i].Data()));
        if (cutNumberBaseEff[i].CompareTo("bla") != 0){
            FileNameEffBaseEtaPrime[i]                           = Form("%s/%s/EtaPrime_MC_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), cutNumberBaseEff[i].Data());
            fileEffBaseEtaPrime[i]                               = new TFile(FileNameEffBaseEtaPrime[i]);
            if (fileEffBaseEtaPrime[i]->IsZombie()){
                enableTriggerEffEtaPrime[i]                         = kFALSE;
                cout << "Didn't find the effi base file for " << triggerName[i].Data() << endl;
                cout << "ABORTING: as base effi file was requested" << endl;

            } else {
                enableTriggerEffEtaPrime[i]                         = kTRUE;
                enableTriggerEffEtaPrimeAll                         = kTRUE;
            }
            if (enableTriggerEffEtaPrime[i]){
                TString effiNameBase                        = "TrueMesonEffiPt";
                TH1D* histoEffiEtaPrimeTemp                      = (TH1D*)fileCorrectedEtaPrime[i]->Get(effiNameBase.Data());
                TH1D* histoEffiBaseEtaPrimeTemp                  = (TH1D*)fileEffBaseEtaPrime[i]->Get(effiNameBase.Data());
                histoEffBaseEtaPrime[i]                          = (TH1D*)fileEffBaseEtaPrime[i]->Get(nameEfficiency.Data());
                histoEffBaseEtaPrime[i]->SetName(Form("EfficiencyBase_%s",  cutNumber[i].Data()));
                histoTriggerEffEtaPrime[i]                       = (TH1D*)histoEffiEtaPrimeTemp->Clone(Form("TriggerEfficiency_%s", cutNumber[i].Data()));
                histoTriggerEffEtaPrime[i]->Divide(histoTriggerEffEtaPrime[i],histoEffiBaseEtaPrimeTemp,1.,1.,"B");

                // limit trigger efficiency to 1
                if(optionEnergy.CompareTo("8TeV") == 0){
                  for(Int_t j = 1; j<histoTriggerEffEtaPrime[i]->GetNbinsX()+1; j++){
                    Double_t binC = histoTriggerEffEtaPrime[i]->GetBinContent(j);
                    if(binC > 1.){
                      histoEffBaseEtaPrime[i]->SetBinContent(j, histoEffBaseEtaPrime[i]->GetBinContent(j)*binC);
                      histoTriggerEffEtaPrime[i]->SetBinContent(j, 1.);
                    }
                  }
                }

                histoEffTimesAccEtaPrime[i]                      = (TH1D*)histoEffBaseEtaPrime[i]->Clone(Form("EffTimeAcc_%s",  cutNumber[i].Data()));
                if(histoAcceptanceEtaPrimeWOEvtWeights[i]){
                  histoAcceptanceEtaPrimeWOEvtWeights[i]->Sumw2();
                  histoEffTimesAccEtaPrime[i]->Multiply(histoAcceptanceEtaPrimeWOEvtWeights[i]);
                  histoEfficiencyEtaPrime[i]->Multiply(histoAcceptanceEtaPrime[i]);
                  histoEfficiencyEtaPrime[i]->Divide(histoEfficiencyEtaPrime[i],histoAcceptanceEtaPrimeWOEvtWeights[i],1.,1.,"B");
                  histoAcceptanceEtaPrime[i] = histoAcceptanceEtaPrimeWOEvtWeights[i];
                  histoAcceptanceEtaPrime[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
                }else histoEffTimesAccEtaPrime[i]->Multiply(histoAcceptanceEtaPrime[i]);
                histoEffTimesAccEtaPrime[i]->Scale(deltaRapid[i]*2*TMath::Pi());
            } else {
                histoEffBaseEtaPrime[i]                          = NULL;
                histoTriggerEffEtaPrime[i]                       = NULL;
            }
        } else {
            histoEffTimesAccEtaPrime[i]                      = (TH1D*)histoEfficiencyEtaPrime[i]->Clone(Form("EffTimeAcc_%s",  cutNumber[i].Data()));
            if(histoAcceptanceEtaPrimeWOEvtWeights[i]){
              histoAcceptanceEtaPrimeWOEvtWeights[i]->Sumw2();
              histoEffTimesAccEtaPrime[i]->Multiply(histoAcceptanceEtaPrimeWOEvtWeights[i]);
              histoEfficiencyEtaPrime[i]->Multiply(histoAcceptanceEtaPrime[i]);
              histoEfficiencyEtaPrime[i]->Divide(histoEfficiencyEtaPrime[i],histoAcceptanceEtaPrimeWOEvtWeights[i],1.,1.,"B");
              histoAcceptanceEtaPrime[i] = histoAcceptanceEtaPrimeWOEvtWeights[i];
              histoAcceptanceEtaPrime[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
            }else histoEffTimesAccEtaPrime[i]->Multiply(histoAcceptanceEtaPrime[i]);
            histoEffTimesAccEtaPrime[i]->Scale(deltaRapid[i]*2*TMath::Pi());
        }

        //Scale spectrum to MBOR
        if (optionEnergy.CompareTo("2.76TeV")==0 &&
            (triggerName[i].Contains("INT7")|| triggerName[i].Contains("EMC7") || triggerName[i].Contains("EG1") || triggerName[i].Contains("EG2")) &&
            isMC.CompareTo("data") == 0){
            histoCorrectedYieldEtaPrime[i]->Scale(0.8613) ;
            histoRawYieldEtaPrime[i]->Scale(0.8613) ;
        }
        if (triggerName[i].CompareTo("INT7") != 0 && triggerName[i].CompareTo("MB") != 0 && triggerName[i].CompareTo("INT1") != 0){
            nRealTriggers++;
        }

        histoMCInputEtaPrime[i]                                  = (TH1D*)fileCorrectedEtaPrime[i]->Get(nameMCYield.Data());
        histoMCInputEtaPrime[i]->SetName(Form("EtaPrime_Input_Reweighted_%s",cutNumber[i].Data()));
        histoMCInputEtaPrime[i]->Sumw2();
        histoMCInputEtaPrime[i]->Rebin(4);
        histoMCInputEtaPrime[i]->Scale(1./4);
        //***************************************************************************************************************
        //****************************** Calculate trigger rejection factors ********************************************
        //***************************************************************************************************************
        if ((modeNormal == 0 || modeNormal == 2 || modeNormal == 3 || modeNormal == 4 || modeNormal == 5 ) && hasClusterOutput ){
            histoRawClusterPt[i]                        = (TH1D*)fileUnCorrectedEtaPrime[i]->Get("ClusterPtPerEvent");
            histoRawClusterE[i]                         = (TH1D*)fileUnCorrectedEtaPrime[i]->Get("ClusterEPerEvent");
            if (!histoRawClusterPt[i]){
                cout << "INFO: couldn't find cluster input, disabeling it!" << endl;
                hasClusterOutput                        = kFALSE;
                triggRejecFac[i][trigSteps[i][0]]       = 1;
                triggRejecFacErr[i][trigSteps[i][0]]    = 0;
            } else {
                histoRawClusterPt[i]->SetName(Form("ClusterPtPerEvent_%s",cutNumber[i].Data()));
                histoRatioRawClusterPt[i]                   = (TH1D*)histoRawClusterPt[i]->Clone(Form("RatioCluster_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()));
                histoRatioRawClusterPt[i]->Divide(histoRatioRawClusterPt[i],histoRawClusterPt[trigSteps[i][0]],1.,1.,"");

                if (histoRawClusterE[i]){
                    histoRawClusterE[i]->SetName(Form("ClusterEPerEvent_%s",cutNumber[i].Data()));
                    histoRatioRawClusterE[i]                   = (TH1D*)histoRawClusterE[i]->Clone(Form("RatioCluster_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()));
                    histoRatioRawClusterE[i]->Divide(histoRatioRawClusterE[i],histoRawClusterE[trigSteps[i][0]],1.,1.,"");
                }

                Int_t binMinTrigg                           = histoRatioRawClusterPt[i]->FindBin(minPt[i]);
                minPt[i]                                    = histoRatioRawClusterPt[i]->GetBinCenter(binMinTrigg);

                TF1* pol0                                   = new TF1("pol0","[0]",minPt[i],maxPt[i]); //
                histoRatioRawClusterPt[i]->Fit(pol0,"QNRMEX0+","",minPt[i],maxPt[i]);

                histoTriggerRejection[i]                    = new TH1D (Form("triggRejectMean_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),
                                                                        Form("triggRejectMean_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),
                                                                        10, 0.5,10.5 );

                for (Int_t k = 0; k<11; k++ ){
                    histoRatioRawClusterPt[i]->Fit(pol0,"QNRMEX0+","",histoRatioRawClusterPt[i]->GetBinCenter(binMinTrigg-2+k),maxPt[i]);
                    histoTriggerRejection[i]->SetBinContent(k+1,pol0->GetParameter(0));
                    histoTriggerRejection[i]->SetBinError(k+1,pol0->GetParError(0));
                }

                TF1* pol0_2                                 = new TF1("pol0_2","[0]",0.5,10.5); //
                histoTriggerRejection[i]->Fit(pol0_2,"QNRMEX0+","",0.5,10.5);
                triggRejecFac[i][trigSteps[i][0]]           = pol0_2->GetParameter(0);
                Double_t largestDev = 0;
                for (Int_t  k = 0; k<11; k++ ){
                    Double_t diffToFit                      = abs(histoTriggerRejection[i]->GetBinContent(k+1)-triggRejecFac[i][trigSteps[i][0]])+histoTriggerRejection[i]->GetBinError(k+1);
                    if (diffToFit > largestDev) largestDev  = diffToFit;
                }
                triggRejecFacErr[i][trigSteps[i][0]] = largestDev;


                delete pol0;
                delete pol0_2;
                cout << "trigger rejection factor " << triggerName[i].Data() << "/" << triggerName[trigSteps[i][0]].Data() << ": " << triggRejecFac[i][trigSteps[i][0]]
                    << "+-" << triggRejecFacErr[i][trigSteps[i][0]] << endl;

                if (enableTriggerRejecCompMC){
                    histoMCRawClusterPt[i]                  = (TH1D*)fileUnCorrectedMCEtaPrime[i]->Get("ClusterPtPerEvent");
                    histoMCRawClusterPt[i]->SetName(Form("MCClusterPtPerEvent_%s",cutNumber[i].Data()));
                    histoMCRatioRawClusterPt[i]             = (TH1D*)histoMCRawClusterPt[i]->Clone(Form("MCRatioCluster_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()));
                    histoMCRatioRawClusterPt[i]->Sumw2();
                    histoMCRatioRawClusterPt[i]->Divide(histoMCRatioRawClusterPt[i],histoMCRawClusterPt[trigSteps[i][0]],1.,1.,"");

                    TF1* pol0_MC                            = new TF1("pol0_MC","[0]",minPt[i],maxPt[i]); //
                    histoMCRatioRawClusterPt[i]->Fit(pol0_MC,"QNRMEX0+","",minPt[i],maxPt[i]);

                    Double_t scaleFactorMC                    = triggRejecFac[i][trigSteps[i][0]]/pol0_MC->GetParameter(0);
                    histoMCRatioRawClusterPt[i]->Scale(scaleFactorMC);

                    cout << "data: "<<triggRejecFac[i][trigSteps[i][0]] << "\t MC: " << pol0_MC->GetParameter(0) << "\t scale factor: " << scaleFactorMC<< endl;
                }
            }
        } else {
            triggRejecFac[i][trigSteps[i][0]]       = 1;
            triggRejecFacErr[i][trigSteps[i][0]]    = 0;

        }
    }

    // figure out which trigger we are using and set correct x-sections
    Int_t isV0AND           = 0;
    if (histoEventQualtity[0]->GetNbinsX() > 7){
        if (histoEventQualtity[0]->GetBinContent(9) > 0){
            isV0AND         = 1;
        }
    }
    if (optionEnergy.CompareTo("5TeV") == 0 || optionEnergy.CompareTo("5TeV2017") == 0 || optionEnergy.CompareTo("8TeV") == 0 ){
        isV0AND             = 1;
    }
    Double_t xSection       = ReturnCorrectXSection( optionEnergy, isV0AND);
    //***************************************************************************************************************
    //*******************************Plotting trigger rejection factors = fits log scale all in one *****************
    //***************************************************************************************************************
    if (hasClusterOutput){
        Size_t textSizeSpectra2         = 0.0415;
        Int_t textPixelPP               = textSizeSpectra2*1100;
        TCanvas* canvasTriggerReject    = new TCanvas("canvasTriggerReject","",0,0,1500,1100);// gives the page size
        DrawGammaCanvasSettings( canvasTriggerReject, 0.076, 0.015, 0.015, 0.085);
        canvasTriggerReject->SetLogy(1);

        Double_t minTriggReject = 0.1;
        Double_t maxTriggReject = 4200;
        if (optionEnergy.CompareTo("2.76TeV") == 0)
            maxTriggReject = 8200;
        else if (modeNormal == 4 && optionEnergy.CompareTo("pPb_5.023TeV") == 0)
            maxTriggReject = 200;
        else if (optionEnergy.CompareTo("13TeV") == 0)
            maxTriggReject = 4200;

        TH2F * histo2DTriggReject;
        histo2DTriggReject = new TH2F("histo2DTriggReject","histo2DTriggReject",1000,0., maxPtGlobalCluster,10000,minTriggReject, maxTriggReject);
        SetStyleHistoTH2ForGraphs(histo2DTriggReject, "#it{p}_{T} (GeV/#it{c})","#it{R}", //"#frac{N_{clus,trig A}/N_{Evt, trig A}}{N_{clus,trig B}/N_{Evt,trig B}}",
                                0.85*textSizeSpectra2,textSizeSpectra2, 0.85*textSizeSpectra2,textSizeSpectra2, 0.85,0.85);
        histo2DTriggReject->DrawCopy();

        TLegend* legendTriggReject = GetAndSetLegend2(0.33, 0.12, 0.92, 0.12+(0.9*(nrOfTrigToBeComb-2+1)*textSizeSpectra2),textPixelPP);
        legendTriggReject->SetMargin(0.02);
        legendTriggReject->SetNColumns(3);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (i != trigSteps[i][0] ){
                for (Int_t j = 1; j < histoRatioRawClusterPt[i]->GetNbinsX()+1; j++){
                    if (histoRatioRawClusterPt[i]->GetBinError(j)/histoRatioRawClusterPt[i]->GetBinContent(j) > 1){
                            histoRatioRawClusterPt[i]->SetBinContent(j,1e6);
                            histoRatioRawClusterPt[i]->SetBinError(j,0);
                    }
                }

                DrawGammaSetMarker(histoRatioRawClusterPt[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoRatioRawClusterPt[i]->DrawCopy("e,same");
                legendTriggReject->AddEntry(histoRatioRawClusterPt[i],Form("   %s/%s",triggerNameLabel[i].Data(),triggerNameLabel[trigSteps[i][0]].Data() ),"p");
                legendTriggReject->AddEntry((TObject*)0,Form("  %3.1f < #it{p}_{T} < %3.1f",minPt[i],maxPt[i]),"");
                if (triggerName[i].Contains("EMC1"))
                    legendTriggReject->AddEntry((TObject*)0,Form("     %3.0f #pm %3.0f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else if (triggerName[i].Contains("EMC7"))
                    legendTriggReject->AddEntry((TObject*)0,Form("     %3.1f #pm %3.1f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else if (triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA"))
                    legendTriggReject->AddEntry((TObject*)0,Form("     %3.2f #pm %3.2f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else if (triggerName[i].Contains("EG1") )
                    legendTriggReject->AddEntry((TObject*)0,Form("     %3.3f #pm %3.3f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else
                    legendTriggReject->AddEntry((TObject*)0,Form("     %3.2f #pm %3.2f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");


                TF1* pol0 = new TF1("pol0","[0]",minPt[i],maxPt[i]); //
                histoRatioRawClusterPt[i]->Fit(pol0,"NRME+","",minPt[i],maxPt[i]);

                fileFitsOutput << triggerNameLabel[i].Data() << "-" << triggerNameLabel[trigSteps[i][0]].Data() << endl;
                fileFitsOutput << WriteParameterToFile(pol0) << endl;

                TH1D* triggRejecCLPol0 = (TH1D*)histoRatioRawClusterPt[i]->Clone(Form("CL_%i",i));
                for (Int_t j = 1; j < triggRejecCLPol0->GetNbinsX()+1; j++){
                    triggRejecCLPol0->SetBinContent(j,triggRejecFac[i][trigSteps[i][0]]);
                    triggRejecCLPol0->SetBinError(j,triggRejecFacErr[i][trigSteps[i][0]]);
                }
                triggRejecCLPol0->SetStats(kFALSE);
                triggRejecCLPol0->SetFillColor(colorTriggShade[i]);
                triggRejecCLPol0->SetMarkerSize(0);
                triggRejecCLPol0->Draw("e3,same");

                pol0->SetParameter(0,triggRejecFac[i][trigSteps[i][0]]);
                pol0->SetParError(0,triggRejecFacErr[i][trigSteps[i][0]]);
                pol0->SetLineColor(colorTrigg[i]-1);
                pol0->SetLineStyle(7);
                pol0->SetRange(minPt[i],maxPt[i]);
                pol0->Draw("same");
                histoRatioRawClusterPt[i]->DrawCopy("e,same");

                TF1* pol1 = new TF1("pol1","[0]+[1]*x",minPt[i],maxPt[i]); //
                histoRatioRawClusterPt[i]->Fit(pol1,"NRME0+","",minPt[i],maxPt[i]);
                fileFitsOutput << WriteParameterToFile(pol1) << endl;

            }
        }
        legendTriggReject->Draw();
        histo2DTriggReject->Draw("same,axis");
        TLatex *labelPerfTriggRejec  = NULL;
        if (isMC.CompareTo("MC") == 0)
            labelPerfTriggRejec  = new TLatex(0.11, 0.925,"ALICE simulation");
        else
            labelPerfTriggRejec  = new TLatex(0.11, 0.925,"ALICE performance");

        SetStyleTLatex( labelPerfTriggRejec, textSizeSpectra2,4);
        labelPerfTriggRejec->Draw();

        TLatex *labelEnergyTriggRejec  = new TLatex(0.95, 0.925,collisionSystem.Data());
        SetStyleTLatex( labelEnergyTriggRejec, textSizeSpectra2,4, 1, 42, kTRUE, 31);
        labelEnergyTriggRejec->Draw();

        TLatex *labelPerfTriggFitRange = new TLatex(0.523, 0.12+(0.9*(nrOfTrigToBeComb-2+1)*textSizeSpectra2)+0.01, "Fit range (GeV/#it{c})");
        SetStyleTLatex( labelPerfTriggFitRange, textSizeSpectra2,4);
        labelPerfTriggFitRange->Draw();

        TLatex *labelPerfTriggRejecFac = new TLatex(0.753, 0.12+(0.9*(nrOfTrigToBeComb-2+1)*textSizeSpectra2)+0.01, "Trigger rejection");
        SetStyleTLatex( labelPerfTriggRejecFac, textSizeSpectra2,4);
        labelPerfTriggRejecFac->Draw();

        canvasTriggerReject->Update();
        canvasTriggerReject->SaveAs(Form("%s/%s_TriggerRejectionFactors.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        if (histoRatioRawClusterE[0]){
            //***************************************************************************************************************
            //*******************************Plotting trigger rejection factors = fits log scale all in one vs E ************
            //***************************************************************************************************************

            TH2F * histo2DTriggRejectE;
            histo2DTriggRejectE = new TH2F("histo2DTriggRejectE","histo2DTriggRejectE",1000,0., maxPtGlobalCluster,10000,minTriggReject, maxTriggReject);
            SetStyleHistoTH2ForGraphs(histo2DTriggRejectE, "#it{E} (GeV)","#it{R}", //"#frac{N_{clus,trig A}/N_{Evt, trig A}}{N_{clus,trig B}/N_{Evt,trig B}}",
                                    0.85*textSizeSpectra2,textSizeSpectra2, 0.85*textSizeSpectra2,textSizeSpectra2, 0.85,0.85);
            histo2DTriggRejectE->DrawCopy();

            TLegend* legendTriggRejectE = GetAndSetLegend2(0.33, 0.12, 0.92, 0.12+(0.9*(nrOfTrigToBeComb-2+1)*textSizeSpectra2),textPixelPP);
            legendTriggRejectE->SetMargin(0.02);
            legendTriggRejectE->SetNColumns(3);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if (i != trigSteps[i][0] ){
                    for (Int_t j = 1; j < histoRatioRawClusterE[i]->GetNbinsX()+1; j++){
                        if (histoRatioRawClusterE[i]->GetBinError(j)/histoRatioRawClusterE[i]->GetBinContent(j) > 1){
                                histoRatioRawClusterE[i]->SetBinContent(j,1e6);
                                histoRatioRawClusterE[i]->SetBinError(j,0);
                        }
                    }

                    DrawGammaSetMarker(histoRatioRawClusterE[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoRatioRawClusterE[i]->DrawCopy("e,same");

                    TF1* pol0 = new TF1("pol0","[0]",minPt[i],maxPt[i]); //
                    histoRatioRawClusterE[i]->Fit(pol0,"NRME+","",minPt[i],maxPt[i]);
                    legendTriggRejectE->AddEntry(histoRatioRawClusterPt[i],Form("   %s/%s",triggerNameLabel[i].Data(),triggerNameLabel[trigSteps[i][0]].Data() ),"p");
                    legendTriggRejectE->AddEntry((TObject*)0,Form("  %3.1f < #it{E} < %3.1f",minPt[i],maxPt[i]),"");
                    if (triggerName[i].Contains("EMC1"))
                        legendTriggRejectE->AddEntry((TObject*)0,Form("     %3.0f #pm %3.0f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                    else if (triggerName[i].Contains("EMC7"))
                        legendTriggRejectE->AddEntry((TObject*)0,Form("     %3.1f #pm %3.1f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                    else if (triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA"))
                        legendTriggRejectE->AddEntry((TObject*)0,Form("     %3.2f #pm %3.2f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                    else if (triggerName[i].Contains("EG1") )
                        legendTriggRejectE->AddEntry((TObject*)0,Form("     %3.3f #pm %3.3f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                    else
                        legendTriggRejectE->AddEntry((TObject*)0,Form("     %3.2f #pm %3.2f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");

                    TH1D* triggRejecCLPol0 = (TH1D*)histoRatioRawClusterE[i]->Clone(Form("CL_%i",i));
                    for (Int_t j = 1; j < triggRejecCLPol0->GetNbinsX()+1; j++){
                        triggRejecCLPol0->SetBinContent(j,triggRejecFac[i][trigSteps[i][0]]);
                        triggRejecCLPol0->SetBinError(j,triggRejecFacErr[i][trigSteps[i][0]]);
                    }
                    triggRejecCLPol0->SetStats(kFALSE);
                    triggRejecCLPol0->SetFillColor(colorTriggShade[i]);
                    triggRejecCLPol0->SetMarkerSize(0);
                    triggRejecCLPol0->Draw("e3,same");

                    pol0->SetParameter(0,triggRejecFac[i][trigSteps[i][0]]);
                    pol0->SetParError(0,triggRejecFacErr[i][trigSteps[i][0]]);
                    pol0->SetLineColor(colorTrigg[i]-1);
                    pol0->SetLineStyle(7);
                    pol0->SetRange(minPt[i],maxPt[i]);
                    pol0->Draw("same");
                    histoRatioRawClusterE[i]->DrawCopy("e,same");
                }
            }
            legendTriggRejectE->Draw();
            histo2DTriggReject->Draw("same,axis");
            labelPerfTriggRejec->Draw();
            labelEnergyTriggRejec->Draw();

            TLatex *labelPerfTriggFitRangeE = new TLatex(0.523, 0.12+(0.9*(nrOfTrigToBeComb-2+1)*textSizeSpectra2)+0.01, "Fit range (GeV)");
            SetStyleTLatex( labelPerfTriggFitRangeE, textSizeSpectra2,4);
            labelPerfTriggFitRangeE->Draw();
            labelPerfTriggRejecFac->Draw();

            canvasTriggerReject->Update();
            canvasTriggerReject->SaveAs(Form("%s/%s_TriggerRejectionFactors_E.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        delete canvasTriggerReject;

        //***************************************************************************************************************
        //**************************Plotting trigger rejection factors = fits linear scale all in one *******************
        //***************************************************************************************************************

        TCanvas* canvasTriggerRejectLinear = new TCanvas("canvasTriggerReject","",0,0,1500,1100);// gives the page size
        DrawGammaCanvasSettings( canvasTriggerRejectLinear, 0.076, 0.015, textSizeSpectra2, 0.08);
        canvasTriggerRejectLinear->SetLogy(0);

        Double_t minTriggRejectLin = 0;
        Double_t maxTriggRejectLin = 2000;
        if (isMC.CompareTo("MC")== 0) maxTriggRejectLin = 2500;
        if (modeNormal == 4 && optionEnergy.CompareTo("pPb_5.023TeV") == 0)
            maxTriggRejectLin = 100;
        if (modeNormal == 4 && optionEnergy.CompareTo("2.76TeV") == 0)
            maxTriggRejectLin = 2000;

        if( optionEnergy.CompareTo("8TeV")==0 ){
            if (modeNormal == 2 || modeNormal == 4  )
                maxTriggRejectLin = 310;
        } else if( optionEnergy.CompareTo("13TeV")==0 ){
            if (modeNormal == 2 || modeNormal == 4  )
                maxTriggRejectLin = 1005;
        }
        TH2F * histo2DTriggRejectLinear;
        histo2DTriggRejectLinear = new TH2F("histo2DTriggRejectLinear","histo2DTriggRejectLinear",1000,0., maxPtGlobalCluster,15000,minTriggRejectLin, maxTriggRejectLin);
        SetStyleHistoTH2ForGraphs(histo2DTriggRejectLinear, "#it{p}_{T} (GeV/#it{c})","#it{R}", //"#frac{N_{clus,trig A}/N_{Evt, trig A}}{N_{clus,trig B}/N_{Evt,trig B}}",
                                0.85*textSizeSpectra2,textSizeSpectra2, 0.85*textSizeSpectra2,textSizeSpectra2, 0.85,0.85);
        histo2DTriggRejectLinear->DrawCopy();

        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (i != trigSteps[i][0] ){
                DrawGammaSetMarker(histoRatioRawClusterPt[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoRatioRawClusterPt[i]->DrawCopy("e1,same");
                TF1* pol0 = new TF1("pol0","[0]",minPt[i],maxPt[i]); //

                histoRatioRawClusterPt[i]->Fit(pol0,"NRME+","",minPt[i],maxPt[i]);
                TH1D* triggRejecCLPol0 = (TH1D*)histoRatioRawClusterPt[i]->Clone(Form("CL_%i",i));
                triggRejecCLPol0->SetStats(kFALSE);
                triggRejecCLPol0->SetFillColor(colorTriggShade[i]);
                triggRejecCLPol0->SetMarkerSize(0);
                for (Int_t j = 1; j < triggRejecCLPol0->GetNbinsX()+1; j++){
                    triggRejecCLPol0->SetBinContent(j,triggRejecFac[i][trigSteps[i][0]]);
                    triggRejecCLPol0->SetBinError(j,triggRejecFacErr[i][trigSteps[i][0]]);
                }
                triggRejecCLPol0->Draw("e3,same");

                pol0->SetParameter(0,triggRejecFac[i][trigSteps[i][0]]);
                pol0->SetParError(0,triggRejecFacErr[i][trigSteps[i][0]]);
                pol0->SetLineColor(colorTrigg[i]-1);
                pol0->SetLineStyle(7);
                pol0->SetRange(minPt[i],maxPt[i]);
                pol0->Draw("same");
                histoRatioRawClusterPt[i]->DrawCopy("e1,same");
            }
        }
        legendTriggReject->Draw();
        labelPerfTriggFitRange->Draw();
        labelPerfTriggRejecFac->Draw();

        histo2DTriggRejectLinear->Draw("same,axis");

        canvasTriggerRejectLinear->Update();
        canvasTriggerRejectLinear->SaveAs(Form("%s/%s_TriggerRejectionFactorsLinY.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        //***************************************************************************************************************
        //********************Plotting trigger rejection factors = fits, linear scale 1 trigger at a time ***************
        //***************************************************************************************************************

        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (i != trigSteps[i][0] ){
                histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,histoRatioRawClusterPt[i]->GetMaximum()*1.5);
                if (triggerName[i].Contains("EMC7"))
                    histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,250);
                if (triggerName[i].Contains("EG2"))
                    histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,30);
                if (triggerName[i].Contains("EG1"))
                    histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,8);
                if (optionEnergy.CompareTo("8TeV")==0){
                    if (triggerName[i].Contains("EMC7"))
                        histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,150);
                    if (triggerName[i].Contains("EGA"))
                        histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,300);
                } else if (optionEnergy.CompareTo("13TeV")==0){
                    if (triggerName[i].Contains("EG2"))
                        histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,1000);
                    if (triggerName[i].Contains("EG1"))
                        histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,30);
                }


                histo2DTriggRejectLinear->DrawCopy();

                DrawGammaSetMarker(histoRatioRawClusterPt[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoRatioRawClusterPt[i]->DrawCopy("e1,same");
                TF1* pol0 = new TF1("pol0","[0]",minPt[i],maxPt[i]); //

                histoRatioRawClusterPt[i]->Fit(pol0,"NRME+","",minPt[i],maxPt[i]);
                TH1D* triggRejecCLPol0 = (TH1D*)histoRatioRawClusterPt[i]->Clone(Form("CL_%i",i));
                for (Int_t j = 1; j < triggRejecCLPol0->GetNbinsX()+1; j++){
                    triggRejecCLPol0->SetBinContent(j,triggRejecFac[i][trigSteps[i][0]]);
                    triggRejecCLPol0->SetBinError(j,triggRejecFacErr[i][trigSteps[i][0]]);
                }

                triggRejecCLPol0->SetStats(kFALSE);
                triggRejecCLPol0->SetFillColor(colorTriggShade[i]);
                triggRejecCLPol0->SetMarkerSize(0);
                triggRejecCLPol0->Draw("e3,same");

                pol0->SetParameter(0,triggRejecFac[i][trigSteps[i][0]]);
                pol0->SetParError(0,triggRejecFacErr[i][trigSteps[i][0]]);
                pol0->SetLineColor(colorTrigg[i]-1);
                pol0->SetLineStyle(7);
                pol0->SetRange(minPt[i],maxPt[i]);
                pol0->Draw("same");
                histoRatioRawClusterPt[i]->DrawCopy("e1,same");

                TF1* pol1 = new TF1("pol1","[0]+[1]*x",minPt[i],maxPt[i]); //
                histoRatioRawClusterPt[i]->Fit(pol1,"NRME0+","",minPt[i],maxPt[i]);
                pol1->SetLineColor(colorTrigg[i]-2);
                pol1->SetLineStyle(9);
                pol1->SetRange(minPt[i],maxPt[i]);
                pol1->Draw("same");

                TLegend* legendTriggRejectSingle = GetAndSetLegend2(0.33, 0.12, 0.92, 0.12+(1.05*(2)*textSizeSpectra2),textPixelPP);
                legendTriggRejectSingle->SetMargin(0.02);
                legendTriggRejectSingle->SetNColumns(3);
                legendTriggRejectSingle->AddEntry(histoRatioRawClusterPt[i],Form("   %s/%s",triggerNameLabel[i].Data(),triggerNameLabel[trigSteps[i][0]].Data() ),"p");

                legendTriggRejectSingle->AddEntry((TObject*)0,Form("%3.1f< #it{p}_{T} < %3.1f",minPt[i],maxPt[i]),"");
                if (triggerName[i].Contains("EMC1"))
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.0f #pm %3.0f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else if (triggerName[i].Contains("EMC7"))
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.1f #pm %3.1f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else if (triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA"))
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.2f #pm %3.2f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else if (triggerName[i].Contains("EG1") )
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.3f #pm %3.3f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.2f #pm %3.2f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");

                legendTriggRejectSingle->Draw();

                histo2DTriggRejectLinear->Draw("same,axis");

                canvasTriggerRejectLinear->Update();
                canvasTriggerRejectLinear->SaveAs(Form("%s/%s_TriggerRejectionFactorsLinY_Single_%s_%s.%s",outputDir.Data(), isMC.Data(), triggerName[i].Data(),
                                                    triggerName[trigSteps[i][0]].Data(), suffix.Data()));

                if (enableTriggerRejecCompMC) {
                    histo2DTriggRejectLinear->DrawCopy();

                    triggRejecCLPol0->Draw("e3,same");
                    pol0->Draw("same");

                    DrawGammaSetMarker(histoRatioRawClusterPt[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoRatioRawClusterPt[i]->DrawCopy("e1,same");
                    DrawGammaSetMarker(histoMCRatioRawClusterPt[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i]+2, colorTriggShade[i]+2);
                    histoMCRatioRawClusterPt[i]->DrawCopy("e1,same");

                    legendTriggRejectSingle->AddEntry(histoMCRatioRawClusterPt[i],"   scaled MC","p");
                    legendTriggRejectSingle->Draw();

                    histo2DTriggRejectLinear->Draw("same,axis");

                    canvasTriggerRejectLinear->Update();
                    canvasTriggerRejectLinear->SaveAs(Form("%s/%s_TriggerRejectionFactorsLinY_SingleWithMC_%s_%s.%s",outputDir.Data(), isMC.Data(), triggerName[i].Data(),
                                                        triggerName[trigSteps[i][0]].Data(), suffix.Data()));

                }

            }
        }

        //***************************************************************************************************************
        //********************Plotting trigger rejection factors (ENERGY!) = fits, linear scale 1 trigger at a time ***************
        //***************************************************************************************************************

        histo2DTriggRejectLinear->SetXTitle("#it{E} (GeV)");
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (i != trigSteps[i][0] ){
                histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,histoRatioRawClusterE[i]->GetMaximum()*1.5);
                if (triggerName[i].Contains("EMC7"))
                    histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,250);
                if (triggerName[i].Contains("EG2"))
                    histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,30);
                if (triggerName[i].Contains("EG1"))
                    histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,8);
                if (optionEnergy.CompareTo("8TeV")==0){
                    if (triggerName[i].Contains("EMC7"))
                        histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,150);
                    if (triggerName[i].Contains("EGA"))
                        histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,300);
                } else if (optionEnergy.CompareTo("13TeV")==0){
                    if (triggerName[i].Contains("EG2"))
                        histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,1000);
                    if (triggerName[i].Contains("EG1"))
                        histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,30);
                }

                histo2DTriggRejectLinear->DrawCopy();

                DrawGammaSetMarker(histoRatioRawClusterE[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoRatioRawClusterE[i]->DrawCopy("e1,same");
                TF1* pol0 = new TF1("pol0","[0]",minPt[i],maxPt[i]); //

                histoRatioRawClusterE[i]->Fit(pol0,"NRME+","",minPt[i],maxPt[i]);
                TH1D* triggRejecCLPol0 = (TH1D*)histoRatioRawClusterE[i]->Clone(Form("CL_%i",i));
                for (Int_t j = 1; j < triggRejecCLPol0->GetNbinsX()+1; j++){
                    triggRejecCLPol0->SetBinContent(j,triggRejecFac[i][trigSteps[i][0]]);
                    triggRejecCLPol0->SetBinError(j,triggRejecFacErr[i][trigSteps[i][0]]);
                }

                triggRejecCLPol0->SetStats(kFALSE);
                triggRejecCLPol0->SetFillColor(colorTriggShade[i]);
                triggRejecCLPol0->SetMarkerSize(0);
                triggRejecCLPol0->Draw("e3,same");

                pol0->SetParameter(0,triggRejecFac[i][trigSteps[i][0]]);
                pol0->SetParError(0,triggRejecFacErr[i][trigSteps[i][0]]);
                pol0->SetLineColor(colorTrigg[i]-1);
                pol0->SetLineStyle(7);
                pol0->SetRange(minPt[i],maxPt[i]);
                pol0->Draw("same");
                histoRatioRawClusterE[i]->DrawCopy("e1,same");

                TF1* pol1 = new TF1("pol1","[0]+[1]*x",minPt[i],maxPt[i]); //
                histoRatioRawClusterE[i]->Fit(pol1,"NRME0+","",minPt[i],maxPt[i]);
                pol1->SetLineColor(colorTrigg[i]-2);
                pol1->SetLineStyle(9);
                pol1->SetRange(minPt[i],maxPt[i]);
                pol1->Draw("same");

                TLegend* legendTriggRejectSingle = GetAndSetLegend2(0.33, 0.12, 0.92, 0.12+(1.05*(2)*textSizeSpectra2),textPixelPP);
                legendTriggRejectSingle->SetMargin(0.02);
                legendTriggRejectSingle->SetNColumns(3);
                legendTriggRejectSingle->AddEntry(histoRatioRawClusterE[i],Form("   %s/%s",triggerNameLabel[i].Data(),triggerNameLabel[trigSteps[i][0]].Data() ),"p");

                legendTriggRejectSingle->AddEntry((TObject*)0,Form("%3.1f< #it{p}_{T} < %3.1f",minPt[i],maxPt[i]),"");
                if (triggerName[i].Contains("EMC1"))
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.0f #pm %3.0f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else if (triggerName[i].Contains("EMC7"))
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.1f #pm %3.1f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else if (triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA"))
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.2f #pm %3.2f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else if (triggerName[i].Contains("EG1") )
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.3f #pm %3.3f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");
                else
                    legendTriggRejectSingle->AddEntry((TObject*)0,Form("     %3.2f #pm %3.2f",triggRejecFac[i][trigSteps[i][0]], triggRejecFacErr[i][trigSteps[i][0]]),"");

                legendTriggRejectSingle->Draw();

                histo2DTriggRejectLinear->Draw("same,axis");

                canvasTriggerRejectLinear->Update();
                canvasTriggerRejectLinear->SaveAs(Form("%s/%s_TriggerRejectionFactorsLinY_E_Single_%s_%s.%s",outputDir.Data(), isMC.Data(), triggerName[i].Data(),
                                                    triggerName[trigSteps[i][0]].Data(), suffix.Data()));

                if (enableTriggerRejecCompMC) {
                    histo2DTriggRejectLinear->DrawCopy();

                    triggRejecCLPol0->Draw("e3,same");
                    pol0->Draw("same");

                    DrawGammaSetMarker(histoRatioRawClusterPt[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoRatioRawClusterPt[i]->DrawCopy("e1,same");
                    DrawGammaSetMarker(histoMCRatioRawClusterPt[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i]+2, colorTriggShade[i]+2);
                    histoMCRatioRawClusterPt[i]->DrawCopy("e1,same");

                    legendTriggRejectSingle->AddEntry(histoMCRatioRawClusterPt[i],"   scaled MC","p");
                    legendTriggRejectSingle->Draw();

                    histo2DTriggRejectLinear->Draw("same,axis");

                    canvasTriggerRejectLinear->Update();
                    canvasTriggerRejectLinear->SaveAs(Form("%s/%s_TriggerRejectionFactorsLinY_E_SingleWithMC_%s_%s.%s",outputDir.Data(), isMC.Data(), triggerName[i].Data(),
                                                        triggerName[trigSteps[i][0]].Data(), suffix.Data()));

                }

            }
        }
        delete canvasTriggerRejectLinear;

        //***************************************************************************************************************
        //************************************Plotting trigger rejection trials *************************************************
        //***************************************************************************************************************
        TCanvas* canvasTriggerRejecTrial = new TCanvas("canvasTriggerRejecTrial","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasTriggerRejecTrial, 0.09, 0.017, textSizeSpectra, 0.08);
        canvasTriggerRejecTrial->SetLogy(0);

        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if ( triggerName[i].CompareTo(triggerName[trigSteps[i][0]].Data()) != 0){
                histoTriggerRejection[i]->SetTitle("");
                DrawGammaSetMarker(histoTriggerRejection[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoTriggerRejection[i]->Draw("e1,p");

                TF1* pol0_2 = new TF1("pol0_2","[0]",0.5,10.5); //
                histoTriggerRejection[i]->Fit(pol0_2,"QNRMEX0+","",0.5,10.5);
                pol0_2->SetLineColor(colorTrigg[i]-1);
                pol0_2->SetLineStyle(7);
                pol0_2->Draw("same");

                canvasTriggerRejecTrial->Update();
                canvasTriggerRejecTrial->SaveAs(Form("%s/TriggRejectMean_%s_%s_%s.%s",outputDir.Data(), isMC.Data(), triggerName[i].Data(), triggerName[trigSteps[i][0]].Data(),suffix.Data()));

                delete pol0_2;
            }
        }
        delete canvasTriggerRejecTrial;



        //***************************************************************************************************************
        //************************************Plotting unscaled invariant raw-yield EtaPrime *********************************
        //***************************************************************************************************************

        TCanvas* canvasClusterYield = new TCanvas("canvasClusterYield","",0,0,1000,1350);// gives the page size
        DrawGammaCanvasSettings( canvasClusterYield, 0.16, 0.02, 0.015, 0.07);
        canvasClusterYield->SetLogy();

        Double_t minClusYieldUnscaled    = 7e-9;
        Double_t maxClusYieldUnscaled    = 5;
        if (modeNormal == 4) {
            minClusYieldUnscaled         = 7e-10;
            maxClusYieldUnscaled         = 5;
        }
        if(optionEnergy.CompareTo("8TeV")==0){
        if (modeNormal == 2) {
            minClusYieldUnscaled         = 7e-10;
            maxClusYieldUnscaled         = 5;
        } else if (modeNormal == 4) {
            minClusYieldUnscaled         = 7e-10;
            maxClusYieldUnscaled         = 5;
        }
        }
        TH2F * histo2DClusUnscaled       = new TH2F("histo2DClusUnscaled", "histo2DClusUnscaled", 1000, 0., maxPtGlobalCluster, 10000, minClusYieldUnscaled, maxClusYieldUnscaled);
        SetStyleHistoTH2ForGraphs(histo2DClusUnscaled, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{#gamma, raw}}{#it{N}_{evt}d#it{p}_{T}} (#it{c}/GeV)^{2}",
                                0.85*textSizeSpectra,0.04, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.7);
        histo2DClusUnscaled->GetXaxis()->SetLabelOffset(-0.005);
        histo2DClusUnscaled->DrawCopy();

        Double_t minXLegendClus = 0.62;
        Int_t nColumnsClus      = 2;
        if (nrOfTrigToBeComb > 6){
            minXLegendClus = 0.3;
            nColumnsClus   = 3;
        }
        Double_t maxYLegendClus = 0.95;
        Double_t minYLegendClus = 0.95-(nrOfTrigToBeComb/nColumnsClus*0.75*textSizeSpectra);
        if(optionEnergy.CompareTo("8TeV")==0) minYLegendClus-=0.03;

        TLegend* legendClusUnscaled      = GetAndSetLegend2(minXLegendClus, minYLegendClus, 0.95, maxYLegendClus,0.85*textSizePixelSpectra);
        legendClusUnscaled->SetNColumns(nColumnsClus);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            DrawGammaSetMarker(histoRawClusterPt[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoRawClusterPt[i]->DrawCopy("e1,same");
            legendClusUnscaled->AddEntry(histoRawClusterPt[i],triggerNameLabel[i].Data(),"p");
        }
        legendClusUnscaled->Draw();


        TLatex *labelEnergyClusUnscaled  = new TLatex(0.2, 0.12+(2*textSizeSpectra*0.85*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyClusUnscaled, 0.85*textSizeSpectra,4);
        labelEnergyClusUnscaled->Draw();

        TLatex *labelClusterUnscaled     = new TLatex(0.2, 0.12+textSizeSpectra*0.85,"#gamma_{cand}");
        SetStyleTLatex( labelClusterUnscaled, 0.85*textSizeSpectra,4);
        labelClusterUnscaled->Draw();

        TLatex *labelDetProcClus = NULL;
        if (modeNormal == 4 || modeNormal == 2){
            labelDetProcClus = new TLatex(0.2, 0.12,detectionProcessClus.Data());
        } else {
            labelDetProcClus = new TLatex(0.2, 0.12,"rec. with EMCal");
        }

        SetStyleTLatex( labelDetProcClus, 0.85*textSizeSpectra,4);
        labelDetProcClus->Draw();

        canvasClusterYield->Update();
        canvasClusterYield->SaveAs(Form("%s/Cluster_%s_YieldUnscaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        if (histoRawClusterE[0]){
            TH2F * histo2DClusEUnscaled       = new TH2F("histo2DClusEUnscaled", "histo2DClusEUnscaled", 1000, 0., maxPtGlobalCluster, 10000, minClusYieldUnscaled, maxClusYieldUnscaled);
            SetStyleHistoTH2ForGraphs(histo2DClusEUnscaled, "#it{E} (GeV)","#frac{d#it{N}_{#gamma, raw}}{#it{N}_{evt}d#it{E}}} (1/GeV)^{2}",
                                    0.85*textSizeSpectra,0.04, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.7);
            histo2DClusEUnscaled->GetXaxis()->SetLabelOffset(-0.005);
            histo2DClusEUnscaled->DrawCopy();

            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                DrawGammaSetMarker(histoRawClusterE[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoRawClusterE[i]->DrawCopy("e1,same");
            }
            legendClusUnscaled->Draw();


            labelEnergyClusUnscaled->Draw();
            labelClusterUnscaled->Draw();
            labelDetProcClus->Draw();

            canvasClusterYield->Update();
            canvasClusterYield->SaveAs(Form("%s/Cluster_E_%s_YieldUnscaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        delete canvasClusterYield;
    }

    //***************************************************************************************************************
    //************************************Plotting efficiencies EtaPrime *************************************************
    //***************************************************************************************************************
    TCanvas* canvasEffi = new TCanvas("canvasEffi","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasEffi, 0.09, 0.017, 0.015, 0.08);
    canvasEffi->SetLogy(1);

    Double_t minEffiEtaPrime = 1e-4;
    Double_t maxEffiEtaPrime = 1e-1;
    if (modeNormal == 4){
        maxEffiEtaPrime      = 8e-1;
    } else if (modeNormal == 2){
        if(optionEnergy.CompareTo("8TeV")==0)
            minEffiEtaPrime  = 5e-5;
    } else if (modeNormal == 0){
        maxEffiEtaPrime      = 8e-3;
        minEffiEtaPrime      = 5e-5;
    } else if (modeNormal == 5){
      maxEffiEtaPrime        = 8e-1;
      minEffiEtaPrime        = 5e-2;
    }

    TH2F * histo2DEffiEtaPrime;
    histo2DEffiEtaPrime = new TH2F("histo2DEffiEtaPrime","histo2DEffiEtaPrime",1000,0., maxPtGlobalEtaPrime,10000,minEffiEtaPrime, maxEffiEtaPrime);
    SetStyleHistoTH2ForGraphs(histo2DEffiEtaPrime, "#it{p}_{T} (GeV/#it{c})","#it{#varepsilon}_{#eta'}#upoint#it{#kappa}_{trigg}",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.05);
    histo2DEffiEtaPrime->DrawCopy();
    histo2DEffiEtaPrime->SetYTitle("#it{#varepsilon}_{#eta'}");

    Double_t minXLegendEffi = 0.62;
    Int_t nColumnsEffi      = 2;
    if (nrOfTrigToBeComb > 6){
        minXLegendEffi = 0.4;
        nColumnsEffi   = 3;
    }
    Double_t minYLegendEffi = 0.13;
    Double_t maxYLegendEffi = minYLegendEffi+(1.05*nrOfTrigToBeComb/nColumnsEffi*0.85*textSizeSpectra);

    TLegend* legendEffiEtaPrime = GetAndSetLegend2(minXLegendEffi, minYLegendEffi, 0.95, maxYLegendEffi ,28);
    legendEffiEtaPrime->SetNColumns(nColumnsEffi);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoEfficiencyEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoEfficiencyEtaPrime[i]->DrawCopy("e1,same");
        legendEffiEtaPrime->AddEntry(histoEfficiencyEtaPrime[i],triggerNameLabel[i].Data(),"p");
    }
    legendEffiEtaPrime->Draw();

    TLatex *labelEnergyEffi = new TLatex(0.95, maxYLegendEffi+0.02+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
    SetStyleTLatex( labelEnergyEffi, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
    labelEnergyEffi->Draw();

    TLatex *labelEtaPrimeEffi = new TLatex(0.95, maxYLegendEffi+0.02+0.99*textSizeSpectra*0.85,"#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelEtaPrimeEffi, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
    labelEtaPrimeEffi->Draw();

    TLatex *labelDetProcEffi = new TLatex(0.95, maxYLegendEffi+0.02,detectionProcess.Data());
    SetStyleTLatex( labelDetProcEffi, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
    labelDetProcEffi->Draw();

    canvasEffi->Update();
    canvasEffi->SaveAs(Form("%s/EtaPrime_Efficiency.%s",outputDir.Data(),suffix.Data()));

    //***************************************************************************************************************
    //************************************ Plotting trigger efficiencies EtaPrime ****************************************
    //***************************************************************************************************************
    if (enableTriggerEffEtaPrimeAll){
        TCanvas* canvasTriggerEffi = new TCanvas("canvasTriggerEffi","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasTriggerEffi, 0.09, 0.017, 0.015, 0.08);
        canvasTriggerEffi->SetLogy(0);

        Double_t minEffiTrigEtaPrime         = 0;
        Double_t maxEffiTrigEtaPrime         = 1.1;
        //if(optionEnergy.CompareTo("8TeV") == 0 && modeNormal == 2) maxEffiTrigEtaPrime = 1.5;

        TH2F * histo2DTriggerEffiEtaPrime;
        histo2DTriggerEffiEtaPrime = new TH2F("histo2DTriggerEffiEtaPrime","histo2DTriggerEffiEtaPrime",1000,0., maxPtGlobalEtaPrime,10000,minEffiTrigEtaPrime, maxEffiTrigEtaPrime);
        SetStyleHistoTH2ForGraphs(histo2DTriggerEffiEtaPrime, "#it{p}_{T} (GeV/#it{c})","#it{#kappa}_{trigg, #eta'}",
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.05);
        histo2DTriggerEffiEtaPrime->DrawCopy();

        TLegend* legendTriggerEffiEtaPrime = GetAndSetLegend2(0.62, 0.165, 0.95, 0.165+(1.05*nRealTriggers/2*0.85*textSizeSpectra),28);
        legendTriggerEffiEtaPrime->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (enableTriggerEffEtaPrime[i]){
                DrawGammaSetMarker(histoTriggerEffEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoTriggerEffEtaPrime[i]->DrawCopy("e1,same");
                legendTriggerEffiEtaPrime->AddEntry(histoTriggerEffEtaPrime[i],triggerNameLabel[i].Data(),"p");
            }
        }
        legendTriggerEffiEtaPrime->Draw();

        labelEnergyEffi->Draw();
        labelEtaPrimeEffi->Draw();
        labelDetProcEffi->Draw();

        canvasTriggerEffi->Update();
        canvasTriggerEffi->SaveAs(Form("%s/EtaPrime_TriggerEfficiency.%s",outputDir.Data(),suffix.Data()));
        delete canvasTriggerEffi;

        //***************************************************************************************************************
        //******************* Plotting efficiencies EtaPrime without trigger efficiency folded in  ***************************
        //***************************************************************************************************************
        canvasEffi->cd();
        histo2DEffiEtaPrime->DrawCopy();

        TGraphErrors* graphEffBaseEtaPrime[12]   = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
        TLegend* legendEffiEtaPrimeW0TriggEff    = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
        legendEffiEtaPrimeW0TriggEff->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if ( triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1") ){
                DrawGammaSetMarker(histoEfficiencyEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoEfficiencyEtaPrime[i]->DrawCopy("e1,same");
                legendEffiEtaPrimeW0TriggEff->AddEntry(histoEfficiencyEtaPrime[i],triggerNameLabel[i].Data(),"p");
            } else {
                if (enableTriggerEffEtaPrime[i]){
                    graphEffBaseEtaPrime[i]      = new TGraphErrors(histoEffBaseEtaPrime[i]);
                    cout << "trying to plot trigg eff" << endl;
                    DrawGammaSetMarkerTGraph(graphEffBaseEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    graphEffBaseEtaPrime[i]->Draw("p,e1,same");
                    legendEffiEtaPrimeW0TriggEff->AddEntry(graphEffBaseEtaPrime[i],triggerNameLabel[i].Data(),"p");
                }
            }
        }
        legendEffiEtaPrimeW0TriggEff->Draw();

        labelEnergyEffi->Draw();
        labelEtaPrimeEffi->Draw();
        labelDetProcEffi->Draw();

        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/EtaPrime_EfficiencyW0TriggEff.%s",outputDir.Data(),suffix.Data()));
    }

    // delete canvasEffi;

    //***************************************************************************************************************
    //************************************Plotting acceptance EtaPrime *************************************************
    //***************************************************************************************************************
    TCanvas* canvasAcc = new TCanvas("canvasAcc","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasAcc, 0.1, 0.017, 0.015, 0.08);
    canvasAcc->SetLogy(0);

    Double_t minAccEtaPrime = 0.15;
    Double_t maxAccEtaPrime = 0.3;
    if (modeNormal == 0){
        maxAccEtaPrime       = 1.05;
        minAccEtaPrime       = 0.4;
    } else if (modeNormal == 2){
        minAccEtaPrime       = 0.15;
        maxAccEtaPrime       = 0.3;
    } else if (modeNormal == 3){
        minAccEtaPrime       = 0.;
        maxAccEtaPrime       = 0.1;
    } else if (modeNormal == 4){
        minAccEtaPrime       = 0.15;
        maxAccEtaPrime       = 0.3;
    } else if (modeNormal == 5){
        minAccEtaPrime       = 0.;
        maxAccEtaPrime       = 0.05;
    }


    TH2F * histo2DAccEtaPrime;
    histo2DAccEtaPrime = new TH2F("histo2DAccEtaPrime","histo2DAccEtaPrime",1000,0., maxPtGlobalEtaPrime,10000,minAccEtaPrime, maxAccEtaPrime);
    SetStyleHistoTH2ForGraphs(histo2DAccEtaPrime, "#it{p}_{T} (GeV/#it{c})","A_{#eta'}",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.25);
    histo2DAccEtaPrime->DrawCopy();

        Double_t minXLegendAcc = 0.62;
        Int_t nColumnsAcc      = 2;
        if (nrOfTrigToBeComb > 6){
            minXLegendAcc = 0.4;
            nColumnsAcc   = 3;
        }
        Double_t minYLegendAcc = 0.13;
        Double_t maxYLegendAcc = minYLegendAcc+(1.05*nrOfTrigToBeComb/nColumnsAcc*0.85*textSizeSpectra);

        TLegend* legendAccEtaPrime = GetAndSetLegend2(minXLegendAcc, minYLegendAcc, 0.95,maxYLegendAcc,28);
        legendAccEtaPrime->SetNColumns(nColumnsAcc);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            DrawGammaSetMarker(histoAcceptanceEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoAcceptanceEtaPrime[i]->DrawCopy("e1,same");
            legendAccEtaPrime->AddEntry(histoAcceptanceEtaPrime[i],triggerNameLabel[i].Data(),"p");
        }
        legendAccEtaPrime->Draw();

        labelEnergyEffi->Draw();
        labelEtaPrimeEffi->Draw();
        labelDetProcEffi->Draw();

    canvasAcc->Update();
    canvasAcc->SaveAs(Form("%s/EtaPrime_Acceptance.%s",outputDir.Data(),suffix.Data()));

    //***************************************************************************************************************
    //**************************** Mass and Width general plotting definitions **************************************
    //***************************************************************************************************************
    TCanvas* canvasMass         = new TCanvas("canvasMass","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasMass, 0.11, 0.017, 0.015, 0.08);

    Double_t minMassEtaPrime         = 0.8;
    Double_t maxMassEtaPrime         = 1.2;

    TH2F * histo2DMassEtaPrime       = new TH2F("histo2DMassEtaPrime","histo2DMassEtaPrime",1000,0., maxPtGlobalEtaPrime,10000,minMassEtaPrime, maxMassEtaPrime);
    SetStyleHistoTH2ForGraphs(histo2DMassEtaPrime, "#it{p}_{T} (GeV/#it{c})","#it{M}_{#eta'} (Mev/#it{c}^{2})",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.4);
    TLegend* legendMassEtaPrime      = GetAndSetLegend2(0.52, 0.88, 0.95, 0.88+(1.05*4/2*0.85*textSizeSpectra),28);
    TLegend* legendMassRedEtaPrime   = GetAndSetLegend2(0.52, 0.88, 0.95, 0.88+(1.05*4/2*0.85*textSizeSpectra),28);
    TLatex *labelEnergyMass     = new TLatex(0.14, 0.85+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
    TLatex *labelEtaPrimeMass        = new TLatex(0.14, 0.85+0.99*textSizeSpectra*0.85,"#eta' #rightarrow #gamma#gamma");
    TLatex *labelDetProcMass    = new TLatex(0.14, 0.85,detectionProcess.Data());
    TLegend* legendMassEtaPrime2     = GetAndSetLegend2(0.46, 0.81, 0.80, 0.81+(1.05*8/2*0.85*textSizeSpectra),28);
    TLegend* legendMassRedEtaPrime2  = GetAndSetLegend2(0.46, 0.81, 0.80, 0.81+(1.05*8/2*0.85*textSizeSpectra),28);

    TCanvas* canvasWidth        = new TCanvas("canvasWidth","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasWidth, 0.09, 0.017, 0.035, 0.08);
    Double_t minWidthEtaPrime        = 0.0;
    Double_t maxWidthEtaPrime        = 0.4;
    TH2F * histo2DWidthEtaPrime      = new TH2F("histo2DWidthEtaPrime","histo2DWidthEtaPrime",1000,0., maxPtGlobalEtaPrime,10000,minWidthEtaPrime, maxWidthEtaPrime);
    SetStyleHistoTH2ForGraphs(histo2DWidthEtaPrime, "#it{p}_{T} (GeV/#it{c})","#sigma_{#eta'} (Mev/#it{c}^{2})",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
    TLegend* legendWidthEtaPrime     = GetAndSetLegend2(0.52, 0.87, 0.95, 0.87+(1.05*4/2*0.85*textSizeSpectra),28);
    TLegend* legendWidthRedEtaPrime  = GetAndSetLegend2(0.52, 0.87, 0.95, 0.87+(1.05*4/2*0.85*textSizeSpectra),28);
    TLatex* labelEnergyWidth    = new TLatex(0.14, 0.84+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
    TLatex* labelEtaPrimeWidth       = new TLatex(0.14, 0.84+0.99*textSizeSpectra*0.85,"#eta' #rightarrow #gamma#gamma");
    TLatex* labelDetProcWidth   = new TLatex(0.14, 0.84,detectionProcess.Data());
    TLegend* legendWidthEtaPrime2    = GetAndSetLegend2(0.46, 0.80, 0.80, 0.80+(1.05*8/2*0.85*textSizeSpectra),28);
    TLegend* legendWidthRedEtaPrime2 = GetAndSetLegend2(0.46, 0.80, 0.80, 0.80+(1.05*8/2*0.85*textSizeSpectra),28);

    //***************************************************************************************************************
    //************************************Plotting Mass EtaPrime *********************************************************
    //***************************************************************************************************************
    canvasMass->cd();
    histo2DMassEtaPrime->DrawCopy();

    legendMassEtaPrime->SetNColumns(2);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if((optionEnergy.CompareTo("2.76TeV")==0 ) ||
            (optionEnergy.CompareTo("8TeV")==0) ||
            (optionEnergy.CompareTo("pPb_5.023TeV")==0) ||
            (optionEnergy.CompareTo("XeXe_5.44TeV")==0)
            ){
            DrawGammaSetMarker(histoMassEtaPrimeData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoMassEtaPrimeData[i]->DrawCopy("e1,same");
            legendMassEtaPrime->AddEntry(histoMassEtaPrimeData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
            DrawGammaSetMarker(histoMassEtaPrimeMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
            histoMassEtaPrimeMC[i]->DrawCopy("e1,same");
            legendMassEtaPrime->AddEntry(histoMassEtaPrimeMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
        }
    }
    legendMassEtaPrime->Draw();

    SetStyleTLatex( labelEnergyMass, 0.85*textSizeSpectra,4);
    labelEnergyMass->Draw();


    SetStyleTLatex( labelEtaPrimeMass, 0.85*textSizeSpectra,4);
    labelEtaPrimeMass->Draw();

    SetStyleTLatex( labelDetProcMass, 0.85*textSizeSpectra,4);
    labelDetProcMass->Draw();

    canvasMass->Update();
    canvasMass->SaveAs(Form("%s/EtaPrime_%s_Mass.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    histo2DMassEtaPrime->DrawCopy();

    //***************************************************************************************************************
    //************************************Plotting Width EtaPrime *********************************************************
    //***************************************************************************************************************
    canvasWidth->cd();
    histo2DWidthEtaPrime->DrawCopy();

    legendWidthEtaPrime->SetNColumns(2);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if((optionEnergy.CompareTo("2.76TeV")==0 ) ||
            (optionEnergy.CompareTo("8TeV")==0) ||
            (optionEnergy.CompareTo("pPb_5.023TeV")==0) ||
            (optionEnergy.CompareTo("XeXe_5.44TeV")==0)
            ){
            DrawGammaSetMarker(histoWidthEtaPrimeData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoWidthEtaPrimeData[i]->DrawCopy("e1,same");
            legendWidthEtaPrime->AddEntry(histoWidthEtaPrimeData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
            DrawGammaSetMarker(histoWidthEtaPrimeMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
            histoWidthEtaPrimeMC[i]->DrawCopy("e1,same");
            legendWidthEtaPrime->AddEntry(histoWidthEtaPrimeMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
        }
    }
    legendWidthEtaPrime->Draw();

    SetStyleTLatex( labelEnergyWidth, 0.85*textSizeSpectra,4);
    labelEnergyWidth->Draw();

    SetStyleTLatex( labelEtaPrimeWidth, 0.85*textSizeSpectra,4);
    labelEtaPrimeWidth->Draw();

    SetStyleTLatex( labelDetProcWidth, 0.85*textSizeSpectra,4);
    labelDetProcWidth->Draw();

    canvasWidth->Update();
    canvasWidth->SaveAs(Form("%s/EtaPrime_%s_Width.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    histo2DWidthEtaPrime->DrawCopy();

    //***************************************************************************************************************
    //************************************Plotting unscaled invariant raw-yield EtaPrime *********************************
    //***************************************************************************************************************
    textSizePixelSpectra = textSizeSpectra*1000;

    TCanvas* canvasRawUnscaled = new TCanvas("canvasRawUnscaled","",0,0,1000,1350);// gives the page size
    DrawGammaCanvasSettings( canvasRawUnscaled, 0.16, 0.02, 0.015, 0.07);
    canvasRawUnscaled->SetLogy();

    Double_t minCorrYieldRawUnscaled    = 7e-8;
    Double_t maxCorrYieldRawUnscaled    = 4e-2;
    if (modeNormal == 4) {
        minCorrYieldRawUnscaled         = 7e-8;
        maxCorrYieldRawUnscaled         = 1;
    }

    if(optionEnergy.CompareTo("8TeV")==0){
      if(modeNormal == 2){
        minCorrYieldRawUnscaled         = 2e-8;
        maxCorrYieldRawUnscaled         = 8e-3;
      } else if(modeNormal == 4){
        minCorrYieldRawUnscaled         = 1e-7;
        maxCorrYieldRawUnscaled         = 2e-1;
      } else if (modeNormal == 0){
        minCorrYieldRawUnscaled         = 2e-8;
        maxCorrYieldRawUnscaled         = 4e-3;
      }
    }

    if(optionEnergy.CompareTo("pPb_5.023TeV")==0){
      if(modeNormal == 2){
        minCorrYieldRawUnscaled         = 2e-8;
        maxCorrYieldRawUnscaled         = 8e-3;
      }
    }

    TH2F * histo2DRawUnscaled       = new TH2F("histo2DRawUnscaled", "histo2DRawUnscaled", 1000, 0., maxPtGlobalEtaPrime, 10000, minCorrYieldRawUnscaled, maxCorrYieldRawUnscaled);
    SetStyleHistoTH2ForGraphs(histo2DRawUnscaled, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{#pi_{0}, raw}}{#it{N}_{evt}d#it{p}_{T}} (#it{c}/GeV)^{2}",
                            0.85*textSizeSpectra,0.04, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.7);
    histo2DRawUnscaled->GetXaxis()->SetLabelOffset(-0.005);
    histo2DRawUnscaled->DrawCopy();

    Double_t minXLegendRaw = 0.62;
    Int_t nColumnsRaw      = 2;
    if (nrOfTrigToBeComb > 6){
        minXLegendRaw = 0.3;
        nColumnsRaw   = 3;
    }
    Double_t maxYLegendRaw = 0.95;
    Double_t minYLegendRaw = 0.95-(nrOfTrigToBeComb/nColumnsRaw*0.75*textSizeSpectra);
    if(optionEnergy.CompareTo("8TeV")==0 || optionEnergy.CompareTo("pPb_5.023TeV")==0 ) minYLegendRaw-=0.03;

    TLegend* legendRawUnscaled      = GetAndSetLegend2(minXLegendRaw, minYLegendRaw, 0.95, maxYLegendRaw,0.85*textSizePixelSpectra);
    legendRawUnscaled->SetNColumns(nColumnsRaw);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoRawYieldEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoRawYieldEtaPrime[i]->DrawCopy("e1,same");
        legendRawUnscaled->AddEntry(histoRawYieldEtaPrime[i],triggerNameLabel[i].Data(),"p");
    }
    legendRawUnscaled->Draw();


    TLatex *labelEnergyRawUnscaled  = new TLatex(0.2, 0.12+(2*textSizeSpectra*0.85*0.75),collisionSystem.Data());
    SetStyleTLatex( labelEnergyRawUnscaled, 0.85*textSizeSpectra,4);
    labelEnergyRawUnscaled->Draw();

    TLatex *labelEtaPrimeRawUnscaled     = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelEtaPrimeRawUnscaled, 0.85*textSizeSpectra,4);
    labelEtaPrimeRawUnscaled->Draw();

    TLatex *labelDetProcRawUnscaled = new TLatex(0.2, 0.12,detectionProcess.Data());
    SetStyleTLatex( labelDetProcRawUnscaled, 0.85*textSizeSpectra,4);
    labelDetProcRawUnscaled->Draw();

    canvasRawUnscaled->Update();
    canvasRawUnscaled->SaveAs(Form("%s/EtaPrime_%s_RawYieldUnscaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
    delete canvasRawUnscaled;

    //***************************************************************************************************************
    //************************************Plotting unscaled invariant yield EtaPrime *************************************
    //***************************************************************************************************************
    TCanvas* canvasCorrUnscaled = new TCanvas("canvasCorrUnscaled","",0,0,1000,1350);// gives the page size
    DrawGammaCanvasSettings( canvasCorrUnscaled, 0.15, 0.017, 0.015, 0.07);
    canvasCorrUnscaled->SetLogy();

    Double_t minCorrYieldUnscaled   = 2e-10;
    Double_t maxCorrYieldUnscaled   = 1e2;
    if (modeNormal == 0){
        minCorrYieldUnscaled        = 1e-8;
        maxCorrYieldUnscaled        = 1e2;
    }

    if(optionEnergy.CompareTo("8TeV")==0){
      if(modeNormal == 2){
        minCorrYieldUnscaled        = 1e-8;
        maxCorrYieldUnscaled        = 1;
      }else if(modeNormal == 4){
        minCorrYieldUnscaled        = 2e-8;
        maxCorrYieldUnscaled        = 0.2;
      }
    }

    if(optionEnergy.CompareTo("pPb_5.023TeV")==0){
      if(modeNormal == 2){
        minCorrYieldUnscaled        = 1e-8;
        maxCorrYieldUnscaled        = 1;
      }else if(modeNormal == 4){
        minCorrYieldUnscaled        = 2e-8;
        maxCorrYieldUnscaled        = 0.2;
      }
    }

    TH2F * histo2DInvYieldUnscaled = new TH2F("histo2DInvYieldUnscaled","histo2DInvYieldUnscaled",1000,0., maxPtGlobalEtaPrime,10000,minCorrYieldUnscaled,maxCorrYieldUnscaled);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldUnscaled, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                            0.85*textSizeSpectra,0.04, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.55);
    histo2DInvYieldUnscaled->DrawCopy();

    TLegend* legendUnscaled = GetAndSetLegend2(minXLegendRaw, minYLegendRaw, 0.95, maxYLegendRaw,0.85*textSizePixelSpectra);
    legendUnscaled->SetNColumns(nColumnsRaw);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoCorrectedYieldEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoCorrectedYieldEtaPrime[i]->DrawCopy("e1,same");
        legendUnscaled->AddEntry(histoCorrectedYieldEtaPrime[i],triggerNameLabel[i].Data(),"p");
    }
    legendUnscaled->Draw();


    TLatex *labelEnergyUnscaled = new TLatex(0.2, 0.12+(2*textSizeSpectra*0.85*0.75),collisionSystem.Data());
    SetStyleTLatex( labelEnergyUnscaled, 0.85*textSizeSpectra,4);
    labelEnergyUnscaled->Draw();

    TLatex *labelEtaPrimeUnscaled = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelEtaPrimeUnscaled, 0.85*textSizeSpectra,4);
    labelEtaPrimeUnscaled->Draw();

    TLatex *labelDetProcUnscaled = new TLatex(0.2, 0.12,detectionProcess.Data());
    SetStyleTLatex( labelDetProcUnscaled, 0.85*textSizeSpectra,4);
    labelDetProcUnscaled->Draw();

    canvasCorrUnscaled->Update();
    canvasCorrUnscaled->SaveAs(Form("%s/EtaPrime_%s_CorrectedYieldUnscaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
    delete canvasCorrUnscaled;

    //***************************************************************************************************************
    //******************************* Scaling corrected yield by trigger rejection factors **************************
    //***************************************************************************************************************

    TH1D*     histoCorrectedYieldEtaPrimeScaled                  [MaxNumberOfFiles];
    TH1D*     histoCorrectedYieldEtaPrimeScaledMasked            [MaxNumberOfFiles];
    histoCorrectedYieldEtaPrimeScaled[0]                         = (TH1D*)histoCorrectedYieldEtaPrime[0]->Clone(Form("CorrectedYieldEtaPrimeScaled_%s", triggerName[0].Data()));
    for (Int_t i = 1; i< nrOfTrigToBeComb; i++){
        fileFitsOutput << triggerName[i].Data() << endl;
        histoCorrectedYieldEtaPrimeScaled[i] = (TH1D*)histoCorrectedYieldEtaPrime[i]->Clone(Form("CorrectedYieldEtaPrimeScaled_%s", triggerName[i].Data()));
        histoCorrectedYieldEtaPrimeScaled[i]->Sumw2();
        histoCorrectedYieldEtaPrimeScaled[i]->Scale(1./triggRejecFac[i][trigSteps[i][0]]);
        fileFitsOutput << trigSteps[i][0] << "\t" << trigSteps[i][1] << "\t" << trigSteps[i][2] << endl;
        if (trigSteps[i][1]!= trigSteps[i][0]){
            fileFitsOutput << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << endl;
            histoCorrectedYieldEtaPrimeScaled[i]->Scale(1./triggRejecFac[trigSteps[i][0]][trigSteps[i][1]]);
        }
        if (trigSteps[i][2]!= trigSteps[i][1]){
            fileFitsOutput << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << "\t"<< triggRejecFac[trigSteps[i][1]][trigSteps[i][2]] << endl;
            histoCorrectedYieldEtaPrimeScaled[i]->Scale(1./triggRejecFac[trigSteps[i][1]][trigSteps[i][2]]);
        }
    }

    // initialize all vectors for sytstematics and general creation of graphs, we can have a maximu of 100 data points at the moment
    Double_t xValueFinalEtaPrime                                 [400];
    Double_t xErrorLowFinalEtaPrime                              [400];
    Double_t xErrorHighFinalEtaPrime                             [400];
    Double_t yValueFinalEtaPrime                                 [400];
    Double_t yErrorLowFinalEtaPrime                              [400];
    Double_t yErrorHighFinalEtaPrime                             [400];
    Int_t nPointFinalEtaPrime                                         = 0;

    Double_t yErrorSysLowFinalEtaPrime                           [400];
    Double_t yErrorSysHighFinalEtaPrime                          [400];

    Double_t ptSysRelEtaPrime                                    [MaxNumberOfFiles][400];
    Double_t yErrorSysLowRelEtaPrime                             [MaxNumberOfFiles][400];
    Double_t yErrorSysHighRelEtaPrime                            [MaxNumberOfFiles][400];
    Bool_t sysAvailEtaPrime                                      [MaxNumberOfFiles];

    Bool_t sysAvailSingleEtaPrime                                [MaxNumberOfFiles];
    Int_t numberBinsSysAvailSingleEtaPrime                       [MaxNumberOfFiles];

    // graphs for easier reduction of measurements to desired range and systematic errors
    TGraphAsymmErrors* graphsCorrectedYieldShrunkEtaPrime        [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldSysShrunkEtaPrime     [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldRemoved0EtaPrime      [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldSysRemoved0EtaPrime   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphMassEtaPrimeData                     [MaxNumberOfFiles];
    TGraphAsymmErrors* graphMassEtaPrimeMC                       [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedMassEtaPrimeData              [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedMassEtaPrimeMC                [MaxNumberOfFiles];
    TGraphAsymmErrors* graphWidthEtaPrimeData                    [MaxNumberOfFiles];
    TGraphAsymmErrors* graphWidthEtaPrimeMC                      [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedWidthEtaPrimeData             [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedWidthEtaPrimeMC               [MaxNumberOfFiles];
    TGraphAsymmErrors* graphAcceptanceEtaPrime                   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedAcceptanceEtaPrime            [MaxNumberOfFiles];
    TGraphAsymmErrors* graphEfficiencyEtaPrime                   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEfficiencyEtaPrime            [MaxNumberOfFiles];
    TGraphAsymmErrors* graphEffTimesAccEtaPrime                  [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEffTimesAccEtaPrime           [MaxNumberOfFiles];

    Int_t nRelSysErrEtaPrimeSources          = 0;
    TGraphAsymmErrors* graphRelSysErrEtaPrimeSource              [30][MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedRelSysErrEtaPrimeSource       [30][MaxNumberOfFiles];
    for (Int_t j = 0; j< 30; j++){
        for (Int_t k = 0; k< MaxNumberOfFiles; k++){
            graphRelSysErrEtaPrimeSource[j][k]           = NULL;
            graphOrderedRelSysErrEtaPrimeSource[j][k]    = NULL;
        }
    }
    // definition of predefined arrays for trigger correlation filling
    TH1D*               histoStatEtaPrime    [12];
    TGraphAsymmErrors*  graphSystEtaPrime    [12];
    TH1D*               histoRelStatEtaPrime [12];
    TGraphAsymmErrors*  graphRelSystEtaPrime [12];

    Int_t offSetsEtaPrime[12]        =   { 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0 };
    Int_t offSetsEtaPrimeSys[12]     =   { 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0 };

    if(optionEnergy.CompareTo("8TeV")==0){
      if(modeNormal == 2){
        offSetsEtaPrime[1] = 3; //INT7
        offSetsEtaPrime[3] = 0; //EMC7
        offSetsEtaPrime[4] = 3; //EGA
      }else if(modeNormal == 4){
        offSetsEtaPrime[1] = 0; //INT7
        offSetsEtaPrime[3] = 0; //EMC7
        offSetsEtaPrime[4] = 4; //EGA
      }
    }

    // set all graphs to NULL first
    for (Int_t j = 0; j<12; j++){
        histoStatEtaPrime[j]                 = NULL;
        graphSystEtaPrime[j]                 = NULL;
        histoRelStatEtaPrime[j]              = NULL;
        graphRelSystEtaPrime[j]              = NULL;
        graphOrderedMassEtaPrimeData[j]      = NULL;
        graphOrderedMassEtaPrimeMC[j]        = NULL;
        graphOrderedWidthEtaPrimeData[j]     = NULL;
        graphOrderedWidthEtaPrimeMC[j]       = NULL;
        graphOrderedAcceptanceEtaPrime[j]    = NULL;
        graphOrderedEfficiencyEtaPrime[j]    = NULL;
        graphOrderedEffTimesAccEtaPrime[j]   = NULL;
    }


    //****************************************************************************************************************
    //************* Processing of each individual trigger, reducing ranges & adding systematics **********************
    //****************************************************************************************************************
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        // read systematics, if fileName is set to "bla" no action has to be performed and systematics will be disabled for the rest of the analysis
        if (sysFileEtaPrime[i].CompareTo("bla") != 0){
            sysAvailEtaPrime[i]                = kTRUE;
            ifstream  fileSysErrEtaPrime;
            fileSysErrEtaPrime.open(sysFileEtaPrime[i].Data(),ios_base::in);
            cout << sysFileEtaPrime[i].Data() << endl;
            gSystem->Exec(Form("cp %s %s/SystematicErrorAveraged_EtaPrime_%s.txt", sysFileEtaPrime[i].Data(), outputDir.Data(),triggerName[i].Data()));


            Int_t iPtBin = 0;
            cout << "reading sys file summed" << endl;
            while(!fileSysErrEtaPrime.eof() && iPtBin < 400){
                Double_t garbage = 0;
                fileSysErrEtaPrime >>ptSysRelEtaPrime[i][iPtBin] >> yErrorSysLowRelEtaPrime[i][iPtBin] >> yErrorSysHighRelEtaPrime[i][iPtBin]>>    garbage >> garbage;
                cout << iPtBin << "\t"<< ptSysRelEtaPrime[i][iPtBin]<< "\t"  << yErrorSysLowRelEtaPrime[i][iPtBin] << "\t"  <<yErrorSysHighRelEtaPrime[i][iPtBin] << "\t"  << endl;;
                iPtBin++;
            }
            fileSysErrEtaPrime.close();
            // read in detailed systematics
            string sysFileEtaPrimeDet = sysFileEtaPrime[i].Data();

            if(!replace(sysFileEtaPrimeDet, "Averaged", "AveragedSingle")){
                cout << "WARNING: could not find detailed systematics file " << sysFileEtaPrimeDet << ", skipping... " << endl;
                sysAvailSingleEtaPrime[i] = kFALSE;
            }else{
                ifstream fileSysErrDetailedEtaPrime;
                fileSysErrDetailedEtaPrime.open(sysFileEtaPrimeDet,ios_base::in);
                if(fileSysErrDetailedEtaPrime.is_open()) {
                    sysAvailSingleEtaPrime[i] = kTRUE;
                    gSystem->Exec(Form("cp %s %s/SystematicErrorAveragedSingle_EtaPrime_%s.txt", ((TString)sysFileEtaPrimeDet).Data(), outputDir.Data(),triggerName[i].Data()));
                } else{
                    sysAvailSingleEtaPrime[i] = kFALSE;
                    cout << "No single errors were found" << endl;
                }

                if (sysAvailSingleEtaPrime[i]){
                    cout << sysFileEtaPrimeDet << endl;
                    iPtBin = 0;
                    string line;
                    Int_t iPtBinColumn = 0;
                    while (getline(fileSysErrDetailedEtaPrime, line) && iPtBin < 400) {
                        istringstream ss(line);
                        TString temp="";
                        iPtBinColumn = 0;
                        while(ss && iPtBinColumn < 400){
                            ss >> temp;
                            if( !(iPtBin==0 && temp.CompareTo("bin")==0) && !temp.IsNull()){
                                ptSysDetail[i][iPtBin].push_back(temp);
                                iPtBinColumn++;
                            }
                        }
                        if(iPtBin == 0){
                            ptSysDetail[i][iPtBin++].push_back("TotalError");
                            iPtBinColumn++;
                        }else iPtBin++;
                    }
                    numberBinsSysAvailSingleEtaPrime[i] = iPtBin;
                    fileSysErrDetailedEtaPrime.close();
                 }
             }
        } else {
            sysAvailEtaPrime[i]             = kFALSE;
            sysAvailSingleEtaPrime[i]       = kFALSE;
        }
        cout << sysAvailEtaPrime[i] << "\t" << sysAvailSingleEtaPrime[i] << endl;
        // continue;


        // print out input spectrum from statistical histogram
        cout << "step 0" << endl;
        for (Int_t j = 1; j< histoCorrectedYieldEtaPrimeScaled[i]->GetNbinsX()+1; j++ ){
            cout << histoCorrectedYieldEtaPrimeScaled[i]->GetBinCenter(j) << "\t" << histoCorrectedYieldEtaPrimeScaled[i]->GetBinContent(j) << endl;
        }
        //*******************************************************************
        //**** create graphs from original histograms for every quantity ****
        //*******************************************************************
        cout << "step 1" << endl;
        // create correct yield graphs
        graphsCorrectedYieldShrunkEtaPrime[i]        = new TGraphAsymmErrors(histoCorrectedYieldEtaPrimeScaled[i]);
        graphsCorrectedYieldRemoved0EtaPrime[i]      = new TGraphAsymmErrors(histoCorrectedYieldEtaPrimeScaled[i]);
        graphsCorrectedYieldSysShrunkEtaPrime[i]     = new TGraphAsymmErrors(histoCorrectedYieldEtaPrimeScaled[i]);
        graphsCorrectedYieldSysRemoved0EtaPrime[i]   = new TGraphAsymmErrors(histoCorrectedYieldEtaPrimeScaled[i]);
        histoCorrectedYieldEtaPrimeScaledMasked[i]   = (TH1D*)histoCorrectedYieldEtaPrimeScaled[i]->Clone(Form("EtaPrime_ScaledMasked_%s",triggerName[i].Data()));

        // create supporting figure graphs
        graphMassEtaPrimeData[i]                 = new TGraphAsymmErrors(histoMassEtaPrimeData[i]);
        graphMassEtaPrimeMC[i]                   = new TGraphAsymmErrors(histoMassEtaPrimeMC[i]);
        graphWidthEtaPrimeData[i]                = new TGraphAsymmErrors(histoWidthEtaPrimeData[i]);
        graphWidthEtaPrimeMC[i]                  = new TGraphAsymmErrors(histoWidthEtaPrimeMC[i]);
        graphAcceptanceEtaPrime[i]                   = new TGraphAsymmErrors(histoAcceptanceEtaPrime[i]);

        if ( triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1")){
            graphEfficiencyEtaPrime[i]               = new TGraphAsymmErrors(histoEfficiencyEtaPrime[i]);
        } else { // only possible if same file had been processed with pure MB cut on the same MC
            if (enableTriggerEffEtaPrime[i]){ // only possible if same file had been processed with pure MB cut on the same MC
                graphEfficiencyEtaPrime[i]           = new TGraphAsymmErrors(histoEffBaseEtaPrime[i]);
            } else {
                graphEfficiencyEtaPrime[i]           = NULL;
            }
        }
        graphEffTimesAccEtaPrime[i]                  = new TGraphAsymmErrors(histoEffTimesAccEtaPrime[i]);

        // remove 0 bins at beginning
        Int_t binsToMask = 1;
        while (histoCorrectedYieldEtaPrimeScaledMasked[i]->GetBinCenter(binsToMask) < ptFromSpecEtaPrime[i][0] ){
            histoCorrectedYieldEtaPrimeScaledMasked[i]->SetBinContent(binsToMask,0.);
            histoCorrectedYieldEtaPrimeScaledMasked[i]->SetBinError(binsToMask,0.);
            binsToMask++;
        }
        while (graphAcceptanceEtaPrime[i]->GetX()[0] < ptFromSpecEtaPrime[i][0] ){
            graphMassEtaPrimeData[i]->RemovePoint(0);
            graphMassEtaPrimeMC[i]->RemovePoint(0);
            graphWidthEtaPrimeData[i]->RemovePoint(0);
            graphWidthEtaPrimeMC[i]->RemovePoint(0);
            graphAcceptanceEtaPrime[i]->RemovePoint(0);
            graphEffTimesAccEtaPrime[i]->RemovePoint(0);
            if (enableTriggerEffEtaPrime[i] || (triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1"))){
                graphEfficiencyEtaPrime[i]->RemovePoint(0);
            }
        }

        // check if trigger is supposed to be used for combination, otherwise put all graphs to NULL
        if ( maskedFullyEtaPrime[i] ){
            graphMassEtaPrimeData[i]     = NULL;
            graphMassEtaPrimeMC[i]       = NULL;
            graphWidthEtaPrimeData[i]    = NULL;
            graphWidthEtaPrimeMC[i]      = NULL;
            graphAcceptanceEtaPrime[i]   = NULL;
            graphEffTimesAccEtaPrime[i]  = NULL;
            graphEfficiencyEtaPrime[i]   = NULL;
            graphsCorrectedYieldShrunkEtaPrime[i]    = NULL;
            graphsCorrectedYieldRemoved0EtaPrime[i]  = NULL;
            graphsCorrectedYieldSysShrunkEtaPrime[i] = NULL;
            graphsCorrectedYieldSysRemoved0EtaPrime[i]   = NULL;
            sysAvailEtaPrime[i]          = kFALSE;
            for (Int_t f = 1; f < histoCorrectedYieldEtaPrimeScaledMasked[i]->GetNbinsX()+1; f++ ){
                histoCorrectedYieldEtaPrimeScaledMasked[i]->SetBinContent(f,0.);
                histoCorrectedYieldEtaPrimeScaledMasked[i]->SetBinError(f,0.);
            }
            nrOfTrigToBeCombEtaPrimeRed--;
            cout << "trigger " << triggerName[i] << " was masked" << endl;
            continue;
        // if upper boundary > -1, remove all points above
        } else if (ptFromSpecEtaPrime[i][1] > -1) {
            for (Int_t f = histoCorrectedYieldEtaPrimeScaledMasked[i]->GetXaxis()->FindBin(ptFromSpecEtaPrime[i][1]); f < histoCorrectedYieldEtaPrimeScaledMasked[i]->GetNbinsX()+1; f++ ){
                histoCorrectedYieldEtaPrimeScaledMasked[i]->SetBinContent(f,0.);
                histoCorrectedYieldEtaPrimeScaledMasked[i]->SetBinError(f,0.);
            }
            while (graphAcceptanceEtaPrime[i]->GetX()[graphAcceptanceEtaPrime[i]->GetN()-1] > ptFromSpecEtaPrime[i][1] ){
                graphMassEtaPrimeData[i]->RemovePoint(graphMassEtaPrimeData[i]->GetN()-1);
                graphMassEtaPrimeMC[i]->RemovePoint(graphMassEtaPrimeMC[i]->GetN()-1);
                graphWidthEtaPrimeData[i]->RemovePoint(graphWidthEtaPrimeData[i]->GetN()-1);
                graphWidthEtaPrimeMC[i]->RemovePoint(graphWidthEtaPrimeMC[i]->GetN()-1);
                graphAcceptanceEtaPrime[i]->RemovePoint(graphAcceptanceEtaPrime[i]->GetN()-1);
                graphEffTimesAccEtaPrime[i]->RemovePoint(graphEffTimesAccEtaPrime[i]->GetN()-1);
                if (enableTriggerEffEtaPrime[i] || (triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1"))){
                    graphEfficiencyEtaPrime[i]->RemovePoint(graphEfficiencyEtaPrime[i]->GetN()-1);
                }
            }
        }
        // if no points are left in graph put graph to NULL
        if (graphAcceptanceEtaPrime[i]->GetN() == 0){
            graphMassEtaPrimeData[i]     = NULL;
            graphMassEtaPrimeMC[i]       = NULL;
            graphWidthEtaPrimeData[i]    = NULL;
            graphWidthEtaPrimeMC[i]      = NULL;
            graphAcceptanceEtaPrime[i]       = NULL;
            graphEffTimesAccEtaPrime[i]      = NULL;
            graphEfficiencyEtaPrime[i]       = NULL;
        }

        // Remove 0 points at beginning for graphs
        // if (graphsCorrectedYieldShrunkEtaPrime[i])graphsCorrectedYieldShrunkEtaPrime[i]->Print();
        cout << "step 2" << endl;
        while (graphsCorrectedYieldShrunkEtaPrime[i]->GetY()[0] == 0) graphsCorrectedYieldShrunkEtaPrime[i]->RemovePoint(0);
        if (graphsCorrectedYieldShrunkEtaPrime[i]) graphsCorrectedYieldShrunkEtaPrime[i]->Print();
        while (graphsCorrectedYieldRemoved0EtaPrime[i]->GetY()[0] == 0) graphsCorrectedYieldRemoved0EtaPrime[i]->RemovePoint(0);
        if (graphsCorrectedYieldRemoved0EtaPrime[i]) graphsCorrectedYieldRemoved0EtaPrime[i]->Print();
        cout << "sys shrunk" << endl;
        while (graphsCorrectedYieldSysShrunkEtaPrime[i]->GetY()[0] == 0) graphsCorrectedYieldSysShrunkEtaPrime[i]->RemovePoint(0);
        if (graphsCorrectedYieldSysShrunkEtaPrime[i])graphsCorrectedYieldSysShrunkEtaPrime[i]->Print();
        cout << "sys shrunk 2" << endl;
        while (graphsCorrectedYieldSysRemoved0EtaPrime[i]->GetY()[0] == 0) graphsCorrectedYieldSysRemoved0EtaPrime[i]->RemovePoint(0);
        if (graphsCorrectedYieldSysRemoved0EtaPrime[i])graphsCorrectedYieldSysRemoved0EtaPrime[i]->Print();

        // put systematics on graphs
        if (graphsCorrectedYieldSysRemoved0EtaPrime[i]){
            if (sysAvailSingleEtaPrime[i]){
                nRelSysErrEtaPrimeSources                    = (Int_t)ptSysDetail[i][0].size()-1;
                for (Int_t k = 0; k < nRelSysErrEtaPrimeSources; k++ ){
                    graphRelSysErrEtaPrimeSource[k][i]       = (TGraphAsymmErrors*) graphsCorrectedYieldSysRemoved0EtaPrime[i]->Clone(Form("RelSysErrEtaPrimeSource%s_%s",((TString)ptSysDetail[i][0].at(k+1)).Data(), triggerName[i].Data()));
                    cout << Form("RelSysErrSource%s_%s",((TString)ptSysDetail[i][0].at(k+1)).Data(), triggerName[i].Data()) << endl;
                }
            }
            for (Int_t j = 0; j< graphsCorrectedYieldSysRemoved0EtaPrime[i]->GetN(); j++){
                if (sysAvailEtaPrime[i]){
                    Int_t counter = 0;
                    while(counter < 400 && TMath::Abs(graphsCorrectedYieldSysRemoved0EtaPrime[i]->GetX()[j] - ptSysRelEtaPrime[i][counter])> 0.001) counter++;
                    if (counter < 400){
                        cout << ptSysRelEtaPrime[i][counter]<< "\t found it" << endl;
                        Double_t yErrorSysLowDummy  = TMath::Abs(yErrorSysLowRelEtaPrime[i][counter]/100*graphsCorrectedYieldSysRemoved0EtaPrime[i]->GetY()[j]);
                        Double_t yErrorSysHighDummy = yErrorSysHighRelEtaPrime[i][counter]/100*graphsCorrectedYieldSysRemoved0EtaPrime[i]->GetY()[j];
                        graphsCorrectedYieldSysRemoved0EtaPrime[i]->SetPointEYlow(j,yErrorSysLowDummy);
                        graphsCorrectedYieldSysRemoved0EtaPrime[i]->SetPointEYhigh(j,yErrorSysHighDummy);
                        if (sysAvailSingleEtaPrime[i]){
                            for (Int_t k = 0; k < nRelSysErrEtaPrimeSources; k++ ){
                                graphRelSysErrEtaPrimeSource[k][i]->SetPoint(j, graphsCorrectedYieldSysRemoved0EtaPrime[i]->GetX()[j] ,((TString)ptSysDetail[i][counter+1].at(k+1)).Atof());
                                graphRelSysErrEtaPrimeSource[k][i]->SetPointEYhigh(j,0);
                                graphRelSysErrEtaPrimeSource[k][i]->SetPointEYlow(j,0);
                            }
                        }
                    } else {
                        graphsCorrectedYieldSysRemoved0EtaPrime[i]->SetPointEYlow(j,0);
                        graphsCorrectedYieldSysRemoved0EtaPrime[i]->SetPointEYhigh(j,0);
                        if (sysAvailSingleEtaPrime[i]){
                            for (Int_t k = 0; k < nRelSysErrEtaPrimeSources; k++ ){
                                graphRelSysErrEtaPrimeSource[k][i]->SetPoint(j, graphsCorrectedYieldSysRemoved0EtaPrime[i]->GetX()[j] ,0);
                                graphRelSysErrEtaPrimeSource[k][i]->SetPointEYlow(j,0);
                                graphRelSysErrEtaPrimeSource[k][i]->SetPointEYhigh(j,0);
                            }
                        }
                    }
                } else {
                    graphsCorrectedYieldSysRemoved0EtaPrime[i]->SetPointEYlow(j,0);
                    graphsCorrectedYieldSysRemoved0EtaPrime[i]->SetPointEYhigh(j,0);
                    averagedEtaPrime = kFALSE;
                    if (sysAvailSingleEtaPrime[i]){
                        for (Int_t k = 0; k < nRelSysErrEtaPrimeSources; k++ ){
                            graphRelSysErrEtaPrimeSource[k][i]->SetPoint(j, graphsCorrectedYieldSysRemoved0EtaPrime[i]->GetX()[j] ,0);
                            graphRelSysErrEtaPrimeSource[k][i]->SetPointEYlow(j,0);
                            graphRelSysErrEtaPrimeSource[k][i]->SetPointEYhigh(j,0);
                        }
                    }
                }
            }
        }

        cout << "step 3" << endl;
        cout << "range to be accepted: " << ptFromSpecEtaPrime[i][0] << "\t-\t" << ptFromSpecEtaPrime[i][1] << endl;
        // remove points at beginning according to ranges set for individual triggers
        while(graphsCorrectedYieldShrunkEtaPrime[i]->GetX()[0] < ptFromSpecEtaPrime[i][0])
            graphsCorrectedYieldShrunkEtaPrime[i]->RemovePoint(0);
        while(graphsCorrectedYieldSysShrunkEtaPrime[i]->GetX()[0] < ptFromSpecEtaPrime[i][0])
            graphsCorrectedYieldSysShrunkEtaPrime[i]->RemovePoint(0);
        // remove points at the end according to ranges set for individual triggers
        while (graphsCorrectedYieldShrunkEtaPrime[i]->GetX()[graphsCorrectedYieldShrunkEtaPrime[i]->GetN()-1] > ptFromSpecEtaPrime[i][1])
            graphsCorrectedYieldShrunkEtaPrime[i]->RemovePoint(graphsCorrectedYieldShrunkEtaPrime[i]->GetN()-1);
        graphsCorrectedYieldShrunkEtaPrime[i]->Print();
        while (graphsCorrectedYieldSysShrunkEtaPrime[i]->GetX()[graphsCorrectedYieldSysShrunkEtaPrime[i]->GetN()-1] > ptFromSpecEtaPrime[i][1])
            graphsCorrectedYieldSysShrunkEtaPrime[i]->RemovePoint(graphsCorrectedYieldSysShrunkEtaPrime[i]->GetN()-1);


        // put systematics on shrunk graphs
        for (Int_t j = 0; j< graphsCorrectedYieldShrunkEtaPrime[i]->GetN(); j++){
            xValueFinalEtaPrime[nPointFinalEtaPrime]      = graphsCorrectedYieldShrunkEtaPrime[i]->GetX()[j];
            xErrorHighFinalEtaPrime[nPointFinalEtaPrime]  = graphsCorrectedYieldShrunkEtaPrime[i]->GetEXhigh()[j];
            xErrorLowFinalEtaPrime[nPointFinalEtaPrime]   = graphsCorrectedYieldShrunkEtaPrime[i]->GetEXlow()[j];
            yValueFinalEtaPrime[nPointFinalEtaPrime]      = graphsCorrectedYieldShrunkEtaPrime[i]->GetY()[j];
            yErrorHighFinalEtaPrime[nPointFinalEtaPrime]  = graphsCorrectedYieldShrunkEtaPrime[i]->GetEYhigh()[j];
            yErrorLowFinalEtaPrime[nPointFinalEtaPrime]   = graphsCorrectedYieldShrunkEtaPrime[i]->GetEYlow()[j];
            if (sysAvailEtaPrime[i]){
                Int_t counter = 0;
                while(counter < 400 && TMath::Abs(xValueFinalEtaPrime[nPointFinalEtaPrime] - ptSysRelEtaPrime[i][counter])> 0.001) counter++;
                if (counter < 400){
                    cout << ptSysRelEtaPrime[i][counter]<< "\t found it" << endl;
                    yErrorSysLowFinalEtaPrime[nPointFinalEtaPrime] = TMath::Abs(yErrorSysLowRelEtaPrime[i][counter]/100*graphsCorrectedYieldShrunkEtaPrime[i]->GetY()[j]);
                    yErrorSysHighFinalEtaPrime[nPointFinalEtaPrime] = yErrorSysHighRelEtaPrime[i][counter]/100*graphsCorrectedYieldShrunkEtaPrime[i]->GetY()[j];

                } else {
                    yErrorSysLowFinalEtaPrime[nPointFinalEtaPrime] = 0;
                    yErrorSysHighFinalEtaPrime[nPointFinalEtaPrime] = 0;

                }
            } else {
                yErrorSysLowFinalEtaPrime[nPointFinalEtaPrime] = 0;
                yErrorSysHighFinalEtaPrime[nPointFinalEtaPrime] = 0;
            }
            graphsCorrectedYieldSysShrunkEtaPrime[i]->SetPointEYlow(j,yErrorSysLowFinalEtaPrime[nPointFinalEtaPrime]);
            graphsCorrectedYieldSysShrunkEtaPrime[i]->SetPointEYhigh(j,yErrorSysHighFinalEtaPrime[nPointFinalEtaPrime]);
            nPointFinalEtaPrime++;
        }

        // Set correct trigger order for combination function
        Int_t nCorrOrder    = GetOrderedTrigger(triggerName[i]);
        if (nCorrOrder == -1){
            cout << "ERROR: trigger name not defined" << endl;
            return;
        }

        if ( graphsCorrectedYieldShrunkEtaPrime[i]){
            histoStatEtaPrime[nCorrOrder]    = histoCorrectedYieldEtaPrimeScaledMasked[i];
            graphSystEtaPrime[nCorrOrder]    = graphsCorrectedYieldSysShrunkEtaPrime[i];
            offSetsEtaPrimeSys[nCorrOrder]   = histoStatEtaPrime[nCorrOrder]->GetXaxis()->FindBin(graphSystEtaPrime[nCorrOrder]->GetX()[0])-1;
            if (graphMassEtaPrimeData[i])
                graphOrderedMassEtaPrimeData[nCorrOrder]         = graphMassEtaPrimeData[i];
            if (graphMassEtaPrimeMC[i])
                graphOrderedMassEtaPrimeMC[nCorrOrder]           = graphMassEtaPrimeMC[i];
            if (graphWidthEtaPrimeData[i])
                graphOrderedWidthEtaPrimeData[nCorrOrder]        = graphWidthEtaPrimeData[i];
            if (graphWidthEtaPrimeMC[i])
                graphOrderedWidthEtaPrimeMC[nCorrOrder]          = graphWidthEtaPrimeMC[i];
            if (graphAcceptanceEtaPrime[i])
                graphOrderedAcceptanceEtaPrime[nCorrOrder]       = graphAcceptanceEtaPrime[i];
            if (graphEfficiencyEtaPrime[i])
                graphOrderedEfficiencyEtaPrime[nCorrOrder]       = graphEfficiencyEtaPrime[i];
            if (graphEffTimesAccEtaPrime[i])
                graphOrderedEffTimesAccEtaPrime[nCorrOrder]      = graphEffTimesAccEtaPrime[i];

            if (sysAvailSingleEtaPrime[i]){
                for (Int_t k = 0; k < nRelSysErrEtaPrimeSources; k++ ){
                    if (graphRelSysErrEtaPrimeSource[k][i])
                        graphOrderedRelSysErrEtaPrimeSource[k][nCorrOrder]   = graphRelSysErrEtaPrimeSource[k][i];
                }
            }
        }
        if (triggerName[i].Contains("INT7") && optionEnergy.CompareTo("8TeV")==0 && modeNormal == 2)
            offSetsEtaPrimeSys[1]+=3;
        if ((triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA")) && optionEnergy.CompareTo("8TeV")==0 && modeNormal == 4 )
            offSetsEtaPrimeSys[4]+=4;
        if ((triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA")) && optionEnergy.CompareTo("8TeV")==0 && modeNormal == 2 )
            offSetsEtaPrimeSys[4]+=3;
    }


    // create weighted graphs for spectra and supporting graphs
    TString nameWeightsLogFileEtaPrime =     Form("%s/weightsEtaPrime_%s.dat",outputDir.Data(),isMC.Data());
    TGraphAsymmErrors* graphCorrectedYieldWeightedAverageEtaPrimeStat    = NULL;
    TGraphAsymmErrors* graphCorrectedYieldWeightedAverageEtaPrimeSys     = NULL;
    TGraphAsymmErrors* graphCorrectedYieldWeightedAverageEtaPrimeTot     = NULL;
    TGraphAsymmErrors* graphMassEtaPrimeDataWeighted                     = NULL;
    TGraphAsymmErrors* graphMassEtaPrimeMCWeighted                       = NULL;
    TGraphAsymmErrors* graphWidthEtaPrimeDataWeighted                    = NULL;
    TGraphAsymmErrors* graphWidthEtaPrimeMCWeighted                      = NULL;
    TGraphAsymmErrors* graphAcceptanceEtaPrimeWeighted                   = NULL;
    TGraphAsymmErrors* graphEfficiencyEtaPrimeWeighted                   = NULL;
    TGraphAsymmErrors* graphEffTimesAccEtaPrimeWeighted                  = NULL;
    TGraphAsymmErrors* graphRelSysErrEtaPrimeSourceWeighted[30];

    for (Int_t k = 0; k < 30; k++){
        graphRelSysErrEtaPrimeSourceWeighted[k]                          = NULL;
    }
    // Calculate averaged etaprime spectrum & respective supporting graphs according to statistical and systematic errors taking correctly into account the cross correlations
    if (averagedEtaPrime){
        cout << maxNAllowedEtaPrime << endl;
        // Calculate average etaprime spectrum
        graphCorrectedYieldWeightedAverageEtaPrimeTot        = CombinePtPointsSpectraTriggerCorrMat(    histoStatEtaPrime, graphSystEtaPrime,
                                                                                                   binningEtaPrime,  maxNAllowedEtaPrime,
                                                                                                   offSetsEtaPrime ,offSetsEtaPrimeSys,
                                                                                                   graphCorrectedYieldWeightedAverageEtaPrimeStat, graphCorrectedYieldWeightedAverageEtaPrimeSys,
                                                                                                   nameWeightsLogFileEtaPrime.Data(),
                                                                                                   modeNormal, optionEnergy, "EtaPrime", "",
                                                                                                   fileInputCorrFactors
                                                                                               );

        // preparations for weight readout
        Double_t xValuesReadEtaPrime[400];
        Double_t weightsReadEtaPrime[12][400];
        Int_t availableMeasEtaPrime[12]       = {-1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1};
        Int_t nMeasSetEtaPrime               = nrOfTrigToBeCombEtaPrimeRed;
        Int_t nPtBinsReadEtaPrime            = 0;

        // labeling and plotting settings
        TString nameTriggerWeighted[12]  = {"INT1", "INT7", "EMC1", "EMC7", strEG2_A.Data(), "EG1",
                                            "INT1_NLM1", "INT7_NLM1", "EMC1_NLM1", "EMC7_NLM1", "EG2_NLM1", "EG1_NLM1"};
        Color_t colorTriggWeighted[12]   = {GetDefaultTriggerColorName(nameTriggerWeighted[0], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[1], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[2], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[3], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[4], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[5], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[6], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[7], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[8], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[9], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[10], 0),
                                           GetDefaultTriggerColorName(nameTriggerWeighted[11], 0)
                                          };
        Marker_t markerTriggWeighted[12] = {GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[0], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[1], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[2], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[3], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[4], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[5], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[6], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[7], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[8], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[9], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[10], 0),
                                           GetDefaultTriggerMarkerStyleName(nameTriggerWeighted[11], 0)
                                          };

        // Reading weights from output file for plotting
        ifstream fileWeightsEtaPrime;
        fileWeightsEtaPrime.open(nameWeightsLogFileEtaPrime,ios_base::in);
        cout << "reading" << nameWeightsLogFileEtaPrime << endl;

        while(!fileWeightsEtaPrime.eof() && nPtBinsReadEtaPrime < 400){
            TString garbage = "";
            if (nPtBinsReadEtaPrime == 0){
                fileWeightsEtaPrime >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
                for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                    fileWeightsEtaPrime >> availableMeasEtaPrime[i] ;
                }
                cout << "read following measurements: ";
                for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                    cout << availableMeasEtaPrime[i] << "\t" ;
                }
                cout << endl;
            } else {
                fileWeightsEtaPrime >> xValuesReadEtaPrime[nPtBinsReadEtaPrime-1];
                for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                    fileWeightsEtaPrime >> weightsReadEtaPrime[availableMeasEtaPrime[i]][nPtBinsReadEtaPrime-1] ;
                }
                cout << "read: "<<  nPtBinsReadEtaPrime << "\t"<< xValuesReadEtaPrime[nPtBinsReadEtaPrime-1] << "\t" ;
                for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                    cout << weightsReadEtaPrime[availableMeasEtaPrime[i]][nPtBinsReadEtaPrime-1] << "\t";
                }
                cout << endl;
            }
            nPtBinsReadEtaPrime++;
        }
        nPtBinsReadEtaPrime = nPtBinsReadEtaPrime-2 ;
        fileWeightsEtaPrime.close();

        // creating & filling the weight graphs
        TGraph* graphWeightsEtaPrime[12];
        for (Int_t i = 0; i < 12; i++){
            graphWeightsEtaPrime[i]                    = NULL;
        }
        for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
            cout << i << "\t" << availableMeasEtaPrime[i] << endl;
            graphWeightsEtaPrime[availableMeasEtaPrime[i]]  = new TGraph(nPtBinsReadEtaPrime,xValuesReadEtaPrime,weightsReadEtaPrime[availableMeasEtaPrime[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsReadEtaPrime; n++){
                if (graphWeightsEtaPrime[availableMeasEtaPrime[i]]->GetY()[bin] == 0) graphWeightsEtaPrime[availableMeasEtaPrime[i]]->RemovePoint(bin);
                else bin++;
            }
            graphWeightsEtaPrime[availableMeasEtaPrime[i]]->Print();
        }


        //  **********************************************************************************************************************
        //  ******************************************* Plotting weights EtaPrime *****************************************************
        //  **********************************************************************************************************************
        Int_t textSizeLabelsPixel = 900*0.04;

        TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);// gives the page size
        DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);

        TH2F * histo2DWeights;
        histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,0.,maxPtGlobalEtaPrime,1000,-0.5,1.1);
        SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DWeights->Draw("copy");

            TLegend* legendWeightsEtaPrime = GetAndSetLegend2(0.12, 0.14, 0.55, 0.14+(0.035*nMeasSetEtaPrime/2*1.35), 32);
            legendWeightsEtaPrime->SetNColumns(2);
            for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                DrawGammaSetMarkerTGraph(graphWeightsEtaPrime[availableMeasEtaPrime[i]],markerTriggWeighted[availableMeasEtaPrime[i]], sizeTrigg[availableMeasEtaPrime[i]],
                                         colorTriggWeighted[availableMeasEtaPrime[i]], colorTriggWeighted[availableMeasEtaPrime[i]]);
                graphWeightsEtaPrime[availableMeasEtaPrime[i]]->Draw("p,same,e1");
                legendWeightsEtaPrime->AddEntry(graphWeightsEtaPrime[availableMeasEtaPrime[i]],nameTriggerWeighted[availableMeasEtaPrime[i]],"p");
            }
            legendWeightsEtaPrime->Draw();

            TLatex *labelWeightsEnergy = new TLatex(0.95,0.24,collisionSystem.Data());
            SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4,1,42,kTRUE, 31);
            labelWeightsEnergy->SetTextFont(43);
            labelWeightsEnergy->Draw();
            TLatex *labelWeightsEtaPrime = new TLatex(0.95,0.20,"#eta' #rightarrow #gamma#gamma");
            SetStyleTLatex( labelWeightsEtaPrime, 0.85*textSizeLabelsPixel,4,1,42,kTRUE, 31);
            labelWeightsEtaPrime->SetTextFont(43);
            labelWeightsEtaPrime->Draw();
            TLatex *labelDetProcWeights    = new TLatex(0.95, 0.16,detectionProcess.Data());
            SetStyleTLatex( labelDetProcWeights, 0.85*textSizeLabelsPixel,4,1,42,kTRUE, 31);
            labelDetProcWeights->SetTextFont(43);
            labelDetProcWeights->Draw();


    //      DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
            DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/%s_WeightsEtaPrimeTriggers.%s",outputDir.Data(), isMC.Data(), suffix.Data()));
        delete canvasWeights;

        // Calculating relative error for etaprime
        for (Int_t i = 0; i < 12; i++){
            if (histoStatEtaPrime[i])
                histoRelStatEtaPrime[i]      = CalculateRelErrUpTH1D( histoStatEtaPrime[i], Form("relativeStatErrorEtaPrime_%s", nameTriggerWeighted[i].Data()));
            if (graphSystEtaPrime[i])
                graphRelSystEtaPrime[i]      = CalculateRelErrUpAsymmGraph( graphSystEtaPrime[i], Form("relativeSysErrorEtaPrime_%s", nameTriggerWeighted[i].Data()));
        }

        TGraphAsymmErrors* graphRelErrorEtaPrimeTot        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaPrimeTot, "relativeTotalErrorEtaPrime");
        while (graphRelErrorEtaPrimeTot->GetY()[0] < 0 ) graphRelErrorEtaPrimeTot->RemovePoint(0);

        TGraphAsymmErrors* graphRelErrorEtaPrimeStat       = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaPrimeStat, "relativeStatErrorEtaPrime");
        while (graphRelErrorEtaPrimeStat->GetY()[0] < 0 ) graphRelErrorEtaPrimeStat->RemovePoint(0);

        TGraphAsymmErrors* graphRelErrorEtaPrimeSys        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaPrimeSys, "relativeSysErrorEtaPrime");
        while (graphRelErrorEtaPrimeSys->GetY()[0] < 0 ) graphRelErrorEtaPrimeSys->RemovePoint(0);

        const char *SysErrDatnameMeanSingleErrCheck = Form("%s/SystematicErrorAveragedSingle%s_EtaPrime_%s_Check.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
        fstream SysErrDatAverSingleCheck;
        SysErrDatAverSingleCheck.precision(4);
        cout << SysErrDatnameMeanSingleErrCheck << endl;
        if(sysAvailSingleEtaPrime[0]){
          SysErrDatAverSingleCheck.open(SysErrDatnameMeanSingleErrCheck, ios::out);
          SysErrDatAverSingleCheck << "pt \t Stat err \t sys err \t tot err " << endl;
          for (Int_t i = 0; i < graphRelErrorEtaPrimeTot->GetN(); i++){
              if (graphRelErrorEtaPrimeStat->GetY()[i] > 0) SysErrDatAverSingleCheck << graphRelErrorEtaPrimeStat->GetX()[i] << "\t" << graphRelErrorEtaPrimeStat->GetY()[i] <<"\t" << graphRelErrorEtaPrimeSys->GetY()[i] <<  "\t" << graphRelErrorEtaPrimeTot->GetY()[i] << endl;
          }
          SysErrDatAverSingleCheck << endl;
          SysErrDatAverSingleCheck.close();
        }

        // plot sys relative errors for individual triggers
        TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);

        TH2F * histo2DRelSysErr;
        histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.,maxPtGlobalEtaPrime,1000,0,60.5);
        SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelSysErr->Draw("copy");
            TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.035*(nMeasSetEtaPrime+1)/2), 0.95, 0.92, 32);
            legendRelSysErr->SetNColumns(2);
            for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                cout << "plotting graph: " << availableMeasEtaPrime[i] << "\t" <<graphRelSystEtaPrime[availableMeasEtaPrime[i]]->GetName() << endl;
                DrawGammaSetMarkerTGraph(graphRelSystEtaPrime[availableMeasEtaPrime[i]], markerTriggWeighted[availableMeasEtaPrime[i]], sizeTrigg[availableMeasEtaPrime[i]],
                                         colorTriggWeighted[availableMeasEtaPrime[i]], colorTriggWeighted[availableMeasEtaPrime[i]]);
                graphRelSystEtaPrime[availableMeasEtaPrime[i]]->Draw("p,same,z");
                legendRelSysErr->AddEntry(graphRelSystEtaPrime[availableMeasEtaPrime[i]],nameTriggerWeighted[availableMeasEtaPrime[i]],"p");
            }
            legendRelSysErr->Draw();

            TLatex *labelRelErrEnergy    = new TLatex(0.15,0.89,collisionSystem.Data());
            SetStyleTLatex( labelRelErrEnergy, 0.85*textSizeLabelsPixel,4);
            labelRelErrEnergy->SetTextFont(43);
            labelRelErrEnergy->Draw();
            TLatex *labelRelErrEtaPrime       = new TLatex(0.15,0.85,"#eta' #rightarrow #gamma#gamma");
            SetStyleTLatex( labelRelErrEtaPrime, 0.85*textSizeLabelsPixel,4);
            labelRelErrEtaPrime->SetTextFont(43);
            labelRelErrEtaPrime->Draw();
            TLatex *labelDetProcRelErr    = new TLatex(0.15, 0.81,detectionProcess.Data());
            SetStyleTLatex( labelDetProcRelErr, 0.85*textSizeLabelsPixel,4);
            labelDetProcRelErr->SetTextFont(43);
            labelDetProcRelErr->Draw();

        canvasRelSysErr->SaveAs(Form("%s/EtaPrime_RelSysErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));


            DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaPrimeSys, 24, 1.5, kGray+1 , kGray+1);
            // graphRelErrorEtaPrimeSys->SetLineStyle(7);
            graphRelErrorEtaPrimeSys->Draw("same,pze1");
            legendRelSysErr->AddEntry(graphRelErrorEtaPrimeSys,"average","p");
            legendRelSysErr->Draw();

            labelRelErrEnergy->Draw();
            labelRelErrEtaPrime->Draw();
            labelDetProcRelErr->Draw();

        canvasRelSysErr->SaveAs(Form("%s/EtaPrime_RelSysErrWithAverage_SingleMeas.%s",outputDir.Data(),suffix.Data()));


        // plot stat relative errors for individual triggers
        TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);

        TH2F * histo2DRelStatErr;
        histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.,maxPtGlobalEtaPrime,1000,0,60.5);
        SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelStatErr->Draw("copy");

            TLegend* legendRelStatErr       = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaPrime/2), 0.95, 0.92, 32);
            legendRelStatErr->SetNColumns(2);
            for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                 cout << "plotting graph: " << availableMeasEtaPrime[i] << "\t" <<histoRelStatEtaPrime[availableMeasEtaPrime[i]]->GetName() << endl;
                if (histoRelStatEtaPrime[availableMeasEtaPrime[i]] && modeNormal == 2){
                    TGraphAsymmErrors* dummyGraph = new TGraphAsymmErrors(histoRelStatEtaPrime[availableMeasEtaPrime[i]]);
                    dummyGraph->Print();
                    DrawGammaSetMarkerTGraph(dummyGraph, markerTriggWeighted[availableMeasEtaPrime[i]], sizeTrigg[availableMeasEtaPrime[i]],
                                         colorTriggWeighted[availableMeasEtaPrime[i]], colorTriggWeighted[availableMeasEtaPrime[i]]);
                    dummyGraph->Draw("pX,same");
                    legendRelStatErr->AddEntry(dummyGraph,nameTriggerWeighted[availableMeasEtaPrime[i]],"p");

                     for (Int_t j = 1; j < histoRelStatEtaPrime[availableMeasEtaPrime[i]]->GetNbinsX()+1; j++){
                        cout << j << ": " << histoRelStatEtaPrime[availableMeasEtaPrime[i]]->GetBinContent(j) << endl;
                     }
                } else {
                    DrawGammaSetMarker(histoRelStatEtaPrime[availableMeasEtaPrime[i]],markerTriggWeighted[availableMeasEtaPrime[i]], sizeTrigg[availableMeasEtaPrime[i]],
                                            colorTriggWeighted[availableMeasEtaPrime[i]], colorTriggWeighted[availableMeasEtaPrime[i]]);
                    histoRelStatEtaPrime[availableMeasEtaPrime[i]]->DrawCopy("p,same,z");
                    legendRelStatErr->AddEntry(histoRelStatEtaPrime[availableMeasEtaPrime[i]],nameTriggerWeighted[availableMeasEtaPrime[i]],"p");
                }
            }
            legendRelStatErr->Draw();

            labelRelErrEnergy->Draw();
            labelRelErrEtaPrime->Draw();
            labelDetProcRelErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/EtaPrime_RelStatErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));

        // plot full error for final result decomposed
        TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);

        TH2F * histo2DRelTotErrEtaPrime;
        histo2DRelTotErrEtaPrime                 = new TH2F("histo2DRelTotErrEtaPrime","histo2DRelTotErrEtaPrime",11000,0.,maxPtGlobalEtaPrime,1000,0,40.5);
        SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaPrime, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelTotErrEtaPrime->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaPrimeTot, 20, 1.5, kBlack , kBlack);
            graphRelErrorEtaPrimeTot->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaPrimeStat, 24, 1.5, kGray+2 , kGray+2);
            graphRelErrorEtaPrimeStat->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaPrimeSys, 24, 1.5, kGray+1 , kGray+1);
            graphRelErrorEtaPrimeSys->SetLineStyle(7);
            graphRelErrorEtaPrimeSys->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr2       = GetAndSetLegend2(0.72, 0.92-(0.035*3), 0.9, 0.92, 32);
            legendRelTotErr2->AddEntry(graphRelErrorEtaPrimeTot,"tot","p");
            legendRelTotErr2->AddEntry(graphRelErrorEtaPrimeStat,"stat","l");
            legendRelTotErr2->AddEntry(graphRelErrorEtaPrimeSys,"sys","l");
            legendRelTotErr2->Draw();

            labelRelErrEnergy->Draw();
            labelRelErrEtaPrime->Draw();
            labelDetProcRelErr->Draw();

        canvasRelTotErr->SaveAs(Form("%s/EtaPrime_RelErrorsFulldecomp.%s",outputDir.Data(),suffix.Data()));


        // Calculate relative sys error weighted
        if (sysAvailSingleEtaPrime[0]){
            for (Int_t k = 0; k< nRelSysErrEtaPrimeSources ; k++ ){
                graphRelSysErrEtaPrimeSourceWeighted[k]      = CalculateWeightedQuantity(    graphOrderedRelSysErrEtaPrimeSource[k],
                                                                                        graphWeightsEtaPrime,
                                                                                        binningEtaPrime,  maxNAllowedEtaPrime,
                                                                                        MaxNumberOfFiles
                                                                                   );
                if (!graphRelSysErrEtaPrimeSourceWeighted[k]){
                    cout << "Aborted in CalculateWeightedQuantity for " << endl;
                    return;
                } else {
                    graphRelSysErrEtaPrimeSourceWeighted[k]->SetName(Form("RelSysErrEtaPrimeSourceWeighted%s", ((TString)ptSysDetail[0][0].at(k+1)).Data()));
                    while (graphRelSysErrEtaPrimeSourceWeighted[k]->GetY()[0] == -10000 )   graphRelSysErrEtaPrimeSourceWeighted[k]->RemovePoint(0);
                }
            }
        }
        // return;

        // Calculation of averaged supporting plots with weights from spectra
        graphMassEtaPrimeDataWeighted                    = CalculateWeightedQuantity(    graphOrderedMassEtaPrimeData,
                                                                                    graphWeightsEtaPrime,
                                                                                    binningEtaPrime,  maxNAllowedEtaPrime,
                                                                                    MaxNumberOfFiles
                                                                                );
        if (!graphMassEtaPrimeDataWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }
        graphMassEtaPrimeMCWeighted                      = CalculateWeightedQuantity(    graphOrderedMassEtaPrimeMC,
                                                                                    graphWeightsEtaPrime,
                                                                                    binningEtaPrime,  maxNAllowedEtaPrime,
                                                                                    MaxNumberOfFiles
                                                                                );
        if (!graphMassEtaPrimeMCWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }

        graphWidthEtaPrimeDataWeighted                   = CalculateWeightedQuantity(    graphOrderedWidthEtaPrimeData,
                                                                                    graphWeightsEtaPrime,
                                                                                    binningEtaPrime,  maxNAllowedEtaPrime,
                                                                                    MaxNumberOfFiles
                                                                                );
        if (!graphWidthEtaPrimeDataWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }

        graphWidthEtaPrimeMCWeighted                     = CalculateWeightedQuantity(    graphOrderedWidthEtaPrimeMC,
                                                                                    graphWeightsEtaPrime,
                                                                                    binningEtaPrime,  maxNAllowedEtaPrime,
                                                                                    MaxNumberOfFiles
                                                                                );
        if (!graphWidthEtaPrimeMCWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }

        cout << "weighting EtaPrime acceptance" << endl;
        graphAcceptanceEtaPrimeWeighted                      = CalculateWeightedQuantity(    graphOrderedAcceptanceEtaPrime,
                                                                                        graphWeightsEtaPrime,
                                                                                        binningEtaPrime,  maxNAllowedEtaPrime,
                                                                                        MaxNumberOfFiles
                                                                                    );
        if (!graphAcceptanceEtaPrimeWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }

        cout << "weighting EtaPrime efficiency" << endl;
        graphEfficiencyEtaPrimeWeighted                      = CalculateWeightedQuantity(    graphOrderedEfficiencyEtaPrime,
                                                                                        graphWeightsEtaPrime,
                                                                                        binningEtaPrime,  maxNAllowedEtaPrime,
                                                                                        MaxNumberOfFiles
                                                                                    );
        if (!graphEfficiencyEtaPrimeWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }
        cout << "weighting EtaPrime efficiency x acceptance" << endl;
        graphEffTimesAccEtaPrimeWeighted                     = CalculateWeightedQuantity(    graphOrderedEffTimesAccEtaPrime,
                                                                                        graphWeightsEtaPrime,
                                                                                        binningEtaPrime,  maxNAllowedEtaPrime,
                                                                                        MaxNumberOfFiles
                                                                                    );

        if (!graphEffTimesAccEtaPrimeWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }

        // remove points in spectrum which should have been masked
        if (graphMassEtaPrimeDataWeighted)
            while (graphMassEtaPrimeDataWeighted->GetY()[0] == -10000 )   graphMassEtaPrimeDataWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted EtaPrime mass data graph" << endl;
        if (graphMassEtaPrimeMCWeighted)
            while (graphMassEtaPrimeMCWeighted->GetY()[0] == -10000)     graphMassEtaPrimeMCWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted EtaPrime mass MC graph" << endl;
        if (graphWidthEtaPrimeDataWeighted)
            while (graphWidthEtaPrimeDataWeighted->GetY()[0] == -10000)  graphWidthEtaPrimeDataWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted EtaPrime width data graph" << endl;
        if (graphWidthEtaPrimeMCWeighted)
            while (graphWidthEtaPrimeMCWeighted->GetY()[0] == -10000)    graphWidthEtaPrimeMCWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted EtaPrime width MC graph" << endl;

        if (graphAcceptanceEtaPrimeWeighted)
            while (graphAcceptanceEtaPrimeWeighted->GetY()[0] == -10000)     graphAcceptanceEtaPrimeWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted EtaPrime acceptance graph" << endl;

        if (graphEfficiencyEtaPrimeWeighted)
            while (graphEfficiencyEtaPrimeWeighted->GetY()[0] == -10000)     graphEfficiencyEtaPrimeWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted EtaPrime efficiency graph" << endl;

        if (graphEffTimesAccEtaPrimeWeighted)
            while (graphEffTimesAccEtaPrimeWeighted->GetY()[0] == -10000)    graphEffTimesAccEtaPrimeWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted EtaPrime acceptance x efficiency graph" << endl;


        //  **********************************************************************************************************************
        //  **************************************** Combine+write detailed Systematics ******************************************
        //  **********************************************************************************************************************

        const char *SysErrDatnameMeanSingleErr = Form("%s/SystematicErrorAveragedSingle%s_EtaPrime_%s.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
        fstream SysErrDatAverSingle;
        SysErrDatAverSingle.precision(4);
        cout << SysErrDatnameMeanSingleErr << endl;
        if(sysAvailSingleEtaPrime[0] && graphRelSysErrEtaPrimeSourceWeighted[0]){
            SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
            for(Int_t iColumn = 0; iColumn < (Int_t)ptSysDetail[0][0].size(); iColumn++) SysErrDatAverSingle << ptSysDetail[0][0].at(iColumn) << "\t";
            SysErrDatAverSingle << endl;
            for (Int_t i = 1; i < nrOfTrigToBeComb; i++){
                if(!sysAvailSingleEtaPrime[i]) continue;
                for(Int_t iCol = 0; iCol < (Int_t)ptSysDetail[i][0].size(); iCol++){
                    if( ((TString)ptSysDetail[i][0].at(iCol)).CompareTo(((TString)ptSysDetail[0][0].at(iCol))) ){
                        cout << "ERROR: Systematic error type at pos " << iCol << " does not agree for " << availableMeasEtaPrime[i] << " & " << availableMeasEtaPrime[0] << ", returning!" << endl;
                        return;
                    }
                }
            }

            for(Int_t i=0; i<graphRelSysErrEtaPrimeSourceWeighted[0]->GetN(); i++){
                SysErrDatAverSingle << graphRelSysErrEtaPrimeSourceWeighted[0]->GetX()[i] << "\t";
                Int_t nColumns = (Int_t)ptSysDetail[0][0].size();
                for(Int_t iErr=0; iErr<nColumns-1; iErr++)
                    SysErrDatAverSingle << graphRelSysErrEtaPrimeSourceWeighted[iErr]->GetY()[i] << "\t";
                SysErrDatAverSingle << endl;

            }
        }
        SysErrDatAverSingle.close();

        // ***************************************************************************************************
        // ********************* Plot all mean erros separately after smoothing ******************************
        // ***************************************************************************************************
        if(sysAvailSingleEtaPrime[0] && graphRelSysErrEtaPrimeSourceWeighted[0]){
            TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);// gives the page size
            DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);

                // create dummy histo
                TH2D *histo2DNewSysErrMean ;
                histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 100,0.,maxPtGlobalEtaPrime,1000.,-0.5,30.);
                SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                        1,0.9, 510, 510);
                histo2DNewSysErrMean->Draw();

                // Give legend position for plotting
                Double_t minXLegend     = 0.12;
                Double_t maxYLegend     = 0.95;
                Double_t widthLegend    = 0.25;
                if (nRelSysErrEtaPrimeSources> 7)
                    widthLegend         = 0.5;
                Double_t heightLegend   = 1.05* 0.035 * (nRelSysErrEtaPrimeSources+3);
                if (nRelSysErrEtaPrimeSources> 7)
                    heightLegend        = 1.05* 0.035 * (nRelSysErrEtaPrimeSources/2+1);

                // create legend
                TLegend* legendMeanNew = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
                legendMeanNew->SetMargin(0.1);
                if (nRelSysErrEtaPrimeSources> 7) legendMeanNew->SetNColumns(2);

                for(Int_t i = 0;i< nRelSysErrEtaPrimeSources-1 ;i++){
                    DrawGammaSetMarkerTGraphAsym(graphRelSysErrEtaPrimeSourceWeighted[i], GetMarkerStyleSystematics(  ptSysDetail[0][0].at(i+1), modeNormal), 1.,
                                                GetColorSystematics( ptSysDetail[0][0].at(i+1), modeNormal),GetColorSystematics( ptSysDetail[0][0].at(i+1), modeNormal));
                    graphRelSysErrEtaPrimeSourceWeighted[i]->Draw("pX0,csame");
                    legendMeanNew->AddEntry(graphRelSysErrEtaPrimeSourceWeighted[i],GetSystematicsName(ptSysDetail[0][0].at(i+1)),"p");
                }

                DrawGammaSetMarkerTGraphAsym(graphRelSysErrEtaPrimeSourceWeighted[nRelSysErrEtaPrimeSources-1], 20, 1.,kBlack,kBlack);
                graphRelSysErrEtaPrimeSourceWeighted[nRelSysErrEtaPrimeSources-1]->Draw("p,csame");
                legendMeanNew->AddEntry(graphRelSysErrEtaPrimeSourceWeighted[nRelSysErrEtaPrimeSources-1],"quad. sum.","p");
                legendMeanNew->Draw();

                // labeling
                TLatex *labelEnergySysDetailed = new TLatex(0.95, 0.93,collisionSystem.Data());
                SetStyleTLatex( labelEnergySysDetailed, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
                labelEnergySysDetailed->Draw();

                TLatex *labelEtaPrimeSysDetailed     = new TLatex(0.95, 0.93-0.99*textSizeSpectra*0.85,"#eta' #rightarrow #gamma#gamma");
                SetStyleTLatex( labelEtaPrimeSysDetailed, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
                labelEtaPrimeSysDetailed->Draw();

                TLatex *labelDetProcSysDetailed = new TLatex(0.95, 0.93-2*0.99*textSizeSpectra*0.85,detectionProcess.Data());
                SetStyleTLatex( labelDetProcSysDetailed, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
                labelDetProcSysDetailed->Draw();

            canvasNewSysErrMean->Update();
            canvasNewSysErrMean->SaveAs(Form("%s/EtaPrime_SysErrorsSeparatedSourcesReweighted_%s.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }

        for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++){
            for(Int_t iB=0; iB<50; iB++) ptSysDetail[iR][iB].clear();
        }


    // if averaging wasn't enabled pick values according to predefined ranges ("cherry picking points")
    } else {
        graphCorrectedYieldWeightedAverageEtaPrimeStat        = new TGraphAsymmErrors(nPointFinalEtaPrime, xValueFinalEtaPrime, yValueFinalEtaPrime,
                                                                                    xErrorLowFinalEtaPrime, xErrorHighFinalEtaPrime,yErrorLowFinalEtaPrime, yErrorHighFinalEtaPrime);
        graphCorrectedYieldWeightedAverageEtaPrimeSys         = new TGraphAsymmErrors(nPointFinalEtaPrime, xValueFinalEtaPrime, yValueFinalEtaPrime,
                                                                                xErrorLowFinalEtaPrime, xErrorHighFinalEtaPrime,yErrorSysLowFinalEtaPrime, yErrorSysHighFinalEtaPrime);

        graphMassEtaPrimeDataWeighted                    = graphMassEtaPrimeData[0];
        graphMassEtaPrimeMCWeighted                      = graphMassEtaPrimeMC[0];
        graphWidthEtaPrimeDataWeighted                   = graphWidthEtaPrimeData[0];
        graphWidthEtaPrimeMCWeighted                     = graphWidthEtaPrimeMC[0];

        graphAcceptanceEtaPrimeWeighted                      = graphAcceptanceEtaPrime[0];
        graphEfficiencyEtaPrimeWeighted                      = graphEfficiencyEtaPrime[0];
        graphEffTimesAccEtaPrimeWeighted                     = graphEffTimesAccEtaPrime[0];

        TGraphAsymmErrors* graphRelErrorEtaPrimeStat       = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaPrimeStat, "relativeStatErrorEtaPrime");
        while (graphRelErrorEtaPrimeStat->GetY()[0] < 0 ) graphRelErrorEtaPrimeStat->RemovePoint(0);

        TGraphAsymmErrors* graphRelErrorEtaPrimeSys        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaPrimeSys, "relativeSysErrorEtaPrime");
        while (graphRelErrorEtaPrimeSys->GetY()[0] < 0 ) graphRelErrorEtaPrimeSys->RemovePoint(0);

        const char *SysErrDatnameMeanSingleErrCheck = Form("%s/SystematicErrorAveragedSingle%s_EtaPrime_%s_Check.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
        fstream SysErrDatAverSingleCheck;
        SysErrDatAverSingleCheck.precision(4);
        cout << SysErrDatnameMeanSingleErrCheck << endl;

        SysErrDatAverSingleCheck.open(SysErrDatnameMeanSingleErrCheck, ios::out);
        SysErrDatAverSingleCheck << "pt \t Stat err \t sys err \t tot err " << endl;
        for (Int_t i = 0; i < graphRelErrorEtaPrimeStat->GetN(); i++){
            if (graphRelErrorEtaPrimeStat->GetY()[i] > 0) SysErrDatAverSingleCheck << graphRelErrorEtaPrimeStat->GetX()[i] << "\t" << graphRelErrorEtaPrimeStat->GetY()[i] <<"\t" << graphRelErrorEtaPrimeSys->GetY()[i] << endl;
        }
        SysErrDatAverSingleCheck << endl;
        SysErrDatAverSingleCheck.close();

    }
    // print final graphs
    cout << "stat etaprime" << endl;
    graphCorrectedYieldWeightedAverageEtaPrimeStat->Print();
    cout << "sys etaprime" << endl;
    if (graphCorrectedYieldWeightedAverageEtaPrimeSys) graphCorrectedYieldWeightedAverageEtaPrimeSys->Print();

    if (graphCorrectedYieldWeightedAverageEtaPrimeStat){
        while (graphCorrectedYieldWeightedAverageEtaPrimeStat->GetX()[0]< minPtGlobalEtaPrime){
            graphCorrectedYieldWeightedAverageEtaPrimeStat->RemovePoint(0);
        }
    }

    if (graphCorrectedYieldWeightedAverageEtaPrimeSys){
        while (graphCorrectedYieldWeightedAverageEtaPrimeSys->GetX()[0]< minPtGlobalEtaPrime){
            graphCorrectedYieldWeightedAverageEtaPrimeSys->RemovePoint(0);
        }
    }

    //***************************************************************************************************************
    //************************************Plotting Mass EtaPrime reduced range  ******************************************
    //***************************************************************************************************************
    canvasMass->cd();
    histo2DMassEtaPrime->DrawCopy();

    legendMassRedEtaPrime->SetNColumns(2);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if((optionEnergy.CompareTo("2.76TeV")==0 ) ||
            (optionEnergy.CompareTo("8TeV")==0) ||
            (optionEnergy.CompareTo("pPb_5.023TeV")==0) ||
            (optionEnergy.CompareTo("XeXe_5.44TeV")==0)
            ){
            if (graphMassEtaPrimeData[i] && !maskedFullyEtaPrime[i]) {
                DrawGammaSetMarkerTGraphAsym(graphMassEtaPrimeData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                graphMassEtaPrimeData[i]->Draw("p,e1,same");
                legendMassRedEtaPrime->AddEntry(graphMassEtaPrimeData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
            }
            if (graphMassEtaPrimeMC[i] && !maskedFullyEtaPrime[i]){
                DrawGammaSetMarkerTGraphAsym(graphMassEtaPrimeMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                graphMassEtaPrimeMC[i]->Draw("p,e1,same");
                legendMassRedEtaPrime->AddEntry(graphMassEtaPrimeMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
            }
        }
    }
    legendMassRedEtaPrime->Draw();
    labelEnergyMass->Draw();
    labelEtaPrimeMass->Draw();
    labelDetProcMass->Draw();

    canvasMass->Update();
    canvasMass->SaveAs(Form("%s/EtaPrime_%s_Mass_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    //***************************************************************************************************************
    //********************************* EtaPrime Mass weighted ***********************************************************
    //***************************************************************************************************************
    if (graphMassEtaPrimeDataWeighted || graphMassEtaPrimeMCWeighted){
        canvasMass->cd();
        histo2DMassEtaPrime->DrawCopy();
        TLegend* legendMassEtaPrimeWeighted = GetAndSetLegend2(0.52, 0.88, 0.75, 0.88+(1.05*2*0.85*textSizeSpectra),28);
        if (graphMassEtaPrimeDataWeighted){
            DrawGammaSetMarkerTGraphAsym(graphMassEtaPrimeDataWeighted, 20, 1, kBlack, kBlack);
            graphMassEtaPrimeDataWeighted->Draw("p,e1,same");
            legendMassEtaPrimeWeighted->AddEntry(graphMassEtaPrimeDataWeighted, "Data", "p");
        }
        if (graphMassEtaPrimeDataWeighted){
            DrawGammaSetMarkerTGraphAsym(graphMassEtaPrimeMCWeighted, 24, 1, kGray+2, kGray+2);
            graphMassEtaPrimeMCWeighted->Draw("p,e1,same");
            legendMassEtaPrimeWeighted->AddEntry(graphMassEtaPrimeMCWeighted, "MC", "p");
        }
        legendMassEtaPrimeWeighted->Draw();
        labelEnergyMass->Draw();
        labelEtaPrimeMass->Draw();
        labelDetProcMass->Draw();

        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/EtaPrime_%s_Mass_Weighted.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        TGraphAsymmErrors* graphMassDifferenceEtaPrimeDatavsMC       = CalculateAsymGraphDifferenceToGraph(graphMassEtaPrimeDataWeighted,graphMassEtaPrimeMCWeighted);
        graphMassDifferenceEtaPrimeDatavsMC->Print();

        TGraphAsymmErrors* graphMassRelDifferenceEtaPrimeDatavsMC    = CalculateAsymGraphRatioToGraph(graphMassDifferenceEtaPrimeDatavsMC,graphMassEtaPrimeMCWeighted);
        graphMassRelDifferenceEtaPrimeDatavsMC->Print();
        graphMassRelDifferenceEtaPrimeDatavsMC = ScaleGraph(graphMassRelDifferenceEtaPrimeDatavsMC, 100.);

        canvasMass->cd();
        TH2F * histo2DRelMassDiffEtaPrime       = new TH2F("histo2DRelMassDiffEtaPrime","histo2DRelMassDiffEtaPrime",1000,0., maxPtGlobalEtaPrime,10000,-3, 3);
        SetStyleHistoTH2ForGraphs(histo2DRelMassDiffEtaPrime, "#it{p}_{T} (GeV/#it{c})","#it{M}_{#eta', data}-#it{M}_{#eta', MC}/ #it{M}_{#eta', MC} (%)",
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.4);
        histo2DRelMassDiffEtaPrime->DrawCopy();

        if (graphMassRelDifferenceEtaPrimeDatavsMC){
            TF1 *fitEtaPrimeRelMassDiff = new TF1("fitEtaPrimeRelMassDiff","[0]",0.,maxPtGlobalEtaPrime);
            fitEtaPrimeRelMassDiff->SetParameter(0,0.);
            graphMassRelDifferenceEtaPrimeDatavsMC->Fit(fitEtaPrimeRelMassDiff,"QNRMEX0+","",0.,maxPtGlobalEtaPrime);
            DrawGammaSetMarkerTF1( fitEtaPrimeRelMassDiff, 1, 0.5, kRed);
            fitEtaPrimeRelMassDiff->Draw("same");

            TLegend* legendEtaPrimeRelMassDiff = GetAndSetLegend2(0.15, 0.12, 0.6, 0.12+(0.035*1), 0.035, 2, "", 42, 0.15);
            legendEtaPrimeRelMassDiff->AddEntry(fitEtaPrimeRelMassDiff,"fit with constant");
            legendEtaPrimeRelMassDiff->AddEntry((TObject*)0,Form("%0.4f #pm %0.4f", fitEtaPrimeRelMassDiff->GetParameter(0), fitEtaPrimeRelMassDiff->GetParError(0)),"");
            legendEtaPrimeRelMassDiff->Draw();

            DrawGammaSetMarkerTGraphAsym(graphMassRelDifferenceEtaPrimeDatavsMC, 20, 1, kBlack, kBlack);
            graphMassRelDifferenceEtaPrimeDatavsMC->Draw("p,e1,same");

            fileFitsOutput << "average rel mass diff: " << fitEtaPrimeRelMassDiff->GetParameter(0) << "+-"<< fitEtaPrimeRelMassDiff->GetParError(0) << endl;
        }
        DrawGammaLines(0., maxPtGlobalEtaPrime , 0., 0., 1, kGray+2, 7);
            // legendMassEtaPrimeWeighted->Draw();
        labelEnergyMass->Draw();
        labelEtaPrimeMass->Draw();
        labelDetProcMass->Draw();

        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/EtaPrime_%s_RelMassDiff_Weighted.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    }
    //***************************************************************************************************************
    //************************************Plotting Width EtaPrime reduced range  *****************************************
    //***************************************************************************************************************
    canvasWidth->cd();
    histo2DWidthEtaPrime->DrawCopy();

    legendWidthRedEtaPrime->SetNColumns(2);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if((optionEnergy.CompareTo("2.76TeV")==0 ) ||
            (optionEnergy.CompareTo("8TeV")==0) ||
            (optionEnergy.CompareTo("pPb_5.023TeV")==0) ||
            (optionEnergy.CompareTo("XeXe_5.44TeV")==0 )

            ){
            if (graphWidthEtaPrimeData[i] && !maskedFullyEtaPrime[i]) {
                DrawGammaSetMarkerTGraphAsym(graphWidthEtaPrimeData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                graphWidthEtaPrimeData[i]->Draw("p,e1,same");
                legendWidthRedEtaPrime->AddEntry(graphWidthEtaPrimeData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
            }
            if (graphWidthEtaPrimeMC[i] && !maskedFullyEtaPrime[i]){
                DrawGammaSetMarkerTGraphAsym(graphWidthEtaPrimeMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                graphWidthEtaPrimeMC[i]->Draw("p,e1,same");
                legendWidthRedEtaPrime->AddEntry(graphWidthEtaPrimeMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
            }
        }
    }
    legendWidthRedEtaPrime->Draw();
    labelEnergyWidth->Draw();
    labelEtaPrimeWidth->Draw();
    labelDetProcWidth->Draw();

    canvasWidth->Update();
    canvasWidth->SaveAs(Form("%s/EtaPrime_%s_Width_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    //***************************************************************************************************************
    //********************************* EtaPrime Width weighted **********************************************************
    //***************************************************************************************************************
    canvasWidth->cd();
    histo2DWidthEtaPrime->DrawCopy();
    TLegend* legendWidthEtaPrimeWeighted = GetAndSetLegend2(0.52, 0.86, 0.75, 0.86+(1.05*2*0.85*textSizeSpectra),28);
    if (graphWidthEtaPrimeDataWeighted){
        DrawGammaSetMarkerTGraphAsym(graphWidthEtaPrimeDataWeighted, 20, 1, kBlack, kBlack);
        graphWidthEtaPrimeDataWeighted->Draw("p,e1,same");
        legendWidthEtaPrimeWeighted->AddEntry(graphWidthEtaPrimeDataWeighted, "Data", "p");
    }
    if (graphWidthEtaPrimeDataWeighted){
        DrawGammaSetMarkerTGraphAsym(graphWidthEtaPrimeMCWeighted, 24, 1, kGray+2, kGray+2);
        graphWidthEtaPrimeMCWeighted->Draw("p,e1,same");
        legendWidthEtaPrimeWeighted->AddEntry(graphWidthEtaPrimeMCWeighted, "MC", "p");
    }
    legendWidthEtaPrimeWeighted->Draw();
    labelEnergyWidth->Draw();
    labelEtaPrimeWidth->Draw();
    labelDetProcWidth->Draw();

    canvasWidth->Update();
    canvasWidth->SaveAs(Form("%s/EtaPrime_%s_Width_Weighted.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
    delete canvasMass;
    delete canvasWidth;

    //***************************************************************************************************************
    //************************************* Efficiency weighted *****************************************************
    //***************************************************************************************************************
    if (graphEfficiencyEtaPrimeWeighted){
        DrawGammaCanvasSettings( canvasEffi, 0.09, 0.017, 0.015, 0.08);
        canvasEffi->SetLogy(1);
        canvasEffi->cd();
        histo2DEffiEtaPrime->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphEfficiencyEtaPrimeWeighted, 20, 1, kGray+2, kGray+2);
        graphEfficiencyEtaPrimeWeighted->Draw("p,e1,same");

        TLatex *labelEnergyEffiWOTrigg = new TLatex(0.95, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyEffiWOTrigg, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
        labelEnergyEffiWOTrigg->Draw();

        TLatex *labelEtaPrimeEffiWOTrigg = new TLatex(0.95, 0.15+0.99*textSizeSpectra*0.85,"#eta' #rightarrow #gamma#gamma");
        SetStyleTLatex( labelEtaPrimeEffiWOTrigg, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
        labelEtaPrimeEffiWOTrigg->Draw();

        TLatex *labelDetProcEffiWOTrigg = new TLatex(0.95, 0.15,detectionProcess.Data());
        SetStyleTLatex( labelDetProcEffiWOTrigg, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
        labelDetProcEffiWOTrigg->Draw();

        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/EtaPrime_EfficiencyW0TriggEff_Weighted.%s",outputDir.Data(),suffix.Data()));
    }

    //***************************************************************************************************************
    //************************************* Efficiency weighted *****************************************************
    //***************************************************************************************************************
    if (graphEffTimesAccEtaPrimeWeighted){
      DrawGammaCanvasSettings( canvasEffi, 0.09, 0.017, 0.015, 0.08);
      canvasEffi->SetLogy(1);
      canvasEffi->SetLogx(1);
      canvasEffi->cd();
      TH2F * histo2DAccEff;
      histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.23,  maxPtGlobalEtaPrime*2, 1000, 8e-5, 2e-0 );
      SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                                 0.85*textSizeSpectra, textSizeSpectra, 0.85*textSizeSpectra, textSizeSpectra, 0.9, 1.04);//(#times #epsilon_{pur})
      histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
      histo2DAccEff->GetXaxis()->SetLabelOffset(-0.01);
      histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
      histo2DAccEff->DrawCopy();

      DrawGammaSetMarkerTGraphAsym(graphEffTimesAccEtaPrimeWeighted, 20, 1, kGray+2, kGray+2);
      graphEffTimesAccEtaPrimeWeighted->Draw("p,e1,same");

      TLatex *labelEnergyEffiWOTrigg = new TLatex(0.95, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
      SetStyleTLatex( labelEnergyEffiWOTrigg, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
      labelEnergyEffiWOTrigg->Draw();

      TLatex *labelEtaPrimeEffiWOTrigg = new TLatex(0.95, 0.15+0.99*textSizeSpectra*0.85,"#eta' #rightarrow #gamma#gamma");
      SetStyleTLatex( labelEtaPrimeEffiWOTrigg, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
      labelEtaPrimeEffiWOTrigg->Draw();

      TLatex *labelDetProcEffiWOTrigg = new TLatex(0.95, 0.15,detectionProcess.Data());
      SetStyleTLatex( labelDetProcEffiWOTrigg, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
      labelDetProcEffiWOTrigg->Draw();

      canvasEffi->Update();
      canvasEffi->SaveAs(Form("%s/EtaPrime_EfficiencyTimesAcceptanceW0TriggEff_Weighted.%s",outputDir.Data(),suffix.Data()));
      canvasEffi->SetLogx(0);
    }

    //***************************************************************************************************************
    //************************************* Acceptance weighted *****************************************************
    //***************************************************************************************************************
    if (graphAcceptanceEtaPrimeWeighted){
        DrawGammaCanvasSettings( canvasAcc, 0.1, 0.017, 0.015, 0.08);
        canvasAcc->cd();
        canvasAcc->SetLogy(0);
        histo2DAccEtaPrime->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphAcceptanceEtaPrimeWeighted, 20, 1, kGray+2, kGray+2);
        graphAcceptanceEtaPrimeWeighted->Draw("p,e1,same");

        TLatex *labelEnergyAcc = new TLatex(0.95, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyAcc, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
        labelEnergyAcc->Draw();

        TLatex *labelEtaPrimeAcc = new TLatex(0.95, 0.15+0.99*textSizeSpectra*0.85,"#eta' #rightarrow #gamma#gamma");
        SetStyleTLatex( labelEtaPrimeAcc, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
        labelEtaPrimeAcc->Draw();

        TLatex *labelDetProcAcc = new TLatex(0.95, 0.15,detectionProcess.Data());
        SetStyleTLatex( labelDetProcAcc, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
        labelDetProcAcc->Draw();

        canvasAcc->Update();
        canvasAcc->SaveAs(Form("%s/EtaPrime_Acceptance_weighted.%s",outputDir.Data(),suffix.Data()));
    }
    // delete canvasAcc;

    //***************************************************************************************************************
    //************************************Plotting scaled invariant yield *****************************************
    //***************************************************************************************************************
    TCanvas* canvasCorrScaled = new TCanvas("canvasCorrScaled","",0,0,1000,1350);// gives the page size
    DrawGammaCanvasSettings( canvasCorrScaled, 0.15, 0.017, 0.015, 0.07);
    canvasCorrScaled->SetLogy();
    // canvasCorrScaled->SetGridx();
    Double_t minCorrYield       = 2e-10;
    Double_t maxCorrYield       = 1e0;
    if (modeNormal == 0){
        minCorrYield            = 2e-9;
        maxCorrYield            = 1e1;
    }

    if(optionEnergy.CompareTo("8TeV")==0){
      if(modeNormal == 2){
        minCorrYield       = 1e-10;
        maxCorrYield       = 1;
      }else if(modeNormal == 4){
        minCorrYield       = 1e-8;
        maxCorrYield       = 0.2;
      }
    }

    if(optionEnergy.CompareTo("pPb_5.023TeV")==0){
      if(modeNormal == 2){
        minCorrYield       = 1e-8;
        maxCorrYield       = 1;
      }else if(modeNormal == 4){
        minCorrYield       = 1e-8;
        maxCorrYield       = 0.2;
      }else if(modeNormal == 3){
        minCorrYield       = 1e-7;
        maxCorrYield       = 10;
      }
    }



    TH2F * histo2DInvYieldScaled;
    histo2DInvYieldScaled = new TH2F("histo2DInvYieldScaled","histo2DInvYieldScaled",1000,0., maxPtGlobalEtaPrime,10000,minCorrYield,maxCorrYield);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldScaled, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.55);
    histo2DInvYieldScaled->DrawCopy();

    Double_t factorEtaPrime              = 10000;
    if (modeNormal == 2) factorEtaPrime        = 10000;
    if (modeNormal == 4) factorEtaPrime        = 5000;

    // two component modeNormall fit
    Double_t paramTCM[5] = {graphCorrectedYieldWeightedAverageEtaPrimeStat->GetY()[0],0.3,graphCorrectedYieldWeightedAverageEtaPrimeStat->GetY()[0]/factorEtaPrime,0.8,3};
    if(modeNormal == 4 && optionEnergy.CompareTo("8TeV")==0){
      paramTCM[0]=0.008; paramTCM[1]=0.65; paramTCM[2]=1.72; paramTCM[3]=0.47; paramTCM[4]=2.93;
    }
    TF1* fitInvYieldEtaPrime = FitObject("tcm","fitInvYieldEtaPrime","EtaPrime",graphCorrectedYieldWeightedAverageEtaPrimeStat,minPtGlobalEtaPrime,maxPtGlobalEtaPrime,paramTCM,"QNRMEX0+");

    // Tsallis fit
    // Double_t paramGraph[3]                  = {1000, 8., 0.13};
    // TF1* fitInvYieldEtaPrime                     = FitObject("l","fitInvYieldEtaPrime","EtaPrime",graphCorrectedYieldWeightedAverageEtaPrimeStat,minPtGlobalEtaPrime,maxPtGlobalEtaPrime,paramGraph,"QNRME+");

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageEtaPrimeSys, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
    graphCorrectedYieldWeightedAverageEtaPrimeSys->Draw("p,E2,same");

    TLegend* legendScaled = GetAndSetLegend2(0.72, 0.95-(1.15*nrOfTrigToBeCombEtaPrimeRed*0.85*textSizeSpectra), 0.95, 0.95,32);

    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if (graphsCorrectedYieldSysShrunkEtaPrime[i]) DrawGammaSetMarkerTGraphAsym(graphsCorrectedYieldSysShrunkEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
        if (graphsCorrectedYieldSysShrunkEtaPrime[i])DrawGammaSetMarker(histoCorrectedYieldEtaPrimeScaled[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        if ( !maskedFullyEtaPrime[i] )histoCorrectedYieldEtaPrimeScaled[i]->DrawCopy("e1,same");
        if ( !maskedFullyEtaPrime[i] )graphsCorrectedYieldSysShrunkEtaPrime[i]->Draw("p,E2,same");
        if (graphsCorrectedYieldSysShrunkEtaPrime[i] && !maskedFullyEtaPrime[i])legendScaled->AddEntry(histoCorrectedYieldEtaPrimeScaled[i],triggerNameLabel[i].Data(),"p");
    }

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageEtaPrimeStat, 24, 2, kRed , kRed, 1, kTRUE);
    graphCorrectedYieldWeightedAverageEtaPrimeStat->Draw("p,E,same");
    legendScaled->AddEntry(graphCorrectedYieldWeightedAverageEtaPrimeStat,"Final","p");
    legendScaled->Draw();

    DrawGammaSetMarkerTF1( fitInvYieldEtaPrime, 7, 2, kGray+2);
    fitInvYieldEtaPrime->Draw("same");

    labelEnergyUnscaled->Draw();
    labelEtaPrimeUnscaled->Draw();
    labelDetProcUnscaled->Draw();

    canvasCorrScaled->Update();
    canvasCorrScaled->SaveAs(Form("%s/EtaPrime_%s_CorrectedYieldScaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    //***************************************************************************************************************
    //************************************Plotting final invariant yield ********************************************
    //***************************************************************************************************************

    histo2DInvYieldScaled->DrawCopy();
    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageEtaPrimeSys, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
    graphCorrectedYieldWeightedAverageEtaPrimeSys->Draw("p,E2,same");

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageEtaPrimeStat, 24, 2, kBlack , kBlack, 1, kTRUE);
    graphCorrectedYieldWeightedAverageEtaPrimeStat->Draw("p,E,same");

    fitInvYieldEtaPrime->Draw("same");

    labelEnergyUnscaled->Draw();
    labelEtaPrimeUnscaled->Draw();
    labelDetProcUnscaled->Draw();

    canvasCorrScaled->Update();
    canvasCorrScaled->SaveAs(Form("%s/EtaPrime_%s_CorrectedYieldFinal.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    if (fitBinShiftEtaPrimeTCM){
        DrawGammaSetMarkerTF1( fitBinShiftEtaPrimeTCM, 9, 2, kRed+2);
        fitBinShiftEtaPrimeTCM->Draw("same");
        canvasCorrScaled->Update();
        canvasCorrScaled->SaveAs(Form("%s/EtaPrime_%s_CorrectedYieldFinalWithBinShiftFit.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
    }


    histo2DInvYieldScaled->DrawCopy();
    graphCorrectedYieldWeightedAverageEtaPrimeSys->Draw("p,E2,same");
    graphCorrectedYieldWeightedAverageEtaPrimeStat->Draw("p,E,same");

    fitInvYieldEtaPrime->Draw("same");

    DrawGammaSetMarker(histoMCInputEtaPrime[0],  0, 0, kBlue+2, kBlue+2);
    histoMCInputEtaPrime[0]->Draw("same,hist,c");

    labelEnergyUnscaled->Draw();
    labelEtaPrimeUnscaled->Draw();
    labelDetProcUnscaled->Draw();

    canvasCorrScaled->Update();
    canvasCorrScaled->SaveAs(Form("%s/EtaPrime_%s_CorrectedYieldFinal_WithMC.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    delete canvasCorrScaled;

    //***************************************************************************************************************
    //****************************** Ratio to fit for individual spectra full range *********************************
    //***************************************************************************************************************
    TCanvas* canvasRatioSpec = new TCanvas("canvasRatioSpec","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasRatioSpec, 0.09, 0.017, 0.015, 0.08);
    canvasRatioSpec->SetLogy(0);

    TH2F * histo2DRatioToFitEtaPrime;
    histo2DRatioToFitEtaPrime = new TH2F("histo2DRatioToFitEtaPrime","histo2DRatioToFitEtaPrime",1000,0., maxPtGlobalEtaPrime,1000,0.55, 1.85);
    SetStyleHistoTH2ForGraphs(histo2DRatioToFitEtaPrime, "#it{p}_{T} (GeV/#it{c})","Data/Fit",
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
    histo2DRatioToFitEtaPrime->DrawCopy();

    TH1D* histoCorrectedYieldToFitEtaPrime[MaxNumberOfFiles];
    TGraphAsymmErrors* graphCorrectedYieldToFitEtaPrime[MaxNumberOfFiles];
    TLegend* legendRatioSpecEtaPrime = GetAndSetLegend2(0.12, 0.95-(1.05*nrOfTrigToBeCombEtaPrimeRed/2*0.85*textSizeSpectra), 0.5, 0.95,28);
    legendRatioSpecEtaPrime->SetNColumns(2);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if (graphsCorrectedYieldSysRemoved0EtaPrime[i] && !maskedFullyEtaPrime[i]){
          histoCorrectedYieldToFitEtaPrime[i] = CalculateHistoRatioToFit (histoCorrectedYieldEtaPrimeScaled[i], fitInvYieldEtaPrime);
          DrawGammaSetMarker(histoCorrectedYieldToFitEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
          legendRatioSpecEtaPrime->AddEntry(histoCorrectedYieldToFitEtaPrime[i],triggerNameLabel[i].Data(),"p");
          graphCorrectedYieldToFitEtaPrime[i] = CalculateGraphErrRatioToFit(graphsCorrectedYieldSysRemoved0EtaPrime[i], fitInvYieldEtaPrime);
          DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitEtaPrime[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
          graphCorrectedYieldToFitEtaPrime[i]->Draw("p,E2,same");
        }

        DrawGammaLines(0., maxPtGlobalEtaPrime , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0., maxPtGlobalEtaPrime , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0., maxPtGlobalEtaPrime , 0.9, 0.9,0.1, kGray, 7);

        if (graphsCorrectedYieldSysRemoved0EtaPrime[i]){
          histoCorrectedYieldToFitEtaPrime[i]->DrawCopy("e1,same");
        }
    }
    legendRatioSpecEtaPrime->Draw();

    TLatex *labelEnergyRatio = new TLatex(0.95, 0.93, collisionSystem.Data());
    SetStyleTLatex( labelEnergyRatio, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
    labelEnergyRatio->Draw();

    TLatex *labelEtaPrimeRatio = new TLatex(0.95, 0.93-textSizeSpectra*0.85*1.04, "#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelEtaPrimeRatio, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
    labelEtaPrimeRatio->Draw();

    TLatex *labelDetProcRatio = new TLatex(0.95, 0.93-(2*textSizeSpectra*0.85), detectionProcess.Data());
    SetStyleTLatex( labelDetProcRatio, 0.85*textSizeSpectra,4,1,42,kTRUE, 31);
    labelDetProcRatio->Draw();

    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/EtaPrime_%s_RatioSpectraToFit.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

    //***************************************************************************************************************
    //****************************** Ratio to fit for individual spectra used range *********************************
    //***************************************************************************************************************

    histo2DRatioToFitEtaPrime->DrawCopy();

    TGraphAsymmErrors* graphCorrectedYieldToFitEtaPrimeUsed[MaxNumberOfFiles];
    TGraphAsymmErrors* graphCorrectedYieldToFitEtaPrimeSysUsed[MaxNumberOfFiles];
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if (!maskedFullyEtaPrime[i]) {
            if (graphsCorrectedYieldShrunkEtaPrime[i]) graphCorrectedYieldToFitEtaPrimeUsed[i] = CalculateGraphErrRatioToFit (graphsCorrectedYieldShrunkEtaPrime[i], fitInvYieldEtaPrime);
            if (graphsCorrectedYieldSysShrunkEtaPrime[i]) graphCorrectedYieldToFitEtaPrimeSysUsed[i] = CalculateGraphErrRatioToFit(graphsCorrectedYieldSysShrunkEtaPrime[i], fitInvYieldEtaPrime);

            if (graphsCorrectedYieldSysShrunkEtaPrime[i])DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitEtaPrimeSysUsed[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
            if (graphsCorrectedYieldShrunkEtaPrime[i])DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitEtaPrimeUsed[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);

            if (graphsCorrectedYieldSysShrunkEtaPrime[i])graphCorrectedYieldToFitEtaPrimeSysUsed[i]->Draw("p,E2,same");

            DrawGammaLines(0., maxPtGlobalEtaPrime , 1., 1.,0.1, kGray+2);
            DrawGammaLines(0., maxPtGlobalEtaPrime , 1.1, 1.1,0.1, kGray, 7);
            DrawGammaLines(0., maxPtGlobalEtaPrime , 0.9, 0.9,0.1, kGray, 7);

            if (graphsCorrectedYieldShrunkEtaPrime[i])graphCorrectedYieldToFitEtaPrimeUsed[i]->Draw("e1,same");
        }
    }

    legendRatioSpecEtaPrime->Draw();

    labelEnergyRatio->Draw();
    labelEtaPrimeRatio->Draw();
    labelDetProcRatio->Draw();

    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/EtaPrime_%s_RatioSpectraToFitUsed.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

    //***************************************************************************************************************
    //****************************** Ratio to fit for final spectrum ************************************************
    //***************************************************************************************************************

    histo2DRatioToFitEtaPrime->DrawCopy();

    TGraphAsymmErrors* graphCorrectedYieldFinalStatToFitEtaPrime;
    TGraphAsymmErrors* graphCorrectedYieldFinalSysToFitEtaPrime;
    TH1D* histoMCInputToFit;

    graphCorrectedYieldFinalStatToFitEtaPrime    = CalculateGraphErrRatioToFit (graphCorrectedYieldWeightedAverageEtaPrimeStat, fitInvYieldEtaPrime);
    graphCorrectedYieldFinalSysToFitEtaPrime     = CalculateGraphErrRatioToFit(graphCorrectedYieldWeightedAverageEtaPrimeSys, fitInvYieldEtaPrime);
    histoMCInputToFit                       = (TH1D*)histoMCInputEtaPrime[0]->Clone("EtaPrimeMCToFit");
    histoMCInputToFit                       = CalculateHistoRatioToFitNLO (histoMCInputToFit, fitInvYieldEtaPrime, minPtGlobalEtaPrime);

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldFinalSysToFitEtaPrime, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldFinalStatToFitEtaPrime, 24, 2, kBlack , kBlack, 1, kTRUE);

    graphCorrectedYieldFinalSysToFitEtaPrime->Draw("p,E2,same");

    DrawGammaLines(0., maxPtGlobalEtaPrime , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0., maxPtGlobalEtaPrime , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(0., maxPtGlobalEtaPrime , 0.9, 0.9,0.1, kGray, 7);

    graphCorrectedYieldFinalStatToFitEtaPrime->Draw("p,E,same");

    labelEnergyRatio->Draw();
    labelEtaPrimeRatio->Draw();
    labelDetProcRatio->Draw();

    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/EtaPrime_%s_RatioSpectraToFitFinal.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

    TH2F * histo2DRatioToFitEtaPrime2;
    histo2DRatioToFitEtaPrime2 = new TH2F("histo2DRatioToFitEtaPrime2","histo2DRatioToFitEtaPrime2",1000,0., maxPtGlobalEtaPrime,1000,0.25, 1.85);
    SetStyleHistoTH2ForGraphs(histo2DRatioToFitEtaPrime2, "#it{p}_{T} (GeV/#it{c})","Data/Fit",
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
    histo2DRatioToFitEtaPrime2->DrawCopy();

    graphCorrectedYieldFinalSysToFitEtaPrime->Draw("p,E2,same");

    DrawGammaLines(0., maxPtGlobalEtaPrime , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0., maxPtGlobalEtaPrime , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(0., maxPtGlobalEtaPrime , 0.9, 0.9,0.1, kGray, 7);


    graphCorrectedYieldFinalStatToFitEtaPrime->Draw("p,E,same");
    DrawGammaSetMarker(histoMCInputToFit, 20, 1.2, kBlue+1, kBlue+1);
    histoMCInputToFit->Draw("same,pe");
    DrawGammaLines(0., maxPtGlobalEtaPrime , 0.6, 0.6,1, kBlue-7, 7);

    labelEnergyRatio->Draw();
    labelEtaPrimeRatio->Draw();
    labelDetProcRatio->Draw();

    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/EtaPrime_%s_RatioSpectraToFitFinal_withMC.%s",outputDir.Data(),isMC.Data(), suffix.Data()));


    delete canvasRatioSpec;

    //*********************************************************************************************************
    //********************** ComparisonFile Output ************************************************************
    //*********************************************************************************************************
    cout << "starting to write out data" << endl;
    TGraphAsymmErrors* graphInvXSectionWeightedAverageEtaPrimeStat   = NULL;
    TH1D* histoInvXSectionWeightedAverageEtaPrimeStat                = NULL;
    TH1D* histoInvYieldWeightedAverageEtaPrimeStat                   = NULL;
    TGraphAsymmErrors* graphInvXSectionWeightedAverageEtaPrimeSys    = NULL;
    cout << "xSection = " << xSection << endl;

    if (graphCorrectedYieldWeightedAverageEtaPrimeStat){
        histoInvYieldWeightedAverageEtaPrimeStat                  = new TH1D("histoInvYieldWeightedAverageEtaPrimeStat", "", maxNAllowedEtaPrime, binningEtaPrime);
        Int_t firstBinEtaPrime = 1;
        while (histoInvYieldWeightedAverageEtaPrimeStat->GetBinCenter(firstBinEtaPrime) < graphCorrectedYieldWeightedAverageEtaPrimeStat->GetX()[0]){
            histoInvYieldWeightedAverageEtaPrimeStat->SetBinContent(firstBinEtaPrime, 0);
            histoInvYieldWeightedAverageEtaPrimeStat->SetBinError(firstBinEtaPrime, 0);
            firstBinEtaPrime++;
        }
        for (Int_t i = 0; i < graphCorrectedYieldWeightedAverageEtaPrimeStat->GetN(); i++){
            histoInvYieldWeightedAverageEtaPrimeStat->SetBinContent(i+firstBinEtaPrime, graphCorrectedYieldWeightedAverageEtaPrimeStat->GetY()[i]);
            histoInvYieldWeightedAverageEtaPrimeStat->SetBinError(i+firstBinEtaPrime, graphCorrectedYieldWeightedAverageEtaPrimeStat->GetEYlow()[i]);
        }
        if (xSection != 0){
            graphInvXSectionWeightedAverageEtaPrimeStat                  = ScaleGraph(graphCorrectedYieldWeightedAverageEtaPrimeStat,xSection*recalcBarn);
            histoInvXSectionWeightedAverageEtaPrimeStat                  = new TH1D("histoInvXSectionWeightedAverageEtaPrimeStat", "", maxNAllowedEtaPrime, binningEtaPrime);
            for (Int_t i = 0; i < graphInvXSectionWeightedAverageEtaPrimeStat->GetN(); i++){
                histoInvXSectionWeightedAverageEtaPrimeStat->SetBinContent(i+firstBinEtaPrime, graphInvXSectionWeightedAverageEtaPrimeStat->GetY()[i]);
                histoInvXSectionWeightedAverageEtaPrimeStat->SetBinError(i+firstBinEtaPrime, graphInvXSectionWeightedAverageEtaPrimeStat->GetEYlow()[i]);
            }
            graphInvXSectionWeightedAverageEtaPrimeStat->Print();
            if (graphCorrectedYieldWeightedAverageEtaPrimeSys)
                graphInvXSectionWeightedAverageEtaPrimeSys                   = ScaleGraph(graphCorrectedYieldWeightedAverageEtaPrimeSys,xSection*recalcBarn);
            graphInvXSectionWeightedAverageEtaPrimeSys->Print();
        }
    }

    cout << "Writing output file" << endl;
    TString fileNameOutputComp  = Form("%s/%s_%sResultsFullCorrection_PP.root",outputDir.Data(),isMC.Data(),system.Data());
    if (optionEnergy.Contains("pPb"))
        fileNameOutputComp      = Form("%s/%s_%sResultsFullCorrection_pPb.root",outputDirDay.Data(),isMC.Data(),system.Data());
    else if (optionEnergy.Contains("PbPb"))
        fileNameOutputComp      = Form("%s/%s_%sResultsFullCorrection_PbPb.root",outputDirDay.Data(),isMC.Data(),system.Data());
    else if (optionEnergy.Contains("XeXe"))
        fileNameOutputComp      = Form("%s/%s_%sResultsFullCorrection_XeXe.root",outputDirDay.Data(),isMC.Data(),system.Data());

    TFile* fileOutputForComparisonFullyCorrected = new TFile(fileNameOutputComp,"UPDATE");
        for (Int_t i=0; i< nrOfTrigToBeComb; i++){
            if (histoEventQualtity[i])                  histoEventQualtity[i]->Write(Form("NEvents_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (modeNormal == 2 || modeNormal == 3 || modeNormal == 4 || modeNormal == 5 ){
                if (histoTriggerRejection[i])           histoTriggerRejection[i]->Write(Form("TriggRejectMean_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),TObject::kOverwrite);
                if (histoRatioRawClusterPt[i])          histoRatioRawClusterPt[i]->Write(Form("TriggRejectvsPt_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),TObject::kOverwrite);
                if (histoRatioRawClusterE[i])           histoRatioRawClusterE[i]->Write(Form("TriggRejectvsE_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),TObject::kOverwrite);
            }
        }
        fileOutputForComparisonFullyCorrected->mkdir(Form("EtaPrime%s%s",fCent.Data(),optionEnergy.Data()));
        TDirectoryFile* directoryEtaPrime = (TDirectoryFile*)fileOutputForComparisonFullyCorrected->Get(Form("EtaPrime%s%s",fCent.Data(),optionEnergy.Data()));
        fileOutputForComparisonFullyCorrected->cd(Form("EtaPrime%s%s",fCent.Data(),optionEnergy.Data()));
        if (graphCorrectedYieldWeightedAverageEtaPrimeStat)  graphCorrectedYieldWeightedAverageEtaPrimeStat->Write("graphCorrectedYieldEtaPrime",TObject::kOverwrite);
        if (histoInvYieldWeightedAverageEtaPrimeStat)  histoInvYieldWeightedAverageEtaPrimeStat->Write("CorrectedYieldEtaPrime",TObject::kOverwrite);
        if (graphCorrectedYieldWeightedAverageEtaPrimeSys)   graphCorrectedYieldWeightedAverageEtaPrimeSys->Write("EtaPrimeSystError",TObject::kOverwrite);
        if (xSection != 0){
            if (graphInvXSectionWeightedAverageEtaPrimeStat){
                graphInvXSectionWeightedAverageEtaPrimeStat->Write("graphInvCrossSectionEtaPrime",TObject::kOverwrite);
                cout << "graph InvSection etaprime stat" << endl;
                graphInvXSectionWeightedAverageEtaPrimeStat->Print();
            }
            if (histoInvXSectionWeightedAverageEtaPrimeStat)     histoInvXSectionWeightedAverageEtaPrimeStat->Write("InvCrossSectionEtaPrime",TObject::kOverwrite);
            if (graphInvXSectionWeightedAverageEtaPrimeSys){
                graphInvXSectionWeightedAverageEtaPrimeSys->Write("InvCrossSectionEtaPrimeSys",TObject::kOverwrite);
                cout << "graph InvSection etaprime sys" << endl;
                graphInvXSectionWeightedAverageEtaPrimeSys->Print();
            }
        }
        if (graphAcceptanceEtaPrimeWeighted)                 graphAcceptanceEtaPrimeWeighted->Write("AcceptanceEtaPrime", TObject::kOverwrite);
        if (graphEfficiencyEtaPrimeWeighted)                 graphEfficiencyEtaPrimeWeighted->Write("EfficiencyEtaPrime", TObject::kOverwrite);
        if (graphEffTimesAccEtaPrimeWeighted)                graphEffTimesAccEtaPrimeWeighted->Write("EffTimesAccEtaPrime", TObject::kOverwrite);
        if (graphMassEtaPrimeDataWeighted)                   graphMassEtaPrimeDataWeighted->Write("EtaPrime_Mass_data", TObject::kOverwrite);
        if (graphMassEtaPrimeMCWeighted)                     graphMassEtaPrimeMCWeighted->Write("EtaPrime_Mass_MC", TObject::kOverwrite);
        if (graphWidthEtaPrimeDataWeighted)                  graphWidthEtaPrimeDataWeighted->Write("EtaPrime_Width_data", TObject::kOverwrite);
        if (graphWidthEtaPrimeMCWeighted)                    graphWidthEtaPrimeMCWeighted->Write("EtaPrime_Width_MC", TObject::kOverwrite);
        if (sysAvailSingleEtaPrime[0]){
            for (Int_t k = 0; k< nRelSysErrEtaPrimeSources ; k++ ){
                if (graphRelSysErrEtaPrimeSourceWeighted[k]) graphRelSysErrEtaPrimeSourceWeighted[k]->Write(graphRelSysErrEtaPrimeSourceWeighted[k]->GetName(), TObject::kOverwrite);
            }
        }
        for (Int_t i=0; i< nrOfTrigToBeComb; i++){
            cout << "trigger: " << triggerName[i].Data() << endl;
            if (histoAcceptanceEtaPrime[i])                  histoAcceptanceEtaPrime[i]->Write(Form("AcceptanceEtaPrime_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoEfficiencyEtaPrime[i])                  histoEfficiencyEtaPrime[i]->Write(Form("EfficiencyEtaPrime_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoRawYieldEtaPrime[i])                    histoRawYieldEtaPrime[i]->Write(Form("RAWYieldPerEventsEtaPrime_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoMCInputEtaPrime[i])                     histoMCInputEtaPrime[i]->Write(Form("EtaPrime_Input_Reweighted_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoEffTimesAccEtaPrime[i])                 histoEffTimesAccEtaPrime[i]->Write(Form("EffTimesAccEtaPrime_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoMassEtaPrimeData[i])                    histoMassEtaPrimeData[i]->Write(Form("EtaPrime_Mass_data_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoMassEtaPrimeMC[i])                      histoMassEtaPrimeMC[i]->Write(Form("EtaPrime_Mass_MC_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoWidthEtaPrimeData[i])                   histoWidthEtaPrimeData[i]->Write(Form("EtaPrime_Width_data_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoWidthEtaPrimeMC[i])                     histoWidthEtaPrimeMC[i]->Write(Form("EtaPrime_Width_MC_%s",triggerName[i].Data()),TObject::kOverwrite);

            if (histoInvMassSig[i])                     histoInvMassSig[i]->Write(Form("EtaPrime_InvMassSig_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoInvMassSigPlusBG[i])               histoInvMassSigPlusBG[i]->Write(Form("EtaPrime_InvMassSigPlusBG_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoInvMassBG[i])                      histoInvMassBG[i]->Write(Form("EtaPrime_InvMassBG_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (fitInvMassSig[i])                       fitInvMassSig[i]->Write(Form("EtaPrime_InvMassSigFit_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (graphsCorrectedYieldRemoved0EtaPrime[i])     graphsCorrectedYieldRemoved0EtaPrime[i]->Write(Form("CorrectedYieldEtaPrime_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (graphsCorrectedYieldSysRemoved0EtaPrime[i])  graphsCorrectedYieldSysRemoved0EtaPrime[i]->Write(Form("EtaPrimeSystError_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (enableTriggerEffEtaPrime[i]){
                if (histoTriggerEffEtaPrime[i])              histoTriggerEffEtaPrime[i]->Write(Form("TriggerEfficiencyEtaPrime_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoEffBaseEtaPrime[i])                 histoEffBaseEtaPrime[i]->Write(Form("EfficiencyBaseEtaPrime_%s",triggerName[i].Data()),TObject::kOverwrite);
            }
        }

    fileOutputForComparisonFullyCorrected->Write();
    fileOutputForComparisonFullyCorrected->Close();
}



//****************************************************************************
//****** Main function for non-ROOT compilation ******************************
//****************************************************************************
/*
    Tip: add the following functions to your bashrc.sh:
        function compile () {
            g++ "$1" -I$(root-config --incdir) $(root-config --libs --evelibs --glibs) -lMinuit -o "${1%.*}"
        }
        function alicompile () {
            g++ "$1" -I$ALIEN_INCLUDE_DIR -L$ALIEN_LIBRARY -o "${1%.*}"
        }
    You can add the flag "-fsanitize=address -g" to debug and/or search for memory leaks.

    Then compile and run the macro using:
        compile <yourcode.C>
        ./yourcode arg1 arg2 arg3 ...
*/
int main( int argc, char* argv[] )
{
    // Default arguments for ExtractSignal
        TString fileListNameEtaPrime    = "triggerFileListEtaPrime.txt";
        Int_t   mode                    = 4;
        Int_t   numberOfTrigg           = 6;
        TString suffix                  = "eps";
        TString isMC                    = "";
        TString optionEnergy            = "";
        TString period                  = "";
        Bool_t  pileUpApplied           = kTRUE;
        Float_t maxPtGlobalEtaPrime     = 16.;
        Bool_t  averagedEtaPrime        = kFALSE;
        TString nameFileFitsShift       = "";
        Bool_t  hasClusterOutput        = kTRUE;
        TString fileInputCorrFactors    = "";

    // Import main call arguments
        TString import;
        if(argc >  1) fileListNameEtaPrime = argv[1];
        if(argc >  2) { // mode
            istringstream sstr(argv[2]);
            sstr >> mode;
        }
        if(argc >  3) { // numberOfTrigg
            istringstream sstr(argv[3]);
            sstr >> numberOfTrigg;
        }
        if(argc >  4) suffix       = argv[4];
        if(argc >  5) isMC         = argv[5];
        if(argc >  6) optionEnergy = argv[6];
        if(argc >  7) period       = argv[7];
        if(argc >  8) { // pileUpApplied
            import = argv[8];
            if( import.EqualTo("kTRUE")  || import.EqualTo("1") ) pileUpApplied = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("0") ) pileUpApplied = kFALSE;
        }
        if(argc >  9) { // maxPtGlobalEtaPrime
            istringstream sstr(argv[9]);
            sstr >> maxPtGlobalEtaPrime;
        }
        if(argc > 10) { // averagedEtaPrime
            import = argv[10];
            if( import.EqualTo("kTRUE")  || import.EqualTo("1") ) averagedEtaPrime = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("0") ) averagedEtaPrime = kFALSE;
        }
        if(argc > 11) nameFileFitsShift = argv[11];
        if(argc > 12) { // hasClusterOutput
            import = argv[12];
            if( import.EqualTo("kTRUE")  || import.EqualTo("1") ) hasClusterOutput = kTRUE;
            if( import.EqualTo("kFALSE") || import.EqualTo("0") ) hasClusterOutput = kFALSE;
        }
        if(argc > 13) fileInputCorrFactors = argv[13];

    // Function call ExtractSignalV2
        cout << "Executing \"ProduceFinalResultsPatchedTriggersEtaPrime\" with the following arguments:" << endl
            << "  fileListNameEtaPrime: \"" << fileListNameEtaPrime << "\"" << endl
            << "  mode:                 \"" << mode                 << "\"" << endl
            << "  numberOfTrigg:        \"" << numberOfTrigg        << "\"" << endl
            << "  suffix:               \"" << suffix               << "\"" << endl
            << "  isMC:                 \"" << isMC                 << "\"" << endl
            << "  optionEnergy:         \"" << optionEnergy         << "\"" << endl
            << "  period:               \"" << period               << "\"" << endl
            << "  pileUpApplied:        \"" << pileUpApplied        << "\"" << endl
            << "  maxPtGlobalEtaPrime:  \"" << maxPtGlobalEtaPrime  << "\"" << endl
            << "  averagedEtaPrime:     \"" << averagedEtaPrime     << "\"" << endl
            << "  nameFileFitsShift:    \"" << nameFileFitsShift    << "\"" << endl
            << "  hasClusterOutput:     \"" << hasClusterOutput     << "\"" << endl
            << "  fileInputCorrFactors: \"" << fileInputCorrFactors << "\"" << endl
        << endl;
        ProduceFinalResultsPatchedTriggersEtaPrime(
            fileListNameEtaPrime,
            mode,
            numberOfTrigg,
            suffix,
            isMC,
            optionEnergy,
            period,
            pileUpApplied,
            maxPtGlobalEtaPrime,
            averagedEtaPrime,
            nameFileFitsShift,
            hasClusterOutput,
            fileInputCorrFactors
        );

    // Main function return
    return 0;
}