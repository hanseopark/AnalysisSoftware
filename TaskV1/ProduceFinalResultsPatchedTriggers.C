/****************************************************************************************************************************
 ******     provided by Gamma Conversion Group, PWGGA,                                                                  *****
 ******        Friederike Bock, friederike.bock@cern.ch                                                                 *****
 *****************************************************************************************************************************/

#include <Riostream.h>
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
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

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


void  ProduceFinalResultsPatchedTriggers(   TString fileListNamePi0     = "triggerFileListPi0.txt",
                                            Int_t   mode                = 4,
                                            Int_t   numberOfTrigg       = 6,
                                            TString suffix              = "eps", 
                                            TString isMC                = "",
                                            TString optionEnergy        = "",
                                            TString period              = "", 
                                            Bool_t  pileUpApplied       = kTRUE,
                                            Float_t maxPtGlobalPi0      = 16.,
                                            Bool_t  averagedPi0         = kFALSE,
                                            Bool_t  enableEta           = kTRUE,
                                            Float_t maxPtGlobalEta      = 14.,
                                            Bool_t  averagedEta         = kFALSE,
                                            Bool_t  v2ClusterizerMerged = kFALSE,
                                            TString nameFileFitsShift   = "",
                                            Bool_t  hasClusterOutput    = kTRUE,
                                            TString fileInputCorrFactors= ""
                                        ){

    //***************************************************************************************************************
    //************************************ Layouting preparations & general setup ***********************************
    //***************************************************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();
    SetPlotStyle();
    
    TString dateForOutput           = ReturnDateStringForOutput();
    TString collisionSystem         = ReturnFullCollisionsSystem(optionEnergy);
    TString detectionProcess        = ReturnFullTextReconstructionProcess(mode);
    TString detectionProcessClus    = ReturnFullTextReconstructionProcess(mode,2);
    
    if (isMC.CompareTo("MC") == 0) collisionSystem = "MC, "+collisionSystem;
    
    Double_t binningPi0[100];
    Int_t maxNBinsPi0               = GetBinning( binningPi0, "Pi0", optionEnergy, mode );
    Int_t maxNAllowedPi0            = 0;
    Int_t nRealTriggers             = 0;
    while (binningPi0[maxNAllowedPi0] < maxPtGlobalPi0 ) maxNAllowedPi0++;
    for (Int_t i= 0; i< maxNAllowedPi0+1; i++){
        cout << binningPi0[i] << endl;
    }
    
    Double_t binningEta[100];
    Int_t maxNBinsEta               = GetBinning( binningEta, "Eta", optionEnergy, mode );
    Int_t maxNAllowedEta            = 0;
    while (binningEta[maxNAllowedEta] < maxPtGlobalEta ) maxNAllowedEta++;
    for (Int_t i= 0; i< maxNAllowedEta+1; i++){
        cout << binningEta[i] << endl;
    }
    
    Double_t maxPtGlobalCluster     = 25;
    if (optionEnergy.CompareTo("2.76TeV")==0){
        if (mode==2){
            maxPtGlobalCluster          = 25;
        } else if (mode == 4){
            maxPtGlobalCluster          = 30;
        } else if (mode == 10){
            maxPtGlobalCluster          = 50;
        }
    } else if (optionEnergy.CompareTo("8TeV")==0){
      if(mode==2 || mode==4){
        maxPtGlobalCluster          = 50;
      } else if (mode == 10){
        maxPtGlobalCluster          = 70;
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
    
    //***************************************************************************************************************
    //********************************** setting correct histo names ************************************************
    //***************************************************************************************************************
    TString nameCorrectedYield                          = "CorrectedYieldTrueEff";
    TString nameEfficiency                              = "TrueMesonEffiPt";
    TString nameSecFK0sEfficiency                       = "TrueSecFromK0SEffiPt";
    TString nameAcceptance                              = "fMCMesonAccepPt";
    TString nameAcceptanceWOEvtWeights                  = "fMCMesonAccepPtWOEvtWeights";
    TString nameMassMC                                  = "histoTrueMassMeson";
    TString nameWidthMC                                 = "histoTrueFWHMMeson";
    TString nameMCYield                                 = "MCYield_Meson_oldBinWOWeights";
    if ( mode == 4 || mode == 5 ){
        nameCorrectedYield                              = "CorrectedYieldNormEff";
        nameEfficiency                                  = "MesonEffiPt";
        nameMassMC                                      = "histoMassMesonRecMC";
        nameWidthMC                                     = "histoFWHMMesonRecMC";
    } else if (mode == 2 || mode == 3){
        nameMassMC                                      = "histoMassMesonRecMC";
        nameWidthMC                                     = "histoFWHMMesonRecMC";        
    } else if (mode == 10){
        nameCorrectedYield                              = "CorrectedYieldTrueEff";
        nameEfficiency                                  = "PrimaryMesonEfficiency";
        nameAcceptance                                  = "fHistoMCAcceptancePt";
        nameMCYield                                     = "MCYield_Meson_oldBin";
        nameSecFK0sEfficiency                           = "TrueMesonEffiSecFromK0SPt";
    }
    
    TString cutNumber       [MaxNumberOfFiles];
    TString triggerName     [MaxNumberOfFiles];
    TString triggerNameLabel[MaxNumberOfFiles];
    Float_t minPt           [MaxNumberOfFiles];
    Float_t maxPt           [MaxNumberOfFiles];
    Int_t trigSteps         [MaxNumberOfFiles][3];
    Float_t ptFromSpecPi0   [MaxNumberOfFiles][2];
    Float_t ptFromSpecEta   [MaxNumberOfFiles][2];
    Bool_t maskedFullyPi0   [MaxNumberOfFiles];
    Bool_t maskedFullyEta   [MaxNumberOfFiles];
    TString sysFilePi0      [MaxNumberOfFiles];
    TString sysFileEta      [MaxNumberOfFiles];
    TString sysFileEtaToPi0 [MaxNumberOfFiles];
    TString cutNumberBaseEff[MaxNumberOfFiles];
    //***************************************************************************************************************
    //*************************** read setting from configuration file **********************************************
    //***************************************************************************************************************
    ifstream in(fileListNamePi0.Data());
    cout<<"Available Triggers:"<<endl;
    // general number of triggers set
    Int_t nrOfTrigToBeComb      = 0;
    // number of triggers which are really used for the respective analysis
    Int_t nrOfTrigToBeCombPi0Red        = 0;
    Int_t nrOfTrigToBeCombEtaRed        = 0;
    Int_t nrOfTrigToBeCombEtaToPi0Red   = 0;
    while(!in.eof() && nrOfTrigToBeComb<numberOfTrigg ){
        in >> cutNumber[nrOfTrigToBeComb] >> minPt[nrOfTrigToBeComb] >> maxPt[nrOfTrigToBeComb] >> triggerName[nrOfTrigToBeComb] 
        >> trigSteps[nrOfTrigToBeComb][0]  >> trigSteps[nrOfTrigToBeComb][1]  >> trigSteps[nrOfTrigToBeComb][2] >> ptFromSpecPi0[nrOfTrigToBeComb][0] >> ptFromSpecPi0[nrOfTrigToBeComb][1] 
        >> ptFromSpecEta[nrOfTrigToBeComb][0] >> ptFromSpecEta[nrOfTrigToBeComb][1] >> sysFilePi0[nrOfTrigToBeComb] >> sysFileEta[nrOfTrigToBeComb] >> sysFileEtaToPi0[nrOfTrigToBeComb] >>
        cutNumberBaseEff[nrOfTrigToBeComb];
        cout<< cutNumber[nrOfTrigToBeComb]<< "\t"<< triggerName[nrOfTrigToBeComb] << "\t transverse momentum range: " << minPt[nrOfTrigToBeComb]<< "\t to "<< maxPt[nrOfTrigToBeComb] <<endl;
        cout << trigSteps[nrOfTrigToBeComb][0] << "\t" << trigSteps[nrOfTrigToBeComb][1] << "\t"<< trigSteps[nrOfTrigToBeComb][2] << endl;
        nrOfTrigToBeComb++;
    }

    for (Int_t i = 0; i < nrOfTrigToBeComb; i++){
        // figure out which triggers are fully masked for the pi0
        if ((ptFromSpecPi0[i][1] == -1 && ptFromSpecPi0[i][0] == -1 )|| (ptFromSpecPi0[i][0] == 0 && ptFromSpecPi0[i][1] == 0)){
            maskedFullyPi0[i]       = kTRUE;
        } else {
            maskedFullyPi0[i]       = kFALSE;
        }
        // figure out which triggers are fully masked for the eta
        if ((ptFromSpecEta[i][1] == -1 && ptFromSpecEta[i][0] == -1 )|| (ptFromSpecEta[i][0] == 0 && ptFromSpecEta[i][1] == 0)){
            maskedFullyEta[i]       = kTRUE;
        } else {
            maskedFullyEta[i]       = kFALSE;
        }
    }
    
    Double_t minPtGlobalPi0         = ptFromSpecPi0[0][0];
    Double_t minPtGlobalEta         = ptFromSpecEta[0][0];
    if (maskedFullyPi0[0]) 
        minPtGlobalPi0              = 12;
    if (maskedFullyEta[0]) 
        minPtGlobalEta              = 12;
    
    
    for (Int_t j = 1; j < nrOfTrigToBeComb; j++){
        if (minPtGlobalPi0 > ptFromSpecPi0[j][0] && !maskedFullyPi0[j] ) 
            minPtGlobalPi0 = ptFromSpecPi0[j][0];
        if (minPtGlobalEta > ptFromSpecEta[j][0] && !maskedFullyEta[j] ) 
            minPtGlobalEta = ptFromSpecEta[j][0];
    }
    cout << "global minimum pT for pi0: " << minPtGlobalPi0 << endl;
    cout << "global minimum pT for eta: " << minPtGlobalEta << endl;
    
    // set individual triggers to total number of availabel triggers, will be reduced later according to usage
    nrOfTrigToBeCombPi0Red      = nrOfTrigToBeComb;
    nrOfTrigToBeCombEtaRed      = nrOfTrigToBeComb;
    nrOfTrigToBeCombEtaToPi0Red = nrOfTrigToBeComb;
    
    // variables to keep track of NLM 
    TString fNLMString      = "";
    TString fNLMStringOutput= "";
    TString fMergedClusterCutNrExampl   = "";
    Int_t fNLMmin           = 0;
    
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

        // put correct labels if only 1 set of NLM is used for merged ana
        if (mode == 10){
            TString fEventCutSelection                  = "";
            TString fClusterCutSelection                = "";
            TString fClusterMergedCutSelection          = "";
            TString dummyString                         = "";
            TString fMesonCutSelection                  = "";
            ReturnSeparatedCutNumberAdvanced( cutNumber[i].Data(), fEventCutSelection, fClusterCutSelection, fClusterMergedCutSelection, dummyString, fMesonCutSelection, mode);
            if (i == 0){
                fMergedClusterCutNrExampl = fClusterMergedCutSelection;
                fNLMmin = ReturnClusterNLM(fClusterMergedCutSelection);
                if (fNLMmin){
                    fNLMString                          = Form("%i local maximum", fNLMmin);        
                    fNLMStringOutput                    = Form("LM%i", fNLMmin);
                } else {
                    fNLMString                          = Form("%i local maxima", fNLMmin);
                    fNLMStringOutput                    = Form("LM%i", fNLMmin);
                }    
            } else {
                if ( fNLMmin != ReturnClusterNLM(fClusterMergedCutSelection)){
                    fMergedClusterCutNrExampl           = "";
                    fNLMString                          = "";
                    fNLMStringOutput                    = "";
                    fNLMmin                             = 0;
                }    
            }
        }    
    }
    
    // defining output directory
    TString outputDir =    Form("%s/%s/%s/FinalResultsTriggersPatched%s", suffix.Data(),optionEnergy.Data(),dateForOutput.Data(),fNLMStringOutput.Data());
    if(optionEnergy.CompareTo("8TeV") == 0){
      if(mode == 4) outputDir = Form("%s/%s/%s/FinalResultsTriggersPatched_EMCAL%s", suffix.Data(),optionEnergy.Data(),dateForOutput.Data(),fNLMStringOutput.Data());
      if(mode == 2) outputDir = Form("%s/%s/%s/FinalResultsTriggersPatched_PCMEMCAL%s", suffix.Data(),optionEnergy.Data(),dateForOutput.Data(),fNLMStringOutput.Data());
    }
    gSystem->Exec("mkdir -p "+outputDir);
   
    gSystem->Exec(Form("cp %s %s/configurationFile.txt", fileListNamePi0.Data(), outputDir.Data()));
    
    // put correct labels if only 1 set of NLM is used for merged ana
    if (fNLMmin > 0 && mode == 10){
        detectionProcess    = ReturnFullTextReconstructionProcess(mode, 0, "#pi^{0}", fMergedClusterCutNrExampl);
    }    
    
    vector<TString>** ptSysDetail     = new vector<TString>*[MaxNumberOfFiles];
    for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++) ptSysDetail[iR] = new vector<TString>[50];
    TString sysStringComb = "PCM";
    if(mode == 2) sysStringComb = "PCMEMC";
    else if(mode == 3) sysStringComb = "PCMPHOS";
    else if(mode == 4) sysStringComb = "EMCEMC";
    else if(mode == 5) sysStringComb = "PHOS";
    else if(mode == 10) sysStringComb = "EMCm";

    //***************************************************************************************************************
    //******************************** Load Pi0 histograms **********************************************************
    //***************************************************************************************************************
    TString FileNameCorrectedPi0    [MaxNumberOfFiles];
    TFile*  fileCorrectedPi0        [MaxNumberOfFiles];
    TString FileNameUnCorrectedPi0  [MaxNumberOfFiles];
    TFile*  fileUnCorrectedPi0      [MaxNumberOfFiles];
    TString FileNameUnCorrectedMCPi0[MaxNumberOfFiles];
    TFile*  fileUnCorrectedMCPi0    [MaxNumberOfFiles];

    TH1D*   histoCorrectedYieldPi0  [MaxNumberOfFiles];
    TH1D*   histoRawClusterPt       [MaxNumberOfFiles];
    TH1D*   histoEfficiencyPi0      [MaxNumberOfFiles];
    TH1D*   histoAcceptancePi0      [MaxNumberOfFiles];
    
    TH1D*   histoAcceptancePi0WOEvtWeights [MaxNumberOfFiles];
    TH1D*   histoPurityPi0          [MaxNumberOfFiles];
    TH1D*   histoEffTimesAccPi0     [MaxNumberOfFiles];
    TH1D*   histoRawYieldPi0        [MaxNumberOfFiles];
    TH1D*   histoRatioRawClusterPt  [MaxNumberOfFiles];
    TH1D*   histoTriggerRejection   [MaxNumberOfFiles];
    TH1D*   histoMassPi0Data        [MaxNumberOfFiles];
    TH1D*   histoMassPi0MC          [MaxNumberOfFiles];
    TH1D*   histoWidthPi0Data       [MaxNumberOfFiles];
    TH1D*   histoWidthPi0MC         [MaxNumberOfFiles];
    TH1F*   histoEventQualtity      [MaxNumberOfFiles];
    TH1D*   histoMCInputPi0         [MaxNumberOfFiles];
    
    Double_t triggRejecFac          [MaxNumberOfFiles][MaxNumberOfFiles];
    Double_t triggRejecFacErr       [MaxNumberOfFiles][MaxNumberOfFiles];
    
    TString FileNameEffBasePi0      [MaxNumberOfFiles];
    TFile*  fileEffBasePi0          [MaxNumberOfFiles];
    TH1D*   histoEffBasePi0         [MaxNumberOfFiles];
    TH1D*   histoTriggerEffPi0      [MaxNumberOfFiles];
    Bool_t  enableTriggerEffPi0     [MaxNumberOfFiles];
    Bool_t  enableTriggerEffPi0All                          = kFALSE;
    Bool_t  enableTriggerRejecCompMC                        = kFALSE;
    
    TH1D*   histoMCRawClusterPt     [MaxNumberOfFiles];
    TH1D*   histoMCRatioRawClusterPt[MaxNumberOfFiles];
    TString rapidityRange                                   = "";
    Double_t deltaRapid             [MaxNumberOfFiles];
    
    TH1D*   histoInvMassSigPlusBG   [MaxNumberOfFiles];
    TH1D*   histoInvMassSig         [MaxNumberOfFiles];
    TH1D*   histoInvMassBG          [MaxNumberOfFiles];
    TF1*    fitInvMassSig           [MaxNumberOfFiles];

    // sec corr histos
    TH1D*   histoSecEffiPi0FromK0s  [MaxNumberOfFiles];
    TH1D*   histoEffectCorrPi0FromK0s[MaxNumberOfFiles];
    Bool_t hasSecEffi                                       = kFALSE;
    Bool_t hasSecCorrFac                                    = kFALSE;
    
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
        Int_t exampleBin                                    = ReturnSingleInvariantMassBinPlotting ("Pi0", optionEnergy, mode, trigger.Atoi());
        
        FileNameCorrectedPi0[i]                             = Form("%s/%s/Pi0_%s_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(), 
                                                                    cutNumber[i].Data());
        if (mode == 10) FileNameCorrectedPi0[i]             = Form("%s/%s/Pi0_%s_GammaMergedCorrection_%s.root", cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(), 
                                                                    cutNumber[i].Data());
        cout<< FileNameCorrectedPi0[i] << endl;
        fileCorrectedPi0[i]                                 = new TFile(FileNameCorrectedPi0[i]);
        if (fileCorrectedPi0[i]->IsZombie()) return;
        // read uncorrected file
        FileNameUnCorrectedPi0[i]                           = Form("%s/%s/Pi0_%s_GammaConvV1WithoutCorrection_%s.root",cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(), 
                                                                    cutNumber[i].Data());
        if (mode == 10) FileNameUnCorrectedPi0[i]           = Form("%s/%s/Pi0_%s_GammaMergedWithoutCorrection_%s.root",cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(), 
                                                                    cutNumber[i].Data());

        cout<< FileNameUnCorrectedPi0[i] << endl;
        fileUnCorrectedPi0[i]                               = new TFile(FileNameUnCorrectedPi0[i]);
        if (fileUnCorrectedPi0[i]->IsZombie()) return;

        if (isMC.CompareTo("data") == 0){
            FileNameUnCorrectedMCPi0[i]                     = Form("%s/%s/Pi0_MC_GammaConvV1WithoutCorrection_%s.root",cutNumber[i].Data(), optionEnergy.Data(), cutNumber[i].Data());
            if (mode == 10) FileNameUnCorrectedMCPi0[i]     = Form("%s/%s/Pi0_MC_GammaMergedWithoutCorrection_%s.root",cutNumber[i].Data(), optionEnergy.Data(), cutNumber[i].Data());

            cout<< FileNameUnCorrectedMCPi0[i] << endl;
            fileUnCorrectedMCPi0[i]                         = new TFile(FileNameUnCorrectedMCPi0[i]);
            if (fileUnCorrectedMCPi0[i]->IsZombie()) 
                enableTriggerRejecCompMC                    = kFALSE;
            else 
                enableTriggerRejecCompMC                    = kTRUE;
        }
        deltaRapid[i]                                       =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
        cout << "For Trigger: " << triggerName[i] << " using rapidity: " <<  rapidityRange.Data() << endl;
        
        histoEventQualtity[i]                               = (TH1F*)fileCorrectedPi0[i]->Get("NEvents");
        histoCorrectedYieldPi0[i]                           = (TH1D*)fileCorrectedPi0[i]->Get(nameCorrectedYield.Data());
        histoCorrectedYieldPi0[i]->SetName(Form("CorrectedYield_%s",cutNumber[i].Data()));
        histoEfficiencyPi0[i]                               = (TH1D*)fileCorrectedPi0[i]->Get(nameEfficiency.Data());
        histoEfficiencyPi0[i]->SetName(Form("Efficiency_%s",  cutNumber[i].Data()));
        histoAcceptancePi0[i]                               = (TH1D*)fileCorrectedPi0[i]->Get(nameAcceptance.Data());
        histoAcceptancePi0[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
        histoAcceptancePi0WOEvtWeights[i]                   = (TH1D*)fileCorrectedPi0[i]->Get(nameAcceptanceWOEvtWeights.Data());
        
        histoEffectCorrPi0FromK0s[i]                        = NULL;
        histoEffectCorrPi0FromK0s[i]                        = (TH1D*)fileCorrectedPi0[i]->Get("RatioSecYieldFromK0SMesonFromToyToRaw");
        if (!histoEffectCorrPi0FromK0s[i] || isMC.CompareTo("MC") == 0 )
            histoEffectCorrPi0FromK0s[i]                    = (TH1D*)fileCorrectedPi0[i]->Get("RatioSecYieldFromK0SMesonToRaw");
        if (histoEffectCorrPi0FromK0s[i])
            hasSecCorrFac                                   = kTRUE;
        histoSecEffiPi0FromK0s[i]                           = NULL;
        histoSecEffiPi0FromK0s[i]                           = (TH1D*)fileCorrectedPi0[i]->Get(nameSecFK0sEfficiency.Data());
        if (histoSecEffiPi0FromK0s[i])
            hasSecEffi                                      = kTRUE;
        
        if(histoAcceptancePi0WOEvtWeights[i]) histoAcceptancePi0WOEvtWeights[i]->SetName(Form("AcceptanceWOEvtWeights_%s",  cutNumber[i].Data()));
        if (mode == 10){
            histoPurityPi0[i]                               = (TH1D*)fileCorrectedPi0[i]->Get("MesonPurity");
            histoPurityPi0[i]->SetName(Form("Purity_%s",  cutNumber[i].Data()));
        }    

        histoRawYieldPi0[i]                                 = (TH1D*)fileCorrectedPi0[i]->Get("histoYieldMesonPerEvent");
        histoRawYieldPi0[i]->SetName(Form("RAWYieldPerEvent_%s",cutNumber[i].Data()));
        if (mode != 10){
            histoMassPi0Data[i]                                 = (TH1D*)fileCorrectedPi0[i]->Get("histoMassMeson");
            histoMassPi0Data[i]->SetName(Form("Pi0_Mass_data_%s",cutNumber[i].Data()));
            histoMassPi0MC[i]                                   = (TH1D*)fileCorrectedPi0[i]->Get(nameMassMC.Data());
            histoMassPi0MC[i]->SetName(Form("Pi0_Mass_MC_%s",cutNumber[i].Data()));
            histoWidthPi0Data[i]                                = (TH1D*)fileCorrectedPi0[i]->Get("histoFWHMMeson");
            histoWidthPi0Data[i]->SetName(Form("Pi0_Width_data_%s",cutNumber[i].Data()));
            histoWidthPi0MC[i]                                  = (TH1D*)fileCorrectedPi0[i]->Get(nameWidthMC.Data());
            histoWidthPi0MC[i]->SetName(Form("Pi0_Width_MC_%s",cutNumber[i].Data()));
            histoInvMassSig[i]                                  = (TH1D*)fileCorrectedPi0[i]->Get(Form("InvMassSig_PtBin%02d",exampleBin));
            if (histoInvMassSig[i]) histoInvMassSig[i]->SetName(Form("Pi0_InvMassSig_Example_%s",triggerName[i].Data()));
            histoInvMassSigPlusBG[i]                            = (TH1D*)fileCorrectedPi0[i]->Get(Form("InvMassSigPlusBG_PtBin%02d",exampleBin));
            if (histoInvMassSigPlusBG[i]) histoInvMassSigPlusBG[i]->SetName(Form("Pi0_InvMassSigPlusBG_Example_%s",triggerName[i].Data()));
            histoInvMassBG[i]                                   = (TH1D*)fileCorrectedPi0[i]->Get(Form("InvMassBG_PtBin%02d",exampleBin));
            if (histoInvMassBG[i]) histoInvMassBG[i]->SetName(Form("Pi0_InvMassBG_Example_%s",triggerName[i].Data()));
            fitInvMassSig[i]                                    = (TF1*)fileCorrectedPi0[i]->Get(Form("FitInvMassSig_PtBin%02d",exampleBin));
            if (fitInvMassSig[i]) fitInvMassSig[i]->SetName(Form("Pi0_InvMassSigFit_Example_%s",triggerName[i].Data()));
        }
        if (cutNumberBaseEff[i].CompareTo("bla") != 0){
            FileNameEffBasePi0[i]                           = Form("%s/%s/Pi0_MC_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), cutNumberBaseEff[i].Data());
            if (mode == 10) FileNameEffBasePi0[i]           = Form("%s/%s/Pi0_MC_GammaMergedCorrection_%s.root", cutNumber[i].Data(), optionEnergy.Data(), cutNumberBaseEff[i].Data());
            fileEffBasePi0[i]                               = new TFile(FileNameEffBasePi0[i]);
            if (fileEffBasePi0[i]->IsZombie()){
                enableTriggerEffPi0[i]                         = kFALSE;
                cout << "Didn't find the effi base file for " << triggerName[i].Data() << endl;
                cout << "ABORTING: as base effi file was requested" << endl;
                
            } else {
                enableTriggerEffPi0[i]                         = kTRUE;
                enableTriggerEffPi0All                         = kTRUE;
            }    
            if (enableTriggerEffPi0[i]){
                TString effiNameBase                        = "TrueMesonEffiPt";
                if (mode == 10) 
                    effiNameBase                            = nameEfficiency;
                TH1D* histoEffiPi0Temp                      = (TH1D*)fileCorrectedPi0[i]->Get(effiNameBase.Data());
                TH1D* histoEffiBasePi0Temp                  = (TH1D*)fileEffBasePi0[i]->Get(effiNameBase.Data());
                histoEffBasePi0[i]                          = (TH1D*)fileEffBasePi0[i]->Get(nameEfficiency.Data());
                histoEffBasePi0[i]->SetName(Form("EfficiencyBase_%s",  cutNumber[i].Data()));
                histoTriggerEffPi0[i]                       = (TH1D*)histoEffiPi0Temp->Clone(Form("TriggerEfficiency_%s", cutNumber[i].Data()));
                histoTriggerEffPi0[i]->Divide(histoTriggerEffPi0[i],histoEffiBasePi0Temp,1.,1.,"B");
                histoEffTimesAccPi0[i]                      = (TH1D*)histoEffBasePi0[i]->Clone(Form("EffTimeAcc_%s",  cutNumber[i].Data()));
                if(histoAcceptancePi0WOEvtWeights[i]){
                  histoAcceptancePi0WOEvtWeights[i]->Sumw2();
                  histoEffTimesAccPi0[i]->Multiply(histoAcceptancePi0WOEvtWeights[i]);
                  histoEfficiencyPi0[i]->Multiply(histoAcceptancePi0[i]);
                  histoEfficiencyPi0[i]->Divide(histoEfficiencyPi0[i],histoAcceptancePi0WOEvtWeights[i],1.,1.,"B");
                  histoAcceptancePi0[i] = histoAcceptancePi0WOEvtWeights[i];
                  histoAcceptancePi0[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
                }else histoEffTimesAccPi0[i]->Multiply(histoAcceptancePi0[i]);
                histoEffTimesAccPi0[i]->Scale(deltaRapid[i]*2*TMath::Pi());
            } else {
                histoEffBasePi0[i]                          = NULL;
                histoTriggerEffPi0[i]                       = NULL;
            }
        } else {
            histoEffTimesAccPi0[i]                      = (TH1D*)histoEfficiencyPi0[i]->Clone(Form("EffTimeAcc_%s",  cutNumber[i].Data()));
            if(histoAcceptancePi0WOEvtWeights[i]){
              histoAcceptancePi0WOEvtWeights[i]->Sumw2();
              histoEffTimesAccPi0[i]->Multiply(histoAcceptancePi0WOEvtWeights[i]);
              histoEfficiencyPi0[i]->Multiply(histoAcceptancePi0[i]);
              histoEfficiencyPi0[i]->Divide(histoEfficiencyPi0[i],histoAcceptancePi0WOEvtWeights[i],1.,1.,"B");
              histoAcceptancePi0[i] = histoAcceptancePi0WOEvtWeights[i];
              histoAcceptancePi0[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
            }else histoEffTimesAccPi0[i]->Multiply(histoAcceptancePi0[i]);
            histoEffTimesAccPi0[i]->Scale(deltaRapid[i]*2*TMath::Pi());
        }
        
        //Scale spectrum to MBOR
        if (optionEnergy.CompareTo("2.76TeV")==0 && 
            (triggerName[i].Contains("INT7")|| triggerName[i].Contains("EMC7") || triggerName[i].Contains("EG1") || triggerName[i].Contains("EG2")) && 
            isMC.CompareTo("data") == 0){
            histoCorrectedYieldPi0[i]->Scale(0.8613) ;
            histoRawYieldPi0[i]->Scale(0.8613) ;
        }
        if (triggerName[i].CompareTo("INT7") != 0 && triggerName[i].CompareTo("MB") != 0 && triggerName[i].CompareTo("INT1") != 0){
            nRealTriggers++;
        }
        
        histoMCInputPi0[i]                                  = (TH1D*)fileCorrectedPi0[i]->Get(nameMCYield.Data());
        histoMCInputPi0[i]->SetName(Form("Pi0_Input_Reweighted_%s",cutNumber[i].Data()));
        histoMCInputPi0[i]->Sumw2();
        histoMCInputPi0[i]->Rebin(4);
        histoMCInputPi0[i]->Scale(1./4);
        //***************************************************************************************************************
        //****************************** Calculate trigger rejection factors ********************************************
        //***************************************************************************************************************
        if ((mode == 0 || mode == 2 || mode == 3 || mode == 4 || mode == 5 || mode == 10) && hasClusterOutput ){
            histoRawClusterPt[i]                        = (TH1D*)fileUnCorrectedPi0[i]->Get("ClusterPtPerEvent");
            if (!histoRawClusterPt[i]){ 
                cout << "INFO: couldn't find cluster input, disabeling it!" << endl;
                hasClusterOutput                        = kFALSE;
                triggRejecFac[i][trigSteps[i][0]]       = 1;
                triggRejecFacErr[i][trigSteps[i][0]]    = 0;
            } else {    
                histoRawClusterPt[i]->SetName(Form("ClusterPtPerEvent_%s",cutNumber[i].Data()));
                histoRatioRawClusterPt[i]                   = (TH1D*)histoRawClusterPt[i]->Clone(Form("RatioCluster_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()));
                histoRatioRawClusterPt[i]->Divide(histoRatioRawClusterPt[i],histoRawClusterPt[trigSteps[i][0]],1.,1.,"");
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
                    histoMCRawClusterPt[i]                  = (TH1D*)fileUnCorrectedMCPi0[i]->Get("ClusterPtPerEvent");
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
    if (optionEnergy.CompareTo("8TeV") == 0){
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
        if (mode == 4 && optionEnergy.CompareTo("pPb_5.023TeV") == 0)
            maxTriggReject = 200;
        if (mode == 10)
            maxTriggReject = 5200;
        
        TH2F * histo2DTriggReject;
        histo2DTriggReject = new TH2F("histo2DTriggReject","histo2DTriggReject",1000,0., maxPtGlobalCluster,10000,minTriggReject, maxTriggReject);
        SetStyleHistoTH2ForGraphs(histo2DTriggReject, "#it{p}_{T} (GeV/#it{c})","#it{R}_{Trig}", //"#frac{N_{clus,trig A}/N_{Evt, trig A}}{N_{clus,trig B}/N_{Evt,trig B}}", 
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
            }
        }
        legendTriggReject->Draw();
        histo2DTriggReject->Draw("same,axis");
        TLatex *labelPerfTriggRejec  = NULL;
        if (isMC.CompareTo("MC") == 0)
            labelPerfTriggRejec  = new TLatex(0.11, 0.93,"ALICE simulation");
        else 
            labelPerfTriggRejec  = new TLatex(0.11, 0.93,"ALICE performance");
        
        SetStyleTLatex( labelPerfTriggRejec, textSizeSpectra2,4);
        labelPerfTriggRejec->Draw();

        TLatex *labelPerfTriggFitRange = new TLatex(0.523, 0.12+(0.9*(nrOfTrigToBeComb-2+1)*textSizeSpectra2)+0.01, "Fit range (GeV/#it{c})");
        SetStyleTLatex( labelPerfTriggFitRange, textSizeSpectra2,4);
        labelPerfTriggFitRange->Draw();

        TLatex *labelPerfTriggRejecFac = new TLatex(0.753, 0.12+(0.9*(nrOfTrigToBeComb-2+1)*textSizeSpectra2)+0.01, "Trigger rejection");
        SetStyleTLatex( labelPerfTriggRejecFac, textSizeSpectra2,4);
        labelPerfTriggRejecFac->Draw();
        
        canvasTriggerReject->Update();
        canvasTriggerReject->SaveAs(Form("%s/%s_TriggerRejectionFactors.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
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
        if (mode == 4 && optionEnergy.CompareTo("pPb_5.023TeV") == 0)
            maxTriggRejectLin = 100;
        if (mode == 4 && optionEnergy.CompareTo("2.76TeV") == 0)
            maxTriggRejectLin = 2000;
        if (mode == 10)
            maxTriggRejectLin = 3000;

        if( optionEnergy.CompareTo("8TeV")==0 ){
        if (mode == 2 || mode == 4 || mode == 10 )
            maxTriggRejectLin = 310;
        }
        
        TH2F * histo2DTriggRejectLinear;
        histo2DTriggRejectLinear = new TH2F("histo2DTriggRejectLinear","histo2DTriggRejectLinear",1000,0., maxPtGlobalCluster,15000,minTriggRejectLin, maxTriggRejectLin);
        SetStyleHistoTH2ForGraphs(histo2DTriggRejectLinear, "#it{p}_{T} (GeV/#it{c})","#it{R}_{Trig}", //"#frac{N_{clus,trig A}/N_{Evt, trig A}}{N_{clus,trig B}/N_{Evt,trig B}}", 
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
        //************************************Plotting unscaled invariant raw-yield Pi0 *********************************
        //***************************************************************************************************************
        
        TCanvas* canvasClusterYield = new TCanvas("canvasClusterYield","",0,0,1000,1350);// gives the page size
        DrawGammaCanvasSettings( canvasClusterYield, 0.16, 0.02, 0.015, 0.07);
        canvasClusterYield->SetLogy();
        
        Double_t minClusYieldUnscaled    = 7e-9;
        Double_t maxClusYieldUnscaled    = 5;
        if (mode == 10) {
            minClusYieldUnscaled         = 2e-10;
            maxClusYieldUnscaled         = 4;
        } else if (mode == 4) {
            minClusYieldUnscaled         = 7e-10;
            maxClusYieldUnscaled         = 5;
        }
        if(optionEnergy.CompareTo("8TeV")==0){
        if (mode == 2) {
            minClusYieldUnscaled         = 7e-10;
            maxClusYieldUnscaled         = 5;
        } else if (mode == 4) {
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
        if (mode == 4 || mode == 2){
            labelDetProcClus = new TLatex(0.2, 0.12,detectionProcessClus.Data());
        } else if (mode == 10) {
            labelDetProcClus = new TLatex(0.2, 0.12,"rec. with EMCal");
        } else {
            labelDetProcClus = new TLatex(0.2, 0.12,"rec. with EMCal");
        }    
            
        SetStyleTLatex( labelDetProcClus, 0.85*textSizeSpectra,4);
        labelDetProcClus->Draw();

        canvasClusterYield->Update();
        canvasClusterYield->SaveAs(Form("%s/Cluster_%s_YieldUnscaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        if (!enableEta) delete canvasClusterYield;    
    }
    
    //***************************************************************************************************************
    //************************************Plotting efficiencies Pi0 *************************************************
    //***************************************************************************************************************
    TCanvas* canvasEffi = new TCanvas("canvasEffi","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasEffi, 0.09, 0.017, 0.015, 0.08);
    canvasEffi->SetLogy(1);

    Double_t minEffiPi0 = 1e-4;
    Double_t maxEffiPi0 = 1e-1;
    if (mode == 4){
        maxEffiPi0      = 8e-1;
    } else if (mode == 10){
        maxEffiPi0      = 8e-1;
        minEffiPi0      = 1e-2;
    } else if (mode == 2){
        if(optionEnergy.CompareTo("8TeV")==0)
            minEffiPi0  = 5e-5;
    } else if (mode == 0){
        maxEffiPi0      = 8e-3;
        minEffiPi0      = 5e-5;
    }
    
    TH2F * histo2DEffiPi0;
    histo2DEffiPi0 = new TH2F("histo2DEffiPi0","histo2DEffiPi0",1000,0., maxPtGlobalPi0,10000,minEffiPi0, maxEffiPi0);
    SetStyleHistoTH2ForGraphs(histo2DEffiPi0, "#it{p}_{T} (GeV/#it{c})","#epsilon_{#pi^{0}}",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
    histo2DEffiPi0->DrawCopy(); 

    Double_t minXLegendEffi = 0.62;
    Int_t nColumnsEffi      = 2; 
    if (nrOfTrigToBeComb > 6){
        minXLegendEffi = 0.4;
        nColumnsEffi   = 3; 
    }
    Double_t minYLegendEffi = 0.13;
    Double_t maxYLegendEffi = minYLegendEffi+(1.05*nrOfTrigToBeComb/nColumnsEffi*0.85*textSizeSpectra);
    
    TLegend* legendEffiPi0 = GetAndSetLegend2(minXLegendEffi, minYLegendEffi, 0.95, maxYLegendEffi ,28);
    legendEffiPi0->SetNColumns(nColumnsEffi);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoEfficiencyPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoEfficiencyPi0[i]->DrawCopy("e1,same"); 
        legendEffiPi0->AddEntry(histoEfficiencyPi0[i],triggerNameLabel[i].Data(),"p"); 
    }
    legendEffiPi0->Draw();

    TLatex *labelEnergyEffi = new TLatex(0.62, maxYLegendEffi+0.02+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
    SetStyleTLatex( labelEnergyEffi, 0.85*textSizeSpectra,4);
    labelEnergyEffi->Draw();

    TLatex *labelPi0Effi = new TLatex(0.62, maxYLegendEffi+0.02+0.99*textSizeSpectra*0.85,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelPi0Effi, 0.85*textSizeSpectra,4);
    labelPi0Effi->Draw();

    TLatex *labelDetProcEffi = new TLatex(0.62, maxYLegendEffi+0.02,detectionProcess.Data());
    SetStyleTLatex( labelDetProcEffi, 0.85*textSizeSpectra,4);
    labelDetProcEffi->Draw();

    canvasEffi->Update();
    canvasEffi->SaveAs(Form("%s/Pi0_Efficiency.%s",outputDir.Data(),suffix.Data()));

    if (hasSecEffi){
        canvasEffi->cd();
        Double_t minEffiSecPi0      = minEffiPi0;
        Double_t maxEffiSecPi0      = maxEffiPi0;
        if (mode == 10){
            maxEffiSecPi0           = 10*maxEffiPi0;
            minEffiSecPi0           = 5*minEffiPi0;
        }    
        TH2F * histo2DEffiSecPi0;
        histo2DEffiSecPi0 = new TH2F("histo2DEffiSecPi0","histo2DEffiSecPi0",1000,0., maxPtGlobalPi0,10000,minEffiSecPi0, maxEffiSecPi0);
        SetStyleHistoTH2ForGraphs(histo2DEffiSecPi0, "#it{p}_{T} (GeV/#it{c})","#epsilon_{sec #pi^{0} from K^{0}_{s}}",
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
        histo2DEffiSecPi0->DrawCopy(); 
    
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (histoSecEffiPi0FromK0s[i]){
                DrawGammaSetMarker(histoSecEffiPi0FromK0s[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoSecEffiPi0FromK0s[i]->DrawCopy("e1,same"); 
            }    
        }
        legendEffiPi0->Draw();

        labelEnergyEffi->Draw();
        labelPi0Effi->Draw();
        labelDetProcEffi->Draw();

        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/Pi0_SecEfficiencyFromK0s.%s",outputDir.Data(),suffix.Data()));        
    }    

    if (hasSecCorrFac){
        canvasEffi->cd();
        canvasEffi->SetLogy(0);
        canvasEffi->SetLeftMargin(0.1);
        canvasEffi->SetTopMargin(0.04);
        Double_t minXLegendEffecSec = 0.62;
        Int_t nColumnsEffecSec      = 2; 
        if (nrOfTrigToBeComb > 6){
            minXLegendEffecSec = 0.4;
            nColumnsEffecSec   = 3; 
        }
        Double_t maxYEffSecCorr     = 0.3;
        if (mode == 2)
            maxYEffSecCorr          = 0.04;
        else if (mode == 4)
            maxYEffSecCorr          = 0.12;
        else if (mode == 10)
            maxYEffSecCorr          = 0.15;
        
        Double_t maxYLegendEffecSec = 0.84;
        Double_t minYLegendEffecSec = maxYLegendEffecSec-(1.05*nrOfTrigToBeComb/nColumnsEffecSec*0.85*textSizeSpectra);

        TH2F * histo2DEffectiveSecCorr;
        histo2DEffectiveSecCorr = new TH2F("histo2DEffectiveSecCorr","histo2DEffectiveSecCorr",1000,0., maxPtGlobalPi0,10000,0, maxYEffSecCorr);
        SetStyleHistoTH2ForGraphs(histo2DEffectiveSecCorr, "#it{p}_{T} (GeV/#it{c})","r_{sec #pi^{0} from K^{0}_{s}}",
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.17);
        histo2DEffectiveSecCorr->DrawCopy(); 

        TLegend* legendEffectSec = GetAndSetLegend2(minXLegendEffecSec, minYLegendEffecSec, 0.95, maxYLegendEffecSec ,28, nColumnsEffecSec);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (histoEffectCorrPi0FromK0s[i]){
                DrawGammaSetMarker(histoEffectCorrPi0FromK0s[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoEffectCorrPi0FromK0s[i]->DrawCopy("e1,same"); 
                legendEffectSec->AddEntry(histoEffectCorrPi0FromK0s[i],triggerNameLabel[i].Data(),"p"); 
            }    
        }
        legendEffectSec->Draw();

        TLatex *labelEnergyEffSec = new TLatex(0.62, maxYLegendEffecSec+0.02+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyEffSec, 0.85*textSizeSpectra,4);
        labelEnergyEffSec->Draw();

        TLatex *labelPi0EffSec = new TLatex(0.62, maxYLegendEffecSec+0.02+0.99*textSizeSpectra*0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelPi0EffSec, 0.85*textSizeSpectra,4);
        labelPi0EffSec->Draw();

        TLatex *labelDetProcEffSec = new TLatex(0.62, maxYLegendEffecSec+0.02,detectionProcess.Data());
        SetStyleTLatex( labelDetProcEffSec, 0.85*textSizeSpectra,4);
        labelDetProcEffSec->Draw();

        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/Pi0_EffectiveSecCorrFromK0s.%s",outputDir.Data(),suffix.Data()));        
        canvasEffi->SetLogy(1);
        canvasEffi->SetLeftMargin(0.09);
        canvasEffi->SetTopMargin(0.015);
    }    
    
    //***************************************************************************************************************
    //************************************ Plotting trigger efficiencies Pi0 ****************************************
    //***************************************************************************************************************
    if (enableTriggerEffPi0All){
        TCanvas* canvasTriggerEffi = new TCanvas("canvasTriggerEffi","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasTriggerEffi, 0.09, 0.017, 0.015, 0.08);
        canvasTriggerEffi->SetLogy(0);

        Double_t minEffiTrigPi0         = 0;
        Double_t maxEffiTrigPi0         = 1.1;
        if(optionEnergy.CompareTo("8TeV") == 0 && mode == 2) maxEffiTrigPi0 = 1.5;
        
        TH2F * histo2DTriggerEffiPi0;
        histo2DTriggerEffiPi0 = new TH2F("histo2DTriggerEffiPi0","histo2DTriggerEffiPi0",1000,0., maxPtGlobalPi0,10000,minEffiTrigPi0, maxEffiTrigPi0);
        SetStyleHistoTH2ForGraphs(histo2DTriggerEffiPi0, "#it{p}_{T} (GeV/#it{c})","#epsilon_{Trigger, #pi^{0}}", 
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
        histo2DTriggerEffiPi0->DrawCopy(); 

        TLegend* legendTriggerEffiPi0 = GetAndSetLegend2(0.62, 0.165, 0.95, 0.165+(1.05*nRealTriggers/2*0.85*textSizeSpectra),28);
        legendTriggerEffiPi0->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (enableTriggerEffPi0[i]){
                DrawGammaSetMarker(histoTriggerEffPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoTriggerEffPi0[i]->DrawCopy("e1,same"); 
                legendTriggerEffiPi0->AddEntry(histoTriggerEffPi0[i],triggerNameLabel[i].Data(),"p"); 
            }    
        }
        legendTriggerEffiPi0->Draw();

        labelEnergyEffi->Draw();
        labelPi0Effi->Draw();
        labelDetProcEffi->Draw();

        canvasTriggerEffi->Update();
        canvasTriggerEffi->SaveAs(Form("%s/Pi0_TriggerEfficiency.%s",outputDir.Data(),suffix.Data()));
        delete canvasTriggerEffi;
        
        //***************************************************************************************************************
        //******************* Plotting efficiencies Pi0 without trigger efficiency folded in  ***************************
        //***************************************************************************************************************
        canvasEffi->cd();
        histo2DEffiPi0->DrawCopy(); 

        TLegend* legendEffiPi0W0TriggEff = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
        legendEffiPi0W0TriggEff->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if ( triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1") ){
                DrawGammaSetMarker(histoEfficiencyPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoEfficiencyPi0[i]->DrawCopy("e1,same"); 
                legendEffiPi0W0TriggEff->AddEntry(histoEfficiencyPi0[i],triggerNameLabel[i].Data(),"p"); 
            } else {
                if (enableTriggerEffPi0[i]){
                    DrawGammaSetMarker(histoEffBasePi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoEffBasePi0[i]->DrawCopy("e1,same"); 
                    legendEffiPi0W0TriggEff->AddEntry(histoEffBasePi0[i],triggerNameLabel[i].Data(),"p"); 
                }
            }    
        }
        legendEffiPi0W0TriggEff->Draw();

        labelEnergyEffi->Draw();
        labelPi0Effi->Draw();
        labelDetProcEffi->Draw();

        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/Pi0_EfficiencyW0TriggEff.%s",outputDir.Data(),suffix.Data()));
    }
//     if (!enableEta) delete canvasEffi;
    
    //***************************************************************************************************************
    //************************************Plotting acceptance Pi0 *************************************************
    //***************************************************************************************************************
    TCanvas* canvasAcc = new TCanvas("canvasAcc","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasAcc, 0.1, 0.017, 0.015, 0.08);
    canvasAcc->SetLogy(0);

    Double_t minAccPi0 = 0.15;
    Double_t maxAccPi0 = 0.3;
    if (mode == 0){
        maxAccPi0       = 1.05;
        minAccPi0       = 0.7;
    }

    if(optionEnergy.CompareTo("8TeV")==0){
      if (mode == 2){
        minAccPi0 = 0.18;
      }else if(mode == 4){
        maxAccPi0 = 0.26;
      }
    }

    TH2F * histo2DAccPi0;
    histo2DAccPi0 = new TH2F("histo2DAccPi0","histo2DAccPi0",1000,0., maxPtGlobalPi0,10000,minAccPi0, maxAccPi0);
    SetStyleHistoTH2ForGraphs(histo2DAccPi0, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}", 
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.25);
    histo2DAccPi0->DrawCopy(); 

    Double_t minXLegendAcc = 0.62;
    Int_t nColumnsAcc      = 2; 
    if (nrOfTrigToBeComb > 6){
        minXLegendAcc = 0.4;
        nColumnsAcc   = 3; 
    }
    Double_t minYLegendAcc = 0.13;
    Double_t maxYLegendAcc = minYLegendAcc+(1.05*nrOfTrigToBeComb/nColumnsAcc*0.85*textSizeSpectra);
    
    TLegend* legendAccPi0 = GetAndSetLegend2(minXLegendAcc, minYLegendAcc, 0.95,maxYLegendAcc,28);
    legendAccPi0->SetNColumns(nColumnsAcc);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoAcceptancePi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoAcceptancePi0[i]->DrawCopy("e1,same"); 
        legendAccPi0->AddEntry(histoAcceptancePi0[i],triggerNameLabel[i].Data(),"p");
    }
    legendAccPi0->Draw();

    labelEnergyEffi->Draw();
    labelPi0Effi->Draw();
    labelDetProcEffi->Draw();

    
    canvasAcc->Update();
    canvasAcc->SaveAs(Form("%s/Pi0_Acceptance.%s",outputDir.Data(),suffix.Data()));

    if (mode == 10){
        //***************************************************************************************************************
        //************************************Plotting efficiencies Pi0 *************************************************
        //***************************************************************************************************************
        TCanvas* canvasPurity = new TCanvas("canvasPurity","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasPurity, 0.09, 0.017, 0.015, 0.08);
        canvasPurity->SetLogy(0);

        Double_t minPurityPi0 = 0.6;
        Double_t maxPurityPi0 = 1.02;

        TH2F * histo2DPurityPi0;
        histo2DPurityPi0 = new TH2F("histo2DPurityPi0","histo2DPurityPi0",1000,0., maxPtGlobalPi0,10000,minPurityPi0, maxPurityPi0);
        SetStyleHistoTH2ForGraphs(histo2DPurityPi0, "#it{p}_{T} (GeV/#it{c})","#it{P}_{#pi^{0}}", 
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
        histo2DPurityPi0->DrawCopy(); 

        Double_t minXLegendPur = 0.62;
        Int_t nColumnsPur      = 2; 
        if (nrOfTrigToBeComb > 6){
            minXLegendPur = 0.4;
            nColumnsPur   = 3; 
        }
        Double_t minYLegendPur = 0.13;
        Double_t maxYLegendPur = minYLegendPur+(1.05*nrOfTrigToBeComb/nColumnsPur*0.85*textSizeSpectra);
        
        TLegend* legendPurityPi0 = GetAndSetLegend2(minXLegendPur, minYLegendPur, 0.95, maxYLegendPur,28);
        legendPurityPi0->SetNColumns(nColumnsPur);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            DrawGammaSetMarker(histoPurityPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoPurityPi0[i]->DrawCopy("e1,same"); 
            legendPurityPi0->AddEntry(histoPurityPi0[i],triggerNameLabel[i].Data(),"p"); 
        }
        legendPurityPi0->Draw();

        TLatex *labelEnergyPurity = new TLatex(0.62, maxYLegendPur+0.02+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyPurity, 0.85*textSizeSpectra,4);
        labelEnergyPurity->Draw();

        TLatex *labelPi0Purity = new TLatex(0.62, maxYLegendPur+0.02+0.99*textSizeSpectra*0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelPi0Purity, 0.85*textSizeSpectra,4);
        labelPi0Purity->Draw();

        TLatex *labelDetProcPurity = new TLatex(0.62, maxYLegendPur+0.02,detectionProcess.Data());
        SetStyleTLatex( labelDetProcPurity, 0.85*textSizeSpectra,4);
        labelDetProcPurity->Draw();

        canvasPurity->Update();
        canvasPurity->SaveAs(Form("%s/Pi0_Purity.%s",outputDir.Data(),suffix.Data()));
        
    }    
    
    //***************************************************************************************************************
    //**************************** Mass and Width general plotting definitions **************************************
    //***************************************************************************************************************
    TCanvas* canvasMass         = new TCanvas("canvasMass","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasMass, 0.11, 0.017, 0.015, 0.08);
    Double_t minMassPi0         = 0.120;
    Double_t maxMassPi0         = 0.160;
    if(optionEnergy.CompareTo("8TeV")==0){
      if(mode == 2){
        minMassPi0              = 0.123;
        maxMassPi0              = 0.150;
      }else if(mode == 4){
        maxMassPi0              = 0.180;
      }
    }

    if (mode == 0){
        maxMassPi0              = 0.142;
        minMassPi0              = 0.130;
    }    
    
    TH2F * histo2DMassPi0       = new TH2F("histo2DMassPi0","histo2DMassPi0",1000,0., maxPtGlobalPi0,10000,minMassPi0, maxMassPi0);
    SetStyleHistoTH2ForGraphs(histo2DMassPi0, "#it{p}_{T} (GeV/#it{c})","#it{M}_{#pi^{0}} (Mev/#it{c}^{2})", 
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.4);
    TLegend* legendMassPi0      = GetAndSetLegend2(0.52, 0.88, 0.95, 0.88+(1.05*4/2*0.85*textSizeSpectra),28);
    TLegend* legendMassRedPi0   = GetAndSetLegend2(0.52, 0.88, 0.95, 0.88+(1.05*4/2*0.85*textSizeSpectra),28);
    TLatex *labelEnergyMass     = new TLatex(0.14, 0.85+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
    TLatex *labelPi0Mass        = new TLatex(0.14, 0.85+0.99*textSizeSpectra*0.85,"#pi^{0} #rightarrow #gamma#gamma");
    TLatex *labelDetProcMass    = new TLatex(0.14, 0.85,detectionProcess.Data());
    TLegend* legendMassPi02     = GetAndSetLegend2(0.46, 0.81, 0.80, 0.81+(1.05*8/2*0.85*textSizeSpectra),28);
    TLegend* legendMassRedPi02  = GetAndSetLegend2(0.46, 0.81, 0.80, 0.81+(1.05*8/2*0.85*textSizeSpectra),28);
    
    TCanvas* canvasWidth        = new TCanvas("canvasWidth","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasWidth, 0.09, 0.017, 0.035, 0.08);
    Double_t minWidthPi0        = 0.0;
    Double_t maxWidthPi0        = 0.0295;
    if (mode == 0){
        minWidthPi0             = 0.0;
        maxWidthPi0             = 0.0125;
    }    
    TH2F * histo2DWidthPi0      = new TH2F("histo2DWidthPi0","histo2DWidthPi0",1000,0., maxPtGlobalPi0,10000,minWidthPi0, maxWidthPi0);
    SetStyleHistoTH2ForGraphs(histo2DWidthPi0, "#it{p}_{T} (GeV/#it{c})","#sigma_{#pi^{0}} (Mev/#it{c}^{2})", 
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
    TLegend* legendWidthPi0     = GetAndSetLegend2(0.52, 0.87, 0.95, 0.87+(1.05*4/2*0.85*textSizeSpectra),28);
    TLegend* legendWidthRedPi0  = GetAndSetLegend2(0.52, 0.87, 0.95, 0.87+(1.05*4/2*0.85*textSizeSpectra),28);
    TLatex* labelEnergyWidth    = new TLatex(0.14, 0.84+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
    TLatex* labelPi0Width       = new TLatex(0.14, 0.84+0.99*textSizeSpectra*0.85,"#pi^{0} #rightarrow #gamma#gamma");
    TLatex* labelDetProcWidth   = new TLatex(0.14, 0.84,detectionProcess.Data());
    TLegend* legendWidthPi02    = GetAndSetLegend2(0.46, 0.80, 0.80, 0.80+(1.05*8/2*0.85*textSizeSpectra),28);
    TLegend* legendWidthRedPi02 = GetAndSetLegend2(0.46, 0.80, 0.80, 0.80+(1.05*8/2*0.85*textSizeSpectra),28);
    
    if (mode != 10){
        //***************************************************************************************************************
        //************************************Plotting Mass Pi0 *********************************************************
        //***************************************************************************************************************
        canvasMass->cd();
        histo2DMassPi0->DrawCopy(); 

        legendMassPi0->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
               (optionEnergy.CompareTo("8TeV")==0)
               ){
                DrawGammaSetMarker(histoMassPi0Data[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoMassPi0Data[i]->DrawCopy("e1,same"); 
                legendMassPi0->AddEntry(histoMassPi0Data[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                DrawGammaSetMarker(histoMassPi0MC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                histoMassPi0MC[i]->DrawCopy("e1,same"); 
                legendMassPi0->AddEntry(histoMassPi0MC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
            }    
        }
        legendMassPi0->Draw();

        SetStyleTLatex( labelEnergyMass, 0.85*textSizeSpectra,4);
        labelEnergyMass->Draw();


        SetStyleTLatex( labelPi0Mass, 0.85*textSizeSpectra,4);
        labelPi0Mass->Draw();

        SetStyleTLatex( labelDetProcMass, 0.85*textSizeSpectra,4);
        labelDetProcMass->Draw();

        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/Pi0_%s_Mass1.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        histo2DMassPi0->DrawCopy(); 


        if (optionEnergy.CompareTo("2.76TeV")==0){
            legendMassPi02->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if( optionEnergy.CompareTo("2.76TeV")==0 && (i==1 || i==3 || i==4 || i==5)){
                    DrawGammaSetMarker(histoMassPi0Data[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoMassPi0Data[i]->DrawCopy("e1,same"); 
                    legendMassPi02->AddEntry(histoMassPi0Data[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                    DrawGammaSetMarker(histoMassPi0MC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    histoMassPi0MC[i]->DrawCopy("e1,same"); 
                    legendMassPi02->AddEntry(histoMassPi0MC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                }    
            }
            legendMassPi02->Draw();
            labelEnergyMass->Draw();
            labelPi0Mass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Pi0_%s_Mass2.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        //***************************************************************************************************************
        //************************************Plotting Width Pi0 *********************************************************
        //***************************************************************************************************************
        canvasWidth->cd();
        histo2DWidthPi0->DrawCopy(); 

        legendWidthPi0->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0)
             ){
                DrawGammaSetMarker(histoWidthPi0Data[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoWidthPi0Data[i]->DrawCopy("e1,same"); 
                legendWidthPi0->AddEntry(histoWidthPi0Data[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                DrawGammaSetMarker(histoWidthPi0MC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                histoWidthPi0MC[i]->DrawCopy("e1,same"); 
                legendWidthPi0->AddEntry(histoWidthPi0MC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
            }    
        }
        legendWidthPi0->Draw();

        SetStyleTLatex( labelEnergyWidth, 0.85*textSizeSpectra,4);
        labelEnergyWidth->Draw();

        SetStyleTLatex( labelPi0Width, 0.85*textSizeSpectra,4);
        labelPi0Width->Draw();

        SetStyleTLatex( labelDetProcWidth, 0.85*textSizeSpectra,4);
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Pi0_%s_Width1.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        histo2DWidthPi0->DrawCopy(); 

        if (optionEnergy.CompareTo("2.76TeV")==0) {
            legendWidthPi02->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if(i==1 || i==3 || i==4 || i==5 ){
                    DrawGammaSetMarker(histoWidthPi0Data[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoWidthPi0Data[i]->DrawCopy("e1,same"); 
                    legendWidthPi02->AddEntry(histoWidthPi0Data[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                    DrawGammaSetMarker(histoWidthPi0MC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    histoWidthPi0MC[i]->DrawCopy("e1,same"); 
                    legendWidthPi02->AddEntry(histoWidthPi0MC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                }    
            }
            legendWidthPi02->Draw();
            labelEnergyWidth->Draw();
            labelPi0Width->Draw();
            labelDetProcWidth->Draw();

            canvasWidth->Update();
            canvasWidth->SaveAs(Form("%s/Pi0_%s_Width2.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }    
    }

    //***************************************************************************************************************
    //************************************Plotting unscaled invariant raw-yield Pi0 *********************************
    //***************************************************************************************************************
    textSizePixelSpectra = textSizeSpectra*1000;
    
    TCanvas* canvasRawUnscaled = new TCanvas("canvasRawUnscaled","",0,0,1000,1350);// gives the page size
    DrawGammaCanvasSettings( canvasRawUnscaled, 0.16, 0.02, 0.015, 0.07);
    canvasRawUnscaled->SetLogy();
    
    Double_t minCorrYieldRawUnscaled    = 7e-8;
    Double_t maxCorrYieldRawUnscaled    = 4e-2;
    if (mode == 10) {
        minCorrYieldRawUnscaled         = 2e-12;
        maxCorrYieldRawUnscaled         = 1;
    } else if (mode == 4) {
        minCorrYieldRawUnscaled         = 7e-8;
        maxCorrYieldRawUnscaled         = 1;
    }

    if(optionEnergy.CompareTo("8TeV")==0){
      if(mode == 2){
        minCorrYieldRawUnscaled         = 2e-8;
        maxCorrYieldRawUnscaled         = 8e-3;
      } else if(mode == 4){
        minCorrYieldRawUnscaled         = 1e-7;
        maxCorrYieldRawUnscaled         = 2e-1;
      } else if (mode == 0){
        minCorrYieldRawUnscaled         = 2e-8;
        maxCorrYieldRawUnscaled         = 4e-3;          
      }    
    }
    TH2F * histo2DRawUnscaled       = new TH2F("histo2DRawUnscaled", "histo2DRawUnscaled", 1000, 0., maxPtGlobalPi0, 10000, minCorrYieldRawUnscaled, maxCorrYieldRawUnscaled);
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
    if(optionEnergy.CompareTo("8TeV")==0) minYLegendRaw-=0.03;
    
    TLegend* legendRawUnscaled      = GetAndSetLegend2(minXLegendRaw, minYLegendRaw, 0.95, maxYLegendRaw,0.85*textSizePixelSpectra);
    legendRawUnscaled->SetNColumns(nColumnsRaw);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoRawYieldPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoRawYieldPi0[i]->DrawCopy("e1,same"); 
        legendRawUnscaled->AddEntry(histoRawYieldPi0[i],triggerNameLabel[i].Data(),"p");
    }
    legendRawUnscaled->Draw();
    
    
    TLatex *labelEnergyRawUnscaled  = new TLatex(0.2, 0.12+(2*textSizeSpectra*0.85*0.75),collisionSystem.Data());
    SetStyleTLatex( labelEnergyRawUnscaled, 0.85*textSizeSpectra,4);
    labelEnergyRawUnscaled->Draw();
    
    TLatex *labelPi0RawUnscaled     = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelPi0RawUnscaled, 0.85*textSizeSpectra,4);
    labelPi0RawUnscaled->Draw();
    
    TLatex *labelDetProcRawUnscaled = new TLatex(0.2, 0.12,detectionProcess.Data());
    SetStyleTLatex( labelDetProcRawUnscaled, 0.85*textSizeSpectra,4);
    labelDetProcRawUnscaled->Draw();

    canvasRawUnscaled->Update();
    canvasRawUnscaled->SaveAs(Form("%s/Pi0_%s_RawYieldUnscaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
    if (!enableEta) delete canvasRawUnscaled;    
    
    //***************************************************************************************************************
    //************************************Plotting unscaled invariant yield Pi0 *************************************
    //***************************************************************************************************************    
    TCanvas* canvasCorrUnscaled = new TCanvas("canvasCorrUnscaled","",0,0,1000,1350);// gives the page size
    DrawGammaCanvasSettings( canvasCorrUnscaled, 0.15, 0.017, 0.015, 0.07);
    canvasCorrUnscaled->SetLogy();
    
    Double_t minCorrYieldUnscaled   = 2e-10;
    Double_t maxCorrYieldUnscaled   = 1e2;
    if (mode == 10) {
        minCorrYieldUnscaled        = 2e-12;
        maxCorrYieldUnscaled        = 1;
    } else if (mode == 0){
        minCorrYieldUnscaled        = 1e-8;
        maxCorrYieldUnscaled        = 1e2;        
    }

    if(optionEnergy.CompareTo("8TeV")==0){
      if(mode == 2){
        minCorrYieldUnscaled        = 1e-8;
        maxCorrYieldUnscaled        = 1;
      }else if(mode == 4){
        minCorrYieldUnscaled        = 2e-8;
        maxCorrYieldUnscaled        = 0.2;
      }
    }
    
    TH2F * histo2DInvYieldUnscaled = new TH2F("histo2DInvYieldUnscaled","histo2DInvYieldUnscaled",1000,0., maxPtGlobalPi0,10000,minCorrYieldUnscaled,maxCorrYieldUnscaled);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldUnscaled, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 
                            0.85*textSizeSpectra,0.04, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.55);
    histo2DInvYieldUnscaled->DrawCopy(); 

    TLegend* legendUnscaled = GetAndSetLegend2(minXLegendRaw, minYLegendRaw, 0.95, maxYLegendRaw,0.85*textSizePixelSpectra);
    legendUnscaled->SetNColumns(nColumnsRaw);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoCorrectedYieldPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoCorrectedYieldPi0[i]->DrawCopy("e1,same"); 
        legendUnscaled->AddEntry(histoCorrectedYieldPi0[i],triggerNameLabel[i].Data(),"p");
    }
    legendUnscaled->Draw();
    
    
    TLatex *labelEnergyUnscaled = new TLatex(0.2, 0.12+(2*textSizeSpectra*0.85*0.75),collisionSystem.Data());
    SetStyleTLatex( labelEnergyUnscaled, 0.85*textSizeSpectra,4);
    labelEnergyUnscaled->Draw();
    
    TLatex *labelPi0Unscaled = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelPi0Unscaled, 0.85*textSizeSpectra,4);
    labelPi0Unscaled->Draw();
    
    TLatex *labelDetProcUnscaled = new TLatex(0.2, 0.12,detectionProcess.Data());
    SetStyleTLatex( labelDetProcUnscaled, 0.85*textSizeSpectra,4);
    labelDetProcUnscaled->Draw();

    canvasCorrUnscaled->Update();
    canvasCorrUnscaled->SaveAs(Form("%s/Pi0_%s_CorrectedYieldUnscaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
    if (!enableEta) delete canvasCorrUnscaled;
    
    //***************************************************************************************************************
    //******************************* Scaling corrected yield by trigger rejection factors **************************
    //***************************************************************************************************************
    
    TH1D*     histoCorrectedYieldPi0Scaled                  [MaxNumberOfFiles];
    TH1D*     histoCorrectedYieldPi0ScaledMasked            [MaxNumberOfFiles];
    histoCorrectedYieldPi0Scaled[0]                         = (TH1D*)histoCorrectedYieldPi0[0]->Clone(Form("CorrectedYieldPi0Scaled_%s", triggerName[0].Data()));
    for (Int_t i = 1; i< nrOfTrigToBeComb; i++){
        cout << triggerName[i].Data() << endl;
        histoCorrectedYieldPi0Scaled[i] = (TH1D*)histoCorrectedYieldPi0[i]->Clone(Form("CorrectedYieldPi0Scaled_%s", triggerName[i].Data()));
        histoCorrectedYieldPi0Scaled[i]->Sumw2();
        histoCorrectedYieldPi0Scaled[i]->Scale(1./triggRejecFac[i][trigSteps[i][0]]);
        cout << trigSteps[i][0] << "\t" << trigSteps[i][1] << "\t" << trigSteps[i][2] << endl;
        if (trigSteps[i][1]!= trigSteps[i][0]){
            cout << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << endl;
            histoCorrectedYieldPi0Scaled[i]->Scale(1./triggRejecFac[trigSteps[i][0]][trigSteps[i][1]]);
        }
        if (trigSteps[i][2]!= trigSteps[i][1]){
            cout << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << "\t"<< triggRejecFac[trigSteps[i][1]][trigSteps[i][2]] << endl;
            histoCorrectedYieldPi0Scaled[i]->Scale(1./triggRejecFac[trigSteps[i][1]][trigSteps[i][2]]);
        }
    }
    
    // initialize all vectors for sytstematics and general creation of graphs, we can have a maximu of 100 data points at the moment
    Double_t xValueFinalPi0                                 [100];
    Double_t xErrorLowFinalPi0                              [100];
    Double_t xErrorHighFinalPi0                             [100];
    Double_t yValueFinalPi0                                 [100];
    Double_t yErrorLowFinalPi0                              [100];
    Double_t yErrorHighFinalPi0                             [100];
    Int_t nPointFinalPi0                                         = 0;

    Double_t yErrorSysLowFinalPi0                           [100];
    Double_t yErrorSysHighFinalPi0                          [100];

    Double_t ptSysRelPi0                                    [MaxNumberOfFiles][100];
    Double_t yErrorSysLowRelPi0                             [MaxNumberOfFiles][100];
    Double_t yErrorSysHighRelPi0                            [MaxNumberOfFiles][100];
    Bool_t sysAvailPi0                                      [MaxNumberOfFiles];

    Bool_t sysAvailSinglePi0                                [MaxNumberOfFiles];
    Int_t numberBinsSysAvailSinglePi0                       [MaxNumberOfFiles];
    
    // graphs for easier reduction of measurements to desired range and systematic errors
    TGraphAsymmErrors* graphsCorrectedYieldShrunkPi0        [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldSysShrunkPi0     [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldRemoved0Pi0      [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldSysRemoved0Pi0   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphMassPi0Data                     [MaxNumberOfFiles];
    TGraphAsymmErrors* graphMassPi0MC                       [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedMassPi0Data              [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedMassPi0MC                [MaxNumberOfFiles];
    TGraphAsymmErrors* graphWidthPi0Data                    [MaxNumberOfFiles];
    TGraphAsymmErrors* graphWidthPi0MC                      [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedWidthPi0Data             [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedWidthPi0MC               [MaxNumberOfFiles];
    TGraphAsymmErrors* graphAcceptancePi0                   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedAcceptancePi0            [MaxNumberOfFiles];
    TGraphAsymmErrors* graphEfficiencyPi0                   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEfficiencyPi0            [MaxNumberOfFiles];
    TGraphAsymmErrors* graphEffTimesAccPi0                  [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEffTimesAccPi0           [MaxNumberOfFiles];
    TGraphAsymmErrors* graphPurityPi0                       [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedPurityPi0                [MaxNumberOfFiles];
        
    // definition of predefined arrays for trigger correlation filling
    TH1D*               histoStatPi0    [12]; 
    TGraphAsymmErrors*  graphSystPi0    [12];
    TH1D*               histoRelStatPi0 [12]; 
    TGraphAsymmErrors*  graphRelSystPi0 [12];

    Int_t offSetsPi0[12]        =   { 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0 };
    Int_t offSetsPi0Sys[12]     =   { 0, 0, 0, 0, 0, 0, 
                                      0, 0, 0, 0, 0, 0 };

    if(optionEnergy.CompareTo("8TeV")==0){
      if(mode == 2){
        offSetsPi0[1] = 3; //INT7
        offSetsPi0[3] = 0; //EMC7
        offSetsPi0[4] = 4; //EGA
      }else if(mode == 4){
        offSetsPi0[1] = 0; //INT7
        offSetsPi0[3] = 0; //EMC7
        offSetsPi0[4] = 4; //EGA
      }
    }
    
    // set all graphs to NULL first 
    for (Int_t j = 0; j<12; j++){
        histoStatPi0[j]                 = NULL;
        graphSystPi0[j]                 = NULL;
        histoRelStatPi0[j]              = NULL;
        graphRelSystPi0[j]              = NULL;
        graphOrderedMassPi0Data[j]      = NULL;
        graphOrderedMassPi0MC[j]        = NULL;
        graphOrderedWidthPi0Data[j]     = NULL;
        graphOrderedWidthPi0MC[j]       = NULL;
        graphOrderedAcceptancePi0[j]    = NULL;
        graphOrderedEfficiencyPi0[j]    = NULL;
        graphOrderedEffTimesAccPi0[j]   = NULL;
        graphOrderedPurityPi0[j]        = NULL;
    }    
    
    //****************************************************************************************************************
    //************* Processing of each individual trigger, reducing ranges & adding systematics **********************
    //****************************************************************************************************************
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        // read systematics, if fileName is set to "bla" no action has to be performed and systematics will be disabled for the rest of the analysis
        if (sysFilePi0[i].CompareTo("bla") != 0){
            sysAvailPi0[i]                = kTRUE;
            ifstream  fileSysErrPi0;
            fileSysErrPi0.open(sysFilePi0[i].Data(),ios_base::in);
            cout << sysFilePi0[i].Data() << endl;
            Int_t counter = 0;
            cout << "reading sys file summed" << endl;
            while(!fileSysErrPi0.eof() && counter < 100){
                Double_t garbage = 0;
                fileSysErrPi0 >>ptSysRelPi0[i][counter] >> yErrorSysLowRelPi0[i][counter] >> yErrorSysHighRelPi0[i][counter]>>    garbage >> garbage;
                cout << counter << "\t"<< ptSysRelPi0[i][counter]<< "\t"  << yErrorSysLowRelPi0[i][counter] << "\t"  <<yErrorSysHighRelPi0[i][counter] << "\t"  << endl;;
                counter++;
            }
            fileSysErrPi0.close();
         // read in detailed systematics
            string sysFilePi0Det = sysFilePi0[i].Data();

            if(!replace(sysFilePi0Det, "Averaged", "AveragedSingle")){
                cout << "WARNING: could not find detailed systematics file " << sysFilePi0Det << ", skipping... " << endl; sysAvailSinglePi0[i] = kFALSE; 
                continue;
            }
            ifstream fileSysErrDetailedPi0;
            fileSysErrDetailedPi0.open(sysFilePi0Det,ios_base::in);
            if(fileSysErrDetailedPi0.is_open()) 
                sysAvailSinglePi0[i] = kTRUE;
            else{
                sysAvailSinglePi0[i] = kFALSE; 
                cout << "No single errors were found" << endl;
            }
            
            if (sysAvailSinglePi0[i]){
                cout << sysFilePi0Det << endl;
                counter = 0;
                string line;
                Int_t counterColumn = 0;
                while (getline(fileSysErrDetailedPi0, line) && counter < 100) {
                    istringstream ss(line);
                    TString temp="";
                    counterColumn = 0;
                    while(ss && counterColumn < 100){
                        ss >> temp;
                        if( !(counter==0 && temp.CompareTo("bin")==0) && !temp.IsNull()){
                        ptSysDetail[i][counter].push_back(temp);
                        counterColumn++;
                        }
                    }
                    if(counter == 0){
                        ptSysDetail[i][counter++].push_back("TotalError");
                        counterColumn++;
                    }else counter++;
                }
                numberBinsSysAvailSinglePi0[i] = counter;
                fileSysErrDetailedPi0.close();
            }   
        } else {
            sysAvailPi0[i]             = kFALSE;
            sysAvailSinglePi0[i]       = kFALSE;
        }
        cout << sysAvailPi0[i] << "\t" << sysAvailSinglePi0[i] << endl;
//         continue;
        
        
        // print out input spectrum from statistical histogram
        cout << "step 0" << endl;
        for (Int_t j = 1; j< histoCorrectedYieldPi0Scaled[i]->GetNbinsX()+1; j++ ){
            cout << histoCorrectedYieldPi0Scaled[i]->GetBinCenter(j) << "\t" << histoCorrectedYieldPi0Scaled[i]->GetBinContent(j) << endl;
        }
        //*******************************************************************
        //**** create graphs from original histograms for every quantity ****
        //*******************************************************************
        cout << "step 1" << endl;
        // create correct yield graphs 
        graphsCorrectedYieldShrunkPi0[i]        = new TGraphAsymmErrors(histoCorrectedYieldPi0Scaled[i]);
        graphsCorrectedYieldRemoved0Pi0[i]      = new TGraphAsymmErrors(histoCorrectedYieldPi0Scaled[i]);
        graphsCorrectedYieldSysShrunkPi0[i]     = new TGraphAsymmErrors(histoCorrectedYieldPi0Scaled[i]);
        graphsCorrectedYieldSysRemoved0Pi0[i]   = new TGraphAsymmErrors(histoCorrectedYieldPi0Scaled[i]);
        histoCorrectedYieldPi0ScaledMasked[i]   = (TH1D*)histoCorrectedYieldPi0Scaled[i]->Clone(Form("Pi0_ScaledMasked_%s",triggerName[i].Data()));
        
        // create supporting figure graphs
        if (mode != 10){ // these are not available for merged cluster analysis 
            graphMassPi0Data[i]                 = new TGraphAsymmErrors(histoMassPi0Data[i]);
            graphMassPi0MC[i]                   = new TGraphAsymmErrors(histoMassPi0MC[i]);
            graphWidthPi0Data[i]                = new TGraphAsymmErrors(histoWidthPi0Data[i]);
            graphWidthPi0MC[i]                  = new TGraphAsymmErrors(histoWidthPi0MC[i]);
        }
        if (mode == 10){ // these are not available for di-photon analysis
            graphPurityPi0[i]                   = new TGraphAsymmErrors(histoPurityPi0[i]);
        }    
        graphAcceptancePi0[i]                   = new TGraphAsymmErrors(histoAcceptancePi0[i]);
        if ( triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1")){
            graphEfficiencyPi0[i]               = new TGraphAsymmErrors(histoEfficiencyPi0[i]);
        } else if (mode == 10) {
            graphEfficiencyPi0[i]               = new TGraphAsymmErrors(histoEfficiencyPi0[i]);
            
        } else { // only possible if same file had been processed with pure MB cut on the same MC
            if (enableTriggerEffPi0[i]){ // only possible if same file had been processed with pure MB cut on the same MC
                graphEfficiencyPi0[i]           = new TGraphAsymmErrors(histoEffBasePi0[i]);
            } else {
                graphEfficiencyPi0[i]           = NULL;
            }    
        }
        graphEffTimesAccPi0[i]                  = new TGraphAsymmErrors(histoEffTimesAccPi0[i]);
        
        // remove 0 bins at beginning
        Int_t binsToMask = 1;
        while (histoCorrectedYieldPi0ScaledMasked[i]->GetBinCenter(binsToMask) < ptFromSpecPi0[i][0] ){
            histoCorrectedYieldPi0ScaledMasked[i]->SetBinContent(binsToMask,0.);
            histoCorrectedYieldPi0ScaledMasked[i]->SetBinError(binsToMask,0.);
            binsToMask++;
        }    
        while (graphAcceptancePi0[i]->GetX()[0] < ptFromSpecPi0[i][0] ){
            if (mode != 10){   
                graphMassPi0Data[i]->RemovePoint(0);
                graphMassPi0MC[i]->RemovePoint(0);
                graphWidthPi0Data[i]->RemovePoint(0);
                graphWidthPi0MC[i]->RemovePoint(0);
            }    
            if (mode == 10){
                graphPurityPi0[i]->RemovePoint(0);
            }    
            graphAcceptancePi0[i]->RemovePoint(0);
            graphEffTimesAccPi0[i]->RemovePoint(0);
            if (enableTriggerEffPi0[i] || (triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1"))){
                graphEfficiencyPi0[i]->RemovePoint(0);
            } else if (mode == 10) {
                graphEfficiencyPi0[i]->RemovePoint(0);
            }    
        }
        
        // check if trigger is supposed to be used for combination, otherwise put all graphs to NULL
        if ( maskedFullyPi0[i] ){
            graphMassPi0Data[i]     = NULL;
            graphMassPi0MC[i]       = NULL;
            graphWidthPi0Data[i]    = NULL;
            graphWidthPi0MC[i]      = NULL;
            graphAcceptancePi0[i]   = NULL;
            graphEffTimesAccPi0[i]  = NULL;
            graphPurityPi0[i]       = NULL;
            graphEfficiencyPi0[i]   = NULL;
            graphsCorrectedYieldShrunkPi0[i]    = NULL;
            graphsCorrectedYieldRemoved0Pi0[i]  = NULL; 
            graphsCorrectedYieldSysShrunkPi0[i] = NULL; 
            graphsCorrectedYieldSysRemoved0Pi0[i]   = NULL; 
            sysAvailPi0[i]          = kFALSE;
            for (Int_t f = 1; f < histoCorrectedYieldPi0ScaledMasked[i]->GetNbinsX()+1; f++ ){
                histoCorrectedYieldPi0ScaledMasked[i]->SetBinContent(f,0.);
                histoCorrectedYieldPi0ScaledMasked[i]->SetBinError(f,0.);
            }
            nrOfTrigToBeCombPi0Red--;
            cout << "trigger " << triggerName[i] << " was masked" << endl;
            continue;
        // if upper boundary > -1, remove all points above    
        } else if (ptFromSpecPi0[i][1] > -1) {
            for (Int_t f = histoCorrectedYieldPi0ScaledMasked[i]->GetXaxis()->FindBin(ptFromSpecPi0[i][1]); f < histoCorrectedYieldPi0ScaledMasked[i]->GetNbinsX()+1; f++ ){
                histoCorrectedYieldPi0ScaledMasked[i]->SetBinContent(f,0.);
                histoCorrectedYieldPi0ScaledMasked[i]->SetBinError(f,0.);
            }
            while (graphAcceptancePi0[i]->GetX()[graphAcceptancePi0[i]->GetN()-1] > ptFromSpecPi0[i][1] ){
                if (mode != 10){
                    graphMassPi0Data[i]->RemovePoint(graphMassPi0Data[i]->GetN()-1);
                    graphMassPi0MC[i]->RemovePoint(graphMassPi0MC[i]->GetN()-1);
                    graphWidthPi0Data[i]->RemovePoint(graphWidthPi0Data[i]->GetN()-1);
                    graphWidthPi0MC[i]->RemovePoint(graphWidthPi0MC[i]->GetN()-1);        
                }    
                if (mode == 10){
                    graphPurityPi0[i]->RemovePoint(graphPurityPi0[i]->GetN()-1);
                }    
                graphAcceptancePi0[i]->RemovePoint(graphAcceptancePi0[i]->GetN()-1);
                graphEffTimesAccPi0[i]->RemovePoint(graphEffTimesAccPi0[i]->GetN()-1);
                if (enableTriggerEffPi0[i] || (triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1"))){
                    graphEfficiencyPi0[i]->RemovePoint(graphEfficiencyPi0[i]->GetN()-1);
                } else if ( mode == 10) {
                    graphEfficiencyPi0[i]->RemovePoint(graphEfficiencyPi0[i]->GetN()-1);
                }    
            }    
        }    
        // if no points are left in graph put graph to NULL  
        if (graphAcceptancePi0[i]->GetN() == 0){
            if (mode != 10){
                graphMassPi0Data[i]     = NULL;
                graphMassPi0MC[i]       = NULL;
                graphWidthPi0Data[i]    = NULL;
                graphWidthPi0MC[i]      = NULL;
            }    
            graphAcceptancePi0[i]       = NULL;
            graphEffTimesAccPi0[i]      = NULL;
            graphEfficiencyPi0[i]       = NULL;
            if (mode == 10){
                graphPurityPi0[i]       = NULL;
            }    
        }
        
        // Remove 0 points at beginning for graphs
//         if (graphsCorrectedYieldShrunkPi0[i])graphsCorrectedYieldShrunkPi0[i]->Print();
        cout << "step 2" << endl;
        while (graphsCorrectedYieldShrunkPi0[i]->GetY()[0] == 0) graphsCorrectedYieldShrunkPi0[i]->RemovePoint(0);
//         if (graphsCorrectedYieldShrunkPi0[i]) graphsCorrectedYieldShrunkPi0[i]->Print();
        while (graphsCorrectedYieldRemoved0Pi0[i]->GetY()[0] == 0) graphsCorrectedYieldRemoved0Pi0[i]->RemovePoint(0);
//         if (graphsCorrectedYieldRemoved0Pi0[i]) graphsCorrectedYieldRemoved0Pi0[i]->Print();
        cout << "sys shrunk" << endl;
        while (graphsCorrectedYieldSysShrunkPi0[i]->GetY()[0] == 0) graphsCorrectedYieldSysShrunkPi0[i]->RemovePoint(0);
//         if (graphsCorrectedYieldSysShrunkPi0[i])graphsCorrectedYieldSysShrunkPi0[i]->Print();
        cout << "sys shrunk 2" << endl;
        while (graphsCorrectedYieldSysRemoved0Pi0[i]->GetY()[0] == 0) graphsCorrectedYieldSysRemoved0Pi0[i]->RemovePoint(0);
//         if (graphsCorrectedYieldSysRemoved0Pi0[i])graphsCorrectedYieldSysRemoved0Pi0[i]->Print();
        
        // put systematics on graphs
        if (graphsCorrectedYieldSysRemoved0Pi0[i]){
            for (Int_t j = 0; j< graphsCorrectedYieldSysRemoved0Pi0[i]->GetN(); j++){
                if (sysAvailPi0[i]){
                    Int_t counter = 0;
                    while(counter < 100 && TMath::Abs(graphsCorrectedYieldSysRemoved0Pi0[i]->GetX()[j] - ptSysRelPi0[i][counter])> 0.001) counter++;
                    if (counter < 100){
                        cout << ptSysRelPi0[i][counter]<< "\t found it" << endl;
                        Double_t yErrorSysLowDummy = TMath::Abs(yErrorSysLowRelPi0[i][counter]/100*graphsCorrectedYieldSysRemoved0Pi0[i]->GetY()[j]);
                        Double_t yErrorSysHighDummy = yErrorSysHighRelPi0[i][counter]/100*graphsCorrectedYieldSysRemoved0Pi0[i]->GetY()[j];
                        graphsCorrectedYieldSysRemoved0Pi0[i]->SetPointEYlow(j,yErrorSysLowDummy);
                        graphsCorrectedYieldSysRemoved0Pi0[i]->SetPointEYhigh(j,yErrorSysHighDummy);
                    } else {
                        graphsCorrectedYieldSysRemoved0Pi0[i]->SetPointEYlow(j,0);
                        graphsCorrectedYieldSysRemoved0Pi0[i]->SetPointEYhigh(j,0);
                    }
                } else {
                    graphsCorrectedYieldSysRemoved0Pi0[i]->SetPointEYlow(j,0);
                    graphsCorrectedYieldSysRemoved0Pi0[i]->SetPointEYhigh(j,0);
                    averagedPi0 = kFALSE;
                }
            }
        }
      
        cout << "step 3" << endl;
        // remove points at beginning according to ranges set for individual triggers
        while(graphsCorrectedYieldShrunkPi0[i]->GetX()[0] < ptFromSpecPi0[i][0])
            graphsCorrectedYieldShrunkPi0[i]->RemovePoint(0);
        while(graphsCorrectedYieldSysShrunkPi0[i]->GetX()[0] < ptFromSpecPi0[i][0])
            graphsCorrectedYieldSysShrunkPi0[i]->RemovePoint(0);
        // remove points at the end according to ranges set for individual triggers
        while (graphsCorrectedYieldShrunkPi0[i]->GetX()[graphsCorrectedYieldShrunkPi0[i]->GetN()-1] > ptFromSpecPi0[i][1])
            graphsCorrectedYieldShrunkPi0[i]->RemovePoint(graphsCorrectedYieldShrunkPi0[i]->GetN()-1);
//         graphsCorrectedYieldShrunkPi0[i]->Print();
        while (graphsCorrectedYieldSysShrunkPi0[i]->GetX()[graphsCorrectedYieldSysShrunkPi0[i]->GetN()-1] > ptFromSpecPi0[i][1]) 
            graphsCorrectedYieldSysShrunkPi0[i]->RemovePoint(graphsCorrectedYieldSysShrunkPi0[i]->GetN()-1);

        // put systematics on shrunk graphs
        for (Int_t j = 0; j< graphsCorrectedYieldShrunkPi0[i]->GetN(); j++){
            xValueFinalPi0[nPointFinalPi0]      = graphsCorrectedYieldShrunkPi0[i]->GetX()[j];
            xErrorHighFinalPi0[nPointFinalPi0]  = graphsCorrectedYieldShrunkPi0[i]->GetEXhigh()[j];
            xErrorLowFinalPi0[nPointFinalPi0]   = graphsCorrectedYieldShrunkPi0[i]->GetEXlow()[j];
            yValueFinalPi0[nPointFinalPi0]      = graphsCorrectedYieldShrunkPi0[i]->GetY()[j];
            yErrorHighFinalPi0[nPointFinalPi0]  = graphsCorrectedYieldShrunkPi0[i]->GetEYhigh()[j];
            yErrorLowFinalPi0[nPointFinalPi0]   = graphsCorrectedYieldShrunkPi0[i]->GetEYlow()[j];
            if (sysAvailPi0[i]){
                Int_t counter = 0;
                while(counter < 100 && TMath::Abs(xValueFinalPi0[nPointFinalPi0] - ptSysRelPi0[i][counter])> 0.001) counter++;
                if (counter < 100){
                    cout << ptSysRelPi0[i][counter]<< "\t found it" << endl;
                    yErrorSysLowFinalPi0[nPointFinalPi0] = TMath::Abs(yErrorSysLowRelPi0[i][counter]/100*graphsCorrectedYieldShrunkPi0[i]->GetY()[j]);
                    yErrorSysHighFinalPi0[nPointFinalPi0] = yErrorSysHighRelPi0[i][counter]/100*graphsCorrectedYieldShrunkPi0[i]->GetY()[j];
                    
                } else {
                    yErrorSysLowFinalPi0[nPointFinalPi0] = 0;
                    yErrorSysHighFinalPi0[nPointFinalPi0] = 0;
                    
                }
            } else {
                yErrorSysLowFinalPi0[nPointFinalPi0] = 0;
                yErrorSysHighFinalPi0[nPointFinalPi0] = 0;
            }
            graphsCorrectedYieldSysShrunkPi0[i]->SetPointEYlow(j,yErrorSysLowFinalPi0[nPointFinalPi0]);
            graphsCorrectedYieldSysShrunkPi0[i]->SetPointEYhigh(j,yErrorSysHighFinalPi0[nPointFinalPi0]);
            nPointFinalPi0++;
        }
        
        // Set correct trigger order for combination function
        if ((triggerName[i].CompareTo("MB") == 0 || triggerName[i].CompareTo("INT1") == 0  || triggerName[i].CompareTo("MB_NLM2") == 0 || triggerName[i].CompareTo("INT1_NLM2") == 0 ) 
            && graphsCorrectedYieldShrunkPi0[i]){
            cout << "filling MB trigger" << endl;
            histoStatPi0[0]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[0]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[0]    = histoStatPi0[0]->GetXaxis()->FindBin(graphSystPi0[0]->GetX()[0])-1;            
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[0]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[0]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[0]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[0]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[0]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[0]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[0]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[0]        = graphPurityPi0[i];
        } else if ((triggerName[i].CompareTo("INT7") == 0 || triggerName[i].CompareTo("INT7_NLM2") == 0) && graphsCorrectedYieldShrunkPi0[i]){   
            cout << "filling INT7 trigger" << endl;
            histoStatPi0[1]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[1]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[1]    = histoStatPi0[1]->GetXaxis()->FindBin(graphSystPi0[1]->GetX()[0])-1;
            if(optionEnergy.CompareTo("8TeV")==0 && mode == 2) offSetsPi0Sys[1]+=3;
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[1]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[1]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[1]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[1]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[1]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[1]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[1]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[1]        = graphPurityPi0[i];
        } else if ((triggerName[i].CompareTo("EMC1") == 0 || triggerName[i].CompareTo("EMC1_NLM2") == 0 ) && graphsCorrectedYieldShrunkPi0[i]){   
            cout << "filling EMC1 trigger" << endl;
            histoStatPi0[2]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[2]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[2]    = histoStatPi0[2]->GetXaxis()->FindBin(graphSystPi0[2]->GetX()[0])-1;
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[2]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[2]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[2]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[2]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[2]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[2]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[2]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[2]        = graphPurityPi0[i];
        } else if ((triggerName[i].CompareTo("EMC7") == 0 || triggerName[i].CompareTo("EMC7_NLM2") == 0 ) && graphsCorrectedYieldShrunkPi0[i]){   
            cout << "filling EMC7 trigger" << endl;
            histoStatPi0[3]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[3]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[3]    = histoStatPi0[3]->GetXaxis()->FindBin(graphSystPi0[3]->GetX()[0])-1;
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[3]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[3]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[3]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[3]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[3]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[3]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[3]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[3]        = graphPurityPi0[i];
        } else if ((triggerName[i].CompareTo("EG2") == 0 || triggerName[i].CompareTo("EG2_NLM2") == 0 ||  triggerName[i].CompareTo("EGA") == 0) && graphsCorrectedYieldShrunkPi0[i]){
            cout << Form("filling %s trigger",strEG2_A.Data()) << endl;
            histoStatPi0[4]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[4]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[4]    = histoStatPi0[4]->GetXaxis()->FindBin(graphSystPi0[4]->GetX()[0])-1;
            if(optionEnergy.CompareTo("8TeV")==0 && (mode == 4 || mode == 2)) offSetsPi0Sys[4]+=4;

            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[4]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[4]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[4]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[4]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[4]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[4]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[4]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[4]        = graphPurityPi0[i];
        } else if ((triggerName[i].CompareTo("EG1") == 0 || triggerName[i].CompareTo("EG1_NLM2") == 0 ) && graphsCorrectedYieldShrunkPi0[i]){   
            cout << "filling EG1 trigger" << endl;
            histoStatPi0[5]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[5]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[5]    = histoStatPi0[5]->GetXaxis()->FindBin(graphSystPi0[5]->GetX()[0])-1;
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[5]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[5]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[5]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[5]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[5]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[5]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[5]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[5]        = graphPurityPi0[i];
        } else if ((triggerName[i].CompareTo("MB_NLM1") == 0 || triggerName[i].CompareTo("INT1_NLM1") == 0  ) && graphsCorrectedYieldShrunkPi0[i]){
            cout << "filling MB trigger NLM1" << endl;
            histoStatPi0[6]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[6]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[6]    = histoStatPi0[6]->GetXaxis()->FindBin(graphSystPi0[6]->GetX()[0])-1;            
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[6]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[6]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[6]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[6]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[6]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[6]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[6]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[6]        = graphPurityPi0[i];
        } else if (triggerName[i].CompareTo("INT7_NLM1") == 0 && graphsCorrectedYieldShrunkPi0[i]){   
            cout << "filling INT7 trigger NLM1" << endl;
            histoStatPi0[7]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[7]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[7]    = histoStatPi0[7]->GetXaxis()->FindBin(graphSystPi0[7]->GetX()[0])-1;
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[7]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[7]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[7]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[7]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[7]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[7]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[7]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[7]        = graphPurityPi0[i];
        } else if (triggerName[i].CompareTo("EMC1_NLM1") == 0  && graphsCorrectedYieldShrunkPi0[i]){   
            cout << "filling EMC1 trigger NLM1" << endl;
            histoStatPi0[8]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[8]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[8]    = histoStatPi0[8]->GetXaxis()->FindBin(graphSystPi0[8]->GetX()[0])-1;
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[8]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[8]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[8]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[8]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[8]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[8]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[8]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[8]        = graphPurityPi0[i];
        } else if (triggerName[i].CompareTo("EMC7_NLM1") == 0  && graphsCorrectedYieldShrunkPi0[i]){   
            cout << "filling EMC7 trigger NLM1" << endl;
            histoStatPi0[9]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[9]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[9]    = histoStatPi0[9]->GetXaxis()->FindBin(graphSystPi0[9]->GetX()[0])-1;
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[9]      = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[9]        = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[9]     = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[9]       = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[9]    = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[9]    = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[9]   = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[9]        = graphPurityPi0[i];
        } else if (triggerName[i].CompareTo("EG2_NLM1") == 0 && graphsCorrectedYieldShrunkPi0[i]){   
            cout << "filling EG2 trigger NLM1" << endl;
            histoStatPi0[10]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[10]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[10]    = histoStatPi0[10]->GetXaxis()->FindBin(graphSystPi0[10]->GetX()[0])-1;
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[10]     = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[10]       = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[10]    = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[10]      = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[10]   = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[10]   = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[10]  = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[10]       = graphPurityPi0[i];
        } else if (triggerName[i].CompareTo("EG1_NLM1") == 0 && graphsCorrectedYieldShrunkPi0[i]){   
            cout << "filling EG1 trigger NLM1" << endl;
            histoStatPi0[11]     = histoCorrectedYieldPi0ScaledMasked[i];
            graphSystPi0[11]     = graphsCorrectedYieldSysShrunkPi0[i];
            offSetsPi0Sys[11]    = histoStatPi0[11]->GetXaxis()->FindBin(graphSystPi0[11]->GetX()[0])-1;
            if (graphMassPi0Data[i]) 
                graphOrderedMassPi0Data[11]     = graphMassPi0Data[i];
            if (graphMassPi0MC[i]) 
                graphOrderedMassPi0MC[11]       = graphMassPi0MC[i];
            if (graphWidthPi0Data[i]) 
                graphOrderedWidthPi0Data[11]    = graphWidthPi0Data[i];
            if (graphWidthPi0MC[i]) 
                graphOrderedWidthPi0MC[11]      = graphWidthPi0MC[i];
            if (graphAcceptancePi0[i])
                graphOrderedAcceptancePi0[11]   = graphAcceptancePi0[i];
            if (graphEfficiencyPi0[i])
                graphOrderedEfficiencyPi0[11]   = graphEfficiencyPi0[i];
            if (graphEffTimesAccPi0[i])
                graphOrderedEffTimesAccPi0[11]  = graphEffTimesAccPi0[i];
            if (graphPurityPi0[i])
                graphOrderedPurityPi0[11]       = graphPurityPi0[i];
        }
    }

//     return;
    
    // create weighted graphs for spectra and supporting graphs
    TString nameWeightsLogFilePi0 =     Form("%s/weightsPi0_%s.dat",outputDir.Data(),isMC.Data());
    TGraphAsymmErrors* graphCorrectedYieldWeightedAveragePi0Stat    = NULL;
    TGraphAsymmErrors* graphCorrectedYieldWeightedAveragePi0Sys     = NULL; 
    TGraphAsymmErrors* graphCorrectedYieldWeightedAveragePi0Tot     = NULL;
    TGraphAsymmErrors* graphMassPi0DataWeighted                     = NULL;
    TGraphAsymmErrors* graphMassPi0MCWeighted                       = NULL;
    TGraphAsymmErrors* graphWidthPi0DataWeighted                    = NULL;
    TGraphAsymmErrors* graphWidthPi0MCWeighted                      = NULL;
    TGraphAsymmErrors* graphAcceptancePi0Weighted                   = NULL;
    TGraphAsymmErrors* graphEfficiencyPi0Weighted                   = NULL;
    TGraphAsymmErrors* graphEffTimesAccPi0Weighted                  = NULL;
    TGraphAsymmErrors* graphPurityPi0Weighted                       = NULL;
    
    // Calculate averaged pi0 spectrum & respective supporting graphs according to statistical and systematic errors taking correctly into account the cross correlations
    if (averagedPi0){
        cout << maxNAllowedPi0 << endl;
        // Calculate average pi0 spectrum
        graphCorrectedYieldWeightedAveragePi0Tot        = CombinePtPointsSpectraTriggerCorrMat(    histoStatPi0, graphSystPi0,  
                                                                                                   binningPi0,  maxNAllowedPi0,
                                                                                                   offSetsPi0 ,offSetsPi0Sys,
                                                                                                   graphCorrectedYieldWeightedAveragePi0Stat, graphCorrectedYieldWeightedAveragePi0Sys,
                                                                                                   nameWeightsLogFilePi0.Data(),
                                                                                                   mode, optionEnergy, "Pi0", v2ClusterizerMerged,
                                                                                                   fileInputCorrFactors
                                                                                               );
    //return;
        // preparations for weight readout
        Double_t xValuesReadPi0[100];
        Double_t weightsReadPi0[12][100];
        Int_t availableMeasPi0[12]       = {-1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1};
        Int_t nMeasSetPi0               = nrOfTrigToBeCombPi0Red;
        Int_t nPtBinsReadPi0            = 0;
        
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
        ifstream fileWeightsPi0;
        fileWeightsPi0.open(nameWeightsLogFilePi0,ios_base::in);
        cout << "reading" << nameWeightsLogFilePi0 << endl;

        while(!fileWeightsPi0.eof() && nPtBinsReadPi0 < 100){
            TString garbage = "";
            if (nPtBinsReadPi0 == 0){
                fileWeightsPi0 >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
                for (Int_t i = 0; i < nMeasSetPi0; i++){
                    fileWeightsPi0 >> availableMeasPi0[i] ;
                }   
                cout << "read following measurements: "; 
                for (Int_t i = 0; i < nMeasSetPi0; i++){
                    cout << availableMeasPi0[i] << "\t" ;
                }   
                cout << endl;
            } else {
                fileWeightsPi0 >> xValuesReadPi0[nPtBinsReadPi0-1];
                for (Int_t i = 0; i < nMeasSetPi0; i++){
                    fileWeightsPi0 >> weightsReadPi0[availableMeasPi0[i]][nPtBinsReadPi0-1] ;
                }   
                cout << "read: "<<  nPtBinsReadPi0 << "\t"<< xValuesReadPi0[nPtBinsReadPi0-1] << "\t" ;
                for (Int_t i = 0; i < nMeasSetPi0; i++){
                    cout << weightsReadPi0[availableMeasPi0[i]][nPtBinsReadPi0-1] << "\t";
                }
                cout << endl;
            }
            nPtBinsReadPi0++;
        }
        nPtBinsReadPi0 = nPtBinsReadPi0-2 ;
        fileWeightsPi0.close();

        // creating & filling the weight graphs
        TGraph* graphWeightsPi0[12];
        for (Int_t i = 0; i < 12; i++){
            graphWeightsPi0[i]                    = NULL;
        }  
        for (Int_t i = 0; i < nMeasSetPi0; i++){
            cout << i << "\t" << availableMeasPi0[i] << endl;
            graphWeightsPi0[availableMeasPi0[i]]  = new TGraph(nPtBinsReadPi0,xValuesReadPi0,weightsReadPi0[availableMeasPi0[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsReadPi0; n++){
                if (graphWeightsPi0[availableMeasPi0[i]]->GetY()[bin] == 0) graphWeightsPi0[availableMeasPi0[i]]->RemovePoint(bin);
                else bin++;
            }
            graphWeightsPi0[availableMeasPi0[i]]->Print();
        }   

        //  **********************************************************************************************************************
        //  **************************************** Combine+write detailed Systematics ******************************************
        //  **********************************************************************************************************************

        const char *SysErrDatnameMeanSingleErr = Form("%s/SystematicErrorAveragedSingle%s_Pi0_%s.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
        fstream SysErrDatAverSingle;
        SysErrDatAverSingle.precision(4);
        cout << SysErrDatnameMeanSingleErr << endl;
        if(sysAvailSinglePi0[0]){
          SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
          for(Int_t iColumn = 0; iColumn < (Int_t)ptSysDetail[0][0].size(); iColumn++) SysErrDatAverSingle << ptSysDetail[0][0].at(iColumn) << "\t";
          SysErrDatAverSingle << endl;
          for (Int_t i = 1; i < nrOfTrigToBeComb; i++){
            if(!sysAvailSinglePi0[i]) continue;
            for(Int_t iCol = 0; iCol < (Int_t)ptSysDetail[i][0].size(); iCol++){
              if( ((TString)ptSysDetail[i][0].at(iCol)).CompareTo(((TString)ptSysDetail[0][0].at(iCol))) ){
                cout << "ERROR: Systematic error type at pos " << iCol << " does not agree for " << availableMeasPi0[i] << " & " << availableMeasPi0[0] << ", returning!" << endl;
                return;
              }
            }
          }

          for(Int_t i=0; i<nPtBinsReadPi0; i++){
            SysErrDatAverSingle << xValuesReadPi0[i] << "\t";
            Int_t nColumns = (Int_t)ptSysDetail[0][0].size();
            Double_t *errors = new Double_t[nColumns-1];
            for(Int_t iErr=0; iErr<nColumns-1; iErr++) errors[iErr] = 0;
            for(Int_t j=0; j<nrOfTrigToBeComb; j++){
              if(!sysAvailSinglePi0[j]) continue;
              Int_t pos = GetBinPosInVec(ptSysDetail[j],numberBinsSysAvailSinglePi0[j],xValuesReadPi0[i]);
              if(pos>-1){
                for(Int_t iErr=1; iErr<nColumns; iErr++) errors[iErr-1] += weightsReadPi0[availableMeasPi0[j]][i]*((TString)ptSysDetail[j][pos].at(iErr)).Atof();
              }
            }
            for(Int_t iErr=0; iErr<nColumns-1; iErr++) SysErrDatAverSingle << errors[iErr] << "\t";
            SysErrDatAverSingle << endl;
            delete[] errors;
          }
        }
        SysErrDatAverSingle.close();

        for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++){
          for(Int_t iB=0; iB<50; iB++) ptSysDetail[iR][iB].clear();
        }

        //  **********************************************************************************************************************
        //  ******************************************* Plotting weights Pi0 *****************************************************
        //  **********************************************************************************************************************
        Int_t textSizeLabelsPixel = 900*0.04;

        TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);// gives the page size
        DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
//         canvasWeights->SetLogx();
    
        TH2F * histo2DWeights;
        histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,0.,maxPtGlobalPi0,1000,-0.5,1.1);
        SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
//         histo2DWeights->GetXaxis()->SetMoreLogLabels();
//         histo2DWeights->GetXaxis()->SetLabelOffset(-0.01);
    //  histo2DWeights->GetYaxis()->SetRangeUser(-10,10);
        histo2DWeights->Draw("copy");
        
            TLegend* legendWeightsPi0 = GetAndSetLegend2(0.12, 0.14, 0.55, 0.14+(0.035*nMeasSetPi0/2*1.35), 32);
            legendWeightsPi0->SetNColumns(2);
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                DrawGammaSetMarkerTGraph(graphWeightsPi0[availableMeasPi0[i]],markerTriggWeighted[availableMeasPi0[i]], sizeTrigg[availableMeasPi0[i]], 
                                         colorTriggWeighted[availableMeasPi0[i]], colorTriggWeighted[availableMeasPi0[i]]);
                graphWeightsPi0[availableMeasPi0[i]]->Draw("p,same,e1");
                legendWeightsPi0->AddEntry(graphWeightsPi0[availableMeasPi0[i]],nameTriggerWeighted[availableMeasPi0[i]],"p");
            }   
            legendWeightsPi0->Draw();

            TLatex *labelWeightsEnergy = new TLatex(0.7,0.24,collisionSystem.Data());
            SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
            labelWeightsEnergy->SetTextFont(43);
            labelWeightsEnergy->Draw();
            TLatex *labelWeightsPi0 = new TLatex(0.7,0.20,"#pi^{0} #rightarrow #gamma#gamma");
            SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
            labelWeightsPi0->SetTextFont(43);
            labelWeightsPi0->Draw();
            TLatex *labelDetProcWeights    = new TLatex(0.7, 0.16,detectionProcess.Data());
            SetStyleTLatex( labelDetProcWeights, 0.85*textSizeLabelsPixel,4);
            labelDetProcWeights->SetTextFont(43);
            labelDetProcWeights->Draw();
            

    //      DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
            DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
            
        canvasWeights->SaveAs(Form("%s/%s_WeightsPi0Triggers.%s",outputDir.Data(), isMC.Data(), suffix.Data()));
        delete canvasWeights;

        // Calculating relative error for pi0        
        for (Int_t i = 0; i < 12; i++){
            if (histoStatPi0[i]) 
                histoRelStatPi0[i]      = CalculateRelErrUpTH1D( histoStatPi0[i], Form("relativeStatErrorPi0_%s", nameTriggerWeighted[i].Data()));
            if (graphSystPi0[i]) 
                graphRelSystPi0[i]      = CalculateRelErrUpAsymmGraph( graphSystPi0[i], Form("relativeSysErrorPi0_%s", nameTriggerWeighted[i].Data()));
        }
        
        
        TGraphAsymmErrors* graphRelErrorPi0Tot        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAveragePi0Tot, "relativeTotalErrorPi0");
        TGraphAsymmErrors* graphRelErrorPi0Stat       = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAveragePi0Stat, "relativeStatErrorPi0");
        TGraphAsymmErrors* graphRelErrorPi0Sys        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAveragePi0Sys, "relativeSysErrorPi0");

        // plot sys relative errors for individual triggers   
        TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    
        TH2F * histo2DRelSysErr;
        histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.,maxPtGlobalPi0,1000,0,60.5);
        SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelSysErr->Draw("copy");
            TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetPi0/2), 0.95, 0.92, 32);
            legendRelSysErr->SetNColumns(2);
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                cout << "plotting graph: " << availableMeasPi0[i] << "\t" <<graphRelSystPi0[availableMeasPi0[i]]->GetName() << endl;
                DrawGammaSetMarkerTGraph(graphRelSystPi0[availableMeasPi0[i]], markerTriggWeighted[availableMeasPi0[i]], sizeTrigg[availableMeasPi0[i]], 
                                         colorTriggWeighted[availableMeasPi0[i]], colorTriggWeighted[availableMeasPi0[i]]);
                graphRelSystPi0[availableMeasPi0[i]]->Draw("p,same,z");
                legendRelSysErr->AddEntry(graphRelSystPi0[availableMeasPi0[i]],nameTriggerWeighted[availableMeasPi0[i]],"p");
            }    
            legendRelSysErr->Draw();

            TLatex *labelRelErrEnergy    = new TLatex(0.15,0.89,collisionSystem.Data());
            SetStyleTLatex( labelRelErrEnergy, 0.85*textSizeLabelsPixel,4);
            labelRelErrEnergy->SetTextFont(43);
            labelRelErrEnergy->Draw();
            TLatex *labelRelErrPi0       = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
            SetStyleTLatex( labelRelErrPi0, 0.85*textSizeLabelsPixel,4);
            labelRelErrPi0->SetTextFont(43);
            labelRelErrPi0->Draw();
            TLatex *labelDetProcRelErr    = new TLatex(0.15, 0.81,detectionProcess.Data());
            SetStyleTLatex( labelDetProcRelErr, 0.85*textSizeLabelsPixel,4);
            labelDetProcRelErr->SetTextFont(43);
            labelDetProcRelErr->Draw();
            
        canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));
        
        // plot stat relative errors for individual triggers    
        TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    
        TH2F * histo2DRelStatErr;
        histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.,maxPtGlobalPi0,1000,0,60.5);
        SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelStatErr->Draw("copy");
        
            TLegend* legendRelStatErr       = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetPi0/2), 0.95, 0.92, 32);
            legendRelStatErr->SetNColumns(2);
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                 cout << "plotting graph: " << availableMeasPi0[i] << "\t" <<histoRelStatPi0[availableMeasPi0[i]]->GetName() << endl;
                if (histoRelStatPi0[availableMeasPi0[i]] && mode == 2){
                    TGraphAsymmErrors* dummyGraph = new TGraphAsymmErrors(histoRelStatPi0[availableMeasPi0[i]]);
                    dummyGraph->Print();
                    DrawGammaSetMarkerTGraph(dummyGraph, markerTriggWeighted[availableMeasPi0[i]], sizeTrigg[availableMeasPi0[i]],
                                         colorTriggWeighted[availableMeasPi0[i]], colorTriggWeighted[availableMeasPi0[i]]);
                    dummyGraph->Draw("pX,same");
                    legendRelStatErr->AddEntry(dummyGraph,nameTriggerWeighted[availableMeasPi0[i]],"p");
                    
                     for (Int_t j = 1; j < histoRelStatPi0[availableMeasPi0[i]]->GetNbinsX()+1; j++){
                        cout << j << ": " << histoRelStatPi0[availableMeasPi0[i]]->GetBinContent(j) << endl;
                     }
                } else {
                    DrawGammaSetMarker(histoRelStatPi0[availableMeasPi0[i]],markerTriggWeighted[availableMeasPi0[i]], sizeTrigg[availableMeasPi0[i]], 
                                            colorTriggWeighted[availableMeasPi0[i]], colorTriggWeighted[availableMeasPi0[i]]);
                    histoRelStatPi0[availableMeasPi0[i]]->DrawCopy("p,same,z");
                    legendRelStatErr->AddEntry(histoRelStatPi0[availableMeasPi0[i]],nameTriggerWeighted[availableMeasPi0[i]],"p");
                }
            }    
            legendRelStatErr->Draw();

            labelRelErrEnergy->Draw();
            labelRelErrPi0->Draw();
            labelDetProcRelErr->Draw();
            
        canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));

        // plot full error for final result decomposed
        TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    
        TH2F * histo2DRelTotErrPi0;
        histo2DRelTotErrPi0                 = new TH2F("histo2DRelTotErrPi0","histo2DRelTotErrPi0",11000,0.,maxPtGlobalPi0,1000,0,40.5);
        SetStyleHistoTH2ForGraphs(histo2DRelTotErrPi0, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelTotErrPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphRelErrorPi0Tot, 20, 1.5, kBlack , kBlack);
            graphRelErrorPi0Tot->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphRelErrorPi0Stat, 24, 1.5, kGray+2 , kGray+2);
            graphRelErrorPi0Stat->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphRelErrorPi0Sys, 24, 1.5, kGray+1 , kGray+1);
            graphRelErrorPi0Sys->SetLineStyle(7);
            graphRelErrorPi0Sys->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr2       = GetAndSetLegend2(0.72, 0.92-(0.035*3), 0.9, 0.92, 32);
            legendRelTotErr2->AddEntry(graphRelErrorPi0Tot,"tot","p");
            legendRelTotErr2->AddEntry(graphRelErrorPi0Stat,"stat","l");
            legendRelTotErr2->AddEntry(graphRelErrorPi0Sys,"sys","l");
            legendRelTotErr2->Draw();

            labelRelErrEnergy->Draw();
            labelRelErrPi0->Draw();
            labelDetProcRelErr->Draw();
            
        canvasRelTotErr->SaveAs(Form("%s/Pi0_RelErrorsFulldecomp.%s",outputDir.Data(),suffix.Data()));
        
        //return;
        // Calculation of averaged supporting plots with weights from spectra
        if (mode != 10){
            graphMassPi0DataWeighted                    = CalculateWeightedQuantity(    graphOrderedMassPi0Data, 
                                                                                        graphWeightsPi0,
                                                                                        binningPi0,  maxNAllowedPi0,
                                                                                        MaxNumberOfFiles
                                                                                   );
            if (!graphMassPi0DataWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            graphMassPi0MCWeighted                      = CalculateWeightedQuantity(    graphOrderedMassPi0MC, 
                                                                                        graphWeightsPi0,
                                                                                        binningPi0,  maxNAllowedPi0,
                                                                                        MaxNumberOfFiles
                                                                                   );
            if (!graphMassPi0MCWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            
            graphWidthPi0DataWeighted                   = CalculateWeightedQuantity(    graphOrderedWidthPi0Data, 
                                                                                        graphWeightsPi0,
                                                                                        binningPi0,  maxNAllowedPi0,
                                                                                        MaxNumberOfFiles
                                                                                   );
            if (!graphWidthPi0DataWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }

            graphWidthPi0MCWeighted                     = CalculateWeightedQuantity(    graphOrderedWidthPi0MC, 
                                                                                        graphWeightsPi0,
                                                                                        binningPi0,  maxNAllowedPi0,
                                                                                        MaxNumberOfFiles
                                                                                   );
            if (!graphWidthPi0MCWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }            
        }    
        
        cout << "weighting Pi0 acceptance" << endl;
        graphAcceptancePi0Weighted                      = CalculateWeightedQuantity(    graphOrderedAcceptancePi0, 
                                                                                        graphWeightsPi0,
                                                                                        binningPi0,  maxNAllowedPi0,
                                                                                        MaxNumberOfFiles
                                                                                    );
        if (!graphAcceptancePi0Weighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }

        cout << "weighting Pi0 efficiency" << endl;
        graphEfficiencyPi0Weighted                      = CalculateWeightedQuantity(    graphOrderedEfficiencyPi0, 
                                                                                        graphWeightsPi0,
                                                                                        binningPi0,  maxNAllowedPi0,
                                                                                        MaxNumberOfFiles
                                                                                    );
        if (!graphEfficiencyPi0Weighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }
        cout << "weighting Pi0 efficiency x acceptance" << endl;
        graphEffTimesAccPi0Weighted                     = CalculateWeightedQuantity(    graphOrderedEffTimesAccPi0, 
                                                                                        graphWeightsPi0,
                                                                                        binningPi0,  maxNAllowedPi0,
                                                                                        MaxNumberOfFiles
                                                                                    );
        
        if (!graphEffTimesAccPi0Weighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }
        if (mode == 10){
            graphPurityPi0Weighted                      = CalculateWeightedQuantity(    graphOrderedPurityPi0, 
                                                                                        graphWeightsPi0,
                                                                                        binningPi0,  maxNAllowedPi0,
                                                                                        MaxNumberOfFiles
                                                                                    );
            if (!graphPurityPi0Weighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            
        }
        
        // remove points in spectrum which should have been masked
        if (mode != 10){
            if (graphMassPi0DataWeighted)
                while (graphMassPi0DataWeighted->GetY()[0] == -10000 )   graphMassPi0DataWeighted->RemovePoint(0);
            else 
                cout << "I don't have a weighted Pi0 mass data graph" << endl;
            if (graphMassPi0MCWeighted)
                while (graphMassPi0MCWeighted->GetY()[0] == -10000)     graphMassPi0MCWeighted->RemovePoint(0);
            else 
                cout << "I don't have a weighted Pi0 mass MC graph" << endl;
            if (graphWidthPi0DataWeighted)
                while (graphWidthPi0DataWeighted->GetY()[0] == -10000)  graphWidthPi0DataWeighted->RemovePoint(0);
            else
                cout << "I don't have a weighted Pi0 width data graph" << endl;
            if (graphWidthPi0MCWeighted)
                while (graphWidthPi0MCWeighted->GetY()[0] == -10000)    graphWidthPi0MCWeighted->RemovePoint(0);
            else 
                cout << "I don't have a weighted Pi0 width MC graph" << endl;
        }        
        if (mode == 10){
            if (graphPurityPi0Weighted)
                while (graphPurityPi0Weighted->GetY()[0] == -10000)     graphPurityPi0Weighted->RemovePoint(0);
            else 
                cout << "I don't have a weighted Pi0 purity graph" << endl;
                
        }    
        if (graphAcceptancePi0Weighted)
            while (graphAcceptancePi0Weighted->GetY()[0] == -10000)     graphAcceptancePi0Weighted->RemovePoint(0);
        else 
            cout << "I don't have a weighted Pi0 acceptance graph" << endl;
        
        if (graphEfficiencyPi0Weighted)
            while (graphEfficiencyPi0Weighted->GetY()[0] == -10000)     graphEfficiencyPi0Weighted->RemovePoint(0);
        else 
            cout << "I don't have a weighted Pi0 efficiency graph" << endl;

        if (graphEffTimesAccPi0Weighted)    
            while (graphEffTimesAccPi0Weighted->GetY()[0] == -10000)    graphEffTimesAccPi0Weighted->RemovePoint(0);
        else 
            cout << "I don't have a weighted Pi0 acceptance x efficiency graph" << endl;

    // if averaging wasn't enabled pick values according to predefined ranges ("cherry picking points")    
    } else {
       graphCorrectedYieldWeightedAveragePi0Stat        = new TGraphAsymmErrors(nPointFinalPi0, xValueFinalPi0, yValueFinalPi0, 
                                                                                xErrorLowFinalPi0, xErrorHighFinalPi0,yErrorLowFinalPi0, yErrorHighFinalPi0);
       graphCorrectedYieldWeightedAveragePi0Sys         = new TGraphAsymmErrors(nPointFinalPi0, xValueFinalPi0, yValueFinalPi0, 
                                                                                xErrorLowFinalPi0, xErrorHighFinalPi0,yErrorSysLowFinalPi0, yErrorSysHighFinalPi0);
    }    
    // print final graphs
    cout << "stat pi0" << endl; 
    graphCorrectedYieldWeightedAveragePi0Stat->Print();
    cout << "sys pi0" << endl;
    if (graphCorrectedYieldWeightedAveragePi0Sys) graphCorrectedYieldWeightedAveragePi0Sys->Print();

    if (graphCorrectedYieldWeightedAveragePi0Stat){
        while (graphCorrectedYieldWeightedAveragePi0Stat->GetX()[0]< minPtGlobalPi0){
            graphCorrectedYieldWeightedAveragePi0Stat->RemovePoint(0);
        }    
    }
    
    if (graphCorrectedYieldWeightedAveragePi0Sys){
        while (graphCorrectedYieldWeightedAveragePi0Sys->GetX()[0]< minPtGlobalPi0){
            graphCorrectedYieldWeightedAveragePi0Sys->RemovePoint(0);
        }    
    }

    if (mode != 10){
        //***************************************************************************************************************
        //************************************Plotting Mass Pi0 reduced range  ******************************************
        //***************************************************************************************************************
        canvasMass->cd();
        histo2DMassPi0->DrawCopy(); 
        
        legendMassRedPi0->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0)
             ){
                if (graphMassPi0Data[i] && !maskedFullyPi0[i]) {
                    DrawGammaSetMarkerTGraphAsym(graphMassPi0Data[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    graphMassPi0Data[i]->Draw("p,e1,same"); 
                    legendMassRedPi0->AddEntry(graphMassPi0Data[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                }
                if (graphMassPi0MC[i] && !maskedFullyPi0[i]){
                    DrawGammaSetMarkerTGraphAsym(graphMassPi0MC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    graphMassPi0MC[i]->Draw("p,e1,same"); 
                    legendMassRedPi0->AddEntry(graphMassPi0MC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                }    
            }    
        }        
        legendMassRedPi0->Draw();
        labelEnergyMass->Draw();
        labelPi0Mass->Draw();
        labelDetProcMass->Draw();

        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/Pi0_%s_Mass1_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        if ( optionEnergy.CompareTo("2.76TeV")==0 ){
            canvasMass->cd();
            histo2DMassPi0->DrawCopy();    
            legendMassRedPi02->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if( i==1 || i==3 || i==4 || i==5){
                    if (graphMassPi0Data[i] && !maskedFullyPi0[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphMassPi0Data[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                        graphMassPi0Data[i]->Draw("p,e1,same"); 
                        legendMassRedPi02->AddEntry(graphMassPi0Data[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                    }
                    if (graphMassPi0MC[i] && !maskedFullyPi0[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphMassPi0MC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                        graphMassPi0MC[i]->Draw("p,e1,same"); 
                        legendMassRedPi02->AddEntry(graphMassPi0MC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                    }    
                }    
            }
            legendMassRedPi02->Draw();
            labelEnergyMass->Draw();
            labelPi0Mass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Pi0_%s_Mass2_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        //***************************************************************************************************************
        //********************************* Pi0 Mass weighted ***********************************************************
        //***************************************************************************************************************    
        if (graphMassPi0DataWeighted || graphMassPi0MCWeighted){
            canvasMass->cd();
            histo2DMassPi0->DrawCopy();    
            TLegend* legendMassPi0Weighted = GetAndSetLegend2(0.52, 0.88, 0.75, 0.88+(1.05*2*0.85*textSizeSpectra),28);
            if (graphMassPi0DataWeighted){
                DrawGammaSetMarkerTGraphAsym(graphMassPi0DataWeighted, 20, 1, kBlack, kBlack);
                graphMassPi0DataWeighted->Draw("p,e1,same");                 
                legendMassPi0Weighted->AddEntry(graphMassPi0DataWeighted, "Data", "p"); 
            }    
            if (graphMassPi0DataWeighted){
                DrawGammaSetMarkerTGraphAsym(graphMassPi0MCWeighted, 24, 1, kGray+2, kGray+2);
                graphMassPi0MCWeighted->Draw("p,e1,same");                 
                legendMassPi0Weighted->AddEntry(graphMassPi0MCWeighted, "MC", "p"); 
            }    
            legendMassPi0Weighted->Draw();
            labelEnergyMass->Draw();
            labelPi0Mass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Pi0_%s_Mass_Weighted.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        //***************************************************************************************************************
        //************************************Plotting Width Pi0 reduced range  *****************************************
        //***************************************************************************************************************
        canvasWidth->cd();
        histo2DWidthPi0->DrawCopy(); 

        legendWidthRedPi0->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0)
             ){
                if (graphWidthPi0Data[i] && !maskedFullyPi0[i]) {
                    DrawGammaSetMarkerTGraphAsym(graphWidthPi0Data[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    graphWidthPi0Data[i]->Draw("p,e1,same"); 
                    legendWidthRedPi0->AddEntry(graphWidthPi0Data[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                } 
                if (graphWidthPi0MC[i] && !maskedFullyPi0[i]){
                    DrawGammaSetMarkerTGraphAsym(graphWidthPi0MC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    graphWidthPi0MC[i]->Draw("p,e1,same"); 
                    legendWidthRedPi0->AddEntry(graphWidthPi0MC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                }    
            }    
        }
        legendWidthRedPi0->Draw();
        labelEnergyWidth->Draw();
        labelPi0Width->Draw();
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Pi0_%s_Width1_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        if (optionEnergy.CompareTo("2.76TeV")==0){
            canvasWidth->cd();
            histo2DWidthPi0->DrawCopy();    
            legendWidthRedPi02->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if( i==1 || i==3 || i==4 || i==5 ){
                    if (graphWidthPi0Data[i] && !maskedFullyPi0[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphWidthPi0Data[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                        graphWidthPi0Data[i]->Draw("p,e1,same"); 
                        legendWidthRedPi02->AddEntry(graphWidthPi0Data[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                    }    
                    if (graphWidthPi0MC[i] && !maskedFullyPi0[i]) {    
                        DrawGammaSetMarkerTGraphAsym(graphWidthPi0MC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                        graphWidthPi0MC[i]->Draw("p,e1,same"); 
                        legendWidthRedPi02->AddEntry(graphWidthPi0MC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                    }    
                }    
            }
            legendWidthRedPi02->Draw();
            labelEnergyWidth->Draw();
            labelPi0Width->Draw();
            labelDetProcWidth->Draw();

            canvasWidth->Update();
            canvasWidth->SaveAs(Form("%s/Pi0_%s_Width2_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        //***************************************************************************************************************
        //********************************* Pi0 Width weighted **********************************************************
        //***************************************************************************************************************        
        canvasWidth->cd();
        histo2DWidthPi0->DrawCopy();    
        TLegend* legendWidthPi0Weighted = GetAndSetLegend2(0.52, 0.86, 0.75, 0.86+(1.05*2*0.85*textSizeSpectra),28);
        if (graphWidthPi0DataWeighted){
            DrawGammaSetMarkerTGraphAsym(graphWidthPi0DataWeighted, 20, 1, kBlack, kBlack);
            graphWidthPi0DataWeighted->Draw("p,e1,same");                 
            legendWidthPi0Weighted->AddEntry(graphWidthPi0DataWeighted, "Data", "p"); 
        }    
        if (graphWidthPi0DataWeighted){
            DrawGammaSetMarkerTGraphAsym(graphWidthPi0MCWeighted, 24, 1, kGray+2, kGray+2);
            graphWidthPi0MCWeighted->Draw("p,e1,same");                 
            legendWidthPi0Weighted->AddEntry(graphWidthPi0MCWeighted, "MC", "p"); 
        }    
        legendWidthPi0Weighted->Draw();
        labelEnergyWidth->Draw();
        labelPi0Width->Draw();
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Pi0_%s_Width_Weighted.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
    }
    if (!enableEta) delete canvasMass;
    if (!enableEta) delete canvasWidth;
    
    //***************************************************************************************************************
    //************************************* Efficiency weighted *****************************************************
    //***************************************************************************************************************    
    if (graphEfficiencyPi0Weighted){
        DrawGammaCanvasSettings( canvasEffi, 0.09, 0.017, 0.015, 0.08);
        canvasEffi->SetLogy(1);
        canvasEffi->cd();
        histo2DEffiPi0->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphEfficiencyPi0Weighted, 20, 1, kGray+2, kGray+2);
        graphEfficiencyPi0Weighted->Draw("p,e1,same");                 

        TLatex *labelEnergyEffiWOTrigg = new TLatex(0.62, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyEffiWOTrigg, 0.85*textSizeSpectra,4);
        labelEnergyEffiWOTrigg->Draw();

        TLatex *labelPi0EffiWOTrigg = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelPi0EffiWOTrigg, 0.85*textSizeSpectra,4);
        labelPi0EffiWOTrigg->Draw();

        TLatex *labelDetProcEffiWOTrigg = new TLatex(0.62, 0.15,detectionProcess.Data());
        SetStyleTLatex( labelDetProcEffiWOTrigg, 0.85*textSizeSpectra,4);
        labelDetProcEffiWOTrigg->Draw();
        
        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/Pi0_EfficiencyW0TriggEff_Weighted.%s",outputDir.Data(),suffix.Data()));
    }
//     if (!enableEta) delete canvasEffi;

    //***************************************************************************************************************
    //***************************************** Purity weighted *****************************************************
    //***************************************************************************************************************    
    if (graphPurityPi0Weighted){
        TCanvas* canvasPurity = new TCanvas("canvasPurity","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasPurity, 0.09, 0.017, 0.015, 0.08);
        canvasPurity->SetLogy(0);

        Double_t minPurityPi0 = 0.6;
        Double_t maxPurityPi0 = 1.02;

        TH2F * histo2DPurityPi0;
        histo2DPurityPi0 = new TH2F("histo2DPurityPi0","histo2DPurityPi0",1000,0., maxPtGlobalPi0,10000,minPurityPi0, maxPurityPi0);
        SetStyleHistoTH2ForGraphs(histo2DPurityPi0, "#it{p}_{T} (GeV/#it{c})","#it{P}_{#pi^{0}}", 
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
        histo2DPurityPi0->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphPurityPi0Weighted, 20, 1, kGray+2, kGray+2);
        graphPurityPi0Weighted->Draw("p,e1,same");                 

        TLatex *labelEnergyPurityWeighted = new TLatex(0.62, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyPurityWeighted, 0.85*textSizeSpectra,4);
        labelEnergyPurityWeighted->Draw();

        TLatex *labelPi0PurityWeighted = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelPi0PurityWeighted, 0.85*textSizeSpectra,4);
        labelPi0PurityWeighted->Draw();

        TLatex *labelDetProcPurityWeighted = new TLatex(0.62, 0.15,detectionProcess.Data());
        SetStyleTLatex( labelDetProcPurityWeighted, 0.85*textSizeSpectra,4);
        labelDetProcPurityWeighted->Draw();

        
        canvasPurity->Update();
        canvasPurity->SaveAs(Form("%s/Pi0_Purity_Weighted.%s",outputDir.Data(),suffix.Data()));
        
    }    
    
    //***************************************************************************************************************
    //************************************* Acceptance weighted *****************************************************
    //***************************************************************************************************************    
    if (graphAcceptancePi0Weighted){
        DrawGammaCanvasSettings( canvasAcc, 0.1, 0.017, 0.015, 0.08);
        canvasAcc->cd();
        canvasAcc->SetLogy(0);
        histo2DAccPi0->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphAcceptancePi0Weighted, 20, 1, kGray+2, kGray+2);
        graphAcceptancePi0Weighted->Draw("p,e1,same");                 
        
        TLatex *labelEnergyAcc = new TLatex(0.62, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyAcc, 0.85*textSizeSpectra,4);
        labelEnergyAcc->Draw();

        TLatex *labelPi0Acc = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelPi0Acc, 0.85*textSizeSpectra,4);
        labelPi0Acc->Draw();

        TLatex *labelDetProcAcc = new TLatex(0.62, 0.15,detectionProcess.Data());
        SetStyleTLatex( labelDetProcAcc, 0.85*textSizeSpectra,4);
        labelDetProcAcc->Draw();

        canvasAcc->Update();
        canvasAcc->SaveAs(Form("%s/Pi0_Acceptance_weighted.%s",outputDir.Data(),suffix.Data()));
    }    
//     if (!enableEta) delete canvasAcc;

    //***************************************************************************************************************
    //************************************Plotting scaled invariant yield *****************************************
    //***************************************************************************************************************
    TCanvas* canvasCorrScaled = new TCanvas("canvasCorrScaled","",0,0,1000,1350);// gives the page size
    DrawGammaCanvasSettings( canvasCorrScaled, 0.15, 0.017, 0.015, 0.07);
    canvasCorrScaled->SetLogy();
//     canvasCorrScaled->SetGridx();
    Double_t minCorrYield       = 2e-10;
    Double_t maxCorrYield       = 1e0;
    if (mode == 10) {
        minCorrYield            = 2e-12;
        maxCorrYield            = 1e-3;
    } else if (mode == 0){
        minCorrYield            = 2e-9;
        maxCorrYield            = 1e1;        
    }    

    if(optionEnergy.CompareTo("8TeV")==0){
      if(mode == 2){
        minCorrYield       = 1e-10;
        maxCorrYield       = 1;
      }else if(mode == 4){
        minCorrYield       = 1e-8;
        maxCorrYield       = 0.2;
      }
    }
    
    TH2F * histo2DInvYieldScaled;
    histo2DInvYieldScaled = new TH2F("histo2DInvYieldScaled","histo2DInvYieldScaled",1000,0., maxPtGlobalPi0,10000,minCorrYield,maxCorrYield);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldScaled, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.55);
    histo2DInvYieldScaled->DrawCopy(); 

    Double_t factorPi0              = 10000;
    if (mode == 2) factorPi0        = 10000;
    if (mode == 4) factorPi0        = 5000;

    // two component model fit
    Double_t paramTCM[5] = {graphCorrectedYieldWeightedAveragePi0Stat->GetY()[0],0.3,graphCorrectedYieldWeightedAveragePi0Stat->GetY()[0]/factorPi0,0.8,3};
    TF1* fitInvYieldPi0 = FitObject("tcm","fitInvYieldPi0","Pi0",graphCorrectedYieldWeightedAveragePi0Stat,minPtGlobalPi0,maxPtGlobalPi0,paramTCM,"QNRMEX0+");

    // Tsallis fit
//     Double_t paramGraph[3]                  = {1000, 8., 0.13};
//     TF1* fitInvYieldPi0                     = FitObject("l","fitInvYieldPi0","Pi0",graphCorrectedYieldWeightedAveragePi0Stat,minPtGlobalPi0,maxPtGlobalPi0,paramGraph,"QNRME+");

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAveragePi0Sys, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
    graphCorrectedYieldWeightedAveragePi0Sys->Draw("p,E2,same");

    TLegend* legendScaled = GetAndSetLegend2(0.72, 0.95-(1.15*nrOfTrigToBeCombPi0Red*0.85*textSizeSpectra), 0.95, 0.95,32);
    
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){        
        if (graphsCorrectedYieldSysShrunkPi0[i]) DrawGammaSetMarkerTGraphAsym(graphsCorrectedYieldSysShrunkPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
        if (graphsCorrectedYieldSysShrunkPi0[i])DrawGammaSetMarker(histoCorrectedYieldPi0Scaled[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        if ( !maskedFullyPi0[i] )histoCorrectedYieldPi0Scaled[i]->DrawCopy("e1,same"); 
        if ( !maskedFullyPi0[i] )graphsCorrectedYieldSysShrunkPi0[i]->Draw("p,E2,same");
        if (graphsCorrectedYieldSysShrunkPi0[i] && !maskedFullyPi0[i])legendScaled->AddEntry(histoCorrectedYieldPi0Scaled[i],triggerNameLabel[i].Data(),"p");
    }
    
    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAveragePi0Stat, 24, 2, kRed , kRed, 1, kTRUE);
    graphCorrectedYieldWeightedAveragePi0Stat->Draw("p,E,same");
    legendScaled->AddEntry(graphCorrectedYieldWeightedAveragePi0Stat,"Final","p");
    legendScaled->Draw();
    
    fitInvYieldPi0->SetLineColor(kGray+2);
    fitInvYieldPi0->SetLineStyle(7);
    fitInvYieldPi0->SetLineWidth(2);
    fitInvYieldPi0->Draw("same");

    labelEnergyUnscaled->Draw();
    labelPi0Unscaled->Draw();
    labelDetProcUnscaled->Draw();
    
    canvasCorrScaled->Update();
    canvasCorrScaled->SaveAs(Form("%s/Pi0_%s_CorrectedYieldScaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    //***************************************************************************************************************
    //************************************Plotting final invariant yield ********************************************
    //***************************************************************************************************************

    histo2DInvYieldScaled->DrawCopy(); 
    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAveragePi0Sys, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
    graphCorrectedYieldWeightedAveragePi0Sys->Draw("p,E2,same");
    
    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAveragePi0Stat, 24, 2, kBlack , kBlack, 1, kTRUE);
    graphCorrectedYieldWeightedAveragePi0Stat->Draw("p,E,same");
    
    fitInvYieldPi0->Draw("same");
        
    labelEnergyUnscaled->Draw();
    labelPi0Unscaled->Draw();
    labelDetProcUnscaled->Draw();
    
    canvasCorrScaled->Update();
    canvasCorrScaled->SaveAs(Form("%s/Pi0_%s_CorrectedYieldFinal.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

    DrawGammaSetMarker(histoMCInputPi0[0],  0, 0, kBlue+2, kBlue+2);    
    histoMCInputPi0[0]->Draw("same,hist,c");
    
    labelEnergyUnscaled->Draw();
    labelPi0Unscaled->Draw();
    labelDetProcUnscaled->Draw();
    
    canvasCorrScaled->Update();
    canvasCorrScaled->SaveAs(Form("%s/Pi0_%s_CorrectedYieldFinal_WithMC.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
    
    if (!enableEta) delete canvasCorrScaled;

    //***************************************************************************************************************
    //****************************** Ratio to fit for individual spectra full range *********************************
    //***************************************************************************************************************
    TCanvas* canvasRatioSpec = new TCanvas("canvasRatioSpec","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasRatioSpec, 0.09, 0.017, 0.015, 0.08);
    canvasRatioSpec->SetLogy(0);
    
    TH2F * histo2DRatioToFitPi0;
    histo2DRatioToFitPi0 = new TH2F("histo2DRatioToFitPi0","histo2DRatioToFitPi0",1000,0., maxPtGlobalPi0,1000,0.55, 1.85);
    SetStyleHistoTH2ForGraphs(histo2DRatioToFitPi0, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
    histo2DRatioToFitPi0->DrawCopy(); 

    TH1D* histoCorrectedYieldToFitPi0[MaxNumberOfFiles];
    TGraphAsymmErrors* graphCorrectedYieldToFitPi0[MaxNumberOfFiles];
    TLegend* legendRatioSpecPi0 = GetAndSetLegend2(0.12, 0.95-(1.05*nrOfTrigToBeCombPi0Red/2*0.85*textSizeSpectra), 0.5, 0.95,28);
    legendRatioSpecPi0->SetNColumns(2);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if (graphsCorrectedYieldSysRemoved0Pi0[i] && !maskedFullyPi0[i]){ 
          histoCorrectedYieldToFitPi0[i] = CalculateHistoRatioToFit (histoCorrectedYieldPi0Scaled[i], fitInvYieldPi0); 
          DrawGammaSetMarker(histoCorrectedYieldToFitPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
          legendRatioSpecPi0->AddEntry(histoCorrectedYieldToFitPi0[i],triggerNameLabel[i].Data(),"p");
          graphCorrectedYieldToFitPi0[i] = CalculateGraphErrRatioToFit(graphsCorrectedYieldSysRemoved0Pi0[i], fitInvYieldPi0); 
          DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
          graphCorrectedYieldToFitPi0[i]->Draw("p,E2,same");
        }
        
        DrawGammaLines(0., maxPtGlobalPi0 , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0., maxPtGlobalPi0 , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0., maxPtGlobalPi0 , 0.9, 0.9,0.1, kGray, 7);

        if (graphsCorrectedYieldSysRemoved0Pi0[i]){ 
          histoCorrectedYieldToFitPi0[i]->DrawCopy("e1,same"); 
        }  
    }
    legendRatioSpecPi0->Draw();

    TLatex *labelEnergyRatio = new TLatex(0.6, 0.93, collisionSystem.Data());
    SetStyleTLatex( labelEnergyRatio, 0.85*textSizeSpectra,4);
    labelEnergyRatio->Draw();
    
    TLatex *labelPi0Ratio = new TLatex(0.6, 0.93-textSizeSpectra*0.85*1.04, "#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelPi0Ratio, 0.85*textSizeSpectra,4);
    labelPi0Ratio->Draw();
    
    TLatex *labelDetProcRatio = new TLatex(0.6, 0.93-(2*textSizeSpectra*0.85), detectionProcess.Data());
    SetStyleTLatex( labelDetProcRatio, 0.85*textSizeSpectra,4);
    labelDetProcRatio->Draw();
    
    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/Pi0_%s_RatioSpectraToFit.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

    //***************************************************************************************************************
    //****************************** Ratio to fit for individual spectra used range *********************************
    //***************************************************************************************************************

    histo2DRatioToFitPi0->DrawCopy(); 

    TGraphAsymmErrors* graphCorrectedYieldToFitPi0Used[MaxNumberOfFiles];
    TGraphAsymmErrors* graphCorrectedYieldToFitPi0SysUsed[MaxNumberOfFiles];
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if (!maskedFullyPi0[i]) {
            if (graphsCorrectedYieldShrunkPi0[i]) graphCorrectedYieldToFitPi0Used[i] = CalculateGraphErrRatioToFit (graphsCorrectedYieldShrunkPi0[i], fitInvYieldPi0); 
            if (graphsCorrectedYieldSysShrunkPi0[i]) graphCorrectedYieldToFitPi0SysUsed[i] = CalculateGraphErrRatioToFit(graphsCorrectedYieldSysShrunkPi0[i], fitInvYieldPi0); 
            
            if (graphsCorrectedYieldSysShrunkPi0[i])DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitPi0SysUsed[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
            if (graphsCorrectedYieldShrunkPi0[i])DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitPi0Used[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);

            if (graphsCorrectedYieldSysShrunkPi0[i])graphCorrectedYieldToFitPi0SysUsed[i]->Draw("p,E2,same");
            
            DrawGammaLines(0., maxPtGlobalPi0 , 1., 1.,0.1, kGray+2);
            DrawGammaLines(0., maxPtGlobalPi0 , 1.1, 1.1,0.1, kGray, 7);
            DrawGammaLines(0., maxPtGlobalPi0 , 0.9, 0.9,0.1, kGray, 7);

            if (graphsCorrectedYieldShrunkPi0[i])graphCorrectedYieldToFitPi0Used[i]->Draw("e1,same"); 
        }    
    }
    
    legendRatioSpecPi0->Draw();

    labelEnergyRatio->Draw();
    labelPi0Ratio->Draw();
    labelDetProcRatio->Draw();
    
    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/Pi0_%s_RatioSpectraToFitUsed.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

    //***************************************************************************************************************
    //****************************** Ratio to fit for final spectrum ************************************************
    //***************************************************************************************************************
    
    histo2DRatioToFitPi0->DrawCopy(); 
    
    TGraphAsymmErrors* graphCorrectedYieldFinalStatToFitPi0;
    TGraphAsymmErrors* graphCorrectedYieldFinalSysToFitPi0;
    TH1D* histoMCInputToFit;
   
    graphCorrectedYieldFinalStatToFitPi0    = CalculateGraphErrRatioToFit (graphCorrectedYieldWeightedAveragePi0Stat, fitInvYieldPi0); 
    graphCorrectedYieldFinalSysToFitPi0     = CalculateGraphErrRatioToFit(graphCorrectedYieldWeightedAveragePi0Sys, fitInvYieldPi0); 
    histoMCInputToFit                       = (TH1D*)histoMCInputPi0[0]->Clone("Pi0MCToFit");
    histoMCInputToFit                       = CalculateHistoRatioToFitNLO (histoMCInputToFit, fitInvYieldPi0, minPtGlobalPi0);
    
    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldFinalSysToFitPi0, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldFinalStatToFitPi0, 24, 2, kBlack , kBlack, 1, kTRUE);
    
    graphCorrectedYieldFinalSysToFitPi0->Draw("p,E2,same");
    
    DrawGammaLines(0., maxPtGlobalPi0 , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0., maxPtGlobalPi0 , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(0., maxPtGlobalPi0 , 0.9, 0.9,0.1, kGray, 7);

    graphCorrectedYieldFinalStatToFitPi0->Draw("p,E,same");
   
    labelEnergyRatio->Draw();
    labelPi0Ratio->Draw();
    labelDetProcRatio->Draw();
    
    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/Pi0_%s_RatioSpectraToFitFinal.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

    TH2F * histo2DRatioToFitPi02;
    histo2DRatioToFitPi02 = new TH2F("histo2DRatioToFitPi02","histo2DRatioToFitPi02",1000,0., maxPtGlobalPi0,1000,0.25, 1.85);
    SetStyleHistoTH2ForGraphs(histo2DRatioToFitPi02, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
    histo2DRatioToFitPi02->DrawCopy(); 
    
    graphCorrectedYieldFinalSysToFitPi0->Draw("p,E2,same");
    
    DrawGammaLines(0., maxPtGlobalPi0 , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0., maxPtGlobalPi0 , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(0., maxPtGlobalPi0 , 0.9, 0.9,0.1, kGray, 7);
    
    
    graphCorrectedYieldFinalStatToFitPi0->Draw("p,E,same");
    DrawGammaSetMarker(histoMCInputToFit, 20, 1.2, kBlue+1, kBlue+1);    
    histoMCInputToFit->Draw("same,pe");
    DrawGammaLines(0., maxPtGlobalPi0 , 0.6, 0.6,1, kBlue-7, 7);
    
    labelEnergyRatio->Draw();
    labelPi0Ratio->Draw();
    labelDetProcRatio->Draw();
    
    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/Pi0_%s_RatioSpectraToFitFinal_withMC.%s",outputDir.Data(),isMC.Data(), suffix.Data()));
    
    
    if (!enableEta) delete canvasRatioSpec;
    
    
    
    //***************************************************************************************************************
    //************************************Loading eta histograms ****************************************************
    //***************************************************************************************************************    
    TString FileNameCorrectedEta        [MaxNumberOfFiles];
    TFile* fileCorrectedEta             [MaxNumberOfFiles];
    TString FileNameCorrectedPi0EtaBin  [MaxNumberOfFiles];
    TFile* fileCorrectedPi0EtaBin       [MaxNumberOfFiles];
    Bool_t foundPi0EtaBinFile           [MaxNumberOfFiles];
    Bool_t doEtaToPi0                                               = kFALSE;
    Bool_t doBinShiftForEtaToPi0                                    = kFALSE;
    TString addNameBinshift                                         = "";
    
    TString FileNameUnCorrectedEta      [MaxNumberOfFiles];
    TFile* fileUnCorrectedEta           [MaxNumberOfFiles];

    TH1D*   histoCorrectedYieldEta      [MaxNumberOfFiles];
    TH1D*   histoEfficiencyEta          [MaxNumberOfFiles];
    TH1D*   histoAcceptanceEta          [MaxNumberOfFiles];
    TH1D*   histoAcceptanceEtaWOEvtWeights [MaxNumberOfFiles];
    TH1D*   histoEffTimesAccEta         [MaxNumberOfFiles];
    TH1D*   histoRawYieldEta            [MaxNumberOfFiles];
    TH1D*   histoMCInputEta             [MaxNumberOfFiles];
    TH1D*   histoMCInputEtaPi0          [MaxNumberOfFiles];
    TH1D*   histoCorrectedYieldPi0EtaBin[MaxNumberOfFiles];
    TH1D*   histoEtaToPi0               [MaxNumberOfFiles];
    TH1D*   histoMassEtaData            [MaxNumberOfFiles];
    TH1D*   histoMassEtaMC              [MaxNumberOfFiles];
    TH1D*   histoWidthEtaData           [MaxNumberOfFiles];
    TH1D*   histoWidthEtaMC             [MaxNumberOfFiles];
    TH1D*   histoEtaInvMassSigPlusBG    [MaxNumberOfFiles];
    TH1D*   histoEtaInvMassSig          [MaxNumberOfFiles];
    TH1D*   histoEtaInvMassBG           [MaxNumberOfFiles];
    TF1*    fitEtaInvMassSig            [MaxNumberOfFiles];
    
    // create pointers for weighted graphs and supporting figutes
    TGraphAsymmErrors* graphCorrectedYieldWeightedAverageEtaStat    = NULL;
    TGraphAsymmErrors* graphCorrectedYieldWeightedAverageEtaSys     = NULL;
    TGraphAsymmErrors* graphMassEtaDataWeighted                     = NULL;
    TGraphAsymmErrors* graphMassEtaMCWeighted                       = NULL;
    TGraphAsymmErrors* graphWidthEtaDataWeighted                    = NULL;
    TGraphAsymmErrors* graphWidthEtaMCWeighted                      = NULL;
    TGraphAsymmErrors* graphAcceptanceEtaWeighted                   = NULL;
    TGraphAsymmErrors* graphEfficiencyEtaWeighted                   = NULL;
    TGraphAsymmErrors* graphEffTimesAccEtaWeighted                  = NULL;
    TGraphAsymmErrors* graphEtaToPi0WeightedAverageStat             = NULL;
    TGraphAsymmErrors* graphEtaToPi0WeightedAverageSys              = NULL;
    // create pointers for shrunk graphs and supporting figures
    TGraphAsymmErrors* graphsCorrectedYieldRemoved0Eta      [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldSysRemoved0Eta   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsEtaToPi0Removed0               [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsEtaToPi0SysRemoved0            [MaxNumberOfFiles];
    TGraphAsymmErrors* graphMassEtaData                     [MaxNumberOfFiles];
    TGraphAsymmErrors* graphMassEtaMC                       [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedMassEtaData              [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedMassEtaMC                [MaxNumberOfFiles];
    TGraphAsymmErrors* graphWidthEtaData                    [MaxNumberOfFiles];
    TGraphAsymmErrors* graphWidthEtaMC                      [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedWidthEtaData             [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedWidthEtaMC               [MaxNumberOfFiles];
    TGraphAsymmErrors* graphAcceptanceEta                   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedAcceptanceEta            [MaxNumberOfFiles];
    TGraphAsymmErrors* graphEfficiencyEta                   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEfficiencyEta            [MaxNumberOfFiles];
    TGraphAsymmErrors* graphEffTimesAccEta                  [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEffTimesAccEta           [MaxNumberOfFiles];
    
    TString FileNameEffBaseEta      [MaxNumberOfFiles];
    TFile*  fileEffBaseEta          [MaxNumberOfFiles];
    TH1D*   histoEffBaseEta         [MaxNumberOfFiles];
    TH1D*   histoTriggerEffEta      [MaxNumberOfFiles];
    Bool_t  enableTriggerEffEta     [MaxNumberOfFiles];
    Bool_t  enableTriggerEffEtaAll                                  = kFALSE;

    // read eta files if we are not in mode 10
    if (enableEta && mode != 10){        
        for (Int_t i=0; i< nrOfTrigToBeComb; i++){
            // Define CutSelections
            TString fEventCutSelection                      = "";
            TString fGammaCutSelection                      = "";
            TString fClusterCutSelection                    = "";
            TString fElectronCutSelection                   = "";
            TString fMesonCutSelection                      = "";

            // disentangle cut selection
            ReturnSeparatedCutNumberAdvanced(cutNumber[i].Data(),fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
        
            TString trigger                                 = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
            Int_t exampleBin                                = ReturnSingleInvariantMassBinPlotting ("Eta", optionEnergy, mode, trigger.Atoi());

            FileNameCorrectedEta[i]                         = Form("%s/%s/Eta_%s_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(), 
                                                                        cutNumber[i].Data());
            cout<< FileNameCorrectedEta[i] << endl;
            fileCorrectedEta[i]                             = new TFile(FileNameCorrectedEta[i]);
            if (fileCorrectedEta[i]->IsZombie()) return;

            // read uncorrected file
            FileNameUnCorrectedEta[i]                       = Form("%s/%s/Eta_%s_GammaConvV1WithoutCorrection_%s.root",cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(), 
                                                                   cutNumber[i].Data());
            cout<< FileNameUnCorrectedEta[i] << endl;
            fileUnCorrectedEta[i]                           = new TFile(FileNameUnCorrectedEta[i]);
            if (fileUnCorrectedEta[i]->IsZombie()) return;
            
            histoCorrectedYieldEta[i]                       = (TH1D*)fileCorrectedEta[i]->Get(nameCorrectedYield.Data());
            histoCorrectedYieldEta[i]->SetName(Form("CorrectedYield_%s",cutNumber[i].Data()));
            histoEfficiencyEta[i]                           = (TH1D*)fileCorrectedEta[i]->Get(nameEfficiency.Data());
            histoEfficiencyEta[i]->SetName(Form("Efficiency_%s",  cutNumber[i].Data()));
            histoAcceptanceEta[i]                           = (TH1D*)fileCorrectedEta[i]->Get(nameAcceptance.Data());
            histoAcceptanceEta[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
            histoAcceptanceEtaWOEvtWeights[i]               = (TH1D*)fileCorrectedEta[i]->Get(nameAcceptanceWOEvtWeights.Data());
            if(histoAcceptanceEtaWOEvtWeights[i]) histoAcceptanceEtaWOEvtWeights[i]->SetName(Form("AcceptanceWOEvtWeights_%s",  cutNumber[i].Data()));
            histoRawYieldEta[i]                             = (TH1D*)fileUnCorrectedEta[i]->Get("histoYieldMesonPerEvent");
            histoRawYieldEta[i]->SetName(Form("RAWYieldPerEvent_%s",cutNumber[i].Data()));
            histoMCInputEta[i]                              = (TH1D*)fileCorrectedEta[i]->Get("MCYield_Meson_oldBinWOWeights");
            histoMCInputEta[i]->Sumw2();
            histoMCInputEta[i]->Rebin(4);
            histoMCInputEta[i]->Scale(1./4);

            histoMCInputEta[i]->SetName(Form("Eta_Input_Reweighted_%s",cutNumber[i].Data()));
            histoMCInputEtaPi0[i]                           = (TH1D*) histoMCInputEta[i]->Clone(Form("EtaToPi0_Input_Reweighted_%s",cutNumber[i].Data()));
            histoMCInputEtaPi0[i]->Sumw2();
            histoMCInputEtaPi0[i]->Divide(histoMCInputEtaPi0[i],histoMCInputPi0[i],1.,1.,"");
            histoMassEtaData[i]                                 = (TH1D*)fileCorrectedEta[i]->Get("histoMassMeson");
            histoMassEtaData[i]->SetName(Form("Eta_Mass_data_%s",cutNumber[i].Data()));
            histoMassEtaMC[i]                                   = (TH1D*)fileCorrectedEta[i]->Get(nameMassMC.Data());
            histoMassEtaMC[i]->SetName(Form("Eta_Mass_MC_%s",cutNumber[i].Data()));
            histoWidthEtaData[i]                                = (TH1D*)fileCorrectedEta[i]->Get("histoFWHMMeson");
            histoWidthEtaData[i]->SetName(Form("Eta_Width_data_%s",cutNumber[i].Data()));
            histoWidthEtaMC[i]                                  = (TH1D*)fileCorrectedEta[i]->Get(nameWidthMC.Data());
            histoWidthEtaMC[i]->SetName(Form("Eta_Width_MC_%s",cutNumber[i].Data()));
            histoEtaInvMassSig[i]                               = (TH1D*)fileCorrectedEta[i]->Get(Form("InvMassSig_PtBin%02d",exampleBin));
            if (histoEtaInvMassSig[i]) histoEtaInvMassSig[i]->SetName(Form("Eta_InvMassSig_Example_%s",triggerName[i].Data()));
            histoEtaInvMassSigPlusBG[i]                         = (TH1D*)fileCorrectedEta[i]->Get(Form("InvMassSigPlusBG_PtBin%02d",exampleBin));
            if (histoEtaInvMassSigPlusBG[i]) histoEtaInvMassSigPlusBG[i]->SetName(Form("Eta_InvMassSigPlusBG_Example_%s",triggerName[i].Data()));
            histoEtaInvMassBG[i]                                = (TH1D*)fileCorrectedEta[i]->Get(Form("InvMassBG_PtBin%02d",exampleBin));
            if (histoEtaInvMassBG[i]) histoEtaInvMassBG[i]->SetName(Form("Eta_InvMassBG_Example_%s",triggerName[i].Data()));
            fitEtaInvMassSig[i]                                 = (TF1*)fileCorrectedEta[i]->Get(Form("FitInvMassSig_PtBin%02d",exampleBin));
            if (fitEtaInvMassSig[i]) fitEtaInvMassSig[i]->SetName(Form("Eta_InvMassSigFit_Example_%s",triggerName[i].Data()));

            if (cutNumberBaseEff[i].CompareTo("bla") != 0){
                FileNameEffBaseEta[i]                           = Form("%s/%s/Eta_MC_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), cutNumberBaseEff[i].Data());
                fileEffBaseEta[i]                               = new TFile(FileNameEffBaseEta[i]);
                if (fileEffBaseEta[i]->IsZombie()){
                    enableTriggerEffEta[i]                      = kFALSE;
                    cout << "Didn't find the effi base file for " << triggerName[i].Data() << endl;
                } else {
                    enableTriggerEffEta[i]                      = kTRUE;
                    enableTriggerEffEtaAll                      = kTRUE;
                }    
                if (enableTriggerEffEta[i]){
                    TH1D* histoEffiEtaTemp                      = (TH1D*)fileCorrectedEta[i]->Get("TrueMesonEffiPt");
                    TH1D* histoEffiBaseEtaTemp                  = (TH1D*)fileEffBaseEta[i]->Get("TrueMesonEffiPt");
                    histoEffBaseEta[i]                          = (TH1D*)fileEffBaseEta[i]->Get(nameEfficiency.Data());
                    histoEffBaseEta[i]->SetName(Form("EfficiencyBase_%s",  cutNumber[i].Data()));
                    histoTriggerEffEta[i]                       = (TH1D*)histoEffiEtaTemp->Clone(Form("TriggerEfficiency_%s", cutNumber[i].Data()));
                    histoTriggerEffEta[i]->Divide(histoTriggerEffEta[i],histoEffiBaseEtaTemp,1.,1.,"B");
                    histoEffTimesAccEta[i]                      = (TH1D*)histoEffBaseEta[i]->Clone(Form("EffTimeAcc_%s",  cutNumber[i].Data()));
                    if(histoAcceptanceEtaWOEvtWeights[i]){
                      histoAcceptanceEtaWOEvtWeights[i]->Sumw2();
                      histoEffTimesAccEta[i]->Multiply(histoAcceptanceEtaWOEvtWeights[i]);
                      histoEfficiencyEta[i]->Multiply(histoAcceptanceEta[i]);
                      histoEfficiencyEta[i]->Divide(histoEfficiencyEta[i],histoAcceptanceEtaWOEvtWeights[i],1.,1.,"B");
                      histoAcceptanceEta[i] = histoAcceptanceEtaWOEvtWeights[i];
                      histoAcceptanceEta[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
                    }else histoEffTimesAccEta[i]->Multiply(histoAcceptanceEta[i]);
                    histoEffTimesAccEta[i]->Scale(deltaRapid[i]*2*TMath::Pi());
                } else {
                    histoEffBaseEta[i]                          = NULL;
                    histoTriggerEffEta[i]                       = NULL;
                }
            } else {
                histoEffTimesAccEta[i]                          = (TH1D*)histoEfficiencyEta[i]->Clone(Form("EffTimeAcc_%s",  cutNumber[i].Data()));
                if(histoAcceptanceEtaWOEvtWeights[i]){
                  histoAcceptanceEtaWOEvtWeights[i]->Sumw2();
                  histoEffTimesAccEta[i]->Multiply(histoAcceptanceEtaWOEvtWeights[i]);
                  histoEfficiencyEta[i]->Multiply(histoAcceptanceEta[i]);
                  histoEfficiencyEta[i]->Divide(histoEfficiencyEta[i],histoAcceptanceEtaWOEvtWeights[i],1.,1.,"B");
                  histoAcceptanceEta[i] = histoAcceptanceEtaWOEvtWeights[i];
                  histoAcceptanceEta[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
                }else histoEffTimesAccEta[i]->Multiply(histoAcceptanceEta[i]);
                histoEffTimesAccEta[i]->Scale(deltaRapid[i]*2*TMath::Pi());

            }    

            // read files for eta to pi0 ratio
            FileNameCorrectedPi0EtaBin[i]                   = Form("%s/%s/Pi0EtaBinning_%s_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(), 
                                                                        cutNumber[i].Data());
            cout<< FileNameCorrectedPi0EtaBin[i] << endl;
            fileCorrectedPi0EtaBin[i]                       = new TFile(FileNameCorrectedPi0EtaBin[i]);
            if (fileCorrectedPi0EtaBin[i]->IsZombie()) {
                foundPi0EtaBinFile[i]                        = kFALSE;
            } else {
                foundPi0EtaBinFile[i]                        = kTRUE;
                doEtaToPi0                                   = kTRUE;
            }    

            if (foundPi0EtaBinFile[i]){
                // check if fit file for binshifting has been given currently only working for 2.76TeV
                if (nameFileFitsShift.CompareTo("") != 0){
                     doBinShiftForEtaToPi0                  = kTRUE;   
                     addNameBinshift                        = "YShifted";
                }
                
                if (! doBinShiftForEtaToPi0){
                    histoCorrectedYieldPi0EtaBin[i]             = (TH1D*)fileCorrectedPi0EtaBin[i]->Get(nameCorrectedYield.Data());
                    histoCorrectedYieldPi0EtaBin[i]->SetName(Form("CorrectedYieldPi0EtaBin_%s",cutNumber[i].Data()));
                    if(optionEnergy.CompareTo("8TeV")==0 && mode==4){
                      histoEtaToPi0[i]                            = (TH1D*)histoCorrectedYieldPi0EtaBin[i]->Clone(Form("EtaToPi0_%s", cutNumber[i].Data()));
                      for(Int_t iB=1; iB<=histoCorrectedYieldPi0EtaBin[i]->GetNbinsX(); iB++){histoEtaToPi0[i]->SetBinContent(iB,histoCorrectedYieldEta[i]->GetBinContent(iB));}
                      histoEtaToPi0[i]->Divide(histoEtaToPi0[i],histoCorrectedYieldPi0EtaBin[i],1.,1.,"");
                    }else{
                      histoEtaToPi0[i]                            = (TH1D*)histoCorrectedYieldEta[i]->Clone(Form("EtaToPi0_%s", cutNumber[i].Data()));
                      histoEtaToPi0[i]->Divide(histoEtaToPi0[i],histoCorrectedYieldPi0EtaBin[i],1.,1.,"");
                    }
                } else {
                    TFile *fileFitsBinShift                     = new TFile(nameFileFitsShift);
                    TF1* fitBinShiftPi0                         = (TF1*)fileFitsBinShift->Get("TsallisFitPi0");
                    TF1* fitBinShiftEta                         = (TF1*)fileFitsBinShift->Get("TsallisFitEta");
                    histoCorrectedYieldPi0EtaBin[i]             = (TH1D*)fileCorrectedPi0EtaBin[i]->Get(nameCorrectedYield.Data());
                    histoCorrectedYieldPi0EtaBin[i]->SetName(Form("CorrectedYieldPi0EtaBin_%s",cutNumber[i].Data()));
                    cout << "shifting pi0 in eta binning: " <<  cutNumber[i].Data() << endl;
                    TH1D* histoCorrectedYieldPi0EtaBinBinShift  = (TH1D*)histoCorrectedYieldPi0EtaBin[i]->Clone(Form("CorrectedYieldPi0EtaBinBinShifted_%s",cutNumber[i].Data()));
                    histoCorrectedYieldPi0EtaBinBinShift        = ApplyYshiftIndividualSpectra( histoCorrectedYieldPi0EtaBinBinShift, fitBinShiftPi0);
                    cout << "shifting eta: " <<  cutNumber[i].Data() << endl;
                    TH1D* histoCorrectedYieldEtaBinShift        = (TH1D*)histoCorrectedYieldEta[i]->Clone(Form("CorrectedYieldEtaBinShifted_%s",cutNumber[i].Data()));
                    histoCorrectedYieldEtaBinShift              = ApplyYshiftIndividualSpectra( histoCorrectedYieldEtaBinShift, fitBinShiftEta);

                    if(optionEnergy.CompareTo("8TeV")==0 && mode==4){
                      histoEtaToPi0[i]                            = (TH1D*)histoCorrectedYieldPi0EtaBinBinShift->Clone(Form("EtaToPi0_%s", cutNumber[i].Data()));
                      for(Int_t iB=1; iB<=histoCorrectedYieldPi0EtaBinBinShift->GetNbinsX(); iB++){histoEtaToPi0[i]->SetBinContent(iB,histoCorrectedYieldEtaBinShift->GetBinContent(iB));}
                      histoEtaToPi0[i]->Divide(histoEtaToPi0[i],histoCorrectedYieldPi0EtaBinBinShift,1.,1.,"");
                    }else{
                      histoEtaToPi0[i]                            = (TH1D*)histoCorrectedYieldEtaBinShift->Clone(Form("EtaToPi0%s_%s", addNameBinshift.Data(), cutNumber[i].Data()));
                      histoEtaToPi0[i]->Divide(histoEtaToPi0[i],histoCorrectedYieldPi0EtaBinBinShift,1.,1.,"");
                    }
                }    
            }
            
            //Scale spectrum to MBOR
            if (optionEnergy.CompareTo("2.76TeV")==0 && 
                (triggerName[i].Contains("INT7")|| triggerName[i].Contains("EMC7") || triggerName[i].Contains("EG1") || triggerName[i].Contains("EG2")) ){
                // &&
                // isMC.CompareTo("data") == 0){
                histoCorrectedYieldEta[i]->Scale(0.8613) ;
                histoRawYieldEta[i]->Scale(0.8613) ;
            }
        }

        //***************************************************************************************************************
        //************************************Plotting efficiencies Eta *************************************************
        //***************************************************************************************************************
        Double_t minEffiEta     = 1e-5;
        Double_t maxEffiEta     = 1e-1;
        if (mode == 4){
            maxEffiEta          = 5e0;
        } else if (mode == 0) {
            maxEffiEta          = 8e-3;
            minEffiEta          = 1e-4;
        }    
        
        if(optionEnergy.CompareTo("8TeV")==0){
          if(mode == 2){
            minEffiEta        = 1e-3;
            maxEffiEta        = 2e-1;
          }else if(mode == 4){
            minEffiEta        = 8e-3;
            maxEffiEta        = 1;
          }
        }

        canvasEffi->cd();
        TH2F * histo2DEffiEta;
        histo2DEffiEta = new TH2F("histo2DEffiEta","histo2DEffiEta",1000,0., maxPtGlobalEta,10000,minEffiEta, maxEffiEta);
        SetStyleHistoTH2ForGraphs(histo2DEffiEta, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco} #times #epsilon_{trigg}",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
        histo2DEffiEta->DrawCopy(); 
        histo2DEffiEta->SetYTitle("#epsilon_{#eta}");

        TLegend* legendEffiEta = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
        legendEffiEta->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            DrawGammaSetMarker(histoEfficiencyEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoEfficiencyEta[i]->DrawCopy("e1,same"); 
            legendEffiEta->AddEntry(histoEfficiencyEta[i],triggerNameLabel[i].Data(),"p");
        }
        legendEffiEta->Draw();

        labelEnergyEffi->Draw();
        TLatex *labelEtaEffi = new TLatex(0.62, maxYLegendEffi+0.02+0.99*textSizeSpectra*0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelEtaEffi, 0.85*textSizeSpectra,4);
        labelEtaEffi->Draw();
        labelDetProcEffi->Draw();
        
        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/Eta_Efficiency.%s",outputDir.Data(),suffix.Data()));

        //***************************************************************************************************************
        //************************************ Plotting trigger efficiencies Eta ****************************************
        //***************************************************************************************************************
        if (enableTriggerEffEtaAll){
            TCanvas* canvasTriggerEffi = new TCanvas("canvasTriggerEffi","",0,0,1000,900);// gives the page size
            DrawGammaCanvasSettings( canvasTriggerEffi, 0.09, 0.017, 0.015, 0.08);
            canvasTriggerEffi->SetLogy(0);

            Double_t minEffiTrigEta         = 0;
            Double_t maxEffiTrigEta         = 1.1;
            if(optionEnergy.CompareTo("8TeV") == 0 && mode == 2) maxEffiTrigEta = 1.5;
            
            TH2F * histo2DTriggerEffiEta;
            histo2DTriggerEffiEta = new TH2F("histo2DTriggerEffiEta","histo2DTriggerEffiEta",1000,0., maxPtGlobalEta,10000,minEffiTrigEta, maxEffiTrigEta);
            SetStyleHistoTH2ForGraphs(histo2DTriggerEffiEta, "#it{p}_{T} (GeV/#it{c})","#epsilon_{Trigger, #eta}", 
                                        0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
            histo2DTriggerEffiEta->DrawCopy(); 

            TLegend* legendTriggerEffiEta = GetAndSetLegend2(0.62, 0.165, 0.95, 0.165+(1.05*nRealTriggers/2*0.85*textSizeSpectra),28);
            legendTriggerEffiEta->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if (enableTriggerEffEta[i]){
                    DrawGammaSetMarker(histoTriggerEffEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoTriggerEffEta[i]->DrawCopy("e1,same"); 
                    legendTriggerEffiEta->AddEntry(histoTriggerEffEta[i],triggerNameLabel[i].Data(),"p"); 
                }    
            }
            legendTriggerEffiEta->Draw();

            labelEnergyEffi->Draw();
            labelEtaEffi->Draw();
            labelDetProcEffi->Draw();

            canvasTriggerEffi->Update();
            canvasTriggerEffi->SaveAs(Form("%s/Eta_TriggerEfficiency.%s",outputDir.Data(),suffix.Data()));
            delete canvasTriggerEffi;

            //***************************************************************************************************************
            //************************************ Plotting standalone efficiencies for eta *********************************
            //***************************************************************************************************************
        
            canvasEffi->cd();
            histo2DEffiEta->DrawCopy(); 

            TLegend* legendEffiEtaW0TriggEff = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
            legendEffiEtaW0TriggEff->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if ( triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1") ){
                    DrawGammaSetMarker(histoEfficiencyEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoEfficiencyEta[i]->DrawCopy("e1,same"); 
                    legendEffiEtaW0TriggEff->AddEntry(histoEfficiencyEta[i],triggerNameLabel[i].Data(),"p"); 
                } else {
                    if (enableTriggerEffEta[i]){
                        DrawGammaSetMarker(histoEffBaseEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                        histoEffBaseEta[i]->DrawCopy("e1,same"); 
                        legendEffiEtaW0TriggEff->AddEntry(histoEffBaseEta[i],triggerNameLabel[i].Data(),"p"); 
                    }
                }    
            }
            legendEffiEtaW0TriggEff->Draw();

            labelEnergyEffi->Draw();
            labelEtaEffi->Draw();
            labelDetProcEffi->Draw();

            canvasEffi->Update();
            canvasEffi->SaveAs(Form("%s/Eta_EfficiencyW0TriggEff.%s",outputDir.Data(),suffix.Data()));
        }

        
        //***************************************************************************************************************
        //************************************Plotting acceptance Eta *************************************************
        //***************************************************************************************************************            
        canvasAcc->cd();
        Double_t minAccEta = 0.05;
        Double_t maxAccEta = 0.3;
        if (mode == 0){
            maxAccEta       = 1.05;
            minAccEta       = 0.5;
        }


        if(optionEnergy.CompareTo("8TeV")==0){
            if(mode == 2){
              minAccEta   = 0.2;
              maxAccEta   = 0.31;
            }else if(mode == 4){
              minAccEta = 0.05;
              maxAccEta = 0.3;
            }
        }

        
        TH2F * histo2DAccEta;
        histo2DAccEta = new TH2F("histo2DAccEta","histo2DAccEta",1000,0., maxPtGlobalEta,10000,minAccEta, maxAccEta);
        SetStyleHistoTH2ForGraphs(histo2DAccEta, "#it{p}_{T} (GeV/#it{c})","A_{#eta}", 
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
        histo2DAccEta->DrawCopy(); 

        TLegend* legendAccEta = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
        legendAccEta->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            DrawGammaSetMarker(histoAcceptanceEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoAcceptanceEta[i]->DrawCopy("e1,same"); 
            legendAccEta->AddEntry(histoAcceptanceEta[i],triggerNameLabel[i].Data(),"p");
        }
        legendAccEta->Draw();

        labelEnergyEffi->Draw();
        labelEtaEffi->Draw();
        labelDetProcEffi->Draw();
        
        canvasAcc->Update();
        canvasAcc->SaveAs(Form("%s/Eta_Acceptance.%s",outputDir.Data(),suffix.Data()));

        //***************************************************************************************************************
        //************************************Plotting Mass Eta *********************************************************
        //***************************************************************************************************************
        canvasMass->cd();
        
        Double_t minMassEta = 0.500;
        Double_t maxMassEta = 0.600;
    
        if (mode == 0){
            minMassEta      = 0.54;
            maxMassEta      = 0.56;
        }
        
        TH2F * histo2DMassEta;
        histo2DMassEta = new TH2F("histo2DMassEta","histo2DMassEta",1000,0., maxPtGlobalEta,10000,minMassEta, maxMassEta);
        SetStyleHistoTH2ForGraphs(histo2DMassEta, "#it{p}_{T} (GeV/#it{c})","#it{M}_{#eta} (Mev/#it{c}^{2})", 
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.4);
        histo2DMassEta->DrawCopy(); 

        TLegend* legendMassEta = GetAndSetLegend2(0.52, 0.88, 0.95, 0.88+(1.05*4/2*0.85*textSizeSpectra),28);
        legendMassEta->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0)
             ){
                DrawGammaSetMarker(histoMassEtaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoMassEtaData[i]->DrawCopy("e1,same"); 
                legendMassEta->AddEntry(histoMassEtaData[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                DrawGammaSetMarker(histoMassEtaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                histoMassEtaMC[i]->DrawCopy("e1,same"); 
                legendMassEta->AddEntry(histoMassEtaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
            }    
        }
        legendMassEta->Draw();
        labelEnergyMass->Draw();

        TLatex *labelEtaMass = new TLatex(0.14, 0.85+0.99*textSizeSpectra*0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelEtaMass, 0.85*textSizeSpectra,4);
        labelEtaMass->Draw();
        labelDetProcMass->Draw();

        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/Eta_%s_Mass1.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        histo2DMassEta->DrawCopy(); 

        TLegend* legendMassEta2 = GetAndSetLegend2(0.46, 0.81, 0.80, 0.81+(1.05*8/2*0.85*textSizeSpectra),28);
        legendMassEta2->SetNColumns(2);
        
        if (optionEnergy.CompareTo("2.76TeV")==0 ){
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if( i==1 || i==3 || i==4 || i==5 ){
                    DrawGammaSetMarker(histoMassEtaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoMassEtaData[i]->DrawCopy("e1,same"); 
                    legendMassEta2->AddEntry(histoMassEtaData[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                    DrawGammaSetMarker(histoMassEtaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    histoMassEtaMC[i]->DrawCopy("e1,same"); 
                    legendMassEta2->AddEntry(histoMassEtaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                }    
            }
            legendMassEta2->Draw();
            labelEnergyMass->Draw();
            labelEtaMass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Eta_%s_Mass2.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }

        //***************************************************************************************************************
        //************************************Plotting Width Eta *********************************************************
        //***************************************************************************************************************
        canvasWidth->cd();
        
        Double_t minWidthEta    = 0.0;
        Double_t maxWidthEta    = 0.060;
        if (mode == 0){
            minWidthEta         = 0.0;
            maxWidthEta         = 0.012;
        }
        
        TH2F * histo2DWidthEta;
        histo2DWidthEta = new TH2F("histo2DWidthEta","histo2DWidthEta",1000,0., maxPtGlobalEta,10000,minWidthEta, maxWidthEta);
        SetStyleHistoTH2ForGraphs(histo2DWidthEta, "#it{p}_{T} (GeV/#it{c})","#sigma_{#eta} (Mev/#it{c}^{2})", 
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
        histo2DWidthEta->DrawCopy(); 

        TLegend* legendWidthEta = GetAndSetLegend2(0.52, 0.87, 0.95, 0.87+(1.05*4/2*0.85*textSizeSpectra),28);
        legendWidthEta->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0)
             ){
                DrawGammaSetMarker(histoWidthEtaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoWidthEtaData[i]->DrawCopy("e1,same"); 
                legendWidthEta->AddEntry(histoWidthEtaData[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                DrawGammaSetMarker(histoWidthEtaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                histoWidthEtaMC[i]->DrawCopy("e1,same"); 
                legendWidthEta->AddEntry(histoWidthEtaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
            }    
        }
        legendWidthEta->Draw();
        labelEnergyWidth->Draw();

        TLatex *labelEtaWidth = new TLatex(0.14, 0.84+0.99*textSizeSpectra*0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelEtaWidth, 0.85*textSizeSpectra,4);
        labelEtaWidth->Draw();
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Eta_%s_Width1.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        histo2DWidthEta->DrawCopy(); 

        TLegend* legendWidthEta2 = GetAndSetLegend2(0.46, 0.80, 0.80, 0.80+(1.05*8/2*0.85*textSizeSpectra),28);
        legendWidthEta2->SetNColumns(2);
        
        if (optionEnergy.CompareTo("2.76TeV")==0 ){
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if( i==1 || i==3 || i==4 || i==5 ){
                    DrawGammaSetMarker(histoWidthEtaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoWidthEtaData[i]->DrawCopy("e1,same"); 
                    legendWidthEta2->AddEntry(histoWidthEtaData[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                    DrawGammaSetMarker(histoWidthEtaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    histoWidthEtaMC[i]->DrawCopy("e1,same"); 
                    legendWidthEta2->AddEntry(histoWidthEtaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                }    
            }
            legendWidthEta2->Draw();
            labelEnergyWidth->Draw();
            labelEtaWidth->Draw();
            labelDetProcWidth->Draw();

            canvasWidth->Update();
            canvasWidth->SaveAs(Form("%s/Eta_%s_Width2.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        //***************************************************************************************************************
        //************************************Plotting unscaled invariant raw-yield Eta *********************************
        //***************************************************************************************************************
        canvasRawUnscaled->cd();
        
        Double_t minCorrYieldRawUnscaledEta     = 7e-8;
        Double_t maxCorrYieldRawUnscaledEta     = 4e-3;
        if (mode == 10) {
            minCorrYieldRawUnscaledEta          = 2e-12;
            maxCorrYieldRawUnscaledEta          = 1;
        } else if (mode == 4) {
            minCorrYieldRawUnscaledEta          = 7e-8;
            maxCorrYieldRawUnscaledEta          = 1e-1;
        } else if (mode == 0) {
            minCorrYieldRawUnscaledEta          = 7e-8;
            maxCorrYieldRawUnscaledEta          = 1e-4;            
        }

        if(optionEnergy.CompareTo("8TeV")==0){
          if(mode == 2){
            minCorrYieldRawUnscaledEta        = 7e-8;
            maxCorrYieldRawUnscaledEta        = 8e-4;
          }else if(mode == 4){
            minCorrYieldRawUnscaledEta        = 2e-8;
            maxCorrYieldRawUnscaledEta        = 8e-3;
          }
        }
        
        TH2F * histo2DRawUnscaledEta       = new TH2F("histo2DRawUnscaledEta", "histo2DRawUnscaledEta", 1000, 0., maxPtGlobalEta, 10000, minCorrYieldRawUnscaledEta, maxCorrYieldRawUnscaledEta);
        SetStyleHistoTH2ForGraphs(histo2DRawUnscaledEta, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{#eta, raw}}{#it{N}_{evt}d#it{p}_{T}} (#it{c}/GeV)^{2}", 
                                0.85*textSizeSpectra,0.04, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.7);
        histo2DRawUnscaledEta->GetXaxis()->SetLabelOffset(-0.005);
        histo2DRawUnscaledEta->DrawCopy(); 

        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            DrawGammaSetMarker(histoRawYieldEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoRawYieldEta[i]->DrawCopy("e1,same"); 
        }
        legendRawUnscaled->Draw();
        labelEnergyRawUnscaled->Draw();
        TLatex *labelEtaRawUnscaled     = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelEtaRawUnscaled, 0.85*textSizeSpectra,4);
        labelEtaRawUnscaled->Draw();
        labelDetProcRawUnscaled->Draw();

        canvasRawUnscaled->Update();
        canvasRawUnscaled->SaveAs(Form("%s/Eta_%s_RawYieldUnscaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        delete canvasRawUnscaled;
        
        //***************************************************************************************************************
        //************************************Plotting unscaled invariant yield Eta *************************************
        //***************************************************************************************************************
        canvasCorrUnscaled->cd();
        
        Double_t minCorrYieldUnscaledEta    = 2e-10;
        Double_t maxCorrYieldUnscaledEta    = 1e2;
        if (mode == 0){
            minCorrYieldUnscaledEta         = 1e-7;
            maxCorrYieldUnscaledEta         = 5e-1;        
        }    

        if(optionEnergy.CompareTo("8TeV")==0){
          if(mode == 2){
            minCorrYieldUnscaledEta         = 1e-7;
            maxCorrYieldUnscaledEta         = 4e-2;
          }else if(mode == 4){
            minCorrYieldUnscaledEta         = 9e-10;
            maxCorrYieldUnscaledEta         = 5e-2;
          }
        }
        
        TH2F * histo2DInvYieldUnscaledEta;
        histo2DInvYieldUnscaledEta = new TH2F("histo2DInvYieldUnscaledEta","histo2DInvYieldUnscaledEta",1000,0., maxPtGlobalEta,10000,minCorrYieldUnscaledEta,maxCorrYieldUnscaledEta);
        SetStyleHistoTH2ForGraphs(histo2DInvYieldUnscaledEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.55);
        histo2DInvYieldUnscaledEta->DrawCopy(); 

        TLegend* legendUnscaledEta = GetAndSetLegend2(0.72, 0.95-(1.15*nrOfTrigToBeComb*0.85*textSizeSpectra), 0.95, 0.95,32);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            DrawGammaSetMarker(histoCorrectedYieldEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoCorrectedYieldEta[i]->DrawCopy("e1,same"); 
            legendUnscaledEta->AddEntry(histoCorrectedYieldEta[i],triggerNameLabel[i].Data(),"p");
        }
        legendUnscaledEta->Draw();

        labelEnergyUnscaled->Draw();
        TLatex *labelEtaUnscaled = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelEtaUnscaled, 0.85*textSizeSpectra,4);
        labelEtaUnscaled->Draw();
        labelDetProcUnscaled->Draw();

        
        canvasCorrUnscaled->Update();
        canvasCorrUnscaled->SaveAs(Form("%s/Eta_%s_CorrectedYieldUnscaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        delete canvasCorrUnscaled;

        //***************************************************************************************************************
        //******************************* Scaling corrected yield by trigger rejection factors Eta **********************
        //***************************************************************************************************************
        
        cout << "**************************************************************************************************" << endl;
        cout << "******************************* Combining different triggers for eta *****************************" << endl;
        cout << "**************************************************************************************************" << endl;
        
        TH1D*     histoCorrectedYieldEtaScaled                  [MaxNumberOfFiles];
        TH1D*     histoCorrectedYieldEtaScaledMasked            [MaxNumberOfFiles];
        histoCorrectedYieldEtaScaled[0]                         = (TH1D*)histoCorrectedYieldEta[0]->Clone(Form("CorrectedYieldEtaScaled_%s", triggerName[0].Data()));
        for (Int_t i = 1; i< nrOfTrigToBeComb; i++){
            cout << triggerName[i].Data() << endl;
            histoCorrectedYieldEtaScaled[i] = (TH1D*)histoCorrectedYieldEta[i]->Clone(Form("CorrectedYieldEtaScaled_%s", triggerName[i].Data()));
            histoCorrectedYieldEtaScaled[i]->Sumw2();
            histoCorrectedYieldEtaScaled[i]->Scale(1./triggRejecFac[i][trigSteps[i][0]]);
            if (trigSteps[i][1]!= trigSteps[i][0]){
                cout << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << endl;
                histoCorrectedYieldEtaScaled[i]->Scale(1./triggRejecFac[trigSteps[i][0]][trigSteps[i][1]]);
            }
            if (trigSteps[i][2]!= trigSteps[i][1]){
                cout << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << "\t"<< triggRejecFac[trigSteps[i][1]][trigSteps[i][2]] << endl;
                histoCorrectedYieldEtaScaled[i]->Scale(1./triggRejecFac[trigSteps[i][1]][trigSteps[i][2]]);
            }
        }
        
        // prepare arrays for systematics 
        Double_t xValueFinalEta                                 [100];
        Double_t xErrorLowFinalEta                              [100];
        Double_t xErrorHighFinalEta                             [100];
        Double_t yValueFinalEta                                 [100];
        Double_t yErrorLowFinalEta                              [100];
        Double_t yErrorHighFinalEta                             [100];
        Int_t nPointFinalEta                                         = 0;

        Double_t yErrorSysLowFinalEta                           [100];
        Double_t yErrorSysHighFinalEta                          [100];

        Double_t ptSysRelEta                                    [MaxNumberOfFiles][100];
        Double_t yErrorSysLowRelEta                             [MaxNumberOfFiles][100];
        Double_t yErrorSysHighRelEta                            [MaxNumberOfFiles][100];
        Bool_t     sysAvailEta                                  [MaxNumberOfFiles];

        Bool_t sysAvailSingleEta                                [MaxNumberOfFiles];
        Int_t numberBinsSysAvailSingleEta                       [MaxNumberOfFiles];
        
        // create graphs for shrunk spectrum
        TGraphAsymmErrors* graphsCorrectedYieldShrunkEta        [MaxNumberOfFiles];
        TGraphAsymmErrors* graphsCorrectedYieldSysShrunkEta     [MaxNumberOfFiles];

        // create inputs for combination function
        TH1D*               histoStatEta    [12]; 
        TGraphAsymmErrors*  graphSystEta    [12];
        TH1D*               histoRelStatEta [12]; 
        TGraphAsymmErrors*  graphRelSystEta [12];

        Int_t offSetsEta[12]            = { 0, 0, 0, 0, 0, 0, 
                                            0, 0, 0, 0, 0, 0 };
        Int_t offSetsEtaSys[12]         = { 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0 };

        if(optionEnergy.CompareTo("8TeV")==0){
          if(mode == 2){
            offSetsEta[1] = 0; //INT7
            offSetsEta[3] = 0; //EMC7
            offSetsEta[4] = 2; //EGA
          }else if(mode == 4){
            offSetsEta[1] = 0; //INT7
            offSetsEta[3] = 0; //EMC7
            offSetsEta[4] = 2; //EGA
          }
        }

        Bool_t hasSysEta                = kFALSE;
        for (Int_t j = 0; j<12; j++){
            histoStatEta[j]                 = NULL;
            graphSystEta[j]                 = NULL;
            histoRelStatEta[j]              = NULL;
            graphRelSystEta[j]              = NULL;
            graphOrderedMassEtaData[j]      = NULL;
            graphOrderedMassEtaMC[j]        = NULL;
            graphOrderedWidthEtaData[j]     = NULL;
            graphOrderedWidthEtaMC[j]       = NULL;
            graphOrderedAcceptanceEta[j]    = NULL;
            graphOrderedEfficiencyEta[j]    = NULL;
            graphOrderedEffTimesAccEta[j]   = NULL;
        }
        
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            // Read systematics file and fill arrays
            cout << triggerName[i].Data() << endl;
            if (sysFileEta[i].CompareTo("bla") != 0){
                sysAvailEta[i]              = kTRUE;
                ifstream  fileSysErrEta;
                fileSysErrEta.open(sysFileEta[i].Data(),ios_base::in);
                cout << sysFileEta[i].Data() << endl;
                Int_t counter = 0;
                while(!fileSysErrEta.eof() && counter < 100){
                    Double_t garbage = 0;
                    fileSysErrEta >>ptSysRelEta[i][counter] >> yErrorSysLowRelEta[i][counter] >> yErrorSysHighRelEta[i][counter]>>    garbage >> garbage;
                    cout << counter << "\t"<< ptSysRelEta[i][counter]<< "\t"  << yErrorSysLowRelEta[i][counter] << "\t"  <<yErrorSysHighRelEta[i][counter] << "\t"  << endl;;
                    counter++;
                }
                fileSysErrEta.close();
                hasSysEta                   = kTRUE;
             // read in detailed systematics
                string sysFileEtaDet = sysFileEta[i].Data();
                if(!replace(sysFileEtaDet, "Averaged", "AveragedSingle")){cout << "WARNING: could not find detailed systematics file " << sysFileEtaDet << ", skipping... " << endl; sysAvailSingleEta[i] = kFALSE; continue;}
                ifstream fileSysErrDetailedEta;
                fileSysErrDetailedEta.open(sysFileEtaDet,ios_base::in);
                if(fileSysErrDetailedEta.is_open()) 
                    sysAvailSingleEta[i] = kTRUE;
                else{ 
                    sysAvailSingleEta[i] = kFALSE; 
                    cout << "couldn't find single errors for eta, jumping" << endl;
                }
                
                if (sysAvailSingleEta[i]){
                    cout << sysFileEtaDet << endl;
                    counter = 0;
                    string line;
                    Int_t counterColumn = 0;
                    while (getline(fileSysErrDetailedEta, line) && counter < 100) {
                        istringstream ss(line);
                        TString temp="";
                        counterColumn = 0;
                        while(ss && counterColumn < 100){
                            ss >> temp;
                            if( !(counter==0 && temp.CompareTo("bin")==0) && !temp.IsNull()){
                            ptSysDetail[i][counter].push_back(temp);
                            counterColumn++;
                            }
                        }
                        if(counter == 0){
                            ptSysDetail[i][counter++].push_back("TotalError");
                            counterColumn++;
                        }else counter++;
                    }
                    numberBinsSysAvailSingleEta[i] = counter;
                    fileSysErrDetailedEta.close();
                }    
            } else {
                sysAvailEta[i]              = kFALSE;
                sysAvailSingleEta[i]        = kTRUE;
            }
            
            cout << "step 1" << endl;
            // intialize graphs for individual triggers to be shrunk
            graphsCorrectedYieldShrunkEta[i]        = new TGraphAsymmErrors(histoCorrectedYieldEtaScaled[i]);
            graphsCorrectedYieldRemoved0Eta[i]      = new TGraphAsymmErrors(histoCorrectedYieldEtaScaled[i]);
            graphsCorrectedYieldSysShrunkEta[i]     = new TGraphAsymmErrors(histoCorrectedYieldEtaScaled[i]);
            graphsCorrectedYieldSysRemoved0Eta[i]   = new TGraphAsymmErrors(histoCorrectedYieldEtaScaled[i]);
            graphMassEtaData[i]                     = new TGraphAsymmErrors(histoMassEtaData[i]);
            graphMassEtaMC[i]                       = new TGraphAsymmErrors(histoMassEtaMC[i]);
            graphWidthEtaData[i]                    = new TGraphAsymmErrors(histoWidthEtaData[i]);
            graphWidthEtaMC[i]                      = new TGraphAsymmErrors(histoWidthEtaMC[i]);
            graphAcceptanceEta[i]                   = new TGraphAsymmErrors(histoAcceptanceEta[i]);
            graphEffTimesAccEta[i]                  = new TGraphAsymmErrors(histoEffTimesAccEta[i]);
            if ( triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1") ){
                graphEfficiencyEta[i]               = new TGraphAsymmErrors(histoEfficiencyEta[i]);
            } else {
                if (enableTriggerEffEta[i]){
                    graphEfficiencyEta[i]           = new TGraphAsymmErrors(histoEffBaseEta[i]);
                } else {
                    graphEfficiencyEta[i]           = NULL;
                }    
            }

            histoCorrectedYieldEtaScaledMasked[i]   = (TH1D*)histoCorrectedYieldEtaScaled[i]->Clone(Form("Eta_ScaledMasked_%s",triggerName[i].Data()));
            
            // remove points at beginning according to ranges set
            Int_t binsToMask = 1;
            while (histoCorrectedYieldEtaScaledMasked[i]->GetBinCenter(binsToMask) < ptFromSpecEta[i][0] ){
                histoCorrectedYieldEtaScaledMasked[i]->SetBinContent(binsToMask,0.);
                histoCorrectedYieldEtaScaledMasked[i]->SetBinError(binsToMask,0.);
                binsToMask++;
            }    
            while (graphMassEtaData[i]->GetX()[0] < ptFromSpecEta[i][0] ){
                graphMassEtaData[i]->RemovePoint(0);
                graphMassEtaMC[i]->RemovePoint(0);
                graphWidthEtaData[i]->RemovePoint(0);
                graphWidthEtaMC[i]->RemovePoint(0);
                graphAcceptanceEta[i]->RemovePoint(0);
                graphEffTimesAccEta[i]->RemovePoint(0);
                if (enableTriggerEffEta[i] || (triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1"))){
                    graphEfficiencyEta[i]->RemovePoint(0);
                }
            }

            // mask unused triggers completely
            if ( maskedFullyEta[i] ){
              graphMassEtaData[i]     = NULL;
              graphMassEtaMC[i]       = NULL;
              graphWidthEtaData[i]    = NULL;
              graphWidthEtaMC[i]      = NULL;
              graphAcceptanceEta[i]   = NULL;
              graphEffTimesAccEta[i]  = NULL;
              graphEfficiencyEta[i]   = NULL;
              graphsCorrectedYieldShrunkEta[i]    = NULL;
              graphsCorrectedYieldRemoved0Eta[i]  = NULL; 
              graphsCorrectedYieldSysShrunkEta[i] = NULL; 
              graphsCorrectedYieldSysRemoved0Eta[i]   = NULL; 
              sysAvailEta[i]          = kFALSE;
              for (Int_t f = 1; f < histoCorrectedYieldPi0ScaledMasked[i]->GetNbinsX()+1; f++ ){
                  histoCorrectedYieldPi0ScaledMasked[i]->SetBinContent(f,0.);
                  histoCorrectedYieldPi0ScaledMasked[i]->SetBinError(f,0.);
              }
              nrOfTrigToBeCombEtaRed--;
              cout << "trigger " << triggerName[i] << " was masked for Eta" << endl;
              continue;
            // shorten graphs at the end according to range set in ptFromSpecEta
            } else if (ptFromSpecEta[i][1] > -1) {
                for (Int_t f = histoCorrectedYieldEtaScaledMasked[i]->GetXaxis()->FindBin(ptFromSpecEta[i][1]); f < histoCorrectedYieldEtaScaledMasked[i]->GetNbinsX()+1; f++ ){
                    histoCorrectedYieldEtaScaledMasked[i]->SetBinContent(f,0.);
                    histoCorrectedYieldEtaScaledMasked[i]->SetBinError(f,0.);
                }    
                while (graphMassEtaData[i]->GetX()[graphMassEtaData[i]->GetN()-1] > ptFromSpecEta[i][1] ){
                    graphMassEtaData[i]->RemovePoint(graphMassEtaData[i]->GetN()-1);
                    graphMassEtaMC[i]->RemovePoint(graphMassEtaMC[i]->GetN()-1);
                    graphWidthEtaData[i]->RemovePoint(graphWidthEtaData[i]->GetN()-1);
                    graphWidthEtaMC[i]->RemovePoint(graphWidthEtaMC[i]->GetN()-1);      
                    graphAcceptanceEta[i]->RemovePoint(graphAcceptanceEta[i]->GetN()-1);      
                    graphEffTimesAccEta[i]->RemovePoint(graphEffTimesAccEta[i]->GetN()-1);
                    if (enableTriggerEffEta[i] || (triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1"))){
                        graphEfficiencyEta[i]->RemovePoint(graphEfficiencyEta[i]->GetN()-1);
                    }
                }                
            }    
            // check if any points are still in graphs, if not put to NULL
            if (graphMassEtaData[i]->GetN() == 0){
                graphMassEtaData[i]     = NULL;
                graphMassEtaMC[i]       = NULL;
                graphWidthEtaData[i]    = NULL;
                graphWidthEtaMC[i]      = NULL;
                graphAcceptanceEta[i]   = NULL;
                graphEfficiencyEta[i]   = NULL;
                graphEffTimesAccEta[i]  = NULL;
            }
        
            // Remove 0 at beginning of the graphs
//             graphsCorrectedYieldShrunkEta[i]->Print();
            cout << "step 2" << endl;
            while (graphsCorrectedYieldShrunkEta[i]->GetY()[0] == 0) graphsCorrectedYieldShrunkEta[i]->RemovePoint(0);
//             graphsCorrectedYieldShrunkEta[i]->Print();
            while (graphsCorrectedYieldRemoved0Eta[i]->GetY()[0] == 0) graphsCorrectedYieldRemoved0Eta[i]->RemovePoint(0);
//             graphsCorrectedYieldRemoved0Eta[i]->Print();
            cout << "sys shrunk" << endl;
            while (graphsCorrectedYieldSysShrunkEta[i]->GetY()[0] == 0) graphsCorrectedYieldSysShrunkEta[i]->RemovePoint(0);
//             graphsCorrectedYieldSysShrunkEta[i]->Print();
            cout << "sys shrunk 2" << endl;
            while (graphsCorrectedYieldSysRemoved0Eta[i]->GetY()[0] == 0) graphsCorrectedYieldSysRemoved0Eta[i]->RemovePoint(0);
//             graphsCorrectedYieldSysRemoved0Eta[i]->Print();
            
            // shrink systematics graphs at the end & fill systematics
            for (Int_t j = 0; j< graphsCorrectedYieldSysRemoved0Eta[i]->GetN(); j++){
                if (sysAvailEta[i]){
                    Int_t counter = 0;
                    while(counter < 100 && TMath::Abs(graphsCorrectedYieldSysRemoved0Eta[i]->GetX()[j] - ptSysRelEta[i][counter])> 0.001) counter++;
                    if (counter < 100){
                        cout << ptSysRelEta[i][counter]<< "\t found it" << endl;
                        Double_t yErrorSysLowDummy = TMath::Abs(yErrorSysLowRelEta[i][counter]/100*graphsCorrectedYieldSysRemoved0Eta[i]->GetY()[j]);
                        Double_t yErrorSysHighDummy = yErrorSysHighRelEta[i][counter]/100*graphsCorrectedYieldSysRemoved0Eta[i]->GetY()[j];
                        graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYlow(j,yErrorSysLowDummy);
                        graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYhigh(j,yErrorSysHighDummy);
                    } else {
                        graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYlow(j,0);
                        graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYhigh(j,0);
                    }
                } else {
                    graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYlow(j,0);
                    graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYhigh(j,0);
                    averagedEta = kFALSE;
                }
            }
            
            // Shrink graphs to desired range from below
            cout << "step 3" << endl;
            while (graphsCorrectedYieldShrunkEta[i]->GetX()[0] < ptFromSpecEta[i][0]) 
                graphsCorrectedYieldShrunkEta[i]->RemovePoint(0);
            while (graphsCorrectedYieldSysShrunkEta[i]->GetX()[0] < ptFromSpecEta[i][0]) 
                graphsCorrectedYieldSysShrunkEta[i]->RemovePoint(0);
            while (graphsCorrectedYieldShrunkEta[i]->GetX()[graphsCorrectedYieldShrunkEta[i]->GetN()-1] > ptFromSpecEta[i][1]) 
                graphsCorrectedYieldShrunkEta[i]->RemovePoint(graphsCorrectedYieldShrunkEta[i]->GetN()-1);
            graphsCorrectedYieldShrunkEta[i]->Print();
            while (graphsCorrectedYieldSysShrunkEta[i]->GetX()[graphsCorrectedYieldSysShrunkEta[i]->GetN()-1] > ptFromSpecEta[i][1]) 
                graphsCorrectedYieldSysShrunkEta[i]->RemovePoint(graphsCorrectedYieldSysShrunkEta[i]->GetN()-1);

            // Shrink graphs to desired range from above and fill systematics
            for (Int_t j = 0; j< graphsCorrectedYieldShrunkEta[i]->GetN(); j++){
                xValueFinalEta[nPointFinalEta] = graphsCorrectedYieldShrunkEta[i]->GetX()[j];
                xErrorHighFinalEta[nPointFinalEta] = graphsCorrectedYieldShrunkEta[i]->GetEXhigh()[j];
                xErrorLowFinalEta[nPointFinalEta] = graphsCorrectedYieldShrunkEta[i]->GetEXlow()[j];
                yValueFinalEta[nPointFinalEta] = graphsCorrectedYieldShrunkEta[i]->GetY()[j];
                yErrorHighFinalEta[nPointFinalEta] = graphsCorrectedYieldShrunkEta[i]->GetEYhigh()[j];
                yErrorLowFinalEta[nPointFinalEta] = graphsCorrectedYieldShrunkEta[i]->GetEYlow()[j];
                if (sysAvailEta[i]){
                    Int_t counter = 0;
                    while(counter < 100 && TMath::Abs(xValueFinalEta[nPointFinalEta] - ptSysRelEta[i][counter])> 0.001) counter++;
                    if (counter < 100){
                        cout << ptSysRelEta[i][counter]<< "\t found it" << endl;
                        yErrorSysLowFinalEta[nPointFinalEta] = TMath::Abs(yErrorSysLowRelEta[i][counter]/100*graphsCorrectedYieldShrunkEta[i]->GetY()[j]);
                        yErrorSysHighFinalEta[nPointFinalEta] = yErrorSysHighRelEta[i][counter]/100*graphsCorrectedYieldShrunkEta[i]->GetY()[j];
                        
                    } else {
                        yErrorSysLowFinalEta[nPointFinalEta] = 0;
                        yErrorSysHighFinalEta[nPointFinalEta] = 0;
                        
                    }
                } else {
                    yErrorSysLowFinalEta[nPointFinalEta] = 0;
                    yErrorSysHighFinalEta[nPointFinalEta] = 0;
                }
                graphsCorrectedYieldSysShrunkEta[i]->SetPointEYlow(j,yErrorSysLowFinalEta[nPointFinalEta]);
                graphsCorrectedYieldSysShrunkEta[i]->SetPointEYhigh(j,yErrorSysHighFinalEta[nPointFinalEta]);
                nPointFinalEta++;
            }
            
            // fill inputs for combinations function in correct order
            if ( (triggerName[i].Contains("MB") || triggerName[i].Contains("INT1")) && graphsCorrectedYieldShrunkEta[i]){
                cout << "filling MB trigger" << endl;
                histoStatEta[0]     = histoCorrectedYieldEtaScaledMasked[i];
                graphSystEta[0]     = graphsCorrectedYieldSysShrunkEta[i];
                offSetsEtaSys[0]    = histoStatEta[0]->GetXaxis()->FindBin(graphSystEta[0]->GetX()[0])-1;
                if (graphMassEtaData[i]) 
                    graphOrderedMassEtaData[0]      = graphMassEtaData[i];
                if (graphMassEtaMC[i]) 
                    graphOrderedMassEtaMC[0]        = graphMassEtaMC[i];
                if (graphWidthEtaData[i]) 
                    graphOrderedWidthEtaData[0]     = graphWidthEtaData[i];
                if (graphWidthEtaMC[i]) 
                    graphOrderedWidthEtaMC[0]       = graphWidthEtaMC[i];
                if (graphAcceptanceEta[i])
                    graphOrderedAcceptanceEta[0]    = graphAcceptanceEta[i];
                if (graphEfficiencyEta[i])
                    graphOrderedEfficiencyEta[0]    = graphEfficiencyEta[i];
                if (graphEffTimesAccEta[i])
                    graphOrderedEffTimesAccEta[0]   = graphEffTimesAccEta[i];                
            } else if (triggerName[i].Contains("INT7") && graphsCorrectedYieldShrunkEta[i]){   
                cout << "filling INT7 trigger" << endl;
                histoStatEta[1]     = histoCorrectedYieldEtaScaledMasked[i];
                graphSystEta[1]     = graphsCorrectedYieldSysShrunkEta[i];
                offSetsEtaSys[1]    = histoStatEta[1]->GetXaxis()->FindBin(graphSystEta[1]->GetX()[0])-1;
                if (graphMassEtaData[i]) 
                    graphOrderedMassEtaData[1]      = graphMassEtaData[i];
                if (graphMassEtaMC[i]) 
                    graphOrderedMassEtaMC[1]        = graphMassEtaMC[i];
                if (graphWidthEtaData[i]) 
                    graphOrderedWidthEtaData[1]     = graphWidthEtaData[i];
                if (graphWidthEtaMC[i]) 
                    graphOrderedWidthEtaMC[1]       = graphWidthEtaMC[i];
                if (graphAcceptanceEta[i])
                    graphOrderedAcceptanceEta[1]    = graphAcceptanceEta[i];
                if (graphEfficiencyEta[i])
                    graphOrderedEfficiencyEta[1]    = graphEfficiencyEta[i];
                if (graphEffTimesAccEta[i])
                    graphOrderedEffTimesAccEta[1]   = graphEffTimesAccEta[i];                
            } else if (triggerName[i].Contains("EMC1") && graphsCorrectedYieldShrunkEta[i]){   
                cout << "filling EMC1 trigger" << endl;
                histoStatEta[2]     = histoCorrectedYieldEtaScaledMasked[i];
                graphSystEta[2]     = graphsCorrectedYieldSysShrunkEta[i];
                offSetsEtaSys[2]    = histoStatEta[2]->GetXaxis()->FindBin(graphSystEta[2]->GetX()[0])-1;
                if (graphMassEtaData[i]) 
                    graphOrderedMassEtaData[2]      = graphMassEtaData[i];
                if (graphMassEtaMC[i]) 
                    graphOrderedMassEtaMC[2]        = graphMassEtaMC[i];
                if (graphWidthEtaData[i]) 
                    graphOrderedWidthEtaData[2]     = graphWidthEtaData[i];
                if (graphWidthEtaMC[i]) 
                    graphOrderedWidthEtaMC[2]       = graphWidthEtaMC[i];
                if (graphAcceptanceEta[i])
                    graphOrderedAcceptanceEta[2]    = graphAcceptanceEta[i];
                if (graphEfficiencyEta[i])
                    graphOrderedEfficiencyEta[2]    = graphEfficiencyEta[i];
                if (graphEffTimesAccEta[i])
                    graphOrderedEffTimesAccEta[2]   = graphEffTimesAccEta[i];                
            } else if (triggerName[i].Contains("EMC7") && graphsCorrectedYieldShrunkEta[i]){   
                cout << "filling EMC7 trigger" << endl;
                histoStatEta[3]     = histoCorrectedYieldEtaScaledMasked[i];
                graphSystEta[3]     = graphsCorrectedYieldSysShrunkEta[i];
                offSetsEtaSys[3]    = histoStatEta[3]->GetXaxis()->FindBin(graphSystEta[3]->GetX()[0])-1;
                if (graphMassEtaData[i]) 
                    graphOrderedMassEtaData[3]      = graphMassEtaData[i];
                if (graphMassEtaMC[i]) 
                    graphOrderedMassEtaMC[3]        = graphMassEtaMC[i];
                if (graphWidthEtaData[i]) 
                    graphOrderedWidthEtaData[3]     = graphWidthEtaData[i];
                if (graphWidthEtaMC[i]) 
                    graphOrderedWidthEtaMC[3]       = graphWidthEtaMC[i];
                if (graphAcceptanceEta[i])
                    graphOrderedAcceptanceEta[3]    = graphAcceptanceEta[i];
                if (graphEfficiencyEta[i])
                    graphOrderedEfficiencyEta[3]    = graphEfficiencyEta[i];
                if (graphEffTimesAccEta[i])
                    graphOrderedEffTimesAccEta[3]   = graphEffTimesAccEta[i];                
            } else if ((triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA")) && graphsCorrectedYieldShrunkEta[i]){
                cout << Form("filling %s trigger",strEG2_A.Data()) << endl;
                histoStatEta[4]     = histoCorrectedYieldEtaScaledMasked[i];
                graphSystEta[4]     = graphsCorrectedYieldSysShrunkEta[i];
                offSetsEtaSys[4]    = histoStatEta[4]->GetXaxis()->FindBin(graphSystEta[4]->GetX()[0])-1;
                if(optionEnergy.CompareTo("8TeV")==0 && (mode == 4 || mode == 2)) offSetsEtaSys[4]+=2;

                if (graphMassEtaData[i]) 
                    graphOrderedMassEtaData[4]      = graphMassEtaData[i];
                if (graphMassEtaMC[i]) 
                    graphOrderedMassEtaMC[4]        = graphMassEtaMC[i];
                if (graphWidthEtaData[i]) 
                    graphOrderedWidthEtaData[4]     = graphWidthEtaData[i];
                if (graphWidthEtaMC[i]) 
                    graphOrderedWidthEtaMC[4]       = graphWidthEtaMC[i];
                if (graphAcceptanceEta[i])
                    graphOrderedAcceptanceEta[4]    = graphAcceptanceEta[i];
                if (graphEfficiencyEta[i])
                    graphOrderedEfficiencyEta[4]    = graphEfficiencyEta[i];
                if (graphEffTimesAccEta[i])
                    graphOrderedEffTimesAccEta[4]   = graphEffTimesAccEta[i];                
            } else if (triggerName[i].Contains("EG1") && graphsCorrectedYieldShrunkEta[i]){   
                cout << "filling EG1 trigger" << endl;
                histoStatEta[5]     = histoCorrectedYieldEtaScaledMasked[i];
                graphSystEta[5]     = graphsCorrectedYieldSysShrunkEta[i];
                offSetsEtaSys[5]    = histoStatEta[5]->GetXaxis()->FindBin(graphSystEta[5]->GetX()[0])-1;
                if (graphMassEtaData[i]) 
                    graphOrderedMassEtaData[5]      = graphMassEtaData[i];
                if (graphMassEtaMC[i]) 
                    graphOrderedMassEtaMC[5]        = graphMassEtaMC[i];
                if (graphWidthEtaData[i]) 
                    graphOrderedWidthEtaData[5]     = graphWidthEtaData[i];
                if (graphWidthEtaMC[i]) 
                    graphOrderedWidthEtaMC[5]       = graphWidthEtaMC[i];
                if (graphAcceptanceEta[i])
                    graphOrderedAcceptanceEta[5]    = graphAcceptanceEta[i];
                if (graphEfficiencyEta[i])
                    graphOrderedEfficiencyEta[5]    = graphEfficiencyEta[i];
                if (graphEffTimesAccEta[i])
                    graphOrderedEffTimesAccEta[5]   = graphEffTimesAccEta[i];                
            }
        }
        
        TString nameWeightsLogFileEta =     Form("%s/weightsEta_%s.dat",outputDir.Data(),isMC.Data());
        TGraphAsymmErrors* graphCorrectedYieldWeightedAverageEtaTot     = NULL;
        // Calculate averaged eta spectrum & respective supporting graphs according to statistical and systematic errors taking correctly into account the cross correlations
        if (averagedEta){
            // Calculate averaged eta spectrum
            graphCorrectedYieldWeightedAverageEtaTot        = CombinePtPointsSpectraTriggerCorrMat( histoStatEta, graphSystEta,  
                                                                                                    binningEta,  maxNAllowedEta,
                                                                                                    offSetsEta ,offSetsEtaSys,
                                                                                                    graphCorrectedYieldWeightedAverageEtaStat, graphCorrectedYieldWeightedAverageEtaSys,
                                                                                                    nameWeightsLogFileEta.Data(),
                                                                                                    mode, optionEnergy, "Eta", v2ClusterizerMerged,
                                                                                                    fileInputCorrFactors
                                                                                                  );
     //return;
            // Prepare arrays for reading weighting numbers 
            Double_t xValuesReadEta[100];
            Double_t weightsReadEta[12][100];
            Int_t availableMeasEta[12]  = { -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1};
            Int_t nMeasSetEta           = nrOfTrigToBeCombEtaRed;
            Int_t nPtBinsReadEta        = 0;
            
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
            ifstream fileWeightsEta;
            fileWeightsEta.open(nameWeightsLogFileEta,ios_base::in);
            cout << "reading" << nameWeightsLogFileEta << endl;
            
            while(!fileWeightsEta.eof() && nPtBinsReadEta < 100){
                TString garbage = "";
                if (nPtBinsReadEta == 0){
                    fileWeightsEta >> garbage ;
                    for (Int_t i = 0; i < nMeasSetEta; i++){
                        fileWeightsEta >> availableMeasEta[i] ;
                    }   
                    cout << "read following measurements: "; 
                    for (Int_t i = 0; i < nMeasSetEta; i++){
                        cout << availableMeasEta[i] << "\t" ;
                    }   
                    cout << endl;
                } else {
                    fileWeightsEta >> xValuesReadEta[nPtBinsReadEta-1];
                    for (Int_t i = 0; i < nMeasSetEta; i++){
                        fileWeightsEta >> weightsReadEta[availableMeasEta[i]][nPtBinsReadEta-1] ;
                    }   
                    cout << "read: "<<  nPtBinsReadEta << "\t"<< xValuesReadEta[nPtBinsReadEta-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetEta; i++){
                        cout << weightsReadEta[availableMeasEta[i]][nPtBinsReadEta-1] << "\t";
                    }
                    cout << endl;
                }
                nPtBinsReadEta++;
            }
            nPtBinsReadEta = nPtBinsReadEta-2 ;
            fileWeightsEta.close();

            // create and fill weighting graphs
            TGraph* graphWeightsEta[12];
            for (Int_t i = 0; i < numberOfTrigg; i++){
                graphWeightsEta[i]                    = NULL;
            }  
            for (Int_t i = 0; i < nMeasSetEta; i++){
                graphWeightsEta[availableMeasEta[i]]  = new TGraph(nPtBinsReadEta,xValuesReadEta,weightsReadEta[availableMeasEta[i]]);
                Int_t bin = 0;
                for (Int_t n = 0; n< nPtBinsReadEta; n++){
                    if (graphWeightsEta[availableMeasEta[i]]->GetY()[bin] == 0) graphWeightsEta[availableMeasEta[i]]->RemovePoint(bin);
                    else bin++;
                }   
            }   

            //  **********************************************************************************************************************
            //  **************************************** Combine+write detailed Systematics ******************************************
            //  **********************************************************************************************************************

            const char *SysErrDatnameMeanSingleErr = Form("%s/SystematicErrorAveragedSingle%s_Eta_%s.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
            fstream SysErrDatAverSingle;
            SysErrDatAverSingle.precision(4);
            cout << SysErrDatnameMeanSingleErr << endl;
            if(sysAvailSingleEta[0]){
              SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
              for(Int_t iColumn = 0; iColumn < (Int_t)ptSysDetail[0][0].size(); iColumn++) SysErrDatAverSingle << ptSysDetail[0][0].at(iColumn) << "\t";
              SysErrDatAverSingle << endl;
              for (Int_t i = 1; i < nrOfTrigToBeComb; i++){
                if(!sysAvailSingleEta[i]) continue;
                for(Int_t iCol = 0; iCol < (Int_t)ptSysDetail[i][0].size(); iCol++){
                  if( ((TString)ptSysDetail[i][0].at(iCol)).CompareTo(((TString)ptSysDetail[0][0].at(iCol))) ){
                    cout << "ERROR: Systematic error type at pos " << iCol << " does not agree for " << availableMeasEta[i] << " & " << availableMeasEta[0] << ", returning!" << endl;
                    return;
                  }
                }
              }

              for(Int_t i=0; i<nPtBinsReadEta; i++){
                SysErrDatAverSingle << xValuesReadEta[i] << "\t";
                Int_t nColumns = (Int_t)ptSysDetail[0][0].size();
                Double_t *errors = new Double_t[nColumns-1];
                for(Int_t iErr=0; iErr<nColumns-1; iErr++) errors[iErr] = 0;
                for(Int_t j=0; j<nrOfTrigToBeComb; j++){
                  if(!sysAvailSingleEta[j]) continue;
                  Int_t pos = GetBinPosInVec(ptSysDetail[j],numberBinsSysAvailSingleEta[j],xValuesReadEta[i]);
                  if(pos>-1){
                    for(Int_t iErr=1; iErr<nColumns; iErr++) errors[iErr-1] += weightsReadEta[availableMeasEta[j]][i]*((TString)ptSysDetail[j][pos].at(iErr)).Atof();
                  }
                }
                for(Int_t iErr=0; iErr<nColumns-1; iErr++) SysErrDatAverSingle << errors[iErr] << "\t";
                SysErrDatAverSingle << endl;
                delete[] errors;
              }
            }
            SysErrDatAverSingle.close();

            for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++){
              for(Int_t iB=0; iB<50; iB++) ptSysDetail[iR][iB].clear();
            }

            //  **********************************************************************************************************************
            //  ******************************************* Plotting weights method for eta ******************************************
            //  **********************************************************************************************************************
            Int_t textSizeLabelsPixel = 900*0.04;

            TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);// gives the page size
            DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
        
            TH2F * histo2DWeightsEta;
            histo2DWeightsEta = new TH2F("histo2DWeightsEta","histo2DWeightsEta",11000,0.,maxPtGlobalEta,1000,-0.5,1.3);
            SetStyleHistoTH2ForGraphs(histo2DWeightsEta, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
            histo2DWeightsEta->Draw("copy");
            
                TLegend* legendWeightsEta = GetAndSetLegend2(0.12, 0.14, 0.55, 0.14+(0.035*nMeasSetEta/2*1.35), 32);
                legendWeightsEta->SetNColumns(2);
                for (Int_t i = 0; i < nMeasSetEta; i++){
                    DrawGammaSetMarkerTGraph(graphWeightsEta[availableMeasEta[i]],  markerTriggWeighted[availableMeasEta[i]], sizeTrigg[availableMeasEta[i]], 
                                            colorTriggWeighted[availableMeasEta[i]], colorTriggWeighted[availableMeasEta[i]]);
                    graphWeightsEta[availableMeasEta[i]]->Draw("p,same,e1");
                    legendWeightsEta->AddEntry(graphWeightsEta[availableMeasEta[i]],nameTriggerWeighted[availableMeasEta[i]],"p");
                }   
                legendWeightsEta->Draw();

                TLatex *labelWeightsEnergy = new TLatex(0.7,0.24,collisionSystem.Data());
                SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
                labelWeightsEnergy->SetTextFont(43);
                labelWeightsEnergy->Draw();
                TLatex *labelWeightsEta = new TLatex(0.7,0.20,"#eta #rightarrow #gamma#gamma");
                SetStyleTLatex( labelWeightsEta, 0.85*textSizeLabelsPixel,4);
                labelWeightsEta->SetTextFont(43);
                labelWeightsEta->Draw();
                TLatex *labelDetProcWeights    = new TLatex(0.7, 0.16,detectionProcess.Data());
                SetStyleTLatex( labelDetProcWeights, 0.85*textSizeLabelsPixel,4);
                labelDetProcWeights->SetTextFont(43);
                labelDetProcWeights->Draw();
                
        //      DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
                DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
                DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
                DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
                DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
                
            canvasWeights->SaveAs(Form("%s/%s_WeightsEtaTriggers.%s", outputDir.Data(), isMC.Data(), suffix.Data()));
            delete canvasWeights;

            // Calculating relative error for eta        
            for (Int_t i = 0; i < 12; i++){
                if (histoStatEta[i]) 
                    histoRelStatEta[i]      = CalculateRelErrUpTH1D( histoStatEta[i], Form("relativeStatErrorEta_%s", nameTriggerWeighted[i].Data()));
                if (graphSystEta[i]) 
                    graphRelSystEta[i]      = CalculateRelErrUpAsymmGraph( graphSystEta[i], Form("relativeSysErrorEta_%s", nameTriggerWeighted[i].Data()));
            }
            TGraphAsymmErrors* graphRelErrorEtaTot        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaTot, "relativeTotalErrorEta");
            TGraphAsymmErrors* graphRelErrorEtaStat       = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaStat, "relativeStatErrorEta");
            TGraphAsymmErrors* graphRelErrorEtaSys        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaSys, "relativeSysErrorEta");

            // plot sys relative errors for individual triggers   
            TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
        
            TH2F * histo2DRelSysErr;
            histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.,maxPtGlobalEta,1000,0,60.5);
            SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
            histo2DRelSysErr->Draw("copy");
                TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEta/2), 0.95, 0.92, 32);
                legendRelSysErr->SetNColumns(2);
                for (Int_t i = 0; i < nMeasSetEta; i++){
                    DrawGammaSetMarkerTGraph(graphRelSystEta[availableMeasEta[i]], markerTriggWeighted[availableMeasEta[i]], sizeTrigg[availableMeasEta[i]], 
                                            colorTriggWeighted[availableMeasEta[i]], colorTriggWeighted[availableMeasEta[i]]);
                    graphRelSystEta[availableMeasEta[i]]->Draw("p,same,z");
                    legendRelSysErr->AddEntry(graphRelSystEta[availableMeasEta[i]],nameTriggerWeighted[availableMeasEta[i]],"p");
                }    
                legendRelSysErr->Draw();

                TLatex *labelRelErrEnergy    = new TLatex(0.15,0.89,collisionSystem.Data());
                SetStyleTLatex( labelRelErrEnergy, 0.85*textSizeLabelsPixel,4);
                labelRelErrEnergy->SetTextFont(43);
                labelRelErrEnergy->Draw();
                TLatex *labelRelErrEta       = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
                SetStyleTLatex( labelRelErrEta, 0.85*textSizeLabelsPixel,4);
                labelRelErrEta->SetTextFont(43);
                labelRelErrEta->Draw();
                TLatex *labelDetProcRelErr    = new TLatex(0.15, 0.81,detectionProcess.Data());
                SetStyleTLatex( labelDetProcRelErr, 0.85*textSizeLabelsPixel,4);
                labelDetProcRelErr->SetTextFont(43);
                labelDetProcRelErr->Draw();
                
            canvasRelSysErr->SaveAs(Form("%s/Eta_RelSysErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));
            
            // plot stat relative errors for individual triggers    
            TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
        
            TH2F * histo2DRelStatErr;
            histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.,maxPtGlobalEta,1000,0,60.5);
            SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
            histo2DRelStatErr->Draw("copy");
                TLegend* legendRelStatErr       = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEta/2), 0.95, 0.92, 32);
                legendRelStatErr->SetNColumns(2);
                for (Int_t i = 0; i < nMeasSetEta; i++){
    //                 cout << "plotting graph: " << availableMeasEta[i] << "\t" <<histoRelStatEta[availableMeasEta[i]]->GetName() << endl;
                    if (availableMeasEta[i] == 0 && histoRelStatEta[availableMeasEta[i]] && mode == 2){
                        TGraphAsymmErrors* dummyGraph = new TGraphAsymmErrors(histoRelStatEta[availableMeasEta[i]]);
                        dummyGraph->Print();
                        DrawGammaSetMarkerTGraph(dummyGraph, markerTriggWeighted[availableMeasEta[i]], sizeTrigg[availableMeasEta[i]], 
                                            colorTriggWeighted[availableMeasEta[i]], colorTriggWeighted[availableMeasEta[i]]);
                        dummyGraph->Draw("pX,same");
                        legendRelSysErr->AddEntry(dummyGraph,nameTriggerWeighted[availableMeasEta[i]],"p");
                        
    //                     for (Int_t j = 1; j < histoRelStatEta[availableMeasEta[i]]->GetNbinsX()+1; j++){
    //                        cout << j << ": " << histoRelStatEta[availableMeasEta[i]]->GetBinContent(j) << endl;
    //                     }    
                    } else {
                        DrawGammaSetMarker(histoRelStatEta[availableMeasEta[i]],markerTriggWeighted[availableMeasEta[i]], sizeTrigg[availableMeasEta[i]], 
                                            colorTriggWeighted[availableMeasEta[i]], colorTriggWeighted[availableMeasEta[i]]);
                        histoRelStatEta[availableMeasEta[i]]->Draw("p,same,z");
                        legendRelStatErr->AddEntry(histoRelStatEta[availableMeasEta[i]],nameTriggerWeighted[availableMeasEta[i]],"p");
                    }    
                }    
                legendRelStatErr->Draw();

                labelRelErrEnergy->Draw();
                labelRelErrEta->Draw();
                labelDetProcRelErr->Draw();
                
            canvasRelStatErr->SaveAs(Form("%s/Eta_RelStatErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));
            
            // plot full error for final result decomposed
            TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
        
            TH2F * histo2DRelTotErrEta;
            histo2DRelTotErrEta                 = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,0.,maxPtGlobalEta,1000,0,60.5);
            SetStyleHistoTH2ForGraphs(histo2DRelTotErrEta, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
            histo2DRelTotErrEta->Draw("copy");

                DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaTot, 20, 1.5, kBlack , kBlack);
                graphRelErrorEtaTot->Draw("p,same,z");
                DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaStat, 24, 1.5, kGray+2 , kGray+2);
                graphRelErrorEtaStat->Draw("l,x0,same,e1");
                DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaSys, 24, 1.5, kGray+1 , kGray+1);
                graphRelErrorEtaSys->SetLineStyle(7);
                graphRelErrorEtaSys->Draw("l,x0,same,e1");

                TLegend* legendRelTotErr2       = GetAndSetLegend2(0.72, 0.92-(0.035*3), 0.9, 0.92, 32);
                legendRelTotErr2->AddEntry(graphRelErrorEtaTot,"tot","p");
                legendRelTotErr2->AddEntry(graphRelErrorEtaStat,"stat","l");
                legendRelTotErr2->AddEntry(graphRelErrorEtaSys,"sys","l");
                legendRelTotErr2->Draw();

                labelRelErrEnergy->Draw();
                labelRelErrEta->Draw();
                labelDetProcRelErr->Draw();
                
            canvasRelTotErr->SaveAs(Form("%s/Eta_RelErrorsFulldecomp.%s",outputDir.Data(),suffix.Data()));
            

            // create averaged supporting graphs
            graphMassEtaDataWeighted                        = CalculateWeightedQuantity(    graphOrderedMassEtaData, 
                                                                                            graphWeightsEta,
                                                                                            binningEta,  maxNAllowedEta,
                                                                                            MaxNumberOfFiles
                                                                                        );
            
            if (!graphMassEtaDataWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            graphMassEtaMCWeighted                          = CalculateWeightedQuantity(    graphOrderedMassEtaMC, 
                                                                                            graphWeightsEta,
                                                                                            binningEta,  maxNAllowedEta,
                                                                                            MaxNumberOfFiles
                                                                                        );
            if (!graphMassEtaMCWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            graphWidthEtaDataWeighted                       = CalculateWeightedQuantity(    graphOrderedWidthEtaData, 
                                                                                            graphWeightsEta,
                                                                                            binningEta,  maxNAllowedEta,
                                                                                            MaxNumberOfFiles
                                                                                        );
            if (!graphWidthEtaDataWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            graphWidthEtaMCWeighted                         = CalculateWeightedQuantity(    graphOrderedWidthEtaMC, 
                                                                                            graphWeightsEta,
                                                                                            binningEta,  maxNAllowedEta,
                                                                                            MaxNumberOfFiles
                                                                                        );
            if (!graphWidthEtaMCWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            graphAcceptanceEtaWeighted                      = CalculateWeightedQuantity(    graphOrderedAcceptanceEta, 
                                                                                            graphWeightsEta,
                                                                                            binningEta,  maxNAllowedEta,
                                                                                            MaxNumberOfFiles
                                                                                        );
            if (!graphAcceptanceEtaWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            graphEfficiencyEtaWeighted                      = CalculateWeightedQuantity(    graphOrderedEfficiencyEta, 
                                                                                            graphWeightsEta,
                                                                                            binningEta,  maxNAllowedEta,
                                                                                            MaxNumberOfFiles
                                                                                        );
            if (!graphEfficiencyEtaWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            graphEffTimesAccEtaWeighted                     = CalculateWeightedQuantity(    graphOrderedEffTimesAccEta, 
                                                                                            graphWeightsEta,
                                                                                            binningEta,  maxNAllowedEta,
                                                                                            MaxNumberOfFiles
                                                                                        );
            if (!graphEffTimesAccEtaWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            // remove masked bins at beginning
            while (graphMassEtaDataWeighted->GetY()[0] == -10000)    graphMassEtaDataWeighted->RemovePoint(0);
            while (graphMassEtaMCWeighted->GetY()[0] == -10000)      graphMassEtaMCWeighted->RemovePoint(0);
            while (graphWidthEtaDataWeighted->GetY()[0] == -10000)   graphWidthEtaDataWeighted->RemovePoint(0);
            while (graphWidthEtaMCWeighted->GetY()[0] == -10000)     graphWidthEtaMCWeighted->RemovePoint(0);
            while (graphAcceptanceEtaWeighted->GetY()[0] == -10000)  graphAcceptanceEtaWeighted->RemovePoint(0);
            while (graphEfficiencyEtaWeighted->GetY()[0] == -10000)  graphEfficiencyEtaWeighted->RemovePoint(0);
            while (graphEffTimesAccEtaWeighted->GetY()[0] == -10000) graphEffTimesAccEtaWeighted->RemovePoint(0);
            
        // if averaging wasn't enabled pick values according to predefined ranges ("cherry picking points")        
        } else {
            graphCorrectedYieldWeightedAverageEtaStat        = new TGraphAsymmErrors(nPointFinalEta, xValueFinalEta, yValueFinalEta, 
                                                                                        xErrorLowFinalEta, xErrorHighFinalEta,yErrorLowFinalEta, yErrorHighFinalEta);
            graphCorrectedYieldWeightedAverageEtaSys         = new TGraphAsymmErrors(nPointFinalEta, xValueFinalEta, yValueFinalEta, 
                                                                                    xErrorLowFinalEta, xErrorHighFinalEta,yErrorSysLowFinalEta, yErrorSysHighFinalEta);
        }   
        // Printing final graphs
        cout << "stat eta" << endl; 
        if (graphCorrectedYieldWeightedAverageEtaStat) graphCorrectedYieldWeightedAverageEtaStat->Print();
        cout << "sys eta" << endl;
        if (graphCorrectedYieldWeightedAverageEtaSys) graphCorrectedYieldWeightedAverageEtaSys->Print();

        if (graphCorrectedYieldWeightedAverageEtaStat){
            while (graphCorrectedYieldWeightedAverageEtaStat->GetX()[0] < minPtGlobalEta){
                graphCorrectedYieldWeightedAverageEtaStat->RemovePoint(0);
            }    
        }
        
        if (graphCorrectedYieldWeightedAverageEtaSys){
            while (graphCorrectedYieldWeightedAverageEtaSys->GetX()[0]< minPtGlobalEta){
                graphCorrectedYieldWeightedAverageEtaSys->RemovePoint(0);
            }    
        }

        
        //***************************************************************************************************************
        //************************************Plotting Mass Eta reduced range  ******************************************
        //***************************************************************************************************************
        canvasMass->cd();
        histo2DMassEta->DrawCopy(); 

        TLegend* legendMassRedEta = GetAndSetLegend2(0.52, 0.88, 0.95, 0.88+(1.05*4/2*0.85*textSizeSpectra),28);
        legendMassRedEta->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0)
             ){
                if (graphMassEtaData[i] && !maskedFullyEta[i]) {
                    DrawGammaSetMarkerTGraphAsym(graphMassEtaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    graphMassEtaData[i]->Draw("p,e1,same"); 
                    legendMassRedEta->AddEntry(graphMassEtaData[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                }
                if (graphMassEtaMC[i] && !maskedFullyEta[i]){
                    DrawGammaSetMarkerTGraphAsym(graphMassEtaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    graphMassEtaMC[i]->Draw("p,e1,same"); 
                    legendMassRedEta->AddEntry(graphMassEtaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                }    
            }    
        }
        legendMassRedEta->Draw();
        labelEnergyMass->Draw();
        labelEtaMass->Draw();
        labelDetProcMass->Draw();

        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/Eta_%s_Mass1_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        canvasMass->cd();
        histo2DMassEta->DrawCopy();    
        TLegend* legendMassRedEta2 = GetAndSetLegend2(0.46, 0.81, 0.80, 0.81+(1.05*8/2*0.85*textSizeSpectra),28);
        legendMassRedEta2->SetNColumns(2);        
        if (optionEnergy.CompareTo("2.76TeV")==0){
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if(i==1 || i==3 || i==4 || i==5 ){
                    if (graphMassEtaData[i] && !maskedFullyEta[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphMassEtaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                        graphMassEtaData[i]->Draw("p,e1,same"); 
                        legendMassRedEta2->AddEntry(graphMassEtaData[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                    }
                    if (graphMassEtaMC[i] && !maskedFullyEta[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphMassEtaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                        graphMassEtaMC[i]->Draw("p,e1,same"); 
                        legendMassRedEta2->AddEntry(graphMassEtaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                    }    
                }    
            }
            legendMassRedEta2->Draw();
            labelEnergyMass->Draw();
            labelEtaMass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Eta_%s_Mass2_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        //***************************************************************************************************************
        //********************************* Eta Mass weighted ***********************************************************
        //***************************************************************************************************************    
        if (graphMassEtaDataWeighted || graphMassEtaMCWeighted){
            canvasMass->cd();
            histo2DMassEta->DrawCopy();    
            TLegend* legendMassEtaWeighted = GetAndSetLegend2(0.52, 0.88, 0.75, 0.88+(1.05*2*0.85*textSizeSpectra),28);
            if (graphMassEtaDataWeighted){
                DrawGammaSetMarkerTGraphAsym(graphMassEtaDataWeighted, 20, 1, kBlack, kBlack);
                graphMassEtaDataWeighted->Draw("p,e1,same");                 
                legendMassEtaWeighted->AddEntry(graphMassEtaDataWeighted, "Data", "p"); 
            }    
            if (graphMassEtaDataWeighted){
                DrawGammaSetMarkerTGraphAsym(graphMassEtaMCWeighted, 24, 1, kGray+2, kGray+2);
                graphMassEtaMCWeighted->Draw("p,e1,same");                 
                legendMassEtaWeighted->AddEntry(graphMassEtaMCWeighted, "MC", "p"); 
            }    
            legendMassEtaWeighted->Draw();
            labelEnergyMass->Draw();
            labelEtaMass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Eta_%s_Mass_Weighted.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        delete canvasMass;
        //***************************************************************************************************************
        //************************************* Efficiency weighted *****************************************************
        //***************************************************************************************************************    
    
        if (graphEfficiencyEtaWeighted){
            canvasEffi->SetLogy(1);
            canvasEffi->cd();
            histo2DEffiEta->DrawCopy(); 

            DrawGammaSetMarkerTGraphAsym(graphEfficiencyEtaWeighted, 20, 1, kGray+2, kGray+2);
            graphEfficiencyEtaWeighted->Draw("p,e1,same");                 

            TLatex *labelEnergyEffiWOTrigg = new TLatex(0.62, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
            SetStyleTLatex( labelEnergyEffiWOTrigg, 0.85*textSizeSpectra,4);
            labelEnergyEffiWOTrigg->Draw();

            TLatex *labelEtaEffiWOTrigg = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#eta #rightarrow #gamma#gamma");
            SetStyleTLatex( labelEtaEffiWOTrigg, 0.85*textSizeSpectra,4);
            labelEtaEffiWOTrigg->Draw();

            TLatex *labelDetProcEffiWOTrigg = new TLatex(0.62, 0.15,detectionProcess.Data());
            SetStyleTLatex( labelDetProcEffiWOTrigg, 0.85*textSizeSpectra,4);
            labelDetProcEffiWOTrigg->Draw();
            
            canvasEffi->Update();
            canvasEffi->SaveAs(Form("%s/Eta_EfficiencyW0TriggEff_Weighted.%s",outputDir.Data(),suffix.Data()));
        }
        delete histo2DEffiEta;
        delete canvasEffi;
        
        //***************************************************************************************************************
        //************************************* Acceptance weighted *****************************************************
        //***************************************************************************************************************    
        if (graphAcceptanceEtaWeighted){
            canvasAcc->cd();

            histo2DAccEta->DrawCopy(); 

            DrawGammaSetMarkerTGraphAsym(graphAcceptanceEtaWeighted, 20, 1, kGray+2, kGray+2);
            graphAcceptanceEtaWeighted->Draw("p,e1,same");                 
            
            TLatex *labelEnergyAcc = new TLatex(0.62, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
            SetStyleTLatex( labelEnergyAcc, 0.85*textSizeSpectra,4);
            labelEnergyAcc->Draw();

            TLatex *labelEtaAcc = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#eta #rightarrow #gamma#gamma");
            SetStyleTLatex( labelEtaAcc, 0.85*textSizeSpectra,4);
            labelEtaAcc->Draw();

            TLatex *labelDetProcAcc = new TLatex(0.62, 0.15,detectionProcess.Data());
            SetStyleTLatex( labelDetProcAcc, 0.85*textSizeSpectra,4);
            labelDetProcAcc->Draw();

            canvasAcc->Update();
            canvasAcc->SaveAs(Form("%s/Eta_Acceptance_weighted.%s",outputDir.Data(),suffix.Data()));
        }    
        delete canvasAcc;
        

        //***************************************************************************************************************
        //************************************Plotting Width Eta reduced range  *****************************************
        //***************************************************************************************************************
        canvasWidth->cd();
        histo2DWidthEta->DrawCopy(); 

        TLegend* legendWidthRedEta = GetAndSetLegend2(0.52, 0.87, 0.95, 0.87+(1.05*4/2*0.85*textSizeSpectra),28);
        legendWidthRedEta->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0)
             ){
                if (graphWidthEtaData[i] && !maskedFullyEta[i]) {
                    DrawGammaSetMarkerTGraphAsym(graphWidthEtaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    graphWidthEtaData[i]->Draw("p,e1,same"); 
                    legendWidthRedEta->AddEntry(graphWidthEtaData[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                } 
                if (graphWidthEtaMC[i] && !maskedFullyEta[i]){
                    DrawGammaSetMarkerTGraphAsym(graphWidthEtaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    graphWidthEtaMC[i]->Draw("p,e1,same"); 
                    legendWidthRedEta->AddEntry(graphWidthEtaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                }    
            }    
        }
        legendWidthRedEta->Draw();
        labelEnergyWidth->Draw();
        labelEtaWidth->Draw();
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Eta_%s_Width1_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        canvasWidth->cd();
        histo2DWidthEta->DrawCopy();    
        TLegend* legendWidthRedEta2 = GetAndSetLegend2(0.46, 0.80, 0.80, 0.80+(1.05*8/2*0.85*textSizeSpectra),28);
        legendWidthRedEta2->SetNColumns(2);        
        if (optionEnergy.CompareTo("2.76TeV")==0 ){
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if( i==1 || i==3 || i==4 || i==5 ){

                    if (graphWidthEtaData[i] && !maskedFullyEta[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphWidthEtaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                        graphWidthEtaData[i]->Draw("p,e1,same"); 
                        legendWidthRedEta2->AddEntry(graphWidthEtaData[i], Form("%s data",triggerNameLabel[i].Data()), "p"); 
                    }    
                    if (graphWidthEtaMC[i] && !maskedFullyEta[i]) {    
                        DrawGammaSetMarkerTGraphAsym(graphWidthEtaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                        graphWidthEtaMC[i]->Draw("p,e1,same"); 
                        legendWidthRedEta2->AddEntry(graphWidthEtaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p"); 
                    }    
                }    
            }
            legendWidthRedEta2->Draw();
            labelEnergyWidth->Draw();
            labelEtaWidth->Draw();
            labelDetProcWidth->Draw();

            canvasWidth->Update();
            canvasWidth->SaveAs(Form("%s/Eta_%s_Width2_Reduced.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        
        //***************************************************************************************************************
        //********************************* Eta Width weighted **********************************************************
        //***************************************************************************************************************        
        canvasWidth->cd();
        histo2DWidthEta->DrawCopy();    
        TLegend* legendWidthEtaWeighted = GetAndSetLegend2(0.52, 0.86, 0.75, 0.86+(1.05*2*0.85*textSizeSpectra),28);
        if (graphWidthEtaDataWeighted){
            DrawGammaSetMarkerTGraphAsym(graphWidthEtaDataWeighted, 20, 1, kBlack, kBlack);
            graphWidthEtaDataWeighted->Draw("p,e1,same");                 
            legendWidthEtaWeighted->AddEntry(graphWidthEtaDataWeighted, "Data", "p"); 
        }    
        if (graphWidthEtaDataWeighted){
            DrawGammaSetMarkerTGraphAsym(graphWidthEtaMCWeighted, 24, 1, kGray+2, kGray+2);
            graphWidthEtaMCWeighted->Draw("p,e1,same");                 
            legendWidthEtaWeighted->AddEntry(graphWidthEtaMCWeighted, "MC", "p"); 
        }    
        legendWidthEtaWeighted->Draw();
        labelEnergyWidth->Draw();
        labelEtaWidth->Draw();
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Eta_%s_Width_Weighted.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        delete canvasWidth;
        
        //***************************************************************************************************************
        //************************************Plotting scaled invariant yield *****************************************
        //***************************************************************************************************************
        canvasCorrScaled->cd();

        Double_t minCorrYieldEta    = 2e-11;
        Double_t maxCorrYieldEta    = 5e-2;
        if (mode == 0){
            minCorrYieldEta         = 5e-8;
            maxCorrYieldEta         = 1e-1;                
        } 

        if(optionEnergy.CompareTo("8TeV")==0){
          if(mode == 2){
            minCorrYieldEta         = 2e-11;
            maxCorrYieldEta         = 5e-2;
          }else if(mode == 4){
            minCorrYieldEta     = 3e-11;
            maxCorrYieldEta     = 7e-3;
          }
        }
        
        TH2F * histo2DInvYieldScaledEta;
        histo2DInvYieldScaledEta = new TH2F("histo2DInvYieldScaledEta","histo2DInvYieldScaledEta",1000,0., maxPtGlobalEta,10000,minCorrYieldEta,maxCorrYieldEta);
        SetStyleHistoTH2ForGraphs(histo2DInvYieldScaledEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.55);
//         histo2DInvYieldScaledEta->SetRangeUser(2e-11,5e-2);
        histo2DInvYieldScaledEta->DrawCopy(); 

        // two component model fit
        cout << "first point eta" << endl;
        graphCorrectedYieldWeightedAverageEtaStat->Print();
        cout << graphCorrectedYieldWeightedAverageEtaStat->GetY()[0] << endl;
        Double_t paramTCM[5] = {graphCorrectedYieldWeightedAverageEtaStat->GetY()[0],0.3,graphCorrectedYieldWeightedAverageEtaStat->GetY()[0]/1000,0.8,3};
        TF1* fitInvYieldEta = FitObject("tcm","fitInvYieldEta","Eta",graphCorrectedYieldWeightedAverageEtaStat,minPtGlobalEta,maxPtGlobalEta,paramTCM,"QNRMEX0+");
        cout << WriteParameterToFile(fitInvYieldEta) << endl;
        
        // tsallis fit
    //     Double_t paramGraph[3]                  = {1000, 8., 0.13};
    //     TF1* fitInvYieldEta                     = FitObject("l","fitInvYieldEta","Eta",graphCorrectedYieldWeightedAverageEtaStat,minPtGlobalEta,maxPtGlobalEta,paramGraph,"QNRME+");

        DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageEtaSys, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
        graphCorrectedYieldWeightedAverageEtaSys->Draw("p,E2,same");

        TLegend* legendScaled = GetAndSetLegend2(0.72, 0.95-(1.15*nrOfTrigToBeCombEtaRed*0.85*textSizeSpectra), 0.95, 0.95,32);
        
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){        
            if (graphsCorrectedYieldSysShrunkEta[i] && !maskedFullyEta[i]) DrawGammaSetMarkerTGraphAsym(graphsCorrectedYieldSysShrunkEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
            if (graphsCorrectedYieldSysShrunkEta[i] && !maskedFullyEta[i]) DrawGammaSetMarker(histoCorrectedYieldEtaScaled[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            if (graphsCorrectedYieldSysShrunkEta[i] && !maskedFullyEta[i]) histoCorrectedYieldEtaScaled[i]->DrawCopy("e1,same"); 
            if (graphsCorrectedYieldSysShrunkEta[i] && !maskedFullyEta[i]) graphsCorrectedYieldSysShrunkEta[i]->Draw("p,E2,same");
            if (graphsCorrectedYieldSysShrunkEta[i] && !maskedFullyEta[i]) legendScaled->AddEntry(histoCorrectedYieldEtaScaled[i],triggerNameLabel[i].Data(),"p");
        }
        
        DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageEtaStat, 24, 2, kRed , kRed, 1, kTRUE);
        graphCorrectedYieldWeightedAverageEtaStat->Draw("p,E,same");
        legendScaled->AddEntry(graphCorrectedYieldWeightedAverageEtaStat,"Final","p");
        legendScaled->Draw();
        
        fitInvYieldEta->SetLineColor(kGray+2);
        fitInvYieldEta->SetLineStyle(7);
        fitInvYieldEta->SetLineWidth(2);
        fitInvYieldEta->Draw("same");

        labelEnergyUnscaled->Draw();
        labelEtaUnscaled->Draw();
        labelDetProcUnscaled->Draw();
        
        canvasCorrScaled->Update();
        canvasCorrScaled->SaveAs(Form("%s/Eta_%s_CorrectedYieldScaledTrigg.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        //***************************************************************************************************************
        //************************************Plotting final invariant yield *****************************************
        //***************************************************************************************************************

        histo2DInvYieldScaledEta->DrawCopy(); 

        graphCorrectedYieldWeightedAverageEtaSys->Draw("p,E2,same");        
        DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageEtaStat, 24, 2, kBlack , kBlack, 1, kTRUE);
        graphCorrectedYieldWeightedAverageEtaStat->Draw("p,E,same");
        
        fitInvYieldEta->Draw("same");

        labelEnergyUnscaled->Draw();
        labelEtaUnscaled->Draw();
        labelDetProcUnscaled->Draw();
        
        canvasCorrScaled->Update();
        canvasCorrScaled->SaveAs(Form("%s/Eta_%s_CorrectedYieldFinal.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        
        delete canvasCorrScaled;

        //***************************************************************************************************************
        //******************************* Ratio to fit for individual spectra full range ********************************
        //***************************************************************************************************************
        canvasRatioSpec->cd();
        
        TH2F * histo2DRatioToFitEta;
        histo2DRatioToFitEta = new TH2F("histo2DRatioToFitEta","histo2DRatioToFitEta",1000,0., maxPtGlobalEta,1000,0.25, 2.05);
        SetStyleHistoTH2ForGraphs(histo2DRatioToFitEta, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.0);
        histo2DRatioToFitEta->DrawCopy(); 

        TH1D* histoCorrectedYieldToFitEta[MaxNumberOfFiles];
        TGraphAsymmErrors* graphCorrectedYieldToFitEta[MaxNumberOfFiles];
        TLegend* legendRatioSpecEta = GetAndSetLegend2(0.15, 0.95-(1.1*nrOfTrigToBeCombEtaRed/2*0.85*textSizeSpectra), 0.45, 0.95,28);
        legendRatioSpecEta->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (graphsCorrectedYieldSysRemoved0Eta[i] && !maskedFullyEta[i]){
              histoCorrectedYieldToFitEta[i] = CalculateHistoRatioToFit (histoCorrectedYieldEtaScaled[i], fitInvYieldEta); 
              graphCorrectedYieldToFitEta[i] = CalculateGraphErrRatioToFit(graphsCorrectedYieldSysRemoved0Eta[i], fitInvYieldEta); 
              
              DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
              DrawGammaSetMarker(histoCorrectedYieldToFitEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);

              graphCorrectedYieldToFitEta[i]->Draw("p,E2,same");
              legendRatioSpecEta->AddEntry(histoCorrectedYieldToFitEta[i],triggerNameLabel[i].Data(),"p");
            }
            DrawGammaLines(0., maxPtGlobalEta , 1., 1.,0.1, kGray+2);
            DrawGammaLines(0., maxPtGlobalEta , 1.1, 1.1,0.1, kGray, 7);
            DrawGammaLines(0., maxPtGlobalEta , 0.9, 0.9,0.1, kGray, 7);

            if (graphsCorrectedYieldSysRemoved0Eta[i]) histoCorrectedYieldToFitEta[i]->DrawCopy("e1,same"); 
        }
        legendRatioSpecEta->Draw();

        labelEnergyRatio->Draw();
        
        TLatex *labelEtaRatio = new TLatex(0.6, 0.93-textSizeSpectra*0.85*1.04, "#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelEtaRatio, 0.85*textSizeSpectra,4);
        labelEtaRatio->Draw();
        
        labelDetProcRatio->Draw();
        
        canvasRatioSpec->Update();
        canvasRatioSpec->SaveAs(Form("%s/Eta_%s_RatioSpectraToFit.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

        //***************************************************************************************************************
        //******************************* Ratio to fit for individual spectra used range ********************************
        //***************************************************************************************************************
        
        histo2DRatioToFitEta->DrawCopy(); 

        TGraphAsymmErrors* graphCorrectedYieldToFitEtaUsed[MaxNumberOfFiles];
        TGraphAsymmErrors* graphCorrectedYieldToFitEtaSysUsed[MaxNumberOfFiles];
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if( !maskedFullyEta[i]){
                if (graphsCorrectedYieldShrunkEta[i]) graphCorrectedYieldToFitEtaUsed[i] = CalculateGraphErrRatioToFit (graphsCorrectedYieldShrunkEta[i], fitInvYieldEta); 
                if (graphsCorrectedYieldSysShrunkEta[i]) graphCorrectedYieldToFitEtaSysUsed[i] = CalculateGraphErrRatioToFit(graphsCorrectedYieldSysShrunkEta[i], fitInvYieldEta); 
                
                if (graphsCorrectedYieldSysShrunkEta[i]) DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitEtaSysUsed[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
                if (graphsCorrectedYieldShrunkEta[i])DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitEtaUsed[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);

                if (graphsCorrectedYieldSysShrunkEta[i]) graphCorrectedYieldToFitEtaSysUsed[i]->Draw("p,E2,same");
                
                DrawGammaLines(0., maxPtGlobalEta , 1., 1.,0.1, kGray+2);
                DrawGammaLines(0., maxPtGlobalEta , 1.1, 1.1,0.1, kGray, 7);
                DrawGammaLines(0., maxPtGlobalEta , 0.9, 0.9,0.1, kGray, 7);

                if (graphsCorrectedYieldShrunkEta[i])graphCorrectedYieldToFitEtaUsed[i]->Draw("e1,same"); 
            }   
        }
        legendRatioSpecEta->Draw();

        labelEnergyRatio->Draw();
        labelEtaRatio->Draw();
        labelDetProcRatio->Draw();
        
        canvasRatioSpec->Update();
        canvasRatioSpec->SaveAs(Form("%s/Eta_%s_RatioSpectraToFitUsed.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

        
        //***************************************************************************************************************
        //***************************************** Ratio to fit for final spectrum *************************************
        //***************************************************************************************************************
        
        histo2DRatioToFitEta->DrawCopy(); 
        
        TGraphAsymmErrors* graphCorrectedYieldFinalStatToFitEta;
        TGraphAsymmErrors* graphCorrectedYieldFinalSysToFitEta;
        graphCorrectedYieldFinalStatToFitEta    = CalculateGraphErrRatioToFit (graphCorrectedYieldWeightedAverageEtaStat, fitInvYieldEta); 
        graphCorrectedYieldFinalSysToFitEta     = CalculateGraphErrRatioToFit(graphCorrectedYieldWeightedAverageEtaSys, fitInvYieldEta); 

        DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldFinalSysToFitEta, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldFinalStatToFitEta, 24, 2, kBlack , kBlack, 1, kTRUE);
        
        graphCorrectedYieldFinalSysToFitEta->Draw("p,E2,same");
        
        DrawGammaLines(0., maxPtGlobalEta , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0., maxPtGlobalEta , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0., maxPtGlobalEta , 0.9, 0.9,0.1, kGray, 7);

        graphCorrectedYieldFinalStatToFitEta->Draw("p,E,same");

        labelEnergyRatio->Draw();
        labelEtaRatio->Draw();
        labelDetProcRatio->Draw();
        
        canvasRatioSpec->Update();
        canvasRatioSpec->SaveAs(Form("%s/Eta_%s_RatioSpectraToFitFinal.%s",outputDir.Data(),isMC.Data(), suffix.Data()));
        
        delete canvasRatioSpec;

        //***************************************************************************************************************
        //************************************ Eta to Pi0 ratio for different triggers **********************************
        //***************************************************************************************************************
        if (doEtaToPi0){
            cout << "**************************************************************************************************" << endl;
            cout << "************************ Combining different triggers for eta to pi0 *****************************" << endl;
            cout << "**************************************************************************************************" << endl;
            
            // prepare systematics arrays
            Double_t xValueFinalEtaToPi0                    [100];
            Double_t xErrorLowFinalEtaToPi0                 [100];
            Double_t xErrorHighFinalEtaToPi0                [100];
            Double_t yValueFinalEtaToPi0                    [100];
            Double_t yErrorLowFinalEtaToPi0                 [100];
            Double_t yErrorHighFinalEtaToPi0                [100];
            Int_t nPointFinalEtaToPi0                                         = 0;

            Double_t yErrorSysLowFinalEtaToPi0              [100];
            Double_t yErrorSysHighFinalEtaToPi0             [100];

            Double_t ptSysRelEtaToPi0                       [MaxNumberOfFiles][100];
            Double_t yErrorSysLowRelEtaToPi0                [MaxNumberOfFiles][100];
            Double_t yErrorSysHighRelEtaToPi0               [MaxNumberOfFiles][100];
            Bool_t   sysAvailEtaToPi0                       [MaxNumberOfFiles];

            Bool_t sysAvailSingleEtaToPi0                   [MaxNumberOfFiles];
            Int_t numberBinsSysAvailSingleEtaToPi0          [MaxNumberOfFiles];
            
            // create graphs for shrunk individual triggers
            TGraphAsymmErrors* graphsEtaToPi0Shrunk         [MaxNumberOfFiles];
            TGraphAsymmErrors* graphsEtaToPi0SysShrunk      [MaxNumberOfFiles];
            TH1D* histoEtaToPi0Masked                       [MaxNumberOfFiles];

            // create variables for combination functions
            TH1D*               histoStatEtaToPi0    [12]; 
            TGraphAsymmErrors*  graphSystEtaToPi0    [12];
            TH1D*               histoRelStatEtaToPi0 [12]; 
            TGraphAsymmErrors*  graphRelSystEtaToPi0 [12];
            Int_t offSetsEtaToPi0[12]           = { 0, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 0 };
            Int_t offSetsEtaToPi0Sys[12]        = { 0, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 0 };

            if(optionEnergy.CompareTo("8TeV")==0){
              if(mode == 2){
                offSetsEtaToPi0[1] = 0; //INT7
                offSetsEtaToPi0[3] = 0; //EMC7
                offSetsEtaToPi0[4] = 2; //EGA
              }else if(mode == 4){
                offSetsEtaToPi0[1] = 0; //INT7
                offSetsEtaToPi0[3] = 0; //EMC7
                offSetsEtaToPi0[4] = 2; //EGA
              }
            }


            Bool_t hasSysEtaToPi0            = kFALSE;
            for (Int_t j = 0; j<12; j++){
                histoStatEtaToPi0[j]        = NULL;
                graphSystEtaToPi0[j]        = NULL;
                histoRelStatEtaToPi0[j]     = NULL;
                graphRelSystEtaToPi0[j]     = NULL;
            }    
            
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                // read systematics file
                cout << triggerName[i].Data() << endl;
                if (sysFileEtaToPi0[i].CompareTo("bla") != 0){
                    sysAvailEtaToPi0[i]              = kTRUE;
                    ifstream  fileSysErrEtaToPi0;
                    fileSysErrEtaToPi0.open(sysFileEtaToPi0[i].Data(),ios_base::in);
                    cout << sysFileEtaToPi0[i].Data() << endl;
                    Int_t counter = 0;
                    while(!fileSysErrEtaToPi0.eof() && counter < 100){
                        Double_t garbage = 0;
                        fileSysErrEtaToPi0 >>ptSysRelEtaToPi0[i][counter] >> yErrorSysLowRelEtaToPi0[i][counter] >> yErrorSysHighRelEtaToPi0[i][counter]>>    garbage >> garbage;
                        cout << counter << "\t"<< ptSysRelEtaToPi0[i][counter]<< "\t"  << yErrorSysLowRelEtaToPi0[i][counter] << "\t"  <<yErrorSysHighRelEtaToPi0[i][counter] << "\t"  << endl;;
                        counter++;
                    }
                    fileSysErrEtaToPi0.close();
                    hasSysEtaToPi0              = kTRUE;
                 // read in detailed systematics
                    string sysFileEtaToPi0Det = sysFileEtaToPi0[i].Data();
                    if(!replace(sysFileEtaToPi0Det, "Averaged", "AveragedSingle")){cout << "WARNING: could not find detailed systematics file " << sysFileEtaToPi0Det << ", skipping... " << endl; sysAvailSingleEtaToPi0[i] = kFALSE; continue;}
                    ifstream fileSysErrDetailedEtaToPi0;
                    fileSysErrDetailedEtaToPi0.open(sysFileEtaToPi0Det,ios_base::in);
                    if(fileSysErrDetailedEtaToPi0.is_open()) 
                        sysAvailSingleEtaToPi0[i] = kTRUE;
                    else{ 
                        sysAvailSingleEtaToPi0[i] = kFALSE; 
                        cout << "couldn't find single errors for eta/pi0, jumping" << endl;    
                    }
                    
                    if (sysAvailSingleEtaToPi0[i]){
                        cout << sysFileEtaToPi0Det << endl;
                        counter = 0;
                        string line;
                        Int_t counterColumn = 0;
                        while (getline(fileSysErrDetailedEtaToPi0, line) && counter < 100) {
                            istringstream ss(line);
                            TString temp="";
                            counterColumn = 0;
                            while(ss && counterColumn < 100){
                                ss >> temp;
                                if( !(counter==0 && temp.CompareTo("bin")==0) && !temp.IsNull()){
                                ptSysDetail[i][counter].push_back(temp);
                                counterColumn++;
                                }
                            }
                            if(counter == 0){
                                ptSysDetail[i][counter++].push_back("TotalError");
                                counterColumn++;
                            }else counter++;
                        }
                        numberBinsSysAvailSingleEtaToPi0[i] = counter;
                        fileSysErrDetailedEtaToPi0.close();
                    }
                } else {
                    sysAvailEtaToPi0[i]         = kFALSE;
                    sysAvailSingleEtaToPi0[i]   = kTRUE;
                }
                
                // fill graphs to be shrunk later
                cout << "step 1" << endl;
                graphsEtaToPi0Shrunk[i]         = new TGraphAsymmErrors(histoEtaToPi0[i]);
                graphsEtaToPi0Removed0[i]       = new TGraphAsymmErrors(histoEtaToPi0[i]);
                graphsEtaToPi0SysShrunk[i]      = new TGraphAsymmErrors(histoEtaToPi0[i]);
                graphsEtaToPi0SysRemoved0[i]    = new TGraphAsymmErrors(histoEtaToPi0[i]);
                histoEtaToPi0Masked[i]          = (TH1D*)histoEtaToPi0[i]->Clone(Form("EtaToPi0%s_Masked_%s",addNameBinshift.Data(), triggerName[i].Data()));
                
                if(optionEnergy.CompareTo("8TeV")==0 && mode==4 && triggerName[i].Contains("EGA")) ptFromSpecEta[i][0] = 16;
                cout << ptFromSpecEta[i][0] << endl;
                // remove 0 bins at beginning according to ptFromSpecEta[i][0]
                Int_t binsToMask = 1;
                while (histoEtaToPi0Masked[i]->GetBinCenter(binsToMask) < ptFromSpecEta[i][0] ){
                    histoEtaToPi0Masked[i]->SetBinContent(binsToMask,0.);
                    histoEtaToPi0Masked[i]->SetBinError(binsToMask,0.);
                    binsToMask++;
                }
                // check if trigger needs to be masked completely
                if ( maskedFullyEta[i] || maskedFullyPi0[i] ){
                    graphsEtaToPi0Shrunk[i]         = NULL;
                    graphsEtaToPi0Removed0[i]       = NULL;
                    graphsEtaToPi0SysShrunk[i]      = NULL;
                    graphsEtaToPi0SysRemoved0[i]    = NULL;
                    for (Int_t f = 1; f < histoEtaToPi0Masked[i]->GetNbinsX()+1; f++ ){
                        histoEtaToPi0Masked[i]->SetBinContent(f,0.);
                        histoEtaToPi0Masked[i]->SetBinError(f,0.);
                    }
                    nrOfTrigToBeCombEtaToPi0Red--;
                    cout << "trigger " << triggerName[i] << " was masked for Eta/pi0" << endl;
                    continue;
                // remove bins at the end according to set range
                } else if (ptFromSpecEta[i][1] > -1){
                    for (Int_t f = histoEtaToPi0Masked[i]->GetXaxis()->FindBin(ptFromSpecEta[i][1]); f < histoEtaToPi0Masked[i]->GetNbinsX()+1; f++ ){
                        histoEtaToPi0Masked[i]->SetBinContent(f,0.);
                        histoEtaToPi0Masked[i]->SetBinError(f,0.);
                    }
                }
                
                // remove 0 bins at beginning
//                 graphsEtaToPi0Shrunk[i]->Print();
                cout << "step 2" << endl;
                while (graphsEtaToPi0Shrunk[i]->GetY()[0] == 0) graphsEtaToPi0Shrunk[i]->RemovePoint(0);
//                 graphsEtaToPi0Shrunk[i]->Print();
                while (graphsEtaToPi0Removed0[i]->GetY()[0] == 0) graphsEtaToPi0Removed0[i]->RemovePoint(0);
//                 graphsEtaToPi0Removed0[i]->Print();
                cout << "sys shrunk" << endl;
                while (graphsEtaToPi0SysShrunk[i]->GetY()[0] == 0) graphsEtaToPi0SysShrunk[i]->RemovePoint(0);
//                 graphsEtaToPi0SysShrunk[i]->Print();
                cout << "sys shrunk 2" << endl;
                while (graphsEtaToPi0SysRemoved0[i]->GetY()[0] == 0) graphsEtaToPi0SysRemoved0[i]->RemovePoint(0);
//                 graphsEtaToPi0SysRemoved0[i]->Print();
                
                // put systematics on graph
                for (Int_t j = 0; j< graphsEtaToPi0SysRemoved0[i]->GetN(); j++){
                    if (sysAvailEtaToPi0[i]){
                        Int_t counter = 0;
                        while(counter < 100 && TMath::Abs(graphsEtaToPi0SysRemoved0[i]->GetX()[j] - ptSysRelEtaToPi0[i][counter])> 0.001) counter++;
                        if (counter < 100){
                            cout << ptSysRelEtaToPi0[i][counter]<< "\t found it" << endl;
                            Double_t yErrorSysLowDummy  = TMath::Abs(yErrorSysLowRelEtaToPi0[i][counter]/100*graphsEtaToPi0SysRemoved0[i]->GetY()[j]);
                            Double_t yErrorSysHighDummy = yErrorSysHighRelEtaToPi0[i][counter]/100*graphsEtaToPi0SysRemoved0[i]->GetY()[j];
                            graphsEtaToPi0SysRemoved0[i]->SetPointEYlow(j,yErrorSysLowDummy);
                            graphsEtaToPi0SysRemoved0[i]->SetPointEYhigh(j,yErrorSysHighDummy);
                        } else {
                            graphsEtaToPi0SysRemoved0[i]->SetPointEYlow(j,0);
                            graphsEtaToPi0SysRemoved0[i]->SetPointEYhigh(j,0);
                        }
                    } else {
                        graphsEtaToPi0SysRemoved0[i]->SetPointEYlow(j,0);
                        graphsEtaToPi0SysRemoved0[i]->SetPointEYhigh(j,0);
                        averagedEta = kFALSE;
                    }
                }
                
                // remove unused bins at beginning
                cout << "step 3" << endl;
                while (graphsEtaToPi0Shrunk[i]->GetX()[0] < ptFromSpecEta[i][0] || graphsEtaToPi0Shrunk[i]->GetX()[0] < ptFromSpecPi0[i][0]) 
                    graphsEtaToPi0Shrunk[i]->RemovePoint(0);
                while (graphsEtaToPi0SysShrunk[i]->GetX()[0] < ptFromSpecEta[i][0] || graphsEtaToPi0SysShrunk[i]->GetX()[0] < ptFromSpecPi0[i][0]) 
                    graphsEtaToPi0SysShrunk[i]->RemovePoint(0);
                // remove unused bins at end
                while ( graphsEtaToPi0Shrunk[i]->GetX()[graphsEtaToPi0Shrunk[i]->GetN()-1] > ptFromSpecEta[i][1] ||
                        graphsEtaToPi0Shrunk[i]->GetX()[graphsEtaToPi0Shrunk[i]->GetN()-1] > ptFromSpecPi0[i][1] ) 
                    graphsEtaToPi0Shrunk[i]->RemovePoint(graphsEtaToPi0Shrunk[i]->GetN()-1);
                graphsEtaToPi0Shrunk[i]->Print();
                while ( graphsEtaToPi0SysShrunk[i]->GetX()[graphsEtaToPi0SysShrunk[i]->GetN()-1] > ptFromSpecEta[i][1] ||
                        graphsEtaToPi0SysShrunk[i]->GetX()[graphsEtaToPi0SysShrunk[i]->GetN()-1] > ptFromSpecPi0[i][1] ) 
                    graphsEtaToPi0SysShrunk[i]->RemovePoint(graphsEtaToPi0SysShrunk[i]->GetN()-1);

                // put systematics on graph
                for (Int_t j = 0; j< graphsEtaToPi0Shrunk[i]->GetN(); j++){
                    xValueFinalEtaToPi0[nPointFinalEtaToPi0]        = graphsEtaToPi0Shrunk[i]->GetX()[j];
                    xErrorHighFinalEtaToPi0[nPointFinalEtaToPi0]    = graphsEtaToPi0Shrunk[i]->GetEXhigh()[j];
                    xErrorLowFinalEtaToPi0[nPointFinalEtaToPi0]     = graphsEtaToPi0Shrunk[i]->GetEXlow()[j];
                    yValueFinalEtaToPi0[nPointFinalEtaToPi0]        = graphsEtaToPi0Shrunk[i]->GetY()[j];
                    yErrorHighFinalEtaToPi0[nPointFinalEtaToPi0]    = graphsEtaToPi0Shrunk[i]->GetEYhigh()[j];
                    yErrorLowFinalEtaToPi0[nPointFinalEtaToPi0]     = graphsEtaToPi0Shrunk[i]->GetEYlow()[j];
                    if (sysAvailEtaToPi0[i]){
                        Int_t counter = 0;
                        while(counter < 100 && TMath::Abs(xValueFinalEtaToPi0[nPointFinalEtaToPi0] - ptSysRelEtaToPi0[i][counter])> 0.001) counter++;
                        if (counter < 100){
                            cout << ptSysRelEtaToPi0[i][counter]<< "\t found it" << endl;
                            yErrorSysLowFinalEtaToPi0[nPointFinalEtaToPi0] = TMath::Abs(yErrorSysLowRelEtaToPi0[i][counter]/100*graphsEtaToPi0Shrunk[i]->GetY()[j]);
                            yErrorSysHighFinalEtaToPi0[nPointFinalEtaToPi0] = yErrorSysHighRelEtaToPi0[i][counter]/100*graphsEtaToPi0Shrunk[i]->GetY()[j];
                            
                        } else {
                            yErrorSysLowFinalEtaToPi0[nPointFinalEtaToPi0]  = 0;
                            yErrorSysHighFinalEtaToPi0[nPointFinalEtaToPi0] = 0;
                            
                        }
                    } else {
                        yErrorSysLowFinalEtaToPi0[nPointFinalEtaToPi0]  = 0;
                        yErrorSysHighFinalEtaToPi0[nPointFinalEtaToPi0] = 0;
                    }
                    graphsEtaToPi0SysShrunk[i]->SetPointEYlow(j,yErrorSysLowFinalEtaToPi0[nPointFinalEtaToPi0]);
                    graphsEtaToPi0SysShrunk[i]->SetPointEYhigh(j,yErrorSysHighFinalEtaToPi0[nPointFinalEtaToPi0]);
                    nPointFinalEtaToPi0++;
                }
                    
                // fill inputs for combinations function in correct order
                if ( (triggerName[i].Contains("MB") || triggerName[i].Contains("INT1")) && graphsEtaToPi0Shrunk[i]){
                    cout << "filling MB trigger" << endl;
                    histoStatEtaToPi0[0]     = histoEtaToPi0Masked[i];
                    graphSystEtaToPi0[0]     = graphsEtaToPi0SysShrunk[i];
                    offSetsEtaToPi0Sys[0]    = histoStatEtaToPi0[0]->GetXaxis()->FindBin(graphSystEtaToPi0[0]->GetX()[0])-1;
                } else if (triggerName[i].Contains("INT7") && graphsEtaToPi0Shrunk[i]){   
                    cout << "filling INT7 trigger" << endl;
                    histoStatEtaToPi0[1]     = histoEtaToPi0Masked[i];
                    graphSystEtaToPi0[1]     = graphsEtaToPi0SysShrunk[i];
                    offSetsEtaToPi0Sys[1]    = histoStatEtaToPi0[1]->GetXaxis()->FindBin(graphSystEtaToPi0[1]->GetX()[0])-1;
                } else if (triggerName[i].Contains("EMC1") && graphsEtaToPi0Shrunk[i]){   
                    cout << "filling EMC1 trigger" << endl;
                    histoStatEtaToPi0[2]     = histoEtaToPi0Masked[i];
                    graphSystEtaToPi0[2]     = graphsEtaToPi0SysShrunk[i];
                    offSetsEtaToPi0Sys[2]    = histoStatEtaToPi0[2]->GetXaxis()->FindBin(graphSystEtaToPi0[2]->GetX()[0])-1;
                } else if (triggerName[i].Contains("EMC7") && graphsEtaToPi0Shrunk[i]){   
                    cout << "filling EMC7 trigger" << endl;
                    histoStatEtaToPi0[3]     = histoEtaToPi0Masked[i];
                    graphSystEtaToPi0[3]     = graphsEtaToPi0SysShrunk[i];
                    offSetsEtaToPi0Sys[3]    = histoStatEtaToPi0[3]->GetXaxis()->FindBin(graphSystEtaToPi0[3]->GetX()[0])-1;
                } else if ((triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA")) && graphsEtaToPi0Shrunk[i]){
                    cout << Form("filling %s trigger",strEG2_A.Data()) << endl;
                    histoStatEtaToPi0[4]     = histoEtaToPi0Masked[i];
                    graphSystEtaToPi0[4]     = graphsEtaToPi0SysShrunk[i];
                    offSetsEtaToPi0Sys[4]    = histoStatEtaToPi0[4]->GetXaxis()->FindBin(graphSystEtaToPi0[4]->GetX()[0])-1;

                    if(optionEnergy.CompareTo("8TeV")==0 && (mode == 4 || mode == 2)) offSetsEtaToPi0Sys[4]+=2;
                } else if (triggerName[i].Contains("EG1") && graphsEtaToPi0Shrunk[i]){   
                    cout << "filling EG1 trigger" << endl;
                    histoStatEtaToPi0[5]     = histoEtaToPi0Masked[i];
                    graphSystEtaToPi0[5]     = graphsEtaToPi0SysShrunk[i];
                    offSetsEtaToPi0Sys[5]    = histoStatEtaToPi0[5]->GetXaxis()->FindBin(graphSystEtaToPi0[5]->GetX()[0])-1;
                }
            }
            
            TString nameWeightsLogFileEtaToPi0                  = Form("%s/weightsEtaToPi0_%s.dat",outputDir.Data(),isMC.Data());
            TGraphAsymmErrors* graphEtaToPi0WeightedAverageTot  = NULL;
            // Calculate averaged eta/pi0 graphs according to statistical and systematic errors taking correctly into account the cross correlations
            if (averagedEta){
                if(optionEnergy.CompareTo("8TeV")==0 && mode==4){
                  maxNAllowedEta -= 3;
                  maxPtGlobalEta = 20;
                }
                //if(optionEnergy.CompareTo("8TeV")==0 && mode==2) maxNAllowedEta -= 2;
                // calculate averaged eta/pi0 graphs
                graphEtaToPi0WeightedAverageTot         = CombinePtPointsSpectraTriggerCorrMat( histoStatEtaToPi0, graphSystEtaToPi0,  
                                                                                                binningEta,  maxNAllowedEta,
                                                                                                offSetsEtaToPi0 ,offSetsEtaToPi0Sys,
                                                                                                graphEtaToPi0WeightedAverageStat, graphEtaToPi0WeightedAverageSys,
                                                                                                nameWeightsLogFileEtaToPi0.Data(),
                                                                                                mode, optionEnergy, "EtaToPi0", kFALSE,
                                                                                                fileInputCorrFactors
                                                                                              );
                // Reading weights from output file for plotting
                ifstream fileWeightsEtaToPi0;
                fileWeightsEtaToPi0.open(nameWeightsLogFileEtaToPi0,ios_base::in);
                cout << "reading" << nameWeightsLogFileEtaToPi0 << endl;
                Double_t xValuesReadEtaToPi0[100];
                Double_t weightsReadEtaToPi0[12][100];
                Int_t availableMeasEtaToPi0[12]     = { -1, -1, -1, -1, -1, -1,
                                                        -1, -1, -1, -1, -1, -1};
                Int_t nMeasSetEtaToPi0              = numberOfTrigg;
                Int_t nPtBinsReadEtaToPi0           = 0;
            
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

                while(!fileWeightsEtaToPi0.eof() && nPtBinsReadEtaToPi0 < 100){
                    TString garbage = "";
                    if (nPtBinsReadEtaToPi0 == 0){
                        fileWeightsEtaToPi0 >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
                        for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                            fileWeightsEtaToPi0 >> availableMeasEtaToPi0[i] ;
                        }   
                        cout << "read following measurements: "; 
                        for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                            cout << availableMeasEtaToPi0[i] << "\t" ;
                        }   
                        cout << endl;
                    } else {
                        fileWeightsEtaToPi0 >> xValuesReadEtaToPi0[nPtBinsReadEtaToPi0-1];
                        for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                            fileWeightsEtaToPi0 >> weightsReadEtaToPi0[availableMeasEtaToPi0[i]][nPtBinsReadEtaToPi0-1] ;
                        }   
                        cout << "read: "<<  nPtBinsReadEtaToPi0 << "\t"<< xValuesReadEtaToPi0[nPtBinsReadEtaToPi0-1] << "\t" ;
                        for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                            cout << weightsReadEtaToPi0[availableMeasEtaToPi0[i]][nPtBinsReadEtaToPi0-1] << "\t";
                        }
                        cout << endl;
                    }
                    nPtBinsReadEtaToPi0++;
                }
                nPtBinsReadEtaToPi0 = nPtBinsReadEtaToPi0-2 ;
                fileWeightsEtaToPi0.close();

                TGraph* graphWeightsEtaToPi0[nMeasSetEtaToPi0];
                for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                    graphWeightsEtaToPi0[availableMeasEtaToPi0[i]] = new TGraph(nPtBinsReadEtaToPi0,xValuesReadEtaToPi0,weightsReadEtaToPi0[availableMeasEtaToPi0[i]]);
                    Int_t bin = 0;
                    for (Int_t n = 0; n< nPtBinsReadEtaToPi0; n++){
                        if (graphWeightsEtaToPi0[availableMeasEtaToPi0[i]]->GetY()[bin] == 0) graphWeightsEtaToPi0[availableMeasEtaToPi0[i]]->RemovePoint(bin);
                        else bin++;
                    }   
                }   

                //  **********************************************************************************************************************
                //  **************************************** Combine+write detailed Systematics ******************************************
                //  **********************************************************************************************************************

                const char *SysErrDatnameMeanSingleErr = Form("%s/SystematicErrorAveragedSingle%s_Pi0EtaBinning_%s.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
                fstream SysErrDatAverSingle;
                SysErrDatAverSingle.precision(4);
                cout << SysErrDatnameMeanSingleErr << endl;
                if(sysAvailSingleEtaToPi0[0]){
                  SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
                  for(Int_t iColumn = 0; iColumn < (Int_t)ptSysDetail[0][0].size(); iColumn++) SysErrDatAverSingle << ptSysDetail[0][0].at(iColumn) << "\t";
                  SysErrDatAverSingle << endl;
                  for (Int_t i = 1; i < nrOfTrigToBeComb; i++){
                    if(!sysAvailSingleEtaToPi0[i]) continue;
                    for(Int_t iCol = 0; iCol < (Int_t)ptSysDetail[i][0].size(); iCol++){
                      if( ((TString)ptSysDetail[i][0].at(iCol)).CompareTo(((TString)ptSysDetail[0][0].at(iCol))) ){
                        cout << "ERROR: Systematic error type at pos " << iCol << " does not agree for " << availableMeasEtaToPi0[i] << " & " << availableMeasEtaToPi0[0] << ", returning!" << endl;
                        return;
                      }
                    }
                  }

                  for(Int_t i=0; i<nPtBinsReadEtaToPi0; i++){
                    SysErrDatAverSingle << xValuesReadEtaToPi0[i] << "\t";
                    Int_t nColumns = (Int_t)ptSysDetail[0][0].size();
                    Double_t *errors = new Double_t[nColumns-1];
                    for(Int_t iErr=0; iErr<nColumns-1; iErr++) errors[iErr] = 0;
                    for(Int_t j=0; j<nrOfTrigToBeComb; j++){
                      if(!sysAvailSingleEtaToPi0[j]) continue;
                      Int_t pos = GetBinPosInVec(ptSysDetail[j],numberBinsSysAvailSingleEtaToPi0[j],xValuesReadEtaToPi0[i]);
                      if(pos>-1){
                        for(Int_t iErr=1; iErr<nColumns; iErr++) errors[iErr-1] += weightsReadEtaToPi0[availableMeasEtaToPi0[j]][i]*((TString)ptSysDetail[j][pos].at(iErr)).Atof();
                      }
                    }
                    for(Int_t iErr=0; iErr<nColumns-1; iErr++) SysErrDatAverSingle << errors[iErr] << "\t";
                    SysErrDatAverSingle << endl;
                    delete[] errors;
                  }
                }
                SysErrDatAverSingle.close();

                for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++){
                  for(Int_t iB=0; iB<50; iB++) ptSysDetail[iR][iB].clear();
                }

                //  **********************************************************************************************************************
                //  ******************************************* Plotting weights for eta/pi0 *********************************************
                //  **********************************************************************************************************************
                Int_t textSizeLabelsPixel = 900*0.04;

                TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);// gives the page size
                DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
            
                TH2F * histo2DWeightsEtaToPi0;
                histo2DWeightsEtaToPi0 = new TH2F("histo2DWeightsEtaToPi0","histo2DWeightsEtaToPi0",11000,0.,maxPtGlobalEta,1000,-0.5,1.3);
                SetStyleHistoTH2ForGraphs(histo2DWeightsEtaToPi0, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
                histo2DWeightsEtaToPi0->Draw("copy");
                
                    TLegend* legendWeightsEtaToPi0 = GetAndSetLegend2(0.12, 0.14, 0.55, 0.14+(0.035*nMeasSetEtaToPi0/2*1.35), 32);
                    legendWeightsEtaToPi0->SetNColumns(2);
                    for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                        DrawGammaSetMarkerTGraph(graphWeightsEtaToPi0[availableMeasEtaToPi0[i]], markerTriggWeighted[availableMeasEtaToPi0[i]], sizeTrigg[availableMeasEtaToPi0[i]], 
                                                colorTriggWeighted[availableMeasEtaToPi0[i]], colorTriggWeighted[availableMeasEtaToPi0[i]]);
                        graphWeightsEtaToPi0[availableMeasEtaToPi0[i]]->Draw("p,same,e1");
                        legendWeightsEtaToPi0->AddEntry(graphWeightsEtaToPi0[availableMeasEtaToPi0[i]],nameTriggerWeighted[availableMeasEtaToPi0[i]],"p");
                    }   
                    legendWeightsEtaToPi0->Draw();

                    TLatex *labelWeightsEnergy = new TLatex(0.7,0.24,collisionSystem.Data());
                    SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
                    labelWeightsEnergy->SetTextFont(43);
                    labelWeightsEnergy->Draw();
                    TLatex *labelWeightsEtaToPi0 = new TLatex(0.7,0.20,"#eta/#pi^{0}");
                    SetStyleTLatex( labelWeightsEtaToPi0, 0.85*textSizeLabelsPixel,4);
                    labelWeightsEtaToPi0->SetTextFont(43);
                    labelWeightsEtaToPi0->Draw();
                    TLatex *labelDetProcWeights    = new TLatex(0.7, 0.16,detectionProcess.Data());
                    SetStyleTLatex( labelDetProcWeights, 0.85*textSizeLabelsPixel,4);
                    labelDetProcWeights->SetTextFont(43);
                    labelDetProcWeights->Draw();

            //      DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
                    DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
                    DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
                    DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
                    DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
                    
                canvasWeights->SaveAs(Form("%s/%s_WeightsEtaToPi0Triggers.%s",outputDir.Data(), isMC.Data(), suffix.Data()));
                delete canvasWeights;

                // Calculating relative error for eta/pi0
                for (Int_t i = 0; i < 12; i++){
                    if (histoStatEtaToPi0[i]) 
                        histoRelStatEtaToPi0[i]      = CalculateRelErrUpTH1D( histoStatEtaToPi0[i], Form("relativeStatErrorEtaToPi0_%s", nameTriggerWeighted[i].Data()));
                    if (graphSystEtaToPi0[i]) 
                        graphRelSystEtaToPi0[i]      = CalculateRelErrUpAsymmGraph( graphSystEtaToPi0[i], Form("relativeSysErrorEtaToPi0_%s", nameTriggerWeighted[i].Data()));
                }
                TGraphAsymmErrors* graphRelErrorEtaToPi0Tot        = CalculateRelErrUpAsymmGraph( graphEtaToPi0WeightedAverageTot, "relativeTotalErrorEtaToPi0");
                TGraphAsymmErrors* graphRelErrorEtaToPi0Stat       = CalculateRelErrUpAsymmGraph( graphEtaToPi0WeightedAverageStat, "relativeStatErrorEtaToPi0");
                TGraphAsymmErrors* graphRelErrorEtaToPi0Sys        = CalculateRelErrUpAsymmGraph( graphEtaToPi0WeightedAverageSys, "relativeSysErrorEtaToPi0");

                // plot sys relative errors for individual triggers   
                TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
            
                TH2F * histo2DRelSysErr;
                histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.,maxPtGlobalEta,1000,0,60.5);
                SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
                histo2DRelSysErr->Draw("copy");
                    TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaToPi0/2), 0.95, 0.92, 32);
                    legendRelSysErr->SetNColumns(2);
                    for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                        DrawGammaSetMarkerTGraph(graphRelSystEtaToPi0[availableMeasEtaToPi0[i]], markerTriggWeighted[availableMeasEtaToPi0[i]], sizeTrigg[availableMeasEtaToPi0[i]], 
                                                colorTriggWeighted[availableMeasEtaToPi0[i]], colorTriggWeighted[availableMeasEtaToPi0[i]]);
                        graphRelSystEtaToPi0[availableMeasEtaToPi0[i]]->Draw("p,same,z");
                        legendRelSysErr->AddEntry(graphRelSystEtaToPi0[availableMeasEtaToPi0[i]],nameTriggerWeighted[availableMeasEtaToPi0[i]],"p");
                    }    
                    legendRelSysErr->Draw();

                    TLatex *labelRelErrEnergy    = new TLatex(0.15,0.89,collisionSystem.Data());
                    SetStyleTLatex( labelRelErrEnergy, 0.85*textSizeLabelsPixel,4);
                    labelRelErrEnergy->SetTextFont(43);
                    labelRelErrEnergy->Draw();
                    TLatex *labelRelErrEtaToPi0       = new TLatex(0.15,0.85,"#eta/#pi^{0}");
                    SetStyleTLatex( labelRelErrEtaToPi0, 0.85*textSizeLabelsPixel,4);
                    labelRelErrEtaToPi0->SetTextFont(43);
                    labelRelErrEtaToPi0->Draw();
                    TLatex *labelDetProcRelErr    = new TLatex(0.15, 0.81,detectionProcess.Data());
                    SetStyleTLatex( labelDetProcRelErr, 0.85*textSizeLabelsPixel,4);
                    labelDetProcRelErr->SetTextFont(43);
                    labelDetProcRelErr->Draw();
                    
                canvasRelSysErr->SaveAs(Form("%s/EtaToPi0_RelSysErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));
                
                // plot stat relative errors for individual triggers    
                TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
            
                TH2F * histo2DRelStatErr;
                histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.,maxPtGlobalEta,1000,0,60.5);
                SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
                histo2DRelStatErr->Draw("copy");
                    TLegend* legendRelStatErr       = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaToPi0/2), 0.95, 0.92, 32);
                    legendRelStatErr->SetNColumns(2);
                    for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                        if (availableMeasEtaToPi0[i] == 0 && histoRelStatEtaToPi0[availableMeasEtaToPi0[i]] && mode == 2){
                            TGraphAsymmErrors* dummyGraph = new TGraphAsymmErrors(histoRelStatEtaToPi0[availableMeasEtaToPi0[i]]);
                            dummyGraph->Print();
                            DrawGammaSetMarkerTGraph(dummyGraph, markerTriggWeighted[availableMeasEtaToPi0[i]], sizeTrigg[availableMeasEtaToPi0[i]], 
                                                colorTriggWeighted[availableMeasEtaToPi0[i]], colorTriggWeighted[availableMeasEtaToPi0[i]]);
                            dummyGraph->Draw("pX,same");
                            legendRelSysErr->AddEntry(dummyGraph,nameTriggerWeighted[availableMeasEtaToPi0[i]],"p");
                            
        //                     for (Int_t j = 1; j < histoRelStatEtaToPi0[availableMeasEtaToPi0[i]]->GetNbinsX()+1; j++){
        //                        cout << j << ": " << histoRelStatEtaToPi0[availableMeasEtaToPi0[i]]->GetBinContent(j) << endl;
        //                     }    
                        } else {
                            DrawGammaSetMarker(histoRelStatEtaToPi0[availableMeasEtaToPi0[i]],markerTriggWeighted[availableMeasEtaToPi0[i]], sizeTrigg[availableMeasEtaToPi0[i]], 
                                                colorTriggWeighted[availableMeasEtaToPi0[i]], colorTriggWeighted[availableMeasEtaToPi0[i]]);
                            histoRelStatEtaToPi0[availableMeasEtaToPi0[i]]->Draw("p,same,z");
                            legendRelStatErr->AddEntry(histoRelStatEtaToPi0[availableMeasEtaToPi0[i]],nameTriggerWeighted[availableMeasEtaToPi0[i]],"p");
                        }    
                    }    
                    legendRelStatErr->Draw();

                    labelRelErrEnergy->Draw();
                    labelRelErrEtaToPi0->Draw();
                    labelDetProcRelErr->Draw();
                    
                canvasRelStatErr->SaveAs(Form("%s/EtaToPi0_RelStatErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));
                
                // plot full error for final result decomposed
                TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
            
                TH2F * histo2DRelTotErrEtaToPi0;
                histo2DRelTotErrEtaToPi0                 = new TH2F("histo2DRelTotErrEtaToPi0","histo2DRelTotErrEtaToPi0",11000,0.,maxPtGlobalEta,1000,0,60.5);
                SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
                histo2DRelTotErrEtaToPi0->Draw("copy");

                    DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaToPi0Tot, 20, 1.5, kBlack , kBlack);
                    graphRelErrorEtaToPi0Tot->Draw("p,same,z");
                    DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaToPi0Stat, 24, 1.5, kGray+2 , kGray+2);
                    graphRelErrorEtaToPi0Stat->Draw("l,x0,same,e1");
                    DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaToPi0Sys, 24, 1.5, kGray+1 , kGray+1);
                    graphRelErrorEtaToPi0Sys->SetLineStyle(7);
                    graphRelErrorEtaToPi0Sys->Draw("l,x0,same,e1");

                    TLegend* legendRelTotErr2       = GetAndSetLegend2(0.72, 0.92-(0.035*3), 0.9, 0.92, 32);
                    legendRelTotErr2->AddEntry(graphRelErrorEtaToPi0Tot,"tot","p");
                    legendRelTotErr2->AddEntry(graphRelErrorEtaToPi0Stat,"stat","l");
                    legendRelTotErr2->AddEntry(graphRelErrorEtaToPi0Sys,"sys","l");
                    legendRelTotErr2->Draw();

                    labelRelErrEnergy->Draw();
                    labelRelErrEtaToPi0->Draw();
                    labelDetProcRelErr->Draw();
                    
                canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_RelErrorsFulldecomp.%s",outputDir.Data(),suffix.Data()));
                
                if(optionEnergy.CompareTo("8TeV")==0 && mode==4) maxNAllowedEta += 3;
                //if(optionEnergy.CompareTo("8TeV")==0 && mode==2) maxNAllowedEta += 2;
            // if averaging wasn't enabled pick values according to predefined ranges ("cherry picking points")            
            } else {
                graphEtaToPi0WeightedAverageStat        = new TGraphAsymmErrors(nPointFinalEtaToPi0, xValueFinalEtaToPi0, yValueFinalEtaToPi0, 
                                                                                         xErrorLowFinalEtaToPi0, xErrorHighFinalEtaToPi0,yErrorLowFinalEtaToPi0, yErrorHighFinalEtaToPi0);
                graphEtaToPi0WeightedAverageSys         = new TGraphAsymmErrors(nPointFinalEtaToPi0, xValueFinalEtaToPi0, yValueFinalEtaToPi0, 
                                                                                        xErrorLowFinalEtaToPi0, xErrorHighFinalEtaToPi0,yErrorSysLowFinalEtaToPi0, yErrorSysHighFinalEtaToPi0);
            }    
            // printing final eta/pi0 ratios
            cout << "stat eta/pi0" << endl; 
            graphEtaToPi0WeightedAverageStat->Print();
            cout << "sys eta/pi0" << endl;
            if (graphEtaToPi0WeightedAverageSys) graphEtaToPi0WeightedAverageSys->Print();

            if (graphEtaToPi0WeightedAverageStat){
                while (graphEtaToPi0WeightedAverageStat->GetX()[0]< minPtGlobalEta){
                    graphEtaToPi0WeightedAverageStat->RemovePoint(0);
                }    
            }
            
            if (graphEtaToPi0WeightedAverageSys){
                while (graphEtaToPi0WeightedAverageSys->GetX()[0]< minPtGlobalEta){
                    graphEtaToPi0WeightedAverageSys->RemovePoint(0);
                }    
            }

            
            //***************************************************************************************************************
            //***************************Plotting individual eta/pi0 ratios full range **************************************
            //***************************************************************************************************************
            
            Size_t textSizeEtaToPi0 = 0.04;            
            TCanvas* canvasEtatoPi0combo = new TCanvas("canvasEtatoPi0combo","",200,10,1200,1100);
            DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.09, 0.017, 0.02, 0.1);

            TH2F * histo2DEtatoPi0combo = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",11000,0,maxPtGlobalEta,1000,0.01,1.2);
            SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}",0.035,0.04, 0.035,0.04, 1.2,1.);
            histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.01,1.05);

            histo2DEtatoPi0combo->DrawCopy();
            DrawGammaLines(0., maxPtGlobalEta , 0.45, 0.45, 0.3, kGray+2);
            DrawGammaLines(0., maxPtGlobalEta , 0.4, 0.4, 0.3, kGray, 7);
            DrawGammaLines(0., maxPtGlobalEta , 0.5, 0.5, 0.3, kGray, 7);


            TLegend* legendEtaToPi0 = GetAndSetLegend2(0.32, 0.14, 0.6, 0.14+(0.035*(nrOfTrigToBeComb/2+1)*0.9), 32);
            legendEtaToPi0->SetNColumns(2);
            
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if (foundPi0EtaBinFile[i] && !maskedFullyEta[i] && !maskedFullyPi0[i] ){
                    DrawGammaSetMarker(histoEtaToPi0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    if (graphsEtaToPi0SysRemoved0[i]){
                        DrawGammaSetMarkerTGraphAsym(graphsEtaToPi0SysRemoved0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
                        graphsEtaToPi0SysRemoved0[i]->Draw("p,E2,same");
                    }    
                    histoEtaToPi0[i]->Draw("e1,same");
                    legendEtaToPi0->AddEntry(histoEtaToPi0[i],triggerNameLabel[i],"p");
                }
            }
            legendEtaToPi0->Draw();

            TLatex *labelEnergyEtaToPi0 = new TLatex(0.65, 0.15+(2*textSizeEtaToPi0*0.85),collisionSystem.Data());
            SetStyleTLatex( labelEnergyEtaToPi0, 0.85*textSizeEtaToPi0,4);
            labelEnergyEtaToPi0->Draw();
            
            TLatex *labelPi0EtaToPi0 = new TLatex(0.65, 0.15+textSizeEtaToPi0*0.85,"#eta/#pi^{0}");
            SetStyleTLatex( labelPi0EtaToPi0, 0.85*textSizeEtaToPi0,4);
            labelPi0EtaToPi0->Draw();
            
            TLatex *labelDetProcEtaToPi0 = new TLatex(0.65, 0.15,detectionProcess.Data());
            SetStyleTLatex( labelDetProcEtaToPi0, 0.85*textSizeEtaToPi0,4);
            labelDetProcEtaToPi0->Draw();

            histo2DEtatoPi0combo->Draw("axis,same");

            canvasEtatoPi0combo->Update();
            canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0%s_%s_diffTriggers.%s",outputDir.Data(), addNameBinshift.Data(), isMC.Data(), suffix.Data()));

            //***************************************************************************************************************
            //***************************Plotting individual eta/pi0 ratios used range with full ****************************
            //***************************************************************************************************************

            histo2DEtatoPi0combo->DrawCopy();
            DrawGammaLines(0., maxPtGlobalEta , 0.45, 0.45, 0.3, kGray+2);
            DrawGammaLines(0., maxPtGlobalEta , 0.4, 0.4, 0.3, kGray, 7);
            DrawGammaLines(0., maxPtGlobalEta , 0.5, 0.5, 0.3, kGray, 7);

            
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if (foundPi0EtaBinFile[i] && !maskedFullyEta[i] && !maskedFullyPi0[i]){
                    if (graphsEtaToPi0Shrunk[i]) DrawGammaSetMarkerTGraphAsym(graphsEtaToPi0Shrunk[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    if (graphsEtaToPi0SysShrunk[i]){
                        DrawGammaSetMarkerTGraphAsym(graphsEtaToPi0SysShrunk[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
                        graphsEtaToPi0SysShrunk[i]->Draw("p,E2,same");
                    }    
                    if (graphsEtaToPi0Shrunk[i]) graphsEtaToPi0Shrunk[i]->Draw("e1,same");
                }
            }
            
            if (graphEtaToPi0WeightedAverageSys){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0WeightedAverageSys,  24, 2, kRed , kRed, 1, kTRUE);
                graphEtaToPi0WeightedAverageSys->Draw("p,E2,same");        
            }    
            if (graphEtaToPi0WeightedAverageStat){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0WeightedAverageStat,  24, 2, kRed , kRed);
                graphEtaToPi0WeightedAverageStat->Draw("p,E,same");
                legendEtaToPi0->AddEntry(graphEtaToPi0WeightedAverageStat, "final","p");
            }    
            legendEtaToPi0->Draw();

            labelEnergyEtaToPi0->Draw();
            labelPi0EtaToPi0->Draw();
            labelDetProcEtaToPi0->Draw();
            
            histo2DEtatoPi0combo->Draw("axis,same");

            canvasEtatoPi0combo->Update();
            canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0%s_%s_diffTriggers_used.%s",outputDir.Data(), addNameBinshift.Data(), isMC.Data(), suffix.Data()));

            //***************************************************************************************************************
            //**************************************Plotting final eta/pi0 ratio ********************************************
            //***************************************************************************************************************
            
            histo2DEtatoPi0combo->DrawCopy();
            DrawGammaLines(0., maxPtGlobalEta , 0.45, 0.45, 0.3, kGray+2);
            DrawGammaLines(0., maxPtGlobalEta , 0.4, 0.4, 0.3, kGray, 7);
            DrawGammaLines(0., maxPtGlobalEta , 0.5, 0.5, 0.3, kGray, 7);

            if (graphEtaToPi0WeightedAverageSys){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0WeightedAverageSys,  24, 2, kGray+1 , kGray+1, 1, kTRUE);
                graphEtaToPi0WeightedAverageSys->Draw("p,E2,same");        
            }    
            if (graphEtaToPi0WeightedAverageStat){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0WeightedAverageStat,  24, 2, kBlack , kBlack);
                graphEtaToPi0WeightedAverageStat->Draw("p,E,same");
            }    

            labelEnergyEtaToPi0->Draw();
            labelPi0EtaToPi0->Draw();
            labelDetProcEtaToPi0->Draw();
            
            histo2DEtatoPi0combo->Draw("axis,same");

            canvasEtatoPi0combo->Update();
            canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0%s_%s_Final.%s",outputDir.Data(), addNameBinshift.Data(), isMC.Data(), suffix.Data()));   
        }
    }
    
    //*********************************************************************************************************
    //********************** ComparisonFile Output ************************************************************
    //*********************************************************************************************************
    TString system              = "PCM";
    if (mode == 2) system       = "PCM-EMCAL";
    if (mode == 3) system       = "PCM-PHOS";
    if (mode == 4) system       = "EMCAL-EMCAL";
    if (mode == 5) system       = "PHOS-PHOS";
    if (mode == 10) system      = "EMC-merged";
    if (mode == 11) system      = "PHOS-merged";
    
    cout << "starting to write out data" << endl;
    TGraphAsymmErrors* graphInvXSectionWeightedAveragePi0Stat   = 0;
    TH1D* histoInvXSectionWeightedAveragePi0Stat                = 0;
    TGraphAsymmErrors* graphInvXSectionWeightedAveragePi0Sys    = 0;
    TGraphAsymmErrors* graphInvXSectionWeightedAverageEtaStat   = 0;
    TH1D* histoInvXSectionWeightedAverageEtaStat                = 0;
    TGraphAsymmErrors* graphInvXSectionWeightedAverageEtaSys    = 0;
    TH1D* histoEtaToPi0WeightedAverageStat                      = 0;
    graphCorrectedYieldWeightedAveragePi0Stat->Print();
    if (graphCorrectedYieldWeightedAveragePi0Stat){
        graphInvXSectionWeightedAveragePi0Stat                  = ScaleGraph(graphCorrectedYieldWeightedAveragePi0Stat,xSection*recalcBarn);
        histoInvXSectionWeightedAveragePi0Stat                  = new TH1D("histoInvXSectionWeightedAveragePi0Stat", "", maxNAllowedPi0, binningPi0);
        Int_t firstBinPi0 = 1;
        while (histoInvXSectionWeightedAveragePi0Stat->GetBinCenter(firstBinPi0) < graphInvXSectionWeightedAveragePi0Stat->GetX()[0]){
            histoInvXSectionWeightedAveragePi0Stat->SetBinContent(firstBinPi0, 0);
            histoInvXSectionWeightedAveragePi0Stat->SetBinError(firstBinPi0, 0);
            firstBinPi0++;
        }    
        for (Int_t i = 0; i < graphInvXSectionWeightedAveragePi0Stat->GetN(); i++){
            histoInvXSectionWeightedAveragePi0Stat->SetBinContent(i+firstBinPi0, graphInvXSectionWeightedAveragePi0Stat->GetY()[i]);
            histoInvXSectionWeightedAveragePi0Stat->SetBinError(i+firstBinPi0, graphInvXSectionWeightedAveragePi0Stat->GetEYlow()[i]);
        }
    }    
    graphInvXSectionWeightedAveragePi0Stat->Print();
    if (graphCorrectedYieldWeightedAveragePi0Sys) 
        graphInvXSectionWeightedAveragePi0Sys                   = ScaleGraph(graphCorrectedYieldWeightedAveragePi0Sys,xSection*recalcBarn);
    graphInvXSectionWeightedAveragePi0Sys->Print();
    if (enableEta && mode != 10){
        if (graphCorrectedYieldWeightedAverageEtaStat) {
            graphInvXSectionWeightedAverageEtaStat              = ScaleGraph(graphCorrectedYieldWeightedAverageEtaStat,xSection*recalcBarn);
            histoInvXSectionWeightedAverageEtaStat              = new TH1D("histoInvXSectionWeightedAverageEtaStat", "", maxNAllowedEta, binningEta);
            Int_t firstBinEta = 1;
            while (histoInvXSectionWeightedAverageEtaStat->GetBinCenter(firstBinEta) < graphInvXSectionWeightedAverageEtaStat->GetX()[0]){
                histoInvXSectionWeightedAverageEtaStat->SetBinContent(firstBinEta, 0);
                histoInvXSectionWeightedAverageEtaStat->SetBinError(firstBinEta, 0);

                firstBinEta++;
            }
            for (Int_t i = 0; i < graphInvXSectionWeightedAverageEtaStat->GetN(); i++){
                histoInvXSectionWeightedAverageEtaStat->SetBinContent(i+firstBinEta, graphInvXSectionWeightedAverageEtaStat->GetY()[i]);
                histoInvXSectionWeightedAverageEtaStat->SetBinError(i+firstBinEta, graphInvXSectionWeightedAverageEtaStat->GetEYlow()[i]);
            }    

        }    
        if (graphCorrectedYieldWeightedAverageEtaSys){
            graphInvXSectionWeightedAverageEtaSys               = ScaleGraph(graphCorrectedYieldWeightedAverageEtaSys,xSection*recalcBarn);
        }
        
        if (doEtaToPi0 && graphEtaToPi0WeightedAverageStat){
            histoEtaToPi0WeightedAverageStat                    = new TH1D("histoEtaToPi0WeightedAverageStat", "", maxNAllowedEta, binningEta);
            Int_t firstBinEtaToPi0 = 1;
            while (histoEtaToPi0WeightedAverageStat->GetBinCenter(firstBinEtaToPi0) < graphEtaToPi0WeightedAverageStat->GetX()[0]){
                histoEtaToPi0WeightedAverageStat->SetBinContent(firstBinEtaToPi0, 0);
                histoEtaToPi0WeightedAverageStat->SetBinError(firstBinEtaToPi0, 0);

                firstBinEtaToPi0++;
            }
            for (Int_t i = 0; i < graphEtaToPi0WeightedAverageStat->GetN(); i++){
                histoEtaToPi0WeightedAverageStat->SetBinContent(i+firstBinEtaToPi0, graphEtaToPi0WeightedAverageStat->GetY()[i]);
                histoEtaToPi0WeightedAverageStat->SetBinError(i+firstBinEtaToPi0, graphEtaToPi0WeightedAverageStat->GetEYlow()[i]);
            }
        }    
        
    }
    
    const char* fileNameOutputComp = Form("%s/%s_%sResultsFullCorrection_PP.root",outputDir.Data(),isMC.Data(),system.Data());
    TFile* fileOutputForComparisonFullyCorrected = new TFile(fileNameOutputComp,"UPDATE");
        for (Int_t i=0; i< nrOfTrigToBeComb; i++){
            if (histoEventQualtity[i])                  histoEventQualtity[i]->Write(Form("NEvents_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (mode == 2 || mode == 3 || mode == 4 || mode == 5 || mode == 10 || mode == 11 ){
                if (histoTriggerRejection[i])           histoTriggerRejection[i]->Write(Form("TriggRejectMean_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),TObject::kOverwrite);
                if (histoRatioRawClusterPt[i])          histoRatioRawClusterPt[i]->Write(Form("TriggRejectvsPt_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),TObject::kOverwrite);
            }
        }
        fileOutputForComparisonFullyCorrected->mkdir(Form("Pi0%s",optionEnergy.Data()));
        TDirectoryFile* directoryPi0 = (TDirectoryFile*)fileOutputForComparisonFullyCorrected->Get(Form("Pi0%s",optionEnergy.Data())); 
        fileOutputForComparisonFullyCorrected->cd(Form("Pi0%s",optionEnergy.Data()));
        if (graphCorrectedYieldWeightedAveragePi0Stat)  graphCorrectedYieldWeightedAveragePi0Stat->Write("CorrectedYieldPi0",TObject::kOverwrite);
        if (graphCorrectedYieldWeightedAveragePi0Sys)   graphCorrectedYieldWeightedAveragePi0Sys->Write("Pi0SystError",TObject::kOverwrite);
        if (graphInvXSectionWeightedAveragePi0Stat){
            graphInvXSectionWeightedAveragePi0Stat->Write("graphInvCrossSectionPi0",TObject::kOverwrite);
            cout << "graph InvSection pi0 stat" << endl;
            graphInvXSectionWeightedAveragePi0Stat->Print();
        }    
        if (histoInvXSectionWeightedAveragePi0Stat)     histoInvXSectionWeightedAveragePi0Stat->Write("InvCrossSectionPi0",TObject::kOverwrite);
        if (graphInvXSectionWeightedAveragePi0Sys){
            graphInvXSectionWeightedAveragePi0Sys->Write("InvCrossSectionPi0Sys",TObject::kOverwrite);
            cout << "graph InvSection pi0 sys" << endl;
            graphInvXSectionWeightedAveragePi0Sys->Print();
        }    
        if (graphAcceptancePi0Weighted)                 graphAcceptancePi0Weighted->Write("AcceptancePi0", TObject::kOverwrite);
        if (graphEfficiencyPi0Weighted)                 graphEfficiencyPi0Weighted->Write("EfficiencyPi0", TObject::kOverwrite);
        if (graphEffTimesAccPi0Weighted)                graphEffTimesAccPi0Weighted->Write("EffTimesAccPi0", TObject::kOverwrite);
        if (graphPurityPi0Weighted)                     graphPurityPi0Weighted->Write("PurityPi0", TObject::kOverwrite);
        if (mode != 10){
            if (graphMassPi0DataWeighted)                   graphMassPi0DataWeighted->Write("Pi0_Mass_data", TObject::kOverwrite);
            if (graphMassPi0MCWeighted)                     graphMassPi0MCWeighted->Write("Pi0_Mass_MC", TObject::kOverwrite);
            if (graphWidthPi0DataWeighted)                  graphWidthPi0DataWeighted->Write("Pi0_Width_data", TObject::kOverwrite);
            if (graphWidthPi0MCWeighted)                    graphWidthPi0MCWeighted->Write("Pi0_Width_MC", TObject::kOverwrite);            
        }
        for (Int_t i=0; i< nrOfTrigToBeComb; i++){
            cout << "trigger: " << triggerName[i].Data() << endl; 
            if (histoAcceptancePi0[i])                  histoAcceptancePi0[i]->Write(Form("AcceptancePi0_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoEfficiencyPi0[i])                  histoEfficiencyPi0[i]->Write(Form("EfficiencyPi0_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoRawYieldPi0[i])                    histoRawYieldPi0[i]->Write(Form("RAWYieldPerEventsPi0_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoMCInputPi0[i])                     histoMCInputPi0[i]->Write(Form("Pi0_Input_Reweighted_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoEffTimesAccPi0[i])                 histoEffTimesAccPi0[i]->Write(Form("EffTimesAccPi0_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (mode != 10){
                if (histoMassPi0Data[i])                    histoMassPi0Data[i]->Write(Form("Pi0_Mass_data_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoMassPi0MC[i])                      histoMassPi0MC[i]->Write(Form("Pi0_Mass_MC_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoWidthPi0Data[i])                   histoWidthPi0Data[i]->Write(Form("Pi0_Width_data_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoWidthPi0MC[i])                     histoWidthPi0MC[i]->Write(Form("Pi0_Width_MC_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoInvMassSig[i])                     histoInvMassSig[i]->Write(Form("Pi0_InvMassSig_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoInvMassSigPlusBG[i])               histoInvMassSigPlusBG[i]->Write(Form("Pi0_InvMassSigPlusBG_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoInvMassBG[i])                      histoInvMassBG[i]->Write(Form("Pi0_InvMassBG_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (fitInvMassSig[i])                       fitInvMassSig[i]->Write(Form("Pi0_InvMassSigFit_Example_%s",triggerName[i].Data()),TObject::kOverwrite);                
            }
            if (graphsCorrectedYieldRemoved0Pi0[i])     graphsCorrectedYieldRemoved0Pi0[i]->Write(Form("CorrectedYieldPi0_%s",triggerName[i].Data()),TObject::kOverwrite);   
            if (graphsCorrectedYieldSysRemoved0Pi0[i])  graphsCorrectedYieldSysRemoved0Pi0[i]->Write(Form("Pi0SystError_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (enableTriggerEffPi0[i]){
                if (histoTriggerEffPi0[i])              histoTriggerEffPi0[i]->Write(Form("TriggerEfficiencyPi0_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoEffBasePi0[i])                 histoEffBasePi0[i]->Write(Form("EfficiencyBasePi0_%s",triggerName[i].Data()),TObject::kOverwrite);
            }
        }
        
        if (enableEta && mode != 10){
            fileOutputForComparisonFullyCorrected->mkdir(Form("Eta%s",optionEnergy.Data()));
            TDirectoryFile* directoryEta = (TDirectoryFile*)fileOutputForComparisonFullyCorrected->Get(Form("Eta%s",optionEnergy.Data())); 
            fileOutputForComparisonFullyCorrected->cd(Form("Eta%s",optionEnergy.Data()));
            if (graphCorrectedYieldWeightedAverageEtaStat)  graphCorrectedYieldWeightedAverageEtaStat->Write("CorrectedYieldEta",TObject::kOverwrite);
            if (graphCorrectedYieldWeightedAverageEtaSys)   graphCorrectedYieldWeightedAverageEtaSys->Write("EtaSystError",TObject::kOverwrite);
            if (graphInvXSectionWeightedAverageEtaStat){
                graphInvXSectionWeightedAverageEtaStat->Write("graphInvCrossSectionEta",TObject::kOverwrite);
                cout << "graph InvSection eta stat" << endl;
                graphInvXSectionWeightedAverageEtaStat->Print();
            }
            if (histoInvXSectionWeightedAverageEtaStat)     histoInvXSectionWeightedAverageEtaStat->Write("InvCrossSectionEta",TObject::kOverwrite);
            if (graphInvXSectionWeightedAverageEtaSys){
                graphInvXSectionWeightedAverageEtaSys->Write("InvCrossSectionEtaSys",TObject::kOverwrite);
                cout << "graph InvSection eta stat" << endl;
                graphInvXSectionWeightedAverageEtaSys->Print();
            }    
            if (graphAcceptanceEtaWeighted)                 graphAcceptanceEtaWeighted->Write("AcceptanceEta", TObject::kOverwrite);
            if (graphEfficiencyEtaWeighted)                 graphEfficiencyEtaWeighted->Write("EfficiencyEta", TObject::kOverwrite);
            if (graphEffTimesAccEtaWeighted)                graphEffTimesAccEtaWeighted->Write("EffTimesAccEta", TObject::kOverwrite);
            if (graphMassEtaDataWeighted)                   graphMassEtaDataWeighted->Write("Eta_Mass_data", TObject::kOverwrite);
            if (graphMassEtaMCWeighted)                     graphMassEtaMCWeighted->Write("Eta_Mass_MC", TObject::kOverwrite);
            if (graphWidthEtaDataWeighted)                  graphWidthEtaDataWeighted->Write("Eta_Width_data", TObject::kOverwrite);
            if (graphWidthEtaMCWeighted)                    graphWidthEtaMCWeighted->Write("Eta_Width_MC", TObject::kOverwrite);

            for (Int_t i=0; i< nrOfTrigToBeComb; i++){
                if (histoAcceptanceEta[i])                  histoAcceptanceEta[i]->Write(Form("AcceptanceEta_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoEfficiencyEta[i])                  histoEfficiencyEta[i]->Write(Form("EfficiencyEta_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoRawYieldEta[i])                    histoRawYieldEta[i]->Write(Form("RAWYieldPerEventsEta_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoMCInputEta[i])                     histoMCInputEta[i]->Write(Form("Eta_Input_Reweighted_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoEffTimesAccEta[i])                 histoEffTimesAccEta[i]->Write(Form("EffTimesAccEta_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoMassEtaData[i])                    histoMassEtaData[i]->Write(Form("Eta_Mass_data_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoMassEtaMC[i])                      histoMassEtaMC[i]->Write(Form("Eta_Mass_MC_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoWidthEtaData[i])                   histoWidthEtaData[i]->Write(Form("Eta_Width_data_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoWidthEtaMC[i])                     histoWidthEtaMC[i]->Write(Form("Eta_Width_MC_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoEtaInvMassSig[i])                  histoEtaInvMassSig[i]->Write(Form("Eta_InvMassSig_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoEtaInvMassSigPlusBG[i])            histoEtaInvMassSigPlusBG[i]->Write(Form("Eta_InvMassSigPlusBG_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoEtaInvMassBG[i])                   histoEtaInvMassBG[i]->Write(Form("Eta_InvMassBG_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (fitEtaInvMassSig[i])                    fitEtaInvMassSig[i]->Write(Form("Eta_InvMassSigFit_Example_%s",triggerName[i].Data()),TObject::kOverwrite);                

                if (graphsCorrectedYieldRemoved0Eta[i])     graphsCorrectedYieldRemoved0Eta[i]->Write(Form("CorrectedYieldEta_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (graphsCorrectedYieldSysRemoved0Eta[i])  graphsCorrectedYieldSysRemoved0Eta[i]->Write(Form("EtaSystError_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (enableTriggerEffEta[i]){
                    if (histoTriggerEffEta[i])              histoTriggerEffEta[i]->Write(Form("TriggerEfficiencyEta_%s",triggerName[i].Data()),TObject::kOverwrite);
                    if (histoEffBaseEta[i])                 histoEffBaseEta[i]->Write(Form("EfficiencyBaseEta_%s",triggerName[i].Data()),TObject::kOverwrite);
                }
            }
            if (doEtaToPi0){
                if (histoEtaToPi0WeightedAverageStat)       histoEtaToPi0WeightedAverageStat->Write(Form("EtaToPi0%sStatError",addNameBinshift.Data()),TObject::kOverwrite);
                if (graphEtaToPi0WeightedAverageStat)       graphEtaToPi0WeightedAverageStat->Write(Form("graphEtaToPi0%sStatError",addNameBinshift.Data()),TObject::kOverwrite);
                if (graphEtaToPi0WeightedAverageSys )       graphEtaToPi0WeightedAverageSys->Write(Form("EtaToPi0%sSystError",addNameBinshift.Data()),TObject::kOverwrite);
                for (Int_t i=0; i< nrOfTrigToBeComb; i++){
                    if (graphsEtaToPi0Removed0[i])          graphsEtaToPi0Removed0[i]->Write(Form("EtaToPi0%sStatError_%s",addNameBinshift.Data(), triggerName[i].Data()),TObject::kOverwrite);
                    if (graphsEtaToPi0SysRemoved0[i])       graphsEtaToPi0SysRemoved0[i]->Write(Form("EtaToPi0%sSystError_%s",addNameBinshift.Data(), triggerName[i].Data()),TObject::kOverwrite);
                }
            }
        }
        
    fileOutputForComparisonFullyCorrected->Write();
    fileOutputForComparisonFullyCorrected->Close();
    
    cout << WriteParameterToFile(fitInvYieldPi0) << endl;
}
