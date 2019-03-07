/****************************************************************************************************************************
 ******     provided by Gamma Conversion Group, PWGGA,                                                                  *****
 ******        Friederike Bock, friederike.bock@cern.ch                                                                 *****
 ******        Daniel Muehlheim, d.muehlheim@cern.ch                                                                    *****
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
void  ProduceFinalResultsPiPlPiMiPiZero(   TString fileListNameOmega     = "triggerFileListOmega.txt",
                                            Int_t   mode                = 4,
                                            Int_t   numberOfTrigg       = 6,
                                            TString suffix              = "eps",
                                            TString isMC                = "",
                                            TString optionEnergy        = "",
                                            TString period              = "",
                                            Bool_t  pileUpApplied       = kTRUE,
                                            Float_t maxPtGlobalOmega      = 16.,
                                            Bool_t  averagedOmega         = kFALSE,
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

    Double_t binningOmega[100];
    Int_t maxNBinsOmegaAbs          = 0;
    Int_t maxNBinsOmega             = GetBinning( binningOmega, maxNBinsOmegaAbs, "Omega", optionEnergy, mode );
    Int_t nRealTriggers             = 0;
    while (binningOmega[maxNBinsOmega] < maxPtGlobalOmega ) maxNBinsOmega++;
    for (Int_t i= 0; i< maxNBinsOmega+1; i++){
        cout << binningOmega[i] << ", ";
    }
    cout << endl;

    Double_t binningEta[100];
    Int_t maxNBinsEtaAbs            = 0;
    Int_t maxNBinsEta               = GetBinning( binningEta, maxNBinsEtaAbs, "Eta", optionEnergy, mode );
    Int_t maxNAllowedEta            = maxNBinsEta;
    while (binningEta[maxNAllowedEta] < maxPtGlobalEta ) maxNAllowedEta++;
    for (Int_t i= 0; i< maxNAllowedEta+1; i++){
        cout << binningEta[i] << ", ";
    }
    cout << endl;

    Double_t maxPtGlobalCluster     = 25;
    if (optionEnergy.CompareTo("2.76TeV")==0){
        if (mode==2){
            maxPtGlobalCluster          = 26;
        } else if (mode == 4){
            maxPtGlobalCluster          = 31;
        }
    } else if (optionEnergy.CompareTo("8TeV")==0){
      if(mode==2 || mode==4){
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

    //***************************************************************************************************************
    //********************************** setting correct histo names ************************************************
    //***************************************************************************************************************
    // Variable to quickly change which type of Inv mass plot is used
    TString InvMassTypeEnding = "_SubPiZero";

    TString nameCorrectedYield                          = Form("CorrectedYieldTrueEff%s",InvMassTypeEnding.Data());
    TString nameEfficiency                              = Form("TrueMesonEffiPt%s",InvMassTypeEnding.Data());
    TString nameAcceptance                              = "fMCMesonAccepPt";
    TString nameAcceptanceWOEvtWeights                  = "fMCMesonAccepPtWOEvtWeights";
    TString nameMassMC                                  = Form("histoTrueMassMeson%s",InvMassTypeEnding.Data());
    TString nameWidthMC                                 = Form("histoTrueFWHMMeson%s",InvMassTypeEnding.Data());
  //  TString nameMCYield                                 = "MCYield_Meson_oldBinWOWeights";
    TString nameMCYield                                 = "MCYield_Meson";
    if ( mode == 4 || mode == 5 ){
        nameCorrectedYield                              = "CorrectedYieldNormEff";
        nameEfficiency                                  = "MesonEffiPt";
        nameMassMC                                      = "histoMassMesonRecMC";
        nameWidthMC                                     = "histoFWHMMesonRecMC";
    } else if ( (mode == 2 || mode == 3) && !(optionEnergy.CompareTo("8TeV")==0 || optionEnergy.CompareTo("pPb_5.023TeV")==0)){
        cout << "using rec quantities for PCM-EMC/PCM-PHOS" << endl;
        nameMassMC                                      = "histoMassMesonRecMC";
        nameWidthMC                                     = "histoFWHMMesonRecMC";
    }

    //***************************************************************************************************************
    //********************************** settings for secondary pion corrs ******************************************
    //***************************************************************************************************************
    // Note: secondary correction is not used for omega! However, this is left in here in case it is needed later
                                                        //    0      1     2        3     4       5     6     7     8     9     10      11    12    2
    Double_t maxYEffSecCorr[4][14]                      = { { 0.080, 0.3 , 0.03000, 0.04, 0.0600, 0.12, 0.3 , 0.3 , 0.3 , 0.3 , 0.15,   0.15, 0.12, 0.04},
                                                            { 0.001, 0.3 , 0.00010, 0.04, 0.0050, 0.12, 0.3 , 0.3 , 0.3 , 0.3 , 0.002,  0.15, 0.12, 0.04},
                                                            { 0.001, 0.3 , 0.00006, 0.04, 0.0003, 0.12, 0.3 , 0.3 , 0.3 , 0.3 , 0.000003,  0.15, 0.12, 0.04},
                                                            { 0.030, 0.3 , 0.01800, 0.04, 0.0300, 0.12, 0.3 , 0.3 , 0.3 , 0.3 , 0.06,   0.15, 0.12, 0.04} };
    Color_t colorSec[4]                                 = {kRed+2, kCyan+2, 807, kBlue};
    Style_t markerStyleSec[4]                           = {21, 33, 30, 28};
    Size_t markerSizeSec[4]                             = {1.5, 1.75, 2., 1.5};

    TString nameSecOmegaPartRead[4]                       = {"K0S", "K0L", "Lambda", "Rest"};
    TString nameSecOmegaPartLabel[4]                      = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "had. int."};
    TString nameSecEfficiency[4]                        = {"TrueSecFromK0SEffiPt", "TrueSecFromK0LEffiPt", "TrueSecFromLambdaEffiPt", "TrueSecFromRestEffiPt"};
    //***************************************************************************************************************


    TString cutNumber       [MaxNumberOfFiles];
    TString triggerName     [MaxNumberOfFiles];
    TString triggerNameLabel[MaxNumberOfFiles];
    Float_t minPt           [MaxNumberOfFiles];
    Float_t maxPt           [MaxNumberOfFiles];
    Int_t trigSteps         [MaxNumberOfFiles][3];
    Float_t ptFromSpecOmega   [MaxNumberOfFiles][2];
    Float_t ptFromSpecEta   [MaxNumberOfFiles][2];
    Bool_t maskedFullyOmega   [MaxNumberOfFiles];
    Bool_t maskedFullyEta   [MaxNumberOfFiles];
    TString sysFileOmega      [MaxNumberOfFiles];
    TString sysFileEta      [MaxNumberOfFiles];
    TString sysFileEtaToOmega [MaxNumberOfFiles];
    TString cutNumberBaseEff[MaxNumberOfFiles];
    //***************************************************************************************************************
    //*************************** read setting from configuration file **********************************************
    //***************************************************************************************************************
    ifstream in(fileListNameOmega.Data());
    cout<<"Available Triggers:"<<endl;
    // general number of triggers set
    Int_t nrOfTrigToBeComb      = 0;
    // number of triggers which are really used for the respective analysis
    Int_t nrOfTrigToBeCombOmegaRed        = 0;
    Int_t nrOfTrigToBeCombEtaRed        = 0;
    Int_t nrOfTrigToBeCombEtaToOmegaRed   = 0;
    while(!in.eof() && nrOfTrigToBeComb<numberOfTrigg ){
        in >> cutNumber[nrOfTrigToBeComb] >> minPt[nrOfTrigToBeComb] >> maxPt[nrOfTrigToBeComb] >> triggerName[nrOfTrigToBeComb]
        >> trigSteps[nrOfTrigToBeComb][0]  >> trigSteps[nrOfTrigToBeComb][1]  >> trigSteps[nrOfTrigToBeComb][2] >> ptFromSpecOmega[nrOfTrigToBeComb][0] >> ptFromSpecOmega[nrOfTrigToBeComb][1]
        >> ptFromSpecEta[nrOfTrigToBeComb][0] >> ptFromSpecEta[nrOfTrigToBeComb][1] >> sysFileOmega[nrOfTrigToBeComb] >> sysFileEta[nrOfTrigToBeComb] >> sysFileEtaToOmega[nrOfTrigToBeComb] >>
        cutNumberBaseEff[nrOfTrigToBeComb];
        cout<< cutNumber[nrOfTrigToBeComb]<< "\t"<< triggerName[nrOfTrigToBeComb] << "\t transverse momentum range: " << minPt[nrOfTrigToBeComb]<< "\t to "<< maxPt[nrOfTrigToBeComb] <<endl;
        cout << trigSteps[nrOfTrigToBeComb][0] << "\t" << trigSteps[nrOfTrigToBeComb][1] << "\t"<< trigSteps[nrOfTrigToBeComb][2] << endl;
        nrOfTrigToBeComb++;
        cout << cutNumberBaseEff[nrOfTrigToBeComb] << endl;
    }

    for (Int_t i = 0; i < nrOfTrigToBeComb; i++){
        // figure out which triggers are fully masked for the Omega
        if ((ptFromSpecOmega[i][1] == -1 && ptFromSpecOmega[i][0] == -1 )|| (ptFromSpecOmega[i][0] == 0 && ptFromSpecOmega[i][1] == 0)){
            maskedFullyOmega[i]       = kTRUE;
        } else {
            maskedFullyOmega[i]       = kFALSE;
        }
        // figure out which triggers are fully masked for the eta
        if ((ptFromSpecEta[i][1] == -1 && ptFromSpecEta[i][0] == -1 )|| (ptFromSpecEta[i][0] == 0 && ptFromSpecEta[i][1] == 0)){
            maskedFullyEta[i]       = kTRUE;
        } else {
            maskedFullyEta[i]       = kFALSE;
        }
    }

    Double_t minPtGlobalOmega         = ptFromSpecOmega[0][0];
    Double_t minPtGlobalEta         = ptFromSpecEta[0][0];
    if (maskedFullyOmega[0])
        minPtGlobalOmega              = 12;
    if (maskedFullyEta[0])
        minPtGlobalEta              = 12;


    for (Int_t j = 1; j < nrOfTrigToBeComb; j++){
        if (minPtGlobalOmega > ptFromSpecOmega[j][0] && !maskedFullyOmega[j] )
            minPtGlobalOmega = ptFromSpecOmega[j][0];
        if (minPtGlobalEta > ptFromSpecEta[j][0] && !maskedFullyEta[j] )
            minPtGlobalEta = ptFromSpecEta[j][0];
    }
    cout << "global minimum pT for Omega: " << minPtGlobalOmega << endl;
    cout << "global minimum pT for eta: " << minPtGlobalEta << endl;

    // set individual triggers to total number of availabel triggers, will be reduced later according to usage
    nrOfTrigToBeCombOmegaRed      = nrOfTrigToBeComb;
    nrOfTrigToBeCombEtaRed      = nrOfTrigToBeComb;
    nrOfTrigToBeCombEtaToOmegaRed = nrOfTrigToBeComb;

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

    }

    // defining output directory
    TString outputDir =    Form("%s/%s/%s/FinalResultsTriggersPatched%s", suffix.Data(),optionEnergy.Data(),dateForOutput.Data(),fNLMStringOutput.Data());
    if(optionEnergy.CompareTo("900GeV") == 0 || optionEnergy.CompareTo("7TeV") == 0 || optionEnergy.CompareTo("8TeV") == 0){
      if(mode == 4) outputDir = Form("%s/%s/%s/FinalResultsTriggersPatched_EMCAL%s", suffix.Data(),optionEnergy.Data(),dateForOutput.Data(),fNLMStringOutput.Data());
      if(mode == 2) outputDir = Form("%s/%s/%s/FinalResultsTriggersPatched_PCMEMCAL%s", suffix.Data(),optionEnergy.Data(),dateForOutput.Data(),fNLMStringOutput.Data());
    }
    gSystem->Exec("mkdir -p "+outputDir);

    gSystem->Exec(Form("cp %s %s/configurationFile.txt", fileListNameOmega.Data(), outputDir.Data()));

    TString nameFinalResDat                     = Form("%s/FitResults.dat",outputDir.Data());
    fstream  fileFitsOutput;
    fileFitsOutput.open(nameFinalResDat.Data(), ios::out);

    vector<TString>** ptSysDetail     = new vector<TString>*[MaxNumberOfFiles];
    for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++) ptSysDetail[iR] = new vector<TString>[50];
    TString             sysStringComb = "PCM";
    if(mode == 2)       sysStringComb = "PCMEMC";
    else if(mode == 3)  sysStringComb = "PCMPHOS";
    else if(mode == 4)  sysStringComb = "EMCEMC";
    else if(mode == 5)  sysStringComb = "PHOS";
    else if(mode == 40 || mode == 60) sysStringComb = "PCM";
    else if(mode == 41 || mode == 61) sysStringComb = "PCMEMC";
    else if(mode == 42 || mode == 62) sysStringComb = "PCMPHOS";
    else if(mode == 44 || mode == 64) sysStringComb = "EMC";
    else if(mode == 45 || mode == 65) sysStringComb = "PHOS";

    //***************************************************************************************************************
    //******************************** Load Omega histograms **********************************************************
    //***************************************************************************************************************
    TString FileNameCorrectedOmega    [MaxNumberOfFiles];
    TFile*  fileCorrectedOmega        [MaxNumberOfFiles];
    TString FileNameUnCorrectedOmega  [MaxNumberOfFiles];
    TFile*  fileUnCorrectedOmega      [MaxNumberOfFiles];
    TString FileNameUnCorrectedMCOmega[MaxNumberOfFiles];
    TFile*  fileUnCorrectedMCOmega    [MaxNumberOfFiles];

    TH1D*   histoCorrectedYieldOmega  [MaxNumberOfFiles];
    TH1D*   histoRawClusterPt       [MaxNumberOfFiles];
    TH1D*   histoRawClusterE        [MaxNumberOfFiles];
    TH1D*   histoEfficiencyOmega      [MaxNumberOfFiles];
    TH1D*   histoAcceptanceOmega      [MaxNumberOfFiles];

    TH1D*   histoAcceptanceOmegaWOEvtWeights [MaxNumberOfFiles];
    TH1D*   histoPurityOmega          [MaxNumberOfFiles];
    TH1D*   histoEffTimesAccOmega     [MaxNumberOfFiles];
    TH1D*   histoRawYieldOmega        [MaxNumberOfFiles];
    TH1D*   histoRatioRawClusterPt  [MaxNumberOfFiles];
    TH1D*   histoRatioRawClusterE   [MaxNumberOfFiles];
    TH1D*   histoTriggerRejection   [MaxNumberOfFiles];
    TH1D*   histoMassOmegaData        [MaxNumberOfFiles];
    TH1D*   histoMassOmegaMC          [MaxNumberOfFiles];
    TH1D*   histoWidthOmegaData       [MaxNumberOfFiles];
    TH1D*   histoWidthOmegaMC         [MaxNumberOfFiles];
    TH1F*   histoEventQualtity      [MaxNumberOfFiles];
    TH1D*   histoMCInputOmega         [MaxNumberOfFiles];

    Double_t triggRejecFac          [MaxNumberOfFiles][MaxNumberOfFiles];
    Double_t triggRejecFacErr       [MaxNumberOfFiles][MaxNumberOfFiles];

    TString FileNameEffBaseOmega      [MaxNumberOfFiles];
    TFile*  fileEffBaseOmega          [MaxNumberOfFiles];
    TH1D*   histoEffBaseOmega         [MaxNumberOfFiles];
    TH1D*   histoTriggerEffOmega      [MaxNumberOfFiles];
    Bool_t  enableTriggerEffOmega     [MaxNumberOfFiles];
    Bool_t  enableTriggerEffOmegaAll                          = kFALSE;
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
    TH1D*   histoSecEffiOmegaFromX  [4][MaxNumberOfFiles];
    TH1D*   histoEffectCorrOmegaFromX[4][MaxNumberOfFiles];
    Bool_t hasSecEffi[4]                                    = {kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t hasSecCorrFac[4]                                 = {kFALSE, kFALSE, kFALSE, kFALSE};

    TH1D* histoDataM02              [MaxNumberOfFiles];
    TH1D* histoMCrecM02             [MaxNumberOfFiles];
    TH1D* histoTrueOmegaM02           [MaxNumberOfFiles];
    TH1D* histoTrueEtaM02           [MaxNumberOfFiles];
    TH1D* histoTrueGammaM02         [MaxNumberOfFiles];
    TH1D* histoTrueElectronM02      [MaxNumberOfFiles];
    TH1D* histoTrueBGM02            [MaxNumberOfFiles];
    TH1D* histoTrueOmegaPureMerged    [MaxNumberOfFiles];
    TH1D* histoTrueOmegaPConvMerged   [MaxNumberOfFiles];
    TH1D* histoTrueOmegaOneGammaM02   [MaxNumberOfFiles];
    TH1D* histoTrueOmegaOneElectronM02[MaxNumberOfFiles];

    // check if fit file for binshifting has to be adjusted for every energy
    TF1* fitBinShiftOmegaTCM                          = 0x0;
    TF1* fitBinShiftOmega                             = 0x0;
    TF1* fitBinShiftEtaTCM                          = 0x0;
    TF1* fitBinShiftEta                             = 0x0;
    Bool_t doBinShiftForEtaToOmega                    = kFALSE;
    TString addNameBinshift                         = "";
    if (nameFileFitsShift.CompareTo("") != 0){
        doBinShiftForEtaToOmega                       = kTRUE;
        addNameBinshift                             = "YShifted";
    }

    if (doBinShiftForEtaToOmega){
        TFile *fileFitsBinShift                         = new TFile(nameFileFitsShift);
        fitBinShiftOmega                                  = (TF1*)fileFitsBinShift->Get("TsallisFitOmega");
        if(!fitBinShiftOmega || optionEnergy.CompareTo("7TeV")==0){
            fitBinShiftOmega                              = (TF1*)fileFitsBinShift->Get("Omega7TeV/Fits/fitBinShiftingOmega");
            if(!fitBinShiftOmega) fitBinShiftOmega          = (TF1*)fileFitsBinShift->Get("Omega7TeV/TsallisFitOmega");
            fitBinShiftOmegaTCM                           = (TF1*)fileFitsBinShift->Get("Omega7TeV/Fits/fitBinShiftingOmega");
            if(!fitBinShiftOmegaTCM) fitBinShiftOmegaTCM    = (TF1*)fileFitsBinShift->Get("Omega7TeV/TwoComponentModelFitOmega");
        }
        fitBinShiftEta                                  = (TF1*)fileFitsBinShift->Get("TsallisFitEta");
        if(!fitBinShiftEta || optionEnergy.CompareTo("7TeV")==0){
            fitBinShiftEta                              = (TF1*)fileFitsBinShift->Get("Eta7TeV/Fits/fitBinShiftingEta");
            if(!fitBinShiftEta) fitBinShiftEta          = (TF1*)fileFitsBinShift->Get("Eta7TeV/TsallisFitEta");
            fitBinShiftEtaTCM                           = (TF1*)fileFitsBinShift->Get("Eta7TeV/Fits/fitBinShiftingEta");
            if(!fitBinShiftEtaTCM) fitBinShiftEtaTCM    = (TF1*)fileFitsBinShift->Get("Eta7TeV/TwoComponentModelFitEta");
        }
        if(!fitBinShiftOmega || optionEnergy.CompareTo("8TeV")==0){
            fitBinShiftOmega                              = (TF1*)fileFitsBinShift->Get("Omega8TeV/TsallisFitOmega");
            fitBinShiftOmegaTCM                           = (TF1*)fileFitsBinShift->Get("Omega8TeV/TwoComponentModelFitOmega");
        }
        if(!fitBinShiftEta || optionEnergy.CompareTo("8TeV")==0){
            fitBinShiftEta                              = (TF1*)fileFitsBinShift->Get("Eta8TeV/TsallisFitEta");
            fitBinShiftEtaTCM                           = (TF1*)fileFitsBinShift->Get("Eta8TeV/TwoComponentModelFitEta");
        }
        if( optionEnergy.CompareTo("pPb_5.023TeV")==0){
            fitBinShiftOmega                              = (TF1*)fileFitsBinShift->Get("TwoComponentModelFitOmega");
            fitBinShiftOmegaTCM                           = (TF1*)fileFitsBinShift->Get("TwoComponentModelFitOmega");
            fitBinShiftEta                              = (TF1*)fileFitsBinShift->Get("TwoComponentModelFitEta");
        }
        cout << fitBinShiftOmega << " - " << fitBinShiftEta << endl;
        cout << "fits for shifting found " << endl;
    }

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
        Int_t exampleBin                                    = ReturnSingleInvariantMassBinPlotting ("Omega", optionEnergy, mode, trigger.Atoi(), scaleFacSinBin);

        FileNameCorrectedOmega[i]                             = Form("%s/%s/Omega_%s_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(),
                                                                    cutNumber[i].Data());
        cout<< FileNameCorrectedOmega[i] << endl;
        fileCorrectedOmega[i]                                 = new TFile(FileNameCorrectedOmega[i]);
        if (fileCorrectedOmega[i]->IsZombie()) return;
        // read uncorrected file
        FileNameUnCorrectedOmega[i]                           = Form("%s/%s/Omega_%s_GammaConvV1WithoutCorrection_%s.root",cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(),
                                                                    cutNumber[i].Data());

        cout<< FileNameUnCorrectedOmega[i] << endl;
        fileUnCorrectedOmega[i]                               = new TFile(FileNameUnCorrectedOmega[i]);
        if (fileUnCorrectedOmega[i]->IsZombie()) return;

        if (isMC.CompareTo("data") == 0){
            FileNameUnCorrectedMCOmega[i]                     = Form("%s/%s/Omega_MC_GammaConvV1WithoutCorrection_%s.root",cutNumber[i].Data(), optionEnergy.Data(), cutNumber[i].Data());

            cout<< FileNameUnCorrectedMCOmega[i] << endl;
            fileUnCorrectedMCOmega[i]                         = new TFile(FileNameUnCorrectedMCOmega[i]);
            if (fileUnCorrectedMCOmega[i]->IsZombie())
                enableTriggerRejecCompMC                    = kFALSE;
            else
                enableTriggerRejecCompMC                    = kTRUE;
        }
        deltaRapid[i]                                       =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
        cout << "For Trigger: " << triggerName[i] << " using rapidity: " <<  rapidityRange.Data() << endl;

        histoEventQualtity[i]                               = (TH1F*)fileCorrectedOmega[i]->Get("NEvents");
        histoCorrectedYieldOmega[i]                           = (TH1D*)fileCorrectedOmega[i]->Get(nameCorrectedYield.Data());
        histoCorrectedYieldOmega[i]->SetName(Form("CorrectedYield_%s",cutNumber[i].Data()));
        histoEfficiencyOmega[i]                               = (TH1D*)fileCorrectedOmega[i]->Get(nameEfficiency.Data());
        histoEfficiencyOmega[i]->SetName(Form("Efficiency_%s",  cutNumber[i].Data()));
        histoAcceptanceOmega[i]                               = (TH1D*)fileCorrectedOmega[i]->Get(nameAcceptance.Data());
        histoAcceptanceOmega[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
        histoAcceptanceOmegaWOEvtWeights[i]                   = (TH1D*)fileCorrectedOmega[i]->Get(nameAcceptanceWOEvtWeights.Data());

        for (Int_t k = 0; k < 4; k++){
            histoEffectCorrOmegaFromX[k][i]                   = NULL;
            histoEffectCorrOmegaFromX[k][i]                   = (TH1D*)fileCorrectedOmega[i]->Get(Form("RatioSecYieldFrom%sMesonFromCocktailToRaw",nameSecOmegaPartRead[k].Data()));
            if (!histoEffectCorrOmegaFromX[k][i])
            histoEffectCorrOmegaFromX[k][i]                   = (TH1D*)fileCorrectedOmega[i]->Get(Form("RatioSecYieldFrom%sMesonFromToyToRaw",nameSecOmegaPartRead[k].Data()));
            if (!histoEffectCorrOmegaFromX[k][i])
            histoEffectCorrOmegaFromX[k][i]                   = (TH1D*)fileCorrectedOmega[i]->Get(Form("RatioSecYieldFrom%sMesonToRaw",nameSecOmegaPartRead[k].Data()));

            if (!histoEffectCorrOmegaFromX[k][i] || isMC.CompareTo("MC") == 0 )
                histoEffectCorrOmegaFromX[k][i]               = (TH1D*)fileCorrectedOmega[i]->Get(Form("RatioSecYieldFrom%sMesonToRaw",nameSecOmegaPartRead[k].Data()));
            if (histoEffectCorrOmegaFromX[k][i])
                hasSecCorrFac[k]                            = kTRUE;
            histoSecEffiOmegaFromX[k][i]                      = NULL;
            histoSecEffiOmegaFromX[k][i]                      = (TH1D*)fileCorrectedOmega[i]->Get(nameSecEfficiency[k].Data());
            if (histoSecEffiOmegaFromX[k][i])
                hasSecEffi[k]                               = kTRUE;
        }

        if(histoAcceptanceOmegaWOEvtWeights[i]) histoAcceptanceOmegaWOEvtWeights[i]->SetName(Form("AcceptanceWOEvtWeights_%s",  cutNumber[i].Data()));

        histoRawYieldOmega[i]                                 = (TH1D*)fileCorrectedOmega[i]->Get(Form("histoYieldMesonPerEvent%s",InvMassTypeEnding.Data()));
        histoRawYieldOmega[i]->SetName(Form("RAWYieldPerEvent_%s",cutNumber[i].Data()));
        if (mode != 10){
            histoMassOmegaData[i]                                 = (TH1D*)fileCorrectedOmega[i]->Get(Form("histoMassMeson%s",InvMassTypeEnding.Data()));
            histoMassOmegaData[i]->SetName(Form("Omega_Mass_data_%s",cutNumber[i].Data()));
            histoMassOmegaMC[i]                                   = (TH1D*)fileCorrectedOmega[i]->Get(nameMassMC.Data());
            histoMassOmegaMC[i]->SetName(Form("Omega_Mass_MC_%s",cutNumber[i].Data()));
            histoWidthOmegaData[i]                                = (TH1D*)fileCorrectedOmega[i]->Get(Form("histoFWHMMeson%s",InvMassTypeEnding.Data()));
            histoWidthOmegaData[i]->SetName(Form("Omega_Width_data_%s",cutNumber[i].Data()));
            histoWidthOmegaMC[i]                                  = (TH1D*)fileCorrectedOmega[i]->Get(nameWidthMC.Data());
            histoWidthOmegaMC[i]->SetName(Form("Omega_Width_MC_%s",cutNumber[i].Data()));
            histoInvMassSig[i]                                  = (TH1D*)fileCorrectedOmega[i]->Get(Form("InvMassSig_PtBin%02d",exampleBin));
            if (histoInvMassSig[i]) histoInvMassSig[i]->SetName(Form("Omega_InvMassSig_Example_%s",triggerName[i].Data()));
            histoInvMassSigPlusBG[i]                            = (TH1D*)fileCorrectedOmega[i]->Get(Form("InvMassSigPlusBG_PtBin%02d",exampleBin));
            if (histoInvMassSigPlusBG[i]) histoInvMassSigPlusBG[i]->SetName(Form("Omega_InvMassSigPlusBG_Example_%s",triggerName[i].Data()));
            histoInvMassBG[i]                                   = (TH1D*)fileCorrectedOmega[i]->Get(Form("InvMassBG_PtBin%02d",exampleBin));
            if (histoInvMassBG[i]) histoInvMassBG[i]->SetName(Form("Omega_InvMassBG_Example_%s",triggerName[i].Data()));
            fitInvMassSig[i]                                    = (TF1*)fileCorrectedOmega[i]->Get(Form("FitInvMassSig_PtBin%02d",exampleBin));
            if (fitInvMassSig[i]) fitInvMassSig[i]->SetName(Form("Omega_InvMassSigFit_Example_%s",triggerName[i].Data()));
        }
        if (cutNumberBaseEff[i].CompareTo("bla") != 0){
            FileNameEffBaseOmega[i]                           = Form("%s/%s/Omega_MC_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), cutNumberBaseEff[i].Data());
            fileEffBaseOmega[i]                               = new TFile(FileNameEffBaseOmega[i]);
            if (fileEffBaseOmega[i]->IsZombie()){
                enableTriggerEffOmega[i]                         = kFALSE;
                cout << "Didn't find the effi base file for " << triggerName[i].Data() << endl;
                cout << "ABORTING: as base effi file was requested" << endl;

            } else {
                enableTriggerEffOmega[i]                         = kTRUE;
                enableTriggerEffOmegaAll                         = kTRUE;
            }
            if (enableTriggerEffOmega[i]){
                TString effiNameBase                        = "TrueMesonEffiPt";
                TH1D* histoEffiOmegaTemp                      = (TH1D*)fileCorrectedOmega[i]->Get(effiNameBase.Data());
                TH1D* histoEffiBaseOmegaTemp                  = (TH1D*)fileEffBaseOmega[i]->Get(effiNameBase.Data());
                histoEffBaseOmega[i]                          = (TH1D*)fileEffBaseOmega[i]->Get(nameEfficiency.Data());
                histoEffBaseOmega[i]->SetName(Form("EfficiencyBase_%s",  cutNumber[i].Data()));
                histoTriggerEffOmega[i]                       = (TH1D*)histoEffiOmegaTemp->Clone(Form("TriggerEfficiency_%s", cutNumber[i].Data()));
                histoTriggerEffOmega[i]->Divide(histoTriggerEffOmega[i],histoEffiBaseOmegaTemp,1.,1.,"B");

//                //limit trigger efficiency to 1
                if(optionEnergy.CompareTo("8TeV") == 0){
                  for(Int_t j = 1; j<histoTriggerEffOmega[i]->GetNbinsX()+1; j++){
                    Double_t binC = histoTriggerEffOmega[i]->GetBinContent(j);
                    if(binC > 1.){
                      histoEffBaseOmega[i]->SetBinContent(j, histoEffBaseOmega[i]->GetBinContent(j)*binC);
                      histoTriggerEffOmega[i]->SetBinContent(j, 1.);
                    }
                  }
                }

                histoEffTimesAccOmega[i]                      = (TH1D*)histoEffBaseOmega[i]->Clone(Form("EffTimeAcc_%s",  cutNumber[i].Data()));
                if(histoAcceptanceOmegaWOEvtWeights[i]){
                  histoAcceptanceOmegaWOEvtWeights[i]->Sumw2();
                  histoEffTimesAccOmega[i]->Multiply(histoAcceptanceOmegaWOEvtWeights[i]);
                  histoEfficiencyOmega[i]->Multiply(histoAcceptanceOmega[i]);
                  histoEfficiencyOmega[i]->Divide(histoEfficiencyOmega[i],histoAcceptanceOmegaWOEvtWeights[i],1.,1.,"B");
                  histoAcceptanceOmega[i] = histoAcceptanceOmegaWOEvtWeights[i];
                  histoAcceptanceOmega[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
                }else histoEffTimesAccOmega[i]->Multiply(histoAcceptanceOmega[i]);
                histoEffTimesAccOmega[i]->Scale(deltaRapid[i]*2*TMath::Pi());
            } else {
                histoEffBaseOmega[i]                          = NULL;
                histoTriggerEffOmega[i]                       = NULL;
            }
        } else {
            histoEffTimesAccOmega[i]                      = (TH1D*)histoEfficiencyOmega[i]->Clone(Form("EffTimeAcc_%s",  cutNumber[i].Data()));
            if(histoAcceptanceOmegaWOEvtWeights[i]){
              histoAcceptanceOmegaWOEvtWeights[i]->Sumw2();
              histoEffTimesAccOmega[i]->Multiply(histoAcceptanceOmegaWOEvtWeights[i]);
              histoEfficiencyOmega[i]->Multiply(histoAcceptanceOmega[i]);
              histoEfficiencyOmega[i]->Divide(histoEfficiencyOmega[i],histoAcceptanceOmegaWOEvtWeights[i],1.,1.,"B");
              histoAcceptanceOmega[i] = histoAcceptanceOmegaWOEvtWeights[i];
              histoAcceptanceOmega[i]->SetName(Form("Acceptance_%s",  cutNumber[i].Data()));
            }else histoEffTimesAccOmega[i]->Multiply(histoAcceptanceOmega[i]);
            histoEffTimesAccOmega[i]->Scale(deltaRapid[i]*2*TMath::Pi());
        }

        //Scale spectrum to MBOR
        if (optionEnergy.CompareTo("2.76TeV")==0 &&
            (triggerName[i].Contains("INT7")|| triggerName[i].Contains("EMC7") || triggerName[i].Contains("EG1") || triggerName[i].Contains("EG2")) &&
            isMC.CompareTo("data") == 0){
            histoCorrectedYieldOmega[i]->Scale(0.8613) ;
            histoRawYieldOmega[i]->Scale(0.8613) ;
        }
        if (triggerName[i].CompareTo("INT7") != 0 && triggerName[i].CompareTo("MB") != 0 && triggerName[i].CompareTo("INT1") != 0){
            nRealTriggers++;
        }

        histoMCInputOmega[i]                                  = (TH1D*)fileCorrectedOmega[i]->Get(nameMCYield.Data());
        histoMCInputOmega[i]->SetName(Form("Omega_Input_Reweighted_%s",cutNumber[i].Data()));
        histoMCInputOmega[i]->Sumw2();
        histoMCInputOmega[i]->Rebin(4);
        histoMCInputOmega[i]->Scale(1./4);
        //***************************************************************************************************************
        //****************************** Calculate trigger rejection factors ********************************************
        //***************************************************************************************************************
        if ((mode == 0 || mode == 2 || mode == 3 || mode == 4 || mode == 5 || mode == 10) && hasClusterOutput ){
            histoRawClusterPt[i]                        = (TH1D*)fileUnCorrectedOmega[i]->Get("ClusterPtPerEvent");
            histoRawClusterE[i]                         = (TH1D*)fileUnCorrectedOmega[i]->Get("ClusterEPerEvent");
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
                    Double_t diffToFit                      = TMath::Abs(histoTriggerRejection[i]->GetBinContent(k+1)-triggRejecFac[i][trigSteps[i][0]])+histoTriggerRejection[i]->GetBinError(k+1);
                    if (diffToFit > largestDev) largestDev  = diffToFit;
                }
                triggRejecFacErr[i][trigSteps[i][0]] = largestDev;


                delete pol0;
                delete pol0_2;
                cout << "trigger rejection factor " << triggerName[i].Data() << "/" << triggerName[trigSteps[i][0]].Data() << ": " << triggRejecFac[i][trigSteps[i][0]]
                    << "+-" << triggRejecFacErr[i][trigSteps[i][0]] << endl;

                if (enableTriggerRejecCompMC){
                    histoMCRawClusterPt[i]                  = (TH1D*)fileUnCorrectedMCOmega[i]->Get("ClusterPtPerEvent");
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
    if (optionEnergy.CompareTo("5TeV") == 0 || optionEnergy.CompareTo("8TeV") == 0 ){
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
        else if (mode == 4 && optionEnergy.CompareTo("pPb_5.023TeV") == 0)
            maxTriggReject = 200;

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
        if (mode == 4 && optionEnergy.CompareTo("pPb_5.023TeV") == 0)
            maxTriggRejectLin = 100;
        if (mode == 4 && optionEnergy.CompareTo("2.76TeV") == 0)
            maxTriggRejectLin = 2000;

        if( optionEnergy.CompareTo("8TeV")==0 ){
        if (mode == 2 || mode == 4 || mode == 10 )
            maxTriggRejectLin = 310;
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
        //************************************Plotting unscaled invariant raw-yield Omega *********************************
        //***************************************************************************************************************

        TCanvas* canvasClusterYield = new TCanvas("canvasClusterYield","",0,0,1000,1350);// gives the page size
        DrawGammaCanvasSettings( canvasClusterYield, 0.16, 0.02, 0.015, 0.07);
        canvasClusterYield->SetLogy();

        Double_t minClusYieldUnscaled    = 7e-9;
        Double_t maxClusYieldUnscaled    = 5;

        if (mode == 4) {
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
        if (!enableEta) delete canvasClusterYield;
    }

    //***************************************************************************************************************
    //************************************Plotting efficiencies Omega *************************************************
    //***************************************************************************************************************
    TCanvas* canvasEffi = new TCanvas("canvasEffi","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasEffi, 0.09, 0.017, 0.037, 0.08);
    canvasEffi->SetLogy(0);

    Double_t minEffiOmega = 1e-4;
    Double_t maxEffiOmega = 1e-1;
    if (mode == 4){
        maxEffiOmega      = 8e-1;
    } else if (mode == 2){
        if(optionEnergy.CompareTo("8TeV")==0)
            minEffiOmega  = 5e-5;
    } else if (mode == 0){
        maxEffiOmega      = 8e-3;
        minEffiOmega      = 5e-5;
    } else if (mode == 40 || mode == 60){
        minEffiOmega = 0.8e-4;
        maxEffiOmega = 1.5e-3;
    } else{
        minEffiOmega = 1e-4;
        maxEffiOmega = 1e-2;
    }

    TH2F * histo2DEffiOmega;
    histo2DEffiOmega = new TH2F("histo2DEffiOmega","histo2DEffiOmega",1000,0., maxPtGlobalOmega*1.1,10000,minEffiOmega, maxEffiOmega);
    SetStyleHistoTH2ForGraphs(histo2DEffiOmega, "#it{p}_{T} (GeV/#it{c})","#it{#varepsilon}_{#omega}#upoint#it{#kappa}_{trigg}",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.05);
    histo2DEffiOmega->DrawCopy();
    histo2DEffiOmega->SetYTitle("#it{#varepsilon}_{#omega}");

    Double_t minXLegendEffi = 0.62;
    Int_t nColumnsEffi      = 2;
    if (nrOfTrigToBeComb > 6){
        minXLegendEffi = 0.4;
        nColumnsEffi   = 3;
    }
    Double_t minYLegendEffi = 0.13;
    Double_t maxYLegendEffi = minYLegendEffi+(1.05*nrOfTrigToBeComb/nColumnsEffi*0.85*textSizeSpectra);

    TLegend* legendEffiOmega = GetAndSetLegend2(minXLegendEffi, minYLegendEffi, 0.95, maxYLegendEffi ,28);
    legendEffiOmega->SetNColumns(nColumnsEffi);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoEfficiencyOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoEfficiencyOmega[i]->DrawCopy("e1,same");
        legendEffiOmega->AddEntry(histoEfficiencyOmega[i],triggerNameLabel[i].Data(),"p");
    }
    legendEffiOmega->Draw();

    TLatex *labelEnergyEffi = new TLatex(0.62, maxYLegendEffi+0.02+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
    SetStyleTLatex( labelEnergyEffi, 0.85*textSizeSpectra,4);
    labelEnergyEffi->Draw();

    TLatex *labelOmegaEffi = new TLatex(0.62, maxYLegendEffi+0.02+0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
    SetStyleTLatex( labelOmegaEffi, 0.85*textSizeSpectra,4);
    labelOmegaEffi->Draw();

    TLatex *labelDetProcEffi = new TLatex(0.62, maxYLegendEffi+0.02,detectionProcess.Data());
    SetStyleTLatex( labelDetProcEffi, 0.85*textSizeSpectra,4);
    labelDetProcEffi->Draw();

    canvasEffi->Update();
    canvasEffi->SaveAs(Form("%s/Omega_Efficiency_%i.%s",outputDir.Data(),mode,suffix.Data()));

    TH2F* histo2DEffiSecOmega[4]          = { NULL, NULL, NULL, NULL};
    TH2F* histo2DEffectiveSecCorr[4]    = { NULL, NULL, NULL, NULL};
    Double_t minEffiSecOmega              = minEffiOmega;
    Double_t maxEffiSecOmega              = maxEffiOmega;
    if (mode == 0){
        maxEffiSecOmega                   = maxEffiOmega/2;
        minEffiSecOmega                   = minEffiOmega/100;
    } else if (mode == 4){
        maxEffiSecOmega                   = 10*maxEffiOmega;
        minEffiSecOmega                   = 5*minEffiOmega;
    }
    for (Int_t k = 0; k< 4; k++){
        if (hasSecEffi[k]){
            canvasEffi->cd();
            histo2DEffiSecOmega[k] = new TH2F(Form("histo2DEffiSecOmega%s",nameSecOmegaPartRead[k].Data()),Form("histo2DEffiSecOmega%s",nameSecOmegaPartRead[k].Data()),
                                            1000, 0., maxPtGlobalOmega,10000,minEffiSecOmega, maxEffiSecOmega);
            SetStyleHistoTH2ForGraphs(histo2DEffiSecOmega[k], "#it{p}_{T} (GeV/#it{c})",Form("#it{#varepsilon}_{sec #omega from %s}",nameSecOmegaPartLabel[k].Data()),
                                        0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
            histo2DEffiSecOmega[k]->DrawCopy();

            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if (histoSecEffiOmegaFromX[k][i]){
                    DrawGammaSetMarker(histoSecEffiOmegaFromX[k][i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoSecEffiOmegaFromX[k][i]->DrawCopy("e1,same");
                }
            }
            legendEffiOmega->Draw();

            labelEnergyEffi->Draw();
            labelOmegaEffi->Draw();
            labelDetProcEffi->Draw();

            canvasEffi->Update();
            canvasEffi->SaveAs(Form("%s/Omega_SecEfficiencyFrom%s.%s",outputDir.Data(), nameSecOmegaPartRead[k].Data(), suffix.Data()));
        }
        if (hasSecCorrFac[k]){
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

            Double_t maxYLegendEffecSec = 0.84-0.025;
            Double_t minYLegendEffecSec = maxYLegendEffecSec-(1.05*nrOfTrigToBeComb/nColumnsEffecSec*0.85*textSizeSpectra);

            histo2DEffectiveSecCorr[k] = new TH2F(Form("histo2DEffectiveSecCorr%s",nameSecOmegaPartRead[k].Data()), Form("histo2DEffectiveSecCorr%s",nameSecOmegaPartRead[k].Data()),
                                                  1000,0., maxPtGlobalOmega,10000,0, maxYEffSecCorr[k][mode]);
            SetStyleHistoTH2ForGraphs(histo2DEffectiveSecCorr[k], "#it{p}_{T} (GeV/#it{c})",Form("r_{sec #omega from %s}", nameSecOmegaPartLabel[k].Data()),
                                        0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.19);
            histo2DEffectiveSecCorr[k]->DrawCopy();

            TLegend* legendEffectSec = GetAndSetLegend2(minXLegendEffecSec, minYLegendEffecSec, 0.95, maxYLegendEffecSec ,28, nColumnsEffecSec);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if (histoEffectCorrOmegaFromX[k][i]){
                    DrawGammaSetMarker(histoEffectCorrOmegaFromX[k][i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoEffectCorrOmegaFromX[k][i]->DrawCopy("e1,same");
                    legendEffectSec->AddEntry(histoEffectCorrOmegaFromX[k][i],triggerNameLabel[i].Data(),"p");
                }
            }
            legendEffectSec->Draw();

            TLatex *labelEnergyEffSec = new TLatex(0.62, maxYLegendEffecSec+0.02+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
            SetStyleTLatex( labelEnergyEffSec, 0.85*textSizeSpectra,4);
            labelEnergyEffSec->Draw();

            TLatex *labelOmegaEffSec = new TLatex(0.62, maxYLegendEffecSec+0.02+0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
            SetStyleTLatex( labelOmegaEffSec, 0.85*textSizeSpectra,4);
            labelOmegaEffSec->Draw();

            TLatex *labelDetProcEffSec = new TLatex(0.62, maxYLegendEffecSec+0.02,detectionProcess.Data());
            SetStyleTLatex( labelDetProcEffSec, 0.85*textSizeSpectra,4);
            labelDetProcEffSec->Draw();

            canvasEffi->Update();
            canvasEffi->SaveAs(Form("%s/Omega_EffectiveSecCorrFrom%s.%s",outputDir.Data(), nameSecOmegaPartRead[k].Data(), suffix.Data()));
            canvasEffi->SetLogy(1);
            canvasEffi->SetLeftMargin(0.09);
            canvasEffi->SetTopMargin(0.015);
            delete labelEnergyEffSec;
            delete labelOmegaEffSec;
            delete labelDetProcEffSec;
        }
    }
    //***************************************************************************************************************
    //************************************ Plotting trigger efficiencies Omega ****************************************
    //***************************************************************************************************************
    if (enableTriggerEffOmegaAll){
        TCanvas* canvasTriggerEffi = new TCanvas("canvasTriggerEffi","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasTriggerEffi, 0.09, 0.017, 0.015, 0.08);
        canvasTriggerEffi->SetLogy(0);

        Double_t minEffiTrigOmega         = 0;
        Double_t maxEffiTrigOmega         = 1.1;
        //if(optionEnergy.CompareTo("8TeV") == 0 && mode == 2) maxEffiTrigOmega = 1.5;

        TH2F * histo2DTriggerEffiOmega;
        histo2DTriggerEffiOmega = new TH2F("histo2DTriggerEffiOmega","histo2DTriggerEffiOmega",1000,0., maxPtGlobalOmega,10000,minEffiTrigOmega, maxEffiTrigOmega);
        SetStyleHistoTH2ForGraphs(histo2DTriggerEffiOmega, "#it{p}_{T} (GeV/#it{c})","#it{#kappa}_{trigg, #omega}",
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.05);
        histo2DTriggerEffiOmega->DrawCopy();

        TLegend* legendTriggerEffiOmega = GetAndSetLegend2(0.62, 0.165, 0.95, 0.165+(1.05*nRealTriggers/2*0.85*textSizeSpectra),28);
        legendTriggerEffiOmega->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if (enableTriggerEffOmega[i]){
                DrawGammaSetMarker(histoTriggerEffOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoTriggerEffOmega[i]->DrawCopy("e1,same");
                legendTriggerEffiOmega->AddEntry(histoTriggerEffOmega[i],triggerNameLabel[i].Data(),"p");
            }
        }
        legendTriggerEffiOmega->Draw();

        labelEnergyEffi->Draw();
        labelOmegaEffi->Draw();
        labelDetProcEffi->Draw();

        canvasTriggerEffi->Update();
        canvasTriggerEffi->SaveAs(Form("%s/Omega_TriggerEfficiency.%s",outputDir.Data(),suffix.Data()));
        delete canvasTriggerEffi;

        //***************************************************************************************************************
        //******************* Plotting efficiencies Omega without trigger efficiency folded in  ***************************
        //***************************************************************************************************************
        canvasEffi->cd();
        histo2DEffiOmega->DrawCopy();

        TGraphErrors* graphEffBaseOmega[12]   = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
        TLegend* legendEffiOmegaW0TriggEff    = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
        legendEffiOmegaW0TriggEff->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if ( triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1") ){
                DrawGammaSetMarker(histoEfficiencyOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoEfficiencyOmega[i]->DrawCopy("e1,same");
                legendEffiOmegaW0TriggEff->AddEntry(histoEfficiencyOmega[i],triggerNameLabel[i].Data(),"p");
            } else {
                if (enableTriggerEffOmega[i]){
                    graphEffBaseOmega[i]      = new TGraphErrors(histoEffBaseOmega[i]);
                    cout << "trying to plot trigg eff" << endl;
                    DrawGammaSetMarkerTGraph(graphEffBaseOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    graphEffBaseOmega[i]->Draw("p,e1,same");
                    legendEffiOmegaW0TriggEff->AddEntry(graphEffBaseOmega[i],triggerNameLabel[i].Data(),"p");
                }
            }
        }
        legendEffiOmegaW0TriggEff->Draw();

        labelEnergyEffi->Draw();
        labelOmegaEffi->Draw();
        labelDetProcEffi->Draw();

        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/Omega_EfficiencyW0TriggEff_%i.%s",outputDir.Data(),mode,suffix.Data()));
    }

//     if (!enableEta) delete canvasEffi;

    //***************************************************************************************************************
    //************************************Plotting acceptance Omega *************************************************
    //***************************************************************************************************************
    TCanvas* canvasAcc = new TCanvas("canvasAcc","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasAcc, 0.1, 0.017, 0.015, 0.08);
    canvasAcc->SetLogy(0);

    Double_t minAccOmega = 0.15;
    Double_t maxAccOmega = 0.3;
    if (mode == 0){
        maxAccOmega       = 1.05;
        minAccOmega       = 0.7;
    }else if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode ==45 ||
             mode == 60 || mode == 61 || mode == 62 || mode == 64 || mode ==65){
      maxAccOmega= 1.0;
      minAccOmega= 0.5;
    }

    if(optionEnergy.CompareTo("8TeV")==0){
      if (mode == 2){
        minAccOmega = 0.18;
      }else if(mode == 4){
        maxAccOmega = 0.26;
      }
    }

    TH2F * histo2DAccOmega;
    histo2DAccOmega = new TH2F("histo2DAccOmega","histo2DAccOmega",1000,0., maxPtGlobalOmega*1.1,10000,minAccOmega, maxAccOmega);
    SetStyleHistoTH2ForGraphs(histo2DAccOmega, "#it{p}_{T} (GeV/#it{c})","A_{#omega}",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.25);
    histo2DAccOmega->DrawCopy();

    Double_t minXLegendAcc = 0.62;
    Int_t nColumnsAcc      = 2;
    if (nrOfTrigToBeComb > 6){
        minXLegendAcc = 0.4;
        nColumnsAcc   = 3;
    }
    Double_t minYLegendAcc = 0.13;
    Double_t maxYLegendAcc = minYLegendAcc+(1.05*nrOfTrigToBeComb/nColumnsAcc*0.85*textSizeSpectra);

    TLegend* legendAccOmega = GetAndSetLegend2(minXLegendAcc, minYLegendAcc, 0.95,maxYLegendAcc,28);
    legendAccOmega->SetNColumns(nColumnsAcc);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoAcceptanceOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoAcceptanceOmega[i]->DrawCopy("e1,same");
        legendAccOmega->AddEntry(histoAcceptanceOmega[i],triggerNameLabel[i].Data(),"p");
    }
    legendAccOmega->Draw();

    labelEnergyEffi->Draw();
    labelOmegaEffi->Draw();
    labelDetProcEffi->Draw();


    canvasAcc->Update();
    canvasAcc->SaveAs(Form("%s/Omega_Acceptance_%i.%s",outputDir.Data(),mode,suffix.Data()));


    //***************************************************************************************************************
    //**************************** Mass and Width general plotting definitions **************************************
    //***************************************************************************************************************
    TCanvas* canvasMass         = new TCanvas("canvasMass","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasMass, 0.11, 0.017, 0.015, 0.08);
    Double_t minMassOmega         = 0.120;
    Double_t maxMassOmega         = 0.160;
    if(optionEnergy.CompareTo("8TeV")==0){
      if(mode == 2){
        minMassOmega              = 0.123;
        maxMassOmega              = 0.150;
      }else if(mode == 4){
        maxMassOmega              = 0.180;
      }
    }

    if (mode == 0){
        maxMassOmega              = 0.142;
        minMassOmega              = 0.130;
    }else if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45 ||
             mode == 60 || mode == 61 || mode == 62 || mode == 64 || mode == 65){
        maxMassOmega              = 0.84;
        minMassOmega              = 0.74;
    }

    TH2F * histo2DMassOmega       = new TH2F("histo2DMassOmega","histo2DMassOmega",1000,0., maxPtGlobalOmega,10000,minMassOmega, maxMassOmega);
    SetStyleHistoTH2ForGraphs(histo2DMassOmega, "#it{p}_{T} (GeV/#it{c})","#it{M}_{#omega} (Mev/#it{c}^{2})",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.4);
    TLegend* legendMassOmega      = GetAndSetLegend2(0.52, 0.88, 0.95, 0.88+(1.05*4/2*0.85*textSizeSpectra),28);
    TLegend* legendMassRedOmega   = GetAndSetLegend2(0.52, 0.88, 0.95, 0.88+(1.05*4/2*0.85*textSizeSpectra),28);
    TLatex *labelEnergyMass     = new TLatex(0.14, 0.85+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
    TLatex *labelOmegaMass        = new TLatex(0.14, 0.85+0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
    TLatex *labelDetProcMass    = new TLatex(0.14, 0.85,detectionProcess.Data());
    TLegend* legendMassOmega2     = GetAndSetLegend2(0.46, 0.81, 0.80, 0.81+(1.05*8/2*0.85*textSizeSpectra),28);
    TLegend* legendMassRedOmega2  = GetAndSetLegend2(0.46, 0.81, 0.80, 0.81+(1.05*8/2*0.85*textSizeSpectra),28);

    TCanvas* canvasWidth        = new TCanvas("canvasWidth","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasWidth, 0.09, 0.017, 0.035, 0.08);
    Double_t minWidthOmega        = 0.0;
    Double_t maxWidthOmega        = 0.0295;
    if (mode == 0){
        minWidthOmega             = 0.0;
        maxWidthOmega             = 0.0125;
    }
    if (mode == 4 && optionEnergy.CompareTo("8TeV")==0){
      maxWidthOmega        = 0.0395;
    }else if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45 ||
             mode == 60 || mode == 61 || mode == 62 || mode == 64 || mode == 65){
      minWidthOmega             = 0.0;
      maxWidthOmega             = 0.02;
    }
    TH2F * histo2DWidthOmega      = new TH2F("histo2DWidthOmega","histo2DWidthOmega",1000,0., maxPtGlobalOmega,10000,minWidthOmega, maxWidthOmega);
    SetStyleHistoTH2ForGraphs(histo2DWidthOmega, "#it{p}_{T} (GeV/#it{c})","#sigma_{#omega} (Mev/#it{c}^{2})",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
    TLegend* legendWidthOmega     = GetAndSetLegend2(0.52, 0.87, 0.95, 0.87+(1.05*4/2*0.85*textSizeSpectra),28);
    TLegend* legendWidthRedOmega  = GetAndSetLegend2(0.52, 0.87, 0.95, 0.87+(1.05*4/2*0.85*textSizeSpectra),28);
    TLatex* labelEnergyWidth    = new TLatex(0.14, 0.84+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
    TLatex* labelOmegaWidth       = new TLatex(0.14, 0.84+0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
    TLatex* labelDetProcWidth   = new TLatex(0.14, 0.84,detectionProcess.Data());
    TLegend* legendWidthOmega2    = GetAndSetLegend2(0.46, 0.80, 0.80, 0.80+(1.05*8/2*0.85*textSizeSpectra),28);
    TLegend* legendWidthRedOmega2 = GetAndSetLegend2(0.46, 0.80, 0.80, 0.80+(1.05*8/2*0.85*textSizeSpectra),28);

    if (mode != 10){
        //***************************************************************************************************************
        //************************************Plotting Mass Omega *********************************************************
        //***************************************************************************************************************
        canvasMass->cd();
        histo2DMassOmega->DrawCopy();

        legendMassOmega->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
               (optionEnergy.CompareTo("8TeV")==0) ||
               (optionEnergy.CompareTo("pPb_5.023TeV")==0)
               ){
                DrawGammaSetMarker(histoMassOmegaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoMassOmegaData[i]->DrawCopy("e1,same");
                legendMassOmega->AddEntry(histoMassOmegaData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
                DrawGammaSetMarker(histoMassOmegaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                histoMassOmegaMC[i]->DrawCopy("e1,same");
                legendMassOmega->AddEntry(histoMassOmegaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
            }
        }
        legendMassOmega->Draw();

        SetStyleTLatex( labelEnergyMass, 0.85*textSizeSpectra,4);
        labelEnergyMass->Draw();


        SetStyleTLatex( labelOmegaMass, 0.85*textSizeSpectra,4);
        labelOmegaMass->Draw();

        SetStyleTLatex( labelDetProcMass, 0.85*textSizeSpectra,4);
        labelDetProcMass->Draw();

        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/Omega_%s_Mass1.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        histo2DMassOmega->DrawCopy();


        if (optionEnergy.CompareTo("2.76TeV")==0){
            legendMassOmega2->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if( optionEnergy.CompareTo("2.76TeV")==0 && (i==1 || i==3 || i==4 || i==5)){
                    DrawGammaSetMarker(histoMassOmegaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoMassOmegaData[i]->DrawCopy("e1,same");
                    legendMassOmega2->AddEntry(histoMassOmegaData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
                    DrawGammaSetMarker(histoMassOmegaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    histoMassOmegaMC[i]->DrawCopy("e1,same");
                    legendMassOmega2->AddEntry(histoMassOmegaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
                }
            }
            legendMassOmega2->Draw();
            labelEnergyMass->Draw();
            labelOmegaMass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Omega_%s_Mass2.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }
        //***************************************************************************************************************
        //************************************Plotting Width Omega *********************************************************
        //***************************************************************************************************************
        canvasWidth->cd();
        histo2DWidthOmega->DrawCopy();

        legendWidthOmega->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0) ||
             (optionEnergy.CompareTo("pPb_5.023TeV")==0)
             ){
                DrawGammaSetMarker(histoWidthOmegaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                histoWidthOmegaData[i]->DrawCopy("e1,same");
                legendWidthOmega->AddEntry(histoWidthOmegaData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
                DrawGammaSetMarker(histoWidthOmegaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                histoWidthOmegaMC[i]->DrawCopy("e1,same");
                legendWidthOmega->AddEntry(histoWidthOmegaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
            }
        }
        legendWidthOmega->Draw();

        SetStyleTLatex( labelEnergyWidth, 0.85*textSizeSpectra,4);
        labelEnergyWidth->Draw();

        SetStyleTLatex( labelOmegaWidth, 0.85*textSizeSpectra,4);
        labelOmegaWidth->Draw();

        SetStyleTLatex( labelDetProcWidth, 0.85*textSizeSpectra,4);
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Omega_%s_Width1_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

        histo2DWidthOmega->DrawCopy();

        if (optionEnergy.CompareTo("2.76TeV")==0) {
            legendWidthOmega2->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if(i==1 || i==3 || i==4 || i==5 ){
                    DrawGammaSetMarker(histoWidthOmegaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    histoWidthOmegaData[i]->DrawCopy("e1,same");
                    legendWidthOmega2->AddEntry(histoWidthOmegaData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
                    DrawGammaSetMarker(histoWidthOmegaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    histoWidthOmegaMC[i]->DrawCopy("e1,same");
                    legendWidthOmega2->AddEntry(histoWidthOmegaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
                }
            }
            legendWidthOmega2->Draw();
            labelEnergyWidth->Draw();
            labelOmegaWidth->Draw();
            labelDetProcWidth->Draw();

            canvasWidth->Update();
            canvasWidth->SaveAs(Form("%s/Omega_%s_Width2_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
        }
    }

    //***************************************************************************************************************
    //************************************Plotting unscaled invariant raw-yield Omega *********************************
    //***************************************************************************************************************
    textSizePixelSpectra = textSizeSpectra*1000;

    TCanvas* canvasRawUnscaled = new TCanvas("canvasRawUnscaled","",0,0,1000,1350);// gives the page size
    DrawGammaCanvasSettings( canvasRawUnscaled, 0.16, 0.02, 0.015, 0.07);
    canvasRawUnscaled->SetLogy();

    Double_t minCorrYieldRawUnscaled    = 7e-8;
    Double_t maxCorrYieldRawUnscaled    = 4e-2;
    if (mode == 4) {
      minCorrYieldRawUnscaled         = 7e-8;
      maxCorrYieldRawUnscaled         = 1;
    } else if (mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45 ||
               mode == 60 || mode == 61 || mode == 62 || mode == 64 || mode == 65){
      minCorrYieldRawUnscaled         = 2e-8;
      maxCorrYieldRawUnscaled         = 4e-3;
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

    if(optionEnergy.CompareTo("pPb_5.023TeV")==0){
      if(mode == 2){
        minCorrYieldRawUnscaled         = 1e-8;
        maxCorrYieldRawUnscaled         = 1e-4;
      }
    }

    TH2F * histo2DRawUnscaled       = new TH2F("histo2DRawUnscaled", "histo2DRawUnscaled", 1000, 0., maxPtGlobalOmega, 10000, minCorrYieldRawUnscaled, maxCorrYieldRawUnscaled);
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
        DrawGammaSetMarker(histoRawYieldOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoRawYieldOmega[i]->DrawCopy("e1,same");
        legendRawUnscaled->AddEntry(histoRawYieldOmega[i],triggerNameLabel[i].Data(),"p");
    }
    legendRawUnscaled->Draw();


    TLatex *labelEnergyRawUnscaled  = new TLatex(0.2, 0.12+(2*textSizeSpectra*0.85*0.75),collisionSystem.Data());
    SetStyleTLatex( labelEnergyRawUnscaled, 0.85*textSizeSpectra,4);
    labelEnergyRawUnscaled->Draw();

    TLatex *labelOmegaRawUnscaled     = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
    SetStyleTLatex( labelOmegaRawUnscaled, 0.85*textSizeSpectra,4);
    labelOmegaRawUnscaled->Draw();

    TLatex *labelDetProcRawUnscaled = new TLatex(0.2, 0.12,detectionProcess.Data());
    SetStyleTLatex( labelDetProcRawUnscaled, 0.85*textSizeSpectra,4);
    labelDetProcRawUnscaled->Draw();

    canvasRawUnscaled->Update();
    canvasRawUnscaled->SaveAs(Form("%s/Omega_%s_RawYieldUnscaledTrigg_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
    if (!enableEta) delete canvasRawUnscaled;

    //***************************************************************************************************************
    //************************************Plotting unscaled invariant yield Omega *************************************
    //***************************************************************************************************************
    TCanvas* canvasCorrUnscaled = new TCanvas("canvasCorrUnscaled","",0,0,1000,1350);// gives the page size
    DrawGammaCanvasSettings( canvasCorrUnscaled, 0.15, 0.017, 0.015, 0.07);
    canvasCorrUnscaled->SetLogy();

    Double_t minCorrYieldUnscaled   = 2e-10;
    Double_t maxCorrYieldUnscaled   = 1e2;
    if (mode == 0){
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

    if(optionEnergy.CompareTo("pPb_5.023TeV")==0){
      if(mode == 2){
        minCorrYieldUnscaled        = 1e-8;
        maxCorrYieldUnscaled        = 1;
      }else if(mode == 4){
        minCorrYieldUnscaled        = 2e-8;
        maxCorrYieldUnscaled        = 0.2;
      }
    }

    TH2F * histo2DInvYieldUnscaled = new TH2F("histo2DInvYieldUnscaled","histo2DInvYieldUnscaled",1000,0., maxPtGlobalOmega,10000,minCorrYieldUnscaled,maxCorrYieldUnscaled);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldUnscaled, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                            0.85*textSizeSpectra,0.04, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.55);
    histo2DInvYieldUnscaled->DrawCopy();

    TLegend* legendUnscaled = GetAndSetLegend2(minXLegendRaw, minYLegendRaw, 0.95, maxYLegendRaw,0.85*textSizePixelSpectra);
    legendUnscaled->SetNColumns(nColumnsRaw);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        DrawGammaSetMarker(histoCorrectedYieldOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        histoCorrectedYieldOmega[i]->DrawCopy("e1,same");
        legendUnscaled->AddEntry(histoCorrectedYieldOmega[i],triggerNameLabel[i].Data(),"p");
    }
    legendUnscaled->Draw();


    TLatex *labelEnergyUnscaled = new TLatex(0.2, 0.12+(2*textSizeSpectra*0.85*0.75),collisionSystem.Data());
    SetStyleTLatex( labelEnergyUnscaled, 0.85*textSizeSpectra,4);
    labelEnergyUnscaled->Draw();

    TLatex *labelOmegaUnscaled = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
    SetStyleTLatex( labelOmegaUnscaled, 0.85*textSizeSpectra,4);
    labelOmegaUnscaled->Draw();

    TLatex *labelDetProcUnscaled = new TLatex(0.2, 0.12,detectionProcess.Data());
    SetStyleTLatex( labelDetProcUnscaled, 0.85*textSizeSpectra,4);
    labelDetProcUnscaled->Draw();

    canvasCorrUnscaled->Update();
    canvasCorrUnscaled->SaveAs(Form("%s/Omega_%s_CorrectedYieldUnscaledTrigg_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
    if (!enableEta) delete canvasCorrUnscaled;

    //***************************************************************************************************************
    //******************************* Scaling corrected yield by trigger rejection factors **************************
    //***************************************************************************************************************

    TH1D*     histoCorrectedYieldOmegaScaled                  [MaxNumberOfFiles];
    TH1D*     histoCorrectedYieldOmegaScaledMasked            [MaxNumberOfFiles];
    histoCorrectedYieldOmegaScaled[0]                         = (TH1D*)histoCorrectedYieldOmega[0]->Clone(Form("CorrectedYieldOmegaScaled_%s", triggerName[0].Data()));
    for (Int_t i = 1; i< nrOfTrigToBeComb; i++){
        fileFitsOutput << triggerName[i].Data() << endl;
        histoCorrectedYieldOmegaScaled[i] = (TH1D*)histoCorrectedYieldOmega[i]->Clone(Form("CorrectedYieldOmegaScaled_%s", triggerName[i].Data()));
        histoCorrectedYieldOmegaScaled[i]->Sumw2();
        histoCorrectedYieldOmegaScaled[i]->Scale(1./triggRejecFac[i][trigSteps[i][0]]);
        fileFitsOutput << trigSteps[i][0] << "\t" << trigSteps[i][1] << "\t" << trigSteps[i][2] << endl;
        if (trigSteps[i][1]!= trigSteps[i][0]){
            fileFitsOutput << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << endl;
            histoCorrectedYieldOmegaScaled[i]->Scale(1./triggRejecFac[trigSteps[i][0]][trigSteps[i][1]]);
        }
        if (trigSteps[i][2]!= trigSteps[i][1]){
            fileFitsOutput << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << "\t"<< triggRejecFac[trigSteps[i][1]][trigSteps[i][2]] << endl;
            histoCorrectedYieldOmegaScaled[i]->Scale(1./triggRejecFac[trigSteps[i][1]][trigSteps[i][2]]);
        }
    }

    // initialize all vectors for sytstematics and general creation of graphs, we can have a maximu of 100 data points at the moment
    Double_t xValueFinalOmega                                 [100];
    Double_t xErrorLowFinalOmega                              [100];
    Double_t xErrorHighFinalOmega                             [100];
    Double_t yValueFinalOmega                                 [100];
    Double_t yErrorLowFinalOmega                              [100];
    Double_t yErrorHighFinalOmega                             [100];
    Int_t nPointFinalOmega                                         = 0;

    Double_t yErrorSysLowFinalOmega                           [100];
    Double_t yErrorSysHighFinalOmega                          [100];

    Double_t ptSysRelOmega                                    [MaxNumberOfFiles][100];
    Double_t yErrorSysLowRelOmega                             [MaxNumberOfFiles][100];
    Double_t yErrorSysHighRelOmega                            [MaxNumberOfFiles][100];
    Bool_t sysAvailOmega                                      [MaxNumberOfFiles];

    Bool_t sysAvailSingleOmega                                [MaxNumberOfFiles];
    Int_t numberBinsSysAvailSingleOmega                       [MaxNumberOfFiles];

    // graphs for easier reduction of measurements to desired range and systematic errors
    TGraphAsymmErrors* graphsCorrectedYieldShrunkOmega        [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldSysShrunkOmega     [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldRemoved0Omega      [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldSysRemoved0Omega   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphMassOmegaData                     [MaxNumberOfFiles];
    TGraphAsymmErrors* graphMassOmegaMC                       [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedMassOmegaData              [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedMassOmegaMC                [MaxNumberOfFiles];
    TGraphAsymmErrors* graphWidthOmegaData                    [MaxNumberOfFiles];
    TGraphAsymmErrors* graphWidthOmegaMC                      [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedWidthOmegaData             [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedWidthOmegaMC               [MaxNumberOfFiles];
    TGraphAsymmErrors* graphAcceptanceOmega                   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedAcceptanceOmega            [MaxNumberOfFiles];
    TGraphAsymmErrors* graphEfficiencyOmega                   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEfficiencyOmega            [MaxNumberOfFiles];
    TGraphAsymmErrors* graphEffTimesAccOmega                  [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEffTimesAccOmega           [MaxNumberOfFiles];
    TGraphAsymmErrors* graphPurityOmega                       [MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedPurityOmega                [MaxNumberOfFiles];
    TGraphAsymmErrors* graphEffectSecCorrOmega                [4][MaxNumberOfFiles];
    TGraphAsymmErrors* graphEfficiencySecOmega                [4][MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEffectSecCorrOmega         [4][MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedEfficiencySecOmega         [4][MaxNumberOfFiles];

    Int_t nRelSysErrOmegaSources          = 0;
    TGraphAsymmErrors* graphRelSysErrOmegaSource              [30][MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedRelSysErrOmegaSource       [30][MaxNumberOfFiles];
    for (Int_t j = 0; j< 30; j++){
        for (Int_t k = 0; k< MaxNumberOfFiles; k++){
            graphRelSysErrOmegaSource[j][k]           = NULL;
            graphOrderedRelSysErrOmegaSource[j][k]    = NULL;
        }
    }
    // definition of predefined arrays for trigger correlation filling
    TH1D*               histoStatOmega    [12];
    TGraphAsymmErrors*  graphSystOmega    [12];
    TH1D*               histoRelStatOmega [12];
    TGraphAsymmErrors*  graphRelSystOmega [12];

    Int_t offSetsOmega[12]        =   { 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0 };
    Int_t offSetsOmegaSys[12]     =   { 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0 };

    if(optionEnergy.CompareTo("8TeV")==0){
      if(mode == 2){
        offSetsOmega[1] = 3; //INT7
        offSetsOmega[3] = 0; //EMC7
        offSetsOmega[4] = 3; //EGA
      }else if(mode == 4){
        offSetsOmega[1] = 0; //INT7
        offSetsOmega[3] = 0; //EMC7
        offSetsOmega[4] = 4; //EGA
      }
    }

    // set all graphs to NULL first
    for (Int_t j = 0; j<12; j++){
        histoStatOmega[j]                 = NULL;
        graphSystOmega[j]                 = NULL;
        histoRelStatOmega[j]              = NULL;
        graphRelSystOmega[j]              = NULL;
        graphOrderedMassOmegaData[j]      = NULL;
        graphOrderedMassOmegaMC[j]        = NULL;
        graphOrderedWidthOmegaData[j]     = NULL;
        graphOrderedWidthOmegaMC[j]       = NULL;
        graphOrderedAcceptanceOmega[j]    = NULL;
        graphOrderedEfficiencyOmega[j]    = NULL;
        graphOrderedEffTimesAccOmega[j]   = NULL;
        graphOrderedPurityOmega[j]        = NULL;
        for (Int_t k = 0; k< 4; k++){
            graphOrderedEfficiencySecOmega[k][j]  = NULL;
            graphOrderedEffectSecCorrOmega[k][j]  = NULL;
        }
    }


    //****************************************************************************************************************
    //************* Processing of each individual trigger, reducing ranges & adding systematics **********************
    //****************************************************************************************************************
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        // read systematics, if fileName is set to "bla" no action has to be performed and systematics will be disabled for the rest of the analysis
        if (sysFileOmega[i].CompareTo("bla") != 0){
            sysAvailOmega[i]                = kTRUE;
            ifstream  fileSysErrOmega;
            fileSysErrOmega.open(sysFileOmega[i].Data(),ios_base::in);
            cout << sysFileOmega[i].Data() << endl;
            gSystem->Exec(Form("cp %s %s/SystematicErrorAveraged_Omega_%s.txt", sysFileOmega[i].Data(), outputDir.Data(),triggerName[i].Data()));


            Int_t iPtBin = 0;
            cout << "reading sys file summed" << endl;
            while(!fileSysErrOmega.eof() && iPtBin < 100){
                Double_t garbage = 0;
                fileSysErrOmega >>ptSysRelOmega[i][iPtBin] >> yErrorSysLowRelOmega[i][iPtBin] >> yErrorSysHighRelOmega[i][iPtBin]>>    garbage >> garbage;
                cout << iPtBin << "\t"<< ptSysRelOmega[i][iPtBin]<< "\t"  << yErrorSysLowRelOmega[i][iPtBin] << "\t"  <<yErrorSysHighRelOmega[i][iPtBin] << "\t"  << endl;;
                iPtBin++;
            }
            fileSysErrOmega.close();
         // read in detailed systematics
            string sysFileOmegaDet = sysFileOmega[i].Data();

            if(!replace(sysFileOmegaDet, "Averaged", "AveragedSingle")){
                cout << "WARNING: could not find detailed systematics file " << sysFileOmegaDet << ", skipping... " << endl;
                sysAvailSingleOmega[i] = kFALSE;
            }else{
                ifstream fileSysErrDetailedOmega;
                fileSysErrDetailedOmega.open(sysFileOmegaDet,ios_base::in);
                if(fileSysErrDetailedOmega.is_open()) {
                    sysAvailSingleOmega[i] = kTRUE;
                    gSystem->Exec(Form("cp %s %s/SystematicErrorAveragedSingle_Omega_%s.txt", ((TString)sysFileOmegaDet).Data(), outputDir.Data(),triggerName[i].Data()));
                } else{
                    sysAvailSingleOmega[i] = kFALSE;
                    cout << "No single errors were found" << endl;
                }

                if (sysAvailSingleOmega[i]){
                    cout << sysFileOmegaDet << endl;
                    iPtBin = 0;
                    string line;
                    Int_t iPtBinColumn = 0;
                    while (getline(fileSysErrDetailedOmega, line) && iPtBin < 100) {
                        istringstream ss(line);
                        TString temp="";
                        iPtBinColumn = 0;
                        while(ss && iPtBinColumn < 100){
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
                    numberBinsSysAvailSingleOmega[i] = iPtBin;
                    fileSysErrDetailedOmega.close();
                 }
             }
        } else {
            sysAvailOmega[i]             = kFALSE;
            sysAvailSingleOmega[i]       = kFALSE;
        }
        cout << sysAvailOmega[i] << "\t" << sysAvailSingleOmega[i] << endl;
//         continue;


        // print out input spectrum from statistical histogram
        cout << "step 0" << endl;
        for (Int_t j = 1; j< histoCorrectedYieldOmegaScaled[i]->GetNbinsX()+1; j++ ){
            cout << histoCorrectedYieldOmegaScaled[i]->GetBinCenter(j) << "\t" << histoCorrectedYieldOmegaScaled[i]->GetBinContent(j) << endl;
        }
        //*******************************************************************
        //**** create graphs from original histograms for every quantity ****
        //*******************************************************************
        cout << "step 1" << endl;
        // create correct yield graphs
        graphsCorrectedYieldShrunkOmega[i]        = new TGraphAsymmErrors(histoCorrectedYieldOmegaScaled[i]);
        graphsCorrectedYieldRemoved0Omega[i]      = new TGraphAsymmErrors(histoCorrectedYieldOmegaScaled[i]);
        graphsCorrectedYieldSysShrunkOmega[i]     = new TGraphAsymmErrors(histoCorrectedYieldOmegaScaled[i]);
        graphsCorrectedYieldSysRemoved0Omega[i]   = new TGraphAsymmErrors(histoCorrectedYieldOmegaScaled[i]);
        histoCorrectedYieldOmegaScaledMasked[i]   = (TH1D*)histoCorrectedYieldOmegaScaled[i]->Clone(Form("Omega_ScaledMasked_%s",triggerName[i].Data()));

        // create supporting figure graphs
        if (mode != 10){ // these are not available for merged cluster analysis
            graphMassOmegaData[i]                 = new TGraphAsymmErrors(histoMassOmegaData[i]);
            graphMassOmegaMC[i]                   = new TGraphAsymmErrors(histoMassOmegaMC[i]);
            graphWidthOmegaData[i]                = new TGraphAsymmErrors(histoWidthOmegaData[i]);
            graphWidthOmegaMC[i]                  = new TGraphAsymmErrors(histoWidthOmegaMC[i]);
        }

        graphAcceptanceOmega[i]                   = new TGraphAsymmErrors(histoAcceptanceOmega[i]);
        for (Int_t k = 0; k< 4; k++){
            if (hasSecCorrFac[k]){
                graphEffectSecCorrOmega[k][i]     = new TGraphAsymmErrors(histoEffectCorrOmegaFromX[k][i]);
            } else {
                graphEffectSecCorrOmega[k][i]     = NULL;
            }
            if (hasSecEffi[k]){
                graphEfficiencySecOmega[k][i]     = new TGraphAsymmErrors(histoSecEffiOmegaFromX[k][i]);
            } else {
                graphEfficiencySecOmega[k][i]     = NULL;
            }
        }

        if ( triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1")){
            graphEfficiencyOmega[i]               = new TGraphAsymmErrors(histoEfficiencyOmega[i]);

        } else { // only possible if same file had been processed with pure MB cut on the same MC
            if (enableTriggerEffOmega[i]){ // only possible if same file had been processed with pure MB cut on the same MC
                graphEfficiencyOmega[i]           = new TGraphAsymmErrors(histoEffBaseOmega[i]);
            } else {
                graphEfficiencyOmega[i]           = NULL;
            }
        }
        graphEffTimesAccOmega[i]                  = new TGraphAsymmErrors(histoEffTimesAccOmega[i]);

        // remove 0 bins at beginning
        Int_t binsToMask = 1;
        while (histoCorrectedYieldOmegaScaledMasked[i]->GetBinCenter(binsToMask) < ptFromSpecOmega[i][0] ){
            histoCorrectedYieldOmegaScaledMasked[i]->SetBinContent(binsToMask,0.);
            histoCorrectedYieldOmegaScaledMasked[i]->SetBinError(binsToMask,0.);
            binsToMask++;
        }
        while (graphAcceptanceOmega[i]->GetX()[0] < ptFromSpecOmega[i][0] ){
            if (mode != 10){
                graphMassOmegaData[i]->RemovePoint(0);
                graphMassOmegaMC[i]->RemovePoint(0);
                graphWidthOmegaData[i]->RemovePoint(0);
                graphWidthOmegaMC[i]->RemovePoint(0);
            }

            graphAcceptanceOmega[i]->RemovePoint(0);
            graphEffTimesAccOmega[i]->RemovePoint(0);
            for (Int_t k = 0; k< 4; k++){
                if (graphEffectSecCorrOmega[k][i])graphEffectSecCorrOmega[k][i]->RemovePoint(0);
                if (graphEfficiencySecOmega[k][i])graphEfficiencySecOmega[k][i]->RemovePoint(0);
            }
            if (enableTriggerEffOmega[i] || (triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1"))){
                graphEfficiencyOmega[i]->RemovePoint(0);
            }
        }

        // check if trigger is supposed to be used for combination, otherwise put all graphs to NULL
        if ( maskedFullyOmega[i] ){
            graphMassOmegaData[i]     = NULL;
            graphMassOmegaMC[i]       = NULL;
            graphWidthOmegaData[i]    = NULL;
            graphWidthOmegaMC[i]      = NULL;
            graphAcceptanceOmega[i]   = NULL;
            graphEffTimesAccOmega[i]  = NULL;
            graphPurityOmega[i]       = NULL;
            graphEfficiencyOmega[i]   = NULL;
            graphsCorrectedYieldShrunkOmega[i]    = NULL;
            graphsCorrectedYieldRemoved0Omega[i]  = NULL;
            graphsCorrectedYieldSysShrunkOmega[i] = NULL;
            graphsCorrectedYieldSysRemoved0Omega[i]   = NULL;
            for (Int_t k = 0; k< 4; k++){
                graphEffectSecCorrOmega[k][i]         = NULL;
                graphEfficiencySecOmega[k][i]         = NULL;
            }
            sysAvailOmega[i]          = kFALSE;
            for (Int_t f = 1; f < histoCorrectedYieldOmegaScaledMasked[i]->GetNbinsX()+1; f++ ){
                histoCorrectedYieldOmegaScaledMasked[i]->SetBinContent(f,0.);
                histoCorrectedYieldOmegaScaledMasked[i]->SetBinError(f,0.);
            }
            nrOfTrigToBeCombOmegaRed--;
            cout << "trigger " << triggerName[i] << " was masked" << endl;
            continue;
        // if upper boundary > -1, remove all points above
        } else if (ptFromSpecOmega[i][1] > -1) {
            for (Int_t f = histoCorrectedYieldOmegaScaledMasked[i]->GetXaxis()->FindBin(ptFromSpecOmega[i][1]); f < histoCorrectedYieldOmegaScaledMasked[i]->GetNbinsX()+1; f++ ){
                histoCorrectedYieldOmegaScaledMasked[i]->SetBinContent(f,0.);
                histoCorrectedYieldOmegaScaledMasked[i]->SetBinError(f,0.);
            }
            while (graphAcceptanceOmega[i]->GetX()[graphAcceptanceOmega[i]->GetN()-1] > ptFromSpecOmega[i][1] ){
                if (mode != 10){
                    graphMassOmegaData[i]->RemovePoint(graphMassOmegaData[i]->GetN()-1);
                    graphMassOmegaMC[i]->RemovePoint(graphMassOmegaMC[i]->GetN()-1);
                    graphWidthOmegaData[i]->RemovePoint(graphWidthOmegaData[i]->GetN()-1);
                    graphWidthOmegaMC[i]->RemovePoint(graphWidthOmegaMC[i]->GetN()-1);
                }

                graphAcceptanceOmega[i]->RemovePoint(graphAcceptanceOmega[i]->GetN()-1);
                graphEffTimesAccOmega[i]->RemovePoint(graphEffTimesAccOmega[i]->GetN()-1);
                for (Int_t k = 0; k< 4; k++){
                    if (graphEffectSecCorrOmega[k][i])graphEffectSecCorrOmega[k][i]->RemovePoint(graphEffectSecCorrOmega[k][i]->GetN()-1);
                    if (graphEfficiencySecOmega[k][i])graphEfficiencySecOmega[k][i]->RemovePoint(graphEfficiencySecOmega[k][i]->GetN()-1);
                }
                if (enableTriggerEffOmega[i] || (triggerName[i].Contains("INT7") || triggerName[i].Contains("MB") || triggerName[i].Contains("INT1"))){
                    graphEfficiencyOmega[i]->RemovePoint(graphEfficiencyOmega[i]->GetN()-1);
                }
            }
        }
        // if no points are left in graph put graph to NULL
        if (graphAcceptanceOmega[i]->GetN() == 0){
            if (mode != 10){
                graphMassOmegaData[i]     = NULL;
                graphMassOmegaMC[i]       = NULL;
                graphWidthOmegaData[i]    = NULL;
                graphWidthOmegaMC[i]      = NULL;
            }
            graphAcceptanceOmega[i]       = NULL;
            graphEffTimesAccOmega[i]      = NULL;
            graphEfficiencyOmega[i]       = NULL;
            for (Int_t k = 0; k< 4; k++){
                graphEffectSecCorrOmega[k][i] = NULL;
                graphEfficiencySecOmega[k][i] = NULL;
            }

        }

        // Remove 0 points at beginning for graphs
//         if (graphsCorrectedYieldShrunkOmega[i])graphsCorrectedYieldShrunkOmega[i]->Print();
        cout << "step 2" << endl;
        while (graphsCorrectedYieldShrunkOmega[i]->GetY()[0] == 0) graphsCorrectedYieldShrunkOmega[i]->RemovePoint(0);
//         if (graphsCorrectedYieldShrunkOmega[i]) graphsCorrectedYieldShrunkOmega[i]->Print();
        while (graphsCorrectedYieldRemoved0Omega[i]->GetY()[0] == 0) graphsCorrectedYieldRemoved0Omega[i]->RemovePoint(0);
//         if (graphsCorrectedYieldRemoved0Omega[i]) graphsCorrectedYieldRemoved0Omega[i]->Print();
        cout << "sys shrunk" << endl;
        while (graphsCorrectedYieldSysShrunkOmega[i]->GetY()[0] == 0) graphsCorrectedYieldSysShrunkOmega[i]->RemovePoint(0);
//         if (graphsCorrectedYieldSysShrunkOmega[i])graphsCorrectedYieldSysShrunkOmega[i]->Print();
        cout << "sys shrunk 2" << endl;
        while (graphsCorrectedYieldSysRemoved0Omega[i]->GetY()[0] == 0) graphsCorrectedYieldSysRemoved0Omega[i]->RemovePoint(0);
//         if (graphsCorrectedYieldSysRemoved0Omega[i])graphsCorrectedYieldSysRemoved0Omega[i]->Print();

        // put systematics on graphs
        if (graphsCorrectedYieldSysRemoved0Omega[i]){
            if (sysAvailSingleOmega[i]){
                nRelSysErrOmegaSources                    = (Int_t)ptSysDetail[i][0].size()-1;
                for (Int_t k = 0; k < nRelSysErrOmegaSources; k++ ){
                    graphRelSysErrOmegaSource[k][i]       = (TGraphAsymmErrors*) graphsCorrectedYieldSysRemoved0Omega[i]->Clone(Form("RelSysErrOmegaSource%s_%s",((TString)ptSysDetail[i][0].at(k+1)).Data(), triggerName[i].Data()));
                    cout << Form("RelSysErrSource%s_%s",((TString)ptSysDetail[i][0].at(k+1)).Data(), triggerName[i].Data()) << endl;
                }
            }
            for (Int_t j = 0; j< graphsCorrectedYieldSysRemoved0Omega[i]->GetN(); j++){
                if (sysAvailOmega[i]){
                    Int_t counter = 0;
                    while(counter < 100 && TMath::Abs(graphsCorrectedYieldSysRemoved0Omega[i]->GetX()[j] - ptSysRelOmega[i][counter])> 0.001) counter++;
                    if (counter < 100){
                        cout << ptSysRelOmega[i][counter]<< "\t found it" << endl;
                        Double_t yErrorSysLowDummy  = TMath::Abs(yErrorSysLowRelOmega[i][counter]/100*graphsCorrectedYieldSysRemoved0Omega[i]->GetY()[j]);
                        Double_t yErrorSysHighDummy = yErrorSysHighRelOmega[i][counter]/100*graphsCorrectedYieldSysRemoved0Omega[i]->GetY()[j];
                        graphsCorrectedYieldSysRemoved0Omega[i]->SetPointEYlow(j,yErrorSysLowDummy);
                        graphsCorrectedYieldSysRemoved0Omega[i]->SetPointEYhigh(j,yErrorSysHighDummy);
                        if (sysAvailSingleOmega[i]){
                            for (Int_t k = 0; k < nRelSysErrOmegaSources; k++ ){
                                graphRelSysErrOmegaSource[k][i]->SetPoint(j, graphsCorrectedYieldSysRemoved0Omega[i]->GetX()[j] ,((TString)ptSysDetail[i][counter+1].at(k+1)).Atof());
                                graphRelSysErrOmegaSource[k][i]->SetPointEYhigh(j,0);
                                graphRelSysErrOmegaSource[k][i]->SetPointEYlow(j,0);
                            }
                        }
                    } else {
                        graphsCorrectedYieldSysRemoved0Omega[i]->SetPointEYlow(j,0);
                        graphsCorrectedYieldSysRemoved0Omega[i]->SetPointEYhigh(j,0);
                        if (sysAvailSingleOmega[i]){
                            for (Int_t k = 0; k < nRelSysErrOmegaSources; k++ ){
                                graphRelSysErrOmegaSource[k][i]->SetPoint(j, graphsCorrectedYieldSysRemoved0Omega[i]->GetX()[j] ,0);
                                graphRelSysErrOmegaSource[k][i]->SetPointEYlow(j,0);
                                graphRelSysErrOmegaSource[k][i]->SetPointEYhigh(j,0);
                            }
                        }
                    }
                } else {
                    graphsCorrectedYieldSysRemoved0Omega[i]->SetPointEYlow(j,0);
                    graphsCorrectedYieldSysRemoved0Omega[i]->SetPointEYhigh(j,0);
                    averagedOmega = kFALSE;
                    if (sysAvailSingleOmega[i]){
                        for (Int_t k = 0; k < nRelSysErrOmegaSources; k++ ){
                            graphRelSysErrOmegaSource[k][i]->SetPoint(j, graphsCorrectedYieldSysRemoved0Omega[i]->GetX()[j] ,0);
                            graphRelSysErrOmegaSource[k][i]->SetPointEYlow(j,0);
                            graphRelSysErrOmegaSource[k][i]->SetPointEYhigh(j,0);
                        }
                    }
                }
            }
        }

        cout << "step 3" << endl;
        cout << "range to be accepted: " << ptFromSpecOmega[i][0] << "\t-\t" << ptFromSpecOmega[i][1] << endl;
        // remove points at beginning according to ranges set for individual triggers
        while(graphsCorrectedYieldShrunkOmega[i]->GetX()[0] < ptFromSpecOmega[i][0])
            graphsCorrectedYieldShrunkOmega[i]->RemovePoint(0);
        while(graphsCorrectedYieldSysShrunkOmega[i]->GetX()[0] < ptFromSpecOmega[i][0])
            graphsCorrectedYieldSysShrunkOmega[i]->RemovePoint(0);
        // remove points at the end according to ranges set for individual triggers
        while (graphsCorrectedYieldShrunkOmega[i]->GetX()[graphsCorrectedYieldShrunkOmega[i]->GetN()-1] > ptFromSpecOmega[i][1])
            graphsCorrectedYieldShrunkOmega[i]->RemovePoint(graphsCorrectedYieldShrunkOmega[i]->GetN()-1);
//         graphsCorrectedYieldShrunkOmega[i]->Print();
        while (graphsCorrectedYieldSysShrunkOmega[i]->GetX()[graphsCorrectedYieldSysShrunkOmega[i]->GetN()-1] > ptFromSpecOmega[i][1])
            graphsCorrectedYieldSysShrunkOmega[i]->RemovePoint(graphsCorrectedYieldSysShrunkOmega[i]->GetN()-1);

        // put systematics on shrunk graphs
        for (Int_t j = 0; j< graphsCorrectedYieldShrunkOmega[i]->GetN(); j++){
            xValueFinalOmega[nPointFinalOmega]      = graphsCorrectedYieldShrunkOmega[i]->GetX()[j];
            xErrorHighFinalOmega[nPointFinalOmega]  = graphsCorrectedYieldShrunkOmega[i]->GetEXhigh()[j];
            xErrorLowFinalOmega[nPointFinalOmega]   = graphsCorrectedYieldShrunkOmega[i]->GetEXlow()[j];
            yValueFinalOmega[nPointFinalOmega]      = graphsCorrectedYieldShrunkOmega[i]->GetY()[j];
            yErrorHighFinalOmega[nPointFinalOmega]  = graphsCorrectedYieldShrunkOmega[i]->GetEYhigh()[j];
            yErrorLowFinalOmega[nPointFinalOmega]   = graphsCorrectedYieldShrunkOmega[i]->GetEYlow()[j];
            if (sysAvailOmega[i]){
                Int_t counter = 0;
                while(counter < 100 && TMath::Abs(xValueFinalOmega[nPointFinalOmega] - ptSysRelOmega[i][counter])> 0.001) counter++;
                if (counter < 100){
                    cout << ptSysRelOmega[i][counter]<< "\t found it" << endl;
                    yErrorSysLowFinalOmega[nPointFinalOmega] = TMath::Abs(yErrorSysLowRelOmega[i][counter]/100*graphsCorrectedYieldShrunkOmega[i]->GetY()[j]);
                    yErrorSysHighFinalOmega[nPointFinalOmega] = yErrorSysHighRelOmega[i][counter]/100*graphsCorrectedYieldShrunkOmega[i]->GetY()[j];

                } else {
                    yErrorSysLowFinalOmega[nPointFinalOmega] = 0;
                    yErrorSysHighFinalOmega[nPointFinalOmega] = 0;

                }
            } else {
                yErrorSysLowFinalOmega[nPointFinalOmega] = 0;
                yErrorSysHighFinalOmega[nPointFinalOmega] = 0;
            }
            graphsCorrectedYieldSysShrunkOmega[i]->SetPointEYlow(j,yErrorSysLowFinalOmega[nPointFinalOmega]);
            graphsCorrectedYieldSysShrunkOmega[i]->SetPointEYhigh(j,yErrorSysHighFinalOmega[nPointFinalOmega]);
            nPointFinalOmega++;
        }

        // Set correct trigger order for combination function
        Int_t nCorrOrder    = GetOrderedTrigger(triggerName[i]);
        if (nCorrOrder == -1){
            cout << "ERROR: trigger name not defined" << endl;
            return;
        }

        if ( graphsCorrectedYieldShrunkOmega[i]){
            histoStatOmega[nCorrOrder]    = histoCorrectedYieldOmegaScaledMasked[i];
            graphSystOmega[nCorrOrder]    = graphsCorrectedYieldSysShrunkOmega[i];
            offSetsOmegaSys[nCorrOrder]   = histoStatOmega[nCorrOrder]->GetXaxis()->FindBin(graphSystOmega[nCorrOrder]->GetX()[0])-1;
            if (graphMassOmegaData[i])
                graphOrderedMassOmegaData[nCorrOrder]         = graphMassOmegaData[i];
            if (graphMassOmegaMC[i])
                graphOrderedMassOmegaMC[nCorrOrder]           = graphMassOmegaMC[i];
            if (graphWidthOmegaData[i])
                graphOrderedWidthOmegaData[nCorrOrder]        = graphWidthOmegaData[i];
            if (graphWidthOmegaMC[i])
                graphOrderedWidthOmegaMC[nCorrOrder]          = graphWidthOmegaMC[i];
            if (graphAcceptanceOmega[i])
                graphOrderedAcceptanceOmega[nCorrOrder]       = graphAcceptanceOmega[i];
            if (graphEfficiencyOmega[i])
                graphOrderedEfficiencyOmega[nCorrOrder]       = graphEfficiencyOmega[i];
            if (graphEffTimesAccOmega[i])
                graphOrderedEffTimesAccOmega[nCorrOrder]      = graphEffTimesAccOmega[i];
            if (graphPurityOmega[i])
                graphOrderedPurityOmega[nCorrOrder]           = graphPurityOmega[i];
            for (Int_t k = 0; k<4; k++){
                if (graphEffectSecCorrOmega[k][i])
                    graphOrderedEffectSecCorrOmega[k][nCorrOrder] = graphEffectSecCorrOmega[k][i];
                if (graphEfficiencySecOmega[k][i])
                    graphOrderedEfficiencySecOmega[k][nCorrOrder] = graphEfficiencySecOmega[k][i];
            }

            if (sysAvailSingleOmega[i]){
                for (Int_t k = 0; k < nRelSysErrOmegaSources; k++ ){
                    if (graphRelSysErrOmegaSource[k][i])
                        graphOrderedRelSysErrOmegaSource[k][nCorrOrder]   = graphRelSysErrOmegaSource[k][i];
                }
            }
        }
        if (triggerName[i].Contains("INT7") && optionEnergy.CompareTo("8TeV")==0 && mode == 2)
            offSetsOmegaSys[1]+=3;
        if ((triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA")) && optionEnergy.CompareTo("8TeV")==0 && mode == 4 )
            offSetsOmegaSys[4]+=4;
        if ((triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA")) && optionEnergy.CompareTo("8TeV")==0 && mode == 2 )
            offSetsOmegaSys[4]+=3;
    }


    // create weighted graphs for spectra and supporting graphs
    TString nameWeightsLogFileOmega =     Form("%s/weightsOmega_%s.dat",outputDir.Data(),isMC.Data());
    TGraphAsymmErrors* graphCorrectedYieldWeightedAverageOmegaStat    = NULL;
    TGraphAsymmErrors* graphCorrectedYieldWeightedAverageOmegaSys     = NULL;
    TGraphAsymmErrors* graphCorrectedYieldWeightedAverageOmegaTot     = NULL;
    TGraphAsymmErrors* graphMassOmegaDataWeighted                     = NULL;
    TGraphAsymmErrors* graphMassOmegaMCWeighted                       = NULL;
    TGraphAsymmErrors* graphWidthOmegaDataWeighted                    = NULL;
    TGraphAsymmErrors* graphWidthOmegaMCWeighted                      = NULL;
    TGraphAsymmErrors* graphAcceptanceOmegaWeighted                   = NULL;
    TGraphAsymmErrors* graphEfficiencyOmegaWeighted                   = NULL;
    TGraphAsymmErrors* graphEffTimesAccOmegaWeighted                  = NULL;
    TGraphAsymmErrors* graphPurityOmegaWeighted                       = NULL;
    TGraphAsymmErrors* graphRelSysErrOmegaSourceWeighted[30];
    TGraphAsymmErrors* graphEffectSecCorrOmegaWeighted[4]             = {NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphEfficiencySecOmegaWeighted[4]             = {NULL, NULL, NULL, NULL};

    for (Int_t k = 0; k < 30; k++){
        graphRelSysErrOmegaSourceWeighted[k]                          = NULL;
    }
    // Calculate averaged Omega spectrum & respective supporting graphs according to statistical and systematic errors taking correctly into account the cross correlations
    if (averagedOmega){
        cout << maxNBinsOmega << endl;
        // Calculate average Omega spectrum
        graphCorrectedYieldWeightedAverageOmegaTot        = CombinePtPointsSpectraTriggerCorrMat(    histoStatOmega, graphSystOmega,
                                                                                                   binningOmega,  maxNBinsOmega,
                                                                                                   offSetsOmega ,offSetsOmegaSys,
                                                                                                   graphCorrectedYieldWeightedAverageOmegaStat, graphCorrectedYieldWeightedAverageOmegaSys,
                                                                                                   nameWeightsLogFileOmega.Data(),
                                                                                                   mode, optionEnergy, "Omega", v2ClusterizerMerged,
                                                                                                   fileInputCorrFactors
                                                                                               );

        // preparations for weight readout
        Double_t xValuesReadOmega[100];
        Double_t weightsReadOmega[12][100];
        Int_t availableMeasOmega[12]       = {-1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1};
        Int_t nMeasSetOmega               = nrOfTrigToBeCombOmegaRed;
        Int_t nPtBinsReadOmega            = 0;

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
        ifstream fileWeightsOmega;
        fileWeightsOmega.open(nameWeightsLogFileOmega,ios_base::in);
        cout << "reading" << nameWeightsLogFileOmega << endl;

        while(!fileWeightsOmega.eof() && nPtBinsReadOmega < 100){
            TString garbage = "";
            if (nPtBinsReadOmega == 0){
                fileWeightsOmega >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
                for (Int_t i = 0; i < nMeasSetOmega; i++){
                    fileWeightsOmega >> availableMeasOmega[i] ;
                }
                cout << "read following measurements: ";
                for (Int_t i = 0; i < nMeasSetOmega; i++){
                    cout << availableMeasOmega[i] << "\t" ;
                }
                cout << endl;
            } else {
                fileWeightsOmega >> xValuesReadOmega[nPtBinsReadOmega-1];
                for (Int_t i = 0; i < nMeasSetOmega; i++){
                    fileWeightsOmega >> weightsReadOmega[availableMeasOmega[i]][nPtBinsReadOmega-1] ;
                }
                cout << "read: "<<  nPtBinsReadOmega << "\t"<< xValuesReadOmega[nPtBinsReadOmega-1] << "\t" ;
                for (Int_t i = 0; i < nMeasSetOmega; i++){
                    cout << weightsReadOmega[availableMeasOmega[i]][nPtBinsReadOmega-1] << "\t";
                }
                cout << endl;
            }
            nPtBinsReadOmega++;
        }
        nPtBinsReadOmega = nPtBinsReadOmega-2 ;
        fileWeightsOmega.close();

        // creating & filling the weight graphs
        TGraph* graphWeightsOmega[12];
        for (Int_t i = 0; i < 12; i++){
            graphWeightsOmega[i]                    = NULL;
        }
        for (Int_t i = 0; i < nMeasSetOmega; i++){
            cout << i << "\t" << availableMeasOmega[i] << endl;
            graphWeightsOmega[availableMeasOmega[i]]  = new TGraph(nPtBinsReadOmega,xValuesReadOmega,weightsReadOmega[availableMeasOmega[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsReadOmega; n++){
                if (graphWeightsOmega[availableMeasOmega[i]]->GetY()[bin] == 0) graphWeightsOmega[availableMeasOmega[i]]->RemovePoint(bin);
                else bin++;
            }
            graphWeightsOmega[availableMeasOmega[i]]->Print();
        }


        //  **********************************************************************************************************************
        //  ******************************************* Plotting weights Omega *****************************************************
        //  **********************************************************************************************************************
        Int_t textSizeLabelsPixel = 900*0.04;

        TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);// gives the page size
        DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);

        TH2F * histo2DWeights;
        histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,0.,maxPtGlobalOmega,1000,-0.5,1.1);
        SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DWeights->Draw("copy");

            TLegend* legendWeightsOmega = GetAndSetLegend2(0.12, 0.14, 0.55, 0.14+(0.035*nMeasSetOmega/2*1.35), 32);
            legendWeightsOmega->SetNColumns(2);
            for (Int_t i = 0; i < nMeasSetOmega; i++){
                DrawGammaSetMarkerTGraph(graphWeightsOmega[availableMeasOmega[i]],markerTriggWeighted[availableMeasOmega[i]], sizeTrigg[availableMeasOmega[i]],
                                         colorTriggWeighted[availableMeasOmega[i]], colorTriggWeighted[availableMeasOmega[i]]);
                graphWeightsOmega[availableMeasOmega[i]]->Draw("p,same,e1");
                legendWeightsOmega->AddEntry(graphWeightsOmega[availableMeasOmega[i]],nameTriggerWeighted[availableMeasOmega[i]],"p");
            }
            legendWeightsOmega->Draw();

            TLatex *labelWeightsEnergy = new TLatex(0.7,0.24,collisionSystem.Data());
            SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
            labelWeightsEnergy->SetTextFont(43);
            labelWeightsEnergy->Draw();
            TLatex *labelWeightsOmega = new TLatex(0.7,0.20,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
            SetStyleTLatex( labelWeightsOmega, 0.85*textSizeLabelsPixel,4);
            labelWeightsOmega->SetTextFont(43);
            labelWeightsOmega->Draw();
            TLatex *labelDetProcWeights    = new TLatex(0.7, 0.16,detectionProcess.Data());
            SetStyleTLatex( labelDetProcWeights, 0.85*textSizeLabelsPixel,4);
            labelDetProcWeights->SetTextFont(43);
            labelDetProcWeights->Draw();


    //      DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
            DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/%s_WeightsOmegaTriggers.%s",outputDir.Data(), isMC.Data(), suffix.Data()));
        delete canvasWeights;

        // Calculating relative error for Omega
        for (Int_t i = 0; i < 12; i++){
            if (histoStatOmega[i])
                histoRelStatOmega[i]      = CalculateRelErrUpTH1D( histoStatOmega[i], Form("relativeStatErrorOmega_%s", nameTriggerWeighted[i].Data()));
            if (graphSystOmega[i])
                graphRelSystOmega[i]      = CalculateRelErrUpAsymmGraph( graphSystOmega[i], Form("relativeSysErrorOmega_%s", nameTriggerWeighted[i].Data()));
        }

        TGraphAsymmErrors* graphRelErrorOmegaTot        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageOmegaTot, "relativeTotalErrorOmega");
        while (graphRelErrorOmegaTot->GetY()[0] < 0 ) graphRelErrorOmegaTot->RemovePoint(0);

        TGraphAsymmErrors* graphRelErrorOmegaStat       = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageOmegaStat, "relativeStatErrorOmega");
        while (graphRelErrorOmegaStat->GetY()[0] < 0 ) graphRelErrorOmegaStat->RemovePoint(0);

        TGraphAsymmErrors* graphRelErrorOmegaSys        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageOmegaSys, "relativeSysErrorOmega");
        while (graphRelErrorOmegaSys->GetY()[0] < 0 ) graphRelErrorOmegaSys->RemovePoint(0);

        const char *SysErrDatnameMeanSingleErrCheck = Form("%s/SystematicErrorAveragedSingle%s_Omega_%s_Check.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
        fstream SysErrDatAverSingleCheck;
        SysErrDatAverSingleCheck.precision(4);
        cout << SysErrDatnameMeanSingleErrCheck << endl;
        if(sysAvailSingleOmega[0]){
          SysErrDatAverSingleCheck.open(SysErrDatnameMeanSingleErrCheck, ios::out);
          SysErrDatAverSingleCheck << "pt \t Stat err \t sys err \t tot err " << endl;
          for (Int_t i = 0; i < graphRelErrorOmegaTot->GetN(); i++){
              if (graphRelErrorOmegaStat->GetY()[i] > 0) SysErrDatAverSingleCheck << graphRelErrorOmegaStat->GetX()[i] << "\t" << graphRelErrorOmegaStat->GetY()[i] <<"\t" << graphRelErrorOmegaSys->GetY()[i] <<  "\t" << graphRelErrorOmegaTot->GetY()[i] << endl;
          }
          SysErrDatAverSingleCheck << endl;
          SysErrDatAverSingleCheck.close();
        }

        // plot sys relative errors for individual triggers
        TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);

        TH2F * histo2DRelSysErr;
        histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.,maxPtGlobalOmega,1000,0,60.5);
        SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelSysErr->Draw("copy");
            TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.035*(nMeasSetOmega+1)/2), 0.95, 0.92, 32);
            legendRelSysErr->SetNColumns(2);
            for (Int_t i = 0; i < nMeasSetOmega; i++){
                cout << "plotting graph: " << availableMeasOmega[i] << "\t" <<graphRelSystOmega[availableMeasOmega[i]]->GetName() << endl;
                DrawGammaSetMarkerTGraph(graphRelSystOmega[availableMeasOmega[i]], markerTriggWeighted[availableMeasOmega[i]], sizeTrigg[availableMeasOmega[i]],
                                         colorTriggWeighted[availableMeasOmega[i]], colorTriggWeighted[availableMeasOmega[i]]);
                graphRelSystOmega[availableMeasOmega[i]]->Draw("p,same,z");
                legendRelSysErr->AddEntry(graphRelSystOmega[availableMeasOmega[i]],nameTriggerWeighted[availableMeasOmega[i]],"p");
            }
            legendRelSysErr->Draw();

            TLatex *labelRelErrEnergy    = new TLatex(0.15,0.89,collisionSystem.Data());
            SetStyleTLatex( labelRelErrEnergy, 0.85*textSizeLabelsPixel,4);
            labelRelErrEnergy->SetTextFont(43);
            labelRelErrEnergy->Draw();
            TLatex *labelRelErrOmega       = new TLatex(0.15,0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
            SetStyleTLatex( labelRelErrOmega, 0.85*textSizeLabelsPixel,4);
            labelRelErrOmega->SetTextFont(43);
            labelRelErrOmega->Draw();
            TLatex *labelDetProcRelErr    = new TLatex(0.15, 0.81,detectionProcess.Data());
            SetStyleTLatex( labelDetProcRelErr, 0.85*textSizeLabelsPixel,4);
            labelDetProcRelErr->SetTextFont(43);
            labelDetProcRelErr->Draw();

        canvasRelSysErr->SaveAs(Form("%s/Omega_RelSysErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));


            DrawGammaSetMarkerTGraphAsym(graphRelErrorOmegaSys, 24, 1.5, kGray+1 , kGray+1);
//             graphRelErrorOmegaSys->SetLineStyle(7);
            graphRelErrorOmegaSys->Draw("same,pze1");
            legendRelSysErr->AddEntry(graphRelErrorOmegaSys,"average","p");
            legendRelSysErr->Draw();

            labelRelErrEnergy->Draw();
            labelRelErrOmega->Draw();
            labelDetProcRelErr->Draw();

        canvasRelSysErr->SaveAs(Form("%s/Omega_RelSysErrWithAverage_SingleMeas.%s",outputDir.Data(),suffix.Data()));


        // plot stat relative errors for individual triggers
        TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);

        TH2F * histo2DRelStatErr;
        histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.,maxPtGlobalOmega,1000,0,60.5);
        SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelStatErr->Draw("copy");

            TLegend* legendRelStatErr       = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetOmega/2), 0.95, 0.92, 32);
            legendRelStatErr->SetNColumns(2);
            for (Int_t i = 0; i < nMeasSetOmega; i++){
                 cout << "plotting graph: " << availableMeasOmega[i] << "\t" <<histoRelStatOmega[availableMeasOmega[i]]->GetName() << endl;
                if (histoRelStatOmega[availableMeasOmega[i]] && mode == 2){
                    TGraphAsymmErrors* dummyGraph = new TGraphAsymmErrors(histoRelStatOmega[availableMeasOmega[i]]);
                    dummyGraph->Print();
                    DrawGammaSetMarkerTGraph(dummyGraph, markerTriggWeighted[availableMeasOmega[i]], sizeTrigg[availableMeasOmega[i]],
                                         colorTriggWeighted[availableMeasOmega[i]], colorTriggWeighted[availableMeasOmega[i]]);
                    dummyGraph->Draw("pX,same");
                    legendRelStatErr->AddEntry(dummyGraph,nameTriggerWeighted[availableMeasOmega[i]],"p");

                     for (Int_t j = 1; j < histoRelStatOmega[availableMeasOmega[i]]->GetNbinsX()+1; j++){
                        cout << j << ": " << histoRelStatOmega[availableMeasOmega[i]]->GetBinContent(j) << endl;
                     }
                } else {
                    DrawGammaSetMarker(histoRelStatOmega[availableMeasOmega[i]],markerTriggWeighted[availableMeasOmega[i]], sizeTrigg[availableMeasOmega[i]],
                                            colorTriggWeighted[availableMeasOmega[i]], colorTriggWeighted[availableMeasOmega[i]]);
                    histoRelStatOmega[availableMeasOmega[i]]->DrawCopy("p,same,z");
                    legendRelStatErr->AddEntry(histoRelStatOmega[availableMeasOmega[i]],nameTriggerWeighted[availableMeasOmega[i]],"p");
                }
            }
            legendRelStatErr->Draw();

            labelRelErrEnergy->Draw();
            labelRelErrOmega->Draw();
            labelDetProcRelErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/Omega_RelStatErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));

        // plot full error for final result decomposed
        TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);

        TH2F * histo2DRelTotErrOmega;
        histo2DRelTotErrOmega                 = new TH2F("histo2DRelTotErrOmega","histo2DRelTotErrOmega",11000,0.,maxPtGlobalOmega,1000,0,40.5);
        SetStyleHistoTH2ForGraphs(histo2DRelTotErrOmega, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelTotErrOmega->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphRelErrorOmegaTot, 20, 1.5, kBlack , kBlack);
            graphRelErrorOmegaTot->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphRelErrorOmegaStat, 24, 1.5, kGray+2 , kGray+2);
            graphRelErrorOmegaStat->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphRelErrorOmegaSys, 24, 1.5, kGray+1 , kGray+1);
            graphRelErrorOmegaSys->SetLineStyle(7);
            graphRelErrorOmegaSys->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr2       = GetAndSetLegend2(0.72, 0.92-(0.035*3), 0.9, 0.92, 32);
            legendRelTotErr2->AddEntry(graphRelErrorOmegaTot,"tot","p");
            legendRelTotErr2->AddEntry(graphRelErrorOmegaStat,"stat","l");
            legendRelTotErr2->AddEntry(graphRelErrorOmegaSys,"sys","l");
            legendRelTotErr2->Draw();

            labelRelErrEnergy->Draw();
            labelRelErrOmega->Draw();
            labelDetProcRelErr->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Omega_RelErrorsFulldecomp.%s",outputDir.Data(),suffix.Data()));


        // Calculate relative sys error weighted
        if (sysAvailSingleOmega[0]){
            for (Int_t k = 0; k< nRelSysErrOmegaSources ; k++ ){
                graphRelSysErrOmegaSourceWeighted[k]      = CalculateWeightedQuantity(    graphOrderedRelSysErrOmegaSource[k],
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                   );
                if (!graphRelSysErrOmegaSourceWeighted[k]){
                    cout << "Aborted in CalculateWeightedQuantity for " << endl;
                    return;
                } else {
                    graphRelSysErrOmegaSourceWeighted[k]->SetName(Form("RelSysErrOmegaSourceWeighted%s", ((TString)ptSysDetail[0][0].at(k+1)).Data()));
                    while (graphRelSysErrOmegaSourceWeighted[k]->GetY()[0] == -10000 )   graphRelSysErrOmegaSourceWeighted[k]->RemovePoint(0);
                }
            }
        }
//         return;

        // Calculation of averaged supporting plots with weights from spectra
        if (mode != 10){
            graphMassOmegaDataWeighted                    = CalculateWeightedQuantity(    graphOrderedMassOmegaData,
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                   );
            if (!graphMassOmegaDataWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
            graphMassOmegaMCWeighted                      = CalculateWeightedQuantity(    graphOrderedMassOmegaMC,
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                   );
            if (!graphMassOmegaMCWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }

            graphWidthOmegaDataWeighted                   = CalculateWeightedQuantity(    graphOrderedWidthOmegaData,
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                   );
            if (!graphWidthOmegaDataWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }

            graphWidthOmegaMCWeighted                     = CalculateWeightedQuantity(    graphOrderedWidthOmegaMC,
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                   );
            if (!graphWidthOmegaMCWeighted){
                cout << "Aborted in CalculateWeightedQuantity" << endl;
                return;
            }
        }

        cout << "weighting Omega acceptance" << endl;
        graphAcceptanceOmegaWeighted                      = CalculateWeightedQuantity(    graphOrderedAcceptanceOmega,
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                    );
        if (!graphAcceptanceOmegaWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }


        for (Int_t k = 0; k< 4; k++){
            cout << "calculating effective sec corr: " << k << endl;
            if (hasSecCorrFac[k]){
                graphEffectSecCorrOmegaWeighted[k]        = CalculateWeightedQuantity(    graphOrderedEffectSecCorrOmega[k],
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                    );
                if (!graphEffectSecCorrOmegaWeighted[k]){
                    cout << "Aborted in CalculateWeightedQuantity for effective sec cor " << k << endl;
                    return;
                }
            }
//             return;
            cout << "calculating effieciency sec Omega: " << k << endl;
            if (hasSecEffi[k]){
                graphEfficiencySecOmegaWeighted[k]        = CalculateWeightedQuantity(    graphOrderedEfficiencySecOmega[k],
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                    );
                if (!graphEffectSecCorrOmegaWeighted[k]){
                    cout << "Aborted in CalculateWeightedQuantity for sec effi " << k << endl;
                    return;
                }
            }
        }


        cout << "weighting Omega efficiency" << endl;
        graphEfficiencyOmegaWeighted                      = CalculateWeightedQuantity(    graphOrderedEfficiencyOmega,
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                    );
        if (!graphEfficiencyOmegaWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }
        cout << "weighting Omega efficiency x acceptance" << endl;
        graphEffTimesAccOmegaWeighted                     = CalculateWeightedQuantity(    graphOrderedEffTimesAccOmega,
                                                                                        graphWeightsOmega,
                                                                                        binningOmega,  maxNBinsOmega,
                                                                                        MaxNumberOfFiles
                                                                                    );

        if (!graphEffTimesAccOmegaWeighted){
            cout << "Aborted in CalculateWeightedQuantity" << endl;
            return;
        }


        // remove points in spectrum which should have been masked
        if (mode != 10){
            if (graphMassOmegaDataWeighted)
                while (graphMassOmegaDataWeighted->GetY()[0] == -10000 )   graphMassOmegaDataWeighted->RemovePoint(0);
            else
                cout << "I don't have a weighted Omega mass data graph" << endl;
            if (graphMassOmegaMCWeighted)
                while (graphMassOmegaMCWeighted->GetY()[0] == -10000)     graphMassOmegaMCWeighted->RemovePoint(0);
            else
                cout << "I don't have a weighted Omega mass MC graph" << endl;
            if (graphWidthOmegaDataWeighted)
                while (graphWidthOmegaDataWeighted->GetY()[0] == -10000)  graphWidthOmegaDataWeighted->RemovePoint(0);
            else
                cout << "I don't have a weighted Omega width data graph" << endl;
            if (graphWidthOmegaMCWeighted)
                while (graphWidthOmegaMCWeighted->GetY()[0] == -10000)    graphWidthOmegaMCWeighted->RemovePoint(0);
            else
                cout << "I don't have a weighted Omega width MC graph" << endl;
        }

        if (graphAcceptanceOmegaWeighted)
            while (graphAcceptanceOmegaWeighted->GetY()[0] == -10000)     graphAcceptanceOmegaWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted Omega acceptance graph" << endl;

        if (graphEfficiencyOmegaWeighted)
            while (graphEfficiencyOmegaWeighted->GetY()[0] == -10000)     graphEfficiencyOmegaWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted Omega efficiency graph" << endl;

        if (graphEffTimesAccOmegaWeighted)
            while (graphEffTimesAccOmegaWeighted->GetY()[0] == -10000)    graphEffTimesAccOmegaWeighted->RemovePoint(0);
        else
            cout << "I don't have a weighted Omega acceptance x efficiency graph" << endl;

        for (Int_t k = 0; k< 4; k++){
            if (graphEffectSecCorrOmegaWeighted[k])
                while (graphEffectSecCorrOmegaWeighted[k]->GetY()[0] == -10000)    graphEffectSecCorrOmegaWeighted[k]->RemovePoint(0);
            else
                cout << "I don't have a weighted Omega effective sec corr graph for "<< k << endl;
               if (graphEfficiencySecOmegaWeighted[k])
                while (graphEfficiencySecOmegaWeighted[k]->GetY()[0] == -10000)    graphEfficiencySecOmegaWeighted[k]->RemovePoint(0);
            else
                cout << "I don't have a weighted sec Omega efficiency  graph for "<< k << endl;
        }

        //  **********************************************************************************************************************
        //  **************************************** Combine+write detailed Systematics ******************************************
        //  **********************************************************************************************************************

        const char *SysErrDatnameMeanSingleErr = Form("%s/SystematicErrorAveragedSingle%s_Omega_%s.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
        fstream SysErrDatAverSingle;
        SysErrDatAverSingle.precision(4);
        cout << SysErrDatnameMeanSingleErr << endl;
        if(sysAvailSingleOmega[0] && graphRelSysErrOmegaSourceWeighted[0]){
            SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
            for(Int_t iColumn = 0; iColumn < (Int_t)ptSysDetail[0][0].size(); iColumn++) SysErrDatAverSingle << ptSysDetail[0][0].at(iColumn) << "\t";
            SysErrDatAverSingle << endl;
            for (Int_t i = 1; i < nrOfTrigToBeComb; i++){
                if(!sysAvailSingleOmega[i]) continue;
                for(Int_t iCol = 0; iCol < (Int_t)ptSysDetail[i][0].size(); iCol++){
                    if( ((TString)ptSysDetail[i][0].at(iCol)).CompareTo(((TString)ptSysDetail[0][0].at(iCol))) ){
                        cout << "ERROR: Systematic error type at pos " << iCol << " does not agree for " << availableMeasOmega[i] << " & " << availableMeasOmega[0] << ", returning!" << endl;
                        return;
                    }
                }
            }

            for(Int_t i=0; i<graphRelSysErrOmegaSourceWeighted[0]->GetN(); i++){
                SysErrDatAverSingle << graphRelSysErrOmegaSourceWeighted[0]->GetX()[i] << "\t";
                Int_t nColumns = (Int_t)ptSysDetail[0][0].size();
                for(Int_t iErr=0; iErr<nColumns-1; iErr++)
                    SysErrDatAverSingle << graphRelSysErrOmegaSourceWeighted[iErr]->GetY()[i] << "\t";
                SysErrDatAverSingle << endl;

            }
        }
        SysErrDatAverSingle.close();

        // ***************************************************************************************************
        // ********************* Plot all mean erros separately after smoothing ******************************
        // ***************************************************************************************************
        if(sysAvailSingleOmega[0] && graphRelSysErrOmegaSourceWeighted[0]){
            TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);// gives the page size
            DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);

                // create dummy histo
                TH2D *histo2DNewSysErrMean ;
                histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 100,0.,maxPtGlobalOmega,1000.,-0.5,30.);
                SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                        1,0.9, 510, 510);
                histo2DNewSysErrMean->Draw();

                // Give legend position for plotting
                Double_t minXLegend     = 0.12;
                Double_t maxYLegend     = 0.95;
                Double_t widthLegend    = 0.25;
                if (nRelSysErrOmegaSources> 7)
                    widthLegend         = 0.5;
                Double_t heightLegend   = 1.05* 0.035 * (nRelSysErrOmegaSources+3);
                if (nRelSysErrOmegaSources> 7)
                    heightLegend        = 1.05* 0.035 * (nRelSysErrOmegaSources/2+1);

                // create legend
                TLegend* legendMeanNew = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
                legendMeanNew->SetMargin(0.1);
                if (nRelSysErrOmegaSources> 7) legendMeanNew->SetNColumns(2);

                for(Int_t i = 0;i< nRelSysErrOmegaSources-1 ;i++){
                    DrawGammaSetMarkerTGraphAsym(graphRelSysErrOmegaSourceWeighted[i], GetMarkerStyleSystematics(  ptSysDetail[0][0].at(i+1), mode), 1.,
                                                GetColorSystematics( ptSysDetail[0][0].at(i+1), mode),GetColorSystematics( ptSysDetail[0][0].at(i+1), mode));
                    graphRelSysErrOmegaSourceWeighted[i]->Draw("pX0,csame");
                    legendMeanNew->AddEntry(graphRelSysErrOmegaSourceWeighted[i],GetSystematicsName(ptSysDetail[0][0].at(i+1)),"p");
                }

                DrawGammaSetMarkerTGraphAsym(graphRelSysErrOmegaSourceWeighted[nRelSysErrOmegaSources-1], 20, 1.,kBlack,kBlack);
                graphRelSysErrOmegaSourceWeighted[nRelSysErrOmegaSources-1]->Draw("p,csame");
                legendMeanNew->AddEntry(graphRelSysErrOmegaSourceWeighted[nRelSysErrOmegaSources-1],"quad. sum.","p");
                legendMeanNew->Draw();

                // labeling
                TLatex *labelEnergySysDetailed = new TLatex(0.7, 0.93,collisionSystem.Data());
                labelEnergySysDetailed->SetTextAlign(31);
                SetStyleTLatex( labelEnergySysDetailed, 0.85*textSizeSpectra,4);
                labelEnergySysDetailed->Draw();

                TLatex *labelOmegaSysDetailed     = new TLatex(0.7, 0.93-0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
                labelOmegaSysDetailed->SetTextAlign(31);
                SetStyleTLatex( labelOmegaSysDetailed, 0.85*textSizeSpectra,4);
                labelOmegaSysDetailed->Draw();

                TLatex *labelDetProcSysDetailed = new TLatex(0.7, 0.93-2*0.99*textSizeSpectra*0.85,detectionProcess.Data());
                labelDetProcSysDetailed->SetTextAlign(31);
                SetStyleTLatex( labelDetProcSysDetailed, 0.85*textSizeSpectra,4);
                labelDetProcSysDetailed->Draw();

            canvasNewSysErrMean->Update();
            canvasNewSysErrMean->SaveAs(Form("%s/Omega_SysErrorsSeparatedSourcesReweighted_%s.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
        }

        for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++){
            for(Int_t iB=0; iB<50; iB++) ptSysDetail[iR][iB].clear();
        }


    // if averaging wasn't enabled pick values according to predefined ranges ("cherry picking points")
    } else {
       graphCorrectedYieldWeightedAverageOmegaStat        = new TGraphAsymmErrors(nPointFinalOmega, xValueFinalOmega, yValueFinalOmega,
                                                                                xErrorLowFinalOmega, xErrorHighFinalOmega,yErrorLowFinalOmega, yErrorHighFinalOmega);
       graphCorrectedYieldWeightedAverageOmegaSys         = new TGraphAsymmErrors(nPointFinalOmega, xValueFinalOmega, yValueFinalOmega,
                                                                                xErrorLowFinalOmega, xErrorHighFinalOmega,yErrorSysLowFinalOmega, yErrorSysHighFinalOmega);

        if (mode != 10){
            graphMassOmegaDataWeighted                    = graphMassOmegaData[0];
            graphMassOmegaMCWeighted                      = graphMassOmegaMC[0];
            graphWidthOmegaDataWeighted                   = graphWidthOmegaData[0];
            graphWidthOmegaMCWeighted                     = graphWidthOmegaMC[0];
        }

        graphAcceptanceOmegaWeighted                      = graphAcceptanceOmega[0];
        graphEfficiencyOmegaWeighted                      = graphEfficiencyOmega[0];
        graphEffTimesAccOmegaWeighted                     = graphEffTimesAccOmega[0];

        for (Int_t k = 0; k<4; k++){
            if (graphEffectSecCorrOmega[k][0]) graphEffectSecCorrOmegaWeighted[k]  = graphEffectSecCorrOmega[k][0];
            if (graphEfficiencySecOmega[k][0]) graphEfficiencySecOmegaWeighted[k]  = graphEfficiencySecOmega[k][0];
        }

        TGraphAsymmErrors* graphRelErrorOmegaStat       = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageOmegaStat, "relativeStatErrorOmega");
        while (graphRelErrorOmegaStat->GetY()[0] < 0 ) graphRelErrorOmegaStat->RemovePoint(0);

        TGraphAsymmErrors* graphRelErrorOmegaSys        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageOmegaSys, "relativeSysErrorOmega");
        while (graphRelErrorOmegaSys->GetY()[0] < 0 ) graphRelErrorOmegaSys->RemovePoint(0);

        const char *SysErrDatnameMeanSingleErrCheck = Form("%s/SystematicErrorAveragedSingle%s_Omega_%s_Check.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
        fstream SysErrDatAverSingleCheck;
        SysErrDatAverSingleCheck.precision(4);
        cout << SysErrDatnameMeanSingleErrCheck << endl;

        SysErrDatAverSingleCheck.open(SysErrDatnameMeanSingleErrCheck, ios::out);
        SysErrDatAverSingleCheck << "pt \t Stat err \t sys err \t tot err " << endl;
        for (Int_t i = 0; i < graphRelErrorOmegaStat->GetN(); i++){
            if (graphRelErrorOmegaStat->GetY()[i] > 0) SysErrDatAverSingleCheck << graphRelErrorOmegaStat->GetX()[i] << "\t" << graphRelErrorOmegaStat->GetY()[i] <<"\t" << graphRelErrorOmegaSys->GetY()[i] << endl;
        }
        SysErrDatAverSingleCheck << endl;
        SysErrDatAverSingleCheck.close();

    }
    // print final graphs
    cout << "stat Omega" << endl;
    graphCorrectedYieldWeightedAverageOmegaStat->Print();
    cout << "sys Omega" << endl;
    if (graphCorrectedYieldWeightedAverageOmegaSys) graphCorrectedYieldWeightedAverageOmegaSys->Print();

    if (graphCorrectedYieldWeightedAverageOmegaStat){
        while (graphCorrectedYieldWeightedAverageOmegaStat->GetX()[0]< minPtGlobalOmega){
            graphCorrectedYieldWeightedAverageOmegaStat->RemovePoint(0);
        }
    }

    if (graphCorrectedYieldWeightedAverageOmegaSys){
        while (graphCorrectedYieldWeightedAverageOmegaSys->GetX()[0]< minPtGlobalOmega){
            graphCorrectedYieldWeightedAverageOmegaSys->RemovePoint(0);
        }
    }

    if (mode != 10){
        //***************************************************************************************************************
        //************************************Plotting Mass Omega reduced range  ******************************************
        //***************************************************************************************************************
        canvasMass->cd();
        histo2DMassOmega->DrawCopy();

        legendMassRedOmega->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0) ||
             (optionEnergy.CompareTo("pPb_5.023TeV")==0)
             ){
                if (graphMassOmegaData[i] && !maskedFullyOmega[i]) {
                    DrawGammaSetMarkerTGraphAsym(graphMassOmegaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    graphMassOmegaData[i]->Draw("p,e1,same");
                    legendMassRedOmega->AddEntry(graphMassOmegaData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
                }
                if (graphMassOmegaMC[i] && !maskedFullyOmega[i]){
                    DrawGammaSetMarkerTGraphAsym(graphMassOmegaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    graphMassOmegaMC[i]->Draw("p,e1,same");
                    legendMassRedOmega->AddEntry(graphMassOmegaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
                }
            }
        }
        legendMassRedOmega->Draw();
        labelEnergyMass->Draw();
        labelOmegaMass->Draw();
        labelDetProcMass->Draw();

        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/Omega_%s_Mass1_Reduced_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

        if ( optionEnergy.CompareTo("2.76TeV")==0 ){
            canvasMass->cd();
            histo2DMassOmega->DrawCopy();
            legendMassRedOmega2->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if( i==1 || i==3 || i==4 || i==5){
                    if (graphMassOmegaData[i] && !maskedFullyOmega[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphMassOmegaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                        graphMassOmegaData[i]->Draw("p,e1,same");
                        legendMassRedOmega2->AddEntry(graphMassOmegaData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
                    }
                    if (graphMassOmegaMC[i] && !maskedFullyOmega[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphMassOmegaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                        graphMassOmegaMC[i]->Draw("p,e1,same");
                        legendMassRedOmega2->AddEntry(graphMassOmegaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
                    }
                }
            }
            legendMassRedOmega2->Draw();
            labelEnergyMass->Draw();
            labelOmegaMass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Omega_%s_Mass2_Reduced_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
        }
        //***************************************************************************************************************
        //********************************* Omega Mass weighted ***********************************************************
        //***************************************************************************************************************
        if (graphMassOmegaDataWeighted || graphMassOmegaMCWeighted){
            canvasMass->cd();
            histo2DMassOmega->DrawCopy();
            TLegend* legendMassOmegaWeighted = GetAndSetLegend2(0.52, 0.88, 0.75, 0.88+(1.05*2*0.85*textSizeSpectra),28);
            if (graphMassOmegaDataWeighted){
                DrawGammaSetMarkerTGraphAsym(graphMassOmegaDataWeighted, 20, 1, kBlack, kBlack);
                graphMassOmegaDataWeighted->Draw("p,e1,same");
                legendMassOmegaWeighted->AddEntry(graphMassOmegaDataWeighted, "Data", "p");
            }
            if (graphMassOmegaDataWeighted){
                DrawGammaSetMarkerTGraphAsym(graphMassOmegaMCWeighted, 24, 1, kGray+2, kGray+2);
                graphMassOmegaMCWeighted->Draw("p,e1,same");
                legendMassOmegaWeighted->AddEntry(graphMassOmegaMCWeighted, "MC", "p");
            }
            legendMassOmegaWeighted->Draw();
            labelEnergyMass->Draw();
            labelOmegaMass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Omega_%s_Mass_Weighted_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

            TGraphAsymmErrors* graphMassDifferenceOmegaDatavsMC       = CalculateAsymGraphDifferenceToGraph(graphMassOmegaDataWeighted,graphMassOmegaMCWeighted);
            graphMassDifferenceOmegaDatavsMC->Print();

            TGraphAsymmErrors* graphMassRelDifferenceOmegaDatavsMC    = CalculateAsymGraphRatioToGraph(graphMassDifferenceOmegaDatavsMC,graphMassOmegaMCWeighted);
            graphMassRelDifferenceOmegaDatavsMC->Print();
            graphMassRelDifferenceOmegaDatavsMC = ScaleGraph(graphMassRelDifferenceOmegaDatavsMC, 100.);

            canvasMass->cd();
            TH2F * histo2DRelMassDiffOmega       = new TH2F("histo2DRelMassDiffOmega","histo2DRelMassDiffOmega",1000,0., maxPtGlobalOmega,10000,-3, 3);
            SetStyleHistoTH2ForGraphs(histo2DRelMassDiffOmega, "#it{p}_{T} (GeV/#it{c})","#it{M}_{#omega, data}-#it{M}_{#omega, MC}/ #it{M}_{#omega, MC} (%)",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.4);
            histo2DRelMassDiffOmega->DrawCopy();

            if (graphMassRelDifferenceOmegaDatavsMC){
              TF1 *fitOmegaRelMassDiff = new TF1("fitOmegaRelMassDiff","[0]",0.,maxPtGlobalOmega);
              fitOmegaRelMassDiff->SetParameter(0,0.);
              graphMassRelDifferenceOmegaDatavsMC->Fit(fitOmegaRelMassDiff,"QNRMEX0+","",0.,maxPtGlobalOmega);
              DrawGammaSetMarkerTF1( fitOmegaRelMassDiff, 1, 0.5, kRed);
              fitOmegaRelMassDiff->Draw("same");

              TLegend* legendOmegaRelMassDiff = GetAndSetLegend2(0.15, 0.12, 0.6, 0.12+(0.035*1), 0.035, 2, "", 42, 0.15);
              legendOmegaRelMassDiff->AddEntry(fitOmegaRelMassDiff,"fit with constant");
              legendOmegaRelMassDiff->AddEntry((TObject*)0,Form("%0.4f #pm %0.4f", fitOmegaRelMassDiff->GetParameter(0), fitOmegaRelMassDiff->GetParError(0)),"");
              legendOmegaRelMassDiff->Draw();

              DrawGammaSetMarkerTGraphAsym(graphMassRelDifferenceOmegaDatavsMC, 20, 1, kBlack, kBlack);
              graphMassRelDifferenceOmegaDatavsMC->Draw("p,e1,same");

              fileFitsOutput << "average rel mass diff: " << fitOmegaRelMassDiff->GetParameter(0) << "+-"<< fitOmegaRelMassDiff->GetParError(0) << endl;
            }
            DrawGammaLines(0., maxPtGlobalOmega , 0., 0., 1, kGray+2, 7);
//             legendMassOmegaWeighted->Draw();
            labelEnergyMass->Draw();
            labelOmegaMass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Omega_%s_RelMassDiff_Weighted.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

        }
        //***************************************************************************************************************
        //************************************Plotting Width Omega reduced range  *****************************************
        //***************************************************************************************************************
        canvasWidth->cd();
        histo2DWidthOmega->DrawCopy();

        legendWidthRedOmega->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
          if((optionEnergy.CompareTo("2.76TeV")==0 && (i==0 || i==2)) ||
             (optionEnergy.CompareTo("8TeV")==0) ||
             (optionEnergy.CompareTo("pPb_5.023TeV")==0)
             ){
                if (graphWidthOmegaData[i] && !maskedFullyOmega[i]) {
                    DrawGammaSetMarkerTGraphAsym(graphWidthOmegaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    graphWidthOmegaData[i]->Draw("p,e1,same");
                    legendWidthRedOmega->AddEntry(graphWidthOmegaData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
                }
                if (graphWidthOmegaMC[i] && !maskedFullyOmega[i]){
                    DrawGammaSetMarkerTGraphAsym(graphWidthOmegaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                    graphWidthOmegaMC[i]->Draw("p,e1,same");
                    legendWidthRedOmega->AddEntry(graphWidthOmegaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
                }
            }
        }
        legendWidthRedOmega->Draw();
        labelEnergyWidth->Draw();
        labelOmegaWidth->Draw();
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Omega_%s_Width1_Reduced_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

        if (optionEnergy.CompareTo("2.76TeV")==0){
            canvasWidth->cd();
            histo2DWidthOmega->DrawCopy();
            legendWidthRedOmega2->SetNColumns(2);
            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if( i==1 || i==3 || i==4 || i==5 ){
                    if (graphWidthOmegaData[i] && !maskedFullyOmega[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphWidthOmegaData[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                        graphWidthOmegaData[i]->Draw("p,e1,same");
                        legendWidthRedOmega2->AddEntry(graphWidthOmegaData[i], Form("%s data",triggerNameLabel[i].Data()), "p");
                    }
                    if (graphWidthOmegaMC[i] && !maskedFullyOmega[i]) {
                        DrawGammaSetMarkerTGraphAsym(graphWidthOmegaMC[i], markerTriggMC[i], sizeTrigg[i], colorTriggShade[i], colorTriggShade[i]);
                        graphWidthOmegaMC[i]->Draw("p,e1,same");
                        legendWidthRedOmega2->AddEntry(graphWidthOmegaMC[i], Form("%s MC", triggerNameLabel[i].Data()), "p");
                    }
                }
            }
            legendWidthRedOmega2->Draw();
            labelEnergyWidth->Draw();
            labelOmegaWidth->Draw();
            labelDetProcWidth->Draw();

            canvasWidth->Update();
            canvasWidth->SaveAs(Form("%s/Omega_%s_Width2_Reduced_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
        }
        //***************************************************************************************************************
        //********************************* Omega Width weighted **********************************************************
        //***************************************************************************************************************
        canvasWidth->cd();
        histo2DWidthOmega->DrawCopy();
        TLegend* legendWidthOmegaWeighted = GetAndSetLegend2(0.52, 0.86, 0.75, 0.86+(1.05*2*0.85*textSizeSpectra),28);
        if (graphWidthOmegaDataWeighted){
            DrawGammaSetMarkerTGraphAsym(graphWidthOmegaDataWeighted, 20, 1, kBlack, kBlack);
            graphWidthOmegaDataWeighted->Draw("p,e1,same");
            legendWidthOmegaWeighted->AddEntry(graphWidthOmegaDataWeighted, "Data", "p");
        }
        if (graphWidthOmegaDataWeighted){
            DrawGammaSetMarkerTGraphAsym(graphWidthOmegaMCWeighted, 24, 1, kGray+2, kGray+2);
            graphWidthOmegaMCWeighted->Draw("p,e1,same");
            legendWidthOmegaWeighted->AddEntry(graphWidthOmegaMCWeighted, "MC", "p");
        }
        legendWidthOmegaWeighted->Draw();
        labelEnergyWidth->Draw();
        labelOmegaWidth->Draw();
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Omega_%s_Width_Weighted_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
    }
    if (!enableEta) delete canvasMass;
    if (!enableEta) delete canvasWidth;

    //***************************************************************************************************************
    //************************************* Efficiency weighted *****************************************************
    //***************************************************************************************************************
    if (graphEfficiencyOmegaWeighted){
        DrawGammaCanvasSettings( canvasEffi, 0.09, 0.017, 0.037, 0.08);
        canvasEffi->SetLogy(0);
        canvasEffi->cd();
        histo2DEffiOmega->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphEfficiencyOmegaWeighted, 20, 1, kBlack, kBlack);
        graphEfficiencyOmegaWeighted->Draw("p,e1,same");

        TLatex *labelEnergyEffiWOTrigg = new TLatex(0.62, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyEffiWOTrigg, 0.85*textSizeSpectra,4);
        labelEnergyEffiWOTrigg->Draw();

        TLatex *labelOmegaEffiWOTrigg = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
        SetStyleTLatex( labelOmegaEffiWOTrigg, 0.85*textSizeSpectra,4);
        labelOmegaEffiWOTrigg->Draw();

        TLatex *labelDetProcEffiWOTrigg = new TLatex(0.62, 0.15,detectionProcess.Data());
        SetStyleTLatex( labelDetProcEffiWOTrigg, 0.85*textSizeSpectra,4);
        labelDetProcEffiWOTrigg->Draw();

        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/Omega_EfficiencyW0TriggEff_Weighted_%i.%s",outputDir.Data(),mode,suffix.Data()));
    }

    //***************************************************************************************************************
    //***************************** Secondary corr factors weighted all separate ************************************
    //***************************************************************************************************************
    for (Int_t k = 0; k< 4; k++ ){
        if (graphEfficiencySecOmegaWeighted[k]){
            canvasEffi->cd();
            histo2DEffiSecOmega[k]->DrawCopy();

            DrawGammaSetMarkerTGraphAsym(graphEfficiencySecOmegaWeighted[k], 20, 1, kGray+2, kGray+2);
            graphEfficiencySecOmegaWeighted[k]->Draw("p,e1,same");


            TLatex* labelEnergySecEff    = new TLatex(0.14, 0.84+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
            SetStyleTLatex( labelEnergySecEff, 0.85*textSizeSpectra,4);

            TLatex* labelOmegaSecEff       = new TLatex(0.14, 0.84+0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
            SetStyleTLatex( labelOmegaSecEff, 0.85*textSizeSpectra,4);

            TLatex* labelDetProcSecEff   = new TLatex(0.14, 0.84,detectionProcess.Data());
            SetStyleTLatex( labelDetProcSecEff, 0.85*textSizeSpectra,4);


            labelEnergySecEff->Draw();
            labelOmegaSecEff->Draw();
            labelDetProcSecEff->Draw();


            canvasEffi->Update();
            canvasEffi->SaveAs(Form("%s/Omega_SecEfficiencyFrom%s_Weighted.%s",outputDir.Data(), nameSecOmegaPartRead[k].Data(), suffix.Data()));
            delete labelEnergySecEff;
            delete labelOmegaSecEff;
            delete labelDetProcSecEff;
        }
        if (graphEffectSecCorrOmegaWeighted[k]){
            canvasEffi->cd();
            canvasEffi->SetLogy(0);
            canvasEffi->SetLeftMargin(0.1);
            canvasEffi->SetTopMargin(0.04);
            histo2DEffectiveSecCorr[k]->DrawCopy();

            DrawGammaSetMarkerTGraphAsym(graphEffectSecCorrOmegaWeighted[k], 20, 1, kGray+2, kGray+2);
            graphEffectSecCorrOmegaWeighted[k]->Draw("p,e1,same");

            TLatex *labelEnergyEffSec = new TLatex(0.95, 0.93-(0.98*textSizeSpectra*0.85),collisionSystem.Data());
            SetStyleTLatex( labelEnergyEffSec, 0.85*textSizeSpectra,4,1, 42, kTRUE, 31);
            labelEnergyEffSec->Draw();

            TLatex *labelOmegaEffSec = new TLatex(0.95,  0.93-(2*textSizeSpectra*0.85),"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
            SetStyleTLatex( labelOmegaEffSec, 0.85*textSizeSpectra,4,1, 42, kTRUE, 31);
            labelOmegaEffSec->Draw();

            TLatex *labelDetProcEffSec = new TLatex(0.95, 0.93-(3*textSizeSpectra*0.85),detectionProcess.Data());
            SetStyleTLatex( labelDetProcEffSec, 0.85*textSizeSpectra,4, 1, 42, kTRUE, 31);
            labelDetProcEffSec->Draw();

            canvasEffi->Update();
            canvasEffi->SaveAs(Form("%s/Omega_EffectiveSecCorrFrom%s_Weighted.%s",outputDir.Data(), nameSecOmegaPartRead[k].Data(), suffix.Data()));
            canvasEffi->SetLogy(1);
            canvasEffi->SetLeftMargin(0.09);
            canvasEffi->SetTopMargin(0.015);
            delete labelEnergyEffSec;
            delete labelOmegaEffSec;
            delete labelDetProcEffSec;
        }

    }
    //***************************************************************************************************************
    //***************************** Secondary corr factors weighted all separate ************************************
    //***************************************************************************************************************
    Int_t nSecEffis = 0;
    Int_t nSecCorrs = 0;
    for (Int_t k = 0; k<4; k++){
        if (hasSecEffi[k]) nSecEffis++;
        if (hasSecCorrFac[k]) nSecCorrs++;
    }
    if (nSecEffis > 0){
        canvasEffi->cd();
        TH2F* histo2DEffiSecOmega2 = new TH2F("histo2DEffiSecOmega","histo2DEffiSecOmega", 1000, 0., maxPtGlobalOmega,10000,minEffiSecOmega, maxEffiSecOmega);
        SetStyleHistoTH2ForGraphs(histo2DEffiSecOmega2, "#it{p}_{T} (GeV/#it{c})","#it{#varepsilon}_{sec #omega from X}",
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
        histo2DEffiSecOmega2->DrawCopy();

        TLegend* legendSecOmegaEffAll = GetAndSetLegend2(0.15, 0.94-(Int_t((nSecEffis+1)/2)*0.85*textSizeSpectra), 0.50, 0.94,32,2,"",43,0.2);
        for (Int_t k = 0; k < 4; k++){
            if (graphEfficiencySecOmegaWeighted[k]){
                DrawGammaSetMarkerTGraphAsym(graphEfficiencySecOmegaWeighted[k], markerStyleSec[k], markerSizeSec[k], colorSec[k], colorSec[k]);
                graphEfficiencySecOmegaWeighted[k]->Draw("p,e1,same");
                legendSecOmegaEffAll->AddEntry(graphEfficiencySecOmegaWeighted[k], Form("X = %s",nameSecOmegaPartLabel[k].Data()), "p");
            }
        }
        legendSecOmegaEffAll->Draw();
        TLatex *labelEnergyEffSec = new TLatex(0.95, 0.95-(0.98*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyEffSec, 0.85*textSizeSpectra,4,1, 42, kTRUE, 31);
        labelEnergyEffSec->Draw();

        TLatex *labelOmegaEffSec = new TLatex(0.95,  0.95-(2*textSizeSpectra*0.85),"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
        SetStyleTLatex( labelOmegaEffSec, 0.85*textSizeSpectra,4,1, 42, kTRUE, 31);
        labelOmegaEffSec->Draw();

        TLatex *labelDetProcEffSec = new TLatex(0.95, 0.95-(3*textSizeSpectra*0.85),detectionProcess.Data());
        SetStyleTLatex( labelDetProcEffSec, 0.85*textSizeSpectra,4, 1, 42, kTRUE, 31);
        labelDetProcEffSec->Draw();

        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/Omega_SecEfficiency_Weighted.%s",outputDir.Data(), suffix.Data()));
    }

    //***************************************************************************************************************
    //***************************************** Purity weighted *****************************************************
    //***************************************************************************************************************
    if (graphPurityOmegaWeighted){
        TCanvas* canvasPurity = new TCanvas("canvasPurity","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasPurity, 0.09, 0.017, 0.015, 0.08);
        canvasPurity->SetLogy(0);

        Double_t minPurityOmega = 0.6;
        Double_t maxPurityOmega = 1.02;

        TH2F * histo2DPurityOmega;
        histo2DPurityOmega = new TH2F("histo2DPurityOmega","histo2DPurityOmega",1000,0., maxPtGlobalOmega,10000,minPurityOmega, maxPurityOmega);
        SetStyleHistoTH2ForGraphs(histo2DPurityOmega, "#it{p}_{T} (GeV/#it{c})","#it{P}_{#omega}",
                                    0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.1);
        histo2DPurityOmega->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphPurityOmegaWeighted, 20, 1, kGray+2, kGray+2);
        graphPurityOmegaWeighted->Draw("p,e1,same");

        TLatex *labelEnergyPurityWeighted = new TLatex(0.62, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyPurityWeighted, 0.85*textSizeSpectra,4);
        labelEnergyPurityWeighted->Draw();

        TLatex *labelOmegaPurityWeighted = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
        SetStyleTLatex( labelOmegaPurityWeighted, 0.85*textSizeSpectra,4);
        labelOmegaPurityWeighted->Draw();

        TLatex *labelDetProcPurityWeighted = new TLatex(0.62, 0.15,detectionProcess.Data());
        SetStyleTLatex( labelDetProcPurityWeighted, 0.85*textSizeSpectra,4);
        labelDetProcPurityWeighted->Draw();


        canvasPurity->Update();
        canvasPurity->SaveAs(Form("%s/Omega_Purity_Weighted.%s",outputDir.Data(),suffix.Data()));

    }

    //***************************************************************************************************************
    //************************************* Acceptance weighted *****************************************************
    //***************************************************************************************************************
    if (graphAcceptanceOmegaWeighted){
        DrawGammaCanvasSettings( canvasAcc, 0.1, 0.017, 0.015, 0.08);
        canvasAcc->cd();
        canvasAcc->SetLogy(0);
        histo2DAccOmega->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphAcceptanceOmegaWeighted, 20, 1, kGray+2, kGray+2);
        graphAcceptanceOmegaWeighted->Draw("p,e1,same");

        TLatex *labelEnergyAcc = new TLatex(0.62, 0.15+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
        SetStyleTLatex( labelEnergyAcc, 0.85*textSizeSpectra,4);
        labelEnergyAcc->Draw();

        TLatex *labelOmegaAcc = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
        SetStyleTLatex( labelOmegaAcc, 0.85*textSizeSpectra,4);
        labelOmegaAcc->Draw();

        TLatex *labelDetProcAcc = new TLatex(0.62, 0.15,detectionProcess.Data());
        SetStyleTLatex( labelDetProcAcc, 0.85*textSizeSpectra,4);
        labelDetProcAcc->Draw();

        canvasAcc->Update();
        canvasAcc->SaveAs(Form("%s/Omega_Acceptance_weighted_%i.%s",outputDir.Data(),mode,suffix.Data()));
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
    if (mode == 0){
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

    if(optionEnergy.CompareTo("pPb_5.023TeV")==0){
      if(mode == 2){
        minCorrYield       = 1e-8;
        maxCorrYield       = 1;
      }else if(mode == 4){
        minCorrYield       = 1e-8;
        maxCorrYield       = 0.2;
      }else if(mode == 3){
        minCorrYield       = 1e-7;
        maxCorrYield       = 10;
      }
    }



    TH2F * histo2DInvYieldScaled;
    histo2DInvYieldScaled = new TH2F("histo2DInvYieldScaled","histo2DInvYieldScaled",1000,0., maxPtGlobalOmega,10000,minCorrYield,maxCorrYield);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldScaled, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.8,1.55);
    histo2DInvYieldScaled->DrawCopy();

    Double_t factorOmega              = 10000;
    if (mode == 2) factorOmega        = 10000;
    if (mode == 4) factorOmega        = 5000;

    // two component model fit
    Double_t paramTCM[5] = {graphCorrectedYieldWeightedAverageOmegaStat->GetY()[0],0.3,graphCorrectedYieldWeightedAverageOmegaStat->GetY()[0]/factorOmega,0.8,3};
    if(mode == 4 && optionEnergy.CompareTo("8TeV")==0){
      paramTCM[0]=0.008; paramTCM[1]=0.65; paramTCM[2]=1.72; paramTCM[3]=0.47; paramTCM[4]=2.93;
    }
    TF1* fitInvYieldOmega = FitObject("tcm","fitInvYieldOmega","Omega",graphCorrectedYieldWeightedAverageOmegaStat,minPtGlobalOmega,maxPtGlobalOmega,paramTCM,"QNRMEX0+");

    // Tsallis fit
//     Double_t paramGraph[3]                  = {1000, 8., 0.13};
//     TF1* fitInvYieldOmega                     = FitObject("l","fitInvYieldOmega","Omega",graphCorrectedYieldWeightedAverageOmegaStat,minPtGlobalOmega,maxPtGlobalOmega,paramGraph,"QNRME+");

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageOmegaSys, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
    graphCorrectedYieldWeightedAverageOmegaSys->Draw("p,E2,same");

    TLegend* legendScaled = GetAndSetLegend2(0.72, 0.95-(1.15*nrOfTrigToBeCombOmegaRed*0.85*textSizeSpectra), 0.95, 0.95,32);

    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if (graphsCorrectedYieldSysShrunkOmega[i]) DrawGammaSetMarkerTGraphAsym(graphsCorrectedYieldSysShrunkOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
        if (graphsCorrectedYieldSysShrunkOmega[i])DrawGammaSetMarker(histoCorrectedYieldOmegaScaled[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        if ( !maskedFullyOmega[i] )histoCorrectedYieldOmegaScaled[i]->DrawCopy("e1,same");
        if ( !maskedFullyOmega[i] )graphsCorrectedYieldSysShrunkOmega[i]->Draw("p,E2,same");
        if (graphsCorrectedYieldSysShrunkOmega[i] && !maskedFullyOmega[i])legendScaled->AddEntry(histoCorrectedYieldOmegaScaled[i],triggerNameLabel[i].Data(),"p");
    }

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageOmegaStat, 24, 2, kRed , kRed, 1, kTRUE);
    graphCorrectedYieldWeightedAverageOmegaStat->Draw("p,E,same");
    legendScaled->AddEntry(graphCorrectedYieldWeightedAverageOmegaStat,"Final","p");
    legendScaled->Draw();

    DrawGammaSetMarkerTF1( fitInvYieldOmega, 7, 2, kGray+2);
    fitInvYieldOmega->Draw("same");

    labelEnergyUnscaled->Draw();
    labelOmegaUnscaled->Draw();
    labelDetProcUnscaled->Draw();

    canvasCorrScaled->Update();
    canvasCorrScaled->SaveAs(Form("%s/Omega_%s_CorrectedYieldScaledTrigg_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

    //***************************************************************************************************************
    //************************************Plotting final invariant yield ********************************************
    //***************************************************************************************************************

    histo2DInvYieldScaled->DrawCopy();
    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageOmegaSys, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
    graphCorrectedYieldWeightedAverageOmegaSys->Draw("p,E2,same");

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldWeightedAverageOmegaStat, 24, 2, kBlack , kBlack, 1, kTRUE);
    graphCorrectedYieldWeightedAverageOmegaStat->Draw("p,E,same");

    fitInvYieldOmega->Draw("same");

    labelEnergyUnscaled->Draw();
    labelOmegaUnscaled->Draw();
    labelDetProcUnscaled->Draw();

    canvasCorrScaled->Update();
    canvasCorrScaled->SaveAs(Form("%s/Omega_%s_CorrectedYieldFinal_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

    if (fitBinShiftOmegaTCM){
        DrawGammaSetMarkerTF1( fitBinShiftOmegaTCM, 9, 2, kRed+2);
        fitBinShiftOmegaTCM->Draw("same");
        canvasCorrScaled->Update();
        canvasCorrScaled->SaveAs(Form("%s/Omega_%s_CorrectedYieldFinalWithBinShiftFit_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
    }


    histo2DInvYieldScaled->DrawCopy();
    graphCorrectedYieldWeightedAverageOmegaSys->Draw("p,E2,same");
    graphCorrectedYieldWeightedAverageOmegaStat->Draw("p,E,same");

    fitInvYieldOmega->Draw("same");

    DrawGammaSetMarker(histoMCInputOmega[0],  0, 0, kBlue+2, kBlue+2);
    histoMCInputOmega[0]->Draw("same,hist,c");

    labelEnergyUnscaled->Draw();
    labelOmegaUnscaled->Draw();
    labelDetProcUnscaled->Draw();

    canvasCorrScaled->Update();
    canvasCorrScaled->SaveAs(Form("%s/Omega_%s_CorrectedYieldFinal_WithMC_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

    if (!enableEta) delete canvasCorrScaled;

    //***************************************************************************************************************
    //****************************** Ratio to fit for individual spectra full range *********************************
    //***************************************************************************************************************
    TCanvas* canvasRatioSpec = new TCanvas("canvasRatioSpec","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasRatioSpec, 0.09, 0.017, 0.015, 0.08);
    canvasRatioSpec->SetLogy(0);

    TH2F * histo2DRatioToFitOmega;
    histo2DRatioToFitOmega = new TH2F("histo2DRatioToFitOmega","histo2DRatioToFitOmega",1000,0., maxPtGlobalOmega,1000,0.55, 1.85);
    SetStyleHistoTH2ForGraphs(histo2DRatioToFitOmega, "#it{p}_{T} (GeV/#it{c})","Data/Fit",
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
    histo2DRatioToFitOmega->DrawCopy();

    TH1D* histoCorrectedYieldToFitOmega[MaxNumberOfFiles];
    TGraphAsymmErrors* graphCorrectedYieldToFitOmega[MaxNumberOfFiles];
    TLegend* legendRatioSpecOmega = GetAndSetLegend2(0.12, 0.95-(1.05*nrOfTrigToBeCombOmegaRed/2*0.85*textSizeSpectra), 0.5, 0.95,28);
    legendRatioSpecOmega->SetNColumns(2);
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if (graphsCorrectedYieldSysRemoved0Omega[i] && !maskedFullyOmega[i]){
          histoCorrectedYieldToFitOmega[i] = CalculateHistoRatioToFit (histoCorrectedYieldOmegaScaled[i], fitInvYieldOmega);
          DrawGammaSetMarker(histoCorrectedYieldToFitOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
          legendRatioSpecOmega->AddEntry(histoCorrectedYieldToFitOmega[i],triggerNameLabel[i].Data(),"p");
          graphCorrectedYieldToFitOmega[i] = CalculateGraphErrRatioToFit(graphsCorrectedYieldSysRemoved0Omega[i], fitInvYieldOmega);
          DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
          graphCorrectedYieldToFitOmega[i]->Draw("p,E2,same");
        }

        DrawGammaLines(0., maxPtGlobalOmega , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0., maxPtGlobalOmega , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0., maxPtGlobalOmega , 0.9, 0.9,0.1, kGray, 7);

        if (graphsCorrectedYieldSysRemoved0Omega[i]){
          histoCorrectedYieldToFitOmega[i]->DrawCopy("e1,same");
        }
    }
    legendRatioSpecOmega->Draw();

    TLatex *labelEnergyRatio = new TLatex(0.6, 0.93, collisionSystem.Data());
    SetStyleTLatex( labelEnergyRatio, 0.85*textSizeSpectra,4);
    labelEnergyRatio->Draw();

    TLatex *labelOmegaRatio = new TLatex(0.6, 0.93-textSizeSpectra*0.85*1.04, "#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
    SetStyleTLatex( labelOmegaRatio, 0.85*textSizeSpectra,4);
    labelOmegaRatio->Draw();

    TLatex *labelDetProcRatio = new TLatex(0.6, 0.93-(2*textSizeSpectra*0.85), detectionProcess.Data());
    SetStyleTLatex( labelDetProcRatio, 0.85*textSizeSpectra,4);
    labelDetProcRatio->Draw();

    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/Omega_%s_RatioSpectraToFit.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

    //***************************************************************************************************************
    //****************************** Ratio to fit for individual spectra used range *********************************
    //***************************************************************************************************************

    histo2DRatioToFitOmega->DrawCopy();

    TGraphAsymmErrors* graphCorrectedYieldToFitOmegaUsed[MaxNumberOfFiles];
    TGraphAsymmErrors* graphCorrectedYieldToFitOmegaSysUsed[MaxNumberOfFiles];
    for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
        if (!maskedFullyOmega[i]) {
            if (graphsCorrectedYieldShrunkOmega[i]) graphCorrectedYieldToFitOmegaUsed[i] = CalculateGraphErrRatioToFit (graphsCorrectedYieldShrunkOmega[i], fitInvYieldOmega);
            if (graphsCorrectedYieldSysShrunkOmega[i]) graphCorrectedYieldToFitOmegaSysUsed[i] = CalculateGraphErrRatioToFit(graphsCorrectedYieldSysShrunkOmega[i], fitInvYieldOmega);

            if (graphsCorrectedYieldSysShrunkOmega[i])DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitOmegaSysUsed[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
            if (graphsCorrectedYieldShrunkOmega[i])DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldToFitOmegaUsed[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);

            if (graphsCorrectedYieldSysShrunkOmega[i])graphCorrectedYieldToFitOmegaSysUsed[i]->Draw("p,E2,same");

            DrawGammaLines(0., maxPtGlobalOmega , 1., 1.,0.1, kGray+2);
            DrawGammaLines(0., maxPtGlobalOmega , 1.1, 1.1,0.1, kGray, 7);
            DrawGammaLines(0., maxPtGlobalOmega , 0.9, 0.9,0.1, kGray, 7);

            if (graphsCorrectedYieldShrunkOmega[i])graphCorrectedYieldToFitOmegaUsed[i]->Draw("e1,same");
        }
    }

    legendRatioSpecOmega->Draw();

    labelEnergyRatio->Draw();
    labelOmegaRatio->Draw();
    labelDetProcRatio->Draw();

    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/Omega_%s_RatioSpectraToFitUsed.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

    //***************************************************************************************************************
    //****************************** Ratio to fit for final spectrum ************************************************
    //***************************************************************************************************************

    histo2DRatioToFitOmega->DrawCopy();

    TGraphAsymmErrors* graphCorrectedYieldFinalStatToFitOmega;
    TGraphAsymmErrors* graphCorrectedYieldFinalSysToFitOmega;
    TH1D* histoMCInputToFit;

    graphCorrectedYieldFinalStatToFitOmega    = CalculateGraphErrRatioToFit (graphCorrectedYieldWeightedAverageOmegaStat, fitInvYieldOmega);
    graphCorrectedYieldFinalSysToFitOmega     = CalculateGraphErrRatioToFit(graphCorrectedYieldWeightedAverageOmegaSys, fitInvYieldOmega);
    histoMCInputToFit                       = (TH1D*)histoMCInputOmega[0]->Clone("OmegaMCToFit");
    histoMCInputToFit                       = CalculateHistoRatioToFitNLO (histoMCInputToFit, fitInvYieldOmega, minPtGlobalOmega);

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldFinalSysToFitOmega, 24, 2, kGray+1 , kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldFinalStatToFitOmega, 24, 2, kBlack , kBlack, 1, kTRUE);

    graphCorrectedYieldFinalSysToFitOmega->Draw("p,E2,same");

    DrawGammaLines(0., maxPtGlobalOmega , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0., maxPtGlobalOmega , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(0., maxPtGlobalOmega , 0.9, 0.9,0.1, kGray, 7);

    graphCorrectedYieldFinalStatToFitOmega->Draw("p,E,same");

    labelEnergyRatio->Draw();
    labelOmegaRatio->Draw();
    labelDetProcRatio->Draw();

    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/Omega_%s_RatioSpectraToFitFinal.%s",outputDir.Data(),isMC.Data(), suffix.Data()));

    TH2F * histo2DRatioToFitOmega2;
    histo2DRatioToFitOmega2 = new TH2F("histo2DRatioToFitOmega2","histo2DRatioToFitOmega2",1000,0., maxPtGlobalOmega,1000,0.25, 1.85);
    SetStyleHistoTH2ForGraphs(histo2DRatioToFitOmega2, "#it{p}_{T} (GeV/#it{c})","Data/Fit",
                            0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
    histo2DRatioToFitOmega2->DrawCopy();

    graphCorrectedYieldFinalSysToFitOmega->Draw("p,E2,same");

    DrawGammaLines(0., maxPtGlobalOmega , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0., maxPtGlobalOmega , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(0., maxPtGlobalOmega , 0.9, 0.9,0.1, kGray, 7);


    graphCorrectedYieldFinalStatToFitOmega->Draw("p,E,same");
    DrawGammaSetMarker(histoMCInputToFit, 20, 1.2, kBlue+1, kBlue+1);
    histoMCInputToFit->Draw("same,pe");
    DrawGammaLines(0., maxPtGlobalOmega , 0.6, 0.6,1, kBlue-7, 7);

    labelEnergyRatio->Draw();
    labelOmegaRatio->Draw();
    labelDetProcRatio->Draw();

    canvasRatioSpec->Update();
    canvasRatioSpec->SaveAs(Form("%s/Omega_%s_RatioSpectraToFitFinal_withMC.%s",outputDir.Data(),isMC.Data(), suffix.Data()));


    if (!enableEta) delete canvasRatioSpec;



    //***************************************************************************************************************
    //************************************Loading eta histograms ****************************************************
    //***************************************************************************************************************
    TString FileNameCorrectedEta        [MaxNumberOfFiles];
    TFile* fileCorrectedEta             [MaxNumberOfFiles];
    TString FileNameCorrectedOmegaEtaBin  [MaxNumberOfFiles];
    TFile* fileCorrectedOmegaEtaBin       [MaxNumberOfFiles];
    Bool_t foundOmegaEtaBinFile           [MaxNumberOfFiles];
    Bool_t doEtaToOmega                                               = kFALSE;

    TString FileNameUnCorrectedEta      [MaxNumberOfFiles];
    TFile* fileUnCorrectedEta           [MaxNumberOfFiles];

    TH1D*   histoCorrectedYieldEta      [MaxNumberOfFiles];
    TH1D*   histoEfficiencyEta          [MaxNumberOfFiles];
    TH1D*   histoAcceptanceEta          [MaxNumberOfFiles];
    TH1D*   histoAcceptanceEtaWOEvtWeights [MaxNumberOfFiles];
    TH1D*   histoEffTimesAccEta         [MaxNumberOfFiles];
    TH1D*   histoRawYieldEta            [MaxNumberOfFiles];
    TH1D*   histoCorrectedYieldOmegaEtaBin[MaxNumberOfFiles];
    TH1D*   histoEtaToOmega               [MaxNumberOfFiles];
    TH1D*   histoMassEtaData            [MaxNumberOfFiles];
    TH1D*   histoMassEtaMC              [MaxNumberOfFiles];
    TH1D*   histoWidthEtaData           [MaxNumberOfFiles];
    TH1D*   histoWidthEtaMC             [MaxNumberOfFiles];
    TH1D*   histoEtaInvMassSigPlusBG    [MaxNumberOfFiles];
    TH1D*   histoEtaInvMassSig          [MaxNumberOfFiles];
    TH1D*   histoEtaInvMassBG           [MaxNumberOfFiles];
    TF1*    fitEtaInvMassSig            [MaxNumberOfFiles];

    TH1D* histoCorrectedYieldOmegaEtaBinBinShift[MaxNumberOfFiles];
    TH1D* histoCorrectedYieldEtaBinShift[MaxNumberOfFiles];

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
    TGraphAsymmErrors* graphEtaToOmegaWeightedAverageStat             = NULL;
    TGraphAsymmErrors* graphEtaToOmegaWeightedAverageSys              = NULL;
    // create pointers for shrunk graphs and supporting figures
    TGraphAsymmErrors* graphsCorrectedYieldRemoved0Eta      [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsCorrectedYieldSysRemoved0Eta   [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsEtaToOmegaRemoved0               [MaxNumberOfFiles];
    TGraphAsymmErrors* graphsEtaToOmegaSysRemoved0            [MaxNumberOfFiles];
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

    Int_t nRelSysErrEtaSources          = 0;
    TGraphAsymmErrors* graphRelSysErrEtaSource              [30][MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedRelSysErrEtaSource       [30][MaxNumberOfFiles];
    TGraphAsymmErrors* graphRelSysErrEtaSourceWeighted      [30];
    Int_t nRelSysErrEtaToOmegaSources     = 0;
    TGraphAsymmErrors* graphRelSysErrEtaToOmegaSource         [30][MaxNumberOfFiles];
    TGraphAsymmErrors* graphOrderedRelSysErrEtaToOmegaSource  [30][MaxNumberOfFiles];
    TGraphAsymmErrors* graphRelSysErrEtaToOmegaSourceWeighted [30];
    for (Int_t j = 0; j< 30; j++){
        graphRelSysErrEtaSourceWeighted[j]          = NULL;
        graphRelSysErrEtaToOmegaSourceWeighted[j]     = NULL;
        for (Int_t k = 0; k< MaxNumberOfFiles; k++){
            graphRelSysErrEtaSource[j][k]               = NULL;
            graphOrderedRelSysErrEtaSource[j][k]        = NULL;
            graphRelSysErrEtaToOmegaSource[j][k]          = NULL;
            graphOrderedRelSysErrEtaToOmegaSource[j][k]   = NULL;
        }
    }
    Bool_t sysAvailSingleEtaToOmega                           [MaxNumberOfFiles];
    Bool_t sysAvailSingleEta                                [MaxNumberOfFiles];

    if ( mode == 4 && optionEnergy.CompareTo("pPb_5.023TeV") == 0 ){
        nameCorrectedYield                              = "CorrectedYieldTrueEff";
        nameEfficiency                                  = "TrueMesonEffiPt";
        nameMassMC                                      = "histoTrueMassMeson";
        nameWidthMC                                     = "histoTrueFWHMMeson";
    }
    cout << "reading the following Eta hists: " << nameCorrectedYield.Data() << "\t" << nameEfficiency.Data() << "\t" << nameMassMC.Data() << "\t" << nameWidthMC.Data() << endl;

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
            Double_t scaleFacSinBin                         = 1.0;
            Int_t exampleBin                                = ReturnSingleInvariantMassBinPlotting ("Eta", optionEnergy, mode, trigger.Atoi(), scaleFacSinBin);

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
            histoRawYieldEta[i]                             = (TH1D*)fileUnCorrectedEta[i]->Get(Form("histoYieldMesonPerEvent%s",InvMassTypeEnding.Data()));
            histoRawYieldEta[i]->SetName(Form("RAWYieldPerEvent_%s",cutNumber[i].Data()));

            histoMassEtaData[i]                                 = (TH1D*)fileCorrectedEta[i]->Get(Form("histoMassMeson%s",InvMassTypeEnding.Data()));
            histoMassEtaData[i]->SetName(Form("Eta_Mass_data_%s",cutNumber[i].Data()));
            histoMassEtaMC[i]                                   = (TH1D*)fileCorrectedEta[i]->Get(nameMassMC.Data());
            histoMassEtaMC[i]->SetName(Form("Eta_Mass_MC_%s",cutNumber[i].Data()));
            histoWidthEtaData[i]                                = (TH1D*)fileCorrectedEta[i]->Get(Form("histoFWHMMeson%s",InvMassTypeEnding.Data()));
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

//                    //limit trigger efficiency to 1
                    if(optionEnergy.CompareTo("8TeV") == 0){
                      for(Int_t j = 1; j<histoTriggerEffEta[i]->GetNbinsX()+1; j++){
                        Double_t binC = histoTriggerEffEta[i]->GetBinContent(j);
                        if(binC > 1.){
                          histoEffBaseEta[i]->SetBinContent(j, histoEffBaseEta[i]->GetBinContent(j)*binC);
                          histoTriggerEffEta[i]->SetBinContent(j, 1.);
                        }
                      }
                    }

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

            // read files for eta to Omega ratio
            FileNameCorrectedOmegaEtaBin[i]                   = Form("%s/%s/OmegaEtaBinning_%s_GammaConvV1Correction_%s.root", cutNumber[i].Data(), optionEnergy.Data(), isMC.Data(),
                                                                        cutNumber[i].Data());
            cout<< FileNameCorrectedOmegaEtaBin[i] << endl;
            fileCorrectedOmegaEtaBin[i]                       = new TFile(FileNameCorrectedOmegaEtaBin[i]);
            if (fileCorrectedOmegaEtaBin[i]->IsZombie()) {
                foundOmegaEtaBinFile[i]                        = kFALSE;
            } else {
                foundOmegaEtaBinFile[i]                        = kTRUE;
                doEtaToOmega                                   = kTRUE;
            }

            if (foundOmegaEtaBinFile[i]){

                if (! doBinShiftForEtaToOmega){
                    histoCorrectedYieldOmegaEtaBin[i]             = (TH1D*)fileCorrectedOmegaEtaBin[i]->Get(nameCorrectedYield.Data());
                    histoCorrectedYieldOmegaEtaBin[i]->SetName(Form("CorrectedYieldOmegaEtaBin_%s",cutNumber[i].Data()));
                    if(optionEnergy.CompareTo("8TeV")==0 && mode==4){
                      histoEtaToOmega[i]                            = (TH1D*)histoCorrectedYieldOmegaEtaBin[i]->Clone(Form("EtaToOmega_%s", cutNumber[i].Data()));
                      for(Int_t iB=1; iB<=histoCorrectedYieldOmegaEtaBin[i]->GetNbinsX(); iB++){histoEtaToOmega[i]->SetBinContent(iB,histoCorrectedYieldEta[i]->GetBinContent(iB));}
                      histoEtaToOmega[i]->Divide(histoEtaToOmega[i],histoCorrectedYieldOmegaEtaBin[i],1.,1.,"");
                    }else{
                      histoEtaToOmega[i]                            = (TH1D*)histoCorrectedYieldEta[i]->Clone(Form("EtaToOmega_%s", cutNumber[i].Data()));
                      histoEtaToOmega[i]->Divide(histoEtaToOmega[i],histoCorrectedYieldOmegaEtaBin[i],1.,1.,"");
                    }
                } else {
                    cout << fitBinShiftOmega << " - " << fitBinShiftEta << endl;
                    histoCorrectedYieldOmegaEtaBin[i]             = (TH1D*)fileCorrectedOmegaEtaBin[i]->Get(nameCorrectedYield.Data());
                    histoCorrectedYieldOmegaEtaBin[i]->SetName(Form("CorrectedYieldOmegaEtaBin_%s",cutNumber[i].Data()));
                    cout << "shifting Omega in eta binning: " <<  cutNumber[i].Data() << endl;
                    histoCorrectedYieldOmegaEtaBinBinShift[i]     = (TH1D*)histoCorrectedYieldOmegaEtaBin[i]->Clone(Form("CorrectedYieldOmegaEtaBinBinShifted_%s",cutNumber[i].Data()));
                    histoCorrectedYieldOmegaEtaBinBinShift[i]     = ApplyYshiftIndividualSpectra( histoCorrectedYieldOmegaEtaBinBinShift[i], fitBinShiftOmega);
                    cout << "shifting eta: " <<  cutNumber[i].Data() << endl;
                    histoCorrectedYieldEtaBinShift[i]           = (TH1D*)histoCorrectedYieldEta[i]->Clone(Form("CorrectedYieldEtaBinShifted_%s",cutNumber[i].Data()));
                    histoCorrectedYieldEtaBinShift[i]           = ApplyYshiftIndividualSpectra( histoCorrectedYieldEtaBinShift[i], fitBinShiftEta);

                    if(optionEnergy.CompareTo("8TeV")==0 && mode==4){
                      histoEtaToOmega[i]                            = (TH1D*)histoCorrectedYieldOmegaEtaBinBinShift[i]->Clone(Form("EtaToOmega_%s", cutNumber[i].Data()));
                      for(Int_t iB=1; iB<=histoCorrectedYieldOmegaEtaBinBinShift[i]->GetNbinsX(); iB++){histoEtaToOmega[i]->SetBinContent(iB,histoCorrectedYieldEtaBinShift[i]->GetBinContent(iB));}
                      histoEtaToOmega[i]->Divide(histoEtaToOmega[i],histoCorrectedYieldOmegaEtaBinBinShift[i],1.,1.,"");
                    }else{
                      histoEtaToOmega[i]                            = (TH1D*)histoCorrectedYieldEtaBinShift[i]->Clone(Form("EtaToOmega%s_%s", addNameBinshift.Data(), cutNumber[i].Data()));
                      histoEtaToOmega[i]->Divide(histoEtaToOmega[i],histoCorrectedYieldOmegaEtaBinBinShift[i],1.,1.,"");
                    }
                }
            } else {
                doBinShiftForEtaToOmega   = kFALSE;
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
        //************************************Plotting binshift corrections *********************************************
        //***************************************************************************************************************

        if(doBinShiftForEtaToOmega){
          canvasEffi->cd();
          canvasEffi->SetLeftMargin(0.1);
          canvasEffi->SetBottomMargin(0.1);
          canvasEffi->SetLogy(0);
          canvasEffi->SetLogx(1);
          TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0., 100.);
          SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (Y) / no shift",
                                  0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.1, 1.2);
          histoBinShift->GetXaxis()->SetRangeUser(minPtGlobalOmega,maxPtGlobalOmega);
          histoBinShift->GetXaxis()->SetMoreLogLabels();
          histoBinShift->GetYaxis()->SetRangeUser(0.8,1.05);
          histoBinShift->DrawCopy();

          TLegend* legendBinShift = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
          legendBinShift->SetNColumns(2);
          for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
              histoCorrectedYieldOmegaEtaBinBinShift[i]->Divide(histoCorrectedYieldOmegaEtaBinBinShift[i],histoCorrectedYieldOmegaEtaBin[i],1.,1.,"B");
              DrawGammaSetMarker(histoCorrectedYieldOmegaEtaBinBinShift[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
              histoCorrectedYieldOmegaEtaBinBinShift[i]->DrawCopy("hist p same");
              legendBinShift->AddEntry(histoCorrectedYieldOmegaEtaBinBinShift[i],triggerNameLabel[i].Data(),"p");
          }
          legendBinShift->Draw();

          labelEnergyEffi->Draw();
          TLatex *labelBinShiftOmega = new TLatex(0.62, maxYLegendEffi+0.02+0.99*textSizeSpectra*0.85,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
          SetStyleTLatex( labelBinShiftOmega, 0.85*textSizeSpectra,4);
          labelBinShiftOmega->Draw();
          labelDetProcEffi->Draw();

          canvasEffi->RedrawAxis();
          canvasEffi->Update();
          canvasEffi->SaveAs(Form("%s/OmegaEtaBinning_%s_BinShiftCorrection.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

          histoBinShift->GetXaxis()->SetRangeUser(minPtGlobalEta,maxPtGlobalEta);
          histoBinShift->DrawCopy();

          TLegend* legendBinShift2 = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
          legendBinShift2->SetNColumns(2);
          for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
              histoCorrectedYieldEtaBinShift[i]->Divide(histoCorrectedYieldEtaBinShift[i],histoCorrectedYieldEta[i],1.,1.,"B");
              DrawGammaSetMarker(histoCorrectedYieldEtaBinShift[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
              histoCorrectedYieldEtaBinShift[i]->DrawCopy("hist p same");
              legendBinShift2->AddEntry(histoCorrectedYieldEtaBinShift[i],triggerNameLabel[i].Data(),"p");
          }
          legendBinShift2->Draw();

          labelEnergyEffi->Draw();
          TLatex *labelBinShiftEta = new TLatex(0.62, maxYLegendEffi+0.02+0.99*textSizeSpectra*0.85,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
          SetStyleTLatex( labelBinShiftEta, 0.85*textSizeSpectra,4);
          labelBinShiftEta->Draw();
          labelDetProcEffi->Draw();

          canvasEffi->RedrawAxis();
          canvasEffi->Update();
          canvasEffi->SaveAs(Form("%s/Eta_%s_BinShiftCorrection.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

          histoBinShift->GetXaxis()->SetRangeUser(minPtGlobalEta,maxPtGlobalEta);
          if(optionEnergy.CompareTo("8TeV")==0 && mode==4)
              histoBinShift->GetXaxis()->SetRangeUser(minPtGlobalEta,20.);
          histoBinShift->GetYaxis()->SetRangeUser(0.95,1.05);
          if(!optionEnergy.CompareTo("8TeV") && mode==0)
              histoBinShift->GetYaxis()->SetRangeUser(0.8,1.2);
          histoBinShift->DrawCopy();

          TLegend* legendBinShift3 = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
          legendBinShift3->SetNColumns(2);
          for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
              if(optionEnergy.CompareTo("8TeV")==0 && mode==4){
                TH1D* histoCorrectedYieldEtaBinShiftTEMP = (TH1D*)histoCorrectedYieldOmegaEtaBinBinShift[i]->Clone(Form("OmegaEtaBinning_%i", i));
                for(Int_t iB=1; iB<=histoCorrectedYieldOmegaEtaBinBinShift[i]->GetNbinsX(); iB++){histoCorrectedYieldEtaBinShiftTEMP->SetBinContent(iB,histoCorrectedYieldEtaBinShift[i]->GetBinContent(iB));}
                histoCorrectedYieldEtaBinShiftTEMP->Divide(histoCorrectedYieldEtaBinShiftTEMP,histoCorrectedYieldOmegaEtaBinBinShift[i],1.,1.,"B");
                histoCorrectedYieldEtaBinShift[i] = histoCorrectedYieldEtaBinShiftTEMP;
              }else{
                histoCorrectedYieldEtaBinShift[i]->Divide(histoCorrectedYieldEtaBinShift[i],histoCorrectedYieldOmegaEtaBinBinShift[i],1.,1.,"B");
              }
              DrawGammaSetMarker(histoCorrectedYieldEtaBinShift[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
              histoCorrectedYieldEtaBinShift[i]->DrawCopy("hist p same");
              legendBinShift3->AddEntry(histoCorrectedYieldEtaBinShift[i],triggerNameLabel[i].Data(),"p");
          }
          legendBinShift3->Draw();

          labelEnergyEffi->Draw();
          TLatex *labelBinShiftEtaToOmega = new TLatex(0.62, maxYLegendEffi+0.02+0.99*textSizeSpectra*0.85,"#eta/#omega");
          SetStyleTLatex( labelBinShiftEtaToOmega, 0.85*textSizeSpectra,4);
          labelBinShiftEtaToOmega->Draw();
          labelDetProcEffi->Draw();

          canvasEffi->RedrawAxis();
          canvasEffi->Update();
          canvasEffi->SaveAs(Form("%s/EtaToOmega_%s_BinShiftCorrection.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
          canvasEffi->SetLogy(1);
          canvasEffi->SetLogx(0);
          canvasEffi->SetLeftMargin(0.09);
          canvasEffi->SetBottomMargin(0.08);
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
        SetStyleHistoTH2ForGraphs(histo2DEffiEta, "#it{p}_{T} (GeV/#it{c})","#it{#varepsilon}_{#eta}#upoint#it{#kappa}_{trigg}",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.05);
        histo2DEffiEta->DrawCopy();
        histo2DEffiEta->SetYTitle("#it{#varepsilon}_{#eta}");

        TLegend* legendEffiEta = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
        legendEffiEta->SetNColumns(2);
        for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
            DrawGammaSetMarker(histoEfficiencyEta[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
            histoEfficiencyEta[i]->DrawCopy("e1,same");
            legendEffiEta->AddEntry(histoEfficiencyEta[i],triggerNameLabel[i].Data(),"p");
        }
        legendEffiEta->Draw();

        labelEnergyEffi->Draw();
        TLatex *labelEtaEffi = new TLatex(0.62, maxYLegendEffi+0.02+0.99*textSizeSpectra*0.85,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
        SetStyleTLatex( labelEtaEffi, 0.85*textSizeSpectra,4);
        labelEtaEffi->Draw();
        labelDetProcEffi->Draw();

        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/Eta_Efficiency_%i.%s",outputDir.Data(),mode,suffix.Data()));

        //***************************************************************************************************************
        //************************************ Plotting trigger efficiencies Eta ****************************************
        //***************************************************************************************************************
        if (enableTriggerEffEtaAll){
            TCanvas* canvasTriggerEffi = new TCanvas("canvasTriggerEffi","",0,0,1000,900);// gives the page size
            DrawGammaCanvasSettings( canvasTriggerEffi, 0.09, 0.017, 0.015, 0.08);
            canvasTriggerEffi->SetLogy(0);

            Double_t minEffiTrigEta         = 0;
            Double_t maxEffiTrigEta         = 1.1;
            //if(optionEnergy.CompareTo("8TeV") == 0 && mode == 2) maxEffiTrigEta = 1.5;

            TH2F * histo2DTriggerEffiEta;
            histo2DTriggerEffiEta = new TH2F("histo2DTriggerEffiEta","histo2DTriggerEffiEta",1000,0., maxPtGlobalEta,10000,minEffiTrigEta, maxEffiTrigEta);
            SetStyleHistoTH2ForGraphs(histo2DTriggerEffiEta, "#it{p}_{T} (GeV/#it{c})","#it{#kappa}_{trigg, #eta}",
                                        0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.05);
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
            canvasEffi->SaveAs(Form("%s/Eta_EfficiencyW0TriggEff_%i.%s",outputDir.Data(),mode,suffix.Data()));
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
             (optionEnergy.CompareTo("8TeV")==0) ||
             (optionEnergy.CompareTo("pPb_5.023TeV")==0)
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

        TLatex *labelEtaMass = new TLatex(0.14, 0.85+0.99*textSizeSpectra*0.85,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
        SetStyleTLatex( labelEtaMass, 0.85*textSizeSpectra,4);
        labelEtaMass->Draw();
        labelDetProcMass->Draw();

        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/Eta_%s_Mass1_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

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
            canvasMass->SaveAs(Form("%s/Eta_%s_Mass2_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
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
             (optionEnergy.CompareTo("8TeV")==0) ||
             (optionEnergy.CompareTo("pPb_5.023TeV")==0)
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

        TLatex *labelEtaWidth = new TLatex(0.14, 0.84+0.99*textSizeSpectra*0.85,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
        SetStyleTLatex( labelEtaWidth, 0.85*textSizeSpectra,4);
        labelEtaWidth->Draw();
        labelDetProcWidth->Draw();

        canvasWidth->Update();
        canvasWidth->SaveAs(Form("%s/Eta_%s_Width1_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

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
            canvasWidth->SaveAs(Form("%s/Eta_%s_Width2_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
        }
        //***************************************************************************************************************
        //************************************Plotting unscaled invariant raw-yield Eta *********************************
        //***************************************************************************************************************
        canvasRawUnscaled->cd();

        Double_t minCorrYieldRawUnscaledEta     = 7e-8;
        Double_t maxCorrYieldRawUnscaledEta     = 4e-3;
        if (mode == 4) {
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
        if(optionEnergy.CompareTo("pPb_5.023TeV")==0){
            if(mode == 2){
                minCorrYieldRawUnscaledEta    = 2e-8;
                maxCorrYieldRawUnscaledEta    = 8e-4;
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
        TLatex *labelEtaRawUnscaled     = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
        SetStyleTLatex( labelEtaRawUnscaled, 0.85*textSizeSpectra,4);
        labelEtaRawUnscaled->Draw();
        labelDetProcRawUnscaled->Draw();

        canvasRawUnscaled->Update();
        canvasRawUnscaled->SaveAs(Form("%s/Eta_%s_RawYieldUnscaledTrigg_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
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

        if(optionEnergy.CompareTo("pPb_5.023TeV")==0){
            if(mode == 2){
                minCorrYieldUnscaledEta         = 1e-8;
                maxCorrYieldUnscaledEta         = 1e-1;
            }else if(mode == 4){
                minCorrYieldUnscaledEta         = 9e-9;
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
        TLatex *labelEtaUnscaled = new TLatex(0.2, 0.12+textSizeSpectra*0.85*0.75,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
        SetStyleTLatex( labelEtaUnscaled, 0.85*textSizeSpectra,4);
        labelEtaUnscaled->Draw();
        labelDetProcUnscaled->Draw();


        canvasCorrUnscaled->Update();
        canvasCorrUnscaled->SaveAs(Form("%s/Eta_%s_CorrectedYieldUnscaledTrigg_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
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
            fileFitsOutput << triggerName[i].Data() << endl;
            histoCorrectedYieldEtaScaled[i] = (TH1D*)histoCorrectedYieldEta[i]->Clone(Form("CorrectedYieldEtaScaled_%s", triggerName[i].Data()));
            histoCorrectedYieldEtaScaled[i]->Sumw2();
            histoCorrectedYieldEtaScaled[i]->Scale(1./triggRejecFac[i][trigSteps[i][0]]);
            if (trigSteps[i][1]!= trigSteps[i][0]){
                fileFitsOutput << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << endl;
                histoCorrectedYieldEtaScaled[i]->Scale(1./triggRejecFac[trigSteps[i][0]][trigSteps[i][1]]);
            }
            if (trigSteps[i][2]!= trigSteps[i][1]){
                fileFitsOutput << triggRejecFac[i][trigSteps[i][0]] << "\t" << triggRejecFac[trigSteps[i][0]][trigSteps[i][1]] << "\t"<< triggRejecFac[trigSteps[i][1]][trigSteps[i][2]] << endl;
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
            offSetsEta[4] = 3; //EGA
          }else if(mode == 4){
            offSetsEta[1] = 0; //INT7
            offSetsEta[3] = 0; //EMC7
            offSetsEta[4] = 3; //EGA
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
                if(!replace(sysFileEtaDet, "Averaged", "AveragedSingle")){
                  cout << "WARNING: could not find detailed systematics file " << sysFileEtaDet << ", skipping... " << endl;
                  sysAvailSingleEta[i] = kFALSE;
                }else{
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
              for (Int_t f = 1; f < histoCorrectedYieldOmegaScaledMasked[i]->GetNbinsX()+1; f++ ){
                  histoCorrectedYieldOmegaScaledMasked[i]->SetBinContent(f,0.);
                  histoCorrectedYieldOmegaScaledMasked[i]->SetBinError(f,0.);
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

            if (graphsCorrectedYieldSysRemoved0Eta[i]){
                if (sysAvailSingleEta[i]){
                    nRelSysErrEtaSources                    = (Int_t)ptSysDetail[i][0].size()-1;
                    for (Int_t k = 0; k < nRelSysErrEtaSources; k++ ){
                        graphRelSysErrEtaSource[k][i]       = (TGraphAsymmErrors*) graphsCorrectedYieldSysRemoved0Eta[i]->Clone(Form("RelSysErrEtaSource%s_%s",((TString)ptSysDetail[i][0].at(k+1)).Data(), triggerName[i].Data()));
                        cout << Form("RelSysErrSource%s_%s",((TString)ptSysDetail[i][0].at(k+1)).Data(), triggerName[i].Data()) << endl;
                    }
                }
            }
            // shrink systematics graphs at the end & fill systematics
            for (Int_t j = 0; j< graphsCorrectedYieldSysRemoved0Eta[i]->GetN(); j++){
                if (sysAvailEta[i]){
                    Int_t counter2 = 0;
                    while(counter2 < 100 && TMath::Abs(graphsCorrectedYieldSysRemoved0Eta[i]->GetX()[j] - ptSysRelEta[i][counter2])> 0.001) counter2++;
                    if (counter2 < 100){
                        cout << ptSysRelEta[i][counter2]<< "\t found it" << endl;
                        Double_t yErrorSysLowDummy = TMath::Abs(yErrorSysLowRelEta[i][counter2]/100*graphsCorrectedYieldSysRemoved0Eta[i]->GetY()[j]);
                        Double_t yErrorSysHighDummy = yErrorSysHighRelEta[i][counter2]/100*graphsCorrectedYieldSysRemoved0Eta[i]->GetY()[j];
                        graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYlow(j,yErrorSysLowDummy);
                        graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYhigh(j,yErrorSysHighDummy);
                        if (sysAvailSingleEta[i]){
                            for (Int_t k = 0; k < nRelSysErrEtaSources; k++ ){
                                graphRelSysErrEtaSource[k][i]->SetPoint(j, graphsCorrectedYieldSysRemoved0Eta[i]->GetX()[j] ,((TString)ptSysDetail[i][counter2+1].at(k+1)).Atof());
                                graphRelSysErrEtaSource[k][i]->SetPointEYhigh(j,0);
                                graphRelSysErrEtaSource[k][i]->SetPointEYlow(j,0);
                            }
                        }
                    } else {
                        graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYlow(j,0);
                        graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYhigh(j,0);
                        if (sysAvailSingleEta[i]){
                            for (Int_t k = 0; k < nRelSysErrEtaSources; k++ ){
                                graphRelSysErrEtaSource[k][i]->SetPoint(j, graphsCorrectedYieldSysRemoved0Eta[i]->GetX()[j] ,0);
                                graphRelSysErrEtaSource[k][i]->SetPointEYlow(j,0);
                                graphRelSysErrEtaSource[k][i]->SetPointEYhigh(j,0);
                            }
                        }
                    }
                } else {
                    graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYlow(j,0);
                    graphsCorrectedYieldSysRemoved0Eta[i]->SetPointEYhigh(j,0);
                    averagedEta = kFALSE;
                    if (sysAvailSingleEta[i]){
                        for (Int_t k = 0; k < nRelSysErrEtaSources; k++ ){
                            graphRelSysErrEtaSource[k][i]->SetPoint(j, graphsCorrectedYieldSysRemoved0Eta[i]->GetX()[j] ,0);
                            graphRelSysErrEtaSource[k][i]->SetPointEYlow(j,0);
                            graphRelSysErrEtaSource[k][i]->SetPointEYhigh(j,0);
                        }
                    }
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
                    Int_t counter2 = 0;
                    while(counter2 < 100 && TMath::Abs(xValueFinalEta[nPointFinalEta] - ptSysRelEta[i][counter2])> 0.001) counter2++;
                    if (counter2 < 100){
                        cout << ptSysRelEta[i][counter2]<< "\t found it" << endl;
                        yErrorSysLowFinalEta[nPointFinalEta] = TMath::Abs(yErrorSysLowRelEta[i][counter2]/100*graphsCorrectedYieldShrunkEta[i]->GetY()[j]);
                        yErrorSysHighFinalEta[nPointFinalEta] = yErrorSysHighRelEta[i][counter2]/100*graphsCorrectedYieldShrunkEta[i]->GetY()[j];

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

            // Set correct trigger order for combination function
            Int_t nCorrOrder    = GetOrderedTrigger(triggerName[i]);
            if (nCorrOrder == -1){
                cout << "ERROR: trigger name not defined" << endl;
                return;
            }

            // fill inputs for combinations function in correct order
            if ( graphsCorrectedYieldShrunkEta[i]){
                histoStatEta[nCorrOrder]     = histoCorrectedYieldEtaScaledMasked[i];
                graphSystEta[nCorrOrder]     = graphsCorrectedYieldSysShrunkEta[i];
                offSetsEtaSys[nCorrOrder]    = histoStatEta[nCorrOrder]->GetXaxis()->FindBin(graphSystEta[nCorrOrder]->GetX()[0])-1;
                if (graphMassEtaData[i])
                    graphOrderedMassEtaData[nCorrOrder]      = graphMassEtaData[i];
                if (graphMassEtaMC[i])
                    graphOrderedMassEtaMC[nCorrOrder]        = graphMassEtaMC[i];
                if (graphWidthEtaData[i])
                    graphOrderedWidthEtaData[nCorrOrder]     = graphWidthEtaData[i];
                if (graphWidthEtaMC[i])
                    graphOrderedWidthEtaMC[nCorrOrder]       = graphWidthEtaMC[i];
                if (graphAcceptanceEta[i])
                    graphOrderedAcceptanceEta[nCorrOrder]    = graphAcceptanceEta[i];
                if (graphEfficiencyEta[i])
                    graphOrderedEfficiencyEta[nCorrOrder]    = graphEfficiencyEta[i];
                if (graphEffTimesAccEta[i])
                    graphOrderedEffTimesAccEta[nCorrOrder]   = graphEffTimesAccEta[i];
                if (sysAvailSingleEta[i]){
                    for (Int_t k = 0; k < nRelSysErrEtaSources; k++ ){
                        if (graphRelSysErrEtaSource[k][i])
                            graphOrderedRelSysErrEtaSource[k][nCorrOrder]    = graphRelSysErrEtaSource[k][i];
                    }
                }
            }
            if ((triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA")) && optionEnergy.CompareTo("8TeV")==0 && (mode == 4 || mode == 2))
              offSetsEtaSys[1]+=0; //INT7
              offSetsEtaSys[3]+=0; //EMC7
              offSetsEtaSys[4]+=3; //EGA
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
                TLatex *labelWeightsEta = new TLatex(0.7,0.20,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
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

            const char *SysErrDatnameMeanSingleErrCheck = Form("%s/SystematicErrorAveragedSingle%s_Eta_%s_Check.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
            fstream SysErrDatAverSingleCheck;
            SysErrDatAverSingleCheck.precision(4);
            cout << SysErrDatnameMeanSingleErrCheck << endl;
            if(sysAvailSingleEta[0]){
              SysErrDatAverSingleCheck.open(SysErrDatnameMeanSingleErrCheck, ios::out);
              SysErrDatAverSingleCheck << "pt \t Stat err \t sys err \t tot err " << endl;
              for (Int_t i = 0; i < graphRelErrorEtaTot->GetN(); i++){
                  if (graphRelErrorEtaStat->GetY()[i] > 0) SysErrDatAverSingleCheck << graphRelErrorEtaStat->GetX()[i] << "\t" << graphRelErrorEtaStat->GetY()[i] <<"\t" << graphRelErrorEtaSys->GetY()[i] <<  "\t" << graphRelErrorEtaTot->GetY()[i] << endl;
              }
              SysErrDatAverSingleCheck << endl;
              SysErrDatAverSingleCheck.close();
            }

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
                TLatex *labelRelErrEta       = new TLatex(0.15,0.85,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
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

            // Calculate relative sys error weighted
            if (sysAvailSingleEta[0]){
                for (Int_t k = 0; k< nRelSysErrEtaSources ; k++ ){
                    graphRelSysErrEtaSourceWeighted[k]      = CalculateWeightedQuantity(    graphOrderedRelSysErrEtaSource[k],
                                                                                            graphWeightsEta,
                                                                                            binningEta,  maxNAllowedEta,
                                                                                            MaxNumberOfFiles
                                                                                    );
                    if (!graphRelSysErrEtaSourceWeighted[k]){
                        cout << "Aborted in CalculateWeightedQuantity for " << endl;
                        return;
                    } else {
                        graphRelSysErrEtaSourceWeighted[k]->SetName(Form("RelSysErrEtaSourceWeighted%s", ((TString)ptSysDetail[0][0].at(k+1)).Data()));
                    }
                }
            }


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
            if (sysAvailSingleEta[0]){
                for (Int_t k = 0; k < nRelSysErrEtaSources; k++){
                    if (graphRelSysErrEtaSourceWeighted[k])
                        while (graphRelSysErrEtaSourceWeighted[k]->GetY()[0] == -10000 )   graphRelSysErrEtaSourceWeighted[k]->RemovePoint(0);
                    else
                        cout << "I don't have a weighted Eta rel sys err graph for error source: " << k << endl;
                }
            }


            //  **********************************************************************************************************************
            //  **************************************** Combine+write detailed Systematics ******************************************
            //  **********************************************************************************************************************

            const char *SysErrDatnameMeanSingleErr = Form("%s/SystematicErrorAveragedSingle%s_Eta_%s.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
            fstream SysErrDatAverSingle;
            SysErrDatAverSingle.precision(4);
            cout << SysErrDatnameMeanSingleErr << endl;
            if(sysAvailSingleEta[0] && graphRelSysErrEtaSourceWeighted[0]){
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

                for(Int_t i=0; i<graphRelSysErrEtaSourceWeighted[0]->GetN(); i++){
                    SysErrDatAverSingle << graphRelSysErrEtaSourceWeighted[0]->GetX()[i] << "\t";
                    Int_t nColumns = (Int_t)ptSysDetail[0][0].size();
                    for(Int_t iErr=0; iErr<nColumns-1; iErr++)
                        SysErrDatAverSingle << graphRelSysErrEtaSourceWeighted[iErr]->GetY()[i] << "\t";
                    SysErrDatAverSingle << endl;

                }
            }
            SysErrDatAverSingle.close();

            // ***************************************************************************************************
            // ********************* Plot all mean erros separately after smoothing ******************************
            // ***************************************************************************************************
            if(sysAvailSingleEta[0] && graphRelSysErrEtaSourceWeighted[0]){
                TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);// gives the page size
                DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);

                    // create dummy histo
                    TH2D *histo2DNewSysErrMean ;
                    histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 100,0.,maxPtGlobalEta,1000.,-0.5,50.);
                    SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                            1,0.9, 510, 510);
                    histo2DNewSysErrMean->Draw();

                    // Give legend position for plotting
                    Double_t minXLegend     = 0.12;
                    Double_t maxYLegend     = 0.95;
                    Double_t widthLegend    = 0.25;
                    if (nRelSysErrEtaSources> 7)
                        widthLegend         = 0.5;
                    Double_t heightLegend   = 1.05* 0.035 * (nRelSysErrEtaSources+3);
                    if (nRelSysErrEtaSources> 7)
                        heightLegend        = 1.05* 0.035 * (nRelSysErrEtaSources/2+1);

                    // create legend
                    TLegend* legendMeanNew = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
                    legendMeanNew->SetMargin(0.1);
                    if (nRelSysErrEtaSources> 7) legendMeanNew->SetNColumns(2);

                    for(Int_t i = 0;i< nRelSysErrEtaSources-1 ;i++){
                        DrawGammaSetMarkerTGraphAsym(graphRelSysErrEtaSourceWeighted[i], GetMarkerStyleSystematics(  ptSysDetail[0][0].at(i+1), mode), 1.,
                                                    GetColorSystematics( ptSysDetail[0][0].at(i+1), mode),GetColorSystematics( ptSysDetail[0][0].at(i+1), mode));
                        graphRelSysErrEtaSourceWeighted[i]->Draw("pX0,csame");
                        legendMeanNew->AddEntry(graphRelSysErrEtaSourceWeighted[i],GetSystematicsName(ptSysDetail[0][0].at(i+1)),"p");
                    }

                    DrawGammaSetMarkerTGraphAsym(graphRelSysErrEtaSourceWeighted[nRelSysErrEtaSources-1], 20, 1.,kBlack,kBlack);
                    graphRelSysErrEtaSourceWeighted[nRelSysErrEtaSources-1]->Draw("p,csame");
                    legendMeanNew->AddEntry(graphRelSysErrEtaSourceWeighted[nRelSysErrEtaSources-1],"quad. sum.","p");
                    legendMeanNew->Draw();

                    // labeling
                    TLatex *labelEnergySysDetailed = new TLatex(0.7, 0.93,collisionSystem.Data());
                    labelEnergySysDetailed->SetTextAlign(31);
                    SetStyleTLatex( labelEnergySysDetailed, 0.85*textSizeSpectra,4);
                    labelEnergySysDetailed->Draw();

                    TLatex *labelEtaSysDetailed     = new TLatex(0.7, 0.93-0.99*textSizeSpectra*0.85,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
                    labelEtaSysDetailed->SetTextAlign(31);
                    SetStyleTLatex( labelEtaSysDetailed, 0.85*textSizeSpectra,4);
                    labelEtaSysDetailed->Draw();

                    TLatex *labelDetProcSysDetailed = new TLatex(0.7, 0.93-2*0.99*textSizeSpectra*0.85,detectionProcess.Data());
                    labelDetProcSysDetailed->SetTextAlign(31);
                    SetStyleTLatex( labelDetProcSysDetailed, 0.85*textSizeSpectra,4);
                    labelDetProcSysDetailed->Draw();

                canvasNewSysErrMean->Update();
                canvasNewSysErrMean->SaveAs(Form("%s/Eta_SysErrorsSeparatedSourcesReweighted_%s.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
            }

            // delete detailed sys array
            for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++){
                for(Int_t iB=0; iB<50; iB++) ptSysDetail[iR][iB].clear();
            }

        // if averaging wasn't enabled pick values according to predefined ranges ("cherry picking points")
        } else {
            graphCorrectedYieldWeightedAverageEtaStat       = new TGraphAsymmErrors(nPointFinalEta, xValueFinalEta, yValueFinalEta,
                                                                                        xErrorLowFinalEta, xErrorHighFinalEta,yErrorLowFinalEta, yErrorHighFinalEta);
            graphCorrectedYieldWeightedAverageEtaSys        = new TGraphAsymmErrors(nPointFinalEta, xValueFinalEta, yValueFinalEta,
                                                                                    xErrorLowFinalEta, xErrorHighFinalEta,yErrorSysLowFinalEta, yErrorSysHighFinalEta);
            graphMassEtaDataWeighted                        = graphMassEtaData[0];
            graphMassEtaMCWeighted                          = graphMassEtaMC[0];
            graphWidthEtaDataWeighted                       = graphWidthEtaData[0];
            graphWidthEtaMCWeighted                         = graphWidthEtaMC[0];

            graphAcceptanceEtaWeighted                      = graphAcceptanceEta[0];
            graphEfficiencyEtaWeighted                      = graphEfficiencyEta[0];
            graphEffTimesAccEtaWeighted                     = graphEffTimesAccEta[0];

            TGraphAsymmErrors* graphRelErrorEtaStat       = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaStat, "relativeStatErrorEta");
            while (graphRelErrorEtaStat->GetY()[0] < 0 ) graphRelErrorEtaStat->RemovePoint(0);

            TGraphAsymmErrors* graphRelErrorEtaSys        = CalculateRelErrUpAsymmGraph( graphCorrectedYieldWeightedAverageEtaSys, "relativeSysErrorEta");
            while (graphRelErrorEtaSys->GetY()[0] < 0 ) graphRelErrorEtaSys->RemovePoint(0);

            const char *SysErrDatnameMeanSingleErrCheck = Form("%s/SystematicErrorAveragedSingle%s_Eta_%s_Check.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
            fstream SysErrDatAverSingleCheck;
            SysErrDatAverSingleCheck.precision(4);
            cout << SysErrDatnameMeanSingleErrCheck << endl;

            SysErrDatAverSingleCheck.open(SysErrDatnameMeanSingleErrCheck, ios::out);
            SysErrDatAverSingleCheck << "pt \t Stat err \t sys err \t tot err " << endl;
            for (Int_t i = 0; i < graphRelErrorEtaStat->GetN(); i++){
                if (graphRelErrorEtaStat->GetY()[i] > 0) SysErrDatAverSingleCheck << graphRelErrorEtaStat->GetX()[i] << "\t" << graphRelErrorEtaStat->GetY()[i] <<"\t" << graphRelErrorEtaSys->GetY()[i] << endl;
            }
            SysErrDatAverSingleCheck << endl;
            SysErrDatAverSingleCheck.close();

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
             (optionEnergy.CompareTo("8TeV")==0) ||
             (optionEnergy.CompareTo("pPb_5.023TeV")==0)
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
        canvasMass->SaveAs(Form("%s/Eta_%s_Mass1_Reduced_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

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
            canvasMass->SaveAs(Form("%s/Eta_%s_Mass2_Reduced_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
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
            canvasMass->SaveAs(Form("%s/Eta_%s_Mass_Weighted_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

            TGraphAsymmErrors* graphMassDifferenceEtaDatavsMC       = CalculateAsymGraphDifferenceToGraph(graphMassEtaDataWeighted,graphMassEtaMCWeighted);
            TGraphAsymmErrors* graphMassRelDifferenceEtaDatavsMC    = CalculateAsymGraphRatioToGraph(graphMassDifferenceEtaDatavsMC,graphMassEtaMCWeighted);
            graphMassRelDifferenceEtaDatavsMC = ScaleGraph(graphMassRelDifferenceEtaDatavsMC, 100.);

            canvasMass->cd();
            TH2F * histo2DRelMassDiffEta       = new TH2F("histo2DRelMassDiffEta","histo2DRelMassDiffEta",1000,0., maxPtGlobalEta,10000,-5, 5);
            SetStyleHistoTH2ForGraphs(histo2DRelMassDiffEta, "#it{p}_{T} (GeV/#it{c})","#it{M}_{#eta, data}-#it{M}_{#eta, MC}/ #it{M}_{#eta, MC} (%)",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.4);
            histo2DRelMassDiffEta->DrawCopy();

            if (graphMassRelDifferenceEtaDatavsMC){
              TF1 *fitEtaRelMassDiff = new TF1("fitEtaRelMassDiff","[0]",0.,maxPtGlobalEta);
              fitEtaRelMassDiff->SetParameter(0,0.);
              graphMassRelDifferenceEtaDatavsMC->Fit(fitEtaRelMassDiff,"QNRMEX0+","",0.,maxPtGlobalEta);
              fitEtaRelMassDiff->SetLineColor(kRed);
              fitEtaRelMassDiff->SetLineWidth(0.5);
              fitEtaRelMassDiff->Draw("same");

              fileFitsOutput << "average rel mass diff: " << fitEtaRelMassDiff->GetParameter(0) << "+-"<< fitEtaRelMassDiff->GetParError(0) << endl;

              TLegend* legendEtaRelMassDiff = GetAndSetLegend2(0.15, 0.12, 0.6, 0.12+(0.035*1), 0.035, 2, "", 42, 0.15);
              legendEtaRelMassDiff->AddEntry(fitEtaRelMassDiff,"fit with constant");
              legendEtaRelMassDiff->AddEntry((TObject*)0,Form("%0.4f #pm %0.4f", fitEtaRelMassDiff->GetParameter(0), fitEtaRelMassDiff->GetParError(0)),"");
              legendEtaRelMassDiff->Draw();

              DrawGammaSetMarkerTGraphAsym(graphMassRelDifferenceEtaDatavsMC, 20, 1, kBlack, kBlack);
              graphMassRelDifferenceEtaDatavsMC->Draw("p,e1,same");
            }
            DrawGammaLines(0., maxPtGlobalEta , 0., 0., 1, kGray+2, 7);
//             legendMassEtaWeighted->Draw();
            labelEnergyMass->Draw();
            labelEtaMass->Draw();
            labelDetProcMass->Draw();

            canvasMass->Update();
            canvasMass->SaveAs(Form("%s/Eta_%s_RelMassDiff_Weighted.%s",outputDir.Data(),isMC.Data(),suffix.Data()));

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

            TLatex *labelEtaEffiWOTrigg = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
            SetStyleTLatex( labelEtaEffiWOTrigg, 0.85*textSizeSpectra,4);
            labelEtaEffiWOTrigg->Draw();

            TLatex *labelDetProcEffiWOTrigg = new TLatex(0.62, 0.15,detectionProcess.Data());
            SetStyleTLatex( labelDetProcEffiWOTrigg, 0.85*textSizeSpectra,4);
            labelDetProcEffiWOTrigg->Draw();

            canvasEffi->Update();
            canvasEffi->SaveAs(Form("%s/Eta_EfficiencyW0TriggEff_Weighted_%i.%s",outputDir.Data(),mode,suffix.Data()));
        }
        delete histo2DEffiEta;

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

            TLatex *labelEtaAcc = new TLatex(0.62, 0.15+0.99*textSizeSpectra*0.85,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
            SetStyleTLatex( labelEtaAcc, 0.85*textSizeSpectra,4);
            labelEtaAcc->Draw();

            TLatex *labelDetProcAcc = new TLatex(0.62, 0.15,detectionProcess.Data());
            SetStyleTLatex( labelDetProcAcc, 0.85*textSizeSpectra,4);
            labelDetProcAcc->Draw();

            canvasAcc->Update();
            canvasAcc->SaveAs(Form("%s/Eta_Acceptance_weighted_%i.%s",outputDir.Data(),mode,suffix.Data()));
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
             (optionEnergy.CompareTo("8TeV")==0) ||
             (optionEnergy.CompareTo("pPb_5.023TeV")==0)
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
        canvasWidth->SaveAs(Form("%s/Eta_%s_Width1_Reduced_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

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
            canvasWidth->SaveAs(Form("%s/Eta_%s_Width2_Reduced_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
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
        canvasWidth->SaveAs(Form("%s/Eta_%s_Width_Weighted_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));
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
                minCorrYieldEta     = 2e-11;
                maxCorrYieldEta     = 5e-2;
            }else if(mode == 4){
                minCorrYieldEta     = 3e-11;
                maxCorrYieldEta     = 7e-3;
            }
        }
        if(optionEnergy.CompareTo("pPb_5.023TeV")==0){
            if(mode == 2){
                minCorrYieldEta     = 1e-8;
                maxCorrYieldEta     = 1e-1;
            }else if(mode == 4){
                minCorrYieldEta     = 3e-9;
                maxCorrYieldEta     = 7e-2;
            }else if(mode == 3){
                minCorrYieldEta     = 7e-6;
                maxCorrYieldEta     = 0.3;
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
        DrawGammaSetMarkerTF1( fitInvYieldEta, 7, 2, kGray+2);
        fitInvYieldEta->Draw("same");

        labelEnergyUnscaled->Draw();
        labelEtaUnscaled->Draw();
        labelDetProcUnscaled->Draw();

        canvasCorrScaled->Update();
        canvasCorrScaled->SaveAs(Form("%s/Eta_%s_CorrectedYieldScaledTrigg_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

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
        canvasCorrScaled->SaveAs(Form("%s/Eta_%s_CorrectedYieldFinal_%i.%s",outputDir.Data(),isMC.Data(),mode,suffix.Data()));

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

        TLatex *labelEtaRatio = new TLatex(0.6, 0.93-textSizeSpectra*0.85*1.04, "#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
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
        //************************************ Eta to Omega ratio for different triggers **********************************
        //***************************************************************************************************************
        if (doEtaToOmega){
            cout << "**************************************************************************************************" << endl;
            cout << "************************ Combining different triggers for eta to Omega *****************************" << endl;
            cout << "**************************************************************************************************" << endl;

            // prepare systematics arrays
            Double_t xValueFinalEtaToOmega                    [100];
            Double_t xErrorLowFinalEtaToOmega                 [100];
            Double_t xErrorHighFinalEtaToOmega                [100];
            Double_t yValueFinalEtaToOmega                    [100];
            Double_t yErrorLowFinalEtaToOmega                 [100];
            Double_t yErrorHighFinalEtaToOmega                [100];
            Int_t nPointFinalEtaToOmega                                         = 0;

            Double_t yErrorSysLowFinalEtaToOmega              [100];
            Double_t yErrorSysHighFinalEtaToOmega             [100];

            Double_t ptSysRelEtaToOmega                       [MaxNumberOfFiles][100];
            Int_t ptSysRelEtaToOmegaNBins                     = 0;
            Double_t yErrorSysLowRelEtaToOmega                [MaxNumberOfFiles][100];
            Double_t yErrorSysHighRelEtaToOmega               [MaxNumberOfFiles][100];
            Bool_t   sysAvailEtaToOmega                       [MaxNumberOfFiles];

            Int_t numberBinsSysAvailSingleEtaToOmega          [MaxNumberOfFiles];

            TGraphAsymmErrors* graphEtaToOmegaBinShiftWeighted = NULL;
            TGraphAsymmErrors* graphEtaToOmegaBinShift        [MaxNumberOfFiles];


            // create graphs for shrunk individual triggers
            TGraphAsymmErrors* graphsEtaToOmegaShrunk         [MaxNumberOfFiles];
            TGraphAsymmErrors* graphsEtaToOmegaSysShrunk      [MaxNumberOfFiles];
            TH1D* histoEtaToOmegaMasked                       [MaxNumberOfFiles];

            // create variables for combination functions
            TH1D*               histoStatEtaToOmega    [12];
            TGraphAsymmErrors*  graphSystEtaToOmega    [12];
            TH1D*               histoRelStatEtaToOmega [12];
            TGraphAsymmErrors*  graphRelSystEtaToOmega [12];
            Int_t offSetsEtaToOmega[12]           = { 0, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 0 };
            Int_t offSetsEtaToOmegaSys[12]        = { 0, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 0 };

            if(optionEnergy.CompareTo("8TeV")==0){
              if(mode == 2){
                offSetsEtaToOmega[1] = 0; //INT7
                offSetsEtaToOmega[3] = 0; //EMC7
                offSetsEtaToOmega[4] = 3; //EGA
              }else if(mode == 4){
                offSetsEtaToOmega[1] = 0; //INT7
                offSetsEtaToOmega[3] = 0; //EMC7
                offSetsEtaToOmega[4] = 3; //EGA
              }
            }


            Bool_t hasSysEtaToOmega            = kFALSE;
            for (Int_t j = 0; j<12; j++){
                histoStatEtaToOmega[j]        = NULL;
                graphSystEtaToOmega[j]        = NULL;
                histoRelStatEtaToOmega[j]     = NULL;
                graphRelSystEtaToOmega[j]     = NULL;
                graphEtaToOmegaBinShift[j]    = NULL;
            }

            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                // read systematics file
                cout << triggerName[i].Data() << endl;
                if (sysFileEtaToOmega[i].CompareTo("bla") != 0){
                    sysAvailEtaToOmega[i]              = kTRUE;
                    ifstream  fileSysErrEtaToOmega;
                    fileSysErrEtaToOmega.open(sysFileEtaToOmega[i].Data(),ios_base::in);
                    cout << sysFileEtaToOmega[i].Data() << endl;
                    Int_t counter = 0;
                    while(!fileSysErrEtaToOmega.eof() && counter < 100){
                        Double_t garbage = 0;
                        fileSysErrEtaToOmega >>ptSysRelEtaToOmega[i][counter] >> yErrorSysLowRelEtaToOmega[i][counter] >> yErrorSysHighRelEtaToOmega[i][counter]>>    garbage >> garbage;
                        cout << counter << "\t"<< ptSysRelEtaToOmega[i][counter]<< "\t"  << yErrorSysLowRelEtaToOmega[i][counter] << "\t"  <<yErrorSysHighRelEtaToOmega[i][counter] << "\t"  << endl;;
                        counter++;
                    }
                    ptSysRelEtaToOmegaNBins = counter;
                    fileSysErrEtaToOmega.close();
                    hasSysEtaToOmega              = kTRUE;
                 // read in detailed systematics
                    string sysFileEtaToOmegaDet = sysFileEtaToOmega[i].Data();
                    if(!replace(sysFileEtaToOmegaDet, "Averaged", "AveragedSingle")){
                      cout << "WARNING: could not find detailed systematics file " << sysFileEtaToOmegaDet << ", skipping... " << endl;
                      sysAvailSingleEtaToOmega[i] = kFALSE;
                    }else{
                      ifstream fileSysErrDetailedEtaToOmega;
                      fileSysErrDetailedEtaToOmega.open(sysFileEtaToOmegaDet,ios_base::in);
                      if(fileSysErrDetailedEtaToOmega.is_open())
                          sysAvailSingleEtaToOmega[i] = kTRUE;
                      else{
                          sysAvailSingleEtaToOmega[i] = kFALSE;
                          cout << "couldn't find single errors for eta/Omega, jumping" << endl;
                      }

                      if (sysAvailSingleEtaToOmega[i]){
                          cout << sysFileEtaToOmegaDet << endl;
                          counter = 0;
                          string line;
                          Int_t counterColumn = 0;
                          while (getline(fileSysErrDetailedEtaToOmega, line) && counter < 100) {
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
                          numberBinsSysAvailSingleEtaToOmega[i] = counter;
                          fileSysErrDetailedEtaToOmega.close();
                      }
                    }
                } else {
                    sysAvailEtaToOmega[i]         = kFALSE;
                    sysAvailSingleEtaToOmega[i]   = kTRUE;
                }

                // fill graphs to be shrunk later
                cout << "step 1" << endl;
                graphsEtaToOmegaShrunk[i]         = new TGraphAsymmErrors(histoEtaToOmega[i]);
                graphsEtaToOmegaRemoved0[i]       = new TGraphAsymmErrors(histoEtaToOmega[i]);
                graphsEtaToOmegaSysShrunk[i]      = new TGraphAsymmErrors(histoEtaToOmega[i]);
                graphsEtaToOmegaSysRemoved0[i]    = new TGraphAsymmErrors(histoEtaToOmega[i]);
                histoEtaToOmegaMasked[i]          = (TH1D*)histoEtaToOmega[i]->Clone(Form("EtaToOmega%s_Masked_%s",addNameBinshift.Data(), triggerName[i].Data()));

                if(optionEnergy.CompareTo("8TeV")==0 && mode==4 && triggerName[i].Contains("EGA")) ptFromSpecEta[i][0] = 16.0;
                cout << ptFromSpecEta[i][0] << endl;
                // remove 0 bins at beginning according to ptFromSpecEta[i][0]
                Int_t binsToMask = 1;
                while (histoEtaToOmegaMasked[i]->GetBinCenter(binsToMask) < ptFromSpecEta[i][0] ){
                    histoEtaToOmegaMasked[i]->SetBinContent(binsToMask,0.);
                    histoEtaToOmegaMasked[i]->SetBinError(binsToMask,0.);
                    binsToMask++;
                }
                // check if trigger needs to be masked completely
                if ( maskedFullyEta[i] || maskedFullyOmega[i] ){
                    graphsEtaToOmegaShrunk[i]         = NULL;
                    graphsEtaToOmegaRemoved0[i]       = NULL;
                    graphsEtaToOmegaSysShrunk[i]      = NULL;
                    graphsEtaToOmegaSysRemoved0[i]    = NULL;
                    for (Int_t f = 1; f < histoEtaToOmegaMasked[i]->GetNbinsX()+1; f++ ){
                        histoEtaToOmegaMasked[i]->SetBinContent(f,0.);
                        histoEtaToOmegaMasked[i]->SetBinError(f,0.);
                    }
                    nrOfTrigToBeCombEtaToOmegaRed--;
                    cout << "trigger " << triggerName[i] << " was masked for Eta/Omega" << endl;
                    continue;
                // remove bins at the end according to set range
                } else if (ptFromSpecEta[i][1] > -1){
                    for (Int_t f = histoEtaToOmegaMasked[i]->GetXaxis()->FindBin(ptFromSpecEta[i][1]); f < histoEtaToOmegaMasked[i]->GetNbinsX()+1; f++ ){
                        histoEtaToOmegaMasked[i]->SetBinContent(f,0.);
                        histoEtaToOmegaMasked[i]->SetBinError(f,0.);
                    }
                }

                // remove 0 bins at beginning
//                 graphsEtaToOmegaShrunk[i]->Print();
                cout << "step 2" << endl;
                while (graphsEtaToOmegaShrunk[i]->GetY()[0] == 0) graphsEtaToOmegaShrunk[i]->RemovePoint(0);
//                 graphsEtaToOmegaShrunk[i]->Print();
                while (graphsEtaToOmegaRemoved0[i]->GetY()[0] == 0) graphsEtaToOmegaRemoved0[i]->RemovePoint(0);
//                 graphsEtaToOmegaRemoved0[i]->Print();
                cout << "sys shrunk" << endl;
                while (graphsEtaToOmegaSysShrunk[i]->GetY()[0] == 0) graphsEtaToOmegaSysShrunk[i]->RemovePoint(0);
//                 graphsEtaToOmegaSysShrunk[i]->Print();
                cout << "sys shrunk 2" << endl;
                while (graphsEtaToOmegaSysRemoved0[i]->GetY()[0] == 0) graphsEtaToOmegaSysRemoved0[i]->RemovePoint(0);
//                 graphsEtaToOmegaSysRemoved0[i]->Print();

                if (graphsEtaToOmegaSysRemoved0[i]){
                    if (sysAvailSingleEtaToOmega[i]){
                        nRelSysErrEtaToOmegaSources                    = (Int_t)ptSysDetail[i][0].size()-1;
                        for (Int_t k = 0; k < nRelSysErrEtaToOmegaSources; k++ ){
                            graphRelSysErrEtaToOmegaSource[k][i]       = (TGraphAsymmErrors*) graphsEtaToOmegaSysRemoved0[i]->Clone(Form("RelSysErrEtaToOmegaSource%s_%s",((TString)ptSysDetail[i][0].at(k+1)).Data(), triggerName[i].Data()));
                            cout << Form("RelSysErrEtaToOmegaSource%s_%s",((TString)ptSysDetail[i][0].at(k+1)).Data(), triggerName[i].Data()) << endl;
                        }
                    }
                }

                // put systematics on graph
                for (Int_t j = 0; j< graphsEtaToOmegaSysRemoved0[i]->GetN(); j++){
                    if (sysAvailEtaToOmega[i]){
                        Int_t counter2 = 0;
                        while(counter2 < 100 && TMath::Abs(graphsEtaToOmegaSysRemoved0[i]->GetX()[j] - ptSysRelEtaToOmega[i][counter2])> 0.001) counter2++;
                        cout << "counter: " << counter2 << ", NBins read: " << ptSysRelEtaToOmegaNBins << endl;
                        if (counter2 < 100 && counter2 < ptSysRelEtaToOmegaNBins){
                            cout << "counter: " << counter2 << ", " << ptSysRelEtaToOmega[i][counter2]<< "\t found it" << endl;
                            Double_t yErrorSysLowDummy  = TMath::Abs(yErrorSysLowRelEtaToOmega[i][counter2]/100*graphsEtaToOmegaSysRemoved0[i]->GetY()[j]);
                            Double_t yErrorSysHighDummy = yErrorSysHighRelEtaToOmega[i][counter2]/100*graphsEtaToOmegaSysRemoved0[i]->GetY()[j];
                            graphsEtaToOmegaSysRemoved0[i]->SetPointEYlow(j,yErrorSysLowDummy);
                            graphsEtaToOmegaSysRemoved0[i]->SetPointEYhigh(j,yErrorSysHighDummy);
                            if (sysAvailSingleEtaToOmega[i]){
                                for (Int_t k = 0; k < nRelSysErrEtaToOmegaSources; k++ ){
                                    graphRelSysErrEtaToOmegaSource[k][i]->SetPoint(j, graphsEtaToOmegaSysRemoved0[i]->GetX()[j] ,((TString)ptSysDetail[i][counter2+1].at(k+1)).Atof());
                                    graphRelSysErrEtaToOmegaSource[k][i]->SetPointEYhigh(j,0);
                                    graphRelSysErrEtaToOmegaSource[k][i]->SetPointEYlow(j,0);
                                }
                            }
                        } else {
                            graphsEtaToOmegaSysRemoved0[i]->SetPointEYlow(j,0);
                            graphsEtaToOmegaSysRemoved0[i]->SetPointEYhigh(j,0);
                            if (sysAvailSingleEtaToOmega[i]){
                                for (Int_t k = 0; k < nRelSysErrEtaToOmegaSources; k++ ){
                                    graphRelSysErrEtaToOmegaSource[k][i]->SetPoint(j, graphsEtaToOmegaSysRemoved0[i]->GetX()[j] ,0);
                                    graphRelSysErrEtaToOmegaSource[k][i]->SetPointEYhigh(j,0);
                                    graphRelSysErrEtaToOmegaSource[k][i]->SetPointEYlow(j,0);
                                }
                            }
                        }
                    } else {
                        graphsEtaToOmegaSysRemoved0[i]->SetPointEYlow(j,0);
                        graphsEtaToOmegaSysRemoved0[i]->SetPointEYhigh(j,0);
                        averagedEta = kFALSE;
                        if (sysAvailSingleEtaToOmega[i]){
                            for (Int_t k = 0; k < nRelSysErrEtaToOmegaSources; k++ ){
                                graphRelSysErrEtaToOmegaSource[k][i]->SetPoint(j, graphsEtaToOmegaSysRemoved0[i]->GetX()[j] ,0);
                                graphRelSysErrEtaToOmegaSource[k][i]->SetPointEYhigh(j,0);
                                graphRelSysErrEtaToOmegaSource[k][i]->SetPointEYlow(j,0);
                            }
                        }
                    }
                }

                if(optionEnergy.CompareTo("8TeV")==0 && mode == 4) ptFromSpecOmega[i][0]-= 0.5;
                // remove unused bins at beginning
                cout << "step 3" << endl;
                while (graphsEtaToOmegaShrunk[i]->GetX()[0] < ptFromSpecEta[i][0] || graphsEtaToOmegaShrunk[i]->GetX()[0] < ptFromSpecOmega[i][0])
                    graphsEtaToOmegaShrunk[i]->RemovePoint(0);
                while (graphsEtaToOmegaSysShrunk[i]->GetX()[0] < ptFromSpecEta[i][0] || graphsEtaToOmegaSysShrunk[i]->GetX()[0] < ptFromSpecOmega[i][0])
                    graphsEtaToOmegaSysShrunk[i]->RemovePoint(0);
                // remove unused bins at end
                while ( graphsEtaToOmegaShrunk[i]->GetX()[graphsEtaToOmegaShrunk[i]->GetN()-1] > ptFromSpecEta[i][1] ||
                        graphsEtaToOmegaShrunk[i]->GetX()[graphsEtaToOmegaShrunk[i]->GetN()-1] > ptFromSpecOmega[i][1] )
                    graphsEtaToOmegaShrunk[i]->RemovePoint(graphsEtaToOmegaShrunk[i]->GetN()-1);
                graphsEtaToOmegaShrunk[i]->Print();
                while ( graphsEtaToOmegaSysShrunk[i]->GetX()[graphsEtaToOmegaSysShrunk[i]->GetN()-1] > ptFromSpecEta[i][1] ||
                        graphsEtaToOmegaSysShrunk[i]->GetX()[graphsEtaToOmegaSysShrunk[i]->GetN()-1] > ptFromSpecOmega[i][1] )
                    graphsEtaToOmegaSysShrunk[i]->RemovePoint(graphsEtaToOmegaSysShrunk[i]->GetN()-1);

                // put systematics on graph
                for (Int_t j = 0; j< graphsEtaToOmegaShrunk[i]->GetN(); j++){
                    xValueFinalEtaToOmega[nPointFinalEtaToOmega]        = graphsEtaToOmegaShrunk[i]->GetX()[j];
                    xErrorHighFinalEtaToOmega[nPointFinalEtaToOmega]    = graphsEtaToOmegaShrunk[i]->GetEXhigh()[j];
                    xErrorLowFinalEtaToOmega[nPointFinalEtaToOmega]     = graphsEtaToOmegaShrunk[i]->GetEXlow()[j];
                    yValueFinalEtaToOmega[nPointFinalEtaToOmega]        = graphsEtaToOmegaShrunk[i]->GetY()[j];
                    yErrorHighFinalEtaToOmega[nPointFinalEtaToOmega]    = graphsEtaToOmegaShrunk[i]->GetEYhigh()[j];
                    yErrorLowFinalEtaToOmega[nPointFinalEtaToOmega]     = graphsEtaToOmegaShrunk[i]->GetEYlow()[j];
                    if (sysAvailEtaToOmega[i]){
                        Int_t counter = 0;
                        while(counter < 100 && TMath::Abs(xValueFinalEtaToOmega[nPointFinalEtaToOmega] - ptSysRelEtaToOmega[i][counter])> 0.001) counter++;
                        if (counter < 100){
                            cout << ptSysRelEtaToOmega[i][counter]<< "\t found it" << endl;
                            yErrorSysLowFinalEtaToOmega[nPointFinalEtaToOmega] = TMath::Abs(yErrorSysLowRelEtaToOmega[i][counter]/100*graphsEtaToOmegaShrunk[i]->GetY()[j]);
                            yErrorSysHighFinalEtaToOmega[nPointFinalEtaToOmega] = yErrorSysHighRelEtaToOmega[i][counter]/100*graphsEtaToOmegaShrunk[i]->GetY()[j];

                        } else {
                            yErrorSysLowFinalEtaToOmega[nPointFinalEtaToOmega]  = 0;
                            yErrorSysHighFinalEtaToOmega[nPointFinalEtaToOmega] = 0;

                        }
                    } else {
                        yErrorSysLowFinalEtaToOmega[nPointFinalEtaToOmega]  = 0;
                        yErrorSysHighFinalEtaToOmega[nPointFinalEtaToOmega] = 0;
                    }
                    graphsEtaToOmegaSysShrunk[i]->SetPointEYlow(j,yErrorSysLowFinalEtaToOmega[nPointFinalEtaToOmega]);
                    graphsEtaToOmegaSysShrunk[i]->SetPointEYhigh(j,yErrorSysHighFinalEtaToOmega[nPointFinalEtaToOmega]);
                    nPointFinalEtaToOmega++;
                }

                // Set correct trigger order for combination function
                Int_t nCorrOrder    = GetOrderedTrigger(triggerName[i]);
                if (nCorrOrder == -1){
                    cout << "ERROR: trigger name not defined" << endl;
                    return;
                }

                // fill inputs for combinations function in correct order
                if ( graphsEtaToOmegaShrunk[i]){
                    histoStatEtaToOmega[nCorrOrder]     = histoEtaToOmegaMasked[i];
                    graphSystEtaToOmega[nCorrOrder]     = graphsEtaToOmegaSysShrunk[i];
                    offSetsEtaToOmegaSys[nCorrOrder]    = histoStatEtaToOmega[nCorrOrder]->GetXaxis()->FindBin(graphSystEtaToOmega[nCorrOrder]->GetX()[0])-1;
                    if (histoCorrectedYieldEtaBinShift[i])
                        graphEtaToOmegaBinShift[nCorrOrder]   = new TGraphAsymmErrors(histoCorrectedYieldEtaBinShift[i]);
                    if (sysAvailSingleEtaToOmega[i]){
                        for (Int_t k = 0; k < nRelSysErrEtaToOmegaSources; k++ ){
                            if (graphRelSysErrEtaToOmegaSource[k][i])
                                graphOrderedRelSysErrEtaToOmegaSource[k][nCorrOrder]    = graphRelSysErrEtaToOmegaSource[k][i];
                        }
                    }
                }
                if ((triggerName[i].Contains("EG2") || triggerName[i].Contains("EGA")) && optionEnergy.CompareTo("8TeV")==0 && (mode == 4 || mode == 2))
                  offSetsEtaToOmegaSys[4]+=3; //EGA
            }

            TString nameWeightsLogFileEtaToOmega                  = Form("%s/weightsEtaToOmega_%s.dat",outputDir.Data(),isMC.Data());
            TGraphAsymmErrors* graphEtaToOmegaWeightedAverageTot  = NULL;
            // Calculate averaged eta/Omega graphs according to statistical and systematic errors taking correctly into account the cross correlations
            if (averagedEta){
                if(optionEnergy.CompareTo("8TeV")==0 && mode==4){
                  maxNAllowedEta -= 3;
                  maxPtGlobalEta = 20;
                }
                //if(optionEnergy.CompareTo("8TeV")==0 && mode==2) maxNAllowedEta -= 2;
                // calculate averaged eta/Omega graphs
                graphEtaToOmegaWeightedAverageTot         = CombinePtPointsSpectraTriggerCorrMat( histoStatEtaToOmega, graphSystEtaToOmega,
                                                                                                binningEta,  maxNAllowedEta,
                                                                                                offSetsEtaToOmega ,offSetsEtaToOmegaSys,
                                                                                                graphEtaToOmegaWeightedAverageStat, graphEtaToOmegaWeightedAverageSys,
                                                                                                nameWeightsLogFileEtaToOmega.Data(),
                                                                                                mode, optionEnergy, "EtaToOmega", kFALSE,
                                                                                                fileInputCorrFactors
                                                                                              );
                // Reading weights from output file for plotting
                ifstream fileWeightsEtaToOmega;
                fileWeightsEtaToOmega.open(nameWeightsLogFileEtaToOmega,ios_base::in);
                cout << "reading" << nameWeightsLogFileEtaToOmega << endl;
                Double_t xValuesReadEtaToOmega[100];
                Double_t weightsReadEtaToOmega[12][100];
                Int_t availableMeasEtaToOmega[12]     = { -1, -1, -1, -1, -1, -1,
                                                        -1, -1, -1, -1, -1, -1};
                Int_t nMeasSetEtaToOmega              = numberOfTrigg;
                Int_t nPtBinsReadEtaToOmega           = 0;

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

                while(!fileWeightsEtaToOmega.eof() && nPtBinsReadEtaToOmega < 100){
                    TString garbage = "";
                    if (nPtBinsReadEtaToOmega == 0){
                        fileWeightsEtaToOmega >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
                        for (Int_t i = 0; i < nMeasSetEtaToOmega; i++){
                            fileWeightsEtaToOmega >> availableMeasEtaToOmega[i] ;
                        }
                        cout << "read following measurements: ";
                        for (Int_t i = 0; i < nMeasSetEtaToOmega; i++){
                            cout << availableMeasEtaToOmega[i] << "\t" ;
                        }
                        cout << endl;
                    } else {
                        fileWeightsEtaToOmega >> xValuesReadEtaToOmega[nPtBinsReadEtaToOmega-1];
                        for (Int_t i = 0; i < nMeasSetEtaToOmega; i++){
                            fileWeightsEtaToOmega >> weightsReadEtaToOmega[availableMeasEtaToOmega[i]][nPtBinsReadEtaToOmega-1] ;
                        }
                        cout << "read: "<<  nPtBinsReadEtaToOmega << "\t"<< xValuesReadEtaToOmega[nPtBinsReadEtaToOmega-1] << "\t" ;
                        for (Int_t i = 0; i < nMeasSetEtaToOmega; i++){
                            cout << weightsReadEtaToOmega[availableMeasEtaToOmega[i]][nPtBinsReadEtaToOmega-1] << "\t";
                        }
                        cout << endl;
                    }
                    nPtBinsReadEtaToOmega++;
                }
                nPtBinsReadEtaToOmega = nPtBinsReadEtaToOmega-2 ;
                fileWeightsEtaToOmega.close();

                TGraph* graphWeightsEtaToOmega[nMeasSetEtaToOmega];
                for (Int_t i = 0; i < nMeasSetEtaToOmega; i++){
                    graphWeightsEtaToOmega[availableMeasEtaToOmega[i]] = new TGraph(nPtBinsReadEtaToOmega,xValuesReadEtaToOmega,weightsReadEtaToOmega[availableMeasEtaToOmega[i]]);
                    Int_t bin = 0;
                    for (Int_t n = 0; n< nPtBinsReadEtaToOmega; n++){
                        if (graphWeightsEtaToOmega[availableMeasEtaToOmega[i]]->GetY()[bin] == 0) graphWeightsEtaToOmega[availableMeasEtaToOmega[i]]->RemovePoint(bin);
                        else bin++;
                    }
                }


                if(doBinShiftForEtaToOmega){
                    // Calculate relative sys error weighted
                    if (sysAvailSingleEtaToOmega[0]){
                        for (Int_t k = 0; k< nRelSysErrEtaToOmegaSources ; k++ ){
                            graphRelSysErrEtaToOmegaSourceWeighted[k]      = CalculateWeightedQuantity(   graphOrderedRelSysErrEtaToOmegaSource[k],
                                                                                                        graphWeightsEtaToOmega,
                                                                                                        binningEta,  maxNAllowedEta,
                                                                                                        MaxNumberOfFiles
                                                                                            );
                            if (!graphRelSysErrEtaToOmegaSourceWeighted[k]){
                                cout << "Aborted in CalculateWeightedQuantity for " << endl;
                                return;
                            } else {
                                graphRelSysErrEtaToOmegaSourceWeighted[k]->SetName(Form("RelSysErrEtaToOmegaSourceWeighted%s", ((TString)ptSysDetail[0][0].at(k+1)).Data()));
                            }
                        }

                    }


                    cout << __LINE__ << endl;
                    cout << "combine EtaToOmegaBinShift" << endl;
                    graphEtaToOmegaBinShiftWeighted                     = CalculateWeightedQuantity(  graphEtaToOmegaBinShift,
                                                                                                    graphWeightsEtaToOmega,
                                                                                                    binningEta,  maxNAllowedEta,
                                                                                                    MaxNumberOfFiles
                                                                                                  );
                  if (!graphEtaToOmegaBinShiftWeighted){
                      cout << "Aborted in CalculateWeightedQuantity" << endl;
                      return;
                  }
                  while (graphEtaToOmegaBinShiftWeighted->GetY()[0] == -10000) graphEtaToOmegaBinShiftWeighted->RemovePoint(0);
                }

                //***************************************************************************************************************
                //************************************Plotting binshift corrections *********************************************
                //***************************************************************************************************************

                if(doBinShiftForEtaToOmega){
                  canvasEffi->cd();
                  canvasEffi->SetLeftMargin(0.1);
                  canvasEffi->SetBottomMargin(0.1);
                  canvasEffi->SetLogy(0);
                  canvasEffi->SetLogx(1);

                  TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0., 100.);
                  SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (Y) / no shift",
                                          0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.1, 1.2);
                  histoBinShift->GetXaxis()->SetMoreLogLabels();
                  histoBinShift->GetXaxis()->SetRangeUser(minPtGlobalEta,maxPtGlobalEta);
                  if(optionEnergy.CompareTo("8TeV")==0 && mode==4) histoBinShift->GetXaxis()->SetRangeUser(minPtGlobalEta,20.);
                  histoBinShift->GetYaxis()->SetRangeUser(0.95,1.05);
                  histoBinShift->DrawCopy();

                  TLegend* legendBinShift3 = GetAndSetLegend2(0.62, 0.13, 0.95, 0.13+(1.05*nrOfTrigToBeComb/2*0.85*textSizeSpectra),28);
                  legendBinShift3->SetNColumns(2);

                  for (Int_t j = 0; j< graphEtaToOmegaBinShiftWeighted->GetN(); j++){
                    graphEtaToOmegaBinShiftWeighted->GetEYlow()[j]=0;
                    graphEtaToOmegaBinShiftWeighted->GetEYhigh()[j]=0;
                  }
                  DrawGammaSetMarkerTGraphAsym(graphEtaToOmegaBinShiftWeighted, 20, 1, kGray+2, kGray+2);
                  graphEtaToOmegaBinShiftWeighted->Draw("p same");
        //          if(optionEnergy.CompareTo("8TeV")==0 && mode==4){
        //            TH1D* histoCorrectedYieldEtaBinShiftTEMP = (TH1D*)histoCorrectedYieldOmegaEtaBinBinShift[i]->Clone(Form("OmegaEtaBinning_%i", i));
        //            for(Int_t iB=1; iB<=histoCorrectedYieldOmegaEtaBinBinShift[i]->GetNbinsX(); iB++){histoCorrectedYieldEtaBinShiftTEMP->SetBinContent(iB,histoCorrectedYieldEtaBinShift[i]->GetBinContent(iB));}
        //            histoCorrectedYieldEtaBinShiftTEMP->Divide(histoCorrectedYieldEtaBinShiftTEMP,histoCorrectedYieldOmegaEtaBinBinShift[i],1.,1.,"B");
        //            histoCorrectedYieldEtaBinShift[i] = histoCorrectedYieldEtaBinShiftTEMP;
        //          }else{
        //            histoCorrectedYieldEtaBinShift[i]->Divide(histoCorrectedYieldEtaBinShift[i],histoCorrectedYieldOmegaEtaBinBinShift[i],1.,1.,"B");
        //          }
        //          DrawGammaSetMarker(histoCorrectedYieldEtaBinShift[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
        //          histoCorrectedYieldEtaBinShift[i]->DrawCopy("hist p same");
        //          legendBinShift3->AddEntry(histoCorrectedYieldEtaBinShift[i],triggerNameLabel[i].Data(),"p");
        //          legendBinShift3->Draw();

                  labelEnergyEffi->Draw();
                  TLatex *labelBinShiftEtaToOmega = new TLatex(0.62, maxYLegendEffi+0.02+0.99*textSizeSpectra*0.85,"#eta/#omega");
                  SetStyleTLatex( labelBinShiftEtaToOmega, 0.85*textSizeSpectra,4);
                  labelBinShiftEtaToOmega->Draw();
                  labelDetProcEffi->Draw();

                  canvasEffi->RedrawAxis();
                  canvasEffi->Update();
                  canvasEffi->SaveAs(Form("%s/EtaToOmega_%s_CombinedBinShiftCorrection.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
                  canvasEffi->SetLogy(1);
                  canvasEffi->SetLogx(0);
                  canvasEffi->SetLeftMargin(0.09);
                  canvasEffi->SetBottomMargin(0.08);
                }
                delete canvasEffi;

                //  **********************************************************************************************************************
                //  ******************************************* Plotting weights for eta/Omega *********************************************
                //  **********************************************************************************************************************
                Int_t textSizeLabelsPixel = 900*0.04;

                TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);// gives the page size
                DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);

                TH2F * histo2DWeightsEtaToOmega;
                histo2DWeightsEtaToOmega = new TH2F("histo2DWeightsEtaToOmega","histo2DWeightsEtaToOmega",11000,0.,maxPtGlobalEta,1000,-0.5,1.3);
                SetStyleHistoTH2ForGraphs(histo2DWeightsEtaToOmega, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
                histo2DWeightsEtaToOmega->Draw("copy");

                    TLegend* legendWeightsEtaToOmega = GetAndSetLegend2(0.12, 0.14, 0.55, 0.14+(0.035*nMeasSetEtaToOmega/2*1.35), 32);
                    legendWeightsEtaToOmega->SetNColumns(2);
                    for (Int_t i = 0; i < nMeasSetEtaToOmega; i++){
                        DrawGammaSetMarkerTGraph(graphWeightsEtaToOmega[availableMeasEtaToOmega[i]], markerTriggWeighted[availableMeasEtaToOmega[i]], sizeTrigg[availableMeasEtaToOmega[i]],
                                                colorTriggWeighted[availableMeasEtaToOmega[i]], colorTriggWeighted[availableMeasEtaToOmega[i]]);
                        graphWeightsEtaToOmega[availableMeasEtaToOmega[i]]->Draw("p,same,e1");
                        legendWeightsEtaToOmega->AddEntry(graphWeightsEtaToOmega[availableMeasEtaToOmega[i]],nameTriggerWeighted[availableMeasEtaToOmega[i]],"p");
                    }
                    legendWeightsEtaToOmega->Draw();

                    TLatex *labelWeightsEnergy = new TLatex(0.7,0.24,collisionSystem.Data());
                    SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
                    labelWeightsEnergy->SetTextFont(43);
                    labelWeightsEnergy->Draw();
                    TLatex *labelWeightsEtaToOmega = new TLatex(0.7,0.20,"#eta/#omega");
                    SetStyleTLatex( labelWeightsEtaToOmega, 0.85*textSizeLabelsPixel,4);
                    labelWeightsEtaToOmega->SetTextFont(43);
                    labelWeightsEtaToOmega->Draw();
                    TLatex *labelDetProcWeights    = new TLatex(0.7, 0.16,detectionProcess.Data());
                    SetStyleTLatex( labelDetProcWeights, 0.85*textSizeLabelsPixel,4);
                    labelDetProcWeights->SetTextFont(43);
                    labelDetProcWeights->Draw();

            //      DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
                    DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
                    DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
                    DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
                    DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);

                canvasWeights->SaveAs(Form("%s/%s_WeightsEtaToOmegaTriggers.%s",outputDir.Data(), isMC.Data(), suffix.Data()));
                delete canvasWeights;

                // Calculating relative error for eta/Omega
                for (Int_t i = 0; i < 12; i++){
                    if (histoStatEtaToOmega[i])
                        histoRelStatEtaToOmega[i]      = CalculateRelErrUpTH1D( histoStatEtaToOmega[i], Form("relativeStatErrorEtaToOmega_%s", nameTriggerWeighted[i].Data()));
                    if (graphSystEtaToOmega[i])
                        graphRelSystEtaToOmega[i]      = CalculateRelErrUpAsymmGraph( graphSystEtaToOmega[i], Form("relativeSysErrorEtaToOmega_%s", nameTriggerWeighted[i].Data()));
                }
                TGraphAsymmErrors* graphRelErrorEtaToOmegaTot        = CalculateRelErrUpAsymmGraph( graphEtaToOmegaWeightedAverageTot, "relativeTotalErrorEtaToOmega");
                TGraphAsymmErrors* graphRelErrorEtaToOmegaStat       = CalculateRelErrUpAsymmGraph( graphEtaToOmegaWeightedAverageStat, "relativeStatErrorEtaToOmega");
                TGraphAsymmErrors* graphRelErrorEtaToOmegaSys        = CalculateRelErrUpAsymmGraph( graphEtaToOmegaWeightedAverageSys, "relativeSysErrorEtaToOmega");

                const char *SysErrDatnameMeanSingleErrCheck = Form("%s/SystematicErrorAveragedSingle%s_EtaToOmega_%s_Check.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
                fstream SysErrDatAverSingleCheck;
                SysErrDatAverSingleCheck.precision(4);
                cout << SysErrDatnameMeanSingleErrCheck << endl;
                if(sysAvailSingleEtaToOmega[0]){
                  SysErrDatAverSingleCheck.open(SysErrDatnameMeanSingleErrCheck, ios::out);
                  SysErrDatAverSingleCheck << "pt \t Stat err \t sys err \t tot err " << endl;
                  for (Int_t i = 0; i < graphRelErrorEtaToOmegaTot->GetN(); i++){
                      if (graphRelErrorEtaToOmegaStat->GetY()[i] > 0) SysErrDatAverSingleCheck << graphRelErrorEtaToOmegaStat->GetX()[i] << "\t" << graphRelErrorEtaToOmegaStat->GetY()[i] <<"\t" << graphRelErrorEtaToOmegaSys->GetY()[i] <<  "\t" << graphRelErrorEtaToOmegaTot->GetY()[i] << endl;
                  }
                  SysErrDatAverSingleCheck << endl;
                  SysErrDatAverSingleCheck.close();
                }


                // plot sys relative errors for individual triggers
                TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);

                TH2F * histo2DRelSysErr;
                histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.,maxPtGlobalEta,1000,0,60.5);
                SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
                histo2DRelSysErr->Draw("copy");
                    TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaToOmega/2), 0.95, 0.92, 32);
                    legendRelSysErr->SetNColumns(2);
                    for (Int_t i = 0; i < nMeasSetEtaToOmega; i++){
                        DrawGammaSetMarkerTGraph(graphRelSystEtaToOmega[availableMeasEtaToOmega[i]], markerTriggWeighted[availableMeasEtaToOmega[i]], sizeTrigg[availableMeasEtaToOmega[i]],
                                                colorTriggWeighted[availableMeasEtaToOmega[i]], colorTriggWeighted[availableMeasEtaToOmega[i]]);
                        graphRelSystEtaToOmega[availableMeasEtaToOmega[i]]->Draw("p,same,z");
                        legendRelSysErr->AddEntry(graphRelSystEtaToOmega[availableMeasEtaToOmega[i]],nameTriggerWeighted[availableMeasEtaToOmega[i]],"p");
                    }
                    legendRelSysErr->Draw();

                    TLatex *labelRelErrEnergy    = new TLatex(0.15,0.89,collisionSystem.Data());
                    SetStyleTLatex( labelRelErrEnergy, 0.85*textSizeLabelsPixel,4);
                    labelRelErrEnergy->SetTextFont(43);
                    labelRelErrEnergy->Draw();
                    TLatex *labelRelErrEtaToOmega       = new TLatex(0.15,0.85,"#eta/#omega");
                    SetStyleTLatex( labelRelErrEtaToOmega, 0.85*textSizeLabelsPixel,4);
                    labelRelErrEtaToOmega->SetTextFont(43);
                    labelRelErrEtaToOmega->Draw();
                    TLatex *labelDetProcRelErr    = new TLatex(0.15, 0.81,detectionProcess.Data());
                    SetStyleTLatex( labelDetProcRelErr, 0.85*textSizeLabelsPixel,4);
                    labelDetProcRelErr->SetTextFont(43);
                    labelDetProcRelErr->Draw();

                canvasRelSysErr->SaveAs(Form("%s/EtaToOmega_RelSysErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));

                // plot stat relative errors for individual triggers
                TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);

                TH2F * histo2DRelStatErr;
                histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.,maxPtGlobalEta,1000,0,60.5);
                SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
                histo2DRelStatErr->Draw("copy");
                    TLegend* legendRelStatErr       = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaToOmega/2), 0.95, 0.92, 32);
                    legendRelStatErr->SetNColumns(2);
                    for (Int_t i = 0; i < nMeasSetEtaToOmega; i++){
                        if (availableMeasEtaToOmega[i] == 0 && histoRelStatEtaToOmega[availableMeasEtaToOmega[i]] && mode == 2){
                            TGraphAsymmErrors* dummyGraph = new TGraphAsymmErrors(histoRelStatEtaToOmega[availableMeasEtaToOmega[i]]);
                            dummyGraph->Print();
                            DrawGammaSetMarkerTGraph(dummyGraph, markerTriggWeighted[availableMeasEtaToOmega[i]], sizeTrigg[availableMeasEtaToOmega[i]],
                                                colorTriggWeighted[availableMeasEtaToOmega[i]], colorTriggWeighted[availableMeasEtaToOmega[i]]);
                            dummyGraph->Draw("pX,same");
                            legendRelSysErr->AddEntry(dummyGraph,nameTriggerWeighted[availableMeasEtaToOmega[i]],"p");

        //                     for (Int_t j = 1; j < histoRelStatEtaToOmega[availableMeasEtaToOmega[i]]->GetNbinsX()+1; j++){
        //                        cout << j << ": " << histoRelStatEtaToOmega[availableMeasEtaToOmega[i]]->GetBinContent(j) << endl;
        //                     }
                        } else {
                            DrawGammaSetMarker(histoRelStatEtaToOmega[availableMeasEtaToOmega[i]],markerTriggWeighted[availableMeasEtaToOmega[i]], sizeTrigg[availableMeasEtaToOmega[i]],
                                                colorTriggWeighted[availableMeasEtaToOmega[i]], colorTriggWeighted[availableMeasEtaToOmega[i]]);
                            histoRelStatEtaToOmega[availableMeasEtaToOmega[i]]->Draw("p,same,z");
                            legendRelStatErr->AddEntry(histoRelStatEtaToOmega[availableMeasEtaToOmega[i]],nameTriggerWeighted[availableMeasEtaToOmega[i]],"p");
                        }
                    }
                    legendRelStatErr->Draw();

                    labelRelErrEnergy->Draw();
                    labelRelErrEtaToOmega->Draw();
                    labelDetProcRelErr->Draw();

                canvasRelStatErr->SaveAs(Form("%s/EtaToOmega_RelStatErr_SingleMeas.%s",outputDir.Data(),suffix.Data()));

                // plot full error for final result decomposed
                TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);

                TH2F * histo2DRelTotErrEtaToOmega;
                histo2DRelTotErrEtaToOmega                 = new TH2F("histo2DRelTotErrEtaToOmega","histo2DRelTotErrEtaToOmega",11000,0.,maxPtGlobalEta,1000,0,60.5);
                SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaToOmega, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
                histo2DRelTotErrEtaToOmega->Draw("copy");

                    DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaToOmegaTot, 20, 1.5, kBlack , kBlack);
                    graphRelErrorEtaToOmegaTot->Draw("p,same,z");
                    DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaToOmegaStat, 24, 1.5, kGray+2 , kGray+2);
                    graphRelErrorEtaToOmegaStat->Draw("l,x0,same,e1");
                    DrawGammaSetMarkerTGraphAsym(graphRelErrorEtaToOmegaSys, 24, 1.5, kGray+1 , kGray+1);
                    graphRelErrorEtaToOmegaSys->SetLineStyle(7);
                    graphRelErrorEtaToOmegaSys->Draw("l,x0,same,e1");

                    TLegend* legendRelTotErr2       = GetAndSetLegend2(0.72, 0.92-(0.035*3), 0.9, 0.92, 32);
                    legendRelTotErr2->AddEntry(graphRelErrorEtaToOmegaTot,"tot","p");
                    legendRelTotErr2->AddEntry(graphRelErrorEtaToOmegaStat,"stat","l");
                    legendRelTotErr2->AddEntry(graphRelErrorEtaToOmegaSys,"sys","l");
                    legendRelTotErr2->Draw();

                    labelRelErrEnergy->Draw();
                    labelRelErrEtaToOmega->Draw();
                    labelDetProcRelErr->Draw();

                canvasRelTotErr->SaveAs(Form("%s/EtaToOmega_RelErrorsFulldecomp.%s",outputDir.Data(),suffix.Data()));

                if(optionEnergy.CompareTo("8TeV")==0 && mode==4) maxNAllowedEta += 3;




                if (sysAvailSingleEtaToOmega[0]){
                    for (Int_t k = 0; k < nRelSysErrEtaToOmegaSources; k++){
                        if (graphRelSysErrEtaToOmegaSourceWeighted[k])
                            while (graphRelSysErrEtaToOmegaSourceWeighted[k]->GetY()[0] == -10000 )   graphRelSysErrEtaToOmegaSourceWeighted[k]->RemovePoint(0);
                        else
                            cout << "I don't have a weighted EtaToOmega rel sys err graph for error source: " << k << endl;
                    }
                }


                //  **********************************************************************************************************************
                //  **************************************** Combine+write detailed Systematics ******************************************
                //  **********************************************************************************************************************

                const char *SysErrDatnameMeanSingleErr = Form("%s/SystematicErrorAveragedSingle%s_EtaToOmega_%s.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
                fstream SysErrDatAverSingle;
                SysErrDatAverSingle.precision(4);
                cout << SysErrDatnameMeanSingleErr << endl;
                if(sysAvailSingleEtaToOmega[0] && graphRelSysErrEtaToOmegaSourceWeighted[0]){
                    SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
                    for(Int_t iColumn = 0; iColumn < (Int_t)ptSysDetail[0][0].size(); iColumn++) SysErrDatAverSingle << ptSysDetail[0][0].at(iColumn) << "\t";
                    SysErrDatAverSingle << endl;
                    for (Int_t i = 1; i < nrOfTrigToBeComb; i++){
                        if(!sysAvailSingleEtaToOmega[i]) continue;
                        for(Int_t iCol = 0; iCol < (Int_t)ptSysDetail[i][0].size(); iCol++){
                            if( ((TString)ptSysDetail[i][0].at(iCol)).CompareTo(((TString)ptSysDetail[0][0].at(iCol))) ){
                                cout << "ERROR: Systematic error type at pos " << iCol << " does not agree for " << availableMeasEtaToOmega[i] << " & " << availableMeasEtaToOmega[0] << ", returning!" << endl;
                                return;
                            }
                        }
                    }

                    for(Int_t i=0; i<graphRelSysErrEtaToOmegaSourceWeighted[0]->GetN(); i++){
                        SysErrDatAverSingle << graphRelSysErrEtaToOmegaSourceWeighted[0]->GetX()[i] << "\t";
                        Int_t nColumns = (Int_t)ptSysDetail[0][0].size();
                        for(Int_t iErr=0; iErr<nColumns-1; iErr++)
                            SysErrDatAverSingle << graphRelSysErrEtaToOmegaSourceWeighted[iErr]->GetY()[i] << "\t";
                        SysErrDatAverSingle << endl;

                    }
                }
                SysErrDatAverSingle.close();

                // ***************************************************************************************************
                // ********************* Plot all mean erros separately after smoothing ******************************
                // ***************************************************************************************************
                if(sysAvailSingleEtaToOmega[0] && graphRelSysErrEtaToOmegaSourceWeighted[0]){
                  TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);// gives the page size
                  DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);

                      // create dummy histo
                      TH2D *histo2DNewSysErrMean ;
                      histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 100,0.,maxPtGlobalEta,1000.,-0.5,50.);
                      SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                              1,0.9, 510, 510);
                      histo2DNewSysErrMean->Draw();

                      // Give legend position for plotting
                      Double_t minXLegend     = 0.12;
                      Double_t maxYLegend     = 0.95;
                      Double_t widthLegend    = 0.25;
                      if (nRelSysErrEtaToOmegaSources> 7)
                          widthLegend         = 0.5;
                      Double_t heightLegend   = 1.05* 0.035 * (nRelSysErrEtaToOmegaSources+3);
                      if (nRelSysErrEtaToOmegaSources> 7)
                          heightLegend        = 1.05* 0.035 * (nRelSysErrEtaToOmegaSources/2+1);

                      // create legend
                      TLegend* legendMeanNew = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
                      legendMeanNew->SetMargin(0.1);
                      if (nRelSysErrEtaToOmegaSources> 7) legendMeanNew->SetNColumns(2);

                      for(Int_t i = 0;i< nRelSysErrEtaToOmegaSources-1 ;i++){
                          DrawGammaSetMarkerTGraphAsym(graphRelSysErrEtaToOmegaSourceWeighted[i], GetMarkerStyleSystematics(  ptSysDetail[0][0].at(i+1), mode), 1.,
                                                      GetColorSystematics( ptSysDetail[0][0].at(i+1), mode),GetColorSystematics( ptSysDetail[0][0].at(i+1), mode));
                          graphRelSysErrEtaToOmegaSourceWeighted[i]->Draw("pX0,csame");
                          legendMeanNew->AddEntry(graphRelSysErrEtaToOmegaSourceWeighted[i],GetSystematicsName(ptSysDetail[0][0].at(i+1)),"p");
                      }

                      DrawGammaSetMarkerTGraphAsym(graphRelSysErrEtaToOmegaSourceWeighted[nRelSysErrEtaToOmegaSources-1], 20, 1.,kBlack,kBlack);
                      graphRelSysErrEtaToOmegaSourceWeighted[nRelSysErrEtaToOmegaSources-1]->Draw("p,csame");
                      legendMeanNew->AddEntry(graphRelSysErrEtaToOmegaSourceWeighted[nRelSysErrEtaToOmegaSources-1],"quad. sum.","p");
                      legendMeanNew->Draw();

                      // labeling
                      TLatex *labelEnergySysDetailed = new TLatex(0.7, 0.93,collisionSystem.Data());
                      labelEnergySysDetailed->SetTextAlign(31);
                      SetStyleTLatex( labelEnergySysDetailed, 0.85*textSizeSpectra,4);
                      labelEnergySysDetailed->Draw();

                      TLatex *labelEtaToOmegaSysDetailed     = new TLatex(0.7, 0.93-0.99*textSizeSpectra*0.85,"#eta/#omega");
                      labelEtaToOmegaSysDetailed->SetTextAlign(31);
                      SetStyleTLatex( labelEtaToOmegaSysDetailed, 0.85*textSizeSpectra,4);
                      labelEtaToOmegaSysDetailed->Draw();

                      TLatex *labelDetProcSysDetailed = new TLatex(0.7, 0.93-2*0.99*textSizeSpectra*0.85,detectionProcess.Data());
                      labelDetProcSysDetailed->SetTextAlign(31);
                      SetStyleTLatex( labelDetProcSysDetailed, 0.85*textSizeSpectra,4);
                      labelDetProcSysDetailed->Draw();

                  canvasNewSysErrMean->Update();
                  canvasNewSysErrMean->SaveAs(Form("%s/EtaToOmega_SysErrorsSeparatedSourcesReweighted_%s.%s",outputDir.Data(),isMC.Data(),suffix.Data()));
                }

                // delete detailed sys array
                for(Int_t iR=0; iR<nrOfTrigToBeComb; iR++){
                    for(Int_t iB=0; iB<50; iB++) ptSysDetail[iR][iB].clear();
                }

            // if averaging wasn't enabled pick values according to predefined ranges ("cherry picking points")
            } else {
                graphEtaToOmegaWeightedAverageStat        = new TGraphAsymmErrors(nPointFinalEtaToOmega, xValueFinalEtaToOmega, yValueFinalEtaToOmega,
                                                                                         xErrorLowFinalEtaToOmega, xErrorHighFinalEtaToOmega,yErrorLowFinalEtaToOmega, yErrorHighFinalEtaToOmega);
                graphEtaToOmegaWeightedAverageSys         = new TGraphAsymmErrors(nPointFinalEtaToOmega, xValueFinalEtaToOmega, yValueFinalEtaToOmega,
                                                                                        xErrorLowFinalEtaToOmega, xErrorHighFinalEtaToOmega,yErrorSysLowFinalEtaToOmega, yErrorSysHighFinalEtaToOmega);
                TGraphAsymmErrors* graphRelErrorEtaToOmegaStat       = CalculateRelErrUpAsymmGraph( graphEtaToOmegaWeightedAverageStat, "relativeStatErrorEtaToOmega");
                while (graphRelErrorEtaToOmegaStat->GetY()[0] < 0 ) graphRelErrorEtaToOmegaStat->RemovePoint(0);

                TGraphAsymmErrors* graphRelErrorEtaToOmegaSys        = CalculateRelErrUpAsymmGraph( graphEtaToOmegaWeightedAverageSys, "relativeSysErrorEtaToOmega");
                while (graphRelErrorEtaToOmegaSys->GetY()[0] < 0 ) graphRelErrorEtaToOmegaSys->RemovePoint(0);

                const char *SysErrDatnameMeanSingleErrCheck = Form("%s/SystematicErrorAveragedSingle%s_EtaToOmega_%s_Check.dat",outputDir.Data(),sysStringComb.Data(),optionEnergy.Data());
                fstream SysErrDatAverSingleCheck;
                SysErrDatAverSingleCheck.precision(4);
                cout << SysErrDatnameMeanSingleErrCheck << endl;

                SysErrDatAverSingleCheck.open(SysErrDatnameMeanSingleErrCheck, ios::out);
                SysErrDatAverSingleCheck << "pt \t Stat err \t sys err \t tot err " << endl;
                for (Int_t i = 0; i < graphRelErrorEtaToOmegaStat->GetN(); i++){
                    if (graphRelErrorEtaToOmegaStat->GetY()[i] > 0) SysErrDatAverSingleCheck << graphRelErrorEtaToOmegaStat->GetX()[i] << "\t" << graphRelErrorEtaToOmegaStat->GetY()[i] <<"\t" << graphRelErrorEtaToOmegaSys->GetY()[i]  << endl;
                }
                SysErrDatAverSingleCheck << endl;
                SysErrDatAverSingleCheck.close();
            }
            // printing final eta/Omega ratios
            cout << "stat eta/Omega" << endl;
            graphEtaToOmegaWeightedAverageStat->Print();
            cout << "sys eta/Omega" << endl;
            if (graphEtaToOmegaWeightedAverageSys) graphEtaToOmegaWeightedAverageSys->Print();

            if (graphEtaToOmegaWeightedAverageStat){
                while (graphEtaToOmegaWeightedAverageStat->GetX()[0]< minPtGlobalEta){
                    graphEtaToOmegaWeightedAverageStat->RemovePoint(0);
                }
            }

            if (graphEtaToOmegaWeightedAverageSys){
                while (graphEtaToOmegaWeightedAverageSys->GetX()[0]< minPtGlobalEta){
                    graphEtaToOmegaWeightedAverageSys->RemovePoint(0);
                }
            }


            //***************************************************************************************************************
            //***************************Plotting individual eta/Omega ratios full range **************************************
            //***************************************************************************************************************

            Size_t textSizeEtaToOmega = 0.04;
            TCanvas* canvasEtatoOmegacombo = new TCanvas("canvasEtatoOmegacombo","",200,10,1200,1100);
            DrawGammaCanvasSettings( canvasEtatoOmegacombo, 0.09, 0.017, 0.02, 0.1);

            TH2F * histo2DEtatoOmegacombo = new TH2F("histo2DEtatoOmegacombo","histo2DEtatoOmegacombo",11000,0,maxPtGlobalEta,1000,0.01,1.2);
            SetStyleHistoTH2ForGraphs(histo2DEtatoOmegacombo, "#it{p}_{T} (GeV/#it{c})","#eta/#omega",0.035,0.04, 0.035,0.04, 1.2,1.);
            histo2DEtatoOmegacombo->GetYaxis()->SetRangeUser(0.01,1.05);

            histo2DEtatoOmegacombo->DrawCopy();
            DrawGammaLines(0., maxPtGlobalEta , 0.45, 0.45, 0.3, kGray+2);
            DrawGammaLines(0., maxPtGlobalEta , 0.4, 0.4, 0.3, kGray, 7);
            DrawGammaLines(0., maxPtGlobalEta , 0.5, 0.5, 0.3, kGray, 7);


            TLegend* legendEtaToOmega = GetAndSetLegend2(0.32, 0.14, 0.6, 0.14+(0.035*(nrOfTrigToBeComb/2+1)*0.9), 32);
            legendEtaToOmega->SetNColumns(2);

            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if (foundOmegaEtaBinFile[i] && !maskedFullyEta[i] && !maskedFullyOmega[i] ){
                    DrawGammaSetMarker(histoEtaToOmega[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    if (graphsEtaToOmegaSysRemoved0[i]){
                        DrawGammaSetMarkerTGraphAsym(graphsEtaToOmegaSysRemoved0[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
                        graphsEtaToOmegaSysRemoved0[i]->Draw("p,E2,same");
                    }
                    histoEtaToOmega[i]->Draw("e1,same");
                    legendEtaToOmega->AddEntry(histoEtaToOmega[i],triggerNameLabel[i],"p");
                }
            }
            legendEtaToOmega->Draw();

            TLatex *labelEnergyEtaToOmega = new TLatex(0.65, 0.15+(2*textSizeEtaToOmega*0.85),collisionSystem.Data());
            SetStyleTLatex( labelEnergyEtaToOmega, 0.85*textSizeEtaToOmega,4);
            labelEnergyEtaToOmega->Draw();

            TLatex *labelOmegaEtaToOmega = new TLatex(0.65, 0.15+textSizeEtaToOmega*0.85,"#eta/#omega");
            SetStyleTLatex( labelOmegaEtaToOmega, 0.85*textSizeEtaToOmega,4);
            labelOmegaEtaToOmega->Draw();

            TLatex *labelDetProcEtaToOmega = new TLatex(0.65, 0.15,detectionProcess.Data());
            SetStyleTLatex( labelDetProcEtaToOmega, 0.85*textSizeEtaToOmega,4);
            labelDetProcEtaToOmega->Draw();

            histo2DEtatoOmegacombo->Draw("axis,same");

            canvasEtatoOmegacombo->Update();
            canvasEtatoOmegacombo->SaveAs(Form("%s/EtaToOmega%s_%s_diffTriggers.%s",outputDir.Data(), addNameBinshift.Data(), isMC.Data(), suffix.Data()));

            //***************************************************************************************************************
            //***************************Plotting individual eta/Omega ratios used range with full ****************************
            //***************************************************************************************************************

            histo2DEtatoOmegacombo->DrawCopy();
            DrawGammaLines(0., maxPtGlobalEta , 0.45, 0.45, 0.3, kGray+2);
            DrawGammaLines(0., maxPtGlobalEta , 0.4, 0.4, 0.3, kGray, 7);
            DrawGammaLines(0., maxPtGlobalEta , 0.5, 0.5, 0.3, kGray, 7);


            for (Int_t i = 0; i< nrOfTrigToBeComb; i++){
                if (foundOmegaEtaBinFile[i] && !maskedFullyEta[i] && !maskedFullyOmega[i]){
                    if (graphsEtaToOmegaShrunk[i]) DrawGammaSetMarkerTGraphAsym(graphsEtaToOmegaShrunk[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i]);
                    if (graphsEtaToOmegaSysShrunk[i]){
                        DrawGammaSetMarkerTGraphAsym(graphsEtaToOmegaSysShrunk[i], markerTrigg[i], sizeTrigg[i], colorTrigg[i], colorTrigg[i], 1, kTRUE);
                        graphsEtaToOmegaSysShrunk[i]->Draw("p,E2,same");
                    }
                    if (graphsEtaToOmegaShrunk[i]) graphsEtaToOmegaShrunk[i]->Draw("e1,same");
                }
            }

            if (graphEtaToOmegaWeightedAverageSys){
                DrawGammaSetMarkerTGraphAsym(graphEtaToOmegaWeightedAverageSys,  24, 2, kRed , kRed, 1, kTRUE);
                graphEtaToOmegaWeightedAverageSys->Draw("p,E2,same");
            }
            if (graphEtaToOmegaWeightedAverageStat){
                DrawGammaSetMarkerTGraphAsym(graphEtaToOmegaWeightedAverageStat,  24, 2, kRed , kRed);
                graphEtaToOmegaWeightedAverageStat->Draw("p,E,same");
                legendEtaToOmega->AddEntry(graphEtaToOmegaWeightedAverageStat, "final","p");
            }
            legendEtaToOmega->Draw();

            labelEnergyEtaToOmega->Draw();
            labelOmegaEtaToOmega->Draw();
            labelDetProcEtaToOmega->Draw();

            histo2DEtatoOmegacombo->Draw("axis,same");

            canvasEtatoOmegacombo->Update();
            canvasEtatoOmegacombo->SaveAs(Form("%s/EtaToOmega%s_%s_diffTriggers_used.%s",outputDir.Data(), addNameBinshift.Data(), isMC.Data(), suffix.Data()));

            //***************************************************************************************************************
            //**************************************Plotting final eta/Omega ratio ********************************************
            //***************************************************************************************************************

            histo2DEtatoOmegacombo->DrawCopy();
            DrawGammaLines(0., maxPtGlobalEta , 0.45, 0.45, 0.3, kGray+2);
            DrawGammaLines(0., maxPtGlobalEta , 0.4, 0.4, 0.3, kGray, 7);
            DrawGammaLines(0., maxPtGlobalEta , 0.5, 0.5, 0.3, kGray, 7);

            if (graphEtaToOmegaWeightedAverageSys){
                DrawGammaSetMarkerTGraphAsym(graphEtaToOmegaWeightedAverageSys,  24, 2, kGray+1 , kGray+1, 1, kTRUE);
                graphEtaToOmegaWeightedAverageSys->Draw("p,E2,same");
            }
            if (graphEtaToOmegaWeightedAverageStat){
                DrawGammaSetMarkerTGraphAsym(graphEtaToOmegaWeightedAverageStat,  24, 2, kBlack , kBlack);
                graphEtaToOmegaWeightedAverageStat->Draw("p,E,same");
            }

            labelEnergyEtaToOmega->Draw();
            labelOmegaEtaToOmega->Draw();
            labelDetProcEtaToOmega->Draw();

            histo2DEtatoOmegacombo->Draw("axis,same");

            canvasEtatoOmegacombo->Update();
            canvasEtatoOmegacombo->SaveAs(Form("%s/EtaToOmega%s_%s_Final.%s",outputDir.Data(), addNameBinshift.Data(), isMC.Data(), suffix.Data()));
        }
    }


    fileFitsOutput << WriteParameterToFile(fitInvYieldOmega) << endl;
    fileFitsOutput.close();

    //*********************************************************************************************************
    //********************** ComparisonFile Output ************************************************************
    //*********************************************************************************************************
    TString system              = "PCM";
    if (mode == 2) system       = "PCM-EMCAL";
    if (mode == 3) system       = "PCM-PHOS";
    if (mode == 4) system       = "EMCAL-EMCAL";
    if (mode == 5) system       = "PHOS-PHOS";
    if (mode == 11) system      = "PHOS-merged";
    if (mode == 40 || mode == 60) system      = "PCM";
    if (mode == 41 || mode == 61) system      = "PCMEMC";
    if (mode == 42 || mode == 62) system      = "PCMPHOS";
    if (mode == 44 || mode == 64) system      = "EMC";
    if (mode == 45 || mode == 65) system      = "PHOS";

    cout << "starting to write out data" << endl;
    TGraphAsymmErrors* graphInvXSectionWeightedAverageOmegaStat   = NULL;
    TH1D* histoInvXSectionWeightedAverageOmegaStat                = NULL;
    TH1D* histoInvYieldWeightedAverageOmegaStat                   = NULL;
    TGraphAsymmErrors* graphInvXSectionWeightedAverageOmegaSys    = NULL;
    TGraphAsymmErrors* graphInvXSectionWeightedAverageEtaStat   = NULL;
    TH1D* histoInvXSectionWeightedAverageEtaStat                = NULL;
    TH1D* histoInvYieldWeightedAverageEtaStat                   = NULL;
    TGraphAsymmErrors* graphInvXSectionWeightedAverageEtaSys    = NULL;
    TH1D* histoEtaToOmegaWeightedAverageStat                      = NULL;
    TH1D* histoEtaToOmegaExtendedUsingFit                         = NULL;


    cout << "xSection = " << xSection << endl;

    if (graphCorrectedYieldWeightedAverageOmegaStat){
        histoInvYieldWeightedAverageOmegaStat                  = new TH1D("histoInvYieldWeightedAverageOmegaStat", "", maxNBinsOmega, binningOmega);
        Int_t firstBinOmega = 1;
        while (histoInvYieldWeightedAverageOmegaStat->GetBinCenter(firstBinOmega) < graphCorrectedYieldWeightedAverageOmegaStat->GetX()[0]){
            histoInvYieldWeightedAverageOmegaStat->SetBinContent(firstBinOmega, 0);
            histoInvYieldWeightedAverageOmegaStat->SetBinError(firstBinOmega, 0);
            firstBinOmega++;
        }
        for (Int_t i = 0; i < graphCorrectedYieldWeightedAverageOmegaStat->GetN(); i++){
            histoInvYieldWeightedAverageOmegaStat->SetBinContent(i+firstBinOmega, graphCorrectedYieldWeightedAverageOmegaStat->GetY()[i]);
            histoInvYieldWeightedAverageOmegaStat->SetBinError(i+firstBinOmega, graphCorrectedYieldWeightedAverageOmegaStat->GetEYlow()[i]);
        }
        if (xSection != 0){
            graphInvXSectionWeightedAverageOmegaStat                  = ScaleGraph(graphCorrectedYieldWeightedAverageOmegaStat,xSection*recalcBarn);
            histoInvXSectionWeightedAverageOmegaStat                  = new TH1D("histoInvXSectionWeightedAverageOmegaStat", "", maxNBinsOmega, binningOmega);
            for (Int_t i = 0; i < graphInvXSectionWeightedAverageOmegaStat->GetN(); i++){
                histoInvXSectionWeightedAverageOmegaStat->SetBinContent(i+firstBinOmega, graphInvXSectionWeightedAverageOmegaStat->GetY()[i]);
                histoInvXSectionWeightedAverageOmegaStat->SetBinError(i+firstBinOmega, graphInvXSectionWeightedAverageOmegaStat->GetEYlow()[i]);
            }
            graphInvXSectionWeightedAverageOmegaStat->Print();
            if (graphCorrectedYieldWeightedAverageOmegaSys)
                graphInvXSectionWeightedAverageOmegaSys                   = ScaleGraph(graphCorrectedYieldWeightedAverageOmegaSys,xSection*recalcBarn);
            graphInvXSectionWeightedAverageOmegaSys->Print();
        }
    }


    if (enableEta && mode != 10){
        if (graphCorrectedYieldWeightedAverageEtaStat){
            histoInvYieldWeightedAverageEtaStat                  = new TH1D("histoInvYieldWeightedAverageEtaStat", "", maxNAllowedEta, binningEta);
            Int_t firstBinEta = 1;
            while (histoInvYieldWeightedAverageEtaStat->GetBinCenter(firstBinEta) < graphCorrectedYieldWeightedAverageEtaStat->GetX()[0]){
                histoInvYieldWeightedAverageEtaStat->SetBinContent(firstBinEta, 0);
                histoInvYieldWeightedAverageEtaStat->SetBinError(firstBinEta, 0);
                firstBinEta++;
            }
            for (Int_t i = 0; i < graphCorrectedYieldWeightedAverageEtaStat->GetN(); i++){
                histoInvYieldWeightedAverageEtaStat->SetBinContent(i+firstBinEta, graphCorrectedYieldWeightedAverageEtaStat->GetY()[i]);
                histoInvYieldWeightedAverageEtaStat->SetBinError(i+firstBinEta, graphCorrectedYieldWeightedAverageEtaStat->GetEYlow()[i]);
            }
            if (xSection != 0){
                graphInvXSectionWeightedAverageEtaStat                  = ScaleGraph(graphCorrectedYieldWeightedAverageEtaStat,xSection*recalcBarn);
                histoInvXSectionWeightedAverageEtaStat                  = new TH1D("histoInvXSectionWeightedAverageEtaStat", "", maxNAllowedEta, binningEta);
                for (Int_t i = 0; i < graphInvXSectionWeightedAverageEtaStat->GetN(); i++){
                    histoInvXSectionWeightedAverageEtaStat->SetBinContent(i+firstBinEta, graphInvXSectionWeightedAverageEtaStat->GetY()[i]);
                    histoInvXSectionWeightedAverageEtaStat->SetBinError(i+firstBinEta, graphInvXSectionWeightedAverageEtaStat->GetEYlow()[i]);
                }
                graphInvXSectionWeightedAverageEtaStat->Print();
                if (graphCorrectedYieldWeightedAverageEtaSys)
                    graphInvXSectionWeightedAverageEtaSys                   = ScaleGraph(graphCorrectedYieldWeightedAverageEtaSys,xSection*recalcBarn);
                graphInvXSectionWeightedAverageEtaSys->Print();
            }
        }

        if (doEtaToOmega && graphEtaToOmegaWeightedAverageStat){
            histoEtaToOmegaWeightedAverageStat                    = new TH1D("histoEtaToOmegaWeightedAverageStat", "", maxNAllowedEta, binningEta);
            Int_t firstBinEtaToOmega = 1;
            while (histoEtaToOmegaWeightedAverageStat->GetBinCenter(firstBinEtaToOmega) < graphEtaToOmegaWeightedAverageStat->GetX()[0]){
                histoEtaToOmegaWeightedAverageStat->SetBinContent(firstBinEtaToOmega, 0);
                histoEtaToOmegaWeightedAverageStat->SetBinError(firstBinEtaToOmega, 0);

                firstBinEtaToOmega++;
            }
            for (Int_t i = 0; i < graphEtaToOmegaWeightedAverageStat->GetN(); i++){
                histoEtaToOmegaWeightedAverageStat->SetBinContent(i+firstBinEtaToOmega, graphEtaToOmegaWeightedAverageStat->GetY()[i]);
                histoEtaToOmegaWeightedAverageStat->SetBinError(i+firstBinEtaToOmega, graphEtaToOmegaWeightedAverageStat->GetEYlow()[i]);
            }

            if(optionEnergy.CompareTo("8TeV") == 0 ){
              histoEtaToOmegaExtendedUsingFit = new TH1D("EtaToOmega_extendedFit","EtaToOmega_extendedFit",maxNAllowedEta,binningEta);
              Int_t i = 0;
              for (; i < graphEtaToOmegaWeightedAverageStat->GetN(); i++){
              //for (; graphEtaToOmegaWeightedAverageStat->GetX()[i] < 10.; i++){
                histoEtaToOmegaExtendedUsingFit->SetBinContent(i+firstBinEtaToOmega, graphEtaToOmegaWeightedAverageStat->GetY()[i]);
                histoEtaToOmegaExtendedUsingFit->SetBinError(i+firstBinEtaToOmega, graphEtaToOmegaWeightedAverageStat->GetEYlow()[i]);
                cout << i << ", " << firstBinEtaToOmega << ", " << graphEtaToOmegaWeightedAverageStat->GetY()[i] << ", " << graphEtaToOmegaWeightedAverageStat->GetEYlow()[i] << endl;
              }
              for (; i+firstBinEtaToOmega <= maxNAllowedEta; i++){
                histoEtaToOmegaExtendedUsingFit->SetBinContent(i+firstBinEtaToOmega, graphInvXSectionWeightedAverageEtaStat->GetY()[i]/fitBinShiftOmegaTCM->Eval((binningEta[i+firstBinEtaToOmega]+binningEta[i+firstBinEtaToOmega-1])/2));
                histoEtaToOmegaExtendedUsingFit->SetBinError(i+firstBinEtaToOmega, graphInvXSectionWeightedAverageEtaStat->GetEYlow()[i]/graphInvXSectionWeightedAverageEtaStat->GetY()[i]);
                cout << i << ", " << firstBinEtaToOmega << ", " << graphInvXSectionWeightedAverageEtaStat->GetY()[i] << ", " << fitBinShiftOmegaTCM->Eval((binningEta[i+firstBinEtaToOmega]+binningEta[i+firstBinEtaToOmega-1])/2) << endl;
              }
            }
        }

    }

    cout << "Writing output file" << endl;
    TString fileNameOutputComp  = Form("%s/%s_%sResultsFullCorrection_PP.root",outputDir.Data(),isMC.Data(),system.Data());
    if (optionEnergy.Contains("pPb"))
        fileNameOutputComp      = Form("%s/%s_%sResultsFullCorrection_pPb.root",outputDir.Data(),isMC.Data(),system.Data());
    else if (optionEnergy.Contains("PbPb"))
        fileNameOutputComp      = Form("%s/%s_%sResultsFullCorrection_PbPb.root",outputDir.Data(),isMC.Data(),system.Data());

    TFile* fileOutputForComparisonFullyCorrected = new TFile(fileNameOutputComp,"UPDATE");
        for (Int_t i=0; i< nrOfTrigToBeComb; i++){
            if (histoEventQualtity[i])                  histoEventQualtity[i]->Write(Form("NEvents_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (mode == 2 || mode == 3 || mode == 4 || mode == 5 || mode == 10 || mode == 11 ){
                if (histoTriggerRejection[i])           histoTriggerRejection[i]->Write(Form("TriggRejectMean_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),TObject::kOverwrite);
                if (histoRatioRawClusterPt[i])          histoRatioRawClusterPt[i]->Write(Form("TriggRejectvsPt_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),TObject::kOverwrite);
                if (histoRatioRawClusterE[i])           histoRatioRawClusterE[i]->Write(Form("TriggRejectvsE_%s_%s",triggerName[i].Data(), triggerName[trigSteps[i][0]].Data()),TObject::kOverwrite);
            }
        }
        fileOutputForComparisonFullyCorrected->mkdir(Form("Omega%s",optionEnergy.Data()));
        TDirectoryFile* directoryOmega = (TDirectoryFile*)fileOutputForComparisonFullyCorrected->Get(Form("Omega%s",optionEnergy.Data()));
        fileOutputForComparisonFullyCorrected->cd(Form("Omega%s",optionEnergy.Data()));
        if (graphCorrectedYieldWeightedAverageOmegaStat)  graphCorrectedYieldWeightedAverageOmegaStat->Write("graphCorrectedYieldOmega",TObject::kOverwrite);
        if (histoInvYieldWeightedAverageOmegaStat)  histoInvYieldWeightedAverageOmegaStat->Write("CorrectedYieldOmega",TObject::kOverwrite);
        if (graphCorrectedYieldWeightedAverageOmegaSys)   graphCorrectedYieldWeightedAverageOmegaSys->Write("OmegaSystError",TObject::kOverwrite);
        if (xSection != 0){
            if (graphInvXSectionWeightedAverageOmegaStat){
                graphInvXSectionWeightedAverageOmegaStat->Write("graphInvCrossSectionOmega",TObject::kOverwrite);
                cout << "graph InvSection Omega stat" << endl;
                graphInvXSectionWeightedAverageOmegaStat->Print();
            }
            if (histoInvXSectionWeightedAverageOmegaStat)     histoInvXSectionWeightedAverageOmegaStat->Write("InvCrossSectionOmega",TObject::kOverwrite);
            if (graphInvXSectionWeightedAverageOmegaSys){
                graphInvXSectionWeightedAverageOmegaSys->Write("InvCrossSectionOmegaSys",TObject::kOverwrite);
                cout << "graph InvSection Omega sys" << endl;
                graphInvXSectionWeightedAverageOmegaSys->Print();
            }
        }
        if (graphAcceptanceOmegaWeighted)                 graphAcceptanceOmegaWeighted->Write("AcceptanceOmega", TObject::kOverwrite);
        if (graphEfficiencyOmegaWeighted)                 graphEfficiencyOmegaWeighted->Write("EfficiencyOmega", TObject::kOverwrite);
        if (graphEffTimesAccOmegaWeighted)                graphEffTimesAccOmegaWeighted->Write("EffTimesAccOmega", TObject::kOverwrite);
        if (graphPurityOmegaWeighted)                     graphPurityOmegaWeighted->Write("PurityOmega", TObject::kOverwrite);
        if (mode != 10){
            if (graphMassOmegaDataWeighted)                   graphMassOmegaDataWeighted->Write("Omega_Mass_data", TObject::kOverwrite);
            if (graphMassOmegaMCWeighted)                     graphMassOmegaMCWeighted->Write("Omega_Mass_MC", TObject::kOverwrite);
            if (graphWidthOmegaDataWeighted)                  graphWidthOmegaDataWeighted->Write("Omega_Width_data", TObject::kOverwrite);
            if (graphWidthOmegaMCWeighted)                    graphWidthOmegaMCWeighted->Write("Omega_Width_MC", TObject::kOverwrite);
        }
        for (Int_t k = 0; k< 4; k++){
            if (graphEffectSecCorrOmegaWeighted[k])       graphEffectSecCorrOmegaWeighted[k]->Write(Form("EffectiveSecondaryOmegaCorrFrom%s",nameSecOmegaPartRead[k].Data()), TObject::kOverwrite);
            if (graphEfficiencySecOmegaWeighted[k])       graphEfficiencySecOmegaWeighted[k]->Write(Form("SecondaryOmegaEfficiencyFrom%s",nameSecOmegaPartRead[k].Data()), TObject::kOverwrite);
        }

        if (sysAvailSingleOmega[0]){
            for (Int_t k = 0; k< nRelSysErrOmegaSources ; k++ ){
                if (graphRelSysErrOmegaSourceWeighted[k]) graphRelSysErrOmegaSourceWeighted[k]->Write(graphRelSysErrOmegaSourceWeighted[k]->GetName(), TObject::kOverwrite);
            }
        }
        for (Int_t i=0; i< nrOfTrigToBeComb; i++){
            cout << "trigger: " << triggerName[i].Data() << endl;
            if (histoAcceptanceOmega[i])                  histoAcceptanceOmega[i]->Write(Form("AcceptanceOmega_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoEfficiencyOmega[i])                  histoEfficiencyOmega[i]->Write(Form("EfficiencyOmega_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoRawYieldOmega[i])                    histoRawYieldOmega[i]->Write(Form("RAWYieldPerEventsOmega_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoMCInputOmega[i])                     histoMCInputOmega[i]->Write(Form("Omega_Input_Reweighted_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (histoEffTimesAccOmega[i])                 histoEffTimesAccOmega[i]->Write(Form("EffTimesAccOmega_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (mode != 10){
                if (histoMassOmegaData[i])                    histoMassOmegaData[i]->Write(Form("Omega_Mass_data_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoMassOmegaMC[i])                      histoMassOmegaMC[i]->Write(Form("Omega_Mass_MC_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoWidthOmegaData[i])                   histoWidthOmegaData[i]->Write(Form("Omega_Width_data_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoWidthOmegaMC[i])                     histoWidthOmegaMC[i]->Write(Form("Omega_Width_MC_%s",triggerName[i].Data()),TObject::kOverwrite);

                if (histoInvMassSig[i])                     histoInvMassSig[i]->Write(Form("Omega_InvMassSig_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoInvMassSigPlusBG[i])               histoInvMassSigPlusBG[i]->Write(Form("Omega_InvMassSigPlusBG_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoInvMassBG[i])                      histoInvMassBG[i]->Write(Form("Omega_InvMassBG_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (fitInvMassSig[i])                       fitInvMassSig[i]->Write(Form("Omega_InvMassSigFit_Example_%s",triggerName[i].Data()),TObject::kOverwrite);
            }
            if (graphsCorrectedYieldRemoved0Omega[i])     graphsCorrectedYieldRemoved0Omega[i]->Write(Form("CorrectedYieldOmega_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (graphsCorrectedYieldSysRemoved0Omega[i])  graphsCorrectedYieldSysRemoved0Omega[i]->Write(Form("OmegaSystError_%s",triggerName[i].Data()),TObject::kOverwrite);
            if (enableTriggerEffOmega[i]){
                if (histoTriggerEffOmega[i])              histoTriggerEffOmega[i]->Write(Form("TriggerEfficiencyOmega_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoEffBaseOmega[i])                 histoEffBaseOmega[i]->Write(Form("EfficiencyBaseOmega_%s",triggerName[i].Data()),TObject::kOverwrite);
            }
            for (Int_t k = 0; k< 4; k++){
                if (histoEffectCorrOmegaFromX[k][i])      histoEffectCorrOmegaFromX[k][i]->Write(Form("EffectiveSecondaryOmegaCorrFrom%s_%s",nameSecOmegaPartRead[k].Data(),triggerName[i].Data()), TObject::kOverwrite);
                if (histoSecEffiOmegaFromX[k][i])         histoSecEffiOmegaFromX[k][i]->Write(Form("SecondaryOmegaEfficiencyFrom%s_%s",nameSecOmegaPartRead[k].Data(),triggerName[i].Data()), TObject::kOverwrite);
            }
        }

        if (enableEta && mode != 10){
            fileOutputForComparisonFullyCorrected->mkdir(Form("Eta%s",optionEnergy.Data()));
            TDirectoryFile* directoryEta = (TDirectoryFile*)fileOutputForComparisonFullyCorrected->Get(Form("Eta%s",optionEnergy.Data()));
            fileOutputForComparisonFullyCorrected->cd(Form("Eta%s",optionEnergy.Data()));
            if (graphCorrectedYieldWeightedAverageEtaStat)  graphCorrectedYieldWeightedAverageEtaStat->Write("graphCorrectedYieldEta",TObject::kOverwrite);
            if (histoInvYieldWeightedAverageEtaStat)  histoInvYieldWeightedAverageEtaStat->Write("CorrectedYieldEta",TObject::kOverwrite);
            if (graphCorrectedYieldWeightedAverageEtaSys)   graphCorrectedYieldWeightedAverageEtaSys->Write("EtaSystError",TObject::kOverwrite);

            if (xSection != 0){
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
            }
            if (graphAcceptanceEtaWeighted)                 graphAcceptanceEtaWeighted->Write("AcceptanceEta", TObject::kOverwrite);
            if (graphEfficiencyEtaWeighted)                 graphEfficiencyEtaWeighted->Write("EfficiencyEta", TObject::kOverwrite);
            if (graphEffTimesAccEtaWeighted)                graphEffTimesAccEtaWeighted->Write("EffTimesAccEta", TObject::kOverwrite);
            if (graphMassEtaDataWeighted)                   graphMassEtaDataWeighted->Write("Eta_Mass_data", TObject::kOverwrite);
            if (graphMassEtaMCWeighted)                     graphMassEtaMCWeighted->Write("Eta_Mass_MC", TObject::kOverwrite);
            if (graphWidthEtaDataWeighted)                  graphWidthEtaDataWeighted->Write("Eta_Width_data", TObject::kOverwrite);
            if (graphWidthEtaMCWeighted)                    graphWidthEtaMCWeighted->Write("Eta_Width_MC", TObject::kOverwrite);
            if (sysAvailSingleEta[0]){
                for (Int_t k = 0; k< nRelSysErrEtaSources ; k++ ){
                    if (graphRelSysErrEtaSourceWeighted[k]) graphRelSysErrEtaSourceWeighted[k]->Write(graphRelSysErrEtaSourceWeighted[k]->GetName(), TObject::kOverwrite);
                }
            }

            for (Int_t i=0; i< nrOfTrigToBeComb; i++){
                if (histoAcceptanceEta[i])                  histoAcceptanceEta[i]->Write(Form("AcceptanceEta_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoEfficiencyEta[i])                  histoEfficiencyEta[i]->Write(Form("EfficiencyEta_%s",triggerName[i].Data()),TObject::kOverwrite);
                if (histoRawYieldEta[i])                    histoRawYieldEta[i]->Write(Form("RAWYieldPerEventsEta_%s",triggerName[i].Data()),TObject::kOverwrite);
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
            if (doEtaToOmega){
                if (histoEtaToOmegaWeightedAverageStat)       histoEtaToOmegaWeightedAverageStat->Write(Form("EtaToOmega%sStatError",addNameBinshift.Data()),TObject::kOverwrite);
                if (histoEtaToOmegaExtendedUsingFit)          histoEtaToOmegaExtendedUsingFit->Write(Form("EtaToOmega_extendedWithFit%s",addNameBinshift.Data()),TObject::kOverwrite);
                if (graphEtaToOmegaWeightedAverageStat)       graphEtaToOmegaWeightedAverageStat->Write(Form("graphEtaToOmega%sStatError",addNameBinshift.Data()),TObject::kOverwrite);
                if (graphEtaToOmegaWeightedAverageSys )       graphEtaToOmegaWeightedAverageSys->Write(Form("EtaToOmega%sSystError",addNameBinshift.Data()),TObject::kOverwrite);
                for (Int_t i=0; i< nrOfTrigToBeComb; i++){
                    if (graphsEtaToOmegaRemoved0[i])          graphsEtaToOmegaRemoved0[i]->Write(Form("EtaToOmega%sStatError_%s",addNameBinshift.Data(), triggerName[i].Data()),TObject::kOverwrite);
                    if (graphsEtaToOmegaSysRemoved0[i])       graphsEtaToOmegaSysRemoved0[i]->Write(Form("EtaToOmega%sSystError_%s",addNameBinshift.Data(), triggerName[i].Data()),TObject::kOverwrite);
                }
            }
        }

    fileOutputForComparisonFullyCorrected->Write();
    fileOutputForComparisonFullyCorrected->Close();
}
