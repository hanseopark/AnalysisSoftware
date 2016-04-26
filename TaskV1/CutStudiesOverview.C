//***********************************************************************************************
//**************************** CutStudiesOverview ***********************************************
//***********************************************************************************************
/************************************************************************************************
 ******     provided by Gamma Conversion Group, PWG4,                                         *****
 ******        Ana Marin, marin@physi.uni-heidelberg.de                                        *****
 ******        Friederike Bock, friederike.bock@cern.ch                                        *****
 ***********************************************************************************************/

#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
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
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

struct SysErrorConversion {
   Double_t value;
   Double_t error;
   //    TString name;
};


void CutStudiesOverview(TString CombineCutsName = "CombineCuts.dat", 
                        TString suffix = "gif", 
                        TString meson = "", 
                        TString isMC = "", 
                        TString optionMult = "", 
                        TString optionEnergy = "", 
                        TString cutVariationName = "", 
                        Int_t NumberOfCuts = 1, 
                        Bool_t optGammaOn = 1, 
                        TString fDalitz = "", 
                        TString optionPeriod = "No", 
                        Int_t mode = 9,
                        Bool_t doBarlow = kFALSE ){

    // Define global arrays
    TString     cutNumber               [50];
    TString     cutNumberAdv            [50];
    Double_t    nColls                  [50];
    TString     prefix2                                         = "";
    Bool_t      doMassRatio                                     = kTRUE;
    Bool_t      correctionFilesAvail                            = kTRUE;
    
    // Set common default plot style
    StyleSettingsThesis();
    SetPlotStyle();

    // Set cutvariation-name to "" for no explicit name
    if (cutVariationName.CompareTo("None")==0){
        cutVariationName                                        = "";
    }
        
    // Define Output Directory    
    TString outputDir                                           = Form("CutStudies/%s",optionEnergy.Data());
    if (cutVariationName.CompareTo("None")!=0) outputDir        = Form("CutStudies/%s/%s",optionEnergy.Data(),cutVariationName.Data());
    TString outputDirRootFile                                   = Form("CutStudies/%s",optionEnergy.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    
    // Define meson names for plots
    TString textMeson;
    Bool_t isEta                                                = kFALSE;
    if (meson.CompareTo("Pi0")==0 || meson.CompareTo("Pi0EtaBinning")==0){
        textMeson                                               = "#pi^{0}";
    } else {
        textMeson                                               = "#eta";
        isEta                                                   = kTRUE;
    }
    // Define input and output MC/data
    if (isMC.CompareTo("kTRUE") ==0){
        prefix2                                                 = "MC";
    } else {
        prefix2                                                 = "data";
    }
    
    // Set collisions system
    TString collisionSystem     = ReturnFullCollisionsSystem(optionEnergy);   
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;     
    }
    TString detectionProcess    = ReturnFullTextReconstructionProcess(mode);
    TString process             = Form("%s #rightarrow #gamma#gamma", textMeson.Data());
    
    // Define colors for differnt cuts
    Color_t color[20]                                           = { kBlack, kAzure, kGreen+2,kOrange+2,kRed, kViolet,  kBlue-9, kMagenta+4,
                                                                    kCyan+3, kCyan-10, kCyan, kGreen+4, kGreen-9, 
                                                                    kGreen,  kYellow+4, kYellow+3, kSpring+10,
                                                                    kMagenta-8, kGray, kGray+3};

    // Read cuts from CutSelection file                    
    ifstream in(CombineCutsName.Data());
    cout<<"Available Cuts:"<<endl;
    string TempCutNumber;
    Int_t Number = 0;
    while(getline(in, TempCutNumber)){
        TString tempCutNumberAdv                                = TempCutNumber;
        cutNumberAdv[Number]                                    = tempCutNumberAdv;
        cutNumber[Number]                                       = tempCutNumberAdv;
        cout<< cutNumber[Number]<<endl;
        Number++;
    }
    cout<<"=========================="<<endl;

    cout << "analysing " << cutVariationName << " cut variations" << endl;

    cout << " " << endl;

    // Determine number of collisions for PbPb case
    if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0) {
        for (Int_t i=0; i< NumberOfCuts; i++){
            nColls[i]                                             = GetNCollFromCutNumber(cutNumber[i]);
        }
    } else {
        for (Int_t i=0; i< NumberOfCuts; i++) nColls[i]         = 1.;
    }

    // Define necessary histogram/file/string arrays 
    const Int_t ConstNumberOfCuts                               = NumberOfCuts;
    const Int_t MaxNumberOfCuts = 20;
    if(ConstNumberOfCuts > MaxNumberOfCuts){
        cout << "Too many cuts, beware!" << endl;
        return;
    }
    TString FileNameCorrected           [MaxNumberOfCuts];
    TString FileNamePi0EtaBinCorrected  [MaxNumberOfCuts];
    TString FileNameUnCorrected         [MaxNumberOfCuts];
    TFile*  Cutcorrfile                 [ConstNumberOfCuts];
    TFile*  CutcorrPi0EtaBinfile        [ConstNumberOfCuts];
    TFile*  Cutuncorrfile               [ConstNumberOfCuts];
    TH1D*   histoCorrectedYieldCut      [ConstNumberOfCuts];
    TH1D*   histoCorrectedYieldPi0EtaCut[ConstNumberOfCuts];
    TH1D*   histoRawClusterPtCut        [ConstNumberOfCuts];
    TH1D*   histoTrueEffiCut            [ConstNumberOfCuts];
    TH1D*   histoAcceptanceCut          [ConstNumberOfCuts];
    TH1D*   histoRawYieldCut            [ConstNumberOfCuts];
    TH1D*   histoMassRatioCut           [ConstNumberOfCuts];
    TH1D*   histoEtaToPi0Cut            [ConstNumberOfCuts];
    TH1D*   histoSBCut                  [ConstNumberOfCuts];
    TH1D*   histoSBNarrowCut            [ConstNumberOfCuts];
    TH1D*   histoSignCut                [ConstNumberOfCuts];
    TH1D*   histoSignNarrowCut          [ConstNumberOfCuts];
    TH1D*   histoRatioCorrectedYieldCut [ConstNumberOfCuts];
    TH1D*   histoRatioTrueEffiCut       [ConstNumberOfCuts];
    TH1D*   histoRatioAcceptanceCut     [ConstNumberOfCuts];
    TH1D*   histoRatioRawYieldCut       [ConstNumberOfCuts];
    TH1D*   histoRatioMassRatioCut      [ConstNumberOfCuts];
    TH1D*   histoRatioEtaToPi0Cut       [ConstNumberOfCuts];
    TH1D*   histoRatioSBCut             [ConstNumberOfCuts];
    TH1D*   histoRatioSBNarrowCut       [ConstNumberOfCuts];
    TH1D*   histoRatioSignCut           [ConstNumberOfCuts];
    TH1D*   histoRatioSignNarrowCut     [ConstNumberOfCuts];
    TH1D*   histoRatioRawClusterPtCut   [ConstNumberOfCuts];
    TString cutStringsName              [MaxNumberOfCuts];

    // Define yield extraction error graphs
    TGraphAsymmErrors* systErrGraphNegYieldExt                  = NULL;
    TGraphAsymmErrors* systErrGraphPosYieldExt                  = NULL;
    TGraphAsymmErrors* systErrGraphNegYieldExtPi0EtaBinning     = NULL;
    TGraphAsymmErrors* systErrGraphPosYieldExtPi0EtaBinning     = NULL;

    TGraphAsymmErrors* systErrGraphBGEstimate                   = NULL;
    TString            centralityString                         = "";
    Double_t           maxPt                                    = 0;
    Bool_t             kSpecialTrigger                          = kFALSE;
    
    if (cutVariationName.CompareTo("SpecialTrigg") == 0){
            kSpecialTrigger                                     = kTRUE;
    }        
    for (Int_t i=0; i< NumberOfCuts; i++){
        // Define CutSelections
        TString fEventCutSelection                              = "";
        TString fGammaCutSelection                              = "";
        TString fClusterCutSelection                            = "";
        TString fElectronCutSelection                           = "";
        TString fMesonCutSelection                              = "";
        // disentangle cut selection
        if (mode == 9){
            ReturnSeparatedCutNumber(cutNumberAdv[i].Data(), fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
            fEventCutSelection                                = fGammaCutSelection(0,7);
            fGammaCutSelection                                = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
            cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
        } else {
            ReturnSeparatedCutNumberAdvanced(cutNumberAdv[i].Data(),fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
        }
        // Check if there was a special trigger among the cuts
        TString fTrigger = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),1);
        if (fTrigger.Atoi()>3 && cutVariationName.CompareTo("") == 0){
            cutVariationName                                  ="SpecialTrigg";
            kSpecialTrigger                                   = kTRUE;
        }    
        // only read corrected file if "Special trigger was used"

        FileNameCorrected[i]                                = Form( "%s/%s/%s_%s_GammaConvV1%sCorrection_%s.root", cutNumberAdv[i].Data(), optionEnergy.Data(), meson.Data(), prefix2.Data(), 
                                                                    fDalitz.Data(), cutNumber[i].Data());
        cout<< FileNameCorrected[i] << endl;
        Cutcorrfile[i]                                      = new TFile(FileNameCorrected[i]);    
        if (Cutcorrfile[i]->IsZombie()){
            correctionFilesAvail                            = kFALSE;
        }
            
        if (isEta){
            FileNamePi0EtaBinCorrected[i]                   = Form("%s/%s/Pi0EtaBinning_%s_GammaConvV1%sCorrection_%s.root", cutNumberAdv[i].Data(), optionEnergy.Data(), prefix2.Data(), 
                                                                    fDalitz.Data(), cutNumber[i].Data());
            CutcorrPi0EtaBinfile[i]                         = new TFile(FileNamePi0EtaBinCorrected[i]);    
            if (CutcorrPi0EtaBinfile[i]->IsZombie()){
                cout << "Pi0/Eta ratio can't be compared file: " << FileNamePi0EtaBinCorrected[i] << " missing! "<< endl;
                isEta                                       = kFALSE;
            }
        }    

        // read uncorrected file
        FileNameUnCorrected[i]                                  = Form("%s/%s/%s_%s_GammaConvV1%sWithoutCorrection_%s.root",cutNumberAdv[i].Data(), optionEnergy.Data(), meson.Data(), prefix2.Data(), 
                                                                    fDalitz.Data(), cutNumber[i].Data());
        cout<< FileNameUnCorrected[i] << endl;
        Cutuncorrfile[i]                                        = new TFile(FileNameUnCorrected[i]);
        if (Cutuncorrfile[i]->IsZombie()) return;
        
        // put proper cutvariation labeling for plots
        if (cutVariationName.Contains("SpecialTrigg")){
            fTrigger                                            = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
            cutStringsName[i]                                   = AnalyseSpecialTriggerCut(fTrigger.Atoi(), optionPeriod);      
        } else if (cutVariationName.Contains("MultiplicityPP")){
            TString minMult                                     = fEventCutSelection(GetEventCentralityMinCutPosition(),1);
            TString maxMult                                     = fEventCutSelection(GetEventCentralityMaxCutPosition(),1);
            cutStringsName[i]                                   = AnalysePPMultiplicityCut(minMult.Atoi(),maxMult.Atoi());      
        } else if (cutVariationName.Contains("V0Reader")){
            TString fV0Reader                                   = fGammaCutSelection(GetPhotonV0FinderCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseV0ReaderCut(fV0Reader.Atoi());      
        } else if (cutVariationName.Contains("Eta")){
            TString fEtaCut                                     = fGammaCutSelection(GetPhotonEtaCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseEtaCut(fEtaCut.Atoi());
        } else if (cutVariationName.Contains("RCutAndPhotonQuality")){
            TString fRCut                                       = fGammaCutSelection(GetPhotonMinRCutPosition(fGammaCutSelection),1);
            TString fPhotonQuality                              = fGammaCutSelection(GetPhotonSharedElectronCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseRCutAndQuality(fRCut.Atoi(), fPhotonQuality.Atoi());
        } else if (cutVariationName.Contains("RCut")){
            TString fRCut                                       = fGammaCutSelection(GetPhotonMinRCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseRCut(fRCut.Atoi());
        } else if (cutVariationName.Contains("SinglePt")){
            TString fSinglePtCut                                = fGammaCutSelection(GetPhotonSinglePtCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseSinglePtCut(fSinglePtCut.Atoi());
        } else if (cutVariationName.Contains("TPCCluster")){
            TString fClusterCut                                 = fGammaCutSelection(GetPhotonClsTPCCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseTPCClusterCut(fClusterCut.Atoi());            
        } else if (cutVariationName.Contains("dEdxE")){     
            TString fdEdxCut                                    = fGammaCutSelection(GetPhotonEDedxSigmaCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseTPCdEdxCutElectronLine(fdEdxCut.Atoi());
        } else if (cutVariationName.Contains("dEdxPi")){    
            TString fdEdxCut                                    = fGammaCutSelection(GetPhotonPiDedxSigmaCutPosition(fGammaCutSelection),3);
            cutStringsName[i]                                   = AnalyseTPCdEdxCutPionLine(fdEdxCut.Data());
        } else if (cutVariationName.Contains("TOF")){
            TString fTOFelectronPIDCut                          = fGammaCutSelection(GetPhotonTOFelectronPIDCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseTOFelectronPIDCut(fTOFelectronPIDCut.Atoi());
        } else if (cutVariationName.Contains("Qt")){
            TString fQtCut                                      = fGammaCutSelection(GetPhotonQtMaxCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseQtMaxCut(fQtCut.Atoi());
        } else if (cutVariationName.Contains("Chi2")){
            TString fChi2Cut                                    = fGammaCutSelection(GetPhotonChi2GammaCutPosition(fGammaCutSelection),1);
            TString fPsiPairCut                                 = fGammaCutSelection(GetPhotonPsiPairCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseChi2GammaCut(fChi2Cut.Atoi(), fPsiPairCut.Atoi());
        } else if (cutVariationName.Contains("PsiPairAndR")){
            TString fPsiPairCut                                 = fGammaCutSelection(GetPhotonPsiPairCutPosition(fGammaCutSelection),1);
            TString fRCut                                       = fGammaCutSelection(GetPhotonMinRCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalysePsiPairAndR(fPsiPairCut.Atoi(), fRCut.Atoi());      
        } else if (cutVariationName.Contains("PsiPair")){
            TString fPsiPairCut                                 = fGammaCutSelection(GetPhotonPsiPairCutPosition(fGammaCutSelection),1);
            TString fChi2Cut                                    = fGammaCutSelection(GetPhotonChi2GammaCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalysePsiPair(fPsiPairCut.Atoi(), fChi2Cut.Atoi());   
        } else if (cutVariationName.Contains("DCAZPhoton")){   
            TString fDCAZCut                                    = fGammaCutSelection(GetPhotonDcaZPrimVtxCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseDCAZPhotonCut(fDCAZCut.Atoi());
        } else if (cutVariationName.Contains("CosPoint")){
            TString fCosPoint                                   = fGammaCutSelection(GetPhotonCosinePointingAngleCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalyseCosPointCut(fCosPoint.Atoi());                 
        } else if (cutVariationName.Contains("PhotonQuality")){
            TString fPhotonQuality                              = fGammaCutSelection(GetPhotonSharedElectronCutPosition(fGammaCutSelection),1);
            cutStringsName[i]                                   = AnalysePhotonQuality(fPhotonQuality.Atoi());                  
        } else if (cutVariationName.Contains("ConvPhi")){    
            cutStringsName[i]                                   = AnalyseConvPhiExclusionCut(fGammaCutSelection);            
        } else if (cutVariationName.Contains("BG")){
            TString fBGCut                                      = fMesonCutSelection(GetMesonBGSchemeCutPosition(),3);
            cutStringsName[i]                                   = AnalyseBackgroundScheme(fBGCut.Data());   
        } else if (cutVariationName.Contains("Rapidity")){
            TString fRapidityCut                                = fMesonCutSelection(GetMesonRapidityCutPosition(),1);
            cutStringsName[i]                                   = AnalyseRapidityMesonCut(fRapidityCut.Atoi());      
        } else if (cutVariationName.Contains("Alpha")){
            TString fAlphaCut                                   = fMesonCutSelection(GetMesonAlphaCutPosition(),1);
            cutStringsName[i]                                   = AnalyseAlphaMesonCut(fAlphaCut.Atoi());
        } else if (cutVariationName.Contains("OpeningAngle")){
          TString fMesonOpeningAngleCut                         = fMesonCutSelection(GetMesonOpeningAngleCutPosition(),1);
          cutStringsName[i]                                     = AnalyseMesonOpeningAngleCut(fMesonOpeningAngleCut.Atoi());
        } else if (cutVariationName.Contains("Cent")){
            cutStringsName[i]                                   = GetCentralityString(fEventCutSelection.Data());
        } else if (cutVariationName.Contains("DiffRapWindow")){
            TString fRapidityCut                                = fMesonCutSelection(GetMesonRapidityCutPosition(),1);
            cutStringsName[i]                                   = AnalyseRapidityMesonCutpPb(fRapidityCut.Atoi());      
        } else if (cutVariationName.Contains("MCSmearing")){
            TString fMCSmearing                                 = fMesonCutSelection(GetMesonUseMCPSmearingCutPosition(),1);
            cutStringsName[i]                                   = AnalyseMCSmearingCut(fMCSmearing.Atoi());      
        } else if (cutVariationName.Contains("ClusterTrackMatchingCalo")){    
            TString fTrackMatching                              = fClusterCutSelection(GetClusterTrackMatchingCutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseTrackMatchingCaloCut(fTrackMatching.Atoi());
        } else if (cutVariationName.Contains("ClusterTrackMatching")){    
            TString fTrackMatching                              = fClusterCutSelection(GetClusterTrackMatchingCutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseTrackMatchingCut(fTrackMatching.Atoi());
        } else if (cutVariationName.Contains("ClusterMaterialTRD")){    
            TString fMinPhi                                     = fClusterCutSelection(GetClusterPhiMinCutPosition(fClusterCutSelection),1);
            TString fMaxPhi                                     = fClusterCutSelection(GetClusterPhiMaxCutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseAcceptanceCutPhiCluster(fMinPhi.Atoi(), fMaxPhi.Atoi());
        } else if (cutVariationName.Contains("ClusterM02")){    
            TString fMinM02Cut                                  = fClusterCutSelection(GetClusterMinM02CutPosition(fClusterCutSelection),1);
            TString fMaxM02Cut                                  = fClusterCutSelection(GetClusterMaxM02CutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseM02Cut(fMinM02Cut.Atoi(), fMaxM02Cut.Atoi());
        } else if (cutVariationName.Contains("ClusterNCells")){    
            TString fNCellsCut                                  = fClusterCutSelection(GetClusterMinNCellsCutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseNCellsCut(fNCellsCut.Atoi());
        } else if (cutVariationName.Contains("ClusterMinEnergy")){    
            TString fMinEnergyCut                               = fClusterCutSelection(GetClusterMinEnergyCutPosition(fClusterCutSelection),1);
            cout << fMinEnergyCut << "\t" << GetClusterMinEnergyCutPosition(fClusterCutSelection) << "\t"<< fClusterCutSelection.Length()<<endl;
            cutStringsName[i]                                   = AnalyseMinEnergyCut(fMinEnergyCut.Atoi());
        } else if (cutVariationName.Contains("ClusterTiming")){    
            TString fTimingCut                                  = fClusterCutSelection(GetClusterTimingCutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseClusterTimingCut(fTimingCut.Atoi());
        } else if (cutVariationName.Contains("ClusterNonLinearity")){
          TString fClusterNonLinearity                          = fClusterCutSelection(GetClusterNonLinearityCutPosition(fClusterCutSelection),2);
          cutStringsName[i]                                     = AnalyseClusterNonLinearityCut(fClusterNonLinearity.Atoi());
        } else {
            cutStringsName[i]                                   = cutNumberAdv[i].Data();
        }
        
        // read histograms from corrected file
        if (correctionFilesAvail){
            // for first cut read yield extraction errors as well
            if(i == 0){
                centralityString                                = GetCentralityString(fEventCutSelection.Data());
                cout << centralityString.Data()<< endl;
                systErrGraphNegYieldExt                         = (TGraphAsymmErrors*)Cutcorrfile[i]->Get(Form("%s_SystErrorRelNeg_YieldExtraction_%s",meson.Data(), centralityString.Data()));
                Double_t* negErrorYield                         = systErrGraphNegYieldExt->GetY();
                for (Int_t j = 0; j < systErrGraphNegYieldExt->GetN(); j++){
                    negErrorYield[j]                            = -1*negErrorYield[j];
                }
                systErrGraphPosYieldExt                         = (TGraphAsymmErrors*)Cutcorrfile[i]->Get(Form("%s_SystErrorRelPos_YieldExtraction_%s",meson.Data(), centralityString.Data()));
                systErrGraphBGEstimate                          = (TGraphAsymmErrors*)Cutcorrfile[i]->Get(Form("%s_SystErrorRel_BGEstimate_%s",meson.Data(), centralityString.Data()));
                
                if (isEta){
                    systErrGraphNegYieldExtPi0EtaBinning        = (TGraphAsymmErrors*)CutcorrPi0EtaBinfile[i]->Get(Form("Pi0EtaBinning_SystErrorRelNeg_YieldExtraction_%s", centralityString.Data()));
                    Double_t* negErrorYield                     = systErrGraphNegYieldExtPi0EtaBinning->GetY();
                    for (Int_t j = 0; j < systErrGraphNegYieldExtPi0EtaBinning->GetN(); j++){
                        negErrorYield[j]                        = -1*negErrorYield[j];
                    }
                    systErrGraphPosYieldExtPi0EtaBinning        = (TGraphAsymmErrors*)CutcorrPi0EtaBinfile[i]->Get(Form("Pi0EtaBinning_SystErrorRelPos_YieldExtraction_%s", centralityString.Data()));
                }    
            }
            
            TString nameCorrectedYield                          = "CorrectedYieldTrueEff";
            TString nameEfficiency                              = "TrueMesonEffiPt";
            TString nameAcceptance                              = "fMCMesonAccepPt";
            if ( mode == 4 || mode == 5 ){
                nameCorrectedYield                              = "CorrectedYieldNormEff";
                nameEfficiency                                  = "MesonEffiPt";
            }    
//             if ((mode == 2 || mode == 3) && !meson.Contains("Pi0")){
//                 nameCorrectedYield                                 = "CorrectedYieldNormEff";
//                 nameEfficiency                                     = "MesonEffiPt";                
//             }
            histoCorrectedYieldCut[i]                           = (TH1D*)Cutcorrfile[i]->Get(nameCorrectedYield.Data());
            histoCorrectedYieldCut[i]->SetName(Form("%s_%s", nameCorrectedYield.Data(),cutNumber[i].Data()));
            histoTrueEffiCut[i]                                 = (TH1D*)Cutcorrfile[i]->Get(nameEfficiency.Data());
            histoTrueEffiCut[i]->SetName(Form("%s_%s", nameEfficiency.Data(), cutNumber[i].Data()));
            histoAcceptanceCut[i]                               = (TH1D*)Cutcorrfile[i]->Get(nameAcceptance.Data());
            histoAcceptanceCut[i]->SetName(Form("%s_%s", nameAcceptance.Data(), cutNumber[i].Data()));
            histoMassRatioCut[i]                                = (TH1D*)Cutcorrfile[i]->Get("histoRatioRecMass");
            if (histoMassRatioCut[i] == NULL) 
                histoMassRatioCut[i]                            = (TH1D*)Cutcorrfile[i]->Get("histoRatioRecMassGauss");
            if (histoMassRatioCut[i] == NULL )
                doMassRatio                                     = kFALSE;
            if (doMassRatio) histoMassRatioCut[i]->SetName(Form("histoMassRatio_%s", cutNumber[i].Data()));
            
            if (isEta){
                histoCorrectedYieldPi0EtaCut[i]                 = (TH1D*)CutcorrPi0EtaBinfile[i]->Get(nameCorrectedYield.Data());
                histoCorrectedYieldPi0EtaCut[i]->SetName(Form("Pi0EtaBinning_%s_%s", nameCorrectedYield.Data(),cutNumber[i].Data()));
                histoEtaToPi0Cut[i]                             = (TH1D*)histoCorrectedYieldCut[i]->Clone(Form("EtaToPi0_%s", cutNumber[i].Data()));
                histoEtaToPi0Cut[i]->Divide(histoEtaToPi0Cut[i],histoCorrectedYieldPi0EtaCut[i],1.,1.,"");
            }
        }
        // read histograms from uncorrected file
        histoRawYieldCut[i]                                     = (TH1D*)Cutuncorrfile[i]->Get("histoYieldMesonPerEvent");
        histoRawYieldCut[i]->SetName(Form("histoYieldMesonPerEvent_%s",cutNumber[i].Data()));
        histoRawYieldCut[i]->Scale(1./nColls[i]);
        cout << "line " << __LINE__ << endl;
        histoSBCut[i]                                           = (TH1D*)Cutuncorrfile[i]->Get("histoSBdefaultMeson");
        histoSBCut[i]->SetName(Form("histoSBdefaultMeson_%s",cutNumber[i].Data()));
        cout << "line " << __LINE__ << endl;
        histoSBNarrowCut[i]                                     = (TH1D*)Cutuncorrfile[i]->Get("histoSBdefaultNarrowMeson");
        histoSBNarrowCut[i]->SetName(Form("histoSBdefaultNarrowMeson_%s",cutNumber[i].Data()));
        cout << "line " << __LINE__ << endl;
        histoSignCut[i]                                         = (TH1D*)Cutuncorrfile[i]->Get("histoSigndefaultMeson");
        histoSignCut[i]->SetName(Form("histoSigndefaultMeson_%s",cutNumber[i].Data()));
        cout << "line " << __LINE__ << endl;
        histoSignNarrowCut[i]                                   = (TH1D*)Cutuncorrfile[i]->Get("histoSigndefaultNarrowMeson");
        histoSignNarrowCut[i]->SetName(Form("histoSigndefaultNarrowMeson_%s",cutNumber[i].Data()));
        cout << "line " << __LINE__ << endl;
        
        if (mode == 2 || mode == 3 || mode == 4 || mode == 5){
            histoRawClusterPtCut[i]                             = (TH1D*)Cutuncorrfile[i]->Get("ClusterPtPerEvent");
        }    
        cout << "line " << __LINE__ << endl;
        
        // calculate ratios for meson measurements
        if (correctionFilesAvail){
            histoRatioCorrectedYieldCut[i]                      = (TH1D*) histoCorrectedYieldCut[i]->Clone(Form("histoRatioCorrectedYieldCut_%s", cutNumber[i].Data()));
            histoRatioCorrectedYieldCut[i]->Divide(histoRatioCorrectedYieldCut[i],histoCorrectedYieldCut[0],1.,1.,"B");
            if (i > 0){
                maxPt= histoCorrectedYieldCut[i]->GetBinCenter(histoCorrectedYieldCut[i]->GetNbinsX()) + 0.5* histoCorrectedYieldCut[i]->GetBinWidth(histoCorrectedYieldCut[i]->GetNbinsX());
            }
            cout << "line " << __LINE__ << endl;
            histoRatioTrueEffiCut[i]                            = (TH1D*) histoTrueEffiCut[i]->Clone(Form("histoRatioTrueEffiCut_%s", cutNumber[i].Data()));
            histoRatioTrueEffiCut[i]->Divide(histoRatioTrueEffiCut[i],histoTrueEffiCut[0],1.,1.,"B");
            cout << "line " << __LINE__ << endl;
            histoRatioAcceptanceCut[i]                          = (TH1D*) histoAcceptanceCut[i]->Clone(Form("histoRatioAcceptanceCut_%s", cutNumber[i].Data()));
            histoRatioAcceptanceCut[i]->Divide(histoRatioAcceptanceCut[i],histoAcceptanceCut[0],1.,1.,"B");
            cout << "line " << __LINE__ << endl;
            if (doMassRatio){
                histoRatioMassRatioCut[i]                       = (TH1D*) histoMassRatioCut[i]->Clone(Form("histoRatioMassRatio_%s", cutNumber[i].Data()));
                histoRatioMassRatioCut[i]->Divide(histoRatioMassRatioCut[i],histoMassRatioCut[0],1.,1.,"B");
            }    
            if (isEta){
                histoRatioEtaToPi0Cut[i]                        = (TH1D*) histoEtaToPi0Cut[i]->Clone(Form("histoEtaToPi0Ratio_%s", cutNumber[i].Data()));
                histoRatioEtaToPi0Cut[i]->Divide(histoRatioEtaToPi0Cut[i],histoEtaToPi0Cut[0],1.,1.,"B");
            }    
        }    
        histoRatioRawYieldCut[i]                                = (TH1D*) histoRawYieldCut[i]->Clone(Form("histoRatioRawYieldCut_%s", cutNumber[i].Data()));
        if (!kSpecialTrigger){
            histoRatioRawYieldCut[i]->Divide(histoRatioRawYieldCut[i],histoRawYieldCut[0],1.,1.,"B");
        } else {
            histoRatioRawYieldCut[i]->Divide(histoRatioRawYieldCut[i],histoRawYieldCut[0],1.,1.,"");
        }    
        cout << "line " << __LINE__ << endl;
        histoRatioSBCut[i]                                      = (TH1D*) histoSBCut[i]->Clone(Form("histoRatioSBCut_%s", cutNumber[i].Data()));
        histoRatioSBCut[i]->Divide(histoRatioSBCut[i],histoSBCut[0],1.,1.,"B");
        histoRatioSBNarrowCut[i]                                = (TH1D*) histoSBNarrowCut[i]->Clone(Form("histoRatioSBNarrowCut_%s", cutNumber[i].Data()));
        histoRatioSBNarrowCut[i]->Divide(histoRatioSBNarrowCut[i],histoSBNarrowCut[0],1.,1.,"B");
        histoRatioSignCut[i]                                    = (TH1D*) histoSignCut[i]->Clone(Form("histoRatioSignCut_%s", cutNumber[i].Data()));
        histoRatioSignCut[i]->Divide(histoRatioSignCut[i],histoSignCut[0],1.,1.,"B");
        histoRatioSignNarrowCut[i]                              = (TH1D*) histoSignNarrowCut[i]->Clone(Form("histoRatioSignNarrowCut_%s", cutNumber[i].Data()));
        histoRatioSignNarrowCut[i]->Divide(histoRatioSignNarrowCut[i],histoSignNarrowCut[0],1.,1.,"B");
        cout << "line " << __LINE__ << endl;
        if (mode == 2 || mode == 3 || mode == 4 || mode == 5){
            histoRatioRawClusterPtCut[i]                        = (TH1D*) histoRawClusterPtCut[i]->Clone(Form("histoRatioRawClusterPtCut%s", cutNumber[i].Data()));
            if (!kSpecialTrigger){
                histoRatioRawClusterPtCut[i]->Divide(histoRatioRawClusterPtCut[i],histoRawClusterPtCut[0],1.,1.,"B");
            } else {
                histoRatioRawClusterPtCut[i]->Divide(histoRatioRawClusterPtCut[i],histoRawClusterPtCut[0],1.,1.,"");
            }    
        }    
        cout << "line " << __LINE__ << endl;
        
    }

    cout<<"=========================="<<endl;


    //**************************************************************************************
    //********************* Plotting RAW-Yield *********************************************
    //**************************************************************************************

        TCanvas* canvasRawYieldMeson = new TCanvas("canvasRawYieldMeson","",1350,1500);  
        DrawGammaCanvasSettings( canvasRawYieldMeson,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padRawYield = new TPad("padRawYield", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padRawYield, 0.12, 0.02, 0.02, 0.);
        padRawYield->SetLogy();
        padRawYield->Draw();
        // lower pad definition
        TPad* padRawYieldRatios = new TPad("padRawYieldRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padRawYieldRatios, 0.12, 0.02, 0.0, 0.2);
        padRawYieldRatios->Draw();

        padRawYield->cd();
        padRawYield->SetTickx();
        padRawYield->SetTicky();

        // Plot raw yield in uppper panel
        padRawYield->cd();
        TLegend* legendRawMeson = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
        if (cutVariationName.Contains("dEdxPi")){
            legendRawMeson->SetTextSize(0.02);
        }
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
                Double_t scaleFactorRaw = 5.;
                if (kSpecialTrigger) scaleFactorRaw = 100.;
                
                DrawGammaSetMarker(histoRawYieldCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoRawYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s RAW Yield/(#it{N}_{ev} #it{N}_{coll})",textMeson.Data()),
                                        kTRUE, scaleFactorRaw, 10e-10, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                legendRawMeson->AddEntry(histoRawYieldCut[i],Form("standard: %s",cutStringsName[i].Data()));
            }
            else {
                if(i<20){
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRawYieldCut[i]->DrawCopy("same,e1,p");
                legendRawMeson->AddEntry(histoRawYieldCut[i],cutStringsName[i].Data());
            }
            
        }
        legendRawMeson->Draw();
        // Labeling of plot
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);
                                                
        padRawYieldRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.5;
                Double_t maxYRatio = 1.55; //qui
//                 if (cutVariationName.Contains("PhotonQuality")){
//                     minYRatio = 0.001;
//                     maxYRatio = 2;
//                     padRawYieldRatios->SetLogy(1);
//                 }      
                SetStyleHistoTH1ForGraphs(histoRatioRawYieldCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioRawYieldCut[i], 20, 1.,color[0],color[0]);
                histoRatioRawYieldCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                
                for(Int_t b = 0; b< histoRatioRawYieldCut[i]->GetNbinsX(); b++){
                    histoRatioRawYieldCut[i]->SetBinError(b+1,histoRawYieldCut[i]->GetBinError(b+1)/histoRawYieldCut[i]->GetBinContent(b+1));
                }
                histoRatioRawYieldCut[i]->SetFillColor(kGray+2);
                histoRatioRawYieldCut[i]->SetFillStyle(0);
                histoRatioRawYieldCut[i]->DrawCopy("p,e2");  

//                 histoRatioRawYieldCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioRawYieldCut[i]->DrawCopy("same,e1,p");
            }
            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasRawYieldMeson->Update();
        canvasRawYieldMeson->SaveAs(Form("%s/%s_%s_RAWYield.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasRawYieldMeson;

    if (cutVariationName.Contains("SpecialTrigg")){
        //**************************************************************************************
        //********************* Plotting RAW-Yield for special triggers  ***********************
        //**************************************************************************************        

        TCanvas* canvasRawYieldsTrigg = new TCanvas("canvasRawYieldsTrigg","",1000,1000);  // gives the page size
        DrawGammaCanvasSettings( canvasRawYieldsTrigg,  0.1, 0.02, 0.04, 0.08);
        canvasRawYieldsTrigg->SetLogy(1);

        // Set legend
        
        TLegend* legendSpecRawTrigger = GetAndSetLegend2(0.15,0.15,0.3,0.15+1.15*0.032*NumberOfCuts, 1000*0.032); 
        // Draw Raw yield for different triggers
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
                Double_t scaleFactorRaw = 10.;
                DrawGammaSetMarker(histoRawYieldCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoRawYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s RAW Yield/(#it{N}_{ev})",textMeson.Data()),
                                        kTRUE, scaleFactorRaw, 10e-10, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                legendSpecRawTrigger->AddEntry(histoRawYieldCut[i],Form("standard: %s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRawYieldCut[i]->DrawCopy("same,e1,p");
                legendSpecRawTrigger->AddEntry(histoRawYieldCut[i],cutStringsName[i].Data());
            }
        }
        legendSpecRawTrigger->Draw();
        // labeling of the plot
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);
        
        canvasRawYieldsTrigg->Update();
        canvasRawYieldsTrigg->SaveAs(Form("%s/%s_%s_TriggerYieldSpectra.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasRawYieldsTrigg;

        //**************************************************************************************
        //***************** Plotting RAW-Yield ratios for special triggers  ********************
        //**************************************************************************************                
        TCanvas* canvasRatioRawYields = new TCanvas("canvasRatioRawYields","",1000,1000);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioRawYields,  0.1, 0.02, 0.04, 0.08);
        if (mode==2 || mode == 4 )canvasRatioRawYields->SetLogy(1); 
        // create legend
        TLegend* legendRatioRaw = GetAndSetLegend2(0.55,0.15,0.7,0.15+1.15*0.032*NumberOfCuts, 1000*0.032);
        // find min und max
        Double_t maxRatio = 2e3;
        Double_t minRatio = 0;
        if (mode==2 || mode == 4){
            maxRatio = 2e5;
            minRatio = 0.1;
        }
        // plot ratios in canvas
        for(Int_t i = 1; i< NumberOfCuts; i++){
            if(i == 1){
                DrawGammaSetMarker(histoRatioRawYieldCut[i], 20, 1., color[1], color[1]);
                DrawAutoGammaMesonHistos( histoRatioRawYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("triggered Yield/standard"),
                                        kFALSE, 5., 10e-10, kTRUE,
                                        kTRUE, minRatio, maxRatio,
                                        kFALSE, 0., 10.);
                legendRatioRaw->AddEntry(histoRatioRawYieldCut[i],Form("%s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioRawYieldCut[i]->DrawCopy("same,e1,p");
                legendRatioRaw->AddEntry(histoRatioRawYieldCut[i],cutStringsName[i].Data());
            }
        }
        // labeling
        legendRatioRaw->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);
        
        canvasRatioRawYields->Update();
        canvasRatioRawYields->SaveAs(Form("%s/%s_%s_TriggerYieldRatio.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasRatioRawYields;
    }
    
    

    //**************************************************************************************
    //************************ Plotting SB  ************************************************
    //**************************************************************************************

        TCanvas* canvasSBMeson = new TCanvas("canvasSBMeson","",1350,1500);  
        DrawGammaCanvasSettings( canvasSBMeson,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padSB = new TPad("padSB", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padSB, 0.12, 0.02, 0.02, 0.);
        padSB->SetLogy();
        padSB->Draw();
        // lower pad definition
        TPad* padSBRatios = new TPad("padSBRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padSBRatios, 0.12, 0.02, 0.0, 0.2);
        padSBRatios->Draw();

        padSB->cd();
        padSB->SetTickx();
        padSB->SetTicky();

        // Plot SB in uppper panel
        padSB->cd();
        
        TLegend* legendSB = GetAndSetLegend2(0.27,0.02,0.35,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
        if (cutVariationName.Contains("dEdxPi")){
            legendSB->SetTextSize(0.02);
        }
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
            
                DrawGammaSetMarker(histoSBCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoSBCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s S/B",textMeson.Data()),
                                        kTRUE, 0., 1e-4, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                legendSB->AddEntry(histoSBCut[i],Form("standard: %s",cutStringsName[i].Data()));
            }
            else {
                if(i<20){
                    DrawGammaSetMarker(histoSBCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoSBCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoSBCut[i]->DrawCopy("same,e1,p");
                legendSB->AddEntry(histoSBCut[i],cutStringsName[i].Data());
            }
            
        }
        legendSB->Draw();
        // Labeling of plot
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);    
            
        padSBRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55; //qui
                if (cutVariationName.Contains("PhotonQuality")){
                    minYRatio = 0.001;
                    maxYRatio = 2;
                    padSBRatios->SetLogy(1);
                }      
                SetStyleHistoTH1ForGraphs(histoRatioSBCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioSBCut[i], 20, 1.,color[0],color[0]);
                histoRatioSBCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoRatioSBCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioSBCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioSBCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioSBCut[i]->DrawCopy("same,e1,p");
            }
            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasSBMeson->Update();
        canvasSBMeson->SaveAs(Form("%s/%s_%s_SB.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasSBMeson;


    //**************************************************************************************
    //************************ Plotting SBNarrow  ************************************************
    //**************************************************************************************

        TCanvas* canvasSBNarrowMeson = new TCanvas("canvasSBNarrowMeson","",1350,1500);  
        DrawGammaCanvasSettings( canvasSBNarrowMeson,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padSBNarrow = new TPad("padSBNarrow", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padSBNarrow, 0.12, 0.02, 0.02, 0.);
        padSBNarrow->SetLogy();
        padSBNarrow->Draw();
        // lower pad definition
        TPad* padSBNarrowRatios = new TPad("padSBNarrowRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padSBNarrowRatios, 0.12, 0.02, 0.0, 0.2);
        padSBNarrowRatios->Draw();

        padSBNarrow->cd();
        padSBNarrow->SetTickx();
        padSBNarrow->SetTicky();

        // Plot SBNarrow in uppper panel
        padSBNarrow->cd();
        TLegend* legendSBNarrow = GetAndSetLegend2(0.27,0.02,0.35,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
        if (cutVariationName.Contains("dEdxPi")){
            legendSBNarrow->SetTextSize(0.02);
        }
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
            
                DrawGammaSetMarker(histoSBNarrowCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoSBNarrowCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s S/B (narrow range)",textMeson.Data()),
                                        kTRUE, 0., 1e-4, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                legendSBNarrow->AddEntry(histoSBNarrowCut[i],Form("standard: %s",cutStringsName[i].Data()));
            }
            else {
                if(i<20){
                    DrawGammaSetMarker(histoSBNarrowCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoSBNarrowCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoSBNarrowCut[i]->DrawCopy("same,e1,p");
                legendSBNarrow->AddEntry(histoSBNarrowCut[i],cutStringsName[i].Data());
            }
            
        }
        legendSBNarrow->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);
            
        padSBNarrowRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55; //qui
                if (cutVariationName.Contains("PhotonQuality")){
                    minYRatio = 0.001;
                    maxYRatio = 2;
                    padSBNarrowRatios->SetLogy(1);
                }      
                SetStyleHistoTH1ForGraphs(histoRatioSBNarrowCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioSBNarrowCut[i], 20, 1.,color[0],color[0]);
                histoRatioSBNarrowCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoRatioSBNarrowCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioSBNarrowCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioSBNarrowCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioSBNarrowCut[i]->DrawCopy("same,e1,p");
            }
            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasSBNarrowMeson->Update();
        canvasSBNarrowMeson->SaveAs(Form("%s/%s_%s_SBNarrow.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasSBNarrowMeson;

        

    //**************************************************************************************
    //************************ Plotting SignNarrow  ************************************************
    //**************************************************************************************

        TCanvas* canvasSignNarrowMeson = new TCanvas("canvasSignNarrowMeson","",1350,1500);  
        DrawGammaCanvasSettings( canvasSignNarrowMeson,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padSignNarrow = new TPad("padSignNarrow", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padSignNarrow, 0.12, 0.02, 0.02, 0.);
        padSignNarrow->SetLogy();
        padSignNarrow->Draw();
        // lower pad definition
        TPad* padSignNarrowRatios = new TPad("padSignNarrowRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padSignNarrowRatios, 0.12, 0.02, 0.0, 0.2);
        padSignNarrowRatios->Draw();

        padSignNarrow->cd();
        padSignNarrow->SetTickx();
        padSignNarrow->SetTicky();

        // Plot SignNarrow in uppper panel
        padSignNarrow->cd();
        TLegend* legendSignNarrow = GetAndSetLegend2(0.27,0.02,0.35,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
        if (cutVariationName.Contains("dEdxPi")){
            legendSignNarrow->SetTextSize(0.02);
        }
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
            
                DrawGammaSetMarker(histoSignNarrowCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoSignNarrowCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s Significance (narrow range)",textMeson.Data()),
                                        kFALSE, 0., 1e-4, kFALSE,
                                        kTRUE, 1., 100.,
                                        kFALSE, 0., 14.);
                legendSignNarrow->AddEntry(histoSignNarrowCut[i],Form("standard: %s",cutStringsName[i].Data()));
            }
            else {
                if(i<20){
                    DrawGammaSetMarker(histoSignNarrowCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoSignNarrowCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoSignNarrowCut[i]->DrawCopy("same,e1,p");
                legendSignNarrow->AddEntry(histoSignNarrowCut[i],cutStringsName[i].Data());
            }
            
        }
        legendSignNarrow->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);    
            
        padSignNarrowRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55; //qui
                if (cutVariationName.Contains("PhotonQuality")){
                    minYRatio = 0.001;
                    maxYRatio = 2;
                    padSignNarrowRatios->SetLogy(1);
                }      
                SetStyleHistoTH1ForGraphs(histoRatioSignNarrowCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioSignNarrowCut[i], 20, 1.,color[0],color[0]);
                histoRatioSignNarrowCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoRatioSignNarrowCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioSignNarrowCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioSignNarrowCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioSignNarrowCut[i]->DrawCopy("same,e1,p");
            }
            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasSignNarrowMeson->Update();
        canvasSignNarrowMeson->SaveAs(Form("%s/%s_%s_SignNarrow.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasSignNarrowMeson;

        

    //**************************************************************************************
    //************************ Plotting Sign  ************************************************
    //**************************************************************************************

        TCanvas* canvasSignMeson = new TCanvas("canvasSignMeson","",1350,1500);  
        DrawGammaCanvasSettings( canvasSignMeson,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padSign = new TPad("padSign", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padSign, 0.12, 0.02, 0.02, 0.);
        padSign->SetLogy();
        padSign->Draw();
        // lower pad definition
        TPad* padSignRatios = new TPad("padSignRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padSignRatios, 0.12, 0.02, 0.0, 0.2);
        padSignRatios->Draw();

        padSign->cd();
        padSign->SetTickx();
        padSign->SetTicky();

        // Plot Sign in uppper panel
        padSign->cd();
        TLegend* legendSign = GetAndSetLegend2(0.27,0.02,0.35,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
        if (cutVariationName.Contains("dEdxPi")){
            legendSign->SetTextSize(0.02);
        }
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
            
                DrawGammaSetMarker(histoSignCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoSignCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s Significance",textMeson.Data()),
                                        kFALSE, 0., 1e-4, kFALSE,
                                        kTRUE, 1., 100.,
                                        kFALSE, 0., 14.);
                legendSign->AddEntry(histoSignCut[i],Form("standard: %s",cutStringsName[i].Data()));
            }
            else {
                if(i<20){
                    DrawGammaSetMarker(histoSignCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoSignCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoSignCut[i]->DrawCopy("same,e1,p");
                legendSign->AddEntry(histoSignCut[i],cutStringsName[i].Data());
            }
            
        }
        legendSign->Draw();
        // Labeling of plot
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);    
            
        padSignRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55; //qui
                if (cutVariationName.Contains("PhotonQuality")){
                    minYRatio = 0.001;
                    maxYRatio = 2;
                    padSignRatios->SetLogy(1);
                }      
                SetStyleHistoTH1ForGraphs(histoRatioSignCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioSignCut[i], 20, 1.,color[0],color[0]);
                histoRatioSignCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoRatioSignCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioSignCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioSignCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioSignCut[i]->DrawCopy("same,e1,p");
            }
            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasSignMeson->Update();
        canvasSignMeson->SaveAs(Form("%s/%s_%s_Sign.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasSignMeson;

        
    if (mode == 2 || mode == 3 || mode == 4 || mode == 5){
        //**************************************************************************************
        //********************* Plotting RAW-Cluster Yield *********************************************
        //**************************************************************************************

            TCanvas* canvasRawClusterPt = new TCanvas("canvasRawClusterPt","",1350,1500);  
            DrawGammaCanvasSettings( canvasRawClusterPt,  0.13, 0.02, 0.02, 0.09);
            // Upper pad definition
            TPad* padRawClusterPt = new TPad("padRawClusterPt", "", 0., 0.33, 1., 1.,-1, -1, -2);
            DrawGammaPadSettings( padRawClusterPt, 0.12, 0.02, 0.02, 0.);
            padRawClusterPt->SetLogy();
            padRawClusterPt->Draw();
            // lower pad definition
            TPad* padRawClusterPtRatios = new TPad("padRawClusterPtRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
            DrawGammaPadSettings( padRawClusterPtRatios, 0.12, 0.02, 0.0, 0.2);
            padRawClusterPtRatios->Draw();

            padRawClusterPt->cd();
            padRawClusterPt->SetTickx();
            padRawClusterPt->SetTicky();

            // Plot raw yield in uppper panel
            padRawClusterPt->cd();
            TLegend* legendRawCluster = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
            if (cutVariationName.Contains("dEdxPi")){
                legendRawCluster->SetTextSize(0.02);
            }
            for(Int_t i = 0; i< NumberOfCuts; i++){
                if(i == 0){
                    Double_t scaleFactorRaw = 5.;
                    if (kSpecialTrigger) scaleFactorRaw = 5.;
                    
                    DrawGammaSetMarker(histoRawClusterPtCut[i], 20, 1., color[0], color[0]);
                    DrawAutoGammaMesonHistos( histoRawClusterPtCut[i],
                                            "", "#it{p}_{T} (GeV/#it{c})", "cluster RAW Yield/(#it{N}_{ev} #it{N}_{coll})",
                                            kTRUE, scaleFactorRaw, 10e-10, kTRUE,
                                            kFALSE, 0.0, 0.030,
                                            kFALSE, 0., 10.);
                    legendRawCluster->AddEntry(histoRawClusterPtCut[i],Form("standard: %s",cutStringsName[i].Data()));
                }
                else {
                    if(i<20){
                        DrawGammaSetMarker(histoRawClusterPtCut[i], 20+i, 1.,color[i],color[i]);
                    } else {
                        DrawGammaSetMarker(histoRawClusterPtCut[i], 20+i, 1.,color[i-20],color[i-20]);
                    }
                    histoRawClusterPtCut[i]->DrawCopy("same,e1,p");
                    legendRawCluster->AddEntry(histoRawClusterPtCut[i],cutStringsName[i].Data());
                }
                
            }
            legendRawCluster->Draw();
            // Labeling of plot
            PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, "#gamma candidates", detectionProcess, 42, 0.03, optionPeriod); 
                
            padRawClusterPtRatios->cd();
            for(Int_t i = 0; i< NumberOfCuts; i++){
                if(i==0){
                    // Set ratio min and max
                    Double_t minYRatio = 0.45;
                    Double_t maxYRatio = 2.05; //qui
                    if (cutVariationName.Contains("PhotonQuality")){
                        minYRatio = 0.001;
                        maxYRatio = 2;
                        padRawClusterPtRatios->SetLogy(1);
                    }      
                    SetStyleHistoTH1ForGraphs(histoRatioRawClusterPtCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                    DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20, 1.,color[0],color[0]);
                    histoRatioRawClusterPtCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                    histoRatioRawClusterPtCut[i]->DrawCopy("p,e1");
                } else{
                    if(i<20){
                        DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20+i, 1.,color[i],color[i]);
                    } else {
                        DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20+i, 1.,color[i-20],color[i-20]);
                    }
                    histoRatioRawClusterPtCut[i]->DrawCopy("same,e1,p");
                }
                DrawGammaLines(0., maxPt,1., 1.,0.1);
            }

            canvasRawClusterPt->Update();
            canvasRawClusterPt->SaveAs(Form("%s/%s_%s_RAWYieldCluster.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
            delete canvasRawClusterPt;

        if (cutVariationName.Contains("SpecialTrigg") && meson.CompareTo("Pi0") == 0 ){
            //**************************************************************************************
            //********************* Plotting RAW-Yield for special triggers  ***********************
            //**************************************************************************************        

            TCanvas* canvasRawClusterPtTrigg = new TCanvas("canvasRawClusterPtTrigg","",1000,1000);  // gives the page size
            DrawGammaCanvasSettings( canvasRawClusterPtTrigg,  0.1, 0.02, 0.04, 0.08);
            canvasRawClusterPtTrigg->SetLogy(1);

            // Set legend
            TLegend* legendSpecRawClusterTrigger =GetAndSetLegend2(0.15,0.11,0.3,0.11+1.15*0.032*NumberOfCuts, 1000*0.032);
            // Draw Raw yield for different triggers
            for(Int_t i = 0; i< NumberOfCuts; i++){
                if(i == 0){
                    Double_t scaleFactorRaw = 10.;
                    DrawGammaSetMarker(histoRawClusterPtCut[i], 20, 1., color[0], color[0]);
                    DrawAutoGammaMesonHistos( histoRawClusterPtCut[i],
                                            "", "#it{p}_{T,cluster} (GeV/#it{c})", "cluster RAW Yield/(#it{N}_{ev})",
                                            kTRUE, scaleFactorRaw, 10e-10, kTRUE,
                                            kFALSE, 0.0, 0.030,
                                            kFALSE, 0., 10.);
                    legendSpecRawClusterTrigger->AddEntry(histoRawClusterPtCut[i],Form("standard: %s",cutStringsName[i].Data()));
                } else {
                    if(i<20){
                        DrawGammaSetMarker(histoRawClusterPtCut[i], 20+i, 1.,color[i],color[i]);
                    } else {
                        DrawGammaSetMarker(histoRawClusterPtCut[i], 20+i, 1.,color[i-20],color[i-20]);
                    }
                    histoRawClusterPtCut[i]->DrawCopy("same,e1,p");
                    legendSpecRawClusterTrigger->AddEntry(histoRawClusterPtCut[i],cutStringsName[i].Data());
                }
            }
            legendSpecRawClusterTrigger->Draw();
            // labeling of the plot
            PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, "#gamma candidates", detectionProcess, 42, 0.03, optionPeriod);
            
            canvasRawClusterPtTrigg->Update();
            canvasRawClusterPtTrigg->SaveAs(Form("%s/%s_TriggerYieldCluster.%s",outputDir.Data(),prefix2.Data(),suffix.Data()));
            delete canvasRawClusterPtTrigg;

            //**************************************************************************************
            //***************** Plotting RAW-Yield ratios for special triggers  ********************
            //**************************************************************************************                
            TCanvas* canvasRatioRawClusterPt = new TCanvas("canvasRatioRawClusterPt","",1000,1000);  // gives the page size
            DrawGammaCanvasSettings( canvasRatioRawClusterPt,  0.1, 0.02, 0.04, 0.08);
            if (mode==2 || mode == 4 )canvasRatioRawClusterPt->SetLogy(1); 
            // create legend
            TLegend* legendRatioRawClusterTrigger = GetAndSetLegend2(0.55,0.15,0.7,0.15+1.15*0.032*NumberOfCuts, 1000*0.032);
            // find min und max
            Double_t maxRatio = 2e3;
            Double_t minRatio = 0;
            if (mode==2 || mode == 4){
                maxRatio = 3e5;
                minRatio = 1;
            }
            // plot ratios in canvas
            for(Int_t i = 1; i< NumberOfCuts; i++){
                if(i == 1){
                    DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20, 1., color[1], color[1]);
                    DrawAutoGammaMesonHistos( histoRatioRawClusterPtCut[i],
                                            "", "#it{p}_{T,cluster} (GeV/#it{c})", Form("triggered Yield/standard"),
                                            kFALSE, 5., 10e-10, kTRUE,
                                            kTRUE, minRatio, maxRatio,
                                            kFALSE, 0., 10.);
                    legendRatioRawClusterTrigger->AddEntry(histoRatioRawClusterPtCut[i],Form("%s",cutStringsName[i].Data()));
                } else {
                    if(i<20){
                        DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20+i, 1.,color[i],color[i]);
                    } else {
                        DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20+i, 1.,color[i-20],color[i-20]);
                    }
                    histoRatioRawClusterPtCut[i]->DrawCopy("same,e1,p");
                    legendRatioRawClusterTrigger->AddEntry(histoRatioRawClusterPtCut[i],cutStringsName[i].Data());
                }
            }
            // labeling
            legendRatioRawClusterTrigger->Draw();
            PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, "#gamma candidates", detectionProcess, 42, 0.03, optionPeriod);

            canvasRatioRawClusterPt->Update();
            canvasRatioRawClusterPt->SaveAs(Form("%s/%s_TriggerYieldClusterRatio.%s",outputDir.Data(),prefix2.Data(),suffix.Data()));
            delete canvasRatioRawClusterPt;
        }
    }
        
        
    if (correctionFilesAvail){
        //*****************************************************************************************
        //******************* Compare Corrected Yields ********************************************
        //*****************************************************************************************
        // Define canvas
        TCanvas* canvasCorrectedYieldMeson = new TCanvas("canvasCorrectedYieldMeson","",1350,1500);  
        DrawGammaCanvasSettings( canvasCorrectedYieldMeson,  0.13, 0.02, 0.02, 0.09);
        // Define upper panel
        TPad* padCorrectedYield = new TPad("padCorrectedYield", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padCorrectedYield, 0.12, 0.02, 0.02, 0.);
        padCorrectedYield->SetLogy(1);
        padCorrectedYield->Draw();
        // Define lower panel
        TPad* padCorrectedYieldRatios = new TPad("padCorrectedYieldRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padCorrectedYieldRatios, 0.12, 0.02, 0.0, 0.2);
        padCorrectedYieldRatios->SetLogy(0);
        padCorrectedYieldRatios->Draw();

        // Plot corrected yield in upper panel
        padCorrectedYield->cd();            
        TLegend* legendCorrectedYieldMeson = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032); 
        if (cutVariationName.Contains("dEdxPi")){
            legendCorrectedYieldMeson->SetTextSize(0.02);
        }
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
                DrawAutoGammaMesonHistos( histoCorrectedYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)",
                                        kTRUE, 5., 5e-10,kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                DrawGammaSetMarker(histoCorrectedYieldCut[i], 20, 1., color[0], color[0]);
                histoCorrectedYieldCut[i]->DrawCopy("e1,p");
                legendCorrectedYieldMeson->AddEntry(histoCorrectedYieldCut[i], Form("standard: %s",cutStringsName[i].Data()));
            }
            else{
                if(i<20){
                    DrawGammaSetMarker(histoCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoCorrectedYieldCut[i]->DrawCopy("same,e1,p");
                legendCorrectedYieldMeson->AddEntry(histoCorrectedYieldCut[i], cutStringsName[i].Data());
            }
        }
        legendCorrectedYieldMeson->Draw();

        // labeling the plot
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);
        
        // plot ratio of corrected yields in lower panel    
        padCorrectedYieldRatios->cd();
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.8;
                Double_t maxYRatio = 1.2;
//                 if(optionMult.CompareTo("Mult")==0) {
//                     if( optionEnergy.Contains("Pb")){ 
//                         minYRatio = 0.005;        maxYRatio = 2.1;
//                     } else {
//                         minYRatio = -1.05;        maxYRatio = 12.5;
//                     }
//                 } else if (cutVariationName.Contains("Cent")) {
//                     minYRatio = 0.005;      maxYRatio = 3.2;
//                 } else if (cutVariationName.Contains("RCutAndPhotonQuality") && !optionEnergy.Contains("PbPb")) {
//                     minYRatio = 0.65;          maxYRatio = 1.15;
//                 } else if (cutVariationName.Contains("ClusterTrackMatching")){
//                     minYRatio = 0.75;      maxYRatio = 1.25;
//                 }    
//                 if (mode != 0 && mode!= 1 ){
//                     minYRatio = 0.75;      maxYRatio = 1.25;
//                 }
                SetStyleHistoTH1ForGraphs(histoRatioCorrectedYieldCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20, 1.,color[0],color[0]);
                histoRatioCorrectedYieldCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                
                for(Int_t b = 0; b< histoRatioCorrectedYieldCut[i]->GetNbinsX(); b++){
                    histoRatioCorrectedYieldCut[i]->SetBinError(b+1,histoCorrectedYieldCut[i]->GetBinError(b+1)/histoCorrectedYieldCut[i]->GetBinContent(b+1));
                }
                histoRatioCorrectedYieldCut[i]->SetFillColor(kGray+2);
                histoRatioCorrectedYieldCut[i]->SetFillStyle(0);
                histoRatioCorrectedYieldCut[i]->DrawCopy("p,e2");  
    
//                 histoRatioCorrectedYieldCut[i]->DrawCopy("p,e1");                
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioCorrectedYieldCut[i]->DrawCopy("same,e1,p");
            }
        }
        


        DrawGammaLines(0., maxPt,1., 1.,0.1);

        canvasCorrectedYieldMeson->Update();
        canvasCorrectedYieldMeson->SaveAs(Form("%s/%s_%s_CorrectedYield.%s",outputDir.Data(), meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasCorrectedYieldMeson;


        //**************************************************************************************
        //********************* Plotting efficiencies *********************************************
        //**************************************************************************************
        // Define canvas    
        TCanvas* canvasTrueEffiMeson = new TCanvas("canvasTrueEffiMeson","",1350,1500);  // gives the page size
        DrawGammaCanvasSettings( canvasTrueEffiMeson,  0.13, 0.02, 0.02, 0.09);
        // Define upper panel
        TPad* padTrueEffi = new TPad("padTrueEffi", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padTrueEffi, 0.12, 0.02, 0.04, 0.);
        padTrueEffi->Draw();
        // Define lower panel
        TPad* padTrueEffiRatios = new TPad("padTrueEffiRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padTrueEffiRatios, 0.12, 0.02, 0.0, 0.2);
        padTrueEffiRatios->Draw();

        // draw efficiency in upper panel
        padTrueEffi->cd();
//         if (mode == 2 || mode == 3 ) padTrueEffi->SetLogy(1);
//         else padTrueEffi->SetLogy(0);

        TLegend* legendEffiMeson = GetAndSetLegend2(0.15,0.92-1.15*0.032*NumberOfCuts,0.3,0.92, 1500*0.75*0.032); 
        if (cutVariationName.Contains("dEdxPi")){
            legendEffiMeson->SetTextSize(0.02);
        }
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
                DrawGammaSetMarker(histoTrueEffiCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoTrueEffiCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{%s}",textMeson.Data()),
                                        kTRUE, 5., 10e-10,kFALSE,
                                        kTRUE, -0.1, 0.00030,
                                        kFALSE, 0., 10.);
                if (mode == 9 || mode == 0 )histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(0.0,0.003);
                if (mode == 2 || mode == 3 )histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(0,0.1);
                if (mode == 4 || mode == 5 )histoTrueEffiCut[i]->GetYaxis()->SetRangeUser(0,0.6);
                histoTrueEffiCut[i]->DrawCopy("e1,p");
                legendEffiMeson->AddEntry(histoTrueEffiCut[i],Form("standard: %s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoTrueEffiCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoTrueEffiCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoTrueEffiCut[i]->DrawCopy("same,e1,p");
                legendEffiMeson->AddEntry(histoTrueEffiCut[i],cutStringsName[i].Data());
            }
            
        }
        legendEffiMeson->Draw();
    
        // Efficiency plot labeling
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.2, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);
                
        // Draw ratio of efficiencies in lower panel
        padTrueEffiRatios->cd();
        if( optionEnergy.Contains("Pb") ) padTrueEffiRatios->SetLogy(0);
        else padTrueEffiRatios->SetLogy(0);
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){      
                Double_t minYRatio = 0.8;
                Double_t maxYRatio = 1.2;
                if (cutVariationName.Contains("MultiplicityPP")){
                    minYRatio = 0.5;
                    maxYRatio = 2.7;
                    
                }    
//                 if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0){
//                     minYRatio = 0.05;        maxYRatio = 2.4;
//                 } else if (cutVariationName.Contains("PhotonQuality")){
//                     padTrueEffiRatios->SetLogy(1);
//                     minYRatio = 0.001;        maxYRatio = 2;
//                 }
                SetStyleHistoTH1ForGraphs(histoRatioTrueEffiCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioTrueEffiCut[i], 20, 1.,color[0],color[0]);
                histoRatioTrueEffiCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                
                for(Int_t b = 0; b< histoRatioTrueEffiCut[i]->GetNbinsX(); b++){
                    histoRatioTrueEffiCut[i]->SetBinError(b+1,histoTrueEffiCut[i]->GetBinError(b+1)/histoTrueEffiCut[i]->GetBinContent(b+1));
                }
                histoRatioTrueEffiCut[i]->SetFillColor(kGray+2);
                histoRatioTrueEffiCut[i]->SetFillStyle(0);
                histoRatioTrueEffiCut[i]->DrawCopy("p,e2");  

//                 histoRatioTrueEffiCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioTrueEffiCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioTrueEffiCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioTrueEffiCut[i]->DrawCopy("same,e1,p");
            }

            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasTrueEffiMeson->Update();
        canvasTrueEffiMeson->SaveAs(Form("%s/%s_%s_Efficiencies.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasTrueEffiMeson;
        

        //**************************************************************************************
        //********************* Plotting Acceptance *********************************************
        //**************************************************************************************
        // Define canvas    
        TCanvas* canvasAcceptanceMeson = new TCanvas("canvasAcceptanceMeson","",1350,1500);  // gives the page size
        DrawGammaCanvasSettings( canvasAcceptanceMeson,  0.13, 0.02, 0.02, 0.09);
        // Define upper panel
        TPad* padAcceptance = new TPad("padAcceptance", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padAcceptance, 0.12, 0.02, 0.04, 0.);
        padAcceptance->Draw();
        // Define lower panel
        TPad* padAcceptanceRatios = new TPad("padAcceptanceRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padAcceptanceRatios, 0.12, 0.02, 0.0, 0.2);
        padAcceptanceRatios->Draw();

        // draw efficiency in upper panel
        padAcceptance->cd();
//         if (mode == 2 || mode == 3 ) padAcceptance->SetLogy(1);
        padAcceptance->SetLogy(0);

        TLegend* legendAcceptMeson = GetAndSetLegend2(0.15,0.92-1.15*0.032*NumberOfCuts,0.3,0.92, 1500*0.75*0.032);
        if (cutVariationName.Contains("dEdxPi")){
            legendAcceptMeson->SetTextSize(0.02);
        }
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i == 0){
                DrawGammaSetMarker(histoAcceptanceCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoAcceptanceCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("A_{%s}",textMeson.Data()),
                                        kFALSE, 5., 10e-10,kFALSE,
                                        kFALSE, -0.1, 0.00030,
                                        kFALSE, 0., 10.);
                if (mode == 9 || mode == 0 )histoAcceptanceCut[i]->GetYaxis()->SetRangeUser(0.0,1.1);
                if (mode == 2 || mode == 3 )histoAcceptanceCut[i]->GetYaxis()->SetRangeUser(0,0.6);
                if (mode == 4 || mode == 5 )histoAcceptanceCut[i]->GetYaxis()->SetRangeUser(0,0.4);
                histoAcceptanceCut[i]->DrawCopy("e1,p");
                legendAcceptMeson->AddEntry(histoAcceptanceCut[i],Form("standard: %s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoAcceptanceCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoAcceptanceCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoAcceptanceCut[i]->DrawCopy("same,e1,p");
                legendAcceptMeson->AddEntry(histoAcceptanceCut[i],cutStringsName[i].Data());
            }
            
        }
        legendAcceptMeson->Draw();
    
        // Acceptance plot labeling
        PutProcessLabelAndEnergyOnPlot( 0.55, 0.20, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);
                
        // Draw ratio of efficiencies in lower panel
        padAcceptanceRatios->cd();
        if( optionEnergy.Contains("Pb") ) padAcceptanceRatios->SetLogy(0);
        else padAcceptanceRatios->SetLogy(0);
        for(Int_t i = 0; i< NumberOfCuts; i++){
            if(i==0){      
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55;
                if( optionEnergy.Contains("Pb") && optionMult.CompareTo("Mult") == 0){
                    minYRatio = 0.05;        maxYRatio = 2.4;
                } else if (cutVariationName.Contains("PhotonQuality")){
                    padAcceptanceRatios->SetLogy(1);
                    minYRatio = 0.001;        maxYRatio = 2;
                }
                SetStyleHistoTH1ForGraphs(histoRatioAcceptanceCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20, 1.,color[0],color[0]);
                histoRatioAcceptanceCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoRatioAcceptanceCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioAcceptanceCut[i]->DrawCopy("same,e1,p");
            }

            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasAcceptanceMeson->Update();
        canvasAcceptanceMeson->SaveAs(Form("%s/%s_%s_Acceptance.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasAcceptanceMeson;


        if (doMassRatio){
            //**************************************************************************************
            //********************* Plotting MassRatio *********************************************
            //**************************************************************************************
            // Define canvas    
            TCanvas* canvasMassRatioMeson = new TCanvas("canvasMassRatioMeson","",1350,1500);  // gives the page size
            DrawGammaCanvasSettings( canvasMassRatioMeson,  0.13, 0.02, 0.02, 0.09);
            // Define upper panel
            TPad* padMassRatio = new TPad("padMassRatio", "", 0., 0.25, 1., 1.,-1, -1, -2);
            DrawGammaPadSettings( padMassRatio, 0.12, 0.02, 0.04, 0.);
            padMassRatio->Draw();
            // Define lower panel
            TPad* padMassRatioRatios = new TPad("padMassRatioRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
            DrawGammaPadSettings( padMassRatioRatios, 0.12, 0.02, 0.0, 0.2);
            padMassRatioRatios->Draw();

            // draw efficiency in upper panel
            padMassRatio->cd();
    //         if (mode == 2 || mode == 3 ) padMassRatio->SetLogy(1);
            padMassRatio->SetLogy(0);

            TLegend* legendMassRatio = GetAndSetLegend2(0.15,0.92-1.15*0.032*NumberOfCuts,0.3,0.92, 1500*0.75*0.032);
            if (cutVariationName.Contains("dEdxPi")){
                legendMassRatio->SetTextSize(0.02);
            }
            for(Int_t i = 0; i< NumberOfCuts; i++){
                if(i == 0){
                    DrawGammaSetMarker(histoMassRatioCut[i], 20, 1., color[0], color[0]);
                    DrawAutoGammaMesonHistos( histoMassRatioCut[i],
                                            "", "#it{p}_{T} (GeV/#it{c})", Form("m_{%s, MC}/m_{%s, data}",textMeson.Data(),textMeson.Data()),
                                            kFALSE, 5., 10e-10,kFALSE,
                                            kFALSE, -0.1, 0.00030,
                                            kFALSE, 0., 10.);
                    histoMassRatioCut[i]->GetYaxis()->SetRangeUser(0.935,1.065);
                    histoMassRatioCut[i]->DrawCopy("e1,p");
                    legendMassRatio->AddEntry(histoMassRatioCut[i],Form("standard: %s",cutStringsName[i].Data()));
                } else {
                    if(i<20){
                        DrawGammaSetMarker(histoMassRatioCut[i], 20+i, 1.,color[i],color[i]);
                    } else {
                        DrawGammaSetMarker(histoMassRatioCut[i], 20+i, 1.,color[i-20],color[i-20]);
                    }
                    histoMassRatioCut[i]->DrawCopy("same,e1,p");
                    legendMassRatio->AddEntry(histoMassRatioCut[i],cutStringsName[i].Data());
                }
                
            }
            legendMassRatio->Draw();
            DrawGammaLines(0., maxPt,1., 1.,0.1);
            DrawGammaLines(0., maxPt,1.02, 1.02,0.1,kGray+1, 7);
            DrawGammaLines(0., maxPt,0.98, 0.98,0.1,kGray+1, 7);
            
            // MassRatio plot labeling
            PutProcessLabelAndEnergyOnPlot( 0.55, 0.20, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);
                    
            // Draw ratio of efficiencies in lower panel
            padMassRatioRatios->cd();
            if( optionEnergy.Contains("Pb") ) padMassRatioRatios->SetLogy(0);
            else padMassRatioRatios->SetLogy(0);
            for(Int_t i = 0; i< NumberOfCuts; i++){
                if(i==0){      
                    Double_t minYRatio = 0.95;
                    Double_t maxYRatio = 1.05;
                    SetStyleHistoTH1ForGraphs(histoRatioMassRatioCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                    DrawGammaSetMarker(histoRatioMassRatioCut[i], 20, 1.,color[0],color[0]);
                    histoRatioMassRatioCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                    histoRatioMassRatioCut[i]->DrawCopy("p,e1");
                } else{
                    if(i<20){
                        DrawGammaSetMarker(histoRatioMassRatioCut[i], 20+i, 1.,color[i],color[i]);
                    } else {
                        DrawGammaSetMarker(histoRatioMassRatioCut[i], 20+i, 1.,color[i-20],color[i-20]);
                    }
                    histoRatioMassRatioCut[i]->DrawCopy("same,e1,p");
                }

                DrawGammaLines(0., maxPt,1., 1.,0.1);
            }

            canvasMassRatioMeson->Update();
            canvasMassRatioMeson->SaveAs(Form("%s/%s_%s_MassRatio.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
            delete canvasMassRatioMeson;
        }
        
        if (isEta){
            //**************************************************************************************
            //********************* Plotting Eta to pi0 ratio  *************************************
            //**************************************************************************************
            // Define canvas    
            TCanvas* canvasEtaToPi0Meson = new TCanvas("canvasEtaToPi0Meson","",1350,1500);  // gives the page size
            DrawGammaCanvasSettings( canvasEtaToPi0Meson,  0.13, 0.02, 0.02, 0.09);
            // Define upper panel
            TPad* padEtaToPi0 = new TPad("padEtaToPi0", "", 0., 0.25, 1., 1.,-1, -1, -2);
            DrawGammaPadSettings( padEtaToPi0, 0.12, 0.02, 0.04, 0.);
            padEtaToPi0->Draw();
            // Define lower panel
            TPad* padEtaToPi0Ratios = new TPad("padEtaToPi0Ratios", "", 0., 0., 1., 0.25,-1, -1, -2);
            DrawGammaPadSettings( padEtaToPi0Ratios, 0.12, 0.02, 0.0, 0.2);
            padEtaToPi0Ratios->Draw();

            // draw efficiency in upper panel
            padEtaToPi0->cd();
    //         if (mode == 2 || mode == 3 ) padEtaToPi0->SetLogy(1);
            padEtaToPi0->SetLogy(0);

            TLegend* legendEtaToPi0 = GetAndSetLegend2(0.15,0.92-1.15*0.032*NumberOfCuts,0.3,0.92, 1500*0.75*0.032);
            if (cutVariationName.Contains("dEdxPi")){
                legendEtaToPi0->SetTextSize(0.02);
            }
            for(Int_t i = 0; i< NumberOfCuts; i++){
                if(i == 0){
                    DrawGammaSetMarker(histoEtaToPi0Cut[i], 20, 1., color[0], color[0]);
                    DrawAutoGammaMesonHistos( histoEtaToPi0Cut[i],
                                            "", "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}",
                                            kFALSE, 5., 10e-10,kFALSE,
                                            kFALSE, -0.1, 0.00030,
                                            kFALSE, 0., 10.);
                    histoEtaToPi0Cut[i]->GetYaxis()->SetRangeUser(0.,1.5);
                    histoEtaToPi0Cut[i]->DrawCopy("e1,p");
                    legendEtaToPi0->AddEntry(histoEtaToPi0Cut[i],Form("standard: %s",cutStringsName[i].Data()));
                } else {
                    if(i<20){
                        DrawGammaSetMarker(histoEtaToPi0Cut[i], 20+i, 1.,color[i],color[i]);
                    } else {
                        DrawGammaSetMarker(histoEtaToPi0Cut[i], 20+i, 1.,color[i-20],color[i-20]);
                    }
                    histoEtaToPi0Cut[i]->DrawCopy("same,e1,p");
                    legendEtaToPi0->AddEntry(histoEtaToPi0Cut[i],cutStringsName[i].Data());
                }
                
            }
            legendEtaToPi0->Draw();
            DrawGammaLines(0., maxPt,0.45, 0.45,0.1);
            DrawGammaLines(0., maxPt,0.55, 0.55,0.1,kGray+1, 7);
            DrawGammaLines(0., maxPt,0.35, 0.35 ,0.1,kGray+1, 7);
            
            // EtaToPi0 plot labeling
            PutProcessLabelAndEnergyOnPlot( 0.55, 0.20, 0.032, collisionSystem, process, detectionProcess, 42, 0.03, optionPeriod);
                    
            // Draw ratio of efficiencies in lower panel
            padEtaToPi0Ratios->cd();
            if( optionEnergy.Contains("Pb") ) padEtaToPi0Ratios->SetLogy(0);
            else padEtaToPi0Ratios->SetLogy(0);
            for(Int_t i = 0; i< NumberOfCuts; i++){
                if(i==0){      
                    Double_t minYRatio = 0.2;
                    Double_t maxYRatio = 1.8;
                    SetStyleHistoTH1ForGraphs(histoRatioEtaToPi0Cut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                    DrawGammaSetMarker(histoRatioEtaToPi0Cut[i], 20, 1.,color[0],color[0]);
                    histoRatioEtaToPi0Cut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                    histoRatioEtaToPi0Cut[i]->DrawCopy("p,e1");
                } else{
                    if(i<20){
                        DrawGammaSetMarker(histoRatioEtaToPi0Cut[i], 20+i, 1.,color[i],color[i]);
                    } else {
                        DrawGammaSetMarker(histoRatioEtaToPi0Cut[i], 20+i, 1.,color[i-20],color[i-20]);
                    }
                    histoRatioEtaToPi0Cut[i]->DrawCopy("same,e1,p");
                    cout<< cutStringsName[i].Data() << endl;
                    for (Int_t l = 1; l<histoRatioEtaToPi0Cut[i]->GetNbinsX()+1;l++ ){
                       cout << l << "\t"<< histoRatioEtaToPi0Cut[i]->GetBinCenter(l) << "\t" << histoRatioEtaToPi0Cut[i]->GetBinContent(l) << endl;
                    }
                }

                DrawGammaLines(0., maxPt,1., 1.,0.1);
            }

            canvasEtaToPi0Meson->Update();
            canvasEtaToPi0Meson->SaveAs(Form("%s/%s_%s_EtaToPi0.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
            delete canvasEtaToPi0Meson;
            
        }    
        
        //*************************************************************************************************
        //******************** Output of the systematic Error due to Signal extraction for Meson ************
        //*************************************************************************************************
        // Determine number of bins
        Int_t NBinsPt = histoCorrectedYieldCut[0]->GetNbinsX();
        const Int_t NBinstPtConst = NBinsPt+1;
        
        // Create array of bin boundaries
        Double_t  BinsXCenter[NBinstPtConst];
        Double_t  BinsXWidth[NBinstPtConst];
        BinsXCenter[0]                      = 0;
        BinsXWidth[0]                       = 0.;
        for (Int_t i = 1; i < NBinsPt +1; i++){
            BinsXCenter[i]                  = histoCorrectedYieldCut[0]->GetBinCenter(i);
            BinsXWidth[i]                   = histoCorrectedYieldCut[0]->GetBinWidth(i)/2.;
        }

        // Create array of Sys Err Objects and fill them
        SysErrorConversion SysErrCut[ConstNumberOfCuts][NBinstPtConst];
        SysErrorConversion SysErrCutRaw[ConstNumberOfCuts][NBinstPtConst];
        for (Int_t j = 0; j < NumberOfCuts; j++){
            for (Int_t i = 1; i < NBinsPt +1; i++){
                SysErrCut[j][i].value       = histoCorrectedYieldCut[j]->GetBinContent(i);
                SysErrCut[j][i].error       = histoCorrectedYieldCut[j]->GetBinError(i);
                SysErrCutRaw[j][i].value    = histoRawYieldCut[j]->GetBinContent(i);
                SysErrCutRaw[j][i].error    = histoRawYieldCut[j]->GetBinError(i);
            }
        }

        // Create Difference arrays
        Double_t DifferenceCut[ConstNumberOfCuts][NBinstPtConst];
        Double_t DifferenceErrorCut[ConstNumberOfCuts][NBinstPtConst];
        Double_t RelDifferenceCut[ConstNumberOfCuts][NBinstPtConst];
        Double_t RelDifferenceErrorCut[ConstNumberOfCuts][NBinstPtConst];
        Double_t RelDifferenceRawCut[ConstNumberOfCuts][NBinstPtConst];
            
        // Create largest difference array
        Double_t LargestDiffNeg[NBinstPtConst];
        Double_t LargestDiffPos[NBinstPtConst];
        Double_t LargestDiffErrorNeg[NBinstPtConst];
        Double_t LargestDiffErrorPos[NBinstPtConst];
        Double_t LargestDiffRelNeg[NBinstPtConst];
        Double_t LargestDiffRelPos[NBinstPtConst];
        Double_t LargestDiffRelErrorNeg[NBinstPtConst];
        Double_t LargestDiffRelErrorPos[NBinstPtConst];

        // Initialize all differences with 0
        for (Int_t j = 0; j < NumberOfCuts; j++){
            for ( Int_t i = 0; i < NBinstPtConst; i++) {
                DifferenceCut[j][i]=0.;
                DifferenceErrorCut[j][i]=0.;
                LargestDiffNeg[i]=0.;
                LargestDiffPos[i]=0.;
                LargestDiffErrorNeg[i]=0.;
                LargestDiffErrorPos[i]=0.;
                RelDifferenceCut[j][i]=0.;
                RelDifferenceRawCut[j][i]=0.;
                RelDifferenceErrorCut[j][i]=0.;
            }
        }


        // Calculate largest difference among cut variation 
        for(Int_t j = 1; j < NumberOfCuts; j++){
            for (Int_t i = 0; i < NBinsPt +1; i++){
                // Calculate difference (rel/abs) and error for corrected yield
                DifferenceCut[j][i] = SysErrCut[j][i].value - SysErrCut[0][i].value;
                DifferenceErrorCut[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCut[j][i].error,2)-TMath::Power(SysErrCut[0][i].error,2)));
                if(SysErrCut[0][i].value != 0){
                    RelDifferenceCut[j][i] = DifferenceCut[j][i]/SysErrCut[0][i].value*100. ;
                    RelDifferenceErrorCut[j][i] = DifferenceErrorCut[j][i]/SysErrCut[0][i].value*100. ;
                } else {
                    RelDifferenceCut[j][i] = -10000.;
                    RelDifferenceErrorCut[j][i] = 100. ;
                }
                // Calculate relativ difference for raw yield
                if(SysErrCutRaw[0][i].value != 0){
                    RelDifferenceRawCut[j][i] = (SysErrCutRaw[j][i].value - SysErrCutRaw[0][i].value)/SysErrCutRaw[0][i].value*100. ;
                } else {
                    RelDifferenceRawCut[j][i] = -10000.;
                }
                if (i == 0){
                    RelDifferenceRawCut[j][i] = 0.;
                    RelDifferenceCut[j][i] = 0. ;
                    RelDifferenceErrorCut[j][i] = 0. ;
                    DifferenceCut[j][i] = 0.;
                    DifferenceErrorCut[j][i] = 0.;
                }  
                    
                if(doBarlow){ 
                    // !!! => Careful, this is meant to be a cross check. If it has to be used for the syst errors
                    // the syste error macros has to be changed accordingly (the mean of pos and  neg dev cannot be used)
                    
                    // Calculate largest differences in positiv and negative direction
                    if(DifferenceCut[j][i] < 0){ // largest negativ deviation
                        // Take deviation if larger than previous largest deviation 
                        // and relative raw yield loss less than 75%
                        if (TMath::Abs(LargestDiffNeg[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                            if( TMath::Abs(DifferenceCut[j][i])/TMath::Abs(DifferenceErrorCut[j][i]) > 1. ){
                                LargestDiffNeg[i] = DifferenceCut[j][i];
                                LargestDiffErrorNeg[i] = DifferenceErrorCut[j][i];
                            } else if( TMath::Abs(DifferenceCut[j][i])/TMath::Abs(DifferenceErrorCut[j][i]) < 1.){
                                cout << "Largest negative difference not updated" << endl;
                            }
                        }
                    } else { // largest positive deviation
                        // Take deviation if larger than previous largest deviation 
                        // and relative raw yield loss less than 75%
                        if (TMath::Abs(LargestDiffPos[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                            if( TMath::Abs(DifferenceCut[j][i])/TMath::Abs(DifferenceErrorCut[j][i]) > 1. ){
                                LargestDiffPos[i] = DifferenceCut[j][i];
                                LargestDiffErrorPos[i] = DifferenceErrorCut[j][i];
                            } else if( TMath::Abs(DifferenceCut[j][i])/TMath::Abs(DifferenceErrorCut[j][i]) < 1.){
                                cout << "Largest positive difference not updated" << endl;
                            }
                        }
                    }
                    
                } else {
                    
                    // Calculate largest differences in positiv and negative direction
                    if(DifferenceCut[j][i] < 0){ // largest negativ deviation
                        // Take deviation if larger than previous largest deviation 
                        // and relative raw yield loss less than 75%
                        if (TMath::Abs(LargestDiffNeg[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                            LargestDiffNeg[i] = DifferenceCut[j][i];
                            LargestDiffErrorNeg[i] = DifferenceErrorCut[j][i];
                        }
                    } else { // largest positive deviation
                        // Take deviation if larger than previous largest deviation 
                        // and relative raw yield loss less than 75%
                        if (TMath::Abs(LargestDiffPos[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                            LargestDiffPos[i] = DifferenceCut[j][i];
                            LargestDiffErrorPos[i] = DifferenceErrorCut[j][i];
                        }
                    }
                    if (i == 0){
                        LargestDiffPos[i] = 0;
                        LargestDiffErrorPos[i] = 0;
                        LargestDiffNeg[i] = 0;
                        LargestDiffErrorNeg[i] = 0;                      
                    }
                }
            }
        }
        
        if(doBarlow){ 
            
            TString SysErrCheckname = Form("%s/%s_%s_SystematicErrorBarlowCheck.dat",outputDir.Data(),meson.Data(),prefix2.Data());
            fstream SysErrDatCheck;
            SysErrDatCheck.open(SysErrCheckname.Data(), ios::out);
            SysErrDatCheck << "Barlow check for the systematic error" << endl;
            for (Int_t l=0; l< NumberOfCuts; l++){
                if (l == 0) {
                    SysErrDatCheck << endl <<"Bin" << "\t" << cutNumber[l] << "\t" <<endl;
                    for(Int_t i = 1; i < (NBinsPt +1); i++){
                        SysErrDatCheck << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << endl;    
                    }
                } else{
                    for(Int_t i = 1; i < (NBinsPt +1); i++){
                        SysErrDatCheck << endl <<"Cut" << "\t" << cutNumber[l] << "\t" <<endl;
                        SysErrDatCheck << "\t Barlow check for " << BinsXCenter[i] << endl;
                        SysErrDatCheck << "Delta = |a_1 - a_2| \t" << TMath::Abs(DifferenceCut[l][i]) << endl;
                        SysErrDatCheck << "sigma_Delta = sqrt( |sigma_2^2 - sigma_1^2| ) \t" << TMath::Abs(DifferenceErrorCut[l][i]) << endl;
                        SysErrDatCheck << "Check: Delta/sigma_Delta < 1? " << (TMath::Abs(DifferenceCut[l][i]))/(TMath::Abs(DifferenceErrorCut[l][i])) << endl;
                    }
                }
            }            
            SysErrDatCheck.close();
        }                        
                    
        // Write systematic error input to log file
        TString SysErrDatname = Form("%s/%s_%s_SystematicErrorCutStudies.dat",outputDir.Data(),meson.Data(),prefix2.Data());
        fstream SysErrDat;
        SysErrDat.open(SysErrDatname.Data(), ios::out);
        SysErrDat << "Calculation of the systematic error due to the yield cuts" << endl;

        cout << "works" << endl;
        for (Int_t l=0; l< NumberOfCuts; l++){
            if (l == 0) {
                SysErrDat << endl <<"Bin" << "\t" << cutNumber[l] << "\t" <<endl;
                for(Int_t i = 1; i < (NBinsPt +1); i++){
                SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << endl;    
                }
            } else{
                SysErrDat << endl <<"Bin" << "\t" << cutNumber[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
                for(Int_t i = 1; i < (NBinsPt +1); i++){
                    if (RelDifferenceRawCut[l][i] > -75.){
                        SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << "\t" <<  DifferenceCut[l][i] << "\t"<< DifferenceErrorCut[l][i] << "\t" <<
                        RelDifferenceCut[l][i] <<  "\t" << RelDifferenceErrorCut[l][i] <<"\t" << RelDifferenceRawCut[l][i]<< endl;
                    } else {
                        SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << "\t" <<  DifferenceCut[l][i] << "\t"<< DifferenceErrorCut[l][i] << "\t" << 
                        RelDifferenceCut[l][i] <<  "\t" << RelDifferenceErrorCut[l][i] <<"\t" << RelDifferenceRawCut[l][i]  <<"\t not considered in largest dev" <<endl;
                    }
                }
            }
        }
        SysErrDat << endl;
        SysErrDat << endl;
        SysErrDat << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
        for(Int_t i = 1; i < (NBinsPt +1); i++){
            SysErrDat << BinsXCenter[i]  << "\t" << LargestDiffNeg[i] << "\t" <<LargestDiffErrorNeg[i]<< "\t" << LargestDiffPos[i] << "\t" << LargestDiffErrorPos[i]<<endl;
        }
        SysErrDat << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
        // Calculate largest relative deviations
        for(Int_t i = 0; i < (NBinsPt +1); i++){
            if ( SysErrCut[0][i].value != 0.){
                LargestDiffRelNeg[i] = - LargestDiffNeg[i]/SysErrCut[0][i].value*100.;
                LargestDiffRelPos[i] = LargestDiffPos[i]/SysErrCut[0][i].value*100.;
                LargestDiffRelErrorNeg[i] = - LargestDiffErrorNeg[i]/SysErrCut[0][i].value*100.;
                LargestDiffRelErrorPos[i] = LargestDiffErrorPos[i]/SysErrCut[0][i].value*100.;
                if (i > 0){
                    SysErrDat << BinsXCenter[i] << "\t" << LargestDiffNeg[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffErrorNeg[i]/SysErrCut[0][i].value*100. << "\t" <<     
                    LargestDiffPos[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffErrorPos[i]/SysErrCut[0][i].value*100.<<endl;
                } else {
                    LargestDiffRelNeg[i] = 0.;
                    LargestDiffRelPos[i] = 0.;
                    LargestDiffRelErrorNeg[i] = 0.;
                    LargestDiffRelErrorPos[i] = 0.;
                }  
            } else {
                LargestDiffRelNeg[i] = 0.;
                LargestDiffRelPos[i] = 0.;
                LargestDiffRelErrorNeg[i] = 0.;
                LargestDiffRelErrorPos[i] = 0.;
            }
        }
        
        SysErrDat.close();

        // Create sys-err graphs
        TGraphAsymmErrors* SystErrGraphNeg = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffRelNeg, BinsXWidth, BinsXWidth, LargestDiffRelErrorNeg, LargestDiffRelErrorNeg);
        SystErrGraphNeg->SetName(Form("%s_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()));
        cout << "Negative error graph" << endl; 
        SystErrGraphNeg->Print();
        TGraphAsymmErrors* SystErrGraphPos = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffRelPos, BinsXWidth, BinsXWidth, LargestDiffRelErrorPos, LargestDiffRelErrorPos);
        SystErrGraphPos->SetName(Form("%s_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()));
        cout << "positive error graph" << endl;
        SystErrGraphPos->Print();
        
        TGraphAsymmErrors* SystErrGraphNegEtaToPi0 = NULL;
        TGraphAsymmErrors* SystErrGraphPosEtaToPi0 = NULL;
        
        if (isEta){
            //*************************************************************************************************
            //******************** Output of the systematic Error due to Signal extraction for Meson ************
            //*************************************************************************************************
            // Determine number of bins
            Int_t NBinsPtEtaToPi0                        = histoEtaToPi0Cut[0]->GetNbinsX();
            const Int_t NBinstPtConstEtaToPi0            = NBinsPtEtaToPi0+1;
            
            // Create array of bin boundaries
            Double_t  BinsXCenterEtaToPi0[NBinstPtConstEtaToPi0];
            Double_t  BinsXWidthEtaToPi0[NBinstPtConstEtaToPi0];
            BinsXCenterEtaToPi0[0]                      = 0;
            BinsXWidthEtaToPi0[0]                       = 0.;
            for (Int_t i = 1; i < NBinsPtEtaToPi0 +1; i++){
                BinsXCenterEtaToPi0[i]                  = histoEtaToPi0Cut[0]->GetBinCenter(i);
                BinsXWidthEtaToPi0[i]                   = histoEtaToPi0Cut[0]->GetBinWidth(i)/2.;
            }

            // Create array of Sys Err Objects and fill them
            SysErrorConversion SysErrCutEtaToPi0[ConstNumberOfCuts][NBinstPtConstEtaToPi0];
            SysErrorConversion SysErrCutEtaToPi0Raw[ConstNumberOfCuts][NBinstPtConstEtaToPi0];
            for (Int_t j = 0; j < NumberOfCuts; j++){
                for (Int_t i = 1; i < NBinsPtEtaToPi0 +1; i++){
                    SysErrCutEtaToPi0[j][i].value       = histoEtaToPi0Cut[j]->GetBinContent(i);
                    SysErrCutEtaToPi0[j][i].error       = histoEtaToPi0Cut[j]->GetBinError(i);
                    SysErrCutEtaToPi0Raw[j][i].value    = histoRawYieldCut[j]->GetBinContent(i);
                    SysErrCutEtaToPi0Raw[j][i].error    = histoRawYieldCut[j]->GetBinError(i);
                }
            }

            // Create Difference arrays
            Double_t DifferenceEtaToPi0Cut[ConstNumberOfCuts][NBinstPtConstEtaToPi0];
            Double_t DifferenceEtaToPi0ErrorCut[ConstNumberOfCuts][NBinstPtConstEtaToPi0];
            Double_t RelDifferenceEtaToPi0Cut[ConstNumberOfCuts][NBinstPtConstEtaToPi0];
            Double_t RelDifferenceEtaToPi0ErrorCut[ConstNumberOfCuts][NBinstPtConstEtaToPi0];
            Double_t RelDifferenceEtaToPi0RawCut[ConstNumberOfCuts][NBinstPtConstEtaToPi0];
                
            // Create largest difference array
            Double_t LargestDiffEtaToPi0Neg[NBinstPtConstEtaToPi0];
            Double_t LargestDiffEtaToPi0Pos[NBinstPtConstEtaToPi0];
            Double_t LargestDiffEtaToPi0ErrorNeg[NBinstPtConstEtaToPi0];
            Double_t LargestDiffEtaToPi0ErrorPos[NBinstPtConstEtaToPi0];
            Double_t LargestDiffEtaToPi0RelNeg[NBinstPtConstEtaToPi0];
            Double_t LargestDiffEtaToPi0RelPos[NBinstPtConstEtaToPi0];
            Double_t LargestDiffEtaToPi0RelErrorNeg[NBinstPtConstEtaToPi0];
            Double_t LargestDiffEtaToPi0RelErrorPos[NBinstPtConstEtaToPi0];

            // Initialize all differences with 0
            for (Int_t j = 0; j < NumberOfCuts; j++){
                for ( Int_t i = 0; i < NBinstPtConstEtaToPi0; i++) {
                    DifferenceEtaToPi0Cut[j][i]         = 0.;
                    DifferenceEtaToPi0ErrorCut[j][i]    = 0.;
                    LargestDiffEtaToPi0Neg[i]           = 0.;
                    LargestDiffEtaToPi0Pos[i]           = 0.;
                    LargestDiffEtaToPi0ErrorNeg[i]      = 0.;
                    LargestDiffEtaToPi0ErrorPos[i]      = 0.;
                    RelDifferenceEtaToPi0Cut[j][i]      = 0.;
                    RelDifferenceEtaToPi0RawCut[j][i]   = 0.;
                    RelDifferenceEtaToPi0ErrorCut[j][i] = 0.;
                }
            }


            // Calculate largest difference among cut variation 
            for(Int_t j = 1; j < NumberOfCuts; j++){
                for (Int_t i = 0; i < NBinsPtEtaToPi0 +1; i++){
                    // Calculate difference (rel/abs) and error for corrected yield
                    DifferenceEtaToPi0Cut[j][i]         = SysErrCutEtaToPi0[j][i].value - SysErrCutEtaToPi0[0][i].value;
                    DifferenceEtaToPi0ErrorCut[j][i]    = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutEtaToPi0[j][i].error,2)-TMath::Power(SysErrCutEtaToPi0[0][i].error,2)));
                    if(SysErrCutEtaToPi0[0][i].value != 0){
                        RelDifferenceEtaToPi0Cut[j][i]      = DifferenceEtaToPi0Cut[j][i]/SysErrCutEtaToPi0[0][i].value*100. ;
                        RelDifferenceEtaToPi0ErrorCut[j][i] = DifferenceEtaToPi0ErrorCut[j][i]/SysErrCutEtaToPi0[0][i].value*100. ;
                    } else {
                        RelDifferenceEtaToPi0Cut[j][i]      = -10000.;
                        RelDifferenceEtaToPi0ErrorCut[j][i] = 100. ;
                    }
                    // Calculate relativ difference for raw yield
                    if(SysErrCutEtaToPi0Raw[0][i].value != 0){
                        RelDifferenceEtaToPi0RawCut[j][i]   = (SysErrCutEtaToPi0Raw[j][i].value - SysErrCutEtaToPi0Raw[0][i].value)/SysErrCutEtaToPi0Raw[0][i].value*100. ;
                    } else {
                        RelDifferenceEtaToPi0RawCut[j][i]   = -10000.;
                    }
                    
                    if (i == 0){
                        DifferenceEtaToPi0Cut[j][i]         = 0.;
                        DifferenceEtaToPi0ErrorCut[j][i]    = 0;
                        RelDifferenceEtaToPi0Cut[j][i]      = 0.;
                        RelDifferenceEtaToPi0ErrorCut[j][i] = 0;
                        RelDifferenceEtaToPi0RawCut[j][i]   = 0;                      
                    }  
                    
                    if(doBarlow){ 
                        // !!! => Careful, this is meant to be a cross check. If it has to be used for the syst errors
                        // the syste error macros has to be changed accordingly (the mean of pos and  neg dev cannot be used)
                        
                        // Calculate largest differences in positiv and negative direction
                        if(DifferenceEtaToPi0Cut[j][i] < 0){ // largest negativ deviation
                            // Take deviation if larger than previous largest deviation 
                            // and relative raw yield loss less than 75%
                            if (TMath::Abs(LargestDiffEtaToPi0Neg[i]) < TMath::Abs(DifferenceEtaToPi0Cut[j][i]) && RelDifferenceEtaToPi0RawCut[j][i] > -75.){
                                if( TMath::Abs(DifferenceEtaToPi0Cut[j][i])/TMath::Abs(DifferenceEtaToPi0ErrorCut[j][i]) > 1. ){
                                    LargestDiffEtaToPi0Neg[i]       = DifferenceEtaToPi0Cut[j][i];
                                    LargestDiffEtaToPi0ErrorNeg[i]  = DifferenceEtaToPi0ErrorCut[j][i];
                                } else if( TMath::Abs(DifferenceEtaToPi0Cut[j][i])/TMath::Abs(DifferenceEtaToPi0ErrorCut[j][i]) < 1.){
                                    cout << "Largest negative difference not updated" << endl;
                                }
                            }
                        } else { // largest positive deviation
                            // Take deviation if larger than previous largest deviation 
                            // and relative raw yield loss less than 75%
                            if (TMath::Abs(LargestDiffEtaToPi0Pos[i]) < TMath::Abs(DifferenceEtaToPi0Cut[j][i]) && RelDifferenceEtaToPi0RawCut[j][i] > -75.){
                                if( TMath::Abs(DifferenceEtaToPi0Cut[j][i])/TMath::Abs(DifferenceEtaToPi0ErrorCut[j][i]) > 1. ){
                                    LargestDiffEtaToPi0Pos[i]       = DifferenceEtaToPi0Cut[j][i];
                                    LargestDiffEtaToPi0ErrorPos[i]  = DifferenceEtaToPi0ErrorCut[j][i];
                                } else if( TMath::Abs(DifferenceEtaToPi0Cut[j][i])/TMath::Abs(DifferenceEtaToPi0ErrorCut[j][i]) < 1.){
                                    cout << "Largest positive difference not updated" << endl;
                                }
                            }
                        }
                        
                    } else {
                        
                        // Calculate largest differences in positiv and negative direction
                        if(DifferenceEtaToPi0Cut[j][i] < 0){ // largest negativ deviation
                            // Take deviation if larger than previous largest deviation 
                            // and relative raw yield loss less than 75%
                            if (TMath::Abs(LargestDiffEtaToPi0Neg[i]) < TMath::Abs(DifferenceEtaToPi0Cut[j][i]) && RelDifferenceEtaToPi0RawCut[j][i] > -75.){
                                LargestDiffEtaToPi0Neg[i]           = DifferenceEtaToPi0Cut[j][i];
                                LargestDiffEtaToPi0ErrorNeg[i]      = DifferenceEtaToPi0ErrorCut[j][i];
                            }
                        } else { // largest positive deviation
                            // Take deviation if larger than previous largest deviation 
                            // and relative raw yield loss less than 75%
                            if (TMath::Abs(LargestDiffEtaToPi0Pos[i]) < TMath::Abs(DifferenceEtaToPi0Cut[j][i]) && RelDifferenceEtaToPi0RawCut[j][i] > -75.){
                                LargestDiffEtaToPi0Pos[i]           = DifferenceEtaToPi0Cut[j][i];
                                LargestDiffEtaToPi0ErrorPos[i]      = DifferenceEtaToPi0ErrorCut[j][i];
                            }
                        }
                        if (i == 0){
                            LargestDiffEtaToPi0Pos[i]           = 0.;
                            LargestDiffEtaToPi0ErrorPos[i]      = 0.;
                            LargestDiffEtaToPi0Neg[i]           = 0.;
                            LargestDiffEtaToPi0ErrorNeg[i]      = 0.;
                        }  

                    }
                }
            }
            
            if(doBarlow){ 
                
                TString SysErrCheckEtaToPi0name     = Form("%s/EtaToPi0_%s_SystematicErrorBarlowCheck.dat",outputDir.Data(),prefix2.Data());
                fstream SysErrDatEtaToPi0Check;
                SysErrDatEtaToPi0Check.open(SysErrCheckEtaToPi0name.Data(), ios::out);
                SysErrDatEtaToPi0Check << "Barlow check for the systematic error" << endl;
                for (Int_t l=0; l< NumberOfCuts; l++){
                    if (l == 0) {
                        SysErrDatEtaToPi0Check << endl <<"Bin" << "\t" << cutNumber[l] << "\t" <<endl;
                        for(Int_t i = 1; i < (NBinsPtEtaToPi0 +1); i++){
                            SysErrDatEtaToPi0Check << BinsXCenterEtaToPi0[i] << "\t" << SysErrCutEtaToPi0[l][i].value << "\t" << SysErrCutEtaToPi0[l][i].error << endl;    
                        }
                    } else{
                        for(Int_t i = 1; i < (NBinsPtEtaToPi0 +1); i++){
                            SysErrDatEtaToPi0Check << endl <<"Cut" << "\t" << cutNumber[l] << "\t" <<endl;
                            SysErrDatEtaToPi0Check << "\t Barlow check for " << BinsXCenterEtaToPi0[i] << endl;
                            SysErrDatEtaToPi0Check << "Delta = |a_1 - a_2| \t" << TMath::Abs(DifferenceEtaToPi0Cut[l][i]) << endl;
                            SysErrDatEtaToPi0Check << "sigma_Delta = sqrt( |sigma_2^2 - sigma_1^2| ) \t" << TMath::Abs(DifferenceEtaToPi0ErrorCut[l][i]) << endl;
                            SysErrDatEtaToPi0Check << "Check: Delta/sigma_Delta < 1? " << (TMath::Abs(DifferenceEtaToPi0Cut[l][i]))/(TMath::Abs(DifferenceEtaToPi0ErrorCut[l][i])) << endl;
                        }
                    }
                }            
                SysErrDatEtaToPi0Check.close();
            }                        
                        
            // Write systematic error input to log file
            TString SysErrDatEtaToPi0name = Form("%s/EtaToPi0_%s_SystematicErrorCutStudies.dat",outputDir.Data(),prefix2.Data());
            fstream SysErrDatEtaToPi0;
            SysErrDatEtaToPi0.open(SysErrDatEtaToPi0name.Data(), ios::out);
            SysErrDatEtaToPi0 << "Calculation of the systematic error due to the yield cuts" << endl;

            cout << "works" << endl;
            for (Int_t l=0; l< NumberOfCuts; l++){
                if (l == 0) {
                    SysErrDatEtaToPi0 << endl <<"Bin" << "\t" << cutNumber[l] << "\t" <<endl;
                    for(Int_t i = 1; i < (NBinsPtEtaToPi0 +1); i++){
                    SysErrDatEtaToPi0 << BinsXCenterEtaToPi0[i] << "\t" << SysErrCutEtaToPi0[l][i].value << "\t" << SysErrCutEtaToPi0[l][i].error << endl;    
                    }
                } else{
                    SysErrDatEtaToPi0 << endl <<"Bin" << "\t" << cutNumber[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
                    for(Int_t i = 1; i < (NBinsPtEtaToPi0 +1); i++){
                        if (RelDifferenceEtaToPi0RawCut[l][i] > -75.){
                            SysErrDatEtaToPi0 << BinsXCenterEtaToPi0[i] << "\t" << SysErrCutEtaToPi0[l][i].value << "\t" << SysErrCutEtaToPi0[l][i].error << "\t" 
                                              << DifferenceEtaToPi0Cut[l][i] << "\t"<< DifferenceEtaToPi0ErrorCut[l][i] << "\t" 
                                              << RelDifferenceEtaToPi0Cut[l][i] <<  "\t" << RelDifferenceEtaToPi0ErrorCut[l][i] <<"\t" << RelDifferenceEtaToPi0RawCut[l][i]<< endl;
                        } else {
                            SysErrDatEtaToPi0 << BinsXCenterEtaToPi0[i] << "\t" << SysErrCutEtaToPi0[l][i].value << "\t" << SysErrCutEtaToPi0[l][i].error << "\t" 
                                              << DifferenceEtaToPi0Cut[l][i] << "\t"<< DifferenceEtaToPi0ErrorCut[l][i] << "\t" 
                                              << RelDifferenceEtaToPi0Cut[l][i] <<  "\t" << RelDifferenceEtaToPi0ErrorCut[l][i] <<"\t" << RelDifferenceEtaToPi0RawCut[l][i]  
                                              << "\t not considered in largest dev" <<endl;
                        }
                    }
                }
            }
            SysErrDatEtaToPi0 << endl;
            SysErrDatEtaToPi0 << endl;
            SysErrDatEtaToPi0 << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
            for(Int_t i = 1; i < (NBinsPtEtaToPi0 +1); i++){
                SysErrDatEtaToPi0 << BinsXCenterEtaToPi0[i]  << "\t" << LargestDiffEtaToPi0Neg[i] << "\t" <<LargestDiffEtaToPi0ErrorNeg[i]<< "\t" 
                                  << LargestDiffEtaToPi0Pos[i] << "\t" << LargestDiffEtaToPi0ErrorPos[i]<<endl;
            }
            SysErrDatEtaToPi0 << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
            // Calculate largest relative deviations
            for(Int_t i = 1; i < (NBinsPtEtaToPi0 +1); i++){
                if ( SysErrCutEtaToPi0[0][i].value != 0.){
                    LargestDiffEtaToPi0RelNeg[i]        = - LargestDiffEtaToPi0Neg[i]/SysErrCutEtaToPi0[0][i].value*100.;
                    LargestDiffEtaToPi0RelPos[i]        = LargestDiffEtaToPi0Pos[i]/SysErrCutEtaToPi0[0][i].value*100.;
                    LargestDiffEtaToPi0RelErrorNeg[i]   = - LargestDiffEtaToPi0ErrorNeg[i]/SysErrCutEtaToPi0[0][i].value*100.;
                    LargestDiffEtaToPi0RelErrorPos[i]   = LargestDiffEtaToPi0ErrorPos[i]/SysErrCutEtaToPi0[0][i].value*100.;
                    if (i > 0){
                        SysErrDatEtaToPi0 << BinsXCenterEtaToPi0[i] << "\t" << LargestDiffEtaToPi0Neg[i]/SysErrCutEtaToPi0[0][i].value*100. << "\t" 
                                          << LargestDiffEtaToPi0ErrorNeg[i]/SysErrCutEtaToPi0[0][i].value*100. << "\t" 
                                          << LargestDiffEtaToPi0Pos[i]/SysErrCutEtaToPi0[0][i].value*100. << "\t" << LargestDiffEtaToPi0ErrorPos[i]/SysErrCutEtaToPi0[0][i].value*100. << endl;
                    } else {
                        LargestDiffEtaToPi0RelNeg[i]        = 0.;
                        LargestDiffEtaToPi0RelPos[i]        = 0.;
                        LargestDiffEtaToPi0RelErrorNeg[i]   = 0.;
                        LargestDiffEtaToPi0RelErrorPos[i]   = 0.;
                    }
                } else {
                    LargestDiffEtaToPi0RelNeg[i]        = 0.;
                    LargestDiffEtaToPi0RelPos[i]        = 0.;
                    LargestDiffEtaToPi0RelErrorNeg[i]   = 0.;
                    LargestDiffEtaToPi0RelErrorPos[i]   = 0.;
                }
            }
            SysErrDatEtaToPi0.close();

            // Create sys-err graphs
            SystErrGraphNegEtaToPi0 = new TGraphAsymmErrors( NBinsPtEtaToPi0+1, BinsXCenterEtaToPi0, LargestDiffEtaToPi0RelNeg, BinsXWidthEtaToPi0, BinsXWidthEtaToPi0, 
                                                                                LargestDiffEtaToPi0RelErrorNeg, LargestDiffEtaToPi0RelErrorNeg);
            SystErrGraphNeg->SetName(Form("EtaToPi0_SystErrorRelNeg_%s",cutVariationName.Data()));
            SystErrGraphPosEtaToPi0 = new TGraphAsymmErrors( NBinsPtEtaToPi0+1, BinsXCenterEtaToPi0, LargestDiffEtaToPi0RelPos, BinsXWidthEtaToPi0, BinsXWidthEtaToPi0, 
                                                                                LargestDiffEtaToPi0RelErrorPos, LargestDiffEtaToPi0RelErrorPos);
            SystErrGraphPos->SetName(Form("EtaToPi0_SystErrorRelPos_%s",cutVariationName.Data()));
   
        }
        
        // Write sys-err graph to root output file
        TString Outputname = Form("%s/%s_%s_SystematicErrorCuts.root",outputDirRootFile.Data(),meson.Data(),prefix2.Data());
        TFile* SystematicErrorFile = new TFile(Outputname.Data(),"UPDATE");

            SystErrGraphPos->Write(Form("%s_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
            SystErrGraphNeg->Write(Form("%s_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
            if (isEta){
                if (SystErrGraphPosEtaToPi0) SystErrGraphPosEtaToPi0->Write(Form("EtaToPi0_SystErrorRelPos_%s",cutVariationName.Data()),TObject::kOverwrite);
                if (SystErrGraphNegEtaToPi0) SystErrGraphNegEtaToPi0->Write(Form("EtaToPi0_SystErrorRelNeg_%s",cutVariationName.Data()),TObject::kOverwrite);
                if (systErrGraphNegYieldExtPi0EtaBinning) systErrGraphNegYieldExtPi0EtaBinning->Write(Form("Pi0EtaBinning_SystErrorRelNeg_YieldExtraction_%s", centralityString.Data()), 
                                                                                                      TObject::kOverwrite);
                if (systErrGraphPosYieldExtPi0EtaBinning) systErrGraphPosYieldExtPi0EtaBinning->Write(Form("Pi0EtaBinning_SystErrorRelPos_YieldExtraction_%s", centralityString.Data()), 
                                                                                                      TObject::kOverwrite);
            }    
            systErrGraphNegYieldExt->Write(Form("%s_SystErrorRelNeg_YieldExtraction_%s",meson.Data(),centralityString.Data()),TObject::kOverwrite);
            systErrGraphPosYieldExt->Write(Form("%s_SystErrorRelPos_YieldExtraction_%s",meson.Data(),centralityString.Data()),TObject::kOverwrite);
            if (systErrGraphBGEstimate) systErrGraphBGEstimate->Write(Form("%s_SystErrorRel_BGEstimate_%s",meson.Data(),centralityString.Data()),TObject::kOverwrite);
        SystematicErrorFile->Write();
        SystematicErrorFile->Close();
        
        // Write out Ratios for V0finder comparison
        if (cutVariationName.Contains("V0")){
            TString Outputname2 = Form("%s/%s_RatioRawYields.root",outputDir.Data(),meson.Data());
            TFile* RatioRawYieldsFile = new TFile(Outputname2.Data(),"UPDATE");
            histoCorrectedYieldCut[1]->Write(Form("%s_histoCorrectedYield_%s_%s",meson.Data(),cutVariationName.Data(),prefix2.Data()),TObject::kOverwrite);
            histoRatioCorrectedYieldCut[1]->Write(Form("%s_histoRatioCorrectedYield_%s_%s",meson.Data(),cutVariationName.Data(),prefix2.Data()),TObject::kOverwrite);
            histoRatioRawYieldCut[1]->Write(Form("%s_histoRatioRawYield_%s_%s",meson.Data(),cutVariationName.Data(),prefix2.Data()),TObject::kOverwrite);
            RatioRawYieldsFile->Write();
            RatioRawYieldsFile->Close();
        }
        Printf("File: %s/%s_%s_SystematicErrorCuts.root",outputDirRootFile.Data(),meson.Data(),prefix2.Data());
            
        //**************************************************************************************************************************
        //*******************************Start of optional pure photon part ********************************************************
        //**************************************************************************************************************************
        if (optGammaOn){
            
            // Define photon histogram/bool/file arrays
            TH1D*   histoCorrectedYieldGammaCut         [ConstNumberOfCuts];        // photon corrected yields 
            TH1D*   histoGammaEffiCut                   [ConstNumberOfCuts];        // photon efficiencies
            TH1D*   histoRawYieldGammaCut               [ConstNumberOfCuts];        // photon raw yields
            TH1D*   histoRatioCorrectedYieldGammaCut    [ConstNumberOfCuts];        // ratio of corrected photon yields
            TH1D*   histoRatioGammaEffiCut              [ConstNumberOfCuts];        // ratio of photon efficiencies
            TH1D*   histoRatioRawYieldGammaCut          [ConstNumberOfCuts];        // ratio of photon raw yields
            Bool_t  gammaAvail                          [ConstNumberOfCuts];        // flag if photon output is available for certain cut
            TString FileNameCorrectedGamma              [MaxNumberOfCuts];        // file name of corrected photon file
            TFile*  CutcorrfileGamma                    [ConstNumberOfCuts];        // files for different cuts
            
            // Reading gamma histograms
            for (Int_t i=0; i< NumberOfCuts; i++){
                FileNameCorrectedGamma[i] = Form("%s/%s/Gamma_%s_%s_GammaConvV1Correction_%s.root",cutNumberAdv[i].Data(),optionEnergy.Data(), meson.Data(),prefix2.Data(), cutNumber[i].Data());
                CutcorrfileGamma[i] = new TFile(FileNameCorrectedGamma[i]);
                // check if file exists
                if (CutcorrfileGamma[i]->IsZombie()){ 
                    gammaAvail[i] = kFALSE; 
                } else { 
                    cout<< FileNameCorrectedGamma[i] << "\t" << "exists" << endl;
                    gammaAvail[i] = kTRUE;
                }
                // only read the histograms if the file is there
                if (gammaAvail[i]){
                    histoCorrectedYieldGammaCut[i] =(TH1D*)CutcorrfileGamma[i]->Get("GammaUnfold");
                    histoCorrectedYieldGammaCut[i]->SetName(Form("CorrectedYieldGamma_%s",cutNumber[i].Data()));
                    histoGammaEffiCut[i] =(TH1D*)CutcorrfileGamma[i]->Get("MCGammaPrimaryRecoEffMCPt");
                    histoGammaEffiCut[i]->SetName(Form("MCGammaPrimaryRecoEffMCPt_%s",cutNumber[i].Data()));
                    histoRawYieldGammaCut[i] =(TH1D*)CutcorrfileGamma[i]->Get("RawGammaSpectrum");
                    histoRawYieldGammaCut[i]->SetName(Form("RawGammaSpectrum_%s",cutNumber[i].Data()));
                    histoRawYieldGammaCut[i]->Scale(1./nColls[i]);
                    if (gammaAvail[0]){
                        histoRatioCorrectedYieldGammaCut[i] = (TH1D*) histoCorrectedYieldGammaCut[i]->Clone(Form("histoRatioCorrectedYieldGammaCut_%s",cutNumber[i].Data()));
                        histoRatioCorrectedYieldGammaCut[i]->Divide(histoRatioCorrectedYieldGammaCut[i],histoCorrectedYieldGammaCut[0],1.,1.,"B");
                        histoRatioRawYieldGammaCut[i] = (TH1D*) histoRawYieldGammaCut[i]->Clone(Form("histoRawYieldGammaCut_%s",cutNumber[i].Data()));
                        histoRatioRawYieldGammaCut[i]->Divide(histoRatioRawYieldGammaCut[i],histoRawYieldGammaCut[0],1.,1.,"B");
                        histoRatioGammaEffiCut[i] = (TH1D*) histoGammaEffiCut[i]->Clone(Form("MCGammaPrimaryRecoEffMCPt_%s",cutNumber[i].Data()));
                        histoRatioGammaEffiCut[i]->Divide(histoRatioGammaEffiCut[i],histoGammaEffiCut[0],1.,1.,"B");
                    }
                    
                }    
            }
            cout<<"=========================="<<endl;
                
            //*****************************************************************************************
            //******************* Compare Corrected Photon Yields ********************************************
            //*****************************************************************************************
                // Define canvas
                TCanvas* canvasCorrectedYieldGamma = new TCanvas("canvasCorrectedYieldGamma","",1350,1500);  // gives the page size
                DrawGammaCanvasSettings( canvasCorrectedYieldGamma, 0.13, 0.02, 0.02, 0.1);
                canvasCorrectedYieldGamma->SetLogy(0);
                // Define upper panel 
                TPad* padCorrectedYieldGamma = new TPad("padCorrectedYieldGamma", "", 0., 0.33, 1., 1.,-1, -1, -2);
                DrawGammaPadSettings( padCorrectedYieldGamma, 0.12, 0.02, 0.02, 0.);
                padCorrectedYieldGamma->SetLogy(1);
                padCorrectedYieldGamma->Draw();
                // Define lower panel
                TPad* padCorrectedYieldGammaRatios = new TPad("padCorrectedYieldGammaRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
                DrawGammaPadSettings( padCorrectedYieldGammaRatios, 0.12, 0.02, 0., 0.2);
                padCorrectedYieldGammaRatios->Draw();
                padCorrectedYieldGammaRatios->SetLogy(0);
                
                // plot corrected photon yield in upper panel if histograms exist
                padCorrectedYieldGamma->cd();
                TLegend* legendCorrectedYieldGamma = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
                if (cutVariationName.Contains("dEdxPi")){
                    legendCorrectedYieldGamma->SetTextSize(0.02);
                }
                for(Int_t i = 0; i< NumberOfCuts; i++){
                    if(i == 0){
                        if (gammaAvail[i]) DrawAutoGammaMesonHistos( histoCorrectedYieldGammaCut[i],
                                                "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}^{#gamma}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)",
                                                kTRUE, 5., 5e-10,kTRUE,
                                                kFALSE, 0.0, 0.030,
                                                kFALSE, 0., 10.);
                        if (gammaAvail[i])DrawGammaSetMarker(histoCorrectedYieldGammaCut[i], 20, 1., color[0], color[0]);
                        if (gammaAvail[i])histoCorrectedYieldGammaCut[i]->DrawCopy("e1,p");
                        if (gammaAvail[i])legendCorrectedYieldGamma->AddEntry(histoCorrectedYieldGammaCut[i], Form("standard: %s",cutStringsName[i].Data()));
                    }
                    else{
                        if(i<20){
                            if (gammaAvail[i])DrawGammaSetMarker(histoCorrectedYieldGammaCut[i], 20+i, 1.,color[i],color[i]);
                        } else {
                            if (gammaAvail[i])DrawGammaSetMarker(histoCorrectedYieldGammaCut[i], 20+i, 1.,color[i-20],color[i-20]);
                        }
                        if (gammaAvail[i])histoCorrectedYieldGammaCut[i]->DrawCopy("same,e1,p");
                        if (gammaAvail[i])legendCorrectedYieldGamma->AddEntry(histoCorrectedYieldGammaCut[i], cutStringsName[i].Data());
                    }
                }
                legendCorrectedYieldGamma->Draw();

                // labeling of plot
                PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, "#gamma", detectionProcess, 42, 0.03, optionPeriod);
                    
                // plot ratios of corrected gammas    
                padCorrectedYieldGammaRatios->cd();
                for(Int_t i = 0; i< NumberOfCuts; i++){
                    if(i==0){
                        if (gammaAvail[0]){
                            SetStyleHistoTH1ForGraphs(histoRatioCorrectedYieldGammaCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                            DrawGammaSetMarker(histoRatioCorrectedYieldGammaCut[i], 20, 1.,color[0],color[0]);
                            histoRatioCorrectedYieldGammaCut[i]->GetYaxis()->SetRangeUser(0.85,1.15);
                            histoRatioCorrectedYieldGammaCut[i]->DrawCopy("p,e1");
                        }
                    } else {
                        if(i<20){
                        if (gammaAvail[0] && gammaAvail[i])    DrawGammaSetMarker(histoRatioCorrectedYieldGammaCut[i], 20+i, 1.,color[i],color[i]);
                        } else {
                            if (gammaAvail[0] && gammaAvail[i])    DrawGammaSetMarker(histoRatioCorrectedYieldGammaCut[i], 20+i, 1.,color[i-20],color[i-20]);
                        }
                        if (gammaAvail[0] && gammaAvail[i])    histoRatioCorrectedYieldGammaCut[i]->DrawCopy("same,e1,p");
                    }
                }
                DrawGammaLines(0., maxPt,1., 1.,0.1);

                canvasCorrectedYieldGamma->Update();
                canvasCorrectedYieldGamma->SaveAs(Form("%s/Gamma_%s_CorrectedYield.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
                delete canvasCorrectedYieldGamma;

            //*****************************************************************************************
            //******************* compare raw photon yields *******************************************
            //*****************************************************************************************
                // define canvas
                TCanvas* canvasRawYieldGamma = new TCanvas("canvasRawYieldGamma","",1350,1500);  // gives the page size
                DrawGammaCanvasSettings( canvasRawYieldGamma, 0.13, 0.02, 0.02, 0.1);
                canvasRawYieldGamma->SetLogy(0);
                // define upper panel
                TPad* padRawYieldGamma = new TPad("padRawYieldGamma", "", 0., 0.33, 1., 1.,-1, -1, -2);
                DrawGammaPadSettings( padRawYieldGamma, 0.12, 0.02, 0.02, 0.);
                padRawYieldGamma->SetLogy(1);
                padRawYieldGamma->Draw();
                // define lower panel
                TPad* padRawYieldGammaRatios = new TPad("padRawYieldGammaRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
                DrawGammaPadSettings( padRawYieldGammaRatios, 0.12, 0.02, 0., 0.2);
                padRawYieldGammaRatios->Draw();
                padRawYieldGammaRatios->SetLogy(0);
                
                // plot raw yields in upper panel
                padRawYieldGamma->cd();          
                TLegend* legendRawYieldGamma = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*NumberOfCuts, 1500*0.75*0.032);
                if (cutVariationName.Contains("dEdxPi")){
                    legendRawYieldGamma->SetTextSize(0.02);
                }
                for(Int_t i = 0; i< NumberOfCuts; i++){
                    if(i == 0){
                        if (gammaAvail[i]) DrawAutoGammaMesonHistos( histoRawYieldGammaCut[i],
                                                "", "#it{p}_{T} (GeV/#it{c})", "#gamma RAW Yield/(#it{N}_{ev} #it{N}_{coll})",
                                                kTRUE, 5., 5e-10,kTRUE,
                                                kFALSE, 0.0, 0.030,
                                                kFALSE, 0., 10.);
                        if (gammaAvail[i])DrawGammaSetMarker(histoRawYieldGammaCut[i], 20, 1., color[0], color[0]);
                        if (gammaAvail[i])histoRawYieldGammaCut[i]->DrawCopy("e1,p");
                        if (gammaAvail[i])legendRawYieldGamma->AddEntry(histoRawYieldGammaCut[i], Form("standard: %s",cutStringsName[i].Data()));
                    } else {
                        if(i<20){
                            if (gammaAvail[i])DrawGammaSetMarker(histoRawYieldGammaCut[i], 20+i, 1.,color[i],color[i]);
                        } else {
                            if (gammaAvail[i])DrawGammaSetMarker(histoRawYieldGammaCut[i], 20+i, 1.,color[i-20],color[i-20]);
                        }
                        if (gammaAvail[i])histoRawYieldGammaCut[i]->DrawCopy("same,e1,p");
                        if (gammaAvail[i])legendRawYieldGamma->AddEntry(histoRawYieldGammaCut[i], cutStringsName[i].Data());
                    }
                    
                }
                legendRawYieldGamma->Draw();
                
                // labeling  
                PutProcessLabelAndEnergyOnPlot( 0.55, 0.95, 0.032, collisionSystem, "#gamma", detectionProcess, 42, 0.03, optionPeriod);

                // draw ratios in lower panel    
                padRawYieldGammaRatios->cd();
                for(Int_t i = 0; i< NumberOfCuts; i++){
                    if(i==0){
                        if (gammaAvail[0]){
                            SetStyleHistoTH1ForGraphs(histoRatioRawYieldGammaCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                            DrawGammaSetMarker(histoRatioRawYieldGammaCut[i], 20, 1.,color[0],color[0]);
                            histoRatioRawYieldGammaCut[i]->GetYaxis()->SetRangeUser(0.05,1.55);
                            histoRatioRawYieldGammaCut[i]->DrawCopy("p,e1");
                        }
                    } else {
                        if(i<20){
                        if (gammaAvail[0] && gammaAvail[i])    DrawGammaSetMarker(histoRatioRawYieldGammaCut[i], 20+i, 1.,color[i],color[i]);
                        } else {
                            if (gammaAvail[0] && gammaAvail[i])    DrawGammaSetMarker(histoRatioRawYieldGammaCut[i], 20+i, 1.,color[i-20],color[i-20]);
                        }
                        if (gammaAvail[0] && gammaAvail[i])    histoRatioRawYieldGammaCut[i]->DrawCopy("same,e1,p");
                    }
                }
                DrawGammaLines(0., maxPt,1., 1.,0.1);

                canvasRawYieldGamma->Update();
                canvasRawYieldGamma->SaveAs(Form("%s/Gamma_%s_RawYield.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
                delete canvasRawYieldGamma;
                
            //*****************************************************************************************
            //******************* compare raw photon yields *******************************************
            //*****************************************************************************************
                // define canvas
                TCanvas* canvasEffiGamma = new TCanvas("canvasEffiGamma","",1350,1500);  // gives the page size
                DrawGammaCanvasSettings( canvasEffiGamma, 0.13, 0.02, 0.02, 0.1);
                canvasEffiGamma->SetLogy(0);
                // define upper panel
                TPad* padEffiGamma = new TPad("padEffiGamma", "", 0., 0.33, 1., 1.,-1, -1, -2);
                DrawGammaPadSettings( padEffiGamma, 0.12, 0.02, 0.02, 0.);
                padEffiGamma->Draw();
                // define lower panel
                TPad* padEffiGammaRatios = new TPad("padEffiGammaRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
                DrawGammaPadSettings( padEffiGammaRatios, 0.12, 0.02, 0., 0.2);
                padEffiGammaRatios->Draw();
                
                // draw photon efficiencies in upper panel
                padEffiGamma->cd();
                TLegend* legendEffiGamma = GetAndSetLegend2(0.15,0.92-1.15*0.032*NumberOfCuts,0.6,0.92, 1500*0.75*0.032);
                if (cutVariationName.Contains("dEdxPi")){
                    legendEffiGamma->SetTextSize(0.02);
                }
                for(Int_t i = 0; i< NumberOfCuts; i++){
                    if(i == 0){
                        if (gammaAvail[i]) DrawAutoGammaMesonHistos( histoGammaEffiCut[i],
                                                "", "#it{p}_{T, MC} (GeV/#it{c})", "#gamma Efficiency",
                                                kTRUE, 1.2, 0,kTRUE,
                                                kFALSE, 0.0, 0.030,
                                                kFALSE, 0., 10.);
                        if (gammaAvail[i])DrawGammaSetMarker(histoGammaEffiCut[i], 20, 1., color[0], color[0]);
                        if (gammaAvail[i])histoGammaEffiCut[i]->DrawCopy("e1,p");
                        if (gammaAvail[i])legendEffiGamma->AddEntry(histoGammaEffiCut[i], Form("standard: %s",cutStringsName[i].Data()));
                    } else {
                        if(i<20){
                            if (gammaAvail[i])DrawGammaSetMarker(histoGammaEffiCut[i], 20+i, 1.,color[i],color[i]);
                        } else {
                            if (gammaAvail[i])DrawGammaSetMarker(histoGammaEffiCut[i], 20+i, 1.,color[i-20],color[i-20]);
                        }
                        if (gammaAvail[i])histoGammaEffiCut[i]->DrawCopy("same,e1,p");
                        if (gammaAvail[i])legendEffiGamma->AddEntry(histoGammaEffiCut[i], cutStringsName[i].Data());
                    }
                    
                }
                legendEffiGamma->Draw();
                // labeling
                PutProcessLabelAndEnergyOnPlot( 0.55, 0.2, 0.032, collisionSystem, "#gamma", detectionProcess, 42, 0.03, optionPeriod);
                    
                // draw ratio of photon efficiencies in lower panel    
                padEffiGammaRatios->cd();
                for(Int_t i = 0; i< NumberOfCuts; i++){
                    if(i==0){
                        if (gammaAvail[0]){
                            SetStyleHistoTH1ForGraphs(histoRatioGammaEffiCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                            DrawGammaSetMarker(histoRatioGammaEffiCut[i], 20, 1.,color[0],color[0]);
                            histoRatioGammaEffiCut[i]->GetYaxis()->SetRangeUser(0.05,1.55);
                            histoRatioGammaEffiCut[i]->DrawCopy("p,e1");
                        }
                    } else {
                        if(i<20){
                            if (gammaAvail[0] && gammaAvail[i])    DrawGammaSetMarker(histoRatioGammaEffiCut[i], 20+i, 1.,color[i],color[i]);
                        } else {
                            if (gammaAvail[0] && gammaAvail[i])    DrawGammaSetMarker(histoRatioGammaEffiCut[i], 20+i, 1.,color[i-20],color[i-20]);
                        }
                        if (gammaAvail[0] && gammaAvail[i])    histoRatioGammaEffiCut[i]->DrawCopy("same,e1,p");
                    }
                }
                DrawGammaLines(0., maxPt,1., 1.,0.1);

                canvasEffiGamma->Update();
                canvasEffiGamma->SaveAs(Form("%s/Gamma_%s_Efficiencies.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
                delete canvasEffiGamma;        
        }
    }
    return;
}
