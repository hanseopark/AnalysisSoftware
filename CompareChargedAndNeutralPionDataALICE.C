/****************************************************************************************************************************
******        provided by Gamma Conversion Group, PWGGA,                                                                *****
******        Friederike Bock, friederike.bock@cern.ch                                                                  *****
*****************************************************************************************************************************/

#include <Riostream.h>
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
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

void CompareChargedAndNeutralPionDataALICE( TString suffix = "eps", 
                                            TString nameFilePP = "CombinedResultsPaperX.root", 
                                            TString nameFilePbPb = "CombinedResultsPbPb_15_May_2013.root", 
//                                             TString fileNamePCMEMCALpp2760GeVPP = "",
                                            TString fileNamePCMPHOSpp2760GeVPP = ""
                                          ){

    gROOT->Reset();    
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis(suffix);    
    SetPlotStyle();
    
    TString dateForOutput                   = ReturnDateStringForOutput();
    TString outputDir                       = Form("%s/%s/ComparisonNeutralAndChargedPions",suffix.Data(),dateForOutput.Data());
    TString fileNameChargedPionPbPb         = "ExternalInputPbPb/IdentifiedCharged/ChargedPionSpectraPbPb_4_Apr_2014.root";
    TString fileNameChargedPionPP           = "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_4_Apr_2014.root";
    TString fileNameEMCalPion7TeVPP         = "ExternalInput/EMCAL/7TeV/pi0Spectrum2011EMCALAddedSignalsEffic_7TeV_150323_evi_11cd.root";
//     TString fileNameEMCalPion7TeVPP      = "ExternalInput/EMCAL/7TeV/WeightedAve7TeV_coooked_20Mar2015.root";
//     TString fileNameEMCalPion2760GeVPP   = "ExternalInput/EMCAL/2.76TeV/FinalCombinedXsec_pi0276TeV_25Apr2015_11a13g.root";
    TString fileNameDalitz7TeV              = "ExternalInput/Dalitz/data_PCMDalitzResultsFullCorrection_PP_NoBinShiftin.root";
    TString fileNameDalitz2760GeV           = "ExternalInput/Dalitz/data_PCMDalitzResultsFullCorrection_PP_NoBinShifting_25102013.root";
    TString fileNameDalitzPbPb              = "ExternalInputPbPb/data_PCMDalitzResults_PbPb_2.76TeV_28102013.root";
    TString fileNameComb2760GeVPP           = "CombinedResultsPaperPP2760GeV_2016_02_15.root";
    
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputFileNeutralPionPP.root ",nameFilePP.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFileNeutralPionPPUpdated2760GeV.root ",fileNameComb2760GeVPP.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFileNeutralPionPbPb.root ",nameFilePbPb.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFileChargedPionPbPb.root ",fileNameChargedPionPbPb.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFileEMCalNeutralPion7TeVPP.root ",fileNameEMCalPion7TeVPP.Data(),outputDir.Data() ));
//     gSystem->Exec(Form("cp %s %s/InputFileEMCalNeutralPion2760GeVPP.root ",fileNameEMCalPion2760GeVPP.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFileDalitzNeutralPionPbPb.root ",fileNameDalitzPbPb.Data(),outputDir.Data() ));
    Bool_t enablePCMEMCALComp2760GeV    = kFALSE;
//     if (fileNamePCMEMCALpp2760GeVPP.CompareTo("")!= 0){
//         gSystem->Exec(Form("cp %s %s/InputFilePCMEMCalNeutralPion2760GeVPP.root ",fileNameEMCalPion2760GeVPP.Data(),outputDir.Data() ));
//         enablePCMEMCALComp2760GeV         = kTRUE;
//     }
    Bool_t enablePCMPHOSComp2760GeV     = kFALSE;
    if (fileNamePCMPHOSpp2760GeVPP.CompareTo("")!= 0){
        gSystem->Exec(Form("cp %s %s/InputFilePCMPHOSNeutralPion2760GeVPP.root ",fileNamePCMPHOSpp2760GeVPP.Data(),outputDir.Data() ));
        enablePCMPHOSComp2760GeV        = kTRUE;
    }
    
    Double_t xSection2760GeVpp          = 55.416*1e-3;
    Double_t xSection2760GeVErrpp       = 3.9;
    Double_t xSection2760GeVppINEL      = 62.8*1e9;
    Double_t xSection900GeVppINEL       = 52.5*1e9;
    Double_t xSection7TeVppINEL         = 73.2*1e9;    
    Double_t recalcBarn                 = 1e12; //NLO in pbarn!!!!
    
    TString collisionSystemPP           = "pp #sqrt{#it{s}} = 2.76 TeV";        
    TString collisionSystemCent0        = "0-5% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";        
    TString collisionSystemCent1        = "5-10% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";        
    TString collisionSystemCent2        = "10-20% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";        
    TString collisionSystemCent         = "0-20% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";        
    TString collisionSystemSemiCent     = "20-40% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";        
    TString collisionSystemSemiPer      = "40-60% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";        
    TString collisionSystemPer          = "60-80% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";        

    Style_t markerStyleCombLowPt        = GetDefaultMarkerStyleDiffDetectors("Comb", kFALSE);
    Style_t markerStyleCombHighPt       = GetDefaultMarkerStyleDiffDetectors("Comb", kFALSE);
    Style_t markerStylePCMLowPt         = GetDefaultMarkerStyleDiffDetectors("PCM", kFALSE);
    Style_t markerStylePCMHighPt        = GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);      
    Style_t markerStylePHOSLowPt        = GetDefaultMarkerStyleDiffDetectors("PCM", kFALSE);
    Style_t markerStylePHOSHighPt       = GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);      
    Style_t markerStylePCMBGLowPt       = GetDefaultMarkerStyleDiffDetectors("PCM", kTRUE);
    Style_t markerStylePCMBGHighPt      = GetDefaultMarkerStyleDiffDetectors("PHOS", kTRUE);      
    Style_t markerStylePHOSBGLowPt      = GetDefaultMarkerStyleDiffDetectors("PCM", kTRUE);
    Style_t markerStylePHOSBGHighPt     = GetDefaultMarkerStyleDiffDetectors("PHOS", kTRUE);      
    Style_t markerStyleEMCALLowPt       = GetDefaultMarkerStyleDiffDetectors("PCM", kFALSE);
    Style_t markerStyleEMCALHighPt      = GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);      
    Style_t markerStyleEMCALMergedHighPt= GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);      
    Style_t markerStyleDalitzLowPt      = GetDefaultMarkerStyleDiffDetectors("PCM", kFALSE);
    Style_t markerStyleDalitzHighPt     = GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);      
    Style_t markerStylePCMEMCALLowPt    = GetDefaultMarkerStyleDiffDetectors("PCM", kFALSE);
    Style_t markerStylePCMEMCALHighPt   = GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);      
    Style_t markerStylePCMPHOSLowPt     = GetDefaultMarkerStyleDiffDetectors("PCM", kFALSE);
    Style_t markerStylePCMPHOSHighPt    = GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);      
    
    Size_t markerSizeComparison = 1.5;

    Color_t colorCombLowPt              = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
    Color_t colorCombHighPt             = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
    Color_t colorPCMLowPt               = GetDefaultColorDiffDetectors("PCM", kFALSE, kFALSE, kFALSE);
    Color_t colorPCMHighPt              = GetDefaultColorDiffDetectors("PCM", kFALSE, kFALSE, kTRUE);
    Color_t colorPHOSLowPt              = GetDefaultColorDiffDetectors("PHOS", kFALSE, kFALSE, kFALSE);
    Color_t colorPHOSHighPt             = GetDefaultColorDiffDetectors("PHOS", kFALSE, kFALSE, kTRUE);
    Color_t colorEMCALLowPt             = GetDefaultColorDiffDetectors("EMCal", kFALSE, kFALSE, kFALSE);
    Color_t colorEMCALHighPt            = GetDefaultColorDiffDetectors("EMCal", kFALSE, kFALSE, kTRUE);
    Color_t colorEMCALMergedHighPt      = GetDefaultColorDiffDetectors("EMCal merged", kFALSE, kFALSE, kTRUE);
    Color_t colorDalitzLowPt            = GetDefaultColorDiffDetectors("Dalitz", kFALSE, kFALSE, kFALSE);
    Color_t colorDalitzHighPt           = GetDefaultColorDiffDetectors("Dalitz", kFALSE, kFALSE, kTRUE);
    Color_t colorPCMEMCALLowPt          = GetDefaultColorDiffDetectors("PCM-EMCal", kFALSE, kFALSE, kFALSE);
    Color_t colorPCMEMCALHighPt         = GetDefaultColorDiffDetectors("PCM-EMCal", kFALSE, kFALSE, kTRUE);
    Color_t colorPCMPHOSLowPt           = GetDefaultColorDiffDetectors("PCM-PHOS", kFALSE, kFALSE, kFALSE);
    Color_t colorPCMPHOSHighPt          = GetDefaultColorDiffDetectors("PCM-PHOS", kFALSE, kFALSE, kTRUE);
    
    cout << "*************************************************************************"<< endl;  
    cout << "******************************  Pi0 pp **********************************"<< endl;
    cout << "*************************************************************************"<< endl;    
    
    TFile* fileNeutralPionCombDataPP                            = new TFile(nameFilePP.Data());
    TGraphAsymmErrors* graphInvYieldPi0Comb7TeV                 = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeV");
    graphInvYieldPi0Comb7TeV                                    = ScaleGraph(graphInvYieldPi0Comb7TeV,1./xSection7TeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0Comb7TeVStatErr          = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeVStatErr");
    graphInvYieldPi0Comb7TeVStatErr                             = ScaleGraph(graphInvYieldPi0Comb7TeVStatErr,1./xSection7TeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0Comb7TeVSysErr           = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeVSysErr");
    graphInvYieldPi0Comb7TeVSysErr                              = ScaleGraph(graphInvYieldPi0Comb7TeVSysErr,1./xSection7TeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PCM7TeVStatErr           = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMStat7TeV");
    graphInvYieldPi0PCM7TeVStatErr                              = ScaleGraph(graphInvYieldPi0PCM7TeVStatErr,1./xSection7TeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PCM7TeVSysErr            = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMSys7TeV");
    graphInvYieldPi0PCM7TeVSysErr                               = ScaleGraph(graphInvYieldPi0PCM7TeVSysErr,1./xSection7TeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PHOS7TeVStatErr          = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSStat7TeV");
    graphInvYieldPi0PHOS7TeVStatErr                             = ScaleGraph(graphInvYieldPi0PHOS7TeVStatErr,1./xSection7TeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PHOS7TeVSysErr           = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSSys7TeV");
    graphInvYieldPi0PHOS7TeVSysErr                              = ScaleGraph(graphInvYieldPi0PHOS7TeVSysErr,1./xSection7TeVppINEL);

    TGraphAsymmErrors* graphInvYieldPi0Comb2760GeVStatErr       = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeVStatErr");
    graphInvYieldPi0Comb2760GeVStatErr                          = ScaleGraph(graphInvYieldPi0Comb2760GeVStatErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0Comb2760GeVSysErr        = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeVSysErr");
    graphInvYieldPi0Comb2760GeVSysErr                           = ScaleGraph(graphInvYieldPi0Comb2760GeVSysErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PCM2760GeVStatErr        = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
    graphInvYieldPi0PCM2760GeVStatErr                           = ScaleGraph(graphInvYieldPi0PCM2760GeVStatErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PCM2760GeVSysErr         = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
    graphInvYieldPi0PCM2760GeVSysErr                            = ScaleGraph(graphInvYieldPi0PCM2760GeVSysErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PHOS2760GeVStatErr       = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOS2760GeVStatErr");
    graphInvYieldPi0PHOS2760GeVStatErr                          = ScaleGraph(graphInvYieldPi0PHOS2760GeVStatErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PHOS2760GeVSysErr        = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr");
    graphInvYieldPi0PHOS2760GeVSysErr                           = ScaleGraph(graphInvYieldPi0PHOS2760GeVSysErr,1./xSection2760GeVppINEL);

    
    TGraphAsymmErrors* graphInvYieldPi0Comb900GeV               = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeV");
    graphInvYieldPi0Comb900GeV                                  = ScaleGraph(graphInvYieldPi0Comb900GeV,1./xSection900GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0Comb900GeVStatErr        = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeVStatErr");
    graphInvYieldPi0Comb900GeVStatErr                           = ScaleGraph(graphInvYieldPi0Comb900GeVStatErr,1./xSection900GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0Comb900GeVSysErr         = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeVSysErr");
    graphInvYieldPi0Comb900GeVSysErr                            = ScaleGraph(graphInvYieldPi0Comb900GeVSysErr,1./xSection900GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PCM900GeVStatErr         = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMStat900GeV");
    graphInvYieldPi0PCM900GeVStatErr                            = ScaleGraph(graphInvYieldPi0PCM900GeVStatErr,1./xSection900GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PCM900GeVSysErr          = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMSys900GeV");
    graphInvYieldPi0PCM900GeVSysErr                             = ScaleGraph(graphInvYieldPi0PCM900GeVSysErr,1./xSection900GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PHOS900GeVStatErr        = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSStat900GeV");
    graphInvYieldPi0PHOS900GeVStatErr                           = ScaleGraph(graphInvYieldPi0PHOS900GeVStatErr,1./xSection900GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PHOS900GeVSysErr         = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSSys900GeV");
    graphInvYieldPi0PHOS900GeVSysErr                            = ScaleGraph(graphInvYieldPi0PHOS900GeVSysErr,1./xSection900GeVppINEL);

    TFile* fileNeutralPionCombDataPPUpdated                     = new TFile(fileNameComb2760GeVPP.Data());
    TGraphAsymmErrors* graphInvYieldPi0CombUp2760GeVStatErr     = (TGraphAsymmErrors*)fileNeutralPionCombDataPPUpdated->Get("graphInvCrossSectionPi0Comb2760GeVAStatErr");
    graphInvYieldPi0CombUp2760GeVStatErr                        = ScaleGraph(graphInvYieldPi0CombUp2760GeVStatErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0CombUp2760GeVSysErr      = (TGraphAsymmErrors*)fileNeutralPionCombDataPPUpdated->Get("graphInvCrossSectionPi0Comb2760GeVASysErr");
    graphInvYieldPi0CombUp2760GeVSysErr                         = ScaleGraph(graphInvYieldPi0CombUp2760GeVSysErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PCMEMCAL2760GeVStatErr   = (TGraphAsymmErrors*)fileNeutralPionCombDataPPUpdated->Get("graphInvCrossSectionPi0PCMEMCAL2760GeVStatErr");
    graphInvYieldPi0PCMEMCAL2760GeVStatErr                      = ScaleGraph(graphInvYieldPi0PCMEMCAL2760GeVStatErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0PCMEMCAL2760GeVSysErr    = (TGraphAsymmErrors*)fileNeutralPionCombDataPPUpdated->Get("graphInvCrossSectionPi0PCMEMCAL2760GeVSysErr");
    graphInvYieldPi0PCMEMCAL2760GeVSysErr                       = ScaleGraph(graphInvYieldPi0PCMEMCAL2760GeVSysErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0EMCAL2760GeVStatErr      = (TGraphAsymmErrors*)fileNeutralPionCombDataPPUpdated->Get("graphInvCrossSectionPi0EMCAL2760GeVStatErr");
    graphInvYieldPi0EMCAL2760GeVStatErr                         = ScaleGraph(graphInvYieldPi0EMCAL2760GeVStatErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0EMCAL2760GeVSysErr       = (TGraphAsymmErrors*)fileNeutralPionCombDataPPUpdated->Get("graphInvCrossSectionPi0EMCAL2760GeVSysErr");
    graphInvYieldPi0EMCAL2760GeVSysErr                          = ScaleGraph(graphInvYieldPi0EMCAL2760GeVSysErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0EMCALMerged2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPPUpdated->Get("graphInvCrossSectionPi0EMCALMerged2760GeVStatErr");
    graphInvYieldPi0EMCALMerged2760GeVStatErr                   = ScaleGraph(graphInvYieldPi0EMCALMerged2760GeVStatErr,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvYieldPi0EMCALMerged2760GeVSysErr = (TGraphAsymmErrors*)fileNeutralPionCombDataPPUpdated->Get("graphInvCrossSectionPi0EMCALMerged2760GeVSysErr");
    graphInvYieldPi0EMCALMerged2760GeVSysErr                    = ScaleGraph(graphInvYieldPi0EMCALMerged2760GeVSysErr,1./xSection2760GeVppINEL);
    
    
    cout << "*************************************************************************"<< endl;  
    cout << "******************************  charged pion PbPb ************************"<< endl;
    cout << "*************************************************************************"<< endl;    

    TFile* fileChargedPionInputPbPb             = new TFile(fileNameChargedPionPbPb.Data());
    TH1D* histoChargedPionSpecHighPtStat0005    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat0005");
    TH1D* histoChargedPionSpecHighPtSyst0005    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst0005");
    TH1D* histoChargedPionSpecHighPtStat0510    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat0510");
    TH1D* histoChargedPionSpecHighPtSyst0510    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst0510");
    TH1D* histoChargedPionSpecHighPtStat1020    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat1020");
    TH1D* histoChargedPionSpecHighPtSyst1020    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst1020");
    TH1D* histoChargedPionSpecHighPtSyst2040    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst2040");
    TH1D* histoChargedPionSpecHighPtStat2040    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat2040");
    TH1D* histoChargedPionSpecHighPtSyst4060    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst4060");
    TH1D* histoChargedPionSpecHighPtStat4060    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat4060");
    TH1D* histoChargedPionSpecHighPtSyst6080    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst6080");
    TH1D* histoChargedPionSpecHighPtStat6080    = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat6080");
    TH1D* histoChargedPionSpecLowPtStat0005     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat0005");
    TH1D* histoChargedPionSpecLowPtSyst0005     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst0005");
    TH1D* histoChargedPionSpecLowPtStat0510     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat0510");
    TH1D* histoChargedPionSpecLowPtSyst0510     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst0510");
    TH1D* histoChargedPionSpecLowPtStat1020     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat1020");
    TH1D* histoChargedPionSpecLowPtSyst1020     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst1020");
    TH1D* histoChargedPionSpecLowPtStat2040     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat2040");
    TH1D* histoChargedPionSpecLowPtSyst2040     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst2040");
    TH1D* histoChargedPionSpecLowPtStat4060     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat4060");
    TH1D* histoChargedPionSpecLowPtSyst4060     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst4060");
    TH1D* histoChargedPionSpecLowPtStat6080     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat6080");
    TH1D* histoChargedPionSpecLowPtSyst6080     = (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst6080");

    cout << "*************************************************************************"<< endl;  
    cout << "******************************  charged pion pp *************************"<< endl;
    cout << "*************************************************************************"<< endl;    

    TFile* fileChargedPionInputpp               = new TFile(fileNameChargedPionPP.Data());
    TH1D* histoChargedPionSpecHighPtStatPP      = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtStat2760GeV");
    TH1D* histoChargedPionSpecHighPtSystPP      = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtSyst2760GeV");
    
    TGraphAsymmErrors* graphChargedPionNegSpecXiangoStatPP2760GeV   = (TGraphAsymmErrors*)fileChargedPionInputpp->Get("graphChargedPionNegSpecXiangoStat2760GeV");
    TGraphAsymmErrors* graphChargedPionNegSpecXiangoSystPP2760GeV   = (TGraphAsymmErrors*)fileChargedPionInputpp->Get("graphChargedPionNegSpecXiangoSyst2760GeV");
   
    TH1D*    histoChargedPionSpecLowPtStat2760GeVCMS            = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat2760GeVCMS");
    TH1D*    histoChargedPionSpecLowPtSys2760GeVCMS             = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys2760GeVCMS");
    TH1D*    histoChargedPionSpecLowPtStatPP2760GeV             = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStatPP2760GeV");
    TH1D*    histoChargedPionSpecLowPtSysPP2760GeV              = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSysPP2760GeV");

    TH1D*    histoChargedPionSpecLowPtStat7TeVCMS               = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat7TeVCMS");
    TH1D*    histoChargedPionSpecLowPtSys7TeVCMS                = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys7TeVCMS");
    TH1D*    histoChargedPionSpecLowPtStatPP7TeV                = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat7TeVALICE");
    TH1D*    histoChargedPionSpecLowPtSysPP7TeV                 = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys7TeVALICE");
    TH1D*    histoChargedPionSpecHighPtStatPP7TeV               = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtStat7TeVALICE");
    TGraphAsymmErrors* graphChargedPionSpecHighPtSystPP7TeV     = (TGraphAsymmErrors*)fileChargedPionInputpp->Get("graphChargedPionSpecHighPtSys7TeVALICE");
    
    TH1D*    histoChargedPionSpecLowPtStat900GeVCMS     = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat900GeVCMS");
    TH1D*    histoChargedPionSpecLowPtSys900GeVCMS      = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys900GeVCMS");
    TH1D*    histoChargedPionSpecLowPtStatPP900GeV      = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat900GeVALICE");
    TH1D*    histoChargedPionSpecLowPtSysPP900GeV       = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys900GeVALICE");


    cout << "*************************************************************************"<< endl;  
    cout << "******************************  Pi0 EMCAL pp  ***************************"<< endl;
    cout << "*************************************************************************"<< endl;    

// Evi's file    
    TFile* fileEMCalPion7TeVPP         = new TFile(fileNameEMCalPion7TeVPP.Data());
    TH1D* histoEMCalPion7TeVPPStat     = (TH1D*)fileEMCalPion7TeVPP->Get("pi0Stat");
    TH1D* histoEMCalPion7TeVPPSyst     = (TH1D*)fileEMCalPion7TeVPP->Get("pi0Syst");
    histoEMCalPion7TeVPPStat->Scale(1./xSection7TeVppINEL);
    histoEMCalPion7TeVPPSyst->Scale(1./xSection7TeVppINEL);

    // Haitao's file
//     TFile* fileEMCalPion7TeVPP         = new TFile(fileNameEMCalPion7TeVPP.Data());
//     TH1D* histoEMCalPion7TeVPPStat     = (TH1D*)fileEMCalPion7TeVPP->Get("h_xsec_Stat");
//     TH1D* histoEMCalPion7TeVPPSyst     = (TH1D*)fileEMCalPion7TeVPP->Get("h_xsec_Syst");
//     histoEMCalPion7TeVPPStat->Scale(1./xSection7TeVppINEL*0.945*0.945);
//     histoEMCalPion7TeVPPSyst->Scale(1./xSection7TeVppINEL*0.945*0.945);
    
    
//     TFile* fileEMCALPion2760PP                         = new TFile(fileNameEMCalPion2760GeVPP.Data());
//     TH1D*    histoNeutralPionsEMCAL2760GeVStatErr     = (TH1D*)fileEMCALPion2760PP->Get("h_xsec_Stat");
//     histoNeutralPionsEMCAL2760GeVStatErr->Scale(1./xSection2760GeVppINEL);
//     TH1D*    histoNeutralPionsEMCAL2760GeVSystErr     = (TH1D*)fileEMCALPion2760PP->Get("h_xsec_Syst");
//     histoNeutralPionsEMCAL2760GeVSystErr->Scale(1./xSection2760GeVppINEL);
//     TDirectory* directoryEMCALPi02760PP                     = (TDirectory*)fileEMCALPion2760PP->Get("Pi02.76TeV"); 
//     TH1D* histoEMCALPi0InvXSectionStat                  = (TH1D*)directoryEMCALPi02760PP->Get("InvCrossSectionPi0");        
//     TGraphAsymmErrors* graphEMCALPi0InvXSectionSys      = (TGraphAsymmErrors*)directoryEMCALPi02760PP->Get("InvCrossSectionPi0Sys");
    
    
    
    cout << "*************************************************************************"<< endl;  
    cout << "******************************  Pi0 PbPb ********************************"<< endl;
    cout << "*************************************************************************"<< endl;    
    
    TFile* fCombResults= new TFile(nameFilePbPb.Data());
    TGraphAsymmErrors*    graphYieldCombStatPi02760GeV          = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPComb2760GeV_StatErr");
    TGraphAsymmErrors*    graphYieldCombSysPi02760GeV           = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPComb2760GeV_SysErr");
    TGraphAsymmErrors*    graphYieldPCMStatPi02760GeV           = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPCM2760GeV_StatErr");
    TGraphAsymmErrors*    graphYieldPCMSysPi02760GeV            = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPCM2760GeV_SysErr");
    TGraphAsymmErrors*    graphYieldPHOSStatPi02760GeV          = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPHOS2760GeV_StatErr");
    TGraphAsymmErrors*    graphYieldPHOSSysPi02760GeV           = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPHOS2760GeV_SysErr");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb0005StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_0005");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb0005SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_0005");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb0005StatErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_0005");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb0005SysErr        = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_0005");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb0005StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_0005");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb0005SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_0005");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb0510StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_0510");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb0510SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_0510");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb0510StatErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_0510");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb0510SysErr        = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_0510");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb0510StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_0510");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb0510SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_0510");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb1020StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_1020");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb1020SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_1020");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb1020StatErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_1020");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb1020SysErr        = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_1020");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb1020StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_1020");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb1020SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_1020");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb2040StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_2040");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb2040SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_2040");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb2040StatErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_2040");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb2040SysErr        = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_2040");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb2040StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_2040");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb2040SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_2040");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb4060StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_4060");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb4060SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_4060");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb4060StatErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_4060");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb4060SysErr        = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_4060");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb4060StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_4060");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb4060SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_4060");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb6080StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_6080");
    TGraphAsymmErrors*    graphYieldPi0CombPbPb6080SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_6080");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb6080StatErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_6080");
    TGraphAsymmErrors*    graphYieldPi0PCMPbPb6080SysErr        = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_6080");
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb6080StatErr      = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_6080");

    //fixed stat error
    graphYieldPi0PHOSPbPb6080StatErr->RemovePoint(graphYieldPi0PHOSPbPb6080StatErr->GetN()-1);
    graphYieldPi0PHOSPbPb6080StatErr->Print();
    TGraphAsymmErrors*    graphYieldPi0PHOSPbPb6080SysErr       = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_6080");
    graphYieldPi0PHOSPbPb6080SysErr->Print();

    cout << "*************************************************************************"<< endl;  
    cout << "******************************  Dalitz *****************************"<< endl;
    cout << "*************************************************************************"<< endl;
    
    TFile* fileDalitzPbPb                                = new TFile(fileNameDalitzPbPb.Data());
    TDirectory*   directoryPi0DalitzPbPb2040             = (TDirectory*)fileDalitzPbPb->Get("Pi0_PbPb_2.76TeV_20-40%"); 
    TH1D* histoYieldPi0DalitzPbPb2040                    = (TH1D*)directoryPi0DalitzPbPb2040->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphYieldSysPi0DalitzPbPb2040    = (TGraphAsymmErrors*)directoryPi0DalitzPbPb2040->Get("Pi0SystError");
    TDirectory*   directoryPi0DalitzPbPb4060             = (TDirectory*)fileDalitzPbPb->Get("Pi0_PbPb_2.76TeV_40-60%"); 
    TH1D* histoYieldPi0DalitzPbPb4060                    = (TH1D*)directoryPi0DalitzPbPb4060->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphYieldSysPi0DalitzPbPb4060    = (TGraphAsymmErrors*)directoryPi0DalitzPbPb4060->Get("Pi0SystError");
    TDirectory*   directoryPi0DalitzPbPb6080             = (TDirectory*)fileDalitzPbPb->Get("Pi0_PbPb_2.76TeV_60-80%"); 
    TH1D* histoYieldPi0DalitzPbPb6080                    = (TH1D*)directoryPi0DalitzPbPb6080->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphYieldSysPi0DalitzPbPb6080    = (TGraphAsymmErrors*)directoryPi0DalitzPbPb6080->Get("Pi0SystError");

    cout << "*************************************************************************"<< endl;    
    cout << "***************************** Dalitz PP  ********************************"<< endl;
    cout << "*************************************************************************"<< endl;
    
    TFile* fileDalitz                                           = new TFile(fileNameDalitz7TeV.Data());
    TFile* fileDalitz2760GeV                                    = new TFile(fileNameDalitz2760GeV.Data());

    TDirectory* directoryPi0Dalit2760GeV                        = (TDirectory*)fileDalitz2760GeV->Get("Pi0Dalitz2.76TeV"); 
    TH1D* histoInvCrossSectionPi0Dalitz2760GeV                  = (TH1D*)directoryPi0Dalit2760GeV->Get("InvCrossSectionPi0");
    histoInvCrossSectionPi0Dalitz2760GeV->Scale(1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvCrossSectionSysPi0Dalitz2760GeV  = (TGraphAsymmErrors*)directoryPi0Dalit2760GeV->Get("InvCrossSectionPi0Sys");
    graphInvCrossSectionSysPi0Dalitz2760GeV                     = ScaleGraph(graphInvCrossSectionSysPi0Dalitz2760GeV,1./xSection2760GeVppINEL);
    
    TDirectory* directoryPi0Dalit7TeV                           = (TDirectory*)fileDalitz->Get("Pi0Dalitz7TeV"); 
    TH1D* histoInvCrossSectionPi0Dalitz7TeV                     = (TH1D*)directoryPi0Dalit7TeV->Get("InvCrossSectionPi0");
    histoInvCrossSectionPi0Dalitz7TeV->Scale(1./xSection7TeVppINEL);
    TGraphAsymmErrors* graphInvCrossSectionSysPi0Dalitz7TeV     = (TGraphAsymmErrors*)directoryPi0Dalit7TeV->Get("InvCrossSectionPi0Sys");
    graphInvCrossSectionSysPi0Dalitz7TeV                        = ScaleGraph(graphInvCrossSectionSysPi0Dalitz7TeV,1./xSection7TeVppINEL);

    
//     TH1D* histoInvCrossSectionPi0PCMEMCAL2760GeV                 = NULL;
//     TGraphAsymmErrors* graphInvCrossSectionSysPi0PCMEMCAL2760GeV = NULL;
//     if (enablePCMEMCALComp2760GeV){
//         cout << "*************************************************************************"<< endl;    
//         cout << "***************************** PCM-EMCAL PP  ********************************"<< endl;
//         cout << "*************************************************************************"<< endl;
//         
//         TFile* filePCMEMCAL2760GeV                                     = new TFile(fileNamePCMEMCALpp2760GeVPP.Data());
// 
//         TDirectory* directoryPi0PCMEMCAL2760GeV                     = (TDirectory*)filePCMEMCAL2760GeV->Get("Pi02.76TeV"); 
//         histoInvCrossSectionPi0PCMEMCAL2760GeV                         = (TH1D*)directoryPi0PCMEMCAL2760GeV->Get("InvCrossSectionPi0");
//         histoInvCrossSectionPi0PCMEMCAL2760GeV->Scale(1./xSection2760GeVppINEL);
//         graphInvCrossSectionSysPi0PCMEMCAL2760GeV                     = (TGraphAsymmErrors*)directoryPi0PCMEMCAL2760GeV->Get("InvCrossSectionPi0Sys");
//         graphInvCrossSectionSysPi0PCMEMCAL2760GeV                     = ScaleGraph(graphInvCrossSectionSysPi0PCMEMCAL2760GeV,1./xSection2760GeVppINEL);
//     }    
    
    TH1D* histoInvCrossSectionPi0PCMPHOS2760GeV                 = NULL;
    TGraphAsymmErrors* graphInvCrossSectionSysPi0PCMPHOS2760GeV = NULL;
    if (enablePCMPHOSComp2760GeV){
        cout << "*************************************************************************"<< endl;    
        cout << "***************************** PCM-PHOS PP *******************************"<< endl;
        cout << "*************************************************************************"<< endl;
        
        TFile* filePCMPHOS2760GeV                               = new TFile(fileNamePCMPHOSpp2760GeVPP.Data());

        TDirectory* directoryPi0PCMPHOS2760GeV                  = (TDirectory*)filePCMPHOS2760GeV->Get("Pi02.76TeV"); 
        histoInvCrossSectionPi0PCMPHOS2760GeV                   = (TH1D*)directoryPi0PCMPHOS2760GeV->Get("InvCrossSectionPi0");
        histoInvCrossSectionPi0PCMPHOS2760GeV->Scale(1./xSection2760GeVppINEL);
        graphInvCrossSectionSysPi0PCMPHOS2760GeV                = (TGraphAsymmErrors*)directoryPi0PCMPHOS2760GeV->Get("InvCrossSectionPi0Sys");
        graphInvCrossSectionSysPi0PCMPHOS2760GeV                = ScaleGraph(graphInvCrossSectionSysPi0PCMPHOS2760GeV,1./xSection2760GeVppINEL);
    }    
    
    
    cout << "*************************************************************************"<< endl;    
    cout << "******************************  PP 2.76TeV ******************************"<< endl;
    cout << "*************************************************************************"<< endl;
    TGraphAsymmErrors* graphYieldCombStatPi02760GeVCopy         = (TGraphAsymmErrors*) graphInvYieldPi0Comb2760GeVStatErr->Clone("graphYieldCombStatPi02760GeVCopy");
    TGraphAsymmErrors* graphYieldCombSysPi02760GeVCopy          = (TGraphAsymmErrors*) graphInvYieldPi0Comb2760GeVSysErr->Clone("graphYieldCombSysPi02760GeVCopy");
    TGraphAsymmErrors* graphYieldPCMStatPi02760GeVCopy          = (TGraphAsymmErrors*) graphInvYieldPi0PCM2760GeVStatErr->Clone("graphYieldPCMStatPi02760GeVCopy");
    TGraphAsymmErrors* graphYieldPCMSysPi02760GeVCopy           = (TGraphAsymmErrors*) graphInvYieldPi0PCM2760GeVSysErr->Clone("graphYieldPCMSysPi02760GeVCopy");
    TGraphAsymmErrors* graphYieldPHOSStatPi02760GeVCopy         = (TGraphAsymmErrors*) graphInvYieldPi0PHOS2760GeVStatErr->Clone("graphYieldPHOSStatPi02760GeVCopy");
    TGraphAsymmErrors* graphYieldPHOSSysPi02760GeVCopy          = (TGraphAsymmErrors*) graphInvYieldPi0PHOS2760GeVSysErr->Clone("graphYieldPHOSSysPi02760GeVCopy");
    
       cout << "combined Spectrum - high Pt" << endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPPHighPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPPHighPtComb    = NULL;
    TGraphErrors* graphYieldCombStatPi02760GeVRebinnedHighPtComb= NULL;
    TGraphErrors* graphYieldCombSysPi02760GeVRebinnedHighPtComb = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsCombPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombStatPi02760GeVCopy, graphYieldCombSysPi02760GeVCopy, 
                                                                                                        histoChargedPionSpecHighPtStatPP, histoChargedPionSpecHighPtSystPP,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphYieldCombStatPi02760GeVRebinnedHighPtComb, &graphYieldCombSysPi02760GeVRebinnedHighPtComb, 
                                                                                                        &graphChargedPionSpecHighPtStatPPHighPtComb, &graphChargedPionSpecHighPtSystPPHighPtComb )    ;
    graphRatioHighPtChargedPionsCombPP->Print();
   
    cout << "combined Spectrum - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPPLowPtComb      = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPPLowPtComb      = NULL;
    TGraphErrors* graphYieldCombStatPi02760GeVRebinnedLowPtComb = NULL;
    TGraphErrors* graphYieldCombSysPi02760GeVRebinnedLowPtComb  = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsCombPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombStatPi02760GeVCopy, graphYieldCombSysPi02760GeVCopy, 
                                                                                                       histoChargedPionSpecLowPtStatPP2760GeV, histoChargedPionSpecLowPtSysPP2760GeV,  
                                                                                                       kTRUE,  kTRUE, 
                                                                                                       &graphYieldCombStatPi02760GeVRebinnedLowPtComb, &graphYieldCombSysPi02760GeVRebinnedLowPtComb,
                                                                                                       &graphChargedPionSpecLowPtStatPPLowPtComb, &graphChargedPionSpecLowPtSystPPLowPtComb )    ;
    graphRatioLowPtChargedPionsCombPP->Print();
    cout << "combined Spectrum - low Pt CMS" << endl;
   
    TGraphErrors* graphChargedPionSpecCMSStatPPCMSComb          = NULL;
    TGraphErrors* graphChargedPionSpecCMSSystPPCMSComb          = NULL;
    TGraphErrors* graphYieldCombStatPi02760GeVRebinnedCMSComb   = NULL;
    TGraphErrors* graphYieldCombSysPi02760GeVRebinnedCMSComb    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsCombPPCMS = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombStatPi02760GeVCopy, graphYieldCombSysPi02760GeVCopy, 
                                                                                                          histoChargedPionSpecLowPtStat2760GeVCMS, histoChargedPionSpecLowPtSys2760GeVCMS,  
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldCombStatPi02760GeVRebinnedCMSComb, &graphYieldCombSysPi02760GeVRebinnedCMSComb,
                                                                                                          &graphChargedPionSpecCMSStatPPCMSComb, &graphChargedPionSpecCMSSystPPCMSComb ) ;
    graphRatioLowPtChargedPionsCombPPCMS->Print();
       
    cout << "PCM - high Pt" << endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPPHighPtPCM     = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPPHighPtPCM     = NULL;
    TGraphErrors* graphYieldPCMStatPi02760GeVRebinnedHighPtPCM  = NULL;
    TGraphErrors* graphYieldPCMSysPi02760GeVRebinnedHighPtPCM   = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCMPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPCMStatPi02760GeVCopy, graphYieldPCMSysPi02760GeVCopy, 
                                                                                                       histoChargedPionSpecHighPtStatPP, histoChargedPionSpecHighPtSystPP,  
                                                                                                       kTRUE,  kTRUE, 
                                                                                                       &graphYieldPCMStatPi02760GeVRebinnedHighPtPCM, &graphYieldPCMSysPi02760GeVRebinnedHighPtPCM, 
                                                                                                       &graphChargedPionSpecHighPtStatPPHighPtPCM, &graphChargedPionSpecHighPtSystPPHighPtPCM ) ;
    graphRatioHighPtChargedPionsPCMPP->Print();
    
    cout << "PCM - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPPLowPtPCM       = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPPLowPtPCM       = NULL;
    TGraphErrors* graphYieldPCMStatPi02760GeVRebinnedLowPtPCM   = NULL;
    TGraphErrors* graphYieldPCMSysPi02760GeVRebinnedLowPtPCM    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCMPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPCMStatPi02760GeVCopy, graphYieldPCMSysPi02760GeVCopy, 
                                                                                                      histoChargedPionSpecLowPtStatPP2760GeV, histoChargedPionSpecLowPtSysPP2760GeV, 
                                                                                                      kTRUE,  kTRUE, 
                                                                                                      &graphYieldPCMStatPi02760GeVRebinnedLowPtPCM, &graphYieldPCMSysPi02760GeVRebinnedLowPtPCM, 
                                                                                                      &graphChargedPionSpecLowPtStatPPLowPtPCM, &graphChargedPionSpecLowPtSystPPLowPtPCM ) ;
    graphRatioLowPtChargedPionsPCMPP->Print();
    
    cout << "PCM - low Pt CMS" << endl;
    TGraphErrors* graphChargedPionSpecCMSStatPPCMSPCM           = NULL;
    TGraphErrors* graphChargedPionSpecCMSSystPPCMSPCM           = NULL;
    TGraphErrors* graphYieldPCMStatPi02760GeVRebinnedCMSPCM     = NULL;
    TGraphErrors* graphYieldPCMSysPi02760GeVRebinnedCMSPCM      = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCMPPCMS = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPCMStatPi02760GeVCopy, graphYieldPCMSysPi02760GeVCopy,
                                                                                                         histoChargedPionSpecLowPtStat2760GeVCMS, histoChargedPionSpecLowPtSys2760GeVCMS, 
                                                                                                         kTRUE,  kTRUE,
                                                                                                         &graphYieldPCMStatPi02760GeVRebinnedCMSPCM, &graphYieldPCMSysPi02760GeVRebinnedCMSPCM,
                                                                                                         &graphChargedPionSpecCMSStatPPCMSPCM, &graphChargedPionSpecCMSSystPPCMSPCM ) ;
    graphRatioLowPtChargedPionsPCMPPCMS->Print();
    
    cout << "PHOS - high Pt" << endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPPHighPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPPHighPtPHOS    = NULL;
    TGraphErrors* graphYieldPHOSStatPi02760GeVRebinnedHighPtPHOS= NULL;
    TGraphErrors* graphYieldPHOSSysPi02760GeVRebinnedHighPtPHOS = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPHOSPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPHOSStatPi02760GeVCopy, graphYieldPHOSSysPi02760GeVCopy,
                                                                                                        histoChargedPionSpecHighPtStatPP, histoChargedPionSpecHighPtSystPP,  
                                                                                                        kTRUE,  kTRUE,
                                                                                                        &graphYieldPHOSStatPi02760GeVRebinnedHighPtPHOS, &graphYieldPHOSSysPi02760GeVRebinnedHighPtPHOS, 
                                                                                                        &graphChargedPionSpecHighPtStatPPHighPtPHOS, &graphChargedPionSpecHighPtSystPPHighPtPHOS )    ;
    graphRatioHighPtChargedPionsPHOSPP->Print();
    
    cout << "PHOS - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPPLowPtPHOS      = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPPLowPtPHOS      = NULL;
    TGraphErrors* graphYieldPHOSStatPi02760GeVRebinnedLowPtPHOS = NULL;
    TGraphErrors* graphYieldPHOSSysPi02760GeVRebinnedLowPtPHOS  = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPHOSPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPHOSStatPi02760GeVCopy, graphYieldPHOSSysPi02760GeVCopy, 
                                                                                                       histoChargedPionSpecLowPtStatPP2760GeV, histoChargedPionSpecLowPtSysPP2760GeV, 
                                                                                                       kTRUE,  kTRUE, 
                                                                                                       &graphYieldPHOSStatPi02760GeVRebinnedLowPtPHOS, &graphYieldPHOSSysPi02760GeVRebinnedLowPtPHOS, 
                                                                                                       &graphChargedPionSpecLowPtStatPPLowPtPHOS, &graphChargedPionSpecLowPtSystPPLowPtPHOS );
    graphRatioLowPtChargedPionsPHOSPP->Print();


    //***************************** ratios Dalitz 2.76 TeV ************************************************
    TH1D* histoInvCrossSectionPi0Dalitz2760GeVCopy                          = (TH1D*) histoInvCrossSectionPi0Dalitz2760GeV->Clone("histoInvCrossSectionPi0Dalitz2760GeVCopy");
    TGraphAsymmErrors* graphInvCrossSectionSysPi0Dalitz2760GeVCopy          = (TGraphAsymmErrors*) graphInvCrossSectionSysPi0Dalitz2760GeV->Clone("graphInvCrossSectionSysPi0Dalitz2760GeVCopy");
    cout << "Dalitz Spectrum - high Pt" << endl;
    TGraphErrors* graphYieldDalitzStatPi02760GeVRebinnedHighPtDalitz        = NULL;
    TGraphErrors* graphYieldDalitzSysPi02760GeVRebinnedHighPtDalitz         = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat2760GeVHighPtDalitz         = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst2760GeVHighPtDalitz         = NULL;   
    TGraphErrors* graphRatioHighPtChargedPionsDalitzPP2760GeV = CalculateRatioBetweenSpectraWithDifferentBinning(histoInvCrossSectionPi0Dalitz2760GeVCopy, graphInvCrossSectionSysPi0Dalitz2760GeVCopy, 
                                                                                                                 histoChargedPionSpecHighPtStatPP, histoChargedPionSpecHighPtSystPP,  
                                                                                                                 kTRUE,  kTRUE, 
                                                                                                                 &graphYieldDalitzStatPi02760GeVRebinnedHighPtDalitz, &graphYieldDalitzSysPi02760GeVRebinnedHighPtDalitz, 
                                                                                                                 &graphChargedPionSpecHighPtStat2760GeVHighPtDalitz, &graphChargedPionSpecHighPtSyst2760GeVHighPtDalitz);
    graphRatioHighPtChargedPionsDalitzPP2760GeV->Print();
    
    cout << "Dalitz Spectrum - low Pt" << endl;
    TGraphErrors* graphYieldDalitzStatPi02760GeVRebinnedLowPtDalitz         = NULL;
    TGraphErrors* graphYieldDalitzSysPi02760GeVRebinnedLowPtDalitz          = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat2760GeVLowPtDalitz           = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst2760GeVLowPtDalitz           = NULL;   
    TGraphErrors* graphRatioLowPtChargedPionsDalitzPP2760GeV = CalculateRatioBetweenSpectraWithDifferentBinning(histoInvCrossSectionPi0Dalitz2760GeVCopy, graphInvCrossSectionSysPi0Dalitz2760GeVCopy, 
                                                                                                                histoChargedPionSpecLowPtStatPP2760GeV, histoChargedPionSpecLowPtSysPP2760GeV, 
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphYieldDalitzStatPi02760GeVRebinnedLowPtDalitz, &graphYieldDalitzSysPi02760GeVRebinnedLowPtDalitz,
                                                                                                                &graphChargedPionSpecLowPtStat2760GeVLowPtDalitz, &graphChargedPionSpecLowPtSyst2760GeVLowPtDalitz);
    graphRatioLowPtChargedPionsDalitzPP2760GeV->Print();

    //***************************** ratios EMCAL 2.76 TeV ************************************************    
    cout << "EMCAL to low pT charged pions" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPP2760GeVLowPtEMCAL          = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPP2760GeVLowPtEMCAL          = NULL;
    TGraphErrors* graphYieldEMCALStatPi02760GeVRebinnedLowPtEMCAL           = NULL;
    TGraphErrors* graphYieldEMCALSysPi02760GeVRebinnedLowPtEMCAL            = NULL;
     TGraphErrors* graphRatioLowPtChargedPionsEMCALPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphInvYieldPi0EMCAL2760GeVStatErr, graphInvYieldPi0EMCAL2760GeVSysErr,   
                                                                                                        histoChargedPionSpecLowPtStatPP2760GeV, histoChargedPionSpecLowPtSysPP2760GeV, 
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphYieldEMCALStatPi02760GeVRebinnedLowPtEMCAL, &graphYieldEMCALSysPi02760GeVRebinnedLowPtEMCAL,
                                                                                                        &graphChargedPionSpecLowPtStatPP2760GeVLowPtEMCAL, &graphChargedPionSpecLowPtSystPP2760GeVLowPtEMCAL );
    graphRatioLowPtChargedPionsEMCALPP->Print();    

    cout << endl<< endl<< endl<< endl<< "*****************************************************************************************************" << endl;
    cout << "EMCAL to high pT charged pions" << endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPP2760GeVHighPtEMCAL        = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPP2760GeVHighPtEMCAL        = NULL;
    TGraphErrors* graphYieldEMCALStatPi02760GeVRebinnedHighPtEMCAL          = NULL;
    TGraphErrors* graphYieldEMCALSysPi02760GeVRebinnedHighPtEMCAL           = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsEMCALPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphInvYieldPi0EMCAL2760GeVStatErr, graphInvYieldPi0EMCAL2760GeVSysErr,
                                                                                                         histoChargedPionSpecHighPtStatPP, histoChargedPionSpecHighPtSystPP,  
                                                                                                         kTRUE,  kTRUE,
                                                                                                         &graphYieldEMCALStatPi02760GeVRebinnedHighPtEMCAL, &graphYieldEMCALSysPi02760GeVRebinnedHighPtEMCAL,
                                                                                                         &graphChargedPionSpecHighPtStatPP2760GeVHighPtEMCAL, &graphChargedPionSpecHighPtSystPP2760GeVHighPtEMCAL );
    graphRatioHighPtChargedPionsEMCALPP->Print();

    //***************************** ratios EMCALmerged 2.76 TeV ************************************************    
    cout << endl<< endl<< endl<< endl<< "*****************************************************************************************************" << endl;
    cout << "EMCALmerged to high pT charged pions" << endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPP2760GeVHighPtEMCALmerged      = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPP2760GeVHighPtEMCALmerged      = NULL;
    TGraphErrors* graphYieldEMCALmergedStatPi02760GeVRebinnedHighPtEMCALmerged  = NULL;
    TGraphErrors* graphYieldEMCALmergedSysPi02760GeVRebinnedHighPtEMCALmerged   = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsEMCALMergedPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphInvYieldPi0EMCALMerged2760GeVStatErr, graphInvYieldPi0EMCALMerged2760GeVSysErr,
                                                                                                         histoChargedPionSpecHighPtStatPP, histoChargedPionSpecHighPtSystPP,  
                                                                                                         kTRUE,  kTRUE,
                                                                                                         &graphYieldEMCALmergedStatPi02760GeVRebinnedHighPtEMCALmerged, &graphYieldEMCALmergedSysPi02760GeVRebinnedHighPtEMCALmerged,
                                                                                                         &graphChargedPionSpecHighPtStatPP2760GeVHighPtEMCALmerged, &graphChargedPionSpecHighPtSystPP2760GeVHighPtEMCALmerged );
    graphRatioHighPtChargedPionsEMCALMergedPP->Print();
    
    
    //***************************** ratios PCM-EMCAL 2.76 TeV *********************************************
    TGraphErrors* graphChargedPionSpecHighPtStatPPHighPtPCMEMCAL            = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPPHighPtPCMEMCAL            = NULL;
    TGraphErrors* graphYieldPCMEMCALStatPi02760GeVRebinnedHighPtPCMEMCAL    = NULL;
    TGraphErrors* graphYieldPCMEMCALSysPi02760GeVRebinnedHighPtPCMEMCAL     = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCMEMCALPP                    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStatPPLowPtPCMEMCAL              = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPPLowPtPCMEMCAL              = NULL;
    TGraphErrors* graphYieldPCMEMCALStatPi02760GeVRebinnedLowPtPCMEMCAL     = NULL;
    TGraphErrors* graphYieldPCMEMCALSysPi02760GeVRebinnedLowPtPCMEMCAL      = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCMEMCALPP                     = NULL;

//     if (enablePCMEMCALComp2760GeV){
        TGraphAsymmErrors* graphInvYieldPi0PCMEMCAL2760GeVStatErrCopy       = (TGraphAsymmErrors*) graphInvYieldPi0PCMEMCAL2760GeVStatErr->Clone("graphInvYieldPi0PCMEMCAL2760GeVStatErrCopy");
        TGraphAsymmErrors* graphInvYieldPi0PCMEMCAL2760GeVSysErrCopy        = (TGraphAsymmErrors*) graphInvYieldPi0PCMEMCAL2760GeVSysErr->Clone("graphInvYieldPi0PCMEMCAL2760GeVSysErrCopy");
//         graphInvYieldPi0PCMEMCAL2760GeVSysErrCopy->RemovePoint(0);
        
        cout << "PCMEMCAL - high Pt" << endl;
        graphRatioHighPtChargedPionsPCMEMCALPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphInvYieldPi0PCMEMCAL2760GeVStatErrCopy, graphInvYieldPi0PCMEMCAL2760GeVSysErrCopy, 
                                                                                              histoChargedPionSpecHighPtStatPP, histoChargedPionSpecHighPtSystPP,  
                                                                                              kTRUE,  kTRUE, 
                                                                                              &graphYieldPCMEMCALStatPi02760GeVRebinnedHighPtPCMEMCAL, &graphYieldPCMEMCALSysPi02760GeVRebinnedHighPtPCMEMCAL, 
                                                                                              &graphChargedPionSpecHighPtStatPPHighPtPCMEMCAL, &graphChargedPionSpecHighPtSystPPHighPtPCMEMCAL ) ;
        graphRatioHighPtChargedPionsPCMEMCALPP->Print();
    
        cout << "PCMEMCAL - low Pt" << endl;
        graphRatioLowPtChargedPionsPCMEMCALPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphInvYieldPi0PCMEMCAL2760GeVStatErrCopy, graphInvYieldPi0PCMEMCAL2760GeVSysErrCopy, 
                                                                                                histoChargedPionSpecLowPtStatPP2760GeV, histoChargedPionSpecLowPtSysPP2760GeV, 
                                                                                                kTRUE,  kTRUE, 
                                                                                                &graphYieldPCMEMCALStatPi02760GeVRebinnedLowPtPCMEMCAL, &graphYieldPCMEMCALSysPi02760GeVRebinnedLowPtPCMEMCAL, 
                                                                                                &graphChargedPionSpecLowPtStatPPLowPtPCMEMCAL, &graphChargedPionSpecLowPtSystPPLowPtPCMEMCAL ) ;
        graphRatioLowPtChargedPionsPCMPP->Print();
//     }
    

    TGraphErrors* graphChargedPionSpecHighPtStatPPHighPtPCMPHOS             = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPPHighPtPCMPHOS             = NULL;
    TGraphErrors* graphYieldPCMPHOSStatPi02760GeVRebinnedHighPtPCMPHOS      = NULL;
    TGraphErrors* graphYieldPCMPHOSSysPi02760GeVRebinnedHighPtPCMPHOS       = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCMPHOSPP                     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStatPPLowPtPCMPHOS               = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPPLowPtPCMPHOS               = NULL;
    TGraphErrors* graphYieldPCMPHOSStatPi02760GeVRebinnedLowPtPCMPHOS       = NULL;
    TGraphErrors* graphYieldPCMPHOSSysPi02760GeVRebinnedLowPtPCMPHOS        = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCMPHOSPP                      = NULL;

    if (enablePCMPHOSComp2760GeV){
        TH1D* histoInvCrossSectionPi0PCMPHOS2760GeVCopy                     = (TH1D*) histoInvCrossSectionPi0PCMPHOS2760GeV->Clone("histoInvCrossSectionPi0PCMPHOS2760GeVCopy");
        TGraphAsymmErrors* graphInvCrossSectionSysPi0PCMPHOS2760GeVCopy     = (TGraphAsymmErrors*) graphInvCrossSectionSysPi0PCMPHOS2760GeV->Clone("graphInvCrossSectionSysPi0PCMPHOS2760GeVCopy");

        cout << "PCMPHOS - high Pt" << endl;
        graphRatioHighPtChargedPionsPCMPHOSPP = CalculateRatioBetweenSpectraWithDifferentBinning(histoInvCrossSectionPi0PCMPHOS2760GeVCopy, graphInvCrossSectionSysPi0PCMPHOS2760GeVCopy, 
                                                                                              histoChargedPionSpecHighPtStatPP, histoChargedPionSpecHighPtSystPP,  
                                                                                              kTRUE,  kTRUE, 
                                                                                              &graphYieldPCMPHOSStatPi02760GeVRebinnedHighPtPCMPHOS, &graphYieldPCMPHOSSysPi02760GeVRebinnedHighPtPCMPHOS, 
                                                                                              &graphChargedPionSpecHighPtStatPPHighPtPCMPHOS, &graphChargedPionSpecHighPtSystPPHighPtPCMPHOS ) ;
        graphRatioHighPtChargedPionsPCMPHOSPP->Print();
    
        cout << "PCMPHOS - low Pt" << endl;
        graphRatioLowPtChargedPionsPCMPHOSPP = CalculateRatioBetweenSpectraWithDifferentBinning(histoInvCrossSectionPi0PCMPHOS2760GeVCopy, graphInvCrossSectionSysPi0PCMPHOS2760GeVCopy, 
                                                                                                histoChargedPionSpecLowPtStatPP2760GeV, histoChargedPionSpecLowPtSysPP2760GeV, 
                                                                                                kTRUE,  kTRUE, 
                                                                                                &graphYieldPCMPHOSStatPi02760GeVRebinnedLowPtPCMPHOS, &graphYieldPCMPHOSSysPi02760GeVRebinnedLowPtPCMPHOS, 
                                                                                                &graphChargedPionSpecLowPtStatPPLowPtPCMPHOS, &graphChargedPionSpecLowPtSystPPLowPtPCMPHOS ) ;
        graphRatioLowPtChargedPionsPCMPHOSPP->Print();
    }
    
    //*************************************************************************************************************
    //***************************** updated combined spectrum *****************************************************
    //*************************************************************************************************************
    TGraphAsymmErrors* graphYieldCombUpStatPi02760GeVCopy                   = (TGraphAsymmErrors*) graphInvYieldPi0CombUp2760GeVStatErr->Clone("graphYieldCombUpStatPi02760GeVCopy");
    TGraphAsymmErrors* graphYieldCombUpSysPi02760GeVCopy                    = (TGraphAsymmErrors*) graphInvYieldPi0CombUp2760GeVSysErr->Clone("graphYieldCombUpSysPi02760GeVCopy");

       cout << "combined Spectrum - high Pt" << endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPPHighPtCombUp              = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPPHighPtCombUp              = NULL;
    TGraphErrors* graphYieldCombStatPi02760GeVRebinnedHighPtCombUp          = NULL;
    TGraphErrors* graphYieldCombSysPi02760GeVRebinnedHighPtCombUp           = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsCombUpPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombUpStatPi02760GeVCopy, graphYieldCombUpSysPi02760GeVCopy, 
                                                                                                        histoChargedPionSpecHighPtStatPP, histoChargedPionSpecHighPtSystPP,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphYieldCombStatPi02760GeVRebinnedHighPtCombUp, &graphYieldCombSysPi02760GeVRebinnedHighPtCombUp, 
                                                                                                        &graphChargedPionSpecHighPtStatPPHighPtCombUp, &graphChargedPionSpecHighPtSystPPHighPtCombUp )    ;
   
    cout << "combined Spectrum - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPPLowPtCombUp                = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPPLowPtCombUp                = NULL;
    TGraphErrors* graphYieldCombStatPi02760GeVRebinnedLowPtCombUp           = NULL;
    TGraphErrors* graphYieldCombSysPi02760GeVRebinnedLowPtCombUp            = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsCombUpPP = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombUpStatPi02760GeVCopy, graphYieldCombUpSysPi02760GeVCopy, 
                                                                                                       histoChargedPionSpecLowPtStatPP2760GeV, histoChargedPionSpecLowPtSysPP2760GeV,  
                                                                                                       kTRUE,  kTRUE, 
                                                                                                       &graphYieldCombStatPi02760GeVRebinnedLowPtCombUp, &graphYieldCombSysPi02760GeVRebinnedLowPtCombUp,
                                                                                                       &graphChargedPionSpecLowPtStatPPLowPtCombUp, &graphChargedPionSpecLowPtSystPPLowPtCombUp )    ;
    
    cout << "combined Spectrum - low Pt CMS" << endl;
   
    TGraphErrors* graphChargedPionSpecCMSStatPPCMSCombUp                    = NULL;
    TGraphErrors* graphChargedPionSpecCMSSystPPCMSCombUp                    = NULL;
    TGraphErrors* graphYieldCombStatPi02760GeVRebinnedCMSCombUp             = NULL;
    TGraphErrors* graphYieldCombSysPi02760GeVRebinnedCMSCombUp              = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsCombUpPPCMS = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombUpStatPi02760GeVCopy, graphYieldCombUpSysPi02760GeVCopy, 
                                                                                                          histoChargedPionSpecLowPtStat2760GeVCMS, histoChargedPionSpecLowPtSys2760GeVCMS,  
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldCombStatPi02760GeVRebinnedCMSCombUp, &graphYieldCombSysPi02760GeVRebinnedCMSCombUp,
                                                                                                          &graphChargedPionSpecCMSStatPPCMSCombUp, &graphChargedPionSpecCMSSystPPCMSCombUp ) ;
    
    cout << "*************************************************************************"<< endl;    
    cout << "******************************  PP 7TeV *********************************"<< endl;
    cout << "*************************************************************************"<< endl;
    TGraphAsymmErrors* graphYieldCombStatPi07TeVCopy                = (TGraphAsymmErrors*) graphInvYieldPi0Comb7TeVStatErr->Clone("graphYieldCombStatPi07TeVCopy");
    TGraphAsymmErrors* graphYieldCombSysPi07TeVCopy                 = (TGraphAsymmErrors*) graphInvYieldPi0Comb7TeVSysErr->Clone("graphYieldCombSysPi07TeVCopy");
    TGraphAsymmErrors* graphYieldPCMStatPi07TeVCopy                 = (TGraphAsymmErrors*) graphInvYieldPi0PCM7TeVStatErr->Clone("graphYieldPCMStatPi07TeVCopy");
    TGraphAsymmErrors* graphYieldPCMSysPi07TeVCopy                  = (TGraphAsymmErrors*) graphInvYieldPi0PCM7TeVSysErr->Clone("graphYieldPCMSysPi07TeVCopy");
    TGraphAsymmErrors* graphYieldPHOSStatPi07TeVCopy                = (TGraphAsymmErrors*) graphInvYieldPi0PHOS7TeVStatErr->Clone("graphYieldPHOSStatPi07TeVCopy");
    TGraphAsymmErrors* graphYieldPHOSSysPi07TeVCopy                 = (TGraphAsymmErrors*) graphInvYieldPi0PHOS7TeVSysErr->Clone("graphYieldPHOSSysPi07TeVCopy");
    TGraphAsymmErrors* graphChargedPionSpecHighPtSystPP7TeVCopy     = (TGraphAsymmErrors*) graphChargedPionSpecHighPtSystPP7TeV->Clone("graphChargedPionSpecHighPtSystPP7TeVCopy");
    
    cout << "combined Spectrum - high Pt" << endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPP7TeVHighPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPP7TeVHighPtComb    = NULL;
    TGraphErrors* graphYieldCombStatPi07TeVRebinnedHighPtComb       = NULL;
    TGraphErrors* graphYieldCombSysPi07TeVRebinnedHighPtComb        = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsCombPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombStatPi07TeVCopy, graphYieldCombSysPi07TeVCopy, 
                                                                                                            histoChargedPionSpecHighPtStatPP7TeV, graphChargedPionSpecHighPtSystPP7TeV,  
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphYieldCombStatPi07TeVRebinnedHighPtComb, &graphYieldCombSysPi07TeVRebinnedHighPtComb, 
                                                                                                            &graphChargedPionSpecHighPtStatPP7TeVHighPtComb, &graphChargedPionSpecHighPtSystPP7TeVHighPtComb)    ;
    graphRatioHighPtChargedPionsCombPP7TeV->Print();
    
    cout << "combined Spectrum - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPP7TeVLowPtComb      = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPP7TeVLowPtComb      = NULL;
    TGraphErrors* graphYieldCombStatPi07TeVRebinnedLowPtComb        = NULL;
    TGraphErrors* graphYieldCombSysPi07TeVRebinnedLowPtComb         = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsCombPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombStatPi07TeVCopy, graphYieldCombSysPi07TeVCopy, 
                                                                                                           histoChargedPionSpecLowPtStatPP7TeV, histoChargedPionSpecLowPtSysPP7TeV,  
                                                                                                           kTRUE,  kTRUE, 
                                                                                                           &graphYieldCombStatPi07TeVRebinnedLowPtComb, &graphYieldCombSysPi07TeVRebinnedLowPtComb, 
                                                                                                           &graphChargedPionSpecLowPtStatPP7TeVLowPtComb, &graphChargedPionSpecLowPtSystPP7TeVLowPtComb);
    graphRatioLowPtChargedPionsCombPP7TeV->Print();

    cout << "combined Spectrum - low Pt CMS" << endl;
    TGraphErrors* graphChargedPionSpecCMSStatPP7TeVCMSComb          = NULL;
    TGraphErrors* graphChargedPionSpecCMSSystPP7TeVCMSComb          = NULL;
    TGraphErrors* graphYieldCombStatPi07TeVRebinnedCMSComb          = NULL;
    TGraphErrors* graphYieldCombSysPi07TeVRebinnedCMSComb           = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsCombPPCMS7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombStatPi07TeVCopy, graphYieldCombSysPi07TeVCopy, 
                                                                                                              histoChargedPionSpecLowPtStat7TeVCMS, histoChargedPionSpecLowPtSys7TeVCMS,  
                                                                                                              kTRUE,  kTRUE, 
                                                                                                              &graphYieldCombStatPi07TeVRebinnedCMSComb, &graphYieldCombSysPi07TeVRebinnedCMSComb, 
                                                                                                              &graphChargedPionSpecCMSStatPP7TeVCMSComb, &graphChargedPionSpecCMSSystPP7TeVCMSComb);
    graphRatioLowPtChargedPionsCombPPCMS7TeV->Print();
        
    cout << "PCM - high Pt" << endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPP7TeVHighPtPCM     = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPP7TeVHighPtPCM     = NULL;
    TGraphErrors* graphYieldPCMStatPi07TeVRebinnedHighPtPCM         = NULL;
    TGraphErrors* graphYieldPCMSysPi07TeVRebinnedHighPtPCM          = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCMPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPCMStatPi07TeVCopy, graphYieldPCMSysPi07TeVCopy,
                                                                                                           histoChargedPionSpecHighPtStatPP7TeV, graphChargedPionSpecHighPtSystPP7TeV, 
                                                                                                           kTRUE,  kTRUE, 
                                                                                                           &graphYieldPCMStatPi07TeVRebinnedHighPtPCM, &graphChargedPionSpecHighPtSystPP7TeVHighPtPCM,
                                                                                                           &graphChargedPionSpecHighPtStatPP7TeVHighPtPCM, &graphChargedPionSpecHighPtSystPP7TeVHighPtPCM);
    graphRatioHighPtChargedPionsPCMPP7TeV->Print();
    
    cout << "PCM - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPP7TeVLowPtPCM       = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPP7TeVLowPtPCM       = NULL;
    TGraphErrors* graphYieldPCMStatPi07TeVRebinnedLowPtPCM          = NULL;
    TGraphErrors* graphYieldPCMSysPi07TeVRebinnedLowPtPCM           = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCMPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPCMStatPi07TeVCopy, graphYieldPCMSysPi07TeVCopy,
                                                                                                          histoChargedPionSpecLowPtStatPP7TeV, histoChargedPionSpecLowPtSysPP7TeV,  
                                                                                                          kTRUE,  kTRUE,
                                                                                                          &graphYieldPCMStatPi07TeVRebinnedLowPtPCM, &graphYieldPCMSysPi07TeVRebinnedLowPtPCM, 
                                                                                                          &graphChargedPionSpecLowPtStatPP7TeVLowPtPCM, &graphChargedPionSpecLowPtSystPP7TeVLowPtPCM);
    graphRatioLowPtChargedPionsPCMPP7TeV->Print();
    
    cout << "PCM - low Pt CMS" << endl;
    TGraphErrors* graphChargedPionSpecCMSStatPP7TeVCMSPCM           = NULL;
    TGraphErrors* graphChargedPionSpecCMSSystPP7TeVCMSPCM           = NULL;
    TGraphErrors* graphYieldPCMStatPi07TeVRebinnedCMSPCM            = NULL;
    TGraphErrors* graphYieldPCMSysPi07TeVRebinnedCMSPCM             = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCMPPCMS7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPCMStatPi07TeVCopy, graphYieldPCMSysPi07TeVCopy,
                                                                                                             histoChargedPionSpecLowPtStat7TeVCMS, histoChargedPionSpecLowPtSys7TeVCMS, 
                                                                                                             kTRUE,  kTRUE, 
                                                                                                             &graphYieldPCMStatPi07TeVRebinnedCMSPCM, &graphYieldPCMSysPi07TeVRebinnedCMSPCM,
                                                                                                             &graphChargedPionSpecCMSStatPP7TeVCMSPCM, &graphChargedPionSpecCMSSystPP7TeVCMSPCM);
    graphRatioLowPtChargedPionsPCMPPCMS7TeV->Print();
    
    cout << "PHOS - high Pt" << endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPP7TeVHighPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPP7TeVHighPtPHOS    = NULL;
    TGraphErrors* graphYieldPHOSStatPi07TeVRebinnedHighPtPHOS       = NULL;
    TGraphErrors* graphYieldPHOSSysPi07TeVRebinnedHighPtPHOS        = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPHOSPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPHOSStatPi07TeVCopy, graphYieldPHOSSysPi07TeVCopy,
                                                                                                            histoChargedPionSpecHighPtStatPP7TeV, graphChargedPionSpecHighPtSystPP7TeV,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &graphYieldPHOSStatPi07TeVRebinnedHighPtPHOS, &graphChargedPionSpecHighPtSystPP7TeVHighPtPHOS, 
                                                                                                            &graphChargedPionSpecHighPtStatPP7TeVHighPtPHOS, &graphChargedPionSpecHighPtSystPP7TeVHighPtPHOS);
    graphRatioHighPtChargedPionsPHOSPP7TeV->Print();

    cout << "PHOS - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPP7TeVLowPtPHOS      = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPP7TeVLowPtPHOS      = NULL;
    TGraphErrors* graphYieldPHOSStatPi07TeVRebinnedLowPtPHOS        = NULL;
    TGraphErrors* graphYieldPHOSSysPi07TeVRebinnedLowPtPHOS         = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPHOSPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPHOSStatPi07TeVCopy, graphYieldPHOSSysPi07TeVCopy,
                                                                                                           histoChargedPionSpecLowPtStatPP7TeV, histoChargedPionSpecLowPtSysPP7TeV,
                                                                                                           kTRUE,  kTRUE, 
                                                                                                           &graphYieldPHOSStatPi07TeVRebinnedLowPtPHOS, &graphYieldPHOSSysPi07TeVRebinnedLowPtPHOS, 
                                                                                                           &graphChargedPionSpecLowPtStatPP7TeVLowPtPHOS, &graphChargedPionSpecLowPtSystPP7TeVLowPtPHOS);
    graphRatioLowPtChargedPionsPHOSPP7TeV->Print();

    cout<< "EMCal- high Pt"<< endl;
    TGraphErrors* graphChargedPionSpecHighPtStatPP7TeVHighPtEMCal   = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSystPP7TeVHighPtEMCal   = NULL;
    TGraphErrors* graphYieldEMCalStatPi07TeVRebinnedHighPtEMCal     = NULL;
    TGraphErrors* graphYieldEMCalSysPi07TeVRebinnedHighPtEMCal      = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsEMCalPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(histoEMCalPion7TeVPPStat, histoEMCalPion7TeVPPSyst, 
                                                                                                             histoChargedPionSpecHighPtStatPP7TeV, graphChargedPionSpecHighPtSystPP7TeV,  
                                                                                                             kTRUE,  kTRUE, 
                                                                                                             &graphYieldEMCalStatPi07TeVRebinnedHighPtEMCal, &graphChargedPionSpecHighPtSystPP7TeVHighPtEMCal,
                                                                                                             &graphChargedPionSpecHighPtStatPP7TeVHighPtEMCal, &graphChargedPionSpecHighPtSystPP7TeVHighPtEMCal);
    graphRatioHighPtChargedPionsEMCalPP7TeV->Print();

    cout << "EMCal - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPP7TeVLowPtEMCal     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPP7TeVLowPtEMCal     = NULL;
    TGraphErrors* graphYieldEMCalStatPi07TeVRebinnedLowPtEMCal      = NULL;
    TGraphErrors* graphYieldEMCalSysPi07TeVRebinnedLowPtEMCal       = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsEMCalPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(histoEMCalPion7TeVPPStat, histoEMCalPion7TeVPPSyst, 
                                                                                                            histoChargedPionSpecLowPtStatPP7TeV, histoChargedPionSpecLowPtSysPP7TeV, 
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &graphYieldEMCalStatPi07TeVRebinnedLowPtEMCal, &graphYieldEMCalSysPi07TeVRebinnedLowPtEMCal,
                                                                                                            &graphChargedPionSpecLowPtStatPP7TeVLowPtEMCal, &graphChargedPionSpecLowPtSystPP7TeVLowPtEMCal);
    graphRatioLowPtChargedPionsEMCalPP7TeV->Print();

    //    //***************************** ratios Dalitz 7 TeV ************************************************
    TH1D* histoInvCrossSectionPi0Dalitz7TeVCopy                     = (TH1D*) histoInvCrossSectionPi0Dalitz7TeV->Clone("histoInvCrossSectionPi0Dalitz7TeVCopy");
    TGraphAsymmErrors* graphInvCrossSectionSysPi0Dalitz7TeVCopy     = (TGraphAsymmErrors*) graphInvCrossSectionSysPi0Dalitz7TeV->Clone("graphInvCrossSectionSysPi0Dalitz7TeVCopy");
    //    
    cout << "Dalitz Spectrum - high Pt" << endl;
    TGraphErrors* graphYieldDalitzStatPi07TeVRebinnedHighPtDalitz   = NULL;
    TGraphErrors* graphYieldDalitzSysPi07TeVRebinnedHighPtDalitz    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat7TeVHighPtDalitz    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst7TeVHighPtDalitz    = NULL;   
    TGraphErrors* graphRatioHighPtChargedPionsDalitzPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(histoInvCrossSectionPi0Dalitz7TeVCopy, graphInvCrossSectionSysPi0Dalitz7TeVCopy, 
                                                                                                              histoChargedPionSpecHighPtStatPP7TeV, graphChargedPionSpecHighPtSystPP7TeV,  
                                                                                                              kTRUE,  kTRUE, 
                                                                                                              &graphYieldDalitzStatPi07TeVRebinnedHighPtDalitz, &graphYieldDalitzSysPi07TeVRebinnedHighPtDalitz, 
                                                                                                              &graphChargedPionSpecHighPtStat7TeVHighPtDalitz, &graphChargedPionSpecHighPtSyst7TeVHighPtDalitz);
    //    graphRatioHighPtChargedPionsDalitzPP7TeV->Print();

    cout << "Dalitz Spectrum - low Pt" << endl;
    TGraphErrors* graphYieldDalitzStatPi07TeVRebinnedLowPtDalitz    = NULL;
    TGraphErrors* graphYieldDalitzSysPi07TeVRebinnedLowPtDalitz     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat7TeVLowPtDalitz      = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst7TeVLowPtDalitz      = NULL;   

    TGraphErrors* graphRatioLowPtChargedPionsDalitzPP7TeV = CalculateRatioBetweenSpectraWithDifferentBinning(histoInvCrossSectionPi0Dalitz7TeVCopy, graphInvCrossSectionSysPi0Dalitz7TeVCopy, 
                                                                                                             histoChargedPionSpecLowPtStatPP7TeV, histoChargedPionSpecLowPtSysPP7TeV,  
                                                                                                             kTRUE,  kTRUE, 
                                                                                                             &graphYieldDalitzStatPi07TeVRebinnedLowPtDalitz, &graphYieldDalitzSysPi07TeVRebinnedLowPtDalitz,
                                                                                                             &graphChargedPionSpecLowPtStat7TeVLowPtDalitz, &graphChargedPionSpecLowPtSyst7TeVLowPtDalitz);
    //    graphRatioLowPtChargedPionsDalitzPP7TeV->Print();
    
    
    cout << "*************************************************************************"<< endl;    
    cout << "******************************  PP 0.9TeV *******************************"<< endl;
    cout << "*************************************************************************"<< endl;
    TGraphAsymmErrors* graphYieldCombStatPi0900GeVCopy      = (TGraphAsymmErrors*) graphInvYieldPi0Comb900GeVStatErr->Clone("graphYieldCombStatPi0900GeVCopy");
    TGraphAsymmErrors* graphYieldCombSysPi0900GeVCopy       = (TGraphAsymmErrors*) graphInvYieldPi0Comb900GeVSysErr->Clone("graphYieldCombSysPi0900GeVCopy");
    TGraphAsymmErrors* graphYieldPCMStatPi0900GeVCopy       = (TGraphAsymmErrors*) graphInvYieldPi0PCM900GeVStatErr->Clone("graphYieldPCMStatPi0900GeVCopy");
    TGraphAsymmErrors* graphYieldPCMSysPi0900GeVCopy        = (TGraphAsymmErrors*) graphInvYieldPi0PCM900GeVSysErr->Clone("graphYieldPCMSysPi0900GeVCopy");
    TGraphAsymmErrors* graphYieldPHOSStatPi0900GeVCopy      = (TGraphAsymmErrors*) graphInvYieldPi0PHOS900GeVStatErr->Clone("graphYieldPHOSStatPi0900GeVCopy");
    TGraphAsymmErrors* graphYieldPHOSSysPi0900GeVCopy       = (TGraphAsymmErrors*) graphInvYieldPi0PHOS900GeVSysErr->Clone("graphYieldPHOSSysPi0900GeVCopy");
    
    
    cout << "combined Spectrum - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPP900GeVLowPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPP900GeVLowPtComb    = NULL;
    TGraphErrors* graphYieldCombStatPi0900GeVRebinnedLowPtComb      = NULL;
    TGraphErrors* graphYieldCombSysPi0900GeVRebinnedLowPtComb       = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsCombPP900GeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombStatPi0900GeVCopy, graphYieldCombSysPi0900GeVCopy, 
                                                                                                             histoChargedPionSpecLowPtStatPP900GeV, histoChargedPionSpecLowPtSysPP900GeV, 
                                                                                                             kTRUE,  kTRUE, 
                                                                                                             &graphYieldCombStatPi0900GeVRebinnedLowPtComb, &graphYieldCombSysPi0900GeVRebinnedLowPtComb,
                                                                                                             &graphChargedPionSpecLowPtStatPP900GeVLowPtComb, &graphChargedPionSpecLowPtSystPP900GeVLowPtComb);
    graphRatioLowPtChargedPionsCombPP900GeV->Print();

    cout << "combined Spectrum - low Pt CMS" << endl;
    TGraphErrors* graphChargedPionSpecCMSStatPP900GeVCMSComb        = NULL;
    TGraphErrors* graphChargedPionSpecCMSSystPP900GeVCMSComb        = NULL;
    TGraphErrors* graphYieldCombStatPi0900GeVRebinnedCMSComb        = NULL;
    TGraphErrors* graphYieldCombSysPi0900GeVRebinnedCMSComb         = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsCombPPCMS900GeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombStatPi0900GeVCopy, graphYieldCombSysPi0900GeVCopy, 
                                                                                                                histoChargedPionSpecLowPtStat900GeVCMS, histoChargedPionSpecLowPtSys900GeVCMS, 
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphYieldCombStatPi0900GeVRebinnedCMSComb, &graphYieldCombSysPi0900GeVRebinnedCMSComb, 
                                                                                                                &graphChargedPionSpecCMSStatPP900GeVCMSComb, &graphChargedPionSpecCMSSystPP900GeVCMSComb);
    graphRatioLowPtChargedPionsCombPPCMS900GeV->Print();
        
    cout << "PCM - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPP900GeVLowPtPCM     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPP900GeVLowPtPCM     = NULL;
    TGraphErrors* graphYieldPCMStatPi0900GeVRebinnedLowPtPCM        = NULL;
    TGraphErrors* graphYieldPCMSysPi0900GeVRebinnedLowPtPCM         = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCMPP900GeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPCMStatPi0900GeVCopy, graphYieldPCMSysPi0900GeVCopy, 
                                                                                                            histoChargedPionSpecLowPtStatPP900GeV, histoChargedPionSpecLowPtSysPP900GeV,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &graphYieldPCMStatPi0900GeVRebinnedLowPtPCM, &graphYieldPCMSysPi0900GeVRebinnedLowPtPCM, 
                                                                                                            &graphChargedPionSpecLowPtStatPP900GeVLowPtPCM, &graphChargedPionSpecLowPtSystPP900GeVLowPtPCM);
    graphRatioLowPtChargedPionsPCMPP900GeV->Print();
    
    cout << "PCM - low Pt CMS" << endl; 
    TGraphErrors* graphChargedPionSpecCMSStatPP900GeVCMSPCM         = NULL;
    TGraphErrors* graphChargedPionSpecCMSSystPP900GeVCMSPCM         = NULL;
    TGraphErrors* graphYieldPCMStatPi0900GeVRebinnedCMSPCM          = NULL;
    TGraphErrors* graphYieldPCMSysPi0900GeVRebinnedCMSPCM           = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCMPPCMS900GeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPCMStatPi0900GeVCopy, graphYieldPCMSysPi0900GeVCopy, 
                                                                                                               histoChargedPionSpecLowPtStat900GeVCMS, histoChargedPionSpecLowPtSys900GeVCMS,  
                                                                                                               kTRUE,  kTRUE, 
                                                                                                               &graphYieldPCMStatPi0900GeVRebinnedCMSPCM, &graphYieldPCMSysPi0900GeVRebinnedCMSPCM, 
                                                                                                               &graphChargedPionSpecCMSStatPP900GeVCMSPCM, &graphChargedPionSpecCMSSystPP900GeVCMSPCM);
    graphRatioLowPtChargedPionsPCMPPCMS900GeV->Print();
    
    cout << "PHOS - low Pt" << endl;
    TGraphErrors* graphChargedPionSpecLowPtStatPP900GeVLowPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSystPP900GeVLowPtPHOS    = NULL;
    TGraphErrors* graphYieldPHOSStatPi0900GeVRebinnedLowPtPHOS      = NULL;
    TGraphErrors* graphYieldPHOSSysPi0900GeVRebinnedLowPtPHOS       = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPHOSPP900GeV = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPHOSStatPi0900GeVCopy, graphYieldPHOSSysPi0900GeVCopy, 
                                                                                                             histoChargedPionSpecLowPtStatPP900GeV, histoChargedPionSpecLowPtSysPP900GeV,  
                                                                                                             kTRUE, kTRUE, 
                                                                                                             &graphYieldPHOSStatPi0900GeVRebinnedLowPtPHOS, &graphYieldPHOSSysPi0900GeVRebinnedLowPtPHOS, 
                                                                                                             &graphChargedPionSpecLowPtStatPP900GeVLowPtPHOS, &graphChargedPionSpecLowPtSystPP900GeVLowPtPHOS);
    graphRatioLowPtChargedPionsPHOSPP900GeV->Print();
    
        
    cout << "*************************************************************************"<< endl;    
    cout << "******************************  PbPb 0-5% *******************************"<< endl;
    cout << "*************************************************************************"<< endl;
        
    TGraphAsymmErrors* graphYieldPi0CombPbPb0005StatErrCopy         = (TGraphAsymmErrors*) graphYieldPi0CombPbPb0005StatErr->Clone("graphYieldPi0CombPbPb0005StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0CombPbPb0005SysErrCopy          = (TGraphAsymmErrors*) graphYieldPi0CombPbPb0005SysErr->Clone("graphYieldPi0CombPbPb0005SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb0005StatErrCopy          = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb0005StatErr->Clone("graphYieldPi0PCMPbPb0005StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb0005SysErrCopy           = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb0005SysErr->Clone("graphYieldPi0PCMPbPb0005SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb0005StatErrCopy         = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb0005StatErr->Clone("graphYieldPi0PHOSPbPb0005StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb0005SysErrCopy          = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb0005SysErr->Clone("graphYieldPi0PHOSPbPb0005SysErrCopy");
    
    cout << "combined Spectrum" << endl;
    TGraphErrors* graphYieldCombStatPi00005RebinnedHighPtComb       = NULL;
    TGraphErrors* graphYieldCombSysPi00005RebinnedHighPtComb        = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat0005HighPtComb      = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst0005HighPtComb      = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsComb0005 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb0005StatErrCopy, graphYieldPi0CombPbPb0005SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat0005, histoChargedPionSpecHighPtSyst0005, 
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldCombStatPi00005RebinnedHighPtComb, &graphYieldCombSysPi00005RebinnedHighPtComb, 
                                                                                                          &graphChargedPionSpecHighPtStat0005HighPtComb, &graphChargedPionSpecHighPtSyst0005HighPtComb);
    graphRatioHighPtChargedPionsComb0005->Print();
   
    TGraphErrors* graphYieldCombStatPi00005RebinnedLowPtComb        = NULL;
    TGraphErrors* graphYieldCombSysPi00005RebinnedLowPtComb         = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat0005LowPtComb        = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst0005LowPtComb        = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsComb0005 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb0005StatErrCopy, graphYieldPi0CombPbPb0005SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat0005, histoChargedPionSpecLowPtSyst0005,  
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldCombStatPi00005RebinnedLowPtComb, &graphYieldCombSysPi00005RebinnedLowPtComb, 
                                                                                                         &graphChargedPionSpecLowPtStat0005LowPtComb, &graphChargedPionSpecLowPtSyst0005LowPtComb);
    graphRatioLowPtChargedPionsComb0005->Print();
    
    cout << "PCM" << endl;
    TGraphErrors* graphYieldPCMStatPi00005RebinnedHighPtPCM         = NULL;
    TGraphErrors* graphYieldPCMSysPi00005RebinnedHighPtPCM          = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat0005HighPtPCM       = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst0005HighPtPCM       = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCM0005 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb0005StatErrCopy, graphYieldPi0PCMPbPb0005SysErrCopy, 
                                                                                                         histoChargedPionSpecHighPtStat0005, histoChargedPionSpecHighPtSyst0005,  
                                                                                                         kTRUE,  kTRUE,
                                                                                                         &graphYieldPCMStatPi00005RebinnedHighPtPCM, &graphYieldPCMSysPi00005RebinnedHighPtPCM, 
                                                                                                         &graphChargedPionSpecHighPtStat0005HighPtPCM, &graphChargedPionSpecHighPtSyst0005HighPtPCM);
    graphRatioHighPtChargedPionsPCM0005->Print();

    TGraphErrors* graphYieldPCMStatPi00005RebinnedLowPtPCM          = NULL;
    TGraphErrors* graphYieldPCMSysPi00005RebinnedLowPtPCM           = NULL;   
    TGraphErrors* graphChargedPionSpecLowPtStat0005LowPtPCM         = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst0005LowPtPCM         = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCM0005 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb0005StatErrCopy, graphYieldPi0PCMPbPb0005SysErrCopy, 
                                                                                                        histoChargedPionSpecLowPtStat0005, histoChargedPionSpecLowPtSyst0005,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphYieldPCMStatPi00005RebinnedLowPtPCM, &graphYieldPCMSysPi00005RebinnedLowPtPCM, 
                                                                                                        &graphChargedPionSpecLowPtStat0005LowPtPCM, &graphChargedPionSpecLowPtSyst0005LowPtPCM);
    graphRatioLowPtChargedPionsPCM0005->Print();

    graphYieldPi0PHOSPbPb0005StatErrCopy->RemovePoint(0);
    graphYieldPi0PHOSPbPb0005SysErrCopy->RemovePoint(0);
    TGraphErrors* graphYieldPHOSStatPi00005RebinnedHighPtPHOS       = NULL;
    TGraphErrors* graphYieldPHOSSysPi00005RebinnedHighPtPHOS        = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat0005HighPtPHOS      = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst0005HighPtPHOS      = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPHOS0005 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb0005StatErrCopy, graphYieldPi0PHOSPbPb0005SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat0005, histoChargedPionSpecHighPtSyst0005,  
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldPHOSStatPi00005RebinnedHighPtPHOS, &graphYieldPHOSSysPi00005RebinnedHighPtPHOS, 
                                                                                                          &graphChargedPionSpecHighPtStat0005HighPtPHOS, &graphChargedPionSpecHighPtSyst0005HighPtPHOS);
    graphRatioHighPtChargedPionsPHOS0005->Print();
   
    TGraphErrors* graphYieldPHOSStatPi00005RebinnedLowPtPHOS        = NULL;
    TGraphErrors* graphYieldPHOSSysPi00005RebinnedLowPtPHOS         = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat0005LowPtPHOS        = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst0005LowPtPHOS        = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPHOS0005 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb0005StatErrCopy, graphYieldPi0PHOSPbPb0005SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat0005, histoChargedPionSpecLowPtSyst0005, 
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldPHOSStatPi00005RebinnedLowPtPHOS, &graphYieldPHOSSysPi00005RebinnedLowPtPHOS, 
                                                                                                         &graphChargedPionSpecLowPtStat0005LowPtPHOS, &graphChargedPionSpecLowPtSyst0005LowPtPHOS);
    graphRatioLowPtChargedPionsPHOS0005->Print();
    
//     return;
    
    cout << "*************************************************************************"<< endl;    
    cout << "******************************  PbPb 5-10% ******************************"<< endl;
    cout << "*************************************************************************"<< endl;
    
    TGraphAsymmErrors* graphYieldPi0CombPbPb0510StatErrCopy         = (TGraphAsymmErrors*) graphYieldPi0CombPbPb0510StatErr->Clone("graphYieldPi0CombPbPb0510StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0CombPbPb0510SysErrCopy          = (TGraphAsymmErrors*) graphYieldPi0CombPbPb0510SysErr->Clone("graphYieldPi0CombPbPb0510SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb0510StatErrCopy          = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb0510StatErr->Clone("graphYieldPi0PCMPbPb0510StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb0510SysErrCopy           = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb0510SysErr->Clone("graphYieldPi0PCMPbPb0510SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb0510StatErrCopy         = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb0510StatErr->Clone("graphYieldPi0PHOSPbPb0510StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb0510SysErrCopy          = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb0510SysErr->Clone("graphYieldPi0PHOSPbPb0510SysErrCopy");
    
    cout << "combined Spectrum" << endl;
    TGraphErrors* graphYieldCombStatPi00510RebinnedHighPtComb       = NULL;
    TGraphErrors* graphYieldCombSysPi00510RebinnedHighPtComb        = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat0510HighPtComb      = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst0510HighPtComb      = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsComb0510 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb0510StatErrCopy, graphYieldPi0CombPbPb0510SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat0510, histoChargedPionSpecHighPtSyst0510,  
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldCombStatPi00510RebinnedHighPtComb, &graphYieldCombSysPi00510RebinnedHighPtComb, 
                                                                                                          &graphChargedPionSpecHighPtStat0510HighPtComb, &graphChargedPionSpecHighPtSyst0510HighPtComb);
    graphRatioHighPtChargedPionsComb0510->Print();
    
    TGraphErrors* graphYieldCombStatPi00510RebinnedLowPtComb        = NULL;
    TGraphErrors* graphYieldCombSysPi00510RebinnedLowPtComb         = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat0510LowPtComb        = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst0510LowPtComb        = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsComb0510 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb0510StatErrCopy, graphYieldPi0CombPbPb0510SysErrCopy,
                                                                                                         histoChargedPionSpecLowPtStat0510, histoChargedPionSpecLowPtSyst0510,  
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldCombStatPi00510RebinnedLowPtComb, &graphYieldCombSysPi00510RebinnedLowPtComb, 
                                                                                                         &graphChargedPionSpecLowPtStat0510LowPtComb, &graphChargedPionSpecLowPtSyst0510LowPtComb);
    graphRatioLowPtChargedPionsComb0510->Print();
    
    
    cout << "PCM" << endl;
    TGraphErrors* graphYieldPCMStatPi00510RebinnedHighPtPCM         = NULL;
    TGraphErrors* graphYieldPCMSysPi00510RebinnedHighPtPCM          = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat0510HighPtPCM       = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst0510HighPtPCM       = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCM0510 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb0510StatErrCopy, graphYieldPi0PCMPbPb0510SysErrCopy, 
                                                                                                         histoChargedPionSpecHighPtStat0510, histoChargedPionSpecHighPtSyst0510, 
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldPCMStatPi00510RebinnedHighPtPCM, &graphYieldPCMSysPi00510RebinnedHighPtPCM, 
                                                                                                         &graphChargedPionSpecHighPtStat0510HighPtPCM, &graphChargedPionSpecHighPtSyst0510HighPtPCM);
    graphRatioHighPtChargedPionsPCM0510->Print();

    TGraphErrors* graphYieldPCMStatPi00510RebinnedLowPtPCM          = NULL;
    TGraphErrors* graphYieldPCMSysPi00510RebinnedLowPtPCM           = NULL;   
    TGraphErrors* graphChargedPionSpecLowPtStat0510LowPtPCM         = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst0510LowPtPCM         = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCM0510 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb0510StatErrCopy, graphYieldPi0PCMPbPb0510SysErrCopy, 
                                                                                                        histoChargedPionSpecLowPtStat0510, histoChargedPionSpecLowPtSyst0510,  
                                                                                                        kTRUE, kTRUE, 
                                                                                                        &graphYieldPCMStatPi00510RebinnedLowPtPCM, &graphYieldPCMSysPi00510RebinnedLowPtPCM, 
                                                                                                        &graphChargedPionSpecLowPtStat0510LowPtPCM, &graphChargedPionSpecLowPtSyst0510LowPtPCM);
    graphRatioLowPtChargedPionsPCM0510->Print();
    
    cout << endl << endl << endl << endl << endl<< "PHOS" << endl;
    graphYieldPi0PHOSPbPb0510StatErrCopy->RemovePoint(0);
    graphYieldPi0PHOSPbPb0510SysErrCopy->RemovePoint(0);
    TGraphErrors* graphYieldPHOSStatPi00510RebinnedHighPtPHOS       = NULL;
    TGraphErrors* graphYieldPHOSSysPi00510RebinnedHighPtPHOS        = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat0510HighPtPHOS      = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst0510HighPtPHOS      = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPHOS0510 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb0510StatErrCopy, graphYieldPi0PHOSPbPb0510SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat0510, histoChargedPionSpecHighPtSyst0510,  
                                                                                                          kTRUE, kTRUE, 
                                                                                                          &graphYieldPHOSStatPi00510RebinnedHighPtPHOS, &graphYieldPHOSSysPi00510RebinnedHighPtPHOS, 
                                                                                                          &graphChargedPionSpecHighPtStat0510HighPtPHOS, &graphChargedPionSpecHighPtSyst0510HighPtPHOS);
    graphRatioHighPtChargedPionsPHOS0510->Print();
    
    TGraphErrors* graphYieldPHOSStatPi00510RebinnedLowPtPHOS        = NULL;
    TGraphErrors* graphYieldPHOSSysPi00510RebinnedLowPtPHOS         = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat0510LowPtPHOS        = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst0510LowPtPHOS        = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPHOS0510 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb0510StatErrCopy, graphYieldPi0PHOSPbPb0510SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat0510, histoChargedPionSpecLowPtSyst0510, 
                                                                                                         kTRUE, kTRUE, 
                                                                                                         &graphYieldPHOSStatPi00510RebinnedLowPtPHOS, &graphYieldPHOSSysPi00510RebinnedLowPtPHOS, 
                                                                                                         &graphChargedPionSpecLowPtStat0510LowPtPHOS, &graphChargedPionSpecLowPtSyst0510LowPtPHOS);
    graphRatioLowPtChargedPionsPHOS0510->Print();
    
    
    cout << "*************************************************************************"<< endl;    
    cout << "******************************  PbPb 10-20% *****************************"<< endl;
    cout << "*************************************************************************"<< endl;
    
    TGraphAsymmErrors* graphYieldPi0CombPbPb1020StatErrCopy     = (TGraphAsymmErrors*) graphYieldPi0CombPbPb1020StatErr->Clone("graphYieldPi0CombPbPb1020StatErrCopy");
    cout << "comb" << endl;
    graphYieldPi0CombPbPb1020StatErrCopy->Print();
    TGraphAsymmErrors* graphYieldPi0CombPbPb1020SysErrCopy      = (TGraphAsymmErrors*) graphYieldPi0CombPbPb1020SysErr->Clone("graphYieldPi0CombPbPb1020SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb1020StatErrCopy      = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb1020StatErr->Clone("graphYieldPi0PCMPbPb1020StatErrCopy");
    cout << "PCM" << endl;
    graphYieldPi0PCMPbPb1020StatErrCopy->Print();
    TGraphAsymmErrors* graphYieldPi0PCMPbPb1020SysErrCopy       = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb1020SysErr->Clone("graphYieldPi0PCMPbPb1020SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb1020StatErrCopy     = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb1020StatErr->Clone("graphYieldPi0PHOSPbPb1020StatErrCopy");
    cout << "PHOS" << endl;
    graphYieldPi0PHOSPbPb1020StatErrCopy->Print();
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb1020SysErrCopy      = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb1020SysErr->Clone("graphYieldPi0PHOSPbPb1020SysErrCopy");
    
    cout << "combined Spectrum" << endl;
    TGraphErrors* graphYieldCombStatPi01020RebinnedHighPtComb   = NULL;
    TGraphErrors* graphYieldCombSysPi01020RebinnedHighPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat1020HighPtComb  = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst1020HighPtComb  = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsComb1020 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb1020StatErrCopy, graphYieldPi0CombPbPb1020SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat1020, histoChargedPionSpecHighPtSyst1020, 
                                                                                                          kTRUE, kTRUE, 
                                                                                                          &graphYieldCombStatPi01020RebinnedHighPtComb, &graphYieldCombSysPi01020RebinnedHighPtComb, 
                                                                                                          &graphChargedPionSpecHighPtStat1020HighPtComb, &graphChargedPionSpecHighPtSyst1020HighPtComb);
    graphRatioHighPtChargedPionsComb1020->Print();
    
    TGraphErrors* graphYieldCombStatPi01020RebinnedLowPtComb    = NULL;
    TGraphErrors* graphYieldCombSysPi01020RebinnedLowPtComb     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat1020LowPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst1020LowPtComb    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsComb1020 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb1020StatErrCopy, graphYieldPi0CombPbPb1020SysErrCopy,
                                                                                                         histoChargedPionSpecLowPtStat1020, histoChargedPionSpecLowPtSyst1020,  
                                                                                                         kTRUE, kTRUE, 
                                                                                                         &graphYieldCombStatPi01020RebinnedLowPtComb, &graphYieldCombSysPi01020RebinnedLowPtComb, 
                                                                                                         &graphChargedPionSpecLowPtStat1020LowPtComb, &graphChargedPionSpecLowPtSyst1020LowPtComb);
    graphRatioLowPtChargedPionsComb1020->Print();
    
    
    cout << "PCM" << endl;
    TGraphErrors* graphYieldPCMStatPi01020RebinnedHighPtPCM     = NULL;
    TGraphErrors* graphYieldPCMSysPi01020RebinnedHighPtPCM      = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat1020HighPtPCM   = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst1020HighPtPCM   = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCM1020 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb1020StatErrCopy, graphYieldPi0PCMPbPb1020SysErrCopy, 
                                                                                                         histoChargedPionSpecHighPtStat1020, histoChargedPionSpecHighPtSyst1020,  
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldPCMStatPi01020RebinnedHighPtPCM, &graphYieldPCMSysPi01020RebinnedHighPtPCM, 
                                                                                                         &graphChargedPionSpecHighPtStat1020HighPtPCM, &graphChargedPionSpecHighPtSyst1020HighPtPCM);
    graphRatioHighPtChargedPionsPCM1020->Print();

    TGraphErrors* graphYieldPCMStatPi01020RebinnedLowPtPCM      = NULL;
    TGraphErrors* graphYieldPCMSysPi01020RebinnedLowPtPCM       = NULL;   
    TGraphErrors* graphChargedPionSpecLowPtStat1020LowPtPCM     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst1020LowPtPCM     = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCM1020 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb1020StatErrCopy, graphYieldPi0PCMPbPb1020SysErrCopy, 
                                                                                                        histoChargedPionSpecLowPtStat1020, histoChargedPionSpecLowPtSyst1020,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphYieldPCMStatPi01020RebinnedLowPtPCM, &graphYieldPCMSysPi01020RebinnedLowPtPCM, 
                                                                                                        &graphChargedPionSpecLowPtStat1020LowPtPCM, &graphChargedPionSpecLowPtSyst1020LowPtPCM);
    graphRatioLowPtChargedPionsPCM1020->Print();
    
    cout << endl << endl << endl << endl << endl<< "PHOS" << endl;
    graphYieldPi0PHOSPbPb1020StatErrCopy->RemovePoint(0);
    graphYieldPi0PHOSPbPb1020SysErrCopy->RemovePoint(0);
    TGraphErrors* graphYieldPHOSStatPi01020RebinnedHighPtPHOS   = NULL;
    TGraphErrors* graphYieldPHOSSysPi01020RebinnedHighPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat1020HighPtPHOS  = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst1020HighPtPHOS  = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPHOS1020 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb1020StatErrCopy, graphYieldPi0PHOSPbPb1020SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat1020, histoChargedPionSpecHighPtSyst1020, 
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldPHOSStatPi01020RebinnedHighPtPHOS, &graphYieldPHOSSysPi01020RebinnedHighPtPHOS, 
                                                                                                          &graphChargedPionSpecHighPtStat1020HighPtPHOS, &graphChargedPionSpecHighPtSyst1020HighPtPHOS);
    graphRatioHighPtChargedPionsPHOS1020->Print();
    
    TGraphErrors* graphYieldPHOSStatPi01020RebinnedLowPtPHOS    = NULL;
    TGraphErrors* graphYieldPHOSSysPi01020RebinnedLowPtPHOS     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat1020LowPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst1020LowPtPHOS    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPHOS1020 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb1020StatErrCopy, graphYieldPi0PHOSPbPb1020SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat1020, histoChargedPionSpecLowPtSyst1020,  
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldPHOSStatPi01020RebinnedLowPtPHOS, &graphYieldPHOSSysPi01020RebinnedLowPtPHOS, 
                                                                                                         &graphChargedPionSpecLowPtStat1020LowPtPHOS, &graphChargedPionSpecLowPtSyst1020LowPtPHOS);
    graphRatioLowPtChargedPionsPHOS1020->Print();
    
    cout << "*************************************************************************"<< endl;    
    cout << "******************************  PbPb 20-40% *****************************"<< endl;
    cout << "*************************************************************************"<< endl;
    
    TGraphAsymmErrors* graphYieldPi0CombPbPb2040StatErrCopy     = (TGraphAsymmErrors*) graphYieldPi0CombPbPb2040StatErr->Clone("graphYieldPi0CombPbPb2040StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0CombPbPb2040SysErrCopy      = (TGraphAsymmErrors*) graphYieldPi0CombPbPb2040SysErr->Clone("graphYieldPi0CombPbPb2040SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb2040StatErrCopy      = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb2040StatErr->Clone("graphYieldPi0PCMPbPb2040StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb2040SysErrCopy       = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb2040SysErr->Clone("graphYieldPi0PCMPbPb2040SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb2040StatErrCopy     = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb2040StatErr->Clone("graphYieldPi0PHOSPbPb2040StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb2040SysErrCopy      = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb2040SysErr->Clone("graphYieldPi0PHOSPbPb2040SysErrCopy");
    
    cout << "combined Spectrum" << endl;
    TGraphErrors* graphYieldCombStatPi02040RebinnedHighPtComb   = NULL;
    TGraphErrors* graphYieldCombSysPi02040RebinnedHighPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat2040HighPtComb  = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst2040HighPtComb  = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsComb2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb2040StatErrCopy, graphYieldPi0CombPbPb2040SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat2040, histoChargedPionSpecHighPtSyst2040, 
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldCombStatPi02040RebinnedHighPtComb, &graphYieldCombSysPi02040RebinnedHighPtComb, 
                                                                                                          &graphChargedPionSpecHighPtStat2040HighPtComb, &graphChargedPionSpecHighPtSyst2040HighPtComb);
    graphRatioHighPtChargedPionsComb2040->Print();
    
    TGraphErrors* graphYieldCombStatPi02040RebinnedLowPtComb    = NULL;
    TGraphErrors* graphYieldCombSysPi02040RebinnedLowPtComb     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat2040LowPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst2040LowPtComb    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsComb2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb2040StatErrCopy, graphYieldPi0CombPbPb2040SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat2040, histoChargedPionSpecLowPtSyst2040,  
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldCombStatPi02040RebinnedLowPtComb, &graphYieldCombSysPi02040RebinnedLowPtComb, 
                                                                                                         &graphChargedPionSpecLowPtStat2040LowPtComb, &graphChargedPionSpecLowPtSyst2040LowPtComb);
    graphRatioLowPtChargedPionsComb2040->Print();
    
    cout << "PCM" << endl;
    TGraphErrors* graphYieldPCMStatPi02040RebinnedHighPtPCM     = NULL;
    TGraphErrors* graphYieldPCMSysPi02040RebinnedHighPtPCM      = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat2040HighPtPCM   = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst2040HighPtPCM   = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCM2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb2040StatErrCopy, graphYieldPi0PCMPbPb2040SysErrCopy, 
                                                                                                         histoChargedPionSpecHighPtStat2040, histoChargedPionSpecHighPtSyst2040,  
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldPCMStatPi02040RebinnedHighPtPCM, &graphYieldPCMSysPi02040RebinnedHighPtPCM, 
                                                                                                         &graphChargedPionSpecHighPtStat2040HighPtPCM, &graphChargedPionSpecHighPtSyst2040HighPtPCM);
    graphRatioHighPtChargedPionsPCM2040->Print();

    TGraphErrors* graphYieldPCMStatPi02040RebinnedLowPtPCM      = NULL;
    TGraphErrors* graphYieldPCMSysPi02040RebinnedLowPtPCM       = NULL;   
    TGraphErrors* graphChargedPionSpecLowPtStat2040LowPtPCM     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst2040LowPtPCM     = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCM2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb2040StatErrCopy, graphYieldPi0PCMPbPb2040SysErrCopy, 
                                                                                                        histoChargedPionSpecLowPtStat2040, histoChargedPionSpecLowPtSyst2040,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphYieldPCMStatPi02040RebinnedLowPtPCM, &graphYieldPCMSysPi02040RebinnedLowPtPCM, 
                                                                                                        &graphChargedPionSpecLowPtStat2040LowPtPCM, &graphChargedPionSpecLowPtSyst2040LowPtPCM);
    graphRatioLowPtChargedPionsPCM2040->Print();
    
    cout << endl << endl << endl << endl << endl<< "PHOS" << endl;
    graphYieldPi0PHOSPbPb2040StatErrCopy->RemovePoint(0);
    graphYieldPi0PHOSPbPb2040SysErrCopy->RemovePoint(0);
    TGraphErrors* graphYieldPHOSStatPi02040RebinnedHighPtPHOS   = NULL;
    TGraphErrors* graphYieldPHOSSysPi02040RebinnedHighPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat2040HighPtPHOS  = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst2040HighPtPHOS  = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPHOS2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb2040StatErrCopy, graphYieldPi0PHOSPbPb2040SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat2040, histoChargedPionSpecHighPtSyst2040, 
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldPHOSStatPi02040RebinnedHighPtPHOS, &graphYieldPHOSSysPi02040RebinnedHighPtPHOS, 
                                                                                                          &graphChargedPionSpecHighPtStat2040HighPtPHOS, &graphChargedPionSpecHighPtSyst2040HighPtPHOS);
    graphRatioHighPtChargedPionsPHOS2040->Print();
    
    TGraphErrors* graphYieldPHOSStatPi02040RebinnedLowPtPHOS    = NULL;
    TGraphErrors* graphYieldPHOSSysPi02040RebinnedLowPtPHOS     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat2040LowPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst2040LowPtPHOS    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPHOS2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb2040StatErrCopy, graphYieldPi0PHOSPbPb2040SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat2040, histoChargedPionSpecLowPtSyst2040,  
                                                                                                         kTRUE,  kTRUE,
                                                                                                         &graphYieldPHOSStatPi02040RebinnedLowPtPHOS, &graphYieldPHOSSysPi02040RebinnedLowPtPHOS, 
                                                                                                         &graphChargedPionSpecLowPtStat2040LowPtPHOS, &graphChargedPionSpecLowPtSyst2040LowPtPHOS);
    graphRatioLowPtChargedPionsPHOS2040->Print();

    cout << "*************************************************************************"<< endl;    
    cout << "******************************  PbPb 40-60% *****************************"<< endl;
    cout << "*************************************************************************"<< endl;
    
    TGraphAsymmErrors* graphYieldPi0CombPbPb4060StatErrCopy     = (TGraphAsymmErrors*) graphYieldPi0CombPbPb4060StatErr->Clone("graphYieldPi0CombPbPb4060StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0CombPbPb4060SysErrCopy      = (TGraphAsymmErrors*) graphYieldPi0CombPbPb4060SysErr->Clone("graphYieldPi0CombPbPb4060SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb4060StatErrCopy      = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb4060StatErr->Clone("graphYieldPi0PCMPbPb4060StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb4060SysErrCopy       = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb4060SysErr->Clone("graphYieldPi0PCMPbPb4060SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb4060StatErrCopy     = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb4060StatErr->Clone("graphYieldPi0PHOSPbPb4060StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb4060SysErrCopy      = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb4060SysErr->Clone("graphYieldPi0PHOSPbPb4060SysErrCopy");
    
    cout << endl << "*************************************************************************"<< endl;  
    cout << "******************************Comb***************************************"<< endl;  
    cout << "*************************************************************************"<< endl << endl;  
    cout << "combined Spectrum" << endl;
    TGraphErrors* graphYieldCombStatPi04060RebinnedHighPtComb   = NULL;
    TGraphErrors* graphYieldCombSysPi04060RebinnedHighPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat4060HighPtComb  = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst4060HighPtComb  = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsComb4060 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb4060StatErrCopy, graphYieldPi0CombPbPb4060SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat4060, histoChargedPionSpecHighPtSyst4060,  
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldCombStatPi04060RebinnedHighPtComb, &graphYieldCombSysPi04060RebinnedHighPtComb, 
                                                                                                          &graphChargedPionSpecHighPtStat4060HighPtComb, &graphChargedPionSpecHighPtSyst4060HighPtComb);
    graphRatioHighPtChargedPionsComb4060->Print();
    
    TGraphErrors* graphYieldCombStatPi04060RebinnedLowPtComb    = NULL;
    TGraphErrors* graphYieldCombSysPi04060RebinnedLowPtComb     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat4060LowPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst4060LowPtComb    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsComb4060 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb4060StatErrCopy, graphYieldPi0CombPbPb4060SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat4060, histoChargedPionSpecLowPtSyst4060, 
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldCombStatPi04060RebinnedLowPtComb, &graphYieldCombSysPi04060RebinnedLowPtComb, 
                                                                                                         &graphChargedPionSpecLowPtStat4060LowPtComb, &graphChargedPionSpecLowPtSyst4060LowPtComb);
    graphRatioLowPtChargedPionsComb4060->Print();
    
    
    cout << "PCM" << endl;
    TGraphErrors* graphYieldPCMStatPi04060RebinnedHighPtPCM     = NULL;
    TGraphErrors* graphYieldPCMSysPi04060RebinnedHighPtPCM      = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat4060HighPtPCM   = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst4060HighPtPCM   = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCM4060 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb4060StatErrCopy, graphYieldPi0PCMPbPb4060SysErrCopy, 
                                                                                                         histoChargedPionSpecHighPtStat4060, histoChargedPionSpecHighPtSyst4060,  
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldPCMStatPi04060RebinnedHighPtPCM, &graphYieldPCMSysPi04060RebinnedHighPtPCM, 
                                                                                                         &graphChargedPionSpecHighPtStat4060HighPtPCM, &graphChargedPionSpecHighPtSyst4060HighPtPCM);
    graphRatioHighPtChargedPionsPCM4060->Print();

    TGraphErrors* graphYieldPCMStatPi04060RebinnedLowPtPCM      = NULL;
    TGraphErrors* graphYieldPCMSysPi04060RebinnedLowPtPCM       = NULL;   
    TGraphErrors* graphChargedPionSpecLowPtStat4060LowPtPCM     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst4060LowPtPCM     = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCM4060 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb4060StatErrCopy, graphYieldPi0PCMPbPb4060SysErrCopy,
                                                                                                        histoChargedPionSpecLowPtStat4060, histoChargedPionSpecLowPtSyst4060,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphYieldPCMStatPi04060RebinnedLowPtPCM, &graphYieldPCMSysPi04060RebinnedLowPtPCM, 
                                                                                                        &graphChargedPionSpecLowPtStat4060LowPtPCM, &graphChargedPionSpecLowPtSyst4060LowPtPCM);
    graphRatioLowPtChargedPionsPCM4060->Print();
    
    cout << endl << endl << endl << endl << endl<< "PHOS" << endl;
    graphYieldPi0PHOSPbPb4060StatErrCopy->RemovePoint(0);
    graphYieldPi0PHOSPbPb4060SysErrCopy->RemovePoint(0);
    TGraphErrors* graphYieldPHOSStatPi04060RebinnedHighPtPHOS   = NULL;
    TGraphErrors* graphYieldPHOSSysPi04060RebinnedHighPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat4060HighPtPHOS  = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst4060HighPtPHOS  = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPHOS4060 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb4060StatErrCopy, graphYieldPi0PHOSPbPb4060SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat4060, histoChargedPionSpecHighPtSyst4060,  
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldPHOSStatPi04060RebinnedHighPtPHOS, &graphYieldPHOSSysPi04060RebinnedHighPtPHOS, 
                                                                                                          &graphChargedPionSpecHighPtStat4060HighPtPHOS, &graphChargedPionSpecHighPtSyst4060HighPtPHOS);
    graphRatioHighPtChargedPionsPHOS4060->Print();
    
    TGraphErrors* graphYieldPHOSStatPi04060RebinnedLowPtPHOS    = NULL;
    TGraphErrors* graphYieldPHOSSysPi04060RebinnedLowPtPHOS     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat4060LowPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst4060LowPtPHOS    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPHOS4060 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb4060StatErrCopy, graphYieldPi0PHOSPbPb4060SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat4060, histoChargedPionSpecLowPtSyst4060,  
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldPHOSStatPi04060RebinnedLowPtPHOS, &graphYieldPHOSSysPi04060RebinnedLowPtPHOS, 
                                                                                                         &graphChargedPionSpecLowPtStat4060LowPtPHOS, &graphChargedPionSpecLowPtSyst4060LowPtPHOS);
    graphRatioLowPtChargedPionsPHOS4060->Print();
    
    cout << "*************************************************************************"<< endl;    
    cout << "******************************  PbPb 60-80% *****************************"<< endl;
    cout << "*************************************************************************"<< endl;
    
    TGraphAsymmErrors* graphYieldPi0CombPbPb6080StatErrCopy     = (TGraphAsymmErrors*) graphYieldPi0CombPbPb6080StatErr->Clone("graphYieldPi0CombPbPb6080StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0CombPbPb6080SysErrCopy      = (TGraphAsymmErrors*) graphYieldPi0CombPbPb6080SysErr->Clone("graphYieldPi0CombPbPb6080SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb6080StatErrCopy      = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb6080StatErr->Clone("graphYieldPi0PCMPbPb6080StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PCMPbPb6080SysErrCopy       = (TGraphAsymmErrors*) graphYieldPi0PCMPbPb6080SysErr->Clone("graphYieldPi0PCMPbPb6080SysErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb6080StatErrCopy     = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb6080StatErr->Clone("graphYieldPi0PHOSPbPb6080StatErrCopy");
    TGraphAsymmErrors* graphYieldPi0PHOSPbPb6080SysErrCopy      = (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb6080SysErr->Clone("graphYieldPi0PHOSPbPb6080SysErrCopy");
    
    cout << "combined Spectrum" << endl;
    TGraphErrors* graphYieldCombStatPi06080RebinnedHighPtComb   = NULL;
    TGraphErrors* graphYieldCombSysPi06080RebinnedHighPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat6080HighPtComb  = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst6080HighPtComb  = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsComb6080 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb6080StatErrCopy, graphYieldPi0CombPbPb6080SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat6080, histoChargedPionSpecHighPtSyst6080, 
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldCombStatPi06080RebinnedHighPtComb, &graphYieldCombSysPi06080RebinnedHighPtComb, 
                                                                                                          &graphChargedPionSpecHighPtStat6080HighPtComb, &graphChargedPionSpecHighPtSyst6080HighPtComb);
    graphRatioHighPtChargedPionsComb6080->Print();
    
    TGraphErrors* graphYieldCombStatPi06080RebinnedLowPtComb    = NULL;
    TGraphErrors* graphYieldCombSysPi06080RebinnedLowPtComb     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat6080LowPtComb    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst6080LowPtComb    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsComb6080 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb6080StatErrCopy, graphYieldPi0CombPbPb6080SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat6080, histoChargedPionSpecLowPtSyst6080, 
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldCombStatPi06080RebinnedLowPtComb, &graphYieldCombSysPi06080RebinnedLowPtComb, 
                                                                                                         &graphChargedPionSpecLowPtStat6080LowPtComb, &graphChargedPionSpecLowPtSyst6080LowPtComb);
    graphRatioLowPtChargedPionsComb6080->Print();
    
    
    cout << "PCM" << endl;
    TGraphErrors* graphYieldPCMStatPi06080RebinnedHighPtPCM     = NULL;
    TGraphErrors* graphYieldPCMSysPi06080RebinnedHighPtPCM      = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat6080HighPtPCM   = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst6080HighPtPCM   = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPCM6080 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb6080StatErrCopy, graphYieldPi0PCMPbPb6080SysErrCopy, 
                                                                                                         histoChargedPionSpecHighPtStat6080, histoChargedPionSpecHighPtSyst6080, 
                                                                                                         kTRUE,  kTRUE,
                                                                                                         &graphYieldPCMStatPi06080RebinnedHighPtPCM, &graphYieldPCMSysPi06080RebinnedHighPtPCM, 
                                                                                                         &graphChargedPionSpecHighPtStat6080HighPtPCM, &graphChargedPionSpecHighPtSyst6080HighPtPCM);
    graphRatioHighPtChargedPionsPCM6080->Print();

    TGraphErrors* graphYieldPCMStatPi06080RebinnedLowPtPCM      = NULL;
    TGraphErrors* graphYieldPCMSysPi06080RebinnedLowPtPCM       = NULL;   
    TGraphErrors* graphChargedPionSpecLowPtStat6080LowPtPCM     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst6080LowPtPCM     = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPCM6080 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb6080StatErrCopy, graphYieldPi0PCMPbPb6080SysErrCopy, 
                                                                                                        histoChargedPionSpecLowPtStat6080, histoChargedPionSpecLowPtSyst6080, 
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphYieldPCMStatPi06080RebinnedLowPtPCM, &graphYieldPCMSysPi06080RebinnedLowPtPCM, 
                                                                                                        &graphChargedPionSpecLowPtStat6080LowPtPCM, &graphChargedPionSpecLowPtSyst6080LowPtPCM);
    graphRatioLowPtChargedPionsPCM6080->Print();
    
    cout << endl << endl << endl << endl << endl<< "PHOS" << endl;
    graphYieldPi0PHOSPbPb6080StatErrCopy->RemovePoint(0);
    graphYieldPi0PHOSPbPb6080SysErrCopy->RemovePoint(0);
    TGraphErrors* graphYieldPHOSStatPi06080RebinnedHighPtPHOS   = NULL;
    TGraphErrors* graphYieldPHOSSysPi06080RebinnedHighPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat6080HighPtPHOS  = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst6080HighPtPHOS  = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsPHOS6080 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb6080StatErrCopy, graphYieldPi0PHOSPbPb6080SysErrCopy, 
                                                                                                          histoChargedPionSpecHighPtStat6080, histoChargedPionSpecHighPtSyst6080,  
                                                                                                          kTRUE,  kTRUE, 
                                                                                                          &graphYieldPHOSStatPi06080RebinnedHighPtPHOS, &graphYieldPHOSSysPi06080RebinnedHighPtPHOS, 
                                                                                                          &graphChargedPionSpecHighPtStat6080HighPtPHOS, &graphChargedPionSpecHighPtSyst6080HighPtPHOS);
    graphRatioHighPtChargedPionsPHOS6080->Print();
    
    TGraphErrors* graphYieldPHOSStatPi06080RebinnedLowPtPHOS    = NULL;
    TGraphErrors* graphYieldPHOSSysPi06080RebinnedLowPtPHOS     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat6080LowPtPHOS    = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst6080LowPtPHOS    = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsPHOS6080 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb6080StatErrCopy, graphYieldPi0PHOSPbPb6080SysErrCopy, 
                                                                                                         histoChargedPionSpecLowPtStat6080, histoChargedPionSpecLowPtSyst6080,  
                                                                                                         kTRUE,  kTRUE, 
                                                                                                         &graphYieldPHOSStatPi06080RebinnedLowPtPHOS, &graphYieldPHOSSysPi06080RebinnedLowPtPHOS, 
                                                                                                         &graphChargedPionSpecLowPtStat6080LowPtPHOS, &graphChargedPionSpecLowPtSyst6080LowPtPHOS);
    graphRatioLowPtChargedPionsPHOS6080->Print();

    cout << "Dalitz 20-40% PbPb" << endl;   
    //***************************** ratios Dalitz 20-40% ************************************************
    TGraphErrors* graphYieldDalitzStatPi02040RebinnedHighPtDalitz   = NULL;
    TGraphErrors* graphYieldDalitzSysPi02040RebinnedHighPtDalitz    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat2040HighPtDalitz    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst2040HighPtDalitz    = NULL;
    TH1D* histoYieldPi0DalitzPbPb2040Copy                           = (TH1D*) histoYieldPi0DalitzPbPb2040->Clone("histoYieldPi0DalitzPbPb2040Copy");
    TGraphAsymmErrors* graphYieldSysPi0DalitzPbPb2040Copy           = (TGraphAsymmErrors*) graphYieldSysPi0DalitzPbPb2040->Clone("graphYieldSysPi0DalitzPbPb2040Copy");
    TGraphErrors* graphRatioHighPtChargedPionsDalitz2040 = CalculateRatioBetweenSpectraWithDifferentBinning(histoYieldPi0DalitzPbPb2040Copy, graphYieldSysPi0DalitzPbPb2040Copy,
                                                                                                            histoChargedPionSpecHighPtStat2040, histoChargedPionSpecHighPtSyst2040,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &graphYieldDalitzStatPi02040RebinnedHighPtDalitz, &graphYieldDalitzSysPi02040RebinnedHighPtDalitz,
                                                                                                            &graphChargedPionSpecHighPtStat2040HighPtDalitz, &graphChargedPionSpecHighPtSyst2040HighPtDalitz);
    
    TGraphErrors* graphYieldDalitzStatPi02040RebinnedLowPtDalitz    = NULL;
    TGraphErrors* graphYieldDalitzSysPi02040RebinnedLowPtDalitz     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat2040LowPtDalitz      = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst2040LowPtDalitz      = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsDalitz2040 = CalculateRatioBetweenSpectraWithDifferentBinning(histoYieldPi0DalitzPbPb2040Copy, graphYieldSysPi0DalitzPbPb2040Copy, 
                                                                                                           histoChargedPionSpecLowPtStat2040, histoChargedPionSpecLowPtSyst2040,  
                                                                                                           kTRUE,  kTRUE, 
                                                                                                           &graphYieldDalitzStatPi02040RebinnedLowPtDalitz, &graphYieldDalitzSysPi02040RebinnedLowPtDalitz,
                                                                                                           &graphChargedPionSpecLowPtStat2040LowPtDalitz, &graphChargedPionSpecLowPtSyst2040LowPtDalitz);

    cout << "Dalitz 40-60% PbPb" << endl;   
    //***************************** ratios Dalitz 40-60% ************************************************    
    TH1D* histoYieldPi0DalitzPbPb4060Copy                           = (TH1D*) histoYieldPi0DalitzPbPb4060->Clone("histoYieldPi0DalitzPbPb4060Copy");
    TGraphAsymmErrors* graphYieldSysPi0DalitzPbPb4060Copy           = (TGraphAsymmErrors*) graphYieldSysPi0DalitzPbPb4060->Clone("graphYieldSysPi0DalitzPbPb4060Copy");
    TGraphErrors* graphYieldDalitzStatPi04060RebinnedHighPtDalitz   = NULL;
    TGraphErrors* graphYieldDalitzSysPi04060RebinnedHighPtDalitz    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat4060HighPtDalitz    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst4060HighPtDalitz    = NULL;
    TGraphErrors* graphRatioHighPtChargedPionsDalitz4060 = CalculateRatioBetweenSpectraWithDifferentBinning(histoYieldPi0DalitzPbPb4060Copy, graphYieldSysPi0DalitzPbPb4060Copy,
                                                                                                            histoChargedPionSpecHighPtStat4060, histoChargedPionSpecHighPtSyst4060,
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &graphYieldDalitzStatPi04060RebinnedHighPtDalitz, &graphYieldDalitzSysPi04060RebinnedHighPtDalitz,
                                                                                                            &graphChargedPionSpecHighPtStat4060HighPtDalitz, &graphChargedPionSpecHighPtSyst4060HighPtDalitz);

    TGraphErrors* graphYieldDalitzStatPi04060RebinnedLowPtDalitz    = NULL;
    TGraphErrors* graphYieldDalitzSysPi04060RebinnedLowPtDalitz     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat4060LowPtDalitz      = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst4060LowPtDalitz      = NULL;
    TGraphErrors* graphRatioLowPtChargedPionsDalitz4060 = CalculateRatioBetweenSpectraWithDifferentBinning(histoYieldPi0DalitzPbPb4060Copy, graphYieldSysPi0DalitzPbPb4060Copy, 
                                                                                                           histoChargedPionSpecLowPtStat4060, histoChargedPionSpecLowPtSyst4060,  
                                                                                                           kTRUE,  kTRUE, 
                                                                                                           &graphYieldDalitzStatPi04060RebinnedLowPtDalitz, &graphYieldDalitzSysPi04060RebinnedLowPtDalitz,
                                                                                                           &graphChargedPionSpecLowPtStat4060LowPtDalitz, &graphChargedPionSpecLowPtSyst4060LowPtDalitz);
    
    cout << "Dalitz 60-80% PbPb" << endl;   
    //***************************** ratios Dalitz 60-80% ************************************************
    TH1D* histoYieldPi0DalitzPbPb6080Copy                           = (TH1D*) histoYieldPi0DalitzPbPb6080->Clone("histoYieldPi0DalitzPbPb6080Copy");
    TGraphAsymmErrors* graphYieldSysPi0DalitzPbPb6080Copy           = (TGraphAsymmErrors*) graphYieldSysPi0DalitzPbPb6080->Clone("graphYieldSysPi0DalitzPbPb6080Copy");
    TGraphErrors* graphYieldDalitzStatPi06080RebinnedHighPtDalitz   = NULL;
    TGraphErrors* graphYieldDalitzSysPi06080RebinnedHighPtDalitz    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtStat6080HighPtDalitz    = NULL;
    TGraphErrors* graphChargedPionSpecHighPtSyst6080HighPtDalitz    = NULL;   
    TGraphErrors* graphRatioHighPtChargedPionsDalitz6080 = CalculateRatioBetweenSpectraWithDifferentBinning(histoYieldPi0DalitzPbPb6080Copy, graphYieldSysPi0DalitzPbPb6080Copy,
                                                                                                            histoChargedPionSpecHighPtStat6080, histoChargedPionSpecHighPtSyst6080,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &graphYieldDalitzStatPi06080RebinnedHighPtDalitz, &graphYieldDalitzSysPi06080RebinnedHighPtDalitz,
                                                                                                            &graphChargedPionSpecHighPtStat6080HighPtDalitz, &graphChargedPionSpecHighPtSyst6080HighPtDalitz);
    
    TGraphErrors* graphYieldDalitzStatPi06080RebinnedLowPtDalitz    = NULL;
    TGraphErrors* graphYieldDalitzSysPi06080RebinnedLowPtDalitz     = NULL;
    TGraphErrors* graphChargedPionSpecLowPtStat6080LowPtDalitz      = NULL;
    TGraphErrors* graphChargedPionSpecLowPtSyst6080LowPtDalitz      = NULL;   
    TGraphErrors* graphRatioLowPtChargedPionsDalitz6080 = CalculateRatioBetweenSpectraWithDifferentBinning(histoYieldPi0DalitzPbPb6080Copy, graphYieldSysPi0DalitzPbPb6080Copy, 
                                                                                                           histoChargedPionSpecLowPtStat6080, histoChargedPionSpecLowPtSyst6080,  
                                                                                                           kTRUE,  kTRUE, 
                                                                                                           &graphYieldDalitzStatPi06080RebinnedLowPtDalitz, &graphYieldDalitzSysPi06080RebinnedLowPtDalitz, 
                                                                                                           &graphChargedPionSpecLowPtStat6080LowPtDalitz, &graphChargedPionSpecLowPtSyst6080LowPtDalitz);

    
    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 combined  PbPb **********************************
    // ***************************************************************************************************************
    TCanvas * canvas6PartCompChargedPions = new TCanvas("canvas6PartCompChargedPions","",10,10,1834,1000);  // gives the page size        
    canvas6PartCompChargedPions->cd();
    DrawGammaCanvasSettings( canvas6PartCompChargedPions, 0.13, 0.0, 0.02, 0.09);
    
    TPad* pad6PartCompChargedPions1 = new TPad("pad6PartCompChargedPions1", "", 0., 0.52, 0.35, 1.,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedPions1, 0.12, 0.0, 0.02, 0.);
    pad6PartCompChargedPions1->Draw();
    TPad* pad6PartCompChargedPions2 = new TPad("pad6PartCompChargedPions2", "", 0., 0., 0.35, 0.52,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedPions2, 0.12, 0.0, 0., 0.12);
    pad6PartCompChargedPions2->Draw();
    
    TPad* pad6PartCompChargedPions3 = new TPad("pad6PartCompChargedPions3", "", 0.35, 0.52, 0.68, 1.,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedPions3, 0.0, 0.0, 0.02, 0.);
    pad6PartCompChargedPions3->Draw();
    TPad* pad6PartCompChargedPions4 = new TPad("pad6PartCompChargedPions4", "", 0.35, 0., 0.68, 0.52,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedPions4, 0.0, 0.0, 0., 0.12);
    pad6PartCompChargedPions4->Draw();

    TPad* pad6PartCompChargedPions5 = new TPad("pad6PartCompChargedPions5", "", 0.68, 0.52, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedPions5, 0.0, 0.02, 0.02, 0.);
    pad6PartCompChargedPions5->Draw();
    TPad* pad6PartCompChargedPions6 = new TPad("pad6PartCompChargedPions6", "", 0.68, 0., 1., 0.52,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedPions6, 0.0, 0.02, 0., 0.12);
    pad6PartCompChargedPions6->Draw();

    TH2F * histo2DCompCombinedRatio2;
    histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.3,40.,1000,0.2,4.    );
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
    histo2DCompCombinedRatio2->GetXaxis()->SetLabelOffset(-0.015);
    SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
    
    TH2F* histo2DCompCombinedRatio;
    histo2DCompCombinedRatio = new TH2F("histo2DCompCombinedRatio","histo2DCompCombinedRatio",1000,0.3,40.,1000,0.2,4.    );
    histo2DCompCombinedRatio->GetXaxis()->SetLabelOffset(-0.015);
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.05,20.);
    SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 

    pad6PartCompChargedPions1->cd();
    pad6PartCompChargedPions1->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsComb0005, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsComb0005, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioHighPtChargedPionsComb0005->Draw("E1psame");
        graphRatioLowPtChargedPionsComb0005->Draw("E1psame");
        TLatex *labelPi0CompChargedPionsPbPb0005 = new TLatex(0.15,0.9,collisionSystemCent0.Data());
        SetStyleTLatex( labelPi0CompChargedPionsPbPb0005, 0.05,4);
        labelPi0CompChargedPionsPbPb0005->Draw();
                
        TLegend* legendPi0CompChargedPionsPbPb0005 = new TLegend(0.18,0.82,0.9,0.88);
        legendPi0CompChargedPionsPbPb0005->SetFillColor(0);
        legendPi0CompChargedPionsPbPb0005->SetLineColor(0);
        legendPi0CompChargedPionsPbPb0005->SetNColumns(2);
        legendPi0CompChargedPionsPbPb0005->SetTextSize(0.045);
        legendPi0CompChargedPionsPbPb0005->AddEntry(graphRatioLowPtChargedPionsComb0005,"#pi^{0}/#pi^{#pm} low #it{p}_{T}","p");
        legendPi0CompChargedPionsPbPb0005->AddEntry(graphRatioHighPtChargedPionsComb0005,"#pi^{0}/#pi^{#pm} high #it{p}_{T}","p");
        legendPi0CompChargedPionsPbPb0005->Draw();
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);
    
    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedPions1->Update();
    pad6PartCompChargedPions2->cd();
    pad6PartCompChargedPions2->SetLogx();
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
     histo2DCompCombinedRatio2->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsComb2040, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsComb2040, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioHighPtChargedPionsComb2040->Draw("E1psame");
        graphRatioLowPtChargedPionsComb2040->Draw("E1psame");

        TLatex *labelPi0CompChargedPionsPbPb2040 = new TLatex(0.15,0.93,collisionSystemSemiCent.Data());
        SetStyleTLatex( labelPi0CompChargedPionsPbPb2040, 0.05,4);
        labelPi0CompChargedPionsPbPb2040->Draw();
        DrawGammaLines(0., 19.5 , 1, 1 ,1, kGray, 2);    

    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedPions2->Update();
    pad6PartCompChargedPions3->cd();
    pad6PartCompChargedPions3->SetLogx();
    histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.25,20.);
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
            
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsComb0510, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsComb0510, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioHighPtChargedPionsComb0510->Draw("E1psame");
        graphRatioLowPtChargedPionsComb0510->Draw("E1psame");
        
        TLatex *labelPi0CompChargedPionsPbPb0510 = new TLatex(0.03,0.9,collisionSystemCent1.Data());
        SetStyleTLatex( labelPi0CompChargedPionsPbPb0510, 0.05,4);
        labelPi0CompChargedPionsPbPb0510->Draw(); 
        DrawGammaLines(0., 19.5 , 1, 1 ,1, kGray, 2);
    
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedPions3->Update();
    pad6PartCompChargedPions4->cd();
    pad6PartCompChargedPions4->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsComb4060, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsComb4060, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioHighPtChargedPionsComb4060->Draw("E1psame");
        graphRatioLowPtChargedPionsComb4060->Draw("E1psame");
        
        TLatex *labelPi0CompChargedPionsPbPb4060 = new TLatex(0.03,0.93,collisionSystemSemiPer.Data());
        SetStyleTLatex( labelPi0CompChargedPionsPbPb4060, 0.047,4);
        labelPi0CompChargedPionsPbPb4060->Draw();        
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedPions4->Update();
    pad6PartCompChargedPions5->cd();
    pad6PartCompChargedPions5->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsComb1020, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsComb1020, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioHighPtChargedPionsComb1020->Draw("E1psame");
        graphRatioLowPtChargedPionsComb1020->Draw("E1psame");
    
        TLatex *labelPi0CompChargedPionsPbPb1020 = new TLatex(0.03,0.9,collisionSystemCent2.Data());
        SetStyleTLatex( labelPi0CompChargedPionsPbPb1020, 0.05,4);
        labelPi0CompChargedPionsPbPb1020->Draw();        
        DrawGammaLines(0., 19.5 , 1, 1 ,1, kGray, 2);
        
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedPions5->Update();
    pad6PartCompChargedPions6->cd();
    pad6PartCompChargedPions6->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsComb6080, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsComb6080, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioHighPtChargedPionsComb6080->Draw("E1psame");
        graphRatioLowPtChargedPionsComb6080->Draw("E1psame");
        
        TLatex *labelPi0CompChargedPionsPbPb6080 = new TLatex(0.04,0.93,collisionSystemPer.Data());
        SetStyleTLatex( labelPi0CompChargedPionsPbPb6080, 0.047,4);
        labelPi0CompChargedPionsPbPb6080->Draw();    
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedPions6->Update();

    canvas6PartCompChargedPions->Update();    
    canvas6PartCompChargedPions->SaveAs(Form("%s/ComparisonChargedToNeutralCombined_6Parted_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    delete pad6PartCompChargedPions1;    
    delete pad6PartCompChargedPions2;    
    delete pad6PartCompChargedPions3;    
    delete pad6PartCompChargedPions4;    
    delete canvas6PartCompChargedPions;

    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 PCM, PHOS PbPb **********************************
    // ***************************************************************************************************************
    TCanvas * canvas6PartCompChargedIndPions = new TCanvas("canvas6PartCompChargedIndPions","",10,10,1834,1000);  // gives the page size        
    canvas6PartCompChargedIndPions->cd();
    DrawGammaCanvasSettings( canvas6PartCompChargedIndPions, 0.13, 0.0, 0.02, 0.09);
    
    TPad* pad6PartCompChargedIndPions1 = new TPad("pad6PartCompChargedIndPions1", "", 0., 0.52, 0.35, 1.,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedIndPions1, 0.12, 0.0, 0.02, 0.);
    pad6PartCompChargedIndPions1->Draw();
    TPad* pad6PartCompChargedIndPions2 = new TPad("pad6PartCompChargedIndPions2", "", 0., 0., 0.35, 0.52,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedIndPions2, 0.12, 0.0, 0., 0.12);
    pad6PartCompChargedIndPions2->Draw();
    
    TPad* pad6PartCompChargedIndPions3 = new TPad("pad6PartCompChargedIndPions3", "", 0.35, 0.52, 0.68, 1.,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedIndPions3, 0.0, 0.0, 0.02, 0.);
    pad6PartCompChargedIndPions3->Draw();
    TPad* pad6PartCompChargedIndPions4 = new TPad("pad6PartCompChargedIndPions4", "", 0.35, 0., 0.68, 0.52,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedIndPions4, 0.0, 0.0, 0., 0.12);
    pad6PartCompChargedIndPions4->Draw();

    TPad* pad6PartCompChargedIndPions5 = new TPad("pad6PartCompChargedIndPions5", "", 0.68, 0.52, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedIndPions5, 0.0, 0.02, 0.02, 0.);
    pad6PartCompChargedIndPions5->Draw();
    TPad* pad6PartCompChargedIndPions6 = new TPad("pad6PartCompChargedIndPions6", "", 0.68, 0., 1., 0.52,-1, -1, -2);
    DrawGammaPadSettings( pad6PartCompChargedIndPions6, 0.0, 0.02, 0., 0.12);
    pad6PartCompChargedIndPions6->Draw();

    pad6PartCompChargedIndPions1->cd();
    pad6PartCompChargedIndPions1->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0005, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0005, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM0005->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM0005->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS0005, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS0005, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS0005->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS0005->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb0005->Draw();
                
        TLegend* legendPi0CompChargedIndPionsPbPb0005 = new TLegend(0.18,0.76,0.9,0.88);
        legendPi0CompChargedIndPionsPbPb0005->SetFillColor(0);
        legendPi0CompChargedIndPionsPbPb0005->SetLineColor(0);
        legendPi0CompChargedIndPionsPbPb0005->SetNColumns(2);
        legendPi0CompChargedIndPionsPbPb0005->SetTextSize(0.045);
        legendPi0CompChargedIndPionsPbPb0005->AddEntry(graphRatioLowPtChargedPionsPCM0005,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompChargedIndPionsPbPb0005->AddEntry(graphRatioHighPtChargedPionsPCM0005,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
        legendPi0CompChargedIndPionsPbPb0005->AddEntry(graphRatioLowPtChargedPionsPHOS0005,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PHOS)","p");
        legendPi0CompChargedIndPionsPbPb0005->AddEntry(graphRatioHighPtChargedPionsPHOS0005,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PHOS)","p");
        legendPi0CompChargedIndPionsPbPb0005->Draw();
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);
    
    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedIndPions1->Update();
    
    pad6PartCompChargedIndPions2->cd();
    pad6PartCompChargedIndPions2->SetLogx();
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
     histo2DCompCombinedRatio2->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM2040, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM2040, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM2040->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM2040->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS2040, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS2040, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS2040->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS2040->Draw("E1psame");

        
        labelPi0CompChargedPionsPbPb2040->Draw();
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray,2 );    

    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedIndPions2->Update();
    pad6PartCompChargedIndPions3->cd();
    pad6PartCompChargedIndPions3->SetLogx();
    histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.25,15.);
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
            
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0510, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0510, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM0510->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM0510->Draw("E1psame");
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS0510, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS0510, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS0510->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS0510->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb0510->Draw(); 
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);
    
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedIndPions3->Update();
    pad6PartCompChargedIndPions4->cd();
    pad6PartCompChargedIndPions4->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM4060, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM4060, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM4060->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM4060->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS4060, markerStylePHOSHighPt, markerSizeComparison,  colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS4060, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS4060->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS4060->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb4060->Draw();        
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedIndPions4->Update();
    pad6PartCompChargedIndPions5->cd();
    pad6PartCompChargedIndPions5->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM1020, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM1020, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM1020->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM1020->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS1020, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS1020, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS1020->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS1020->Draw("E1psame");
    
        labelPi0CompChargedPionsPbPb1020->Draw();        
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);
        
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedIndPions5->Update();
    pad6PartCompChargedIndPions6->cd();
    pad6PartCompChargedIndPions6->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM6080, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM6080, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM6080->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM6080->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS6080, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS6080, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS6080->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS6080->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb6080->Draw();    
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedIndPions6->Update();

    canvas6PartCompChargedIndPions->Update();    
    canvas6PartCompChargedIndPions->SaveAs(Form("%s/ComparisonChargedToNeutralInd_6Parted_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
   
   
    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 PCM only ****************************************
    // ***************************************************************************************************************
    pad6PartCompChargedIndPions1->cd();
    pad6PartCompChargedIndPions1->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0005, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0005, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM0005->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM0005->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb0005->Draw();
                
        TLegend* legendPi0CompChargedOnlyPCMPionsPbPb0005 = new TLegend(0.18,0.76,0.9,0.88);
        legendPi0CompChargedOnlyPCMPionsPbPb0005->SetFillColor(0);
        legendPi0CompChargedOnlyPCMPionsPbPb0005->SetLineColor(0);
        legendPi0CompChargedOnlyPCMPionsPbPb0005->SetNColumns(2);
        legendPi0CompChargedOnlyPCMPionsPbPb0005->SetTextSize(0.045);
        legendPi0CompChargedOnlyPCMPionsPbPb0005->AddEntry(graphRatioLowPtChargedPionsPCM0005,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompChargedOnlyPCMPionsPbPb0005->AddEntry(graphRatioHighPtChargedPionsPCM0005,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
        legendPi0CompChargedOnlyPCMPionsPbPb0005->Draw();
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray, 2);
    
    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedIndPions1->Update();
    
    pad6PartCompChargedIndPions2->cd();
    pad6PartCompChargedIndPions2->SetLogx();
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM2040, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM2040, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM2040->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM2040->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb2040->Draw();
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2); 

    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedIndPions2->Update();
    pad6PartCompChargedIndPions3->cd();
    pad6PartCompChargedIndPions3->SetLogx();
    histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.25,15.);
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
            
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0510, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0510, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM0510->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM0510->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb0510->Draw(); 
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
    
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedIndPions3->Update();
    pad6PartCompChargedIndPions4->cd();
    pad6PartCompChargedIndPions4->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM4060, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM4060, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM4060->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM4060->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb4060->Draw();    
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedIndPions4->Update();
    pad6PartCompChargedIndPions5->cd();
    pad6PartCompChargedIndPions5->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM1020, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM1020, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM1020->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM1020->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb1020->Draw();    
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
        
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedIndPions5->Update();
    pad6PartCompChargedIndPions6->cd();
    pad6PartCompChargedIndPions6->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM6080, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM6080, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM6080->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM6080->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb6080->Draw(); 
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedIndPions6->Update();

    canvas6PartCompChargedIndPions->Update(); 
    canvas6PartCompChargedIndPions->SaveAs(Form("%s/ComparisonChargedToNeutralOnlyPCM_6Parted_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));


    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 PCM, Dalitz PbPb ********************************
    // ***************************************************************************************************************
    pad6PartCompChargedIndPions1->cd();
    pad6PartCompChargedIndPions1->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0005, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0005, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM0005->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM0005->Draw("E1psame");
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsDalitz2040, markerStyleDalitzHighPt, markerSizeComparison, colorDalitzHighPt , colorDalitzHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsDalitz2040, markerStyleDalitzLowPt, markerSizeComparison,  colorDalitzLowPt, colorDalitzLowPt);
        
        labelPi0CompChargedPionsPbPb0005->Draw();
                
        TLegend* legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005 = new TLegend(0.18,0.76,0.9,0.88);
        legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005->SetFillColor(0);
        legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005->SetLineColor(0);
        legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005->SetNColumns(2);
        legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005->SetTextSize(0.045);
        legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005->AddEntry(graphRatioLowPtChargedPionsPCM0005,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005->AddEntry(graphRatioHighPtChargedPionsPCM0005,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
        legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005->AddEntry(graphRatioLowPtChargedPionsDalitz2040,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (Dalitz)","p");
        legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005->AddEntry(graphRatioHighPtChargedPionsDalitz2040,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (Dalitz)","p");
        legendPi0CompChargedOnlyPCMAndDalitzPionsPbPb0005->Draw();
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
    
    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedIndPions1->Update();
    
    pad6PartCompChargedIndPions2->cd();
    pad6PartCompChargedIndPions2->SetLogx();
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
    

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM2040, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM2040, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM2040->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM2040->Draw("E1psame");

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsDalitz2040, markerStyleDalitzHighPt, markerSizeComparison, colorDalitzHighPt , colorDalitzHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsDalitz2040, markerStyleDalitzLowPt, markerSizeComparison,  colorDalitzLowPt, colorDalitzLowPt);
        graphRatioHighPtChargedPionsDalitz2040->Draw("E1psame");
        graphRatioLowPtChargedPionsDalitz2040->Draw("E1psame");

        
        labelPi0CompChargedPionsPbPb2040->Draw();
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2); 

    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedIndPions2->Update();
    pad6PartCompChargedIndPions3->cd();
    pad6PartCompChargedIndPions3->SetLogx();
    histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.25,15.);
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
            
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0510, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0510, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM0510->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM0510->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb0510->Draw(); 
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
    
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedIndPions3->Update();
    pad6PartCompChargedIndPions4->cd();
    pad6PartCompChargedIndPions4->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM4060, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM4060, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM4060->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM4060->Draw("E1psame");

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsDalitz4060, markerStyleDalitzHighPt, markerSizeComparison, colorDalitzHighPt , colorDalitzHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsDalitz4060, markerStyleDalitzLowPt, markerSizeComparison,  colorDalitzLowPt, colorDalitzLowPt);
        graphRatioHighPtChargedPionsDalitz4060->Draw("E1psame");
        graphRatioLowPtChargedPionsDalitz4060->Draw("E1psame");

        
        labelPi0CompChargedPionsPbPb4060->Draw();    
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedIndPions4->Update();
    pad6PartCompChargedIndPions5->cd();
    pad6PartCompChargedIndPions5->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM1020, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM1020, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM1020->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM1020->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb1020->Draw();    
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
        
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedIndPions5->Update();
    pad6PartCompChargedIndPions6->cd();
    pad6PartCompChargedIndPions6->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM6080, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM6080, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCM6080->Draw("E1psame");
        graphRatioLowPtChargedPionsPCM6080->Draw("E1psame");

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsDalitz6080, markerStyleDalitzHighPt, markerSizeComparison, colorDalitzHighPt , colorDalitzHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsDalitz6080, markerStyleDalitzLowPt, markerSizeComparison,  colorDalitzLowPt, colorDalitzLowPt);
        graphRatioHighPtChargedPionsDalitz6080->Draw("E1psame");
        graphRatioLowPtChargedPionsDalitz6080->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb6080->Draw(); 
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedIndPions6->Update();

    canvas6PartCompChargedIndPions->Update(); 
    canvas6PartCompChargedIndPions->SaveAs(Form("%s/ComparisonChargedToNeutralOnlyPCMAndDalitz_6Parted_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    
    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 PHOS only PbPb **********************************
    // ***************************************************************************************************************

    pad6PartCompChargedIndPions1->cd();
    pad6PartCompChargedIndPions1->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS0005, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS0005, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS0005->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS0005->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb0005->Draw();
                
        TLegend* legendPi0CompChargedOnlyPHOSPionsPbPb0005 = new TLegend(0.18,0.76,0.9,0.88);
        legendPi0CompChargedOnlyPHOSPionsPbPb0005->SetFillColor(0);
        legendPi0CompChargedOnlyPHOSPionsPbPb0005->SetLineColor(0);
        legendPi0CompChargedOnlyPHOSPionsPbPb0005->SetNColumns(2);
        legendPi0CompChargedOnlyPHOSPionsPbPb0005->SetTextSize(0.045);
        legendPi0CompChargedOnlyPHOSPionsPbPb0005->AddEntry(graphRatioLowPtChargedPionsPHOS0005,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PHOS)","p");
        legendPi0CompChargedOnlyPHOSPionsPbPb0005->AddEntry(graphRatioHighPtChargedPionsPHOS0005,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PHOS)","p");
        legendPi0CompChargedOnlyPHOSPionsPbPb0005->Draw();
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
    
    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedIndPions1->Update();
    
    pad6PartCompChargedIndPions2->cd();
    pad6PartCompChargedIndPions2->SetLogx();
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS2040, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS2040, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS2040->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS2040->Draw("E1psame");

        
        labelPi0CompChargedPionsPbPb2040->Draw();
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2); 

    histo2DCompCombinedRatio2->Draw("axis,same");
    pad6PartCompChargedIndPions2->Update();
    pad6PartCompChargedIndPions3->cd();
    pad6PartCompChargedIndPions3->SetLogx();
    histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.25,15.);
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
            
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS0510, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS0510, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS0510->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS0510->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb0510->Draw(); 
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
    
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedIndPions3->Update();
    pad6PartCompChargedIndPions4->cd();
    pad6PartCompChargedIndPions4->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS4060, markerStylePHOSHighPt, markerSizeComparison,  colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS4060, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS4060->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS4060->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb4060->Draw();    
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedIndPions4->Update();
    pad6PartCompChargedIndPions5->cd();
    pad6PartCompChargedIndPions5->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS1020, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS1020, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS1020->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS1020->Draw("E1psame");
    
        labelPi0CompChargedPionsPbPb1020->Draw();    
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
        
    histo2DCompCombinedRatio->Draw("axis,same");
    pad6PartCompChargedIndPions5->Update();
    pad6PartCompChargedIndPions6->cd();
    pad6PartCompChargedIndPions6->SetLogx();
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS6080, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS6080, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOS6080->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOS6080->Draw("E1psame");
        
        labelPi0CompChargedPionsPbPb6080->Draw(); 
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

    histo2DCompCombinedRatio->Draw("axis,same");
    
    pad6PartCompChargedIndPions6->Update();

    canvas6PartCompChargedIndPions->Update(); 
    canvas6PartCompChargedIndPions->SaveAs(Form("%s/ComparisonChargedToNeutralOnlyPHOS_6Parted_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    
    delete pad6PartCompChargedIndPions1;    
    delete pad6PartCompChargedIndPions2;    
    delete pad6PartCompChargedIndPions3;    
    delete pad6PartCompChargedIndPions4;    
    delete canvas6PartCompChargedIndPions;

    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 comb pp 2.76TeV *********************************
    // ***************************************************************************************************************    
    TCanvas* canvasCompYieldPPComb = new TCanvas("canvasCompYieldPPComb","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPPComb,  0.12, 0.02, 0.02, 0.12);
    
     canvasCompYieldPPComb->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsCombPP, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        graphRatioHighPtChargedPionsCombPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombPP, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioLowPtChargedPionsCombPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombPPCMS,22,markerSizeComparison, kRed , kRed);
        graphRatioLowPtChargedPionsCombPPCMS->Draw("E1psame");
        
        TLatex *labelRatioPi02760GeV = new TLatex(0.16,0.91,"pp #sqrt{#it{s}} = 2.76 TeV");
        SetStyleTLatex( labelRatioPi02760GeV, 0.06,4);
        labelRatioPi02760GeV->Draw();

        TLegend* legendPi0CompChargedPionsPP = new TLegend(0.15,0.75,0.9,0.85);
        legendPi0CompChargedPionsPP->SetFillColor(0);
        legendPi0CompChargedPionsPP->SetFillStyle(0);
        legendPi0CompChargedPionsPP->SetLineColor(0);
        legendPi0CompChargedPionsPP->SetNColumns(2);
        legendPi0CompChargedPionsPP->SetTextSize(0.045);
        legendPi0CompChargedPionsPP->SetMargin(0.12);
        legendPi0CompChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsCombPP,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP->AddEntry(graphRatioHighPtChargedPionsCombPP,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsCombPPCMS,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (#pi^{#pm} from CMS)","p");
        legendPi0CompChargedPionsPP->Draw();
    
        DrawGammaLines(0., 15.,1., 1., 0.1,kGray,2);
    
    canvasCompYieldPPComb->Update();
    canvasCompYieldPPComb->Print(Form("%s/ComparisonChargedToNeutralCombined_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 PCM, PHOS pp 2.76TeV ****************************
    // ***************************************************************************************************************    
    
    TCanvas* canvasCompYieldPPInd = new TCanvas("canvasCompYieldPPInd","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPPInd,  0.12, 0.02, 0.02, 0.12);
    
     canvasCompYieldPPInd->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMPP, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMPP, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCMPP->Draw("E1psame");
        graphRatioLowPtChargedPionsPCMPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOSPP, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOSPP, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOSPP->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOSPP->Draw("E1psame");

        labelRatioPi02760GeV->Draw();

        TLegend* legendPi0CompIndChargedPionsPP = new TLegend(0.15,0.7,0.9,0.89);
        legendPi0CompIndChargedPionsPP->SetFillColor(0);
        legendPi0CompIndChargedPionsPP->SetFillStyle(0);
        legendPi0CompIndChargedPionsPP->SetLineColor(0);
        legendPi0CompIndChargedPionsPP->SetLineStyle(0);
        legendPi0CompIndChargedPionsPP->SetNColumns(2);
        legendPi0CompIndChargedPionsPP->SetTextSize(0.045);
        legendPi0CompIndChargedPionsPP->SetMargin(0.1);
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsPCMPP,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioHighPtChargedPionsPCMPP,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsPHOSPP,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PHOS)","p");
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioHighPtChargedPionsPHOSPP,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PHOS)","p");
        legendPi0CompIndChargedPionsPP->Draw();

        legendPi0CompIndChargedPionsPP->Draw();
        DrawGammaLines(0., 15.,1., 1., 1,kGray,2);
    
    
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedToNeutralInd_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ************************** Comparison pi0/pi+-, pi0 updated comb pp 2.76TeV ***********************
    // ***************************************************************************************************************    
    canvasCompYieldPPInd->cd();
    canvasCompYieldPPInd->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();

    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsCombUpPP, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        graphRatioHighPtChargedPionsCombUpPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombUpPP, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioLowPtChargedPionsCombUpPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombUpPPCMS,22,markerSizeComparison, kRed , kRed);
        graphRatioLowPtChargedPionsCombUpPPCMS->Draw("E1psame");
        
        labelRatioPi02760GeV->Draw();
        legendPi0CompChargedPionsPP->Draw();
    
        DrawGammaLines(0., 20 , 1, 1 ,1, kGray, 2);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedToNeutralCombinedUpdated_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ************************** Comparison pi0/pi+-, pi0 updated & old comb pp 2.76TeV ***********************
    // ***************************************************************************************************************    
    canvasCompYieldPPInd->cd();
    canvasCompYieldPPInd->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();

    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsCombUpPP, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        graphRatioHighPtChargedPionsCombUpPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombUpPP, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioLowPtChargedPionsCombUpPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombUpPPCMS,22,markerSizeComparison, kRed , kRed);
        graphRatioLowPtChargedPionsCombUpPPCMS->Draw("E1psame");

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsCombPP, markerStyleCombHighPt+4, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        graphRatioHighPtChargedPionsCombPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombPP, markerStyleCombLowPt+4, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioLowPtChargedPionsCombPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombPPCMS,22+4,markerSizeComparison, kRed , kRed);
        graphRatioLowPtChargedPionsCombPPCMS->Draw("E1psame");
        
        
        labelRatioPi02760GeV->Draw();
        legendPi0CompChargedPionsPP->SetY1NDC(0.65);
        legendPi0CompChargedPionsPP->AddEntry((TObject*)0,"","");
        legendPi0CompChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsCombUpPP,"#pi^{0} up. /#pi^{#pm} low #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP->AddEntry(graphRatioHighPtChargedPionsCombUpPP,"#pi^{0} up./#pi^{#pm} high #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsCombUpPPCMS,"#pi^{0} up. /#pi^{#pm} low #it{p}_{T} (#pi^{#pm} from CMS)","p");
        legendPi0CompChargedPionsPP->Draw();
    
        DrawGammaLines(0., 20 , 1, 1 ,1, kGray, 2);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedToNeutralCombinedUpdatedCompPrev_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    
    // ***************************************************************************************************************
    // ************************** Comparison pi0/pi+-, pi0 PCM, PHOS, Dalitz, EMCAL pp 2.76TeV ***********************
    // ***************************************************************************************************************    
    
    canvasCompYieldPPInd->cd();
    canvasCompYieldPPInd->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMPP, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMPP, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCMPP->Draw("E1psame");
        graphRatioLowPtChargedPionsPCMPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOSPP, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOSPP, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOSPP->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOSPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsDalitzPP2760GeV, markerStyleDalitzHighPt, markerSizeComparison, colorDalitzHighPt , colorDalitzHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsDalitzPP2760GeV, markerStyleDalitzLowPt, markerSizeComparison,  colorDalitzLowPt, colorDalitzLowPt);
        graphRatioHighPtChargedPionsDalitzPP2760GeV->Draw("E1psame");
        graphRatioLowPtChargedPionsDalitzPP2760GeV->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsEMCALPP, markerStyleEMCALHighPt, markerSizeComparison, colorEMCALHighPt , colorEMCALHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsEMCALPP, markerStyleEMCALLowPt, markerSizeComparison,  colorEMCALLowPt ,colorEMCALLowPt );
        graphRatioHighPtChargedPionsEMCALPP->Draw("E1psame");
        graphRatioLowPtChargedPionsEMCALPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsEMCALMergedPP, markerStyleEMCALMergedHighPt, markerSizeComparison, colorEMCALMergedHighPt , colorEMCALMergedHighPt);
        graphRatioHighPtChargedPionsEMCALMergedPP->Draw("E1psame");
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMEMCALPP, markerStylePCMEMCALHighPt, markerSizeComparison, colorPCMEMCALHighPt , colorPCMEMCALHighPt);
        graphRatioHighPtChargedPionsPCMEMCALPP->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMEMCALPP, markerStylePCMEMCALLowPt, markerSizeComparison,  colorPCMEMCALLowPt ,colorPCMEMCALLowPt );
        graphRatioLowPtChargedPionsPCMEMCALPP->Draw("E1psame");
        
        labelRatioPi02760GeV->Draw();
        legendPi0CompIndChargedPionsPP->SetY1NDC(0.65);
        legendPi0CompIndChargedPionsPP->SetX1NDC(0.15);
        legendPi0CompIndChargedPionsPP->SetTextSize(0.04);
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsDalitzPP2760GeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (Dalitz)","p");
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioHighPtChargedPionsDalitzPP2760GeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (Dalitz)","p");
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsEMCALPP,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (EMCal)","p");
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioHighPtChargedPionsEMCALPP,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (EMCal)","p");
        legendPi0CompIndChargedPionsPP->AddEntry((TObject*)0,"","");
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioHighPtChargedPionsEMCALMergedPP,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (EMCal merged)","p");
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsPCMEMCALPP,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM-EMCal)","p");
        legendPi0CompIndChargedPionsPP->AddEntry(graphRatioHighPtChargedPionsPCMEMCALPP,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM-EMCal)","p");
        legendPi0CompIndChargedPionsPP->Draw();

        legendPi0CompIndChargedPionsPP->Draw();
        DrawGammaLines(0., 20 , 1, 1 ,1, kGray, 2);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedToNeutralInd_allMeasurements_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    canvasCompYieldPPInd->cd();
    canvasCompYieldPPInd->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();

        graphRatioHighPtChargedPionsPCMPP->Draw("E1psame");
        graphRatioLowPtChargedPionsPCMPP->Draw("E1psame");
        graphRatioHighPtChargedPionsPHOSPP->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOSPP->Draw("E1psame");
        graphRatioHighPtChargedPionsDalitzPP2760GeV->Draw("E1psame");
        graphRatioLowPtChargedPionsDalitzPP2760GeV->Draw("E1psame");
        graphRatioHighPtChargedPionsEMCALPP->Draw("E1psame");
        graphRatioLowPtChargedPionsEMCALPP->Draw("E1psame");
//         if (enablePCMEMCALComp2760GeV){
            if (graphRatioHighPtChargedPionsPCMEMCALPP){
                DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMEMCALPP, markerStylePCMEMCALHighPt, markerSizeComparison, colorPCMEMCALHighPt , colorPCMEMCALHighPt);
                graphRatioHighPtChargedPionsPCMEMCALPP->Draw("E1psame");
            }
            if (graphRatioLowPtChargedPionsPCMEMCALPP){
                DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMEMCALPP, markerStylePCMEMCALLowPt, markerSizeComparison,  colorPCMEMCALLowPt ,colorPCMEMCALLowPt );
                graphRatioLowPtChargedPionsPCMEMCALPP->Draw("E1psame");
            }    
//         }
        if (enablePCMPHOSComp2760GeV){
            if (graphRatioHighPtChargedPionsPCMPHOSPP){
                graphRatioHighPtChargedPionsPCMPHOSPP->RemovePoint(graphRatioHighPtChargedPionsPCMPHOSPP->GetN()-1);
                graphRatioHighPtChargedPionsPCMPHOSPP->RemovePoint(graphRatioHighPtChargedPionsPCMPHOSPP->GetN()-1);
                DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMPHOSPP, markerStylePCMPHOSHighPt, markerSizeComparison, colorPCMPHOSHighPt , colorPCMPHOSHighPt);
                graphRatioHighPtChargedPionsPCMPHOSPP->Draw("E1psame");
            }
            if (graphRatioLowPtChargedPionsPCMPHOSPP){
                graphRatioLowPtChargedPionsPCMPHOSPP->RemovePoint(0);
                DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMPHOSPP, markerStylePCMPHOSLowPt, markerSizeComparison,  colorPCMPHOSLowPt ,colorPCMPHOSLowPt );
                graphRatioLowPtChargedPionsPCMPHOSPP->Draw("E1psame");
            }    
        }
        
        labelRatioPi02760GeV->Draw();

        legendPi0CompIndChargedPionsPP->SetY1NDC(0.65);
        legendPi0CompIndChargedPionsPP->SetX1NDC(0.15);
        legendPi0CompIndChargedPionsPP->SetTextSize(0.04);
        if (enablePCMPHOSComp2760GeV && (graphRatioLowPtChargedPionsPCMPHOSPP|| graphRatioHighPtChargedPionsPCMPHOSPP) ){
            legendPi0CompIndChargedPionsPP->SetY2NDC(0.89);
            legendPi0CompIndChargedPionsPP->SetY1NDC(0.6);
            legendPi0CompIndChargedPionsPP->SetX1NDC(0.15);
            legendPi0CompIndChargedPionsPP->SetTextSize(0.035);
        }    
        if (enablePCMPHOSComp2760GeV && graphRatioLowPtChargedPionsPCMPHOSPP) 
            legendPi0CompIndChargedPionsPP->AddEntry(graphRatioLowPtChargedPionsPCMPHOSPP,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM-PHOS)","p");
        if (enablePCMPHOSComp2760GeV && graphRatioHighPtChargedPionsPCMPHOSPP) 
            legendPi0CompIndChargedPionsPP->AddEntry(graphRatioHighPtChargedPionsPCMPHOSPP,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM-PHOS)","p");
        legendPi0CompIndChargedPionsPP->Draw();

        legendPi0CompIndChargedPionsPP->Draw();
        DrawGammaLines(0., 20 , 1, 1 ,1, kGray, 2);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedToNeutralInd_reallyAllMeasurements_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    
    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 comb pp 7TeV ************************************
    // ***************************************************************************************************************    
    
    TCanvas* canvasCompYieldPP7TeVComb = new TCanvas("canvasCompYieldPP7TeVComb","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPP7TeVComb,  0.12, 0.02, 0.02, 0.12);
    
     canvasCompYieldPP7TeVComb->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsCombPP7TeV, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        graphRatioHighPtChargedPionsCombPP7TeV->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombPP7TeV, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioLowPtChargedPionsCombPP7TeV->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombPPCMS7TeV,22,markerSizeComparison, kRed , kRed);
        graphRatioLowPtChargedPionsCombPPCMS7TeV->Draw("E1psame");
        
        TLatex *labelRatioPi07TeV = new TLatex(0.16,0.9,"pp #sqrt{#it{s}} = 7 TeV");
        SetStyleTLatex( labelRatioPi07TeV, 0.06,4);
        labelRatioPi07TeV->Draw();

        TLegend* legendPi0CompChargedPionsPP7TeV = new TLegend(0.15,0.75,0.9,0.85);
        legendPi0CompChargedPionsPP7TeV->SetFillColor(0);
        legendPi0CompChargedPionsPP7TeV->SetLineColor(0);
        legendPi0CompChargedPionsPP7TeV->SetNColumns(2);
        legendPi0CompChargedPionsPP7TeV->SetTextSize(0.045);
        legendPi0CompChargedPionsPP7TeV->AddEntry(graphRatioLowPtChargedPionsCombPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP7TeV->AddEntry(graphRatioHighPtChargedPionsCombPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP7TeV->AddEntry(graphRatioLowPtChargedPionsCombPPCMS7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (#pi^{#pm} from CMS)","p");
        legendPi0CompChargedPionsPP7TeV->Draw();
        
        DrawGammaLines(0., 15.,1., 1.,0.1,kGray,2);
        
        
        canvasCompYieldPP7TeVComb->Update();
        canvasCompYieldPP7TeVComb->Print(Form("%s/ComparisonChargedToNeutralCombined_PP7TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 comb without CMS pp 7TeV ************************
    // ***************************************************************************************************************    
        
    canvasCompYieldPP7TeVComb->cd();
    canvasCompYieldPP7TeVComb->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
    
        
        graphRatioHighPtChargedPionsCombPP7TeV->Draw("E1psame");
        
        graphRatioLowPtChargedPionsCombPP7TeV->Draw("E1psame");
        
        labelRatioPi07TeV->Draw();

        TLegend* legendPi0CompChargedPionsPP7TeVWOCMS = new TLegend(0.15,0.75,0.9,0.85);
        legendPi0CompChargedPionsPP7TeVWOCMS->SetFillColor(0);
        legendPi0CompChargedPionsPP7TeVWOCMS->SetLineColor(0);
        legendPi0CompChargedPionsPP7TeVWOCMS->SetNColumns(2);
        legendPi0CompChargedPionsPP7TeVWOCMS->SetTextSize(0.045);
        legendPi0CompChargedPionsPP7TeVWOCMS->AddEntry(graphRatioLowPtChargedPionsCombPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP7TeVWOCMS->AddEntry(graphRatioHighPtChargedPionsCombPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP7TeVWOCMS->Draw();
    
        DrawGammaLines(0., 15.,1., 1.,1,kGray,2);
    
    
    canvasCompYieldPP7TeVComb->Update();
    canvasCompYieldPP7TeVComb->Print(Form("%s/ComparisonChargedToNeutralCombinedWOCMS_PP7TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ************************************ Comparison pi0/pi+-, pi0 PCM, PHOS pp 7TeV *******************************
    // ***************************************************************************************************************    
            
    TCanvas* canvasCompYieldPP7TeVInd = new TCanvas("canvasCompYieldPP7TeVInd","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPP7TeVInd,  0.12, 0.02, 0.02, 0.12);
    
     canvasCompYieldPP7TeVInd->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,30.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMPP7TeV, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMPP7TeV, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCMPP7TeV->Draw("E1psame");
        graphRatioLowPtChargedPionsPCMPP7TeV->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOSPP7TeV, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOSPP7TeV, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOSPP7TeV->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOSPP7TeV->Draw("E1psame");

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsEMCalPP7TeV, markerStyleEMCALHighPt, markerSizeComparison, colorEMCALHighPt , colorEMCALHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsEMCalPP7TeV, markerStyleEMCALLowPt, markerSizeComparison, colorEMCALLowPt , colorEMCALLowPt);
        graphRatioHighPtChargedPionsEMCalPP7TeV->Draw("E1psame");
        graphRatioLowPtChargedPionsEMCalPP7TeV->Draw("E1psame");

        labelRatioPi07TeV->Draw();

        TLegend* legendPi0CompIndChargedPionsPP7TeV = new TLegend(0.15,0.73,0.9,0.85);
        legendPi0CompIndChargedPionsPP7TeV->SetFillColor(0);
        legendPi0CompIndChargedPionsPP7TeV->SetLineColor(0);
        legendPi0CompIndChargedPionsPP7TeV->SetNColumns(2);
        legendPi0CompIndChargedPionsPP7TeV->SetTextSize(0.045);
        legendPi0CompIndChargedPionsPP7TeV->AddEntry(graphRatioLowPtChargedPionsPCMPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeV->AddEntry(graphRatioHighPtChargedPionsPCMPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeV->AddEntry(graphRatioLowPtChargedPionsPHOSPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PHOS)","p");
        legendPi0CompIndChargedPionsPP7TeV->AddEntry(graphRatioHighPtChargedPionsPHOSPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PHOS)","p");
        legendPi0CompIndChargedPionsPP7TeV->AddEntry(graphRatioLowPtChargedPionsEMCalPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (EMCal)","p");
        legendPi0CompIndChargedPionsPP7TeV->AddEntry(graphRatioHighPtChargedPionsEMCalPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (EMCal)","p");
        legendPi0CompIndChargedPionsPP7TeV->Draw();

        legendPi0CompIndChargedPionsPP7TeV->Draw();
        DrawGammaLines(0., 30.,1., 1.,1,kGray,2);
        
    canvasCompYieldPP7TeVInd->Update();
    canvasCompYieldPP7TeVInd->Print(Form("%s/ComparisonChargedToNeutralInd_PP7TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    // **************************************************************************************************************
    // ************************************ PCM additional pileup correction ****************************************
    // **************************************************************************************************************
    TF1* fitCorrectionFactorsHistvsPt                                 = new TF1("fitCorrectionFactorsHistvsPt","[0]/pow(x,[1])+[2]");
    fitCorrectionFactorsHistvsPt->SetParameter(0,2.9737546081);
    fitCorrectionFactorsHistvsPt->SetParameter(1,1.4795520406);
    fitCorrectionFactorsHistvsPt->SetParameter(2,2.2652589579);

    TGraphErrors* graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr   = (TGraphErrors*)graphRatioHighPtChargedPionsPCMPP7TeV->Clone("graphRatioHighPtChargedPionsPCMPP7TeV");
    TGraphErrors* graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr    = (TGraphErrors*)graphRatioLowPtChargedPionsPCMPP7TeV->Clone("graphRatioLowPtChargedPionsPCMPP7TeV");
    Double_t* yLow                                                  = graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr->GetY();
    Double_t* xLow                                                  = graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr->GetX();
    Double_t* yHigh                                                 = graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr->GetY();
    Double_t* xHigh                                                 = graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr->GetX();
    for (Int_t i = 0; i < graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr->GetN(); i++){
        cout << xLow[i] << "\t" << 100-fitCorrectionFactorsHistvsPt->Eval(xLow[i]) << endl;
        cout << yLow[i] << "\t" << yLow[i]*(100-fitCorrectionFactorsHistvsPt->Eval(xLow[i]))/100 << "\t" << yLow[i] << endl;
        yLow[i]     = yLow[i]*(100-fitCorrectionFactorsHistvsPt->Eval(xLow[i]))/100;
    }   
    cout << "*************************" << endl;
    for (Int_t i = 0; i < graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr->GetN(); i++){
        cout << xHigh[i] << "\t" <<100-fitCorrectionFactorsHistvsPt->Eval(xHigh[i]) << endl;
        cout << yHigh[i] << "\t" << yHigh[i]*(100-fitCorrectionFactorsHistvsPt->Eval(xHigh[i]))/100 << "\t" << yHigh[i] << endl;
        yHigh[i]     = yHigh[i]*(100-fitCorrectionFactorsHistvsPt->Eval(xHigh[i]))/100;
    }   
   
    // ***************************************************************************************************************
    // ************************** Comparison pi0/pi+-, pi0 PCM only, + add pileup corr pp 7TeV ***********************
    // ***************************************************************************************************************    

    canvasCompYieldPP7TeVInd->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMPP7TeV, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMPP7TeV, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCMPP7TeV->Draw("E1psame");
        graphRatioLowPtChargedPionsPCMPP7TeV->Draw("E1psame");
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr,25,markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr,24,markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr->Draw("E1psame");
        graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr->Draw("E1psame");
        
        
        TLatex *labelPileupCorr = new TLatex(0.16,0.82,"Open symbols out of bunch pileup corrected");
        SetStyleTLatex( labelPileupCorr, 0.05,4);
        labelPileupCorr->Draw();

        labelRatioPi07TeV->Draw();

        TLegend* legendPi0CompIndChargedPionsPP7TeVPileupCorr = new TLegend(0.13,0.65,0.75,0.8);
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->SetFillColor(0);
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->SetLineColor(0);
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->SetNColumns(2);
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->SetTextSize(0.04);
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->AddEntry(graphRatioLowPtChargedPionsPCMPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->AddEntry(graphRatioHighPtChargedPionsPCMPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->AddEntry(graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->AddEntry(graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->AddEntry((TObject*)0,"pileup corr","");
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->AddEntry((TObject*)0,"pileup corr","");
        legendPi0CompIndChargedPionsPP7TeVPileupCorr->Draw();

        legendPi0CompIndChargedPionsPP7TeVPileupCorr->Draw();
        DrawGammaLines(0., 20.,1., 1.,1,kGray,2);
    
    
    canvasCompYieldPP7TeVInd->Update();
    canvasCompYieldPP7TeVInd->Print(Form("%s/ComparisonCharged7TeVPCMPileupCorrected_PP7TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    delete canvasCompYieldPP7TeVInd;

    // ***************************************************************************************************************
    // ************************** Comparison pi0/pi+-, pi0 PCM add pileup corr, Dalitz pp 7TeV ***********************
    // ***************************************************************************************************************    
    
    TCanvas* canvasCompYieldPP7TeVIndDalitz = new TCanvas("canvasCompYieldPP7TeVIndDalitz","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPP7TeVIndDalitz,  0.12, 0.02, 0.02, 0.12);
    canvasCompYieldPP7TeVIndDalitz->cd();
    canvasCompYieldPP7TeVIndDalitz->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->Draw();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr->Draw("E1psame");
        graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr->Draw("E1psame");
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsDalitzPP7TeV, markerStyleDalitzHighPt, markerSizeComparison, colorDalitzHighPt , colorDalitzHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsDalitzPP7TeV, markerStyleDalitzLowPt, markerSizeComparison,  colorDalitzLowPt, colorDalitzLowPt);
        graphRatioHighPtChargedPionsDalitzPP7TeV->Draw("E1psame");
        graphRatioLowPtChargedPionsDalitzPP7TeV->Draw("E1psame");
        
        labelRatioPi07TeV->Draw();

        TLegend* legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz = new TLegend(0.13,0.70,0.75,0.85);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->SetFillColor(0);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->SetLineColor(0);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->SetNColumns(2);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->SetTextSize(0.04);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->AddEntry(graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->AddEntry(graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->AddEntry((TObject*)0,"pileup corr","");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->AddEntry((TObject*)0,"pileup corr","");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->AddEntry(graphRatioLowPtChargedPionsDalitzPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (Dalitz)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->AddEntry(graphRatioHighPtChargedPionsDalitzPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (Dalitz)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitz->Draw();

         DrawGammaLines(0., 20 , 1, 1 ,1,kGray, 2);    
    
    canvasCompYieldPP7TeVIndDalitz->Update();
    canvasCompYieldPP7TeVIndDalitz->Print(Form("%s/ComparisonCharged7TeVPCMPileupCorrectedDalitz_PP7TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    
    // ***************************************************************************************************************
    // ***************** Comparison pi0/pi+-, pi0 PCM add pileup corr, Dalitz, EMCAL, PHOS pp 7TeV *******************
    // ***************************************************************************************************************    

    TCanvas* canvasCompYieldPP7TeVIndDalitzAll = new TCanvas("canvasCompYieldPP7TeVIndDalitzAll","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPP7TeVIndDalitzAll,  0.12, 0.02, 0.02, 0.12);
    canvasCompYieldPP7TeVIndDalitzAll->cd();
    canvasCompYieldPP7TeVIndDalitzAll->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,30.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->Draw();

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr->Draw("E1psame");
        graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr->Draw("E1psame");

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOSPP7TeV, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOSPP7TeV, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioHighPtChargedPionsPHOSPP7TeV->Draw("E1psame");
        graphRatioLowPtChargedPionsPHOSPP7TeV->Draw("E1psame");
        
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsDalitzPP7TeV, markerStyleDalitzHighPt, markerSizeComparison, colorDalitzHighPt , colorDalitzHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsDalitzPP7TeV, markerStyleDalitzLowPt, markerSizeComparison,  colorDalitzLowPt, colorDalitzLowPt);
        graphRatioHighPtChargedPionsDalitzPP7TeV->Draw("E1psame");
        graphRatioLowPtChargedPionsDalitzPP7TeV->Draw("E1psame");
      

        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsEMCalPP7TeV, markerStyleEMCALHighPt, markerSizeComparison, colorEMCALHighPt , colorEMCALHighPt);
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsEMCalPP7TeV, markerStyleEMCALLowPt, markerSizeComparison, colorEMCALLowPt , colorEMCALLowPt);
        graphRatioHighPtChargedPionsEMCalPP7TeV->Draw("E1psame");
        graphRatioLowPtChargedPionsEMCalPP7TeV->Draw("E1psame");

        labelRatioPi07TeV->Draw();

        TLegend* legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll = new TLegend(0.13,0.69,0.75,0.89);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->SetFillColor(0);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->SetFillStyle(0);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->SetLineColor(0);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->SetNColumns(2);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->SetTextSize(0.04);
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry(graphRatioLowPtChargedPionsPCMPP7TeVPileupCorr,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry(graphRatioHighPtChargedPionsPCMPP7TeVPileupCorr,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry((TObject*)0,"pileup corr","");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry((TObject*)0,"pileup corr","");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry(graphRatioLowPtChargedPionsDalitzPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (Dalitz)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry(graphRatioHighPtChargedPionsDalitzPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (Dalitz)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry(graphRatioLowPtChargedPionsPHOSPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PHOS)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry(graphRatioHighPtChargedPionsPHOSPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PHOS)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry(graphRatioLowPtChargedPionsEMCalPP7TeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (EMCal)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->AddEntry(graphRatioHighPtChargedPionsEMCalPP7TeV,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (EMCal)","p");
        legendPi0CompIndChargedPionsPP7TeVPileupCorrDalitzAll->Draw();

        DrawGammaLines(0., 30 , 1, 1 ,1, kGray, 2);
   
    canvasCompYieldPP7TeVIndDalitzAll->Update();
    canvasCompYieldPP7TeVIndDalitzAll->Print(Form("%s/ComparisonCharged7TeVPCMPileupCorrectedDalitzAll_PP7TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    
    // ***************************************************************************************************************
    // ************************************* Comparison pi0/pi+-, pi0 comb pp 0.9TeV *********************************
    // ***************************************************************************************************************    
   
    TCanvas* canvasCompYieldPP900GeVComb = new TCanvas("canvasCompYieldPP900GeVComb","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPP900GeVComb,  0.12, 0.02, 0.02, 0.12);
    
     canvasCompYieldPP900GeVComb->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();
    
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombPP900GeV, markerStyleCombLowPt, markerSizeComparison, kBlue , kBlue);
        graphRatioLowPtChargedPionsCombPP900GeV->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsCombPPCMS900GeV,22,markerSizeComparison, kRed , kRed);
        graphRatioLowPtChargedPionsCombPPCMS900GeV->Draw("E1psame");
        
        
        TLatex *labelRatioPi0900GeV = new TLatex(0.16,0.9,"pp #sqrt{#it{s}} = 0.9 TeV");
        SetStyleTLatex( labelRatioPi0900GeV, 0.06,4);
        labelRatioPi0900GeV->Draw();

        TLegend* legendPi0CompChargedPionsPP900GeV = new TLegend(0.15,0.75,0.9,0.85);
        legendPi0CompChargedPionsPP900GeV->SetFillColor(0);
        legendPi0CompChargedPionsPP900GeV->SetLineColor(0);
        legendPi0CompChargedPionsPP900GeV->SetNColumns(2);
        legendPi0CompChargedPionsPP900GeV->SetTextSize(0.045);
        legendPi0CompChargedPionsPP900GeV->AddEntry(graphRatioLowPtChargedPionsCombPP900GeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP900GeV->AddEntry(graphRatioLowPtChargedPionsCombPPCMS900GeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (#pi^{#pm} from CMS)","p");
        legendPi0CompChargedPionsPP900GeV->Draw();
    
        DrawGammaLines(0., 15.,1., 1., 1, kGray,2);
    
    
    canvasCompYieldPP900GeVComb->Update();
    canvasCompYieldPP900GeVComb->Print(Form("%s/ComparisonChargedToNeutralCombined_PP900GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ************************************* Comparison pi0/pi+-, pi0 PCM, PHOS pp 0.9TeV ****************************
    // ***************************************************************************************************************    
    
    TCanvas* canvasCompYieldPP900GeVInd = new TCanvas("canvasCompYieldPP900GeVInd","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPP900GeVInd,  0.12, 0.02, 0.02, 0.12);
    
     canvasCompYieldPP900GeVInd->SetLogx();
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMPP900GeV, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
        graphRatioLowPtChargedPionsPCMPP900GeV->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOSPP900GeV, markerStylePHOSLowPt, markerSizeComparison,  colorPHOSLowPt, colorPHOSLowPt);
        graphRatioLowPtChargedPionsPHOSPP900GeV->Draw("E1psame");

        labelRatioPi0900GeV->Draw();

        TLegend* legendPi0CompIndChargedPionsPP900GeV = new TLegend(0.15,0.75,0.9,0.85);
        legendPi0CompIndChargedPionsPP900GeV->SetFillColor(0);
        legendPi0CompIndChargedPionsPP900GeV->SetLineColor(0);
        legendPi0CompIndChargedPionsPP900GeV->SetNColumns(2);
        legendPi0CompIndChargedPionsPP900GeV->SetTextSize(0.045);
        legendPi0CompIndChargedPionsPP900GeV->AddEntry(graphRatioLowPtChargedPionsPCMPP900GeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
        legendPi0CompIndChargedPionsPP900GeV->AddEntry(graphRatioLowPtChargedPionsPHOSPP900GeV,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PHOS)","p");
        legendPi0CompIndChargedPionsPP900GeV->Draw();

        legendPi0CompIndChargedPionsPP900GeV->Draw();
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
    
    DrawGammaLines(0., 15.,1., 1.,0.1,kGray,2);
    
    
    canvasCompYieldPP900GeVInd->Update();
    canvasCompYieldPP900GeVInd->Print(Form("%s/ComparisonChargedToNeutralInd_PP900GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    delete canvasCompYieldPP900GeVInd;
   
    
    TFile fCombResultsOutput(Form("ComparisonResults%s.root",dateForOutput.Data()),"UPDATE");
        graphRatioHighPtChargedPionsCombPP->Write("graphRatioHighPtChargedPionsCombPP2760GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsCombPP->Write("graphRatioLowPtChargedPionsCombPP2760GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsCombPPCMS->Write("graphRatioLowPtChargedPionsCombPPCMS2760GeV",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPCMPP->Write("graphRatioHighPtChargedPionsPCMPP2760GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCMPP->Write("graphRatioLowPtChargedPionsPCMPP2760GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCMPPCMS->Write("graphRatioLowPtChargedPionsPCMPPCMS2760GeV",TObject::kOverwrite);   
        graphRatioHighPtChargedPionsPHOSPP->Write("graphRatioHighPtChargedPionsPHOSPP2760GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPHOSPP->Write("graphRatioLowPtChargedPionsPHOSPP2760GeV",TObject::kOverwrite);

        graphRatioHighPtChargedPionsCombUpPP->Write("graphRatioHighPtChargedPionsCombUpPP2760GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsCombUpPP->Write("graphRatioLowPtChargedPionsCombUpPP2760GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsCombUpPPCMS->Write("graphRatioLowPtChargedPionsCombUpPPCMS2760GeV",TObject::kOverwrite);
        
        graphRatioHighPtChargedPionsCombPP7TeV->Write("graphRatioHighPtChargedPionsCombPP7TeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsCombPP7TeV->Write("graphRatioLowPtChargedPionsCombPP7TeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsCombPPCMS7TeV->Write("graphRatioLowPtChargedPionsCombPPCMS7TeV",TObject::kOverwrite);   
        graphRatioHighPtChargedPionsPCMPP7TeV->Write("graphRatioHighPtChargedPionsPCMPP7TeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCMPP7TeV->Write("graphRatioLowPtChargedPionsPCMPP7TeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCMPPCMS7TeV->Write("graphRatioLowPtChargedPionsPCMPPCMS7TeV",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPHOSPP7TeV->Write("graphRatioHighPtChargedPionsPHOSPP7TeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPHOSPP7TeV->Write("graphRatioLowPtChargedPionsPHOSPP7TeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsCombPP900GeV->Write("graphRatioLowPtChargedPionsCombPP900GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsCombPPCMS900GeV->Write("graphRatioLowPtChargedPionsCombPPCMS900GeV",TObject::kOverwrite);   
        graphRatioLowPtChargedPionsPCMPP900GeV->Write("graphRatioLowPtChargedPionsPCMPP900GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCMPPCMS900GeV->Write("graphRatioLowPtChargedPionsPCMPPCMS900GeV",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPHOSPP900GeV->Write("graphRatioLowPtChargedPionsPHOSPP900GeV",TObject::kOverwrite);
        graphRatioHighPtChargedPionsComb0005->Write("graphRatioHighPtChargedPionsComb0005",TObject::kOverwrite);
        graphRatioLowPtChargedPionsComb0005->Write("graphRatioLowPtChargedPionsComb0005",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPCM0005->Write("graphRatioHighPtChargedPionsPCM0005",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCM0005->Write("graphRatioLowPtChargedPionsPCM0005",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPHOS0005->Write("graphRatioHighPtChargedPionsPHOS0005",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPHOS0005->Write("graphRatioLowPtChargedPionsPHOS0005",TObject::kOverwrite);
        graphRatioHighPtChargedPionsComb0510->Write("graphRatioHighPtChargedPionsComb0510",TObject::kOverwrite);
        graphRatioLowPtChargedPionsComb0510->Write("graphRatioLowPtChargedPionsComb0510",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPCM0510->Write("graphRatioHighPtChargedPionsPCM0510",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCM0510->Write("graphRatioLowPtChargedPionsPCM0510",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPHOS0510->Write("graphRatioHighPtChargedPionsPHOS0510",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPHOS0510->Write("graphRatioLowPtChargedPionsPHOS0510",TObject::kOverwrite);
        graphRatioHighPtChargedPionsComb1020->Write("graphRatioHighPtChargedPionsComb1020",TObject::kOverwrite);
        graphRatioLowPtChargedPionsComb1020->Write("graphRatioLowPtChargedPionsComb1020",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPCM1020->Write("graphRatioHighPtChargedPionsPCM1020",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCM1020->Write("graphRatioLowPtChargedPionsPCM1020",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPHOS1020->Write("graphRatioHighPtChargedPionsPHOS1020",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPHOS1020->Write("graphRatioLowPtChargedPionsPHOS1020",TObject::kOverwrite);
        graphRatioHighPtChargedPionsComb2040->Write("graphRatioHighPtChargedPionsComb2040",TObject::kOverwrite);
        graphRatioLowPtChargedPionsComb2040->Write("graphRatioLowPtChargedPionsComb2040",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPCM2040->Write("graphRatioHighPtChargedPionsPCM2040",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCM2040->Write("graphRatioLowPtChargedPionsPCM2040",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPHOS2040->Write("graphRatioHighPtChargedPionsPHOS2040",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPHOS2040->Write("graphRatioLowPtChargedPionsPHOS2040",TObject::kOverwrite);
        graphRatioHighPtChargedPionsComb4060->Write("graphRatioHighPtChargedPionsComb4060",TObject::kOverwrite);
        graphRatioLowPtChargedPionsComb4060->Write("graphRatioLowPtChargedPionsComb4060",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPCM4060->Write("graphRatioHighPtChargedPionsPCM4060",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCM4060->Write("graphRatioLowPtChargedPionsPCM4060",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPHOS4060->Write("graphRatioHighPtChargedPionsPHOS4060",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPHOS4060->Write("graphRatioLowPtChargedPionsPHOS4060",TObject::kOverwrite);
        graphRatioHighPtChargedPionsComb6080->Write("graphRatioHighPtChargedPionsComb6080",TObject::kOverwrite);
        graphRatioLowPtChargedPionsComb6080->Write("graphRatioLowPtChargedPionsComb6080",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPCM6080->Write("graphRatioHighPtChargedPionsPCM6080",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPCM6080->Write("graphRatioLowPtChargedPionsPCM6080",TObject::kOverwrite);
        graphRatioHighPtChargedPionsPHOS6080->Write("graphRatioHighPtChargedPionsPHOS6080",TObject::kOverwrite);
        graphRatioLowPtChargedPionsPHOS6080->Write("graphRatioLowPtChargedPionsPHOS6080",TObject::kOverwrite);      
    fCombResultsOutput.Close();
    
    TFile fCombResultsOutputRebinned(Form("InputSpectraRebinnedToFitCharged_%s.root",dateForOutput.Data()),"UPDATE");
        graphYieldCombStatPi00005RebinnedHighPtComb->Write("graphYieldCombStatPi00005_RebinnedHighPtComb");
        graphYieldCombSysPi00005RebinnedHighPtComb->Write("graphYieldCombSystPi00005_RebinnedHighPtComb");
        graphYieldCombStatPi00005RebinnedLowPtComb->Write("graphYieldCombStatPi00005_RebinnedLowPtComb");
        graphYieldCombSysPi00005RebinnedLowPtComb->Write("graphYieldCombSystPi00005_RebinnedLowPtComb");
        graphYieldCombStatPi00510RebinnedHighPtComb->Write("graphYieldCombStatPi00510_RebinnedHighPtComb");
        graphYieldCombSysPi00510RebinnedHighPtComb->Write("graphYieldCombSystPi00510_RebinnedHighPtComb");
        graphYieldCombStatPi00510RebinnedLowPtComb->Write("graphYieldCombStatPi00510_RebinnedLowPtComb");
        graphYieldCombSysPi00510RebinnedLowPtComb->Write("graphYieldCombSystPi00510_RebinnedLowPtComb");
        graphYieldCombStatPi01020RebinnedHighPtComb->Write("graphYieldCombStatPi01020_RebinnedHighPtComb");
        graphYieldCombSysPi01020RebinnedHighPtComb->Write("graphYieldCombSystPi01020_RebinnedHighPtComb");
        graphYieldCombStatPi01020RebinnedLowPtComb->Write("graphYieldCombStatPi01020_RebinnedLowPtComb");
        graphYieldCombSysPi01020RebinnedLowPtComb->Write("graphYieldCombSystPi01020_RebinnedLowPtComb");
        graphYieldCombStatPi02040RebinnedHighPtComb->Write("graphYieldCombStatPi02040_RebinnedHighPtComb");
        graphYieldCombSysPi02040RebinnedHighPtComb->Write("graphYieldCombSystPi02040_RebinnedHighPtComb");
        graphYieldCombStatPi02040RebinnedLowPtComb->Write("graphYieldCombStatPi02040_RebinnedLowPtComb");
        graphYieldCombSysPi02040RebinnedLowPtComb->Write("graphYieldCombSystPi02040_RebinnedLowPtComb");
        graphYieldCombStatPi04060RebinnedHighPtComb->Write("graphYieldCombStatPi04060_RebinnedHighPtComb");
        graphYieldCombSysPi04060RebinnedHighPtComb->Write("graphYieldCombSystPi04060_RebinnedHighPtComb");
        graphYieldCombStatPi04060RebinnedLowPtComb->Write("graphYieldCombStatPi04060_RebinnedLowPtComb");
        graphYieldCombSysPi04060RebinnedLowPtComb->Write("graphYieldCombSystPi04060_RebinnedLowPtComb");
        graphYieldCombStatPi06080RebinnedHighPtComb->Write("graphYieldCombStatPi06080_RebinnedHighPtComb");
        graphYieldCombSysPi06080RebinnedHighPtComb->Write("graphYieldCombSystPi06080_RebinnedHighPtComb");
        graphYieldCombStatPi06080RebinnedLowPtComb->Write("graphYieldCombStatPi06080_RebinnedLowPtComb");
        graphYieldCombSysPi06080RebinnedLowPtComb->Write("graphYieldCombSystPi06080_RebinnedLowPtComb");
        graphChargedPionSpecHighPtStat0005HighPtComb->Write("graphChargedPionSpecHighPtStat0005_RebinnedHighPtComb");
        graphChargedPionSpecHighPtSyst0005HighPtComb->Write("graphChargedPionSpecHighPtSyst0005_RebinnedHighPtComb");
        graphChargedPionSpecLowPtStat0005LowPtComb->Write("graphChargedPionSpecLowPtStat0005_RebinnedLowPtComb");
        graphChargedPionSpecLowPtSyst0005LowPtComb->Write("graphChargedPionSpecLowPtSyst0005_RebinnedLowPtComb");
        graphChargedPionSpecHighPtStat0510HighPtComb->Write("graphChargedPionSpecHighPtStat0510_RebinnedHighPtComb");
        graphChargedPionSpecHighPtSyst0510HighPtComb->Write("graphChargedPionSpecHighPtSyst0510_RebinnedHighPtComb");
        graphChargedPionSpecLowPtStat0510LowPtComb->Write("graphChargedPionSpecLowPtStat0510_RebinnedLowPtComb");
        graphChargedPionSpecLowPtSyst0510LowPtComb->Write("graphChargedPionSpecLowPtSyst0510_RebinnedLowPtComb");
        graphChargedPionSpecHighPtStat1020HighPtComb->Write("graphChargedPionSpecHighPtStat1020_RebinnedHighPtComb");
        graphChargedPionSpecHighPtSyst1020HighPtComb->Write("graphChargedPionSpecHighPtSyst1020_RebinnedHighPtComb");
        graphChargedPionSpecLowPtStat1020LowPtComb->Write("graphChargedPionSpecLowPtStat1020_RebinnedLowPtComb");
        graphChargedPionSpecLowPtSyst1020LowPtComb->Write("graphChargedPionSpecLowPtSyst1020_RebinnedLowPtComb");
        graphChargedPionSpecHighPtStat2040HighPtComb->Write("graphChargedPionSpecHighPtStat2040_RebinnedHighPtComb");
        graphChargedPionSpecHighPtSyst2040HighPtComb->Write("graphChargedPionSpecHighPtSyst2040_RebinnedHighPtComb");
        graphChargedPionSpecLowPtStat2040LowPtComb->Write("graphChargedPionSpecLowPtStat2040_RebinnedLowPtComb");
        graphChargedPionSpecLowPtSyst2040LowPtComb->Write("graphChargedPionSpecLowPtSyst2040_RebinnedLowPtComb");
        graphChargedPionSpecHighPtStat4060HighPtComb->Write("graphChargedPionSpecHighPtStat4060_RebinnedHighPtComb");
        graphChargedPionSpecHighPtSyst4060HighPtComb->Write("graphChargedPionSpecHighPtSyst4060_RebinnedHighPtComb");
        graphChargedPionSpecLowPtStat4060LowPtComb->Write("graphChargedPionSpecLowPtStat4060_RebinnedLowPtComb");
        graphChargedPionSpecLowPtSyst4060LowPtComb->Write("graphChargedPionSpecLowPtSyst4060_RebinnedLowPtComb");
        graphChargedPionSpecHighPtStat6080HighPtComb->Write("graphChargedPionSpecHighPtStat6080_RebinnedHighPtComb");
        graphChargedPionSpecHighPtSyst6080HighPtComb->Write("graphChargedPionSpecHighPtSyst6080_RebinnedHighPtComb");
        graphChargedPionSpecLowPtStat6080LowPtComb->Write("graphChargedPionSpecLowPtStat6080_RebinnedLowPtComb");
        graphChargedPionSpecLowPtSyst6080LowPtComb->Write("graphChargedPionSpecLowPtSyst6080_RebinnedLowPtComb");
        graphYieldCombStatPi02760GeVRebinnedHighPtComb->Write("graphYieldCombStatPi0_PP2760GeV_RebinnedHighPtComb");
        graphYieldCombSysPi02760GeVRebinnedHighPtComb->Write("graphYieldCombSystPi0_PP2760GeV_RebinnedHighPtComb");
        graphYieldCombStatPi02760GeVRebinnedLowPtComb->Write("graphYieldCombStatPi0_PP2760GeV_RebinnedLowPtComb");
        graphYieldCombSysPi02760GeVRebinnedLowPtComb->Write("graphYieldCombSystPi0_PP2760GeV_RebinnedLowPtComb");
        graphYieldCombStatPi02760GeVRebinnedCMSComb->Write("graphYieldCombStatPi0_PP2760GeV_RebinnedCMSComb");
        graphYieldCombSysPi02760GeVRebinnedCMSComb->Write("graphYieldCombSystPi0_PP2760GeV_RebinnedCMSComb");
        graphChargedPionSpecHighPtStatPPHighPtComb->Write("graphChargedPionSpecHighPtStat_PP2760GeV_RebinnedHighPtComb");
        graphChargedPionSpecHighPtSystPPHighPtComb->Write("graphChargedPionSpecHighPtSyst_PP2760GeV_RebinnedHighPtComb");
        graphChargedPionSpecLowPtStatPPLowPtComb->Write("graphChargedPionSpecLowPtStat_PP2760GeV_RebinnedLowPtComb");
        graphChargedPionSpecLowPtSystPPLowPtComb->Write("graphChargedPionSpecLowPtSyst_PP2760GeV_RebinnedLowPtComb");
        graphChargedPionSpecCMSStatPPCMSComb->Write("graphChargedPionSpecCMSStat_PP2760GeV_RebinnedCMSComb");
        graphChargedPionSpecCMSSystPPCMSComb->Write("graphChargedPionSpecCMSSyst_PP2760GeV_RebinnedCMSComb");
        
        graphYieldPCMStatPi00005RebinnedHighPtPCM->Write("graphYieldPCMStatPi00005_RebinnedHighPtPCM");
        graphYieldPCMSysPi00005RebinnedHighPtPCM->Write("graphYieldPCMSystPi00005_RebinnedHighPtPCM");
        graphYieldPCMStatPi00005RebinnedLowPtPCM->Write("graphYieldPCMStatPi00005_RebinnedLowPtPCM");
        graphYieldPCMSysPi00005RebinnedLowPtPCM->Write("graphYieldPCMSystPi00005_RebinnedLowPtPCM");
        graphYieldPCMStatPi00510RebinnedHighPtPCM->Write("graphYieldPCMStatPi00510_RebinnedHighPtPCM");
        graphYieldPCMSysPi00510RebinnedHighPtPCM->Write("graphYieldPCMSystPi00510_RebinnedHighPtPCM");
        graphYieldPCMStatPi00510RebinnedLowPtPCM->Write("graphYieldPCMStatPi00510_RebinnedLowPtPCM");
        graphYieldPCMSysPi00510RebinnedLowPtPCM->Write("graphYieldPCMSystPi00510_RebinnedLowPtPCM");
        graphYieldPCMStatPi01020RebinnedHighPtPCM->Write("graphYieldPCMStatPi01020_RebinnedHighPtPCM");
        graphYieldPCMSysPi01020RebinnedHighPtPCM->Write("graphYieldPCMSystPi01020_RebinnedHighPtPCM");
        graphYieldPCMStatPi01020RebinnedLowPtPCM->Write("graphYieldPCMStatPi01020_RebinnedLowPtPCM");
        graphYieldPCMSysPi01020RebinnedLowPtPCM->Write("graphYieldPCMSystPi01020_RebinnedLowPtPCM");
        graphYieldPCMStatPi02040RebinnedHighPtPCM->Write("graphYieldPCMStatPi02040_RebinnedHighPtPCM");
        graphYieldPCMSysPi02040RebinnedHighPtPCM->Write("graphYieldPCMSystPi02040_RebinnedHighPtPCM");
        graphYieldPCMStatPi02040RebinnedLowPtPCM->Write("graphYieldPCMStatPi02040_RebinnedLowPtPCM");
        graphYieldPCMSysPi02040RebinnedLowPtPCM->Write("graphYieldPCMSystPi02040_RebinnedLowPtPCM");
        graphYieldPCMStatPi04060RebinnedHighPtPCM->Write("graphYieldPCMStatPi04060_RebinnedHighPtPCM");
        graphYieldPCMSysPi04060RebinnedHighPtPCM->Write("graphYieldPCMSystPi04060_RebinnedHighPtPCM");
        graphYieldPCMStatPi04060RebinnedLowPtPCM->Write("graphYieldPCMStatPi04060_RebinnedLowPtPCM");
        graphYieldPCMSysPi04060RebinnedLowPtPCM->Write("graphYieldPCMSystPi04060_RebinnedLowPtPCM");
        graphYieldPCMStatPi06080RebinnedHighPtPCM->Write("graphYieldPCMStatPi06080_RebinnedHighPtPCM");
        graphYieldPCMSysPi06080RebinnedHighPtPCM->Write("graphYieldPCMSystPi06080_RebinnedHighPtPCM");
        graphYieldPCMStatPi06080RebinnedLowPtPCM->Write("graphYieldPCMStatPi06080_RebinnedLowPtPCM");
        graphYieldPCMSysPi06080RebinnedLowPtPCM->Write("graphYieldPCMSystPi06080_RebinnedLowPtPCM");
        graphChargedPionSpecHighPtStat0005HighPtPCM->Write("graphChargedPionSpecHighPtStat0005_RebinnedHighPtPCM");
        graphChargedPionSpecHighPtSyst0005HighPtPCM->Write("graphChargedPionSpecHighPtSyst0005_RebinnedHighPtPCM");
        graphChargedPionSpecLowPtStat0005LowPtPCM->Write("graphChargedPionSpecLowPtStat0005_RebinnedLowPtPCM");
        graphChargedPionSpecLowPtSyst0005LowPtPCM->Write("graphChargedPionSpecLowPtSyst0005_RebinnedLowPtPCM");
        graphChargedPionSpecHighPtStat0510HighPtPCM->Write("graphChargedPionSpecHighPtStat0510_RebinnedHighPtPCM");
        graphChargedPionSpecHighPtSyst0510HighPtPCM->Write("graphChargedPionSpecHighPtSyst0510_RebinnedHighPtPCM");
        graphChargedPionSpecLowPtStat0510LowPtPCM->Write("graphChargedPionSpecLowPtStat0510_RebinnedLowPtPCM");
        graphChargedPionSpecLowPtSyst0510LowPtPCM->Write("graphChargedPionSpecLowPtSyst0510_RebinnedLowPtPCM");
        graphChargedPionSpecHighPtStat1020HighPtPCM->Write("graphChargedPionSpecHighPtStat1020_RebinnedHighPtPCM");
        graphChargedPionSpecHighPtSyst1020HighPtPCM->Write("graphChargedPionSpecHighPtSyst1020_RebinnedHighPtPCM");
        graphChargedPionSpecLowPtStat1020LowPtPCM->Write("graphChargedPionSpecLowPtStat1020_RebinnedLowPtPCM");
        graphChargedPionSpecLowPtSyst1020LowPtPCM->Write("graphChargedPionSpecLowPtSyst1020_RebinnedLowPtPCM");
        graphChargedPionSpecHighPtStat2040HighPtPCM->Write("graphChargedPionSpecHighPtStat2040_RebinnedHighPtPCM");
        graphChargedPionSpecHighPtSyst2040HighPtPCM->Write("graphChargedPionSpecHighPtSyst2040_RebinnedHighPtPCM");
        graphChargedPionSpecLowPtStat2040LowPtPCM->Write("graphChargedPionSpecLowPtStat2040_RebinnedLowPtPCM");
        graphChargedPionSpecLowPtSyst2040LowPtPCM->Write("graphChargedPionSpecLowPtSyst2040_RebinnedLowPtPCM");
        graphChargedPionSpecHighPtStat4060HighPtPCM->Write("graphChargedPionSpecHighPtStat4060_RebinnedHighPtPCM");
        graphChargedPionSpecHighPtSyst4060HighPtPCM->Write("graphChargedPionSpecHighPtSyst4060_RebinnedHighPtPCM");
        graphChargedPionSpecLowPtStat4060LowPtPCM->Write("graphChargedPionSpecLowPtStat4060_RebinnedLowPtPCM");
        graphChargedPionSpecLowPtSyst4060LowPtPCM->Write("graphChargedPionSpecLowPtSyst4060_RebinnedLowPtPCM");
        graphChargedPionSpecHighPtStat6080HighPtPCM->Write("graphChargedPionSpecHighPtStat6080_RebinnedHighPtPCM");
        graphChargedPionSpecHighPtSyst6080HighPtPCM->Write("graphChargedPionSpecHighPtSyst6080_RebinnedHighPtPCM");
        graphChargedPionSpecLowPtStat6080LowPtPCM->Write("graphChargedPionSpecLowPtStat6080_RebinnedLowPtPCM");
        graphChargedPionSpecLowPtSyst6080LowPtPCM->Write("graphChargedPionSpecLowPtSyst6080_RebinnedLowPtPCM");
        graphYieldPCMStatPi02760GeVRebinnedHighPtPCM->Write("graphYieldPCMStatPi0_PP2760GeV_RebinnedHighPtPCM");
        graphYieldPCMSysPi02760GeVRebinnedHighPtPCM->Write("graphYieldPCMSystPi0_PP2760GeV_RebinnedHighPtPCM");
        graphYieldPCMStatPi02760GeVRebinnedLowPtPCM->Write("graphYieldPCMStatPi0_PP2760GeV_RebinnedLowPtPCM");
        graphYieldPCMSysPi02760GeVRebinnedLowPtPCM->Write("graphYieldPCMSystPi0_PP2760GeV_RebinnedLowPtPCM");
        graphYieldPCMStatPi02760GeVRebinnedCMSPCM->Write("graphYieldPCMStatPi0_PP2760GeV_RebinnedCMSPCM");
        graphYieldPCMSysPi02760GeVRebinnedCMSPCM->Write("graphYieldPCMSystPi0_PP2760GeV_RebinnedCMSPCM");
        graphChargedPionSpecHighPtStatPPHighPtPCM->Write("graphChargedPionSpecHighPtStat_PP2760GeV_RebinnedHighPtPCM");
        graphChargedPionSpecHighPtSystPPHighPtPCM->Write("graphChargedPionSpecHighPtSyst_PP2760GeV_RebinnedHighPtPCM");
        graphChargedPionSpecLowPtStatPPLowPtPCM->Write("graphChargedPionSpecLowPtStat_PP2760GeV_RebinnedLowPtPCM");
        graphChargedPionSpecLowPtSystPPLowPtPCM->Write("graphChargedPionSpecLowPtSyst_PP2760GeV_RebinnedLowPtPCM");
        graphChargedPionSpecCMSStatPPCMSPCM->Write("graphChargedPionSpecCMSStat_PP2760GeV_RebinnedCMSPCM");
        graphChargedPionSpecCMSSystPPCMSPCM->Write("graphChargedPionSpecCMSSyst_PP2760GeV_RebinnedCMSPCM");

        
        graphYieldPHOSStatPi00005RebinnedHighPtPHOS->Write("graphYieldPHOSStatPi00005_RebinnedHighPtPHOS");
        graphYieldPHOSSysPi00005RebinnedHighPtPHOS->Write("graphYieldPHOSSystPi00005_RebinnedHighPtPHOS");
        graphYieldPHOSStatPi00005RebinnedLowPtPHOS->Write("graphYieldPHOSStatPi00005_RebinnedLowPtPHOS");
        graphYieldPHOSSysPi00005RebinnedLowPtPHOS->Write("graphYieldPHOSSystPi00005_RebinnedLowPtPHOS");
        graphYieldPHOSStatPi00510RebinnedHighPtPHOS->Write("graphYieldPHOSStatPi00510_RebinnedHighPtPHOS");
        graphYieldPHOSSysPi00510RebinnedHighPtPHOS->Write("graphYieldPHOSSystPi00510_RebinnedHighPtPHOS");
        graphYieldPHOSStatPi00510RebinnedLowPtPHOS->Write("graphYieldPHOSStatPi00510_RebinnedLowPtPHOS");
        graphYieldPHOSSysPi00510RebinnedLowPtPHOS->Write("graphYieldPHOSSystPi00510_RebinnedLowPtPHOS");
        graphYieldPHOSStatPi01020RebinnedHighPtPHOS->Write("graphYieldPHOSStatPi01020_RebinnedHighPtPHOS");
        graphYieldPHOSSysPi01020RebinnedHighPtPHOS->Write("graphYieldPHOSSystPi01020_RebinnedHighPtPHOS");
        graphYieldPHOSStatPi01020RebinnedLowPtPHOS->Write("graphYieldPHOSStatPi01020_RebinnedLowPtPHOS");
        graphYieldPHOSSysPi01020RebinnedLowPtPHOS->Write("graphYieldPHOSSystPi01020_RebinnedLowPtPHOS");
        graphYieldPHOSStatPi02040RebinnedHighPtPHOS->Write("graphYieldPHOSStatPi02040_RebinnedHighPtPHOS");
        graphYieldPHOSSysPi02040RebinnedHighPtPHOS->Write("graphYieldPHOSSystPi02040_RebinnedHighPtPHOS");
        graphYieldPHOSStatPi02040RebinnedLowPtPHOS->Write("graphYieldPHOSStatPi02040_RebinnedLowPtPHOS");
        graphYieldPHOSSysPi02040RebinnedLowPtPHOS->Write("graphYieldPHOSSystPi02040_RebinnedLowPtPHOS");
        graphYieldPHOSStatPi04060RebinnedHighPtPHOS->Write("graphYieldPHOSStatPi04060_RebinnedHighPtPHOS");
        graphYieldPHOSSysPi04060RebinnedHighPtPHOS->Write("graphYieldPHOSSystPi04060_RebinnedHighPtPHOS");
        graphYieldPHOSStatPi04060RebinnedLowPtPHOS->Write("graphYieldPHOSStatPi04060_RebinnedLowPtPHOS");
        graphYieldPHOSSysPi04060RebinnedLowPtPHOS->Write("graphYieldPHOSSystPi04060_RebinnedLowPtPHOS");
        graphYieldPHOSStatPi06080RebinnedHighPtPHOS->Write("graphYieldPHOSStatPi06080_RebinnedHighPtPHOS");
        graphYieldPHOSSysPi06080RebinnedHighPtPHOS->Write("graphYieldPHOSSystPi06080_RebinnedHighPtPHOS");
        graphYieldPHOSStatPi06080RebinnedLowPtPHOS->Write("graphYieldPHOSStatPi06080_RebinnedLowPtPHOS");
        graphYieldPHOSSysPi06080RebinnedLowPtPHOS->Write("graphYieldPHOSSystPi06080_RebinnedLowPtPHOS");
        graphChargedPionSpecHighPtStat0005HighPtPHOS->Write("graphChargedPionSpecHighPtStat0005_RebinnedHighPtPHOS");
        graphChargedPionSpecHighPtSyst0005HighPtPHOS->Write("graphChargedPionSpecHighPtSyst0005_RebinnedHighPtPHOS");
        graphChargedPionSpecLowPtStat0005LowPtPHOS->Write("graphChargedPionSpecLowPtStat0005_RebinnedLowPtPHOS");
        graphChargedPionSpecLowPtSyst0005LowPtPHOS->Write("graphChargedPionSpecLowPtSyst0005_RebinnedLowPtPHOS");
        graphChargedPionSpecHighPtStat0510HighPtPHOS->Write("graphChargedPionSpecHighPtStat0510_RebinnedHighPtPHOS");
        graphChargedPionSpecHighPtSyst0510HighPtPHOS->Write("graphChargedPionSpecHighPtSyst0510_RebinnedHighPtPHOS");
        graphChargedPionSpecLowPtStat0510LowPtPHOS->Write("graphChargedPionSpecLowPtStat0510_RebinnedLowPtPHOS");
        graphChargedPionSpecLowPtSyst0510LowPtPHOS->Write("graphChargedPionSpecLowPtSyst0510_RebinnedLowPtPHOS");
        graphChargedPionSpecHighPtStat1020HighPtPHOS->Write("graphChargedPionSpecHighPtStat1020_RebinnedHighPtPHOS");
        graphChargedPionSpecHighPtSyst1020HighPtPHOS->Write("graphChargedPionSpecHighPtSyst1020_RebinnedHighPtPHOS");
        graphChargedPionSpecLowPtStat1020LowPtPHOS->Write("graphChargedPionSpecLowPtStat1020_RebinnedLowPtPHOS");
        graphChargedPionSpecLowPtSyst1020LowPtPHOS->Write("graphChargedPionSpecLowPtSyst1020_RebinnedLowPtPHOS");
        graphChargedPionSpecHighPtStat2040HighPtPHOS->Write("graphChargedPionSpecHighPtStat2040_RebinnedHighPtPHOS");
        graphChargedPionSpecHighPtSyst2040HighPtPHOS->Write("graphChargedPionSpecHighPtSyst2040_RebinnedHighPtPHOS");
        graphChargedPionSpecLowPtStat2040LowPtPHOS->Write("graphChargedPionSpecLowPtStat2040_RebinnedLowPtPHOS");
        graphChargedPionSpecLowPtSyst2040LowPtPHOS->Write("graphChargedPionSpecLowPtSyst2040_RebinnedLowPtPHOS");
        graphChargedPionSpecHighPtStat4060HighPtPHOS->Write("graphChargedPionSpecHighPtStat4060_RebinnedHighPtPHOS");
        graphChargedPionSpecHighPtSyst4060HighPtPHOS->Write("graphChargedPionSpecHighPtSyst4060_RebinnedHighPtPHOS");
        graphChargedPionSpecLowPtStat4060LowPtPHOS->Write("graphChargedPionSpecLowPtStat4060_RebinnedLowPtPHOS");
        graphChargedPionSpecLowPtSyst4060LowPtPHOS->Write("graphChargedPionSpecLowPtSyst4060_RebinnedLowPtPHOS");
        graphChargedPionSpecHighPtStat6080HighPtPHOS->Write("graphChargedPionSpecHighPtStat6080_RebinnedHighPtPHOS");
        graphChargedPionSpecHighPtSyst6080HighPtPHOS->Write("graphChargedPionSpecHighPtSyst6080_RebinnedHighPtPHOS");
        graphChargedPionSpecLowPtStat6080LowPtPHOS->Write("graphChargedPionSpecLowPtStat6080_RebinnedLowPtPHOS");
        graphChargedPionSpecLowPtSyst6080LowPtPHOS->Write("graphChargedPionSpecLowPtSyst6080_RebinnedLowPtPHOS");
        graphYieldPHOSStatPi02760GeVRebinnedHighPtPHOS->Write("graphYieldPHOSStatPi0_PP2760GeV_RebinnedHighPtPHOS");
        graphYieldPHOSSysPi02760GeVRebinnedHighPtPHOS->Write("graphYieldPHOSSystPi0_PP2760GeV_RebinnedHighPtPHOS");
        graphYieldPHOSStatPi02760GeVRebinnedLowPtPHOS->Write("graphYieldPHOSStatPi0_PP2760GeV_RebinnedLowPtPHOS");
        graphYieldPHOSSysPi02760GeVRebinnedLowPtPHOS->Write("graphYieldPHOSSystPi0_PP2760GeV_RebinnedLowPtPHOS");
        graphChargedPionSpecHighPtStatPPHighPtPHOS->Write("graphChargedPionSpecHighPtStat_PP2760GeV_RebinnedHighPtPHOS");
        graphChargedPionSpecHighPtSystPPHighPtPHOS->Write("graphChargedPionSpecHighPtSyst_PP2760GeV_RebinnedHighPtPHOS");
        graphChargedPionSpecLowPtStatPPLowPtPHOS->Write("graphChargedPionSpecLowPtStat_PP2760GeV_RebinnedLowPtPHOS");
        graphChargedPionSpecLowPtSystPPLowPtPHOS->Write("graphChargedPionSpecLowPtSyst_PP2760GeV_RebinnedLowPtPHOS");
        
        
    fCombResultsOutputRebinned.Close();
}
