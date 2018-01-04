/****************************************************************************************************************************
******         provided by Gamma Conversion Group, PWGGA,                                                     *****
******        Friederike Bock, friederike.bock@cern.ch                                                    *****
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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ExtractSignalBinning.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    // TString name;
};

void CompareMergedClusterAnaFinal   (   TString fileNameMergedHaitao            = "", 
                                        TString fileNameMergedV2Cluster         = "", 
                                        TString fileNameMergedV1Cluster         = "",  
                                        TString fileNameMergedV1NLM1Cluster     = "",  
                                        TString fileNameMergedV1NLM2Cluster     = "",  
                                        TString fileNameDiClusterV2             = "",
                                        TString suffix                          = "eps",
                                        TString fileNameMergedV2Cluster30MeV    = ""
                                    ){

    TString date = ReturnDateString();
    
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();
    SetPlotStyle();
    
    gStyle->SetEndErrorSize(0);
    
    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystem2760GeV              = "pp, #sqrt{#it{s}} = 2.76 TeV";    
    TString outputDir                           = Form("%s/%s/CompareMergedClusterAnaFinal",suffix.Data(),dateForOutput.Data());
    
    gSystem->Exec("mkdir -p "+outputDir);
        
    Width_t  widthLinesBoxes                    = 1.4;
    Width_t  widthCommonFit                     = 2;
    
    // Definition of colors, styles and markers sizes    
    Color_t  colorDet[5]                        = {kBlack, kRed+2, kBlue+2, kCyan+2, kAzure+2};   
    Style_t  markerStyleDet[5]                  = {20, 21, 24, 24, 24 };
    Size_t   markerSizeDet[5]                   = {2.,2.,2.,2.,2.};          

    Bool_t enableV2CompAdd                      = kFALSE;
    if (fileNameMergedV2Cluster30MeV.CompareTo("") != 0){
        enableV2CompAdd                         = kTRUE;
    }

    //************************** Read data for EMCAL merged Haitao **************************************************    
    TFile* fileMergedHaitao                                     = new TFile(fileNameMergedHaitao.Data());
    TDirectory* directoryMergedHaitaoPi0                        = (TDirectory*)fileMergedHaitao->Get("Pi02.76TeV"); 
        TH1D* histoMergedHaitaoPi0InvXSectionStat               = (TH1D*)directoryMergedHaitaoPi0->Get("InvCrossSection_With_stat");
        TH1D* histoMergedHaitaoPi0InvXSectionSys                = (TH1D*)directoryMergedHaitaoPi0->Get("InvCrossSection_With_syst");
        TGraphAsymmErrors* graphMergedHaitaoPi0InvXSectionStat  = new TGraphAsymmErrors(histoMergedHaitaoPi0InvXSectionStat);
        TGraphAsymmErrors* graphMergedHaitaoPi0InvXSectionSys   = new TGraphAsymmErrors(histoMergedHaitaoPi0InvXSectionSys);
        for (Int_t i = 0; graphMergedHaitaoPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedHaitaoPi0InvXSectionStat->RemovePoint(0);
            graphMergedHaitaoPi0InvXSectionSys->RemovePoint(0);
        }
        while (graphMergedHaitaoPi0InvXSectionStat->GetX()[graphMergedHaitaoPi0InvXSectionStat->GetN()-1] > 40)
            graphMergedHaitaoPi0InvXSectionStat-> RemovePoint(graphMergedHaitaoPi0InvXSectionStat->GetN()-1);
        while (graphMergedHaitaoPi0InvXSectionSys->GetX()[graphMergedHaitaoPi0InvXSectionSys->GetN()-1] > 40)
            graphMergedHaitaoPi0InvXSectionSys-> RemovePoint(graphMergedHaitaoPi0InvXSectionSys->GetN()-1);
            
    //************************** Read data for EMCAL merged Fredi V2 clusterizer **************************************************            
    TFile* fileMergedV2Clus                                     = new TFile(fileNameMergedV2Cluster.Data());
    TDirectory* directoryMergedV2ClusPi0                        = (TDirectory*)fileMergedV2Clus->Get("Pi02.76TeV"); 
        TH1F* histoTriggRejectEMC1vsINT1v2Clus                  = (TH1F*)fileMergedV2Clus->Get("TriggRejectvsPt_EMC1_INT1");
        TH1F* histoTriggRejectEMC7vsINT7v2Clus                  = (TH1F*)fileMergedV2Clus->Get("TriggRejectvsPt_EMC7_INT7");
        TH1F* histoTriggRejectEG2vsEMC7v2Clus                   = (TH1F*)fileMergedV2Clus->Get("TriggRejectvsPt_EG2_EMC7");
        TH1F* histoTriggRejectEG1vsEG2v2Clus                    = (TH1F*)fileMergedV2Clus->Get("TriggRejectvsPt_EG1_EG2");
        TGraphAsymmErrors* graphMergedV2ClusPi0InvXSectionStat      = (TGraphAsymmErrors*)directoryMergedV2ClusPi0->Get("graphInvCrossSectionPi0");
        TGraphAsymmErrors* graphMergedV2ClusPi0InvXSectionSys       = (TGraphAsymmErrors*)directoryMergedV2ClusPi0->Get("InvCrossSectionPi0Sys");
        TGraphAsymmErrors* graphMergedV2ClusPi0Efficiency           = (TGraphAsymmErrors*)directoryMergedV2ClusPi0->Get("EfficiencyPi0");
        TGraphAsymmErrors* graphMergedV2ClusPi0Purity               = (TGraphAsymmErrors*)directoryMergedV2ClusPi0->Get("PurityPi0");        
        TGraphAsymmErrors* graphMergedV2ClusPi0PurityDivEff         = CalculateGraphErrRatioToOtherTGraphErr(graphMergedV2ClusPi0Purity, graphMergedV2ClusPi0Efficiency,kTRUE);
        TGraphAsymmErrors* graphMergedV2ClusPi0EffDivPur            = CalculateGraphErrRatioToOtherTGraphErr(graphMergedV2ClusPi0Efficiency, graphMergedV2ClusPi0Purity, kTRUE);
        for (Int_t i = 0; graphMergedV2ClusPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedV2ClusPi0InvXSectionStat->RemovePoint(0);
        }
        for (Int_t i = 0; graphMergedV2ClusPi0InvXSectionSys->GetX()[0]< 10; i++){
            graphMergedV2ClusPi0InvXSectionSys->RemovePoint(0);
        }

    //************************** Read data for EMCAL merged Fredi V2 clusterizer **************************************************            
    TFile* fileDiClusV2Clus                                     = new TFile(fileNameDiClusterV2.Data());
    TDirectory* directoryDiClusV2Pi0                            = (TDirectory*)fileDiClusV2Clus->Get("Pi02.76TeV"); 
        TGraphAsymmErrors* graphDiClusV2Pi0Efficiency               = (TGraphAsymmErrors*)directoryDiClusV2Pi0->Get("EfficiencyPi0");
        
    //************************** Read data for EMCAL merged Fredi V1 clusterizer **************************************************                
    TFile* fileMergedV1Clus                                     = new TFile(fileNameMergedV1Cluster.Data());
    TDirectory* directoryMergedV1ClusPi0                        = (TDirectory*)fileMergedV1Clus->Get("Pi02.76TeV"); 
        TGraphAsymmErrors* graphMergedV1ClusPi0InvXSectionStat      = (TGraphAsymmErrors*)directoryMergedV1ClusPi0->Get("graphInvCrossSectionPi0");
        TGraphAsymmErrors* graphMergedV1ClusPi0InvXSectionSys       = (TGraphAsymmErrors*)directoryMergedV1ClusPi0->Get("InvCrossSectionPi0Sys");
        TGraphAsymmErrors* graphMergedV1ClusPi0Efficiency           = (TGraphAsymmErrors*)directoryMergedV1ClusPi0->Get("EfficiencyPi0");
        TGraphAsymmErrors* graphMergedV1ClusPi0Purity               = (TGraphAsymmErrors*)directoryMergedV1ClusPi0->Get("PurityPi0");
        TGraphAsymmErrors* graphMergedV1ClusPi0PurityDivEff         = CalculateGraphErrRatioToOtherTGraphErr(graphMergedV1ClusPi0Purity, graphMergedV1ClusPi0Efficiency,kTRUE);
        TGraphAsymmErrors* graphMergedV1ClusPi0EffDivPur            = CalculateGraphErrRatioToOtherTGraphErr(graphMergedV1ClusPi0Efficiency, graphMergedV1ClusPi0Purity, kTRUE);
        for (Int_t i = 0; graphMergedV1ClusPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedV1ClusPi0InvXSectionStat->RemovePoint(0);
        }
        for (Int_t i = 0; graphMergedV1ClusPi0InvXSectionSys->GetX()[0]< 10; i++){
            graphMergedV1ClusPi0InvXSectionSys->RemovePoint(0);
        }
        
    //************************** Read data for EMCAL merged Fredi V1 clusterizer NLM1 **************************************************                    
    TFile* fileMergedV1NLM1Clus                                 = new TFile(fileNameMergedV1NLM1Cluster.Data());
    TDirectory* directoryMergedV1NLM1ClusPi0                    = (TDirectory*)fileMergedV1NLM1Clus->Get("Pi02.76TeV"); 
    TH1F* histoTriggRejectEMC1vsINT1v1ClusNLM1                  = (TH1F*)fileMergedV1NLM1Clus->Get("TriggRejectvsPt_EMC1_NLM1_INT1_NLM1");
    TH1F* histoTriggRejectEMC7vsINT7v1ClusNLM1                  = (TH1F*)fileMergedV1NLM1Clus->Get("TriggRejectvsPt_EMC7_NLM1_INT7_NLM1");
    TH1F* histoTriggRejectEG2vsEMC7v1ClusNLM1                   = (TH1F*)fileMergedV1NLM1Clus->Get("TriggRejectvsPt_EG2_NLM1_EMC7_NLM1");
    TH1F* histoTriggRejectEG1vsEG2v1ClusNLM1                    = (TH1F*)fileMergedV1NLM1Clus->Get("TriggRejectvsPt_EG1_NLM1_EG2_NLM1");
        TGraphAsymmErrors* graphMergedV1NLM1ClusPi0InvXSectionStat      = (TGraphAsymmErrors*)directoryMergedV1NLM1ClusPi0->Get("graphInvCrossSectionPi0");
        TGraphAsymmErrors* graphMergedV1NLM1ClusPi0InvXSectionSys       = (TGraphAsymmErrors*)directoryMergedV1NLM1ClusPi0->Get("InvCrossSectionPi0Sys");
        TGraphAsymmErrors* graphMergedV1NLM1ClusPi0Efficiency           = (TGraphAsymmErrors*)directoryMergedV1NLM1ClusPi0->Get("EfficiencyPi0");
        TGraphAsymmErrors* graphMergedV1NLM1ClusPi0Purity               = (TGraphAsymmErrors*)directoryMergedV1NLM1ClusPi0->Get("PurityPi0");
        TGraphAsymmErrors* graphMergedV1NLM1ClusPi0PurityDivEff         = CalculateGraphErrRatioToOtherTGraphErr(graphMergedV1NLM1ClusPi0Purity, graphMergedV1NLM1ClusPi0Efficiency,kTRUE);
        TGraphAsymmErrors* graphMergedV1NLM1ClusPi0EffDivPur            = CalculateGraphErrRatioToOtherTGraphErr(graphMergedV1NLM1ClusPi0Efficiency, graphMergedV1NLM1ClusPi0Purity, kTRUE);
        for (Int_t i = 0; graphMergedV1NLM1ClusPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedV1NLM1ClusPi0InvXSectionStat->RemovePoint(0);
        }
        for (Int_t i = 0; graphMergedV1NLM1ClusPi0InvXSectionSys->GetX()[0]< 10; i++){
            graphMergedV1NLM1ClusPi0InvXSectionSys->RemovePoint(0);
        }

    //************************** Read data for EMCAL merged Fredi V1 clusterizer NLM2 **************************************************                    
    TFile* fileMergedV1NLM2Clus                                 = new TFile(fileNameMergedV1NLM2Cluster.Data());
    TDirectory* directoryMergedV1NLM2ClusPi0                    = (TDirectory*)fileMergedV1NLM2Clus->Get("Pi02.76TeV"); 
        TGraphAsymmErrors* graphMergedV1NLM2ClusPi0InvXSectionStat      = (TGraphAsymmErrors*)directoryMergedV1NLM2ClusPi0->Get("graphInvCrossSectionPi0");
        TGraphAsymmErrors* graphMergedV1NLM2ClusPi0InvXSectionSys       = (TGraphAsymmErrors*)directoryMergedV1NLM2ClusPi0->Get("InvCrossSectionPi0Sys");
        TGraphAsymmErrors* graphMergedV1NLM2ClusPi0Efficiency           = (TGraphAsymmErrors*)directoryMergedV1NLM2ClusPi0->Get("EfficiencyPi0");
        TGraphAsymmErrors* graphMergedV1NLM2ClusPi0Purity               = (TGraphAsymmErrors*)directoryMergedV1NLM2ClusPi0->Get("PurityPi0");
        TGraphAsymmErrors* graphMergedV1NLM2ClusPi0PurityDivEff         = CalculateGraphErrRatioToOtherTGraphErr(graphMergedV1NLM2ClusPi0Purity, graphMergedV1NLM2ClusPi0Efficiency,kTRUE);
        TGraphAsymmErrors* graphMergedV1NLM2ClusPi0EffDivPur            = CalculateGraphErrRatioToOtherTGraphErr(graphMergedV1NLM2ClusPi0Efficiency, graphMergedV1NLM2ClusPi0Purity, kTRUE);
        for (Int_t i = 0; graphMergedV1NLM2ClusPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedV1NLM2ClusPi0InvXSectionStat->RemovePoint(0);
        }
        for (Int_t i = 0; graphMergedV1NLM2ClusPi0InvXSectionSys->GetX()[0]< 10; i++){
            graphMergedV1NLM2ClusPi0InvXSectionSys->RemovePoint(0);
        }

    //************************** Read data for EMCAL merged Fredi V2 clusterizer **************************************************            
    TFile* fileMergedV2Clus30MeV                                = NULL;
    TDirectory* directoryMergedV2ClusPi030MeV                   = NULL;
    TGraphAsymmErrors* graphMergedV2Clus30MeVPi0InvXSectionStat = NULL;
    TGraphAsymmErrors* graphMergedV2Clus30MeVPi0InvXSectionSys  = NULL;

    if (enableV2CompAdd){
        fileMergedV2Clus30MeV                                   = new TFile(fileNameMergedV2Cluster30MeV.Data());
        directoryMergedV2ClusPi030MeV                           = (TDirectory*)fileMergedV2Clus30MeV->Get("Pi02.76TeV"); 
        graphMergedV2Clus30MeVPi0InvXSectionStat                = (TGraphAsymmErrors*)directoryMergedV2ClusPi030MeV->Get("graphInvCrossSectionPi0");
        graphMergedV2Clus30MeVPi0InvXSectionSys                 = (TGraphAsymmErrors*)directoryMergedV2ClusPi030MeV->Get("InvCrossSectionPi0Sys");
        for (Int_t i = 0; graphMergedV2Clus30MeVPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedV2Clus30MeVPi0InvXSectionStat->RemovePoint(0);
        }
        for (Int_t i = 0; graphMergedV2Clus30MeVPi0InvXSectionSys->GetX()[0]< 10; i++){
            graphMergedV2Clus30MeVPi0InvXSectionSys->RemovePoint(0);
        }
    }   
        
    // **********************************************************************************************************************
    // ******************************************* Compare Haitao to all others *********************************************
    // **********************************************************************************************************************
    TGraphErrors*  dummyA                               = NULL;
    TGraphErrors*  dummyB                               = NULL;
    TGraphErrors*  dummyC                               = NULL;
    TGraphErrors*  dummyD                               = NULL;
        
    TGraphErrors*  graphRatioPi0HaitaoDivV1Tot          = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedHaitaoPi0InvXSectionStat, graphMergedHaitaoPi0InvXSectionSys,  
                                                                                                            graphMergedV1ClusPi0InvXSectionStat, graphMergedV1ClusPi0InvXSectionSys, 
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0HaitaoDivV2Tot          = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedHaitaoPi0InvXSectionStat, graphMergedHaitaoPi0InvXSectionSys,  
                                                                                                            graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionSys, 
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0HaitaoDivV1NLM1Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedHaitaoPi0InvXSectionStat, graphMergedHaitaoPi0InvXSectionSys,  
                                                                                                            graphMergedV1NLM1ClusPi0InvXSectionStat, graphMergedV1NLM1ClusPi0InvXSectionSys, 
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0HaitaoDivV1NLM2Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedHaitaoPi0InvXSectionStat, graphMergedHaitaoPi0InvXSectionSys,  
                                                                                                            graphMergedV1NLM2ClusPi0InvXSectionStat, graphMergedV1NLM2ClusPi0InvXSectionSys, 
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    
    // **********************************************************************************************************************
    // ******************************************* Compare V1 clusterizer (Fredi) to all others *****************************
    // **********************************************************************************************************************
    TGraphErrors*  graphRatioPi0V1ClusDivHaitaoTot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1ClusPi0InvXSectionStat, graphMergedV1ClusPi0InvXSectionSys, 
                                                                                                            graphMergedHaitaoPi0InvXSectionStat, graphMergedHaitaoPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V1ClusDivV1NLM1Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1ClusPi0InvXSectionStat, graphMergedV1ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV1NLM1ClusPi0InvXSectionStat, graphMergedV1NLM1ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V1ClusDivV1NLM2Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1ClusPi0InvXSectionStat, graphMergedV1ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV1NLM2ClusPi0InvXSectionStat, graphMergedV1NLM2ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V1ClusDivV2Tot          = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1ClusPi0InvXSectionStat, graphMergedV1ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    // **********************************************************************************************************************
    // ******************************************* Compare V2 clusterizer (Fredi) to all others *****************************
    // **********************************************************************************************************************
    TGraphErrors*  graphRatioPi0V2ClusDivHaitaoTot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionSys, 
                                                                                                            graphMergedHaitaoPi0InvXSectionStat, graphMergedHaitaoPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V2ClusDivV1Tot          = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV1ClusPi0InvXSectionStat, graphMergedV1ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V2ClusDivV1NLM1Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV1NLM1ClusPi0InvXSectionStat, graphMergedV1NLM1ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V2ClusDivV1NLM2Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV1NLM2ClusPi0InvXSectionStat, graphMergedV1NLM2ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    // **********************************************************************************************************************
    // ******************************************* Compare V1 clusterizer NLM=1 (Fredi) to all others ***********************
    // **********************************************************************************************************************
    TGraphErrors*  graphRatioPi0V1NLM1ClusDivHaitaoTot  = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1NLM1ClusPi0InvXSectionStat, graphMergedV1NLM1ClusPi0InvXSectionSys, 
                                                                                                            graphMergedHaitaoPi0InvXSectionStat, graphMergedHaitaoPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V1NLM1ClusDivV1NLM2Tot  = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1NLM1ClusPi0InvXSectionStat, graphMergedV1NLM1ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV1NLM2ClusPi0InvXSectionStat, graphMergedV1NLM2ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V1NLM1ClusDivV1Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1NLM1ClusPi0InvXSectionStat, graphMergedV1NLM1ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV1ClusPi0InvXSectionStat, graphMergedV1ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V1NLM1ClusDivV2Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1NLM1ClusPi0InvXSectionStat, graphMergedV1NLM1ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    // **********************************************************************************************************************
    // ******************************************* Compare V1 clusterizer NLM=2 (Fredi) to all others ***********************
    // **********************************************************************************************************************
    TGraphErrors*  graphRatioPi0V1NLM2ClusDivHaitaoTot  = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1NLM2ClusPi0InvXSectionStat, graphMergedV1NLM2ClusPi0InvXSectionSys, 
                                                                                                            graphMergedHaitaoPi0InvXSectionStat, graphMergedHaitaoPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V1NLM2ClusDivV1NLM1Tot  = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1NLM2ClusPi0InvXSectionStat, graphMergedV1NLM2ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV1NLM1ClusPi0InvXSectionStat, graphMergedV1NLM1ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V1NLM2ClusDivV1Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1NLM2ClusPi0InvXSectionStat, graphMergedV1NLM2ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV1ClusPi0InvXSectionStat, graphMergedV1ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors*  graphRatioPi0V1NLM2ClusDivV2Tot      = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV1NLM2ClusPi0InvXSectionStat, graphMergedV1NLM2ClusPi0InvXSectionSys, 
                                                                                                            graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    TGraphErrors* graphRatioPi0V2Clus30MeVDivV2Tot      = NULL;
    if (enableV2CompAdd){
        graphRatioPi0V2Clus30MeVDivV2Tot                = CalculateRatioBetweenSpectraWithDifferentBinning( graphMergedV2Clus30MeVPi0InvXSectionStat, graphMergedV2Clus30MeVPi0InvXSectionSys, 
                                                                                                            graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionSys,  
                                                                                                            kTRUE,  kTRUE, 
                                                                                                            &dummyA, &dummyB, 
                                                                                                            &dummyC, &dummyD)    ;
    }    
                                                                                                            
    // **********************************************************************************************************************
    // ************************* Calculate graphs with relativ error around 1 for all 5 measurements ************************
    // **********************************************************************************************************************
    TGraphAsymmErrors* relStatErrHaitao                 = CalculateGraphErrRatioToOtherTGraphErr( graphMergedHaitaoPi0InvXSectionStat, graphMergedHaitaoPi0InvXSectionStat);
    TGraphAsymmErrors* relSysErrHaitao                  = CalculateGraphErrRatioToOtherTGraphErr( graphMergedHaitaoPi0InvXSectionSys, graphMergedHaitaoPi0InvXSectionSys);
    TGraphAsymmErrors* relStatErrV2Clus                 = CalculateGraphErrRatioToOtherTGraphErr( graphMergedV2ClusPi0InvXSectionStat, graphMergedV2ClusPi0InvXSectionStat);
    TGraphAsymmErrors* relSysErrV2Clus                  = CalculateGraphErrRatioToOtherTGraphErr( graphMergedV2ClusPi0InvXSectionSys, graphMergedV2ClusPi0InvXSectionSys);
    TGraphAsymmErrors* relStatErrV1Clus                 = CalculateGraphErrRatioToOtherTGraphErr( graphMergedV1ClusPi0InvXSectionStat, graphMergedV1ClusPi0InvXSectionStat);
    TGraphAsymmErrors* relSysErrV1Clus                  = CalculateGraphErrRatioToOtherTGraphErr( graphMergedV1ClusPi0InvXSectionSys, graphMergedV1ClusPi0InvXSectionSys);
    TGraphAsymmErrors* relStatErrV1NLM1Clus             = CalculateGraphErrRatioToOtherTGraphErr( graphMergedV1NLM1ClusPi0InvXSectionStat, graphMergedV1NLM1ClusPi0InvXSectionStat);
    TGraphAsymmErrors* relSysErrV1NLM1Clus              = CalculateGraphErrRatioToOtherTGraphErr( graphMergedV1NLM1ClusPi0InvXSectionSys, graphMergedV1NLM1ClusPi0InvXSectionSys);
    TGraphAsymmErrors* relStatErrV1NLM2Clus             = CalculateGraphErrRatioToOtherTGraphErr( graphMergedV1NLM2ClusPi0InvXSectionStat, graphMergedV1NLM2ClusPi0InvXSectionStat);
    TGraphAsymmErrors* relSysErrV1NLM2Clus              = CalculateGraphErrRatioToOtherTGraphErr( graphMergedV1NLM2ClusPi0InvXSectionSys, graphMergedV1NLM2ClusPi0InvXSectionSys);

    // **********************************************************************************************************************
    // ******************************************* Ratio of New to Haitao ***************************************************
    // **********************************************************************************************************************
    Int_t textSizeLabelsPixel           = 900*0.04;
   
    TCanvas* canvasRatioMeas   = new TCanvas("canvasRatioMeas","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioMeas, 0.1, 0.02, 0.035, 0.09);
    canvasRatioMeas->SetLogx();
   
    TH2F * histo2DRatioMeas;
    histo2DRatioMeas           = new TH2F("histo2DRatioMeas","histo2DRatioMeas",11000,7,60.,1000,0.2,1.8);
    SetStyleHistoTH2ForGraphs(histo2DRatioMeas, "#it{p}_{T} (GeV/#it{c})","#frac{Meas. X}{Haitao}",0.035,0.04, 0.035,0.04, 1.,1.,510,505);
    histo2DRatioMeas->GetXaxis()->SetMoreLogLabels();
    histo2DRatioMeas->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRatioMeas->GetYaxis()->SetRangeUser(-10,10);
    histo2DRatioMeas->Draw("copy");
// 
//  
        DrawGammaSetMarkerTGraphAsym(relSysErrHaitao, 20, 1, kGray+1 , kGray+1, widthLinesBoxes, kTRUE, kGray+1);
        relSysErrHaitao->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrHaitao, 20, 1, kGray+2 , kGray+2, widthLinesBoxes, kTRUE);
        relStatErrHaitao->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrV2Clus, 20, 1, kRed-6 , kRed-6, widthLinesBoxes, kTRUE);
        relStatErrV2Clus->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrV1Clus, 20, 1, kBlue-6 , kBlue-6, widthLinesBoxes, kTRUE);
        relStatErrV1Clus->Draw("E2same");
    
    
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V2ClusDivHaitaoTot, markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphRatioPi0V2ClusDivHaitaoTot->Draw("p,same,xz");

        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1ClusDivHaitaoTot, markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
        graphRatioPi0V1ClusDivHaitaoTot->Draw("p,same,xz");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1NLM1ClusDivHaitaoTot, markerStyleDet[3], markerSizeDet[3], colorDet[3] , colorDet[3]);
        graphRatioPi0V1NLM1ClusDivHaitaoTot->Draw("p,same,xz");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1NLM2ClusDivHaitaoTot, markerStyleDet[4], markerSizeDet[4], colorDet[4] , colorDet[4]);
        graphRatioPi0V1NLM2ClusDivHaitaoTot->Draw("p,same,xz");
        
        DrawGammaLines(7, 60., 1., 1.,0.1, kGray+2);
        DrawGammaLines(7, 60., 1.2, 1.2,0.1, kGray, 7);
        DrawGammaLines(7, 60., 0.8, 0.8,0.1, kGray, 7);
        DrawGammaLines(7, 60., 1.1, 1.1,0.1, kGray, 3);
        DrawGammaLines(7, 60., 0.9, 0.9,0.1, kGray, 3);
//         
        TLegend* legendRatioHaitao       = GetAndSetLegend2(0.27, 0.29-(0.035*4), 0.75, 0.29, 32);
        legendRatioHaitao->SetMargin(0.15);
        legendRatioHaitao->SetNColumns(2);
        legendRatioHaitao->AddEntry(relStatErrHaitao,"Haitao, stat Err","f");
        legendRatioHaitao->AddEntry(graphRatioPi0V2ClusDivHaitaoTot,"X: V2 clus","p");
        legendRatioHaitao->AddEntry(relSysErrHaitao,"Haitao, sys Err","f");
        legendRatioHaitao->AddEntry(graphRatioPi0V1ClusDivHaitaoTot,"X: V1 clus","p");
        legendRatioHaitao->AddEntry(relStatErrV2Clus,"V2 clus, stat Err","f");
        legendRatioHaitao->AddEntry(graphRatioPi0V1NLM1ClusDivHaitaoTot,"X: V1 clus LM=1","p");
        legendRatioHaitao->AddEntry(relStatErrV1Clus,"V1 clus, stat Err","f");
        legendRatioHaitao->AddEntry(graphRatioPi0V1NLM2ClusDivHaitaoTot,"X: V1 clus LM=2","p");
        legendRatioHaitao->Draw();
// 
        TLatex *labelRatioToOldEnergy   = new TLatex(0.15,0.89,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioToOldEnergy, 0.85*textSizeLabelsPixel,4);
        labelRatioToOldEnergy->SetTextFont(43);
        labelRatioToOldEnergy->Draw();
        TLatex *labelRatioToOldPi0      = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToOldPi0, 0.85*textSizeLabelsPixel,4);
        labelRatioToOldPi0->SetTextFont(43);
        labelRatioToOldPi0->Draw();
        
    canvasRatioMeas->SaveAs(Form("%s/Pi0_RatioOfXToHaitao_PP2760GeV.%s",outputDir.Data(),suffix.Data()));

    histo2DRatioMeas->Draw("copy");
        DrawGammaSetMarkerTGraphAsym(relSysErrHaitao, 20, 1, kGray+1 , kGray+1, widthLinesBoxes, kTRUE, kGray+1);
        relSysErrHaitao->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrHaitao, 20, 1, kGray+2 , kGray+2, widthLinesBoxes, kTRUE);
        relStatErrHaitao->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrV1Clus, 20, 1, kBlue-6 , kBlue-6, widthLinesBoxes, kTRUE);
        relStatErrV1Clus->Draw("E2same");
        
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1ClusDivHaitaoTot, markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
        graphRatioPi0V1ClusDivHaitaoTot->Draw("p,same,xz");
        
        DrawGammaLines(7, 60., 1., 1.,0.1, kGray+2);
        DrawGammaLines(7, 60., 1.2, 1.2,0.1, kGray, 7);
        DrawGammaLines(7, 60., 0.8, 0.8,0.1, kGray, 7);
        DrawGammaLines(7, 60., 1.1, 1.1,0.1, kGray, 3);
        DrawGammaLines(7, 60., 0.9, 0.9,0.1, kGray, 3);

        TLegend* legendRatioOnlyHaitao       = GetAndSetLegend2(0.27, 0.25-(0.035*2.2), 0.75, 0.25, 32);
        legendRatioOnlyHaitao->SetMargin(0.15);
        legendRatioOnlyHaitao->SetNColumns(2);
        legendRatioOnlyHaitao->AddEntry(relStatErrHaitao,"Haitao, stat Err","f");
        legendRatioOnlyHaitao->AddEntry(relSysErrHaitao,"Haitao, sys Err","f");
        legendRatioOnlyHaitao->AddEntry(relStatErrV1Clus,"V1 clus, stat Err","f");
        legendRatioOnlyHaitao->AddEntry(graphRatioPi0V1ClusDivHaitaoTot,"X: V1 clus","p");
        legendRatioOnlyHaitao->Draw();

        labelRatioToOldEnergy->Draw();
        labelRatioToOldPi0->Draw();
        
    canvasRatioMeas->SaveAs(Form("%s/Pi0_RatioOfV1ToHaitaoOnly_PP2760GeV.%s",outputDir.Data(),suffix.Data()));
   
    histo2DRatioMeas->GetYaxis()->SetTitle("#frac{Meas. X}{V1 clus}");
    histo2DRatioMeas->Draw("copy");
        DrawGammaSetMarkerTGraphAsym(relSysErrV1Clus, 20, 1, kBlue-8 , kBlue-8, widthLinesBoxes, kTRUE, kBlue-8);
        relSysErrV1Clus->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrV1Clus, 20, 1, kBlue-6 , kBlue-6, widthLinesBoxes, kTRUE);
        relStatErrV1Clus->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrV2Clus, 20, 1, kRed-6 , kRed-6, widthLinesBoxes, kTRUE);
        relStatErrV2Clus->Draw("E2same");
//         DrawGammaSetMarkerTGraphAsym(relStatErrHaitao, 20, 1, kGray+2 , kGray+2, widthLinesBoxes, kTRUE);
//         relStatErrHaitao->Draw("E2same");
        
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V2ClusDivV1Tot, markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphRatioPi0V2ClusDivV1Tot->Draw("p,same,xz");

//         DrawGammaSetMarkerTGraphErr(graphRatioPi0HaitaoDivV1Tot, markerStyleDet[0], markerSizeDet[0], colorDet[0] , colorDet[0]);
//         graphRatioPi0HaitaoDivV1Tot->Draw("p,same,xz");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1NLM1ClusDivV1Tot, markerStyleDet[3], markerSizeDet[3], colorDet[3] , colorDet[3]);
        graphRatioPi0V1NLM1ClusDivV1Tot->Draw("p,same,xz");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1NLM2ClusDivV1Tot, markerStyleDet[4], markerSizeDet[4], colorDet[4] , colorDet[4]);
        graphRatioPi0V1NLM2ClusDivV1Tot->Draw("p,same,xz");
        
        DrawGammaLines(7, 60., 1., 1.,0.1, kGray+2);
        DrawGammaLines(7, 60., 1.2, 1.2,0.1, kGray, 7);
        DrawGammaLines(7, 60., 0.8, 0.8,0.1, kGray, 7);
        DrawGammaLines(7, 60., 1.1, 1.1,0.1, kGray, 3);
        DrawGammaLines(7, 60., 0.9, 0.9,0.1, kGray, 3);
         
        TLegend* legendRatioV1Clus       = GetAndSetLegend2(0.27, 0.29-(0.035*4), 0.75, 0.29, 32);
        legendRatioV1Clus->SetMargin(0.15);
        legendRatioV1Clus->SetNColumns(2);
        legendRatioV1Clus->AddEntry(relStatErrV1Clus,"V1 clus, stat Err","f");
        legendRatioV1Clus->AddEntry(graphRatioPi0V2ClusDivV1Tot,"X: V2 clus","p");
        legendRatioV1Clus->AddEntry(relSysErrV1Clus,"V1 clus, sys Err","f");
//         legendRatioV1Clus->AddEntry(graphRatioPi0HaitaoDivV1Tot,"X: Haitao","p");
        
        legendRatioV1Clus->AddEntry(graphRatioPi0V1NLM1ClusDivV1Tot,"X: V1 clus LM=1","p");
        legendRatioV1Clus->AddEntry(relStatErrV2Clus,"V2 clus, stat Err","f");
//         legendRatioV1Clus->AddEntry(relStatErrHaitao,"Haitao, stat Err","f");
        legendRatioV1Clus->AddEntry(graphRatioPi0V1NLM2ClusDivV1Tot,"X: V1 clus LM=2","p");
        legendRatioV1Clus->Draw();
 
        labelRatioToOldEnergy->Draw();
        labelRatioToOldPi0->Draw();
        
    canvasRatioMeas->SaveAs(Form("%s/Pi0_RatioOfXToV1Clus_PP2760GeV.%s",outputDir.Data(),suffix.Data()));

    histo2DRatioMeas->GetYaxis()->SetTitle("#frac{Meas. X}{V1 clus LM=1}");
    histo2DRatioMeas->Draw("copy");
        DrawGammaSetMarkerTGraphAsym(relSysErrV1NLM1Clus, 20, 1, kCyan-8 , kCyan-8, widthLinesBoxes, kTRUE, kCyan-8);
        relSysErrV1NLM1Clus->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrV1NLM1Clus, 20, 1, kCyan-2 , kCyan-2, widthLinesBoxes, kTRUE);
        relStatErrV1NLM1Clus->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrV2Clus, 20, 1, kRed-6 , kRed-6, widthLinesBoxes, kTRUE);
        relStatErrV2Clus->Draw("E2same");
//         DrawGammaSetMarkerTGraphAsym(relStatErrHaitao, 20, 1, kGray+2 , kGray+2, widthLinesBoxes, kTRUE);
//         relStatErrHaitao->Draw("E2same");
        
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V2ClusDivV1NLM1Tot, markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphRatioPi0V2ClusDivV1NLM1Tot->Draw("p,same,xz");

//         DrawGammaSetMarkerTGraphErr(graphRatioPi0HaitaoDivV1Tot, markerStyleDet[0], markerSizeDet[0], colorDet[0] , colorDet[0]);
//         graphRatioPi0HaitaoDivV1Tot->Draw("p,same,xz");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1ClusDivV1NLM1Tot, markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
        graphRatioPi0V1ClusDivV1NLM1Tot->Draw("p,same,xz");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1NLM2ClusDivV1NLM1Tot, markerStyleDet[4], markerSizeDet[4], colorDet[4] , colorDet[4]);
        graphRatioPi0V1NLM2ClusDivV1NLM1Tot->Draw("p,same,xz");
        
        DrawGammaLines(7, 60., 1., 1.,0.1, kGray+2);
        DrawGammaLines(7, 60., 1.2, 1.2,0.1, kGray, 7);
        DrawGammaLines(7, 60., 0.8, 0.8,0.1, kGray, 7);
        DrawGammaLines(7, 60., 1.1, 1.1,0.1, kGray, 3);
        DrawGammaLines(7, 60., 0.9, 0.9,0.1, kGray, 3);
         
        TLegend* legendRatioV1NLM1Clus       = GetAndSetLegend2(0.27, 0.29-(0.035*4), 0.75, 0.29, 32);
        legendRatioV1NLM1Clus->SetMargin(0.15);
        legendRatioV1NLM1Clus->SetNColumns(2);
        legendRatioV1NLM1Clus->AddEntry(relStatErrV1NLM1Clus,"V1 clus LM=1, stat Err","f");
        legendRatioV1NLM1Clus->AddEntry(graphRatioPi0V2ClusDivV1NLM1Tot,"X: V2 clus","p");
        legendRatioV1NLM1Clus->AddEntry(relSysErrV1NLM1Clus,"V1 clus LM=1, sys Err","f");
//         legendRatioV1NLM1Clus->AddEntry(graphRatioPi0HaitaoDivV1Tot,"X: Haitao","p");
        legendRatioV1NLM1Clus->AddEntry(graphRatioPi0V1ClusDivV1NLM1Tot,"X: V1 ","p");
        legendRatioV1NLM1Clus->AddEntry(relStatErrV2Clus,"V2 clus, stat Err","f");
//         legendRatioV1NLM1Clus->AddEntry(relStatErrHaitao,"Haitao, stat Err","f");
        legendRatioV1NLM1Clus->AddEntry(graphRatioPi0V1NLM2ClusDivV1NLM1Tot,"X: V1 clus LM=2","p");
        legendRatioV1NLM1Clus->Draw();
 
        labelRatioToOldEnergy->Draw();
        labelRatioToOldPi0->Draw();
        
    canvasRatioMeas->SaveAs(Form("%s/Pi0_RatioOfXToV1NLM1Clus_PP2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    histo2DRatioMeas->GetYaxis()->SetTitle("#frac{Meas. X}{V2 clus}");
    histo2DRatioMeas->GetYaxis()->SetRangeUser(0.2,1.55);
    histo2DRatioMeas->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(relSysErrV2Clus, 20, 1, kRed-8, kRed-8, widthLinesBoxes, kTRUE, kRed-8);
        relSysErrV2Clus->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrV1NLM1Clus, 20, 1, kCyan-2 , kCyan-2, widthLinesBoxes, kTRUE);
        relStatErrV1NLM1Clus->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(relStatErrV2Clus, 20, 1, kRed-6 , kRed-6, widthLinesBoxes, kTRUE);
        relStatErrV2Clus->Draw("E2same");
        
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1NLM1ClusDivV2Tot, markerStyleDet[3], markerSizeDet[3], colorDet[3] , colorDet[3]);
        graphRatioPi0V1NLM1ClusDivV2Tot->Draw("p,same,xz");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1ClusDivV2Tot, markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
        graphRatioPi0V1ClusDivV2Tot->Draw("p,same,xz");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0V1NLM2ClusDivV2Tot, markerStyleDet[4], markerSizeDet[4], colorDet[4] , colorDet[4]);
        graphRatioPi0V1NLM2ClusDivV2Tot->Draw("p,same,xz");
        
        DrawGammaLines(7, 60., 1., 1.,0.1, kGray+2);
        DrawGammaLines(7, 60., 1.2, 1.2,0.1, kGray, 7);
        DrawGammaLines(7, 60., 0.8, 0.8,0.1, kGray, 7);
        DrawGammaLines(7, 60., 1.1, 1.1,0.1, kGray, 3);
        DrawGammaLines(7, 60., 0.9, 0.9,0.1, kGray, 3);
         
        TLegend* legendRatioV2Clus       = GetAndSetLegend2(0.27, 0.29-(0.035*4), 0.75, 0.29, 32);
        legendRatioV2Clus->SetMargin(0.15);
        legendRatioV2Clus->SetNColumns(2);
        legendRatioV2Clus->AddEntry(relStatErrV2Clus,"V2 clus, stat Err","f");
        legendRatioV2Clus->AddEntry(graphRatioPi0V1NLM1ClusDivV2Tot,"X: V1 clus LM=1","p");
        legendRatioV2Clus->AddEntry(relSysErrV2Clus,"V2 clus, sys Err","f");
        legendRatioV2Clus->AddEntry(graphRatioPi0V1ClusDivV2Tot,"X: V1 ","p");
        legendRatioV2Clus->AddEntry(relStatErrV1NLM1Clus,"V1 clus LM=1, stat Err","f");
        legendRatioV2Clus->AddEntry(graphRatioPi0V1NLM2ClusDivV2Tot,"X: V1 clus LM=2","p");
        legendRatioV2Clus->Draw();
 
        labelRatioToOldEnergy->Draw();
        labelRatioToOldPi0->Draw();
        
    canvasRatioMeas->SaveAs(Form("%s/Pi0_RatioOfXToV2Clus_PP2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    if (enableV2CompAdd){
        histo2DRatioMeas->GetYaxis()->SetTitle("#frac{Meas. X}{V2 clus}");
        histo2DRatioMeas->GetYaxis()->SetRangeUser(0.5,1.55);
        histo2DRatioMeas->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(relSysErrV2Clus, 20, 1, kRed-8, kRed-8, widthLinesBoxes, kTRUE, kRed-8);
            relSysErrV2Clus->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(relStatErrV2Clus, 20, 1, kRed-6 , kRed-6, widthLinesBoxes, kTRUE);
            relStatErrV2Clus->Draw("E2same");
            
//             TF1* pol0 = new TF1("pol0","[0]",10,40);//
//             graphRatioPi0V2Clus30MeVDivV2Tot->Fit(pol0,"NRMX0+","",10,40);
            
            DrawGammaSetMarkerTGraphErr(graphRatioPi0V2Clus30MeVDivV2Tot, markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
            graphRatioPi0V2Clus30MeVDivV2Tot->Draw("p,same,xz");
//             pol0->Draw("same");
            
            DrawGammaLines(7, 60., 1., 1.,0.1, kGray+2);
            DrawGammaLines(7, 60., 1.2, 1.2,0.1, kGray, 7);
            DrawGammaLines(7, 60., 0.8, 0.8,0.1, kGray, 7);
            DrawGammaLines(7, 60., 1.1, 1.1,0.1, kGray, 3);
            DrawGammaLines(7, 60., 0.9, 0.9,0.1, kGray, 3);
            
            TLegend* legendRatioV2Clus2       = GetAndSetLegend2(0.27, 0.25-(0.035*3), 0.75, 0.25, 32);
            legendRatioV2Clus2->SetMargin(0.15);
            legendRatioV2Clus2->SetNColumns(2);
            legendRatioV2Clus2->AddEntry(relStatErrV2Clus,"V2 clus, stat Err","f");
            legendRatioV2Clus2->AddEntry(relSysErrV2Clus,"V2 clus, sys Err","f");
            legendRatioV2Clus2->AddEntry(graphRatioPi0V2Clus30MeVDivV2Tot,"X: V2, cutoff 30MeV","p");
            legendRatioV2Clus2->Draw();
    
            labelRatioToOldEnergy->Draw();
            labelRatioToOldPi0->Draw();
            
        canvasRatioMeas->SaveAs(Form("%s/Pi0_RatioOfV2Clus30MeVToV2Clus_PP2760GeV.%s",outputDir.Data(),suffix.Data()));        
    }
    
    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for pi0 single measurement 2.76TeV **************************
    // **********************************************************************************************************************
    textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;
    cout << textSizeLabelsRel << endl;
    
    TCanvas* canvasAcceptanceTimesEff       = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
    canvasAcceptanceTimesEff->SetLogy(1);
    canvasAcceptanceTimesEff->SetLogx(1);
    
        TH2F * hist2DEff;
        hist2DEff                = new TH2F("hist2DEff", "hist2DEff",1000, 5.,  52, 1000, 1e-3, 1e-0 );
        SetStyleHistoTH2ForGraphs( hist2DEff, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff} ",  
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1);//(#times #epsilon_{pur})
        hist2DEff->GetYaxis()->SetLabelOffset(0.001);
        hist2DEff->GetXaxis()->SetLabelOffset(-0.01);
        hist2DEff->GetXaxis()->SetMoreLogLabels(kTRUE);
        hist2DEff->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphDiClusV2Pi0Efficiency,25, markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphDiClusV2Pi0Efficiency->Draw("p,same,e");        
        DrawGammaSetMarkerTGraphAsym(graphMergedV2ClusPi0Efficiency,markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphMergedV2ClusPi0Efficiency->Draw("p,same,e");
//         DrawGammaSetMarkerTGraphAsym(graphMergedV1ClusPi0Efficiency,markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
//         graphMergedV1ClusPi0Efficiency->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM1ClusPi0Efficiency,markerStyleDet[3], markerSizeDet[3], colorDet[3] , colorDet[3]);
        graphMergedV1NLM1ClusPi0Efficiency->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM2ClusPi0Efficiency,markerStyleDet[4], markerSizeDet[4], colorDet[4] , colorDet[4]);
        graphMergedV1NLM2ClusPi0Efficiency->Draw("p,same,e");

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi0->AddEntry(graphMergedV2ClusPi0Efficiency,"V2 clus","p");
        legendEffiAccPi0->AddEntry(graphDiClusV2Pi0Efficiency,"V2 di-clus","p");
//         legendEffiAccPi0->AddEntry(graphMergedV1ClusPi0Efficiency,"V1 clus","p");
        legendEffiAccPi0->AddEntry(graphMergedV1NLM1ClusPi0Efficiency,"V1 clus, LM=1","p");
        legendEffiAccPi0->AddEntry(graphMergedV1NLM2ClusPi0Efficiency,"V1 clus, LM=2","p");        
        legendEffiAccPi0->Draw();

//         TLatex *labelPerfEffi               = new TLatex(0.15,0.92,"ALICE performance");
//         SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
//         labelPerfEffi->Draw();
        TLatex *labelEnergyEffi             = new TLatex(0.15,0.92,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
        labelEnergyEffi->Draw();
        TLatex *labelDetSysEffiPi0          = new TLatex(0.15,0.87,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysEffiPi0, textSizeLabelsRel,4);
        labelDetSysEffiPi0->Draw();

        
    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_Efficiency.%s",outputDir.Data(),suffix.Data()));
    canvasAcceptanceTimesEff->SetLogy(0);
       TH2F * hist2DPur;
        hist2DPur                = new TH2F("hist2DPur", "hist2DPur",1000, 5.,  52, 1000, 0.3, 1.1 );
        SetStyleHistoTH2ForGraphs( hist2DPur, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{pur} ",  
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1);//(#times #epsilon_{pur})
        hist2DPur->GetYaxis()->SetLabelOffset(0.001);
        hist2DPur->GetXaxis()->SetLabelOffset(-0.01);
        hist2DPur->GetXaxis()->SetMoreLogLabels(kTRUE);
        hist2DPur->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphMergedV2ClusPi0Purity,markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphMergedV2ClusPi0Purity->Draw("p,same,e");
//         DrawGammaSetMarkerTGraphAsym(graphMergedV1ClusPi0Purity,markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
//         graphMergedV1ClusPi0Purity->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM1ClusPi0Purity,markerStyleDet[3], markerSizeDet[3], colorDet[3] , colorDet[3]);
        graphMergedV1NLM1ClusPi0Purity->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM2ClusPi0Purity,markerStyleDet[4], markerSizeDet[4], colorDet[4] , colorDet[4]);
        graphMergedV1NLM2ClusPi0Purity->Draw("p,same,e");

        TLegend* legendPurityPi0           = GetAndSetLegend2(0.14, 0.13, 0.35, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
        legendPurityPi0->AddEntry(graphMergedV2ClusPi0Purity,"V2 clus","p");
//         legendPurityPi0->AddEntry(graphMergedV1ClusPi0Purity,"V1 clus","p");
        legendPurityPi0->AddEntry(graphMergedV1NLM1ClusPi0Purity,"V1 clus, LM=1","p");
        legendPurityPi0->AddEntry(graphMergedV1NLM2ClusPi0Purity,"V1 clus, LM=2","p");        
        legendPurityPi0->Draw();

//         TLatex *labelPerfEffi               = new TLatex(0.15,0.92,"ALICE performance");
//         SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
//         labelPerfEffi->Draw();
        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->Draw();
        
    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_Purity.%s",outputDir.Data(),suffix.Data()));
   
    canvasAcceptanceTimesEff->SetLogy(1);
        TH2F * hist2DPurDivEff;
        hist2DPurDivEff                = new TH2F("hist2DPurDivEff", "hist2DPurDivEff",1000, 5.,  52, 1000, 0.7, 1000 );
        SetStyleHistoTH2ForGraphs( hist2DPurDivEff, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{pur}/#epsilon_{eff} ",  
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1);//(#times #epsilon_{pur})
        hist2DPurDivEff->GetYaxis()->SetLabelOffset(0.001);
        hist2DPurDivEff->GetXaxis()->SetLabelOffset(-0.01);
        hist2DPurDivEff->GetXaxis()->SetMoreLogLabels(kTRUE);
        hist2DPurDivEff->DrawCopy(); 

        cout << "effi" << endl;
        graphMergedV2ClusPi0Efficiency->Print();
        cout << "pur" << endl;
        graphMergedV2ClusPi0Purity->Print();
        cout << "pur/effi" << endl;
        graphMergedV2ClusPi0PurityDivEff->Print();
        DrawGammaSetMarkerTGraphAsym(graphMergedV2ClusPi0PurityDivEff,markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphMergedV2ClusPi0PurityDivEff->Draw("p,same,e");
//         DrawGammaSetMarkerTGraphAsym(graphMergedV1ClusPi0PurityDivEff,markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
//         graphMergedV1ClusPi0PurityDivEff->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM1ClusPi0PurityDivEff,markerStyleDet[3], markerSizeDet[3], colorDet[3] , colorDet[3]);
        graphMergedV1NLM1ClusPi0PurityDivEff->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM2ClusPi0PurityDivEff,markerStyleDet[4], markerSizeDet[4], colorDet[4] , colorDet[4]);
        graphMergedV1NLM2ClusPi0PurityDivEff->Draw("p,same,e");

        TLegend* legendPurDifEffPi0           = GetAndSetLegend2(0.62, 0.95-(3*textSizeLabelsRel) , 0.85, 0.95,textSizeLabelsPixel);
        legendPurDifEffPi0->AddEntry(graphMergedV2ClusPi0PurityDivEff,"V2 clus","p");
//         legendPurDifEffPi0->AddEntry(graphMergedV1ClusPi0PurityDivEff,"V1 clus","p");
        legendPurDifEffPi0->AddEntry(graphMergedV1NLM1ClusPi0PurityDivEff,"V1 clus, LM=1","p");
        legendPurDifEffPi0->AddEntry(graphMergedV1NLM2ClusPi0PurityDivEff,"V1 clus, LM=2","p");        
        legendPurDifEffPi0->Draw();

//         TLatex *labelPerfEffi               = new TLatex(0.15,0.92,"ALICE performance");
//         SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
//         labelPerfEffi->Draw();
        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->Draw();
        
    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_PurityDiffEfficiency.%s",outputDir.Data(),suffix.Data()));
    
        TH2F * hist2DEffDivPur;
        hist2DEffDivPur                = new TH2F("hist2DEffDivPur", "hist2DEffDivPur",1000, 5.,  52, 1000, 0.001, 1 );
        SetStyleHistoTH2ForGraphs( hist2DEffDivPur, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}/#epsilon_{pur}",  
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1);//(#times #epsilon_{pur})
        hist2DEffDivPur->GetYaxis()->SetLabelOffset(0.001);
        hist2DEffDivPur->GetXaxis()->SetLabelOffset(-0.01);
        hist2DEffDivPur->GetXaxis()->SetMoreLogLabels(kTRUE);
        hist2DEffDivPur->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphDiClusV2Pi0Efficiency,25, markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphDiClusV2Pi0Efficiency->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV2ClusPi0EffDivPur,markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphMergedV2ClusPi0EffDivPur->Draw("p,same,e");
//         DrawGammaSetMarkerTGraphAsym(graphMergedV1ClusPi0EffDivPur,markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
//         graphMergedV1ClusPi0EffDivPur->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM1ClusPi0EffDivPur,markerStyleDet[3], markerSizeDet[3], colorDet[3] , colorDet[3]);
        graphMergedV1NLM1ClusPi0EffDivPur->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM2ClusPi0EffDivPur,markerStyleDet[4], markerSizeDet[4], colorDet[4] , colorDet[4]);
        graphMergedV1NLM2ClusPi0EffDivPur->Draw("p,same,e");

        TLegend* legendEffDivPurPi0           = GetAndSetLegend2(0.62, 0.13, 0.85, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffDivPurPi0->AddEntry(graphMergedV2ClusPi0PurityDivEff,"V2 clus","p");
        legendEffDivPurPi0->AddEntry(graphDiClusV2Pi0Efficiency,"V2 di-clus","p");
//         legendEffDivPurPi0->AddEntry(graphMergedV1ClusPi0PurityDivEff,"V1 clus","p");
        legendEffDivPurPi0->AddEntry(graphMergedV1NLM1ClusPi0PurityDivEff,"V1 clus, LM=1","p");
        legendEffDivPurPi0->AddEntry(graphMergedV1NLM2ClusPi0PurityDivEff,"V1 clus, LM=2","p");        
        legendEffDivPurPi0->Draw();

        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->Draw();
        
    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_EfficiencyDiffPur.%s",outputDir.Data(),suffix.Data()));
    
    canvasAcceptanceTimesEff->SetLogy(0);
        hist2DEffDivPur->GetYaxis()->SetRangeUser(0,0.8);
        hist2DEffDivPur->DrawCopy(); 
        
        DrawGammaSetMarkerTGraphAsym(graphDiClusV2Pi0Efficiency,25, markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphDiClusV2Pi0Efficiency->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV2ClusPi0EffDivPur,markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        graphMergedV2ClusPi0EffDivPur->Draw("p,same,e");
//         DrawGammaSetMarkerTGraphAsym(graphMergedV1ClusPi0EffDivPur,markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
//         graphMergedV1ClusPi0EffDivPur->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM1ClusPi0EffDivPur,markerStyleDet[3], markerSizeDet[3], colorDet[3] , colorDet[3]);
        graphMergedV1NLM1ClusPi0EffDivPur->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphMergedV1NLM2ClusPi0EffDivPur,markerStyleDet[4], markerSizeDet[4], colorDet[4] , colorDet[4]);
        graphMergedV1NLM2ClusPi0EffDivPur->Draw("p,same,e");

        TLegend* legendEffDiffPurLinPi0           = GetAndSetLegend2(0.62, 0.95-(4*textSizeLabelsRel) , 0.85, 0.95,textSizeLabelsPixel);
        legendEffDiffPurLinPi0->AddEntry(graphMergedV2ClusPi0PurityDivEff,"V2 clus","p");
        legendEffDiffPurLinPi0->AddEntry(graphDiClusV2Pi0Efficiency,"V2 di-clus","p");
//         legendEffDiffPurLinPi0->AddEntry(graphMergedV1ClusPi0PurityDivEff,"V1 clus","p");
        legendEffDiffPurLinPi0->AddEntry(graphMergedV1NLM1ClusPi0PurityDivEff,"V1 clus, LM=1","p");
        legendEffDiffPurLinPi0->AddEntry(graphMergedV1NLM2ClusPi0PurityDivEff,"V1 clus, LM=2","p");        
        legendEffDiffPurLinPi0->Draw();


//         TLatex *labelPerfEffi               = new TLatex(0.15,0.92,"ALICE performance");
//         SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
//         labelPerfEffi->Draw();
        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->Draw();
        
    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_EfficiencyDiffPur_lin.%s",outputDir.Data(),suffix.Data()));
     
    //***************************************************************************************************************
    //**************************Plotting trigger rejection factors = fits linear scale all in one *******************
    //***************************************************************************************************************
    Size_t textSizeSpectra2         = 0.0415;
    Int_t textPixelPP               = textSizeSpectra2*1100;

    TCanvas* canvasTriggerRejectLinear = new TCanvas("canvasTriggerReject","",0,0,1500,1100);// gives the page size
    DrawGammaCanvasSettings( canvasTriggerRejectLinear, 0.118, 0.015, textSizeSpectra2, 0.08);
    canvasTriggerRejectLinear->SetLogy(0);
    
    Double_t minTriggRejectLin = 0;
    Double_t maxTriggRejectLin = 3100;
    
    TH2F * histo2DTriggRejectLinear;
    histo2DTriggRejectLinear = new TH2F("histo2DTriggRejectLinear","histo2DTriggRejectLinear",1000,0., 50,15000,minTriggRejectLin, maxTriggRejectLin);
    SetStyleHistoTH2ForGraphs(histo2DTriggRejectLinear, "#it{p}_{T} (GeV/#it{c})","#frac{N_{clus,trig EMC1} /N_{Evt, trig EMC1}}{N_{clus,trig INT1} /N_{Evt,trig INT1}}", 
                            0.85*textSizeSpectra2,textSizeSpectra2, 0.85*textSizeSpectra2,textSizeSpectra2, 0.85,1.18);
    histo2DTriggRejectLinear->GetXaxis()->SetRangeUser(0,25);
    histo2DTriggRejectLinear->DrawCopy(); 

        Double_t minPtEMC1  = 4.5; 
        Double_t maxPtEMC1  = 20.0; 

        TF1* pol0V1FitEMC1INT1 = new TF1("pol0V1FitEMC1INT1","[0]",minPtEMC1,maxPtEMC1); //
        histoTriggRejectEMC1vsINT1v1ClusNLM1->Fit(pol0V1FitEMC1INT1,"NRME0+","",minPtEMC1,maxPtEMC1);
        TF1* pol0V2FitEMC1INT1 = new TF1("pol0V2FitEMC1INT1","[0]",minPtEMC1,maxPtEMC1); //
        histoTriggRejectEMC1vsINT1v2Clus->Fit(pol0V2FitEMC1INT1,"NRME0+","",minPtEMC1,maxPtEMC1);
        
        DrawGammaSetMarker(histoTriggRejectEMC1vsINT1v1ClusNLM1,markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
        histoTriggRejectEMC1vsINT1v1ClusNLM1->Scale(pol0V2FitEMC1INT1->GetParameter(0)/pol0V1FitEMC1INT1->GetParameter(0));
        histoTriggRejectEMC1vsINT1v1ClusNLM1->Draw("p,same,e");
        DrawGammaSetMarker(histoTriggRejectEMC1vsINT1v2Clus,markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        histoTriggRejectEMC1vsINT1v2Clus->Draw("p,same,e");

        
        TF1* pol1V1FitEMC1INT1 = new TF1("pol1V1FitEMC1INT1","[0]+[1]*x",minPtEMC1,maxPtEMC1); //
        histoTriggRejectEMC1vsINT1v1ClusNLM1->Fit(pol1V1FitEMC1INT1,"NRME0+","",minPtEMC1,maxPtEMC1);
        pol1V1FitEMC1INT1->SetLineColor(colorDet[2]);
        pol1V1FitEMC1INT1->SetLineStyle(9);
        pol1V1FitEMC1INT1->SetRange(minPtEMC1,maxPtEMC1);
        pol1V1FitEMC1INT1->Draw("same");
        TF1* pol1V2FitEMC1INT1 = new TF1("pol1V2FitEMC1INT1","[0]+[1]*x",minPtEMC1,maxPtEMC1); //
        histoTriggRejectEMC1vsINT1v2Clus->Fit(pol1V2FitEMC1INT1,"NRME0+","",minPtEMC1,maxPtEMC1);
        pol1V2FitEMC1INT1->SetLineColor(colorDet[1]);
        pol1V2FitEMC1INT1->SetLineStyle(9);
        pol1V2FitEMC1INT1->SetRange(minPtEMC1,maxPtEMC1);
        pol1V2FitEMC1INT1->Draw("same");
        
        TLegend* legendTriggerRejectEMC1           = GetAndSetLegend2(0.72, 0.13, 0.85, 0.13+(2*textSizeLabelsRel),textSizeLabelsPixel);
        legendTriggerRejectEMC1->AddEntry(histoTriggRejectEMC1vsINT1v2Clus,"V2 clus","p");
        legendTriggerRejectEMC1->AddEntry(histoTriggRejectEMC1vsINT1v1ClusNLM1,"V1 clus","p");
        legendTriggerRejectEMC1->Draw();

    TLatex *labelEnergyTriggerReject             = new TLatex(0.15,0.90,collisionSystem2760GeV.Data());
    SetStyleTLatex( labelEnergyTriggerReject, textSizeLabelsRel,4);
    labelEnergyTriggerReject->Draw();
    TLatex *labelTriggerRejectEMC1             = new TLatex(0.15,0.85,"EMC1/INT1");
    SetStyleTLatex( labelTriggerRejectEMC1, textSizeLabelsRel,4);
    labelTriggerRejectEMC1->Draw();
        
    histo2DTriggRejectLinear->Draw("same,axis");
    
    canvasTriggerRejectLinear->Update();
    canvasTriggerRejectLinear->SaveAs(Form("%s/TriggerRejectionFactorsLinY_EMC1vsINT1.%s",outputDir.Data(),suffix.Data()));
    
    histo2DTriggRejectLinear->GetYaxis()->SetTitle("#frac{N_{clus,trig EMC7} /N_{Evt, trig EMC7}}{N_{clus,trig INT7} /N_{Evt,trig INT7}}");
    histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,320);
    histo2DTriggRejectLinear->GetXaxis()->SetRangeUser(0,25);
    histo2DTriggRejectLinear->DrawCopy(); 

        Double_t minPtEMC7  = 2.8; 
        Double_t maxPtEMC7  = 15.0; 

        TF1* pol0V1FitEMC7INT7 = new TF1("pol0V1FitEMC7INT7","[0]",minPtEMC7,maxPtEMC7); //
        histoTriggRejectEMC7vsINT7v1ClusNLM1->Fit(pol0V1FitEMC7INT7,"NRME0+","",minPtEMC7,maxPtEMC7);
        TF1* pol0V2FitEMC7INT7 = new TF1("pol0V2FitEMC7INT7","[0]",minPtEMC7,maxPtEMC7); //
        histoTriggRejectEMC7vsINT7v2Clus->Fit(pol0V2FitEMC7INT7,"NRME0+","",minPtEMC7,maxPtEMC7);

        DrawGammaSetMarker(histoTriggRejectEMC7vsINT7v1ClusNLM1,markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
        histoTriggRejectEMC7vsINT7v1ClusNLM1->Scale(pol0V2FitEMC7INT7->GetParameter(0)/pol0V1FitEMC7INT7->GetParameter(0));
        histoTriggRejectEMC7vsINT7v1ClusNLM1->Draw("p,same,e");
        DrawGammaSetMarker(histoTriggRejectEMC7vsINT7v2Clus,markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        histoTriggRejectEMC7vsINT7v2Clus->Draw("p,same,e");

        
        TF1* pol1V1FitEMC7INT7 = new TF1("pol1V1FitEMC7INT7","[0]+[1]*x",minPtEMC7,maxPtEMC7); //
        histoTriggRejectEMC7vsINT7v1ClusNLM1->Fit(pol1V1FitEMC7INT7,"NRME0+","",minPtEMC7,maxPtEMC7);
        pol1V1FitEMC7INT7->SetLineColor(colorDet[2]);
        pol1V1FitEMC7INT7->SetLineStyle(9);
        pol1V1FitEMC7INT7->SetRange(minPtEMC7,maxPtEMC7);
        pol1V1FitEMC7INT7->Draw("same");
        TF1* pol1V2FitEMC7INT7 = new TF1("pol1V2FitEMC7INT7","[0]+[1]*x",minPtEMC7,maxPtEMC7); //
        histoTriggRejectEMC7vsINT7v2Clus->Fit(pol1V2FitEMC7INT7,"NRME0+","",minPtEMC7,maxPtEMC7);
        pol1V2FitEMC7INT7->SetLineColor(colorDet[1]);
        pol1V2FitEMC7INT7->SetLineStyle(9);
        pol1V2FitEMC7INT7->SetRange(minPtEMC7,maxPtEMC7);
        pol1V2FitEMC7INT7->Draw("same");

        legendTriggerRejectEMC1->Draw();

    labelEnergyTriggerReject->Draw();
    TLatex *labelTriggerRejectEMC7             = new TLatex(0.15,0.85,"EMC7/INT7");
    SetStyleTLatex( labelTriggerRejectEMC7, textSizeLabelsRel,4);
    labelTriggerRejectEMC7->Draw();

        
    histo2DTriggRejectLinear->Draw("same,axis");
    
    canvasTriggerRejectLinear->Update();
    canvasTriggerRejectLinear->SaveAs(Form("%s/TriggerRejectionFactorsLinY_EMC7vsINT7.%s",outputDir.Data(),suffix.Data()));
    
    histo2DTriggRejectLinear->GetYaxis()->SetTitle("#frac{N_{clus,trig EG2} /N_{Evt, trig EG2}}{N_{clus,trig EMC7} /N_{Evt,trig EMC7}}");
    histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,25);
    histo2DTriggRejectLinear->GetXaxis()->SetRangeUser(0,30);
    histo2DTriggRejectLinear->DrawCopy(); 

        Double_t minPtEG2   = 4.5; 
        Double_t maxPtEG2   = 18.0; 
            
        TF1* pol0V1FitEG2EMC7 = new TF1("pol0V1FitEG2EMC7","[0]",minPtEG2,maxPtEG2); //
        histoTriggRejectEG2vsEMC7v1ClusNLM1->Fit(pol0V1FitEG2EMC7,"NRME0+","",minPtEG2,maxPtEG2);
        TF1* pol0V2FitEG2EMC7 = new TF1("pol0V2FitEG2EMC7","[0]",minPtEG2,maxPtEG2); //
        histoTriggRejectEG2vsEMC7v2Clus->Fit(pol0V2FitEG2EMC7,"NRME0+","",minPtEG2,maxPtEG2);
    
        DrawGammaSetMarker(histoTriggRejectEG2vsEMC7v1ClusNLM1,markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
        histoTriggRejectEG2vsEMC7v1ClusNLM1->Scale(pol0V2FitEG2EMC7->GetParameter(0)/pol0V1FitEG2EMC7->GetParameter(0));
        histoTriggRejectEG2vsEMC7v1ClusNLM1->Draw("p,same,e");
        DrawGammaSetMarker(histoTriggRejectEG2vsEMC7v2Clus,markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        histoTriggRejectEG2vsEMC7v2Clus->Draw("p,same,e");

        TF1* pol1V1FitEG2EMC7 = new TF1("pol1V1FitEG2EMC7","[0]+[1]*x",minPtEG2,maxPtEG2); //
        histoTriggRejectEG2vsEMC7v1ClusNLM1->Fit(pol1V1FitEG2EMC7,"NRME0+","",minPtEG2,maxPtEG2);
        pol1V1FitEG2EMC7->SetLineColor(colorDet[2]);
        pol1V1FitEG2EMC7->SetLineStyle(9);
        pol1V1FitEG2EMC7->SetRange(minPtEG2,maxPtEG2);
        pol1V1FitEG2EMC7->Draw("same");
        TF1* pol1V2FitEG2EMC7 = new TF1("pol1V2FitEG2EMC7","[0]+[1]*x",minPtEG2,maxPtEG2); //
        histoTriggRejectEG2vsEMC7v2Clus->Fit(pol1V2FitEG2EMC7,"NRME0+","",minPtEG2,maxPtEG2);
        pol1V2FitEG2EMC7->SetLineColor(colorDet[1]);
        pol1V2FitEG2EMC7->SetLineStyle(9);
        pol1V2FitEG2EMC7->SetRange(minPtEG2,maxPtEG2);
        pol1V2FitEG2EMC7->Draw("same");
        
        legendTriggerRejectEMC1->Draw();

    labelEnergyTriggerReject->Draw();
    TLatex *labelTriggerRejectEG2             = new TLatex(0.15,0.85,"EG2/EMC7");
    SetStyleTLatex( labelTriggerRejectEG2, textSizeLabelsRel,4);
    labelTriggerRejectEG2->Draw();

        
    histo2DTriggRejectLinear->Draw("same,axis");
    
    canvasTriggerRejectLinear->Update();
    canvasTriggerRejectLinear->SaveAs(Form("%s/TriggerRejectionFactorsLinY_EG2vsEMC7.%s",outputDir.Data(),suffix.Data()));
    
    histo2DTriggRejectLinear->GetYaxis()->SetTitle("#frac{N_{clus,trig EG1} /N_{Evt, trig EG1}}{N_{clus,trig EG2} /N_{Evt,trig EG2}}");
    histo2DTriggRejectLinear->GetYaxis()->SetRangeUser(0,10);
    histo2DTriggRejectLinear->GetXaxis()->SetRangeUser(0,50);
    histo2DTriggRejectLinear->DrawCopy(); 

        Double_t minPtEG1   = 6.5; 
        Double_t maxPtEG1   = 40.0; 
    
        TF1* pol0V1FitEG1EG2 = new TF1("pol0V1FitEG1EG2","[0]",minPtEG1,maxPtEG1); //
        histoTriggRejectEG1vsEG2v1ClusNLM1->Fit(pol0V1FitEG1EG2,"NRME0+","",minPtEG1,maxPtEG1);
        TF1* pol0V2FitEG1EG2 = new TF1("pol0V2FitEG1EG2","[0]",minPtEG1,maxPtEG1); //
        histoTriggRejectEG1vsEG2v2Clus->Fit(pol0V2FitEG1EG2,"NRME0+","",minPtEG1,maxPtEG1);

        DrawGammaSetMarker(histoTriggRejectEG1vsEG2v1ClusNLM1,markerStyleDet[2], markerSizeDet[2], colorDet[2] , colorDet[2]);
        histoTriggRejectEG1vsEG2v1ClusNLM1->Scale(pol0V2FitEG1EG2->GetParameter(0)/pol0V1FitEG1EG2->GetParameter(0));
        histoTriggRejectEG1vsEG2v1ClusNLM1->Draw("p,same,e");
        DrawGammaSetMarker(histoTriggRejectEG1vsEG2v2Clus,markerStyleDet[1], markerSizeDet[1], colorDet[1] , colorDet[1]);
        histoTriggRejectEG1vsEG2v2Clus->Draw("p,same,e");
        
        TF1* pol1V1FitEG1EG2 = new TF1("pol1V1FitEG1EG2","[0]+[1]*x",minPtEG1,maxPtEG1); //
        histoTriggRejectEG1vsEG2v1ClusNLM1->Fit(pol1V1FitEG1EG2,"NRME0+","",minPtEG1,maxPtEG1);
        pol1V1FitEG1EG2->SetLineColor(colorDet[2]);
        pol1V1FitEG1EG2->SetLineStyle(9);
        pol1V1FitEG1EG2->SetRange(minPtEG1,maxPtEG1);
        pol1V1FitEG1EG2->Draw("same");
        TF1* pol1V2FitEG1EG2 = new TF1("pol1V2FitEG1EG2","[0]+[1]*x",minPtEG1,maxPtEG1); //
        histoTriggRejectEG1vsEG2v2Clus->Fit(pol1V2FitEG1EG2,"NRME0+","",minPtEG1,maxPtEG1);
        pol1V2FitEG1EG2->SetLineColor(colorDet[1]);
        pol1V2FitEG1EG2->SetLineStyle(9);
        pol1V2FitEG1EG2->SetRange(minPtEG1,maxPtEG1);
        pol1V2FitEG1EG2->Draw("same");
        
        legendTriggerRejectEMC1->Draw();

    labelEnergyTriggerReject->Draw();
    TLatex *labelTriggerRejectEG1             = new TLatex(0.15,0.85,"EG1/EG2");
    SetStyleTLatex( labelTriggerRejectEG1, textSizeLabelsRel,4);
    labelTriggerRejectEG1->Draw();

        
    histo2DTriggRejectLinear->Draw("same,axis");
    
    canvasTriggerRejectLinear->Update();
    canvasTriggerRejectLinear->SaveAs(Form("%s/TriggerRejectionFactorsLinY_EG1vsEG2.%s",outputDir.Data(),suffix.Data()));
    
}
    
