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
                                        TString suffix                          = "eps"
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
    TFile* fileMergedV2Clus                                  = new TFile(fileNameMergedV2Cluster.Data());
    TDirectory* directoryMergedV2ClusPi0                     = (TDirectory*)fileMergedV2Clus->Get("Pi02.76TeV"); 
        TGraphAsymmErrors* graphMergedV2ClusPi0InvXSectionStat   = (TGraphAsymmErrors*)directoryMergedV2ClusPi0->Get("graphInvCrossSectionPi0");
        TGraphAsymmErrors* graphMergedV2ClusPi0InvXSectionSys    = (TGraphAsymmErrors*)directoryMergedV2ClusPi0->Get("InvCrossSectionPi0Sys");
        for (Int_t i = 0; graphMergedV2ClusPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedV2ClusPi0InvXSectionStat->RemovePoint(0);
        }
        for (Int_t i = 0; graphMergedV2ClusPi0InvXSectionSys->GetX()[0]< 10; i++){
            graphMergedV2ClusPi0InvXSectionSys->RemovePoint(0);
        }
        
    //************************** Read data for EMCAL merged Fredi V1 clusterizer **************************************************                
    TFile* fileMergedV1Clus                                  = new TFile(fileNameMergedV1Cluster.Data());
    TDirectory* directoryMergedV1ClusPi0                     = (TDirectory*)fileMergedV1Clus->Get("Pi02.76TeV"); 
        TGraphAsymmErrors* graphMergedV1ClusPi0InvXSectionStat   = (TGraphAsymmErrors*)directoryMergedV1ClusPi0->Get("graphInvCrossSectionPi0");
        TGraphAsymmErrors* graphMergedV1ClusPi0InvXSectionSys    = (TGraphAsymmErrors*)directoryMergedV1ClusPi0->Get("InvCrossSectionPi0Sys");
        for (Int_t i = 0; graphMergedV1ClusPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedV1ClusPi0InvXSectionStat->RemovePoint(0);
        }
        for (Int_t i = 0; graphMergedV1ClusPi0InvXSectionSys->GetX()[0]< 10; i++){
            graphMergedV1ClusPi0InvXSectionSys->RemovePoint(0);
        }
        
    //************************** Read data for EMCAL merged Fredi V1 clusterizer NLM1 **************************************************                    
    TFile* fileMergedV1NLM1Clus                                  = new TFile(fileNameMergedV1NLM1Cluster.Data());
    TDirectory* directoryMergedV1NLM1ClusPi0                     = (TDirectory*)fileMergedV1NLM1Clus->Get("Pi02.76TeV"); 
        TGraphAsymmErrors* graphMergedV1NLM1ClusPi0InvXSectionStat   = (TGraphAsymmErrors*)directoryMergedV1NLM1ClusPi0->Get("graphInvCrossSectionPi0");
        TGraphAsymmErrors* graphMergedV1NLM1ClusPi0InvXSectionSys    = (TGraphAsymmErrors*)directoryMergedV1NLM1ClusPi0->Get("InvCrossSectionPi0Sys");
        for (Int_t i = 0; graphMergedV1NLM1ClusPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedV1NLM1ClusPi0InvXSectionStat->RemovePoint(0);
        }
        for (Int_t i = 0; graphMergedV1NLM1ClusPi0InvXSectionSys->GetX()[0]< 10; i++){
            graphMergedV1NLM1ClusPi0InvXSectionSys->RemovePoint(0);
        }

    //************************** Read data for EMCAL merged Fredi V1 clusterizer NLM2 **************************************************                    
    TFile* fileMergedV1NLM2Clus                                  = new TFile(fileNameMergedV1NLM2Cluster.Data());
    TDirectory* directoryMergedV1NLM2ClusPi0                     = (TDirectory*)fileMergedV1NLM2Clus->Get("Pi02.76TeV"); 
        TGraphAsymmErrors* graphMergedV1NLM2ClusPi0InvXSectionStat   = (TGraphAsymmErrors*)directoryMergedV1NLM2ClusPi0->Get("graphInvCrossSectionPi0");
        TGraphAsymmErrors* graphMergedV1NLM2ClusPi0InvXSectionSys    = (TGraphAsymmErrors*)directoryMergedV1NLM2ClusPi0->Get("InvCrossSectionPi0Sys");
        for (Int_t i = 0; graphMergedV1NLM2ClusPi0InvXSectionStat->GetX()[0]< 10; i++){
            graphMergedV1NLM2ClusPi0InvXSectionStat->RemovePoint(0);
        }
        for (Int_t i = 0; graphMergedV1NLM2ClusPi0InvXSectionSys->GetX()[0]< 10; i++){
            graphMergedV1NLM2ClusPi0InvXSectionSys->RemovePoint(0);
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
    
}
    
