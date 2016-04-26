
/****************************************************************************************************************************
******      provided by Gamma Conversion Group, PWG4,                                        *****
******      Ana Marin, marin@physi.uni-heidelberg.de                                      *****
******         Kathrin Koch, kkoch@physi.uni-heidelberg.de                                      *****
******      Friederike Bock, friederike.bock@cern.ch                                      *****
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

extern TRandom*   gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*   gSystem;
extern TMinuit*   gMinuit;

void CompareChargedAndNeutralPionDataALICEpPb_FullpT(TString suffix = "pdf"){

    //*************************************** General style settings *******************************************
    gROOT->Reset();   
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();  
    SetPlotStyle();

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurementspPb",suffix.Data(),dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    TString fileNamePCM             = "ExternalInputpPb/PCM/data_PCMResults_pPb_20150804_standard_CatErrors.root";
    TString fileNamePCMDalitz       = "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20150806.root";
    TString fileNamePHOS            = "ExternalInputpPb/PHOS/data_PHOSResults_pPb_20160208.root";
    TString fileNameEMCAL           = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_160213_pPb.root";
    
    Width_t  widthLinesBoxes        = 1;

    TString collisionSystempPb      = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     

    Size_t markerSizeComparison     = 1.5;
    TString nameHistoPCM            = "CorrectedYieldPi0";
    TString nameGraphPCM            = "Pi0SystError";

    
    //************************************ Setting colors and markers *******************************************
    TString nameMeasGlobal[4]                   = { "PCM", "PHOS", "EMCal", "PCM-Dalitz"};
    TString dateAna[4]                          = {"24.06.2015", "23.10.2014", "13.2.2016", "12.9.2014"};

    Color_t colorDet[4];
    Marker_t markerStyleDet[4];
    Size_t markerSizeDet[4];
    for (Int_t i = 0; i < 4; i++){
        colorDet[i]                             = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
        markerStyleDet[i]                       = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerSizeDet[i]                        = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*0.6;
    }    
    
    
    //******************************************* reading PCM file ***********************************************
    TFile* filePCMpPb                           = new TFile(fileNamePCM.Data());
    TDirectory* directoryPCMPi0pPb                      = (TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
    TH1D* histoPCMYieldPi0pPb                           = (TH1D*)directoryPCMPi0pPb->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb        = (TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystError"); 

    
    //******************************************* reading PHOS file ***********************************************
    TFile* filePHOSpPb                          = new TFile(fileNamePHOS.Data());
//     TDirectory *directoryPHOSPi0pPb             = (TDirectory*)filePHOSpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
    TH1D* histoPHOSYieldPi0pPb                          = (TH1D*)filePHOSpPb->Get("hCor_stat");      
    TH1D* histoPHOSYieldPi0pPbSys                       = (TH1D*)filePHOSpPb->Get("hCor_syst");      

    //******************************************* reading EMCAL file ***********************************************
    TFile* fileEMCALpPb                         = new TFile(fileNameEMCAL.Data());
    TDirectory* directoryEMCALPi0pPb                    = (TDirectory*)fileEMCALpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
    TH1D* histoEMCALYieldPi0pPb                         = (TH1D*)directoryEMCALPi0pPb->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphEMCALYieldPi0SysErrpPb      = (TGraphAsymmErrors*)directoryEMCALPi0pPb->Get("Pi0SystError"); 
    
    //******************************************* reading PCM-Dalitz file *****************************************
    TFile* filePCMDalitzpPb                     = new TFile(fileNamePCMDalitz.Data());
    TDirectory* directoryPCMDalitzPi0pPb                = (TDirectory*)filePCMDalitzpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
    TH1D* histoPCMDalitzYieldPi0pPb                     = (TH1D*)directoryPCMDalitzPi0pPb->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMDalitzYieldPi0SysErrpPb  = (TGraphAsymmErrors*)directoryPCMDalitzPi0pPb->Get("Pi0SystError"); 

    
    //*************************************** reading charged pion file *******************************************
    cout << "*************************************************************************"<< endl;  
    cout << "****************** charged pions Full pT ********************************"<< endl;
    cout << "*************************************************************************"<< endl;

    TFile* fileChargedPionsFullpT       = new TFile("ExternalInputpPb/ChargedPionSpectrapPb_15_Feb_2016.root");
    TH1D *histoChargedPionSpecSys       = (TH1D*)fileChargedPionsFullpT->Get("histoChargedPionSpecFullPtSyspPb"); //0-5% 
    TH1D *histoChargedPionSpecStat      = (TH1D*)fileChargedPionsFullpT->Get("histoChargedPionSpecFullPtStatpPb");

    //*************************************** Calculating ratio to charged pions **********************************
    cout << "*************************************************************************"<< endl;  
    cout << "******************************  pPb Min bias ****************************"<< endl;
    cout << "*************************************************************************"<< endl;

    TGraphAsymmErrors* graphPCMYieldPi0SysErrpPbCopy = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrpPb->Clone("graphPCMYieldPi0SysErrpPbCopy");
    TGraphAsymmErrors* graphEMCALYieldPi0SysErrpPbCopy = (TGraphAsymmErrors*) graphEMCALYieldPi0SysErrpPb->Clone("graphEMCALYieldPi0SysErrpPbCopy");
    TGraphAsymmErrors* graphPCMDalitzYieldPi0SysErrpPbCopy = (TGraphAsymmErrors*) graphPCMDalitzYieldPi0SysErrpPb->Clone("graphPCMDalitzYieldPi0SysErrpPbCopy");

    TGraphErrors* bla1 = NULL;
    TGraphErrors* bla2 = NULL;
    TGraphErrors* bla3 = NULL;
    TGraphErrors* bla4 = NULL;
    cout << "*************************************************************************"<< endl;  
    cout << "PCM Spectrum full PT" << endl;
    cout << "*************************************************************************"<< endl;  
    TGraphErrors* graphRatioChargedPionsPCMpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMYieldPi0pPb, graphPCMYieldPi0SysErrpPbCopy, histoChargedPionSpecStat, histoChargedPionSpecSys,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
    graphRatioChargedPionsPCMpPb->Print();
    
    cout << "*************************************************************************"<< endl;  
    cout << "EMCAL Spectrum full PT" << endl;
    cout << "*************************************************************************"<< endl;  
    TGraphErrors* graphRatioChargedPionsEMCALpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoEMCALYieldPi0pPb, graphEMCALYieldPi0SysErrpPbCopy, histoChargedPionSpecStat, histoChargedPionSpecSys,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
    graphRatioChargedPionsEMCALpPb->Print();
    
    cout << "*************************************************************************"<< endl;  
    cout << "PHOS Spectrum full Pt" << endl;
    cout << "*************************************************************************"<< endl;  
    TGraphErrors* graphRatioChargedPionsPHOSpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoPHOSYieldPi0pPb, histoPHOSYieldPi0pPb, histoChargedPionSpecStat, histoChargedPionSpecSys,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
    graphRatioChargedPionsPHOSpPb->Print();

    cout << "*************************************************************************"<< endl;  
    cout << "PCM-Dalitz Spectrum full PT" << endl;
    cout << "*************************************************************************"<< endl;  
    TGraphErrors* graphRatioChargedPionsPCMDalitzpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMDalitzYieldPi0pPb, graphPCMDalitzYieldPi0SysErrpPbCopy, histoChargedPionSpecStat, histoChargedPionSpecSys,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
    graphRatioChargedPionsPCMDalitzpPb->Print();

    //************************************************************************************************************
    //******************************  plotting just minBias individual measurements ******************************
    //************************************************************************************************************

    TCanvas* canvasCompYieldpPbInd = new TCanvas("canvasCompYieldpPbInd","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldpPbInd,  0.09, 0.02, 0.02, 0.12);

    canvasCompYieldpPbInd->SetLogx();
    TH2F * histo2DCompCombinedRatio2;
    histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.3,25.,1000,0.2,4.   );
    SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.05,0.064, 0.05,0.06, 0.8,0.75, 512, 505);
    histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,18.);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.1,2.1);
    histo2DCompCombinedRatio2->GetXaxis()->SetLabelOffset(-0.015);
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPCMpPb, markerStyleDet[0] , markerSizeDet[0], colorDet[0], colorDet[0]);
        graphRatioChargedPionsPCMpPb->Draw("E1psame");

        DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPHOSpPb, markerStyleDet[1], markerSizeDet[1],  colorDet[1], colorDet[1]);
        graphRatioChargedPionsPHOSpPb->Draw("E1psame");

        DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsEMCALpPb, markerStyleDet[2], markerSizeDet[2],  colorDet[2], colorDet[2]);
        graphRatioChargedPionsEMCALpPb->Draw("E1psame");
        
        DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPCMDalitzpPb, markerStyleDet[3], markerSizeDet[3],  colorDet[3], colorDet[3]);
        graphRatioChargedPionsPCMDalitzpPb->Draw("E1psame");
        
        TLatex *labelRatioPi0pPb = new TLatex(0.13,0.91,collisionSystempPb.Data());
        SetStyleTLatex( labelRatioPi0pPb, 0.05,4);
        labelRatioPi0pPb->Draw();
        
        TLegend* legendPi0CompIndChargedPionspPb = GetAndSetLegend2(0.13, 0.15, 0.96, 0.15+(0.05*4/2), 25,2);
        legendPi0CompIndChargedPionspPb->SetMargin(0.14);
        legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsPCMpPb,Form("%s (%s)",nameMeasGlobal[0].Data(),dateAna[0].Data()),"p");
        legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsPHOSpPb,Form("%s (%s)",nameMeasGlobal[1].Data(),dateAna[1].Data()),"p");
        legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsEMCALpPb,Form("%s (%s)",nameMeasGlobal[2].Data(),dateAna[2].Data()),"p");
        legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsPCMDalitzpPb,Form("%s (%s)",nameMeasGlobal[3].Data(),dateAna[3].Data()),"p");
        legendPi0CompIndChargedPionspPb->Draw();

        legendPi0CompIndChargedPionspPb->Draw();
        DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

    DrawGammaLines(0., 20.,1., 1.,0.1,kGray,2);


    canvasCompYieldpPbInd->Update();
    canvasCompYieldpPbInd->Print(Form("%s/ComparisonChargedToNeutralInd_pPb.%s",outputDir.Data(),suffix.Data()));


}
