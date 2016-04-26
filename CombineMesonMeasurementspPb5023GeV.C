/****************************************************************************************************************************
******         provided by Gamma Conversion Group, PWG4,                                                     *****
******        Ana Marin, marin@physi.uni-heidelberg.de                                                    *****
******           Kathrin Koch, kkoch@physi.uni-heidelberg.de                                                     *****
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
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;


void DrawPad(Int_t iPad,Double_t a,Double_t b, Double_t c, Double_t d){
  
  Double_t bottom=1. ;
  TPad * sp;
  if(iPad==0)  sp = new TPad(Form("spectrum%d",iPad),"",0.,0,1.,bottom*0.22*(iPad+1)+0.1);
  else  sp = new TPad(Form("spectrum%d",iPad),"",0.,bottom*0.22*iPad+0.1,1.,bottom*0.22*(iPad+1)+0.1);
 //  TPad * sp = new TPad(Form("spectrum%d",iPad),"",0.,bottom*0.25*iPad,1.,bottom*0.25*(iPad+1));
  sp->Draw();
  sp->cd();
  sp->Range(0,0,1,1);
  sp->SetFillColor(0);
  sp->SetFillStyle(0);
  sp->SetBorderSize(1);
  sp->SetTopMargin(a);
  sp->SetBottomMargin(b);
  sp->SetLeftMargin(c);
  sp->SetRightMargin(d);
  sp->SetLogx();
  sp->SetGridy();
  sp->cd() ;
}  

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //    TString name;
};

const Int_t  Ntotal = 32;
const Int_t  nPtLimits = Ntotal+1;

void CombineMesonMeasurementspPb5023GeV(     /*TString fileNamePCM = "", 
                                        TString fileNamePCMEMCAL = "", 
                                        TString fileNameEMCALLow = "",  
                                        TString fileNameEMCALFull = "",  
                                        TString suffix = "eps", 
                                        TString isMC= "", 
                                        TString thesisPlots = "", 
                                        TString bWCorrection="X"*/){    

    TString date                                        = ReturnDateString();
    
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();
    SetPlotStyle();
    
    TString suffix                                      = "pdf";
    
    TString dateForOutput                               = ReturnDateStringForOutput();
    TString outputDir                                   = Form("CombinepPbSpectra/%s/%s",suffix.Data(),dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    
    
    //cout << dateForOutput.Data() << endl;
    //___________________________________ Declaration of files _____________________________________________
    TString bWCorrection                                = "X";
    
    //TString fileNameNeutralPionDalitz                 = "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_2015-06-28.root";
    TString fileNameNeutralPionDalitz                   = "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20150806.root";
    //TString fileNameNeutralPionPCM                    = "ExternalInputpPb/PCM/data_PCMResults_pPb_20150624_standard_dc4.root";
    TString fileNameNeutralPionPCM                      = "ExternalInputpPb/PCM/data_PCMResults_pPb_20151111_standard_CatErrors.root";
							   
    //TString fileNameNeutralPionPHOS                     = "ExternalInputpPb/PHOS/data_PHOSResults_pPb_20150623.root";
    
    TString fileNameNeutralPionPHOS                     = "ExternalInputpPb/PHOS/data_PHOSResults_pPb_20160208.root";
    
    
    TString nameHistoPHOS                               = "hCor_stat";
    TString nameHistoPHOSSysErrors                      = "hCor_syst";
    TString fileNameNeutralPionEMCAL                    = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_160218_pPb.root"; 
    TString nameHistoEMCal                              = "CorrectedYieldPi0";
    TString nameHistoEMCalSysErrors                     = "Pi0SystError";
    
    
    TString collisionSystempPb                          = "p-Pb #sqrt{s_{NN}} = 5.02 TeV"; 
    
    
    Double_t pTLimits[nPtLimits]                        = { 0.3, 0.4, 0.5, 0.6, 0.7, 
                                                            0.8, 1.0, 1.2, 1.4, 1.6, 
                                                            1.8, 2.0, 2.2, 2.4, 2.6,
                                                            2.8, 3.0, 3.2, 3.4, 3.6,
                                                            3.8, 4.0, 4.5, 5.0, 5.5,
                                                            6.0, 7.0, 8.0, 10.0,12.0,
							    14.0,16.0, 20.0};
                                
    
    Int_t offSets[11]                                   =  { -1, 5, -1, 0, 0,  2, 0, 0, 0,  0, 0};
    Int_t offSetsSys[11]                                =  {  0, 6,  6, 0, 0, 3, 0, 0, 0, 0, 0};
    
    
    TH1D* statErrorCollection[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorCollection[i]                          = NULL;
    }    
    
    TGraphAsymmErrors* sysErrorCollection[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollection[i]                           = NULL;
    }    
    
    
    // **************************************************************************************
    // ****************************** Reading Dalitz ****************************************
    // **************************************************************************************
    TFile* fileNeutralPionsDalitz                       = new TFile(fileNameNeutralPionDalitz.Data());
    TDirectory* directoryDalitzPi0pPb                   = (TDirectory*) fileNeutralPionsDalitz->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
    if( ! directoryDalitzPi0pPb ){
      cout<<"Dalitz: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
      return;
    }

    TH1D* histoDalitzYieldPi0pPb                        = (TH1D*)directoryDalitzPi0pPb->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphDalitzYieldPi0pPbSystErr    = (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystError");
    
    TGraphAsymmErrors* graphRatioCombDalitzSys          = (TGraphAsymmErrors*)graphDalitzYieldPi0pPbSystErr->Clone();
    TGraphAsymmErrors* temp02                           = new TGraphAsymmErrors(histoDalitzYieldPi0pPb);

    cout<<"Dalitz systematic"<<endl;
    graphRatioCombDalitzSys->Print();
    cout<<"Dalitz statistic"<<endl;
    temp02->Print();
    
    // **************************************************************************************
    // ****************************** Reading PCM *******************************************
    // **************************************************************************************    
    TFile* fileNeutralPionPCM                           = new TFile(fileNameNeutralPionPCM.Data());
    TDirectory* directoryPCMPi0pPb                      = (TDirectory*) fileNeutralPionPCM->GetDirectory("Pi0_pPb_5.023TeV_0-100%");    
    if( ! directoryPCMPi0pPb ){
      cout<<"PCM: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
      return; 
    }
    
    TH1D* histoPCMYieldPi0pPb                           = (TH1D*)directoryPCMPi0pPb->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMYieldPi0pPbSystErr       = (TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystError");
    TGraphAsymmErrors* temp01                           = new TGraphAsymmErrors(histoPCMYieldPi0pPb);
    
    cout<<"PCM systematic"<<endl;
    graphPCMYieldPi0pPbSystErr->Print();
    cout<<"PCM statistic"<<endl;
    temp01->Print();
    
    TGraphAsymmErrors* graphRatioCombPCMSys             = (TGraphAsymmErrors*)graphPCMYieldPi0pPbSystErr->Clone();;
    
    
    // **************************************************************************************
    // ******************************* Reading PHOS *****************************************
    // **************************************************************************************    
    TFile* fileNeutralPionPHOS                          = new TFile(fileNameNeutralPionPHOS);
    //TDirectory*    directoryPHOSPi0pPb =                 (TDirectory*)fileNeutralPionPHOS->Get("Pi0_pPb_5.023TeV_0-100%"); 
    if( ! fileNeutralPionPHOS ){
      cout<<"PHOS: The file fileNeutralPionPHOS does not exist "<<endl; 
      return;  
    }
   
    TH1D* histoPHOSYieldPi0pPbStat                      = (TH1D*)fileNeutralPionPHOS->Get(nameHistoPHOS.Data());
    TH1D* histoPHOSYieldPi0pPbSyst                      = (TH1D*)fileNeutralPionPHOS->Get(nameHistoPHOSSysErrors.Data());
    TGraphAsymmErrors* graphPHOSYieldPi0pPbSystErr      = new TGraphAsymmErrors(histoPHOSYieldPi0pPbSyst);
    graphPHOSYieldPi0pPbSystErr->RemovePoint(0);
    TGraphAsymmErrors* temp03                           = new TGraphAsymmErrors(histoPHOSYieldPi0pPbStat);
    
    cout<<"Systematic PHOS "<<endl;
    graphPHOSYieldPi0pPbSystErr->Print();
    cout<<"Statistic PHOS  "<<endl;
    temp03->Print();
    TH1D * histoRatioCombPHOS                           = (TH1D*)histoPHOSYieldPi0pPbSyst->Clone();
    
    // **************************************************************************************
    // ******************************** Reading EMCal ***************************************
    // **************************************************************************************
    TFile* fileNeutralPionEMCAL                         = new TFile(fileNameNeutralPionEMCAL);
    TDirectory* directoryEMCalPi0pPb                    = (TDirectory*)fileNeutralPionEMCAL->Get("Pi0_pPb_5.023TeV_0-100%");
    if( ! directoryEMCalPi0pPb ){
      cout<<"EMCAL: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
      return;  
    }
    
    TH1D* histoEMCALYieldPi0pPbStat                     = (TH1D*)directoryEMCalPi0pPb->Get(nameHistoEMCal.Data());
    TGraphAsymmErrors* graphEMCALYieldPi0pPbSystErr     = (TGraphAsymmErrors*)directoryEMCalPi0pPb->Get(nameHistoEMCalSysErrors.Data());
    
    TGraphAsymmErrors* graphRatioCombEMCALSys           = (TGraphAsymmErrors*)graphEMCALYieldPi0pPbSystErr->Clone();;
    
    // **************************************************************************************
    // ********************************* Combine spectra ************************************
    // **************************************************************************************
    
    TString fileNameOutputWeightingOld                  = Form("%s/WeightingOld.dat",outputDir.Data());

    statErrorCollection[0]          = (TH1D*)histoPCMYieldPi0pPb->Clone("statErrPCMPi0");
    statErrorCollection[1]          = (TH1D*)histoPHOSYieldPi0pPbStat->Clone("statErrPHOSPi0");
    statErrorCollection[2]          = (TH1D*)histoEMCALYieldPi0pPbStat->Clone("statErrEMCALPi0");
    statErrorCollection[5]          = (TH1D*)histoDalitzYieldPi0pPb->Clone("statErrDalitzPi0");
    
    sysErrorCollection[0]           = (TGraphAsymmErrors*)graphPCMYieldPi0pPbSystErr->Clone("sysErrPCMPi0");
    sysErrorCollection[1]           = (TGraphAsymmErrors*)graphPHOSYieldPi0pPbSystErr->Clone("sysErrPHOSPi0");
    sysErrorCollection[2]           = (TGraphAsymmErrors*)graphEMCALYieldPi0pPbSystErr->Clone("sysErrEMCALPi0");
    sysErrorCollection[5]           = (TGraphAsymmErrors*)graphDalitzYieldPi0pPbSystErr->Clone("sysErrDalitzPi0");
    
    
    
    
    TGraphAsymmErrors* graphCombPi0InvCrossSectionStatpPb5023GeV= NULL;
    TGraphAsymmErrors* graphCombPi0InvCrossSectionSyspPb5023GeV = NULL;
    
    TGraphAsymmErrors* graphCombPi0InvCrossSectionTotpPb5023GeV = CombinePtPointsSpectraFullCorrMat(    statErrorCollection,    sysErrorCollection,     
                                                                                                        pTLimits, Ntotal,
                                                                                                        offSets, offSetsSys,
                                                                                                        graphCombPi0InvCrossSectionStatpPb5023GeV, graphCombPi0InvCrossSectionSyspPb5023GeV,
                                                                                                        fileNameOutputWeightingOld,1
                                                                                                    );
    graphCombPi0InvCrossSectionStatpPb5023GeV->Print();
    
    TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVStaClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStatpPb5023GeV->Clone();
    TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVSysClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSyspPb5023GeV->Clone();
    TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVTotClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionTotpPb5023GeV->Clone();
    
    TGraphAsymmErrors* graphRatioCombCombFit                    = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVTotClone ->Clone(); 
    TGraphAsymmErrors* graphRatioCombCombFitSta                 = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVStaClone ->Clone();
    TGraphAsymmErrors* graphRatioCombCombFitSys                 = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVSysClone ->Clone(); 

    
    
    // **************************************************************************************
    // ************************* Plotting yield of different systems ************************
    // **************************************************************************************
    TCanvas* canvasCompYieldpPbInd = new TCanvas("canvasCompYieldpPbInd","",200,10,500,400);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldpPbInd, 0.15, 0.02, 0.02, 0.12);
    
    canvasCompYieldpPbInd->SetLogx();
    canvasCompYieldpPbInd->SetLogy();
    TH2F * histo2DCompCombined;
    histo2DCompCombined = new TH2F("histo2DCompCombined","histo2DCompCombined",1000,0.3,30.,1000,1.2e-9,30.   );
     SetStyleHistoTH2ForGraphs(histo2DCompCombined, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.03,0.04, 0.03,0.04, 1.0,1.4, 512, 508);
     // SetStyleHistoTH2ForGraphs(histo2DCompCombined, "#it{p}_{T} (GeV/#it{c})","arbitrary units", 0.03,0.04, 0.03,0.04, 1.0,1.1, 512, 508);
  //   histo2DCompCombined->GetXaxis()->SetRangeUser(0.,30.);
     //  histo2DCompCombined->GetYaxis()->SetLabelColor(0);
    //  histo2DCompCombined->GetYaxis()->SetRangeUser(0.1,2.1);
    histo2DCompCombined->DrawCopy();
    
    //     TF1 * CurrentFit = new TF1();
    TF1* fitBylinkinCombPi0pPb5023GeVPt = FitObject("tcm","fitInvCrossSectionPi0","Pi0");    
    Double_t *ParameterspPb  = new Double_t[5];// = { 7.4e+10, 0.3, 1e+09,0.3,8};
        
    ParameterspPb[0] = 0.01;//7.4e+1;
    ParameterspPb[1] = 0.3;
    ParameterspPb[2] = .06;
    ParameterspPb[3] = 0.07;
    ParameterspPb[4] = 400.8;
          
    Double_t minPt = 0.3;
    Double_t maxPt = 23.0;

    fitBylinkinCombPi0pPb5023GeVPt->SetRange(minPt,maxPt);
    fitBylinkinCombPi0pPb5023GeVPt->SetParameters(ParameterspPb[0],ParameterspPb[1],ParameterspPb[2],ParameterspPb[3],ParameterspPb[4]); // standard

    graphInvYieldPi0CombpPb5023GeVTotClone->Fit(fitBylinkinCombPi0pPb5023GeVPt,"QVNRMEX0+","",minPt,maxPt);
    graphInvYieldPi0CombpPb5023GeVTotClone->Fit(fitBylinkinCombPi0pPb5023GeVPt,"QVNRMEX0+","",minPt,maxPt);
    
    histo2DCompCombined->DrawCopy();
    
     
    
    //DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVSysClone,20,1, kBlue, kBlue, 1, kTRUE, kBlue);  
     DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVSysClone,20,0.7, 4, 4, 1, kTRUE);  
      
    graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");

    DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVStaClone,20,0.7, kBlue, kBlue);  
    graphInvYieldPi0CombpPb5023GeVStaClone->Draw("pe,same");

     
    TLegend* legendRpPbCombine = new TLegend(0.23,0.30,0.65,0.40);
    legendRpPbCombine->SetFillColor(0);
    legendRpPbCombine->SetLineColor(0);
    legendRpPbCombine->SetTextSize(0.03);
    legendRpPbCombine->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
    legendRpPbCombine->AddEntry(fitBylinkinCombPi0pPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
    legendRpPbCombine->Draw();
	
    
    
    
    fitBylinkinCombPi0pPb5023GeVPt->Draw("same");
    
    canvasCompYieldpPbInd->Print(Form("%s/Comb_pPb.pdf",outputDir.Data()));

    // **************************************************************************************
    // ************************* Creating ratios of difffernt systems ************************
    // **************************************************************************************

    graphRatioCombPCMSys        = CalculateGraphErrRatioToFit(graphRatioCombPCMSys, fitBylinkinCombPi0pPb5023GeVPt);
    graphRatioCombDalitzSys     = CalculateGraphErrRatioToFit(graphRatioCombDalitzSys, fitBylinkinCombPi0pPb5023GeVPt);
    graphRatioCombEMCALSys      = CalculateGraphErrRatioToFit(graphRatioCombEMCALSys, fitBylinkinCombPi0pPb5023GeVPt);
    
    graphRatioCombCombFit       = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVTotClone,fitBylinkinCombPi0pPb5023GeVPt); 
    graphRatioCombCombFitSta    = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVStaClone,fitBylinkinCombPi0pPb5023GeVPt); 
    graphRatioCombCombFitSys    = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVSysClone,fitBylinkinCombPi0pPb5023GeVPt); 
    //graphRatioCombEMCALSys      = CalculateGraphErrRatioToFit(graphRatioCombEMCALSys,fitBylinkinCombPi0pPb5023GeVPt);

    histoRatioCombPHOS          = CalculateHistoRatioToFit(histoRatioCombPHOS,fitBylinkinCombPi0pPb5023GeVPt);
    

    TH1D* histoRatioCombPCM     = CalculateHistoRatioToFit(histoPCMYieldPi0pPb,fitBylinkinCombPi0pPb5023GeVPt);
    TH1D* histoRatioCombDalitz  = CalculateHistoRatioToFit(histoDalitzYieldPi0pPb,fitBylinkinCombPi0pPb5023GeVPt);
    TH1D* histoRatioCombEMCAL   = CalculateHistoRatioToFit(histoEMCALYieldPi0pPbStat,fitBylinkinCombPi0pPb5023GeVPt);
    TH1D* histoRatioCombPHOSsys = CalculateHistoRatioToFit(histoPHOSYieldPi0pPbSyst,fitBylinkinCombPi0pPb5023GeVPt);

    // **************************************************************************************
    // ************************* Plotting ratio of difffernt systems ************************
    // **************************************************************************************

    TCanvas* canvasRatioCompYieldpPbInd = new TCanvas("canvasRatioCompYieldpPbInd","",200,10,900,700);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioCompYieldpPbInd,  0.12, 0.02, 0.02, 0.12);
    canvasRatioCompYieldpPbInd->SetLogx();
//     canvasRatioCompYieldpPbInd->SetGridx();
//     canvasRatioCompYieldpPbInd->SetGridy();


    DrawPad(0,0.0,0.3,0.12,0.05);
    TH2F * ratio2DInvXSectionPi0;
    ratio2DInvXSectionPi0 = new TH2F("ratio2DInvXSectionPi0","ratio2DInvXSectionPi0",1000,0.3,25.,1000,0.51,1.49);
    
    //SetStyleHistoTH1ForGraphs(ratio2DInvXSectionPi0, "#it{p}_{T} (GeV/#it{c})", "Data/Fit",12,14,12,13,4.85, 1.5, 505,505);
    
    //#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}
    
    SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.09,0.11,0.07,0.1,1.,0.45,505,505);
    ratio2DInvXSectionPi0->GetYaxis()->CenterTitle();
	
	    
    /* 
    ratio2DInvXSectionPi0->GetYaxis()->SetLabelSize(0.055) ;
    ratio2DInvXSectionPi0->GetYaxis()->SetTitleSize(0.1) ;
    ratio2DInvXSectionPi0->GetYaxis()->SetTitleOffset(0.3) ;
    ratio2DInvXSectionPi0->SetYTitle("Data/Fit");
    ratio2DInvXSectionPi0->GetXaxis()->SetLabelSize(0.055) ;
    ratio2DInvXSectionPi0->GetXaxis()->SetTitleSize(0.1) ;
    ratio2DInvXSectionPi0->SetXTitle("#it{p}_{T}(GeV/#it{c})");
    */

    ratio2DInvXSectionPi0->DrawCopy();  
    graphRatioCombCombFitSys->SetPointEYlow(30,0.137146);//attribute EMCal Sys error to 
    graphRatioCombCombFitSys->SetPointEYhigh(30,0.137146);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys,20,1, kBlue, kBlue, 1, kTRUE, kBlue-9);  
        graphRatioCombCombFitSys->Draw("E2,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit,20,1, kBlue, kBlue);  
        graphRatioCombCombFit->Draw("same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPCMSys,24,1, 1, 1, 1, kTRUE);  
        graphRatioCombPCMSys->Draw("E2same");

        DrawGammaSetMarker(histoRatioCombPCM, 24,1, 1 , 1);
        histoRatioCombPCM->Draw("same") ;

        TLatex * lt2 = new TLatex(2.,1.25,"PCM") ;
        lt2->SetTextColor(kBlack) ;
	lt2->SetTextSize(0.11) ;
        lt2->Draw() ;

    //=============Dalitz=====================  
    canvasRatioCompYieldpPbInd->cd() ;
    DrawPad(1,0.0,0.,0.12,0.05);
    SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.1,0.15,0.1,0.15,0.8,0.3,505,505);
        ratio2DInvXSectionPi0->DrawCopy(); 
        graphRatioCombCombFitSys->Draw("E2,same") ;

        graphRatioCombCombFit->Draw("pe,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombDalitzSys,24,1, kCyan+2, kCyan+2, 1, kTRUE);  
        graphRatioCombDalitzSys->Draw("E2same");
        DrawGammaSetMarker(histoRatioCombDalitz, 24,1 ,kCyan+2 ,kCyan+2);
        histoRatioCombDalitz->Draw("same") ; 
	TLatex * lt = new TLatex(2.,0.7,"PCM"); 
	lt->SetTextSize(0.16) ;
      lt->SetTextColor(kCyan+2) ;
        lt->DrawText(2.,1.25,"Dalitz") ;
 
       // ===========EMCAL======================  
      canvasRatioCompYieldpPbInd->cd() ;
      DrawPad(2,0.0,0.,0.12,0.05);
      ratio2DInvXSectionPi0->DrawCopy(); 
  
      graphRatioCombCombFitSys->Draw("E2") ;
  
      graphRatioCombCombFit->Draw("pe,same");

            
      DrawGammaSetMarkerTGraphAsym(graphRatioCombEMCALSys,24,1,kGreen+2 , kGreen+2, 1, kTRUE);  
      graphRatioCombEMCALSys->Draw("E2same");
      DrawGammaSetMarker(histoRatioCombEMCAL, 24,1 ,kGreen+2 ,kGreen+2);
      histoRatioCombEMCAL->Draw("same") ;
      cout<<"EMCal"<<endl;
      graphRatioCombEMCALSys->Print();
      lt->SetTextColor(kGreen+2) ;
      lt->DrawText(2.,1.25,"EMCal") ;

      //PHOS  
      canvasRatioCompYieldpPbInd->cd() ;
      DrawPad(3,0.,0.,0.12,0.05);
  
      ratio2DInvXSectionPi0->DrawCopy(); 
      graphRatioCombCombFitSys->Draw("E2") ;
      graphRatioCombCombFit->Draw("pe,same");
      graphRatioCombCombFit->Print();
      graphRatioCombCombFitSys->Print();
      histoRatioCombPHOSsys->SetLineColor(kRed+1);
      histoRatioCombPHOSsys->SetFillStyle(0) ;
      histoRatioCombPHOSsys->Draw("E2same") ;
  
      histoRatioCombPHOS->SetMarkerColor(kRed+1);
      histoRatioCombPHOS->SetLineColor(kRed+1);
      histoRatioCombPHOS->Draw("same");


      lt->SetTextColor(kRed+1) ;
      lt->DrawText(2.,1.25,"PHOS") ;

      canvasRatioCompYieldpPbInd->Update();
      canvasRatioCompYieldpPbInd->Print(Form("%s/RatioCompYieldpPb.pdf",outputDir.Data()));
    // **************************************************************************************
    // ************************* Plotting ratio combined to fit ************************
    // **************************************************************************************

    TCanvas* canvasRatioCompYieldpPb = new TCanvas("canvasRatioCompYieldpPb","",200,10,900,700);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioCompYieldpPb,  0.12, 0.02, 0.02, 0.12);
    canvasRatioCompYieldpPb->SetLogx();
//     canvasRatioCompYieldpPb->SetGridx();
//     canvasRatioCompYieldpPb->SetGridy();


    TH2F * ratio2DInvXSectionPi0a;
    ratio2DInvXSectionPi0a = new TH2F("ratio2DInvXSectionPi0a","ratio2DInvXSectionPi0a",1000,0.3,25.,1000,0.21,2.39);

    //SetStyleHistoTH1ForGraphs(ratio2DInvXSectionPi0a, "#it{p}_{T} (GeV/#it{c})", "Data/Fit",12,14,12,13,4.85, 1.5, 505,505);
    
    //#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}
    
    SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0a, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
    //  ratio2DInvXSectionPi0a->GetYaxis()->CenterTitle();
	
	    
    /* 
    ratio2DInvXSectionPi0a->GetYaxis()->SetLabelSize(0.055) ;
    ratio2DInvXSectionPi0a->GetYaxis()->SetTitleSize(0.1) ;
    ratio2DInvXSectionPi0a->GetYaxis()->SetTitleOffset(0.3) ;
    ratio2DInvXSectionPi0a->SetYTitle("Data/Fit");
    ratio2DInvXSectionPi0a->GetXaxis()->SetLabelSize(0.055) ;
    ratio2DInvXSectionPi0a->GetXaxis()->SetTitleSize(0.1) ;
    ratio2DInvXSectionPi0a->SetXTitle("#it{p}_{T}(GeV/#it{c})");
    */

    ratio2DInvXSectionPi0a->DrawCopy();  

	TLine *lineA=new TLine(0.,1.,25.,1.);
	lineA->SetLineColor(kGray+1);
	lineA->Draw();
        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys,20,1., kBlue, kBlue, 1.5, kTRUE, 0);  
        graphRatioCombCombFitSys->Draw("E2,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit,20,1., kBlue, kBlue,1.5);  
        graphRatioCombCombFit->Draw("p,same");

        TLatex * lt3 = new TLatex(3.,2.2,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
        lt3->SetTextColor(kBlack) ;
	lt3->SetTextSize(0.05) ;
	lt3->DrawLatex(3.,2.0,"ALICE");
	lt3->DrawLatex(3.,1.8,"#pi^{0} #rightarrow #gamma#gamma");
        lt3->Draw() ;

 

      canvasRatioCompYieldpPb->Update();
      canvasRatioCompYieldpPb->Print(Form("%s/RatioCombYieldFitpPb.pdf",outputDir.Data()));


   // **************************************************************************************
    // ************************* Plotting ratio individual spectra to fit ************************
    // **************************************************************************************

    TCanvas* canvasRatioIndYieldpPb = new TCanvas("canvasRatioIndYieldpPb","",200,10,900,700);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioIndYieldpPb,  0.12, 0.02, 0.02, 0.12);
    canvasRatioIndYieldpPb->SetLogx();
//     canvasRatioIndYieldpPb->SetGridx();
//     canvasRatioIndYieldpPb->SetGridy();


    TH2F * ratio2DInvXSectionPi0b;
    ratio2DInvXSectionPi0b = new TH2F("ratio2DInvXSectionPi0b","ratio2DInvXSectionPi0b",1000,0.3,25.,1000,0.21,2.39);

    //SetStyleHistoTH1ForGraphs(ratio2DInvXSectionPi0b, "#it{p}_{T} (GeV/#it{c})", "Data/Fit",12,14,12,13,4.85, 1.5, 505,505);
    
    //#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}
    
    SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0b, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
    //  ratio2DInvXSectionPi0b->GetYaxis()->CenterTitle();
	
	    
    /* 
    ratio2DInvXSectionPi0b->GetYaxis()->SetLabelSize(0.055) ;
    ratio2DInvXSectionPi0b->GetYaxis()->SetTitleSize(0.1) ;
    ratio2DInvXSectionPi0b->GetYaxis()->SetTitleOffset(0.3) ;
    ratio2DInvXSectionPi0b->SetYTitle("Data/Fit");
    ratio2DInvXSectionPi0b->GetXaxis()->SetLabelSize(0.055) ;
    ratio2DInvXSectionPi0b->GetXaxis()->SetTitleSize(0.1) ;
    ratio2DInvXSectionPi0b->SetXTitle("#it{p}_{T}(GeV/#it{c})");
    */

    ratio2DInvXSectionPi0b->DrawCopy();  

	TLine *lineB=new TLine(0.,1.,25.,1.);
	lineB->SetLineColor(kGray+1);
	lineB->Draw();


	DrawGammaSetMarkerTGraphAsym(graphRatioCombPCMSys,20,1, 1, 1, 1, kTRUE);  
        graphRatioCombPCMSys->Draw("E2same");
        DrawGammaSetMarker(histoRatioCombPCM, 20,1, 1 , 1);
        histoRatioCombPCM->Draw("same") ;

        DrawGammaSetMarkerTGraphAsym(graphRatioCombDalitzSys,29,1, kCyan+2, kCyan+2, 1, kTRUE);  
        graphRatioCombDalitzSys->Draw("E2same");
        DrawGammaSetMarker(histoRatioCombDalitz, 29,1.5 ,kCyan+2 ,kCyan+2);
        histoRatioCombDalitz->Draw("same") ; 
 
	DrawGammaSetMarkerTGraphAsym(graphRatioCombEMCALSys,33,1,kGreen+2 , kGreen+2, 1, kTRUE);  
	graphRatioCombEMCALSys->Draw("E2same");
	DrawGammaSetMarker(histoRatioCombEMCAL, 33,1.5 ,kGreen+2 ,kGreen+2);
	histoRatioCombEMCAL->Draw("same") ;
	
	histoRatioCombPHOSsys->Draw("E2same") ;
	
	histoRatioCombPHOS->SetMarkerStyle(21);
	histoRatioCombPHOS->SetMarkerSize(1);
	histoRatioCombPHOS->Draw("same");
	
        TLatex * lt4 = new TLatex(3.,2.2,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
        lt4->SetTextColor(kBlack) ;
	lt4->SetTextSize(0.05) ;
	lt4->DrawLatex(3.,2.,"ALICE");
	lt4->DrawLatex(3.,1.8,"#pi^{0} #rightarrow #gamma#gamma");
        lt4->Draw() ;

 	TLegend* legendSpectraDiffDetMinBiasStrip = new TLegend(0.18,0.7,0.5,0.95);
	legendSpectraDiffDetMinBiasStrip->SetFillColor(0);
	legendSpectraDiffDetMinBiasStrip->SetLineColor(0);
	legendSpectraDiffDetMinBiasStrip->SetTextFont(42);
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombPCM,Form("PCM"),"pf");
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombDalitz,Form("Dalitz"),"pf");
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombEMCAL,Form("EMCal"),"pf");
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombPHOS,Form("PHOS"),"pf");
	legendSpectraDiffDetMinBiasStrip->Draw();

      canvasRatioIndYieldpPb->Update();
      canvasRatioIndYieldpPb->Print(Form("%s/RatioIndYieldFitpPb.pdf",outputDir.Data()));


    TFile fResults(Form("%s/ResultspPb_%s.root",outputDir.Data(), dateForOutput.Data()),"RECREATE");
    fitBylinkinCombPi0pPb5023GeVPt->Write("FitCombinedpPbSpectrum");
      graphCombPi0InvCrossSectionStatpPb5023GeV->Write("CombinedpPbSpectrumStatErr");
      graphCombPi0InvCrossSectionSyspPb5023GeV->Write("CombinedpPbSpectrumSysErr");
      graphCombPi0InvCrossSectionTotpPb5023GeV->Write("CombinedpPbSpectrumTotErr");



}
    
