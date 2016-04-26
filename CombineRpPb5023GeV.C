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

void CombineRpPb5023GeV(){    

    TString date                                        = ReturnDateString();
    
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();
    SetPlotStyle();
    
    TString suffix                                      = "eps";
    
    TString dateForOutput                               = ReturnDateStringForOutput();
    TString outputDir                                   = Form("CombineRpPb/%s/%s",suffix.Data(),dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    
    
    //cout << dateForOutput.Data() << endl;
    //___________________________________ Declaration of files _____________________________________________

    
 
    TString fileNameRpPb                   = "ExternalInputpPb/ResultsRpPb_2016_04_06.root";
    TFile* fileNeutralPionRpPb                           = new TFile(fileNameRpPb.Data());
    
    TString nameHistoPCM                               = "Pi0_RpPb_PCM_StatErr";
    TString nameHistoPCMSysErrors                      = "Pi0_RpPb_PCM_SystErr";
    TString nameHistoDalitz                            = "Pi0_RpPb_Dalitz_StatErr";
    TString nameHistoDalitzSysErrors                   = "Pi0_RpPb_Dalitz_SystErr";
    TString nameHistoPHOS                              = "Pi0_RpPb_PHOS_StatErr";
    TString nameHistoPHOSSysErrors                     = "Pi0_RpPb_PHOS_SystErr";
    TString nameHistoEMCal                             = "Pi0_RpPb_EMCal_StatErr";
    TString nameHistoEMCalSysErrors                    = "Pi0_RpPb_EMCal_StatErr";
    
    
    TString collisionSystempPb                          = "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"; 
    
    
    Double_t pTLimits[nPtLimits]                        = { 0.3, 0.4, 0.5, 0.6, 0.7, 
                                                            0.8, 1.0, 1.2, 1.4, 1.6, 
                                                            1.8, 2.0, 2.2, 2.4, 2.6,
                                                            2.8, 3.0, 3.2, 3.4, 3.6,
                                                            3.8, 4.0, 4.5, 5.0, 5.5,
                                                            6.0, 7.0, 8.0, 10.0,12.0,
							    14.0,16.0, 20.0};
                                
    
    Int_t offSets[11]                                   =  { 0, 6, -1, 0, 0,  3, 0, 0, 0,  0, 0};
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

    TGraphAsymmErrors* graphDalitzYieldPi0pPb           = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get(nameHistoDalitz.Data());
    TH1D* histoDalitzYieldPi0pPb           = (TH1D*)GraphAsymErrorsToHist_withErrors(graphDalitzYieldPi0pPb,nameHistoDalitz.Data());
    TGraphAsymmErrors* graphDalitzYieldPi0pPbSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get(nameHistoDalitzSysErrors.Data());

    // **************************************************************************************
    // ****************************** Reading PCM *******************************************
    // **************************************************************************************    
    TGraphAsymmErrors* graphPCMYieldPi0pPb           = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get(nameHistoPCM.Data());
    TH1D* histoPCMYieldPi0pPb           = (TH1D*)GraphAsymErrorsToHist_withErrors(graphPCMYieldPi0pPb,nameHistoPCM.Data());
    TGraphAsymmErrors* graphPCMYieldPi0pPbSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get(nameHistoPCMSysErrors.Data());
    
    // **************************************************************************************
    // ******************************* Reading PHOS *****************************************
    // **************************************************************************************    
  
    TGraphAsymmErrors* graphPHOSYieldPi0pPb                      = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get(nameHistoPHOS.Data()); 
        TH1D* histoPHOSYieldPi0pPb           = GraphAsymErrorsToHist_withErrors(graphPHOSYieldPi0pPb,nameHistoPHOS.Data());
    TGraphAsymmErrors* graphPHOSYieldPi0pPbSystErr                   = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get(nameHistoPHOSSysErrors.Data());

    // **************************************************************************************
    // ******************************** Reading EMCal ***************************************
    // **************************************************************************************
    
    //   TGraphAsymmErrors* graphEMCALYieldPi0pPb       = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get(nameHistoEMCal.Data());
    //    TH1D* histoEMCALYieldPi0pPb           = GraphAsymErrorsToHist_withErrors(graphEMCALYieldPi0pPb,20.,nameHistoEMCal.Data());
    //   TGraphAsymmErrors* graphEMCALYieldPi0pPbSystErr     = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get(nameHistoEMCalSysErrors.Data());

    // **************************************************************************************
    // ******************************** Reading Model calculations **************************
    // **************************************************************************************
    TGraphAsymmErrors*  graphAsymmErrorsPi0DSS5000    = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_fDSS_errors"); 
    TGraph*  graphPi0CGC                     = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("ColorGlasCondensate"); 
    TGraph*  graphPi0DSS5000                 = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_fDSS_NLO"); 
    TGraph*  graphPi0ESP09sPi0AKK            = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_AKK_NLO"); 
    TGraph*  graphPi0ESP09sPi0KKP            = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_KKP_NLO"); 


	graphAsymmErrorsPi0DSS5000->Draw("same,E3");
	graphPi0CGC->Draw("same,p");
	graphPi0DSS5000->Draw("same,p,l");
	graphPi0ESP09sPi0AKK->Draw("same,p,l");
	graphPi0ESP09sPi0KKP->Draw("same,p,l");
    // **************************************************************************************
    // ******************************** Reading Charged Particles and Pions **************************
    // **************************************************************************************
    TGraphAsymmErrors*  graphRpPbChargedParticlesStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("RpPb_ChargedParticles_StatErr"); 
    TGraphAsymmErrors*  graphRpPbChargedParticlesSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("RpPb_ChargedParticles_SystErr"); 
    TGraphAsymmErrors*  graphRpPbChargedPionsStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("RpPb_ChargedPions_StatErr"); 
    TGraphAsymmErrors*  graphRpPbChargedPionsSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("RpPb_ChargedPions_SystErr"); 
     // **************************************************************************************
    // ********************************* Combine spectra ************************************
    // **************************************************************************************
    
    TString fileNameOutputWeightingOld                  = Form("%s/WeightingOld.dat",outputDir.Data());

    statErrorCollection[0]          = (TH1D*)histoPCMYieldPi0pPb->Clone("statErrPCMPi0");
    statErrorCollection[1]          = (TH1D*)histoPHOSYieldPi0pPb->Clone("statErrPHOSPi0");
    //   statErrorCollection[2]          = (TH1D*)histoEMCALYieldPi0pPb->Clone("statErrEMCALPi0");
    statErrorCollection[5]          = (TH1D*)histoDalitzYieldPi0pPb->Clone("statErrDalitzPi0");
    
    sysErrorCollection[0]           = (TGraphAsymmErrors*)graphPCMYieldPi0pPbSystErr->Clone("sysErrPCMPi0");
    sysErrorCollection[1]           = (TGraphAsymmErrors*)graphPHOSYieldPi0pPbSystErr->Clone("sysErrPHOSPi0");
    //    sysErrorCollection[2]           = (TGraphAsymmErrors*)graphEMCALYieldPi0pPbSystErr->Clone("sysErrEMCALPi0");
    sysErrorCollection[5]           = (TGraphAsymmErrors*)graphDalitzYieldPi0pPbSystErr->Clone("sysErrDalitzPi0");
    cout << "Sys Error PCM:" <<endl;
    graphPCMYieldPi0pPbSystErr->Print();    cout << "Sys Error PHOS:" <<endl;
    graphPHOSYieldPi0pPbSystErr->Print();  cout << "Sys Error Dalitz:" <<endl;
    graphDalitzYieldPi0pPbSystErr->Print();
    
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
    // ************************* Plotting Combined R_pPb ************************
    // **************************************************************************************
    TCanvas* canvasCombRpPb = new TCanvas("canvasCombRpPb","",200,10,1200,700);  // gives the page size
    DrawGammaCanvasSettings( canvasCombRpPb, 0.08, 0.02, 0.02, 0.13);
    
    //  canvasCombRpPb->SetLogx();
    // canvasCombRpPb->SetLogy();
    TH2F * histo2DCombined;
    histo2DCombined = new TH2F("histo2DCombined","histo2DCombined",1000,0.,20.,1000,0.3,1.5);
     SetStyleHistoTH2ForGraphs(histo2DCombined, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
    histo2DCombined->DrawCopy();
    

     
    

     DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVSysClone,20,1.5, 4, 4, 1, kTRUE);  
      
    graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");

    DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVStaClone,20,1.5, kBlue, kBlue);  
    graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");

     
    TLegend* legendRpPbCombine = new TLegend(0.1,0.85,0.55,0.95);
    legendRpPbCombine->SetFillColor(0);
    legendRpPbCombine->SetLineColor(0);
    legendRpPbCombine->SetTextSize(0.03);
    legendRpPbCombine->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");

    
    legendRpPbCombine->Draw();
  
    TLine* line =new TLine(0.,1.,20.,1.);
    line->Draw("same");
 
    
     canvasCombRpPb->Print(Form("%s/Comb_RpPb.%s",outputDir.Data(),suffix.Data()));
    
    
    // **************************************************************************************
    // ************************* Plotting All GA R_pPb ************************
    // **************************************************************************************
    TCanvas* canvasAllGARpPb = new TCanvas("canvasAllGARpPb","",200,10,1200,700);  // gives the page size
    DrawGammaCanvasSettings( canvasAllGARpPb, 0.08, 0.02, 0.02, 0.13);
    
    //  canvasAllGARpPb->SetLogx();
    // canvasAllGARpPb->SetLogy();
    TH2F * histo2DAllGA;
    histo2DAllGA = new TH2F("histo2DAllGA","histo2DAllGA",1000,0.,20.,1000,0.3,1.5);
     SetStyleHistoTH2ForGraphs(histo2DAllGA, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
    histo2DAllGA->DrawCopy();
    
    graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
    graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");

    graphDalitzYieldPi0pPbSystErr->Draw("E2,same");  
    graphDalitzYieldPi0pPb->Draw("p,same");  

    graphPCMYieldPi0pPbSystErr->Draw("E2,same");  
    graphPCMYieldPi0pPb->Draw("p,same");  

    graphPHOSYieldPi0pPbSystErr->Draw("E2,same");  
    graphPHOSYieldPi0pPb->Draw("p,same");  



    TLegend* legendRpPbAllGA = new TLegend(0.2,0.15,0.5,0.35);
    legendRpPbAllGA->SetFillColor(0);
    legendRpPbAllGA->SetLineColor(0);
    legendRpPbAllGA->SetTextSize(0.03);
    legendRpPbAllGA->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"Combined","pef");
    legendRpPbAllGA->AddEntry(graphPCMYieldPi0pPbSystErr,"PCM","pef");
    legendRpPbAllGA->AddEntry(graphDalitzYieldPi0pPbSystErr,"Dalitz","pef");
    legendRpPbAllGA->AddEntry(graphPHOSYieldPi0pPbSystErr,"PHOS","pef");

    
    legendRpPbAllGA->Draw();
  

    line->Draw("same");
 
    
     canvasAllGARpPb->Print(Form("%s/AllGA_RpPb.%s",outputDir.Data(),suffix.Data()));   

 
    // **************************************************************************************
    // ************************* Plotting Combined R_pPb with models ************************
    // **************************************************************************************
    TCanvas* canvasCombinedWithModelsRpPb = new TCanvas("canvasCombinedWithModelsRpPb","",200,10,1200,700);  // gives the page size
    DrawGammaCanvasSettings( canvasCombinedWithModelsRpPb, 0.08, 0.02, 0.02, 0.13);
    
    //  canvasCombinedWithModelsRpPb->SetLogx();
    // canvasCombinedWithModelsRpPb->SetLogy();
    TH2F * histo2DCombinedWithModels;
    histo2DCombinedWithModels = new TH2F("histo2DCombinedWithModels","histo2DCombinedWithModels",1000,0.,20.,1000,0.3,1.5);
     SetStyleHistoTH2ForGraphs(histo2DCombinedWithModels, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
    histo2DCombinedWithModels->DrawCopy();
    
    graphAsymmErrorsPi0DSS5000->Draw("same,E3");
	graphPi0CGC->Draw("same,p");
	graphPi0DSS5000->Draw("same,p,l");
	graphPi0ESP09sPi0AKK->Draw("same,p,l");
	graphPi0ESP09sPi0KKP->Draw("same,p,l");
    graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
    graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");
        Float_t offsetYMaxLegendPCM = -0.15;
	Float_t xMinTheoryPCM = 0.30;

	TLegend* legendRpPbESP09sPCM = new TLegend(0.2,0.15,0.45,0.43);
	legendRpPbESP09sPCM->SetFillColor(0);
	legendRpPbESP09sPCM->SetLineColor(0);
	legendRpPbESP09sPCM->SetTextSize(0.04);
	legendRpPbESP09sPCM->SetTextFont(42);
	legendRpPbESP09sPCM->AddEntry(graphPi0ESP09sPi0KKP,"EPS09s KKP NLO");
	legendRpPbESP09sPCM->AddEntry(graphPi0ESP09sPi0AKK,"EPS09s AKK NLO");
	legendRpPbESP09sPCM->AddEntry(graphPi0DSS5000,"EPS09s fDSS NLO");
	legendRpPbESP09sPCM->AddEntry(graphAsymmErrorsPi0DSS5000,"EPS09s fDSS errors","pef");
	
	
	TLegend* legendRpPbCGCPCM = new TLegend(0.5,0.15,0.85,0.25);
	legendRpPbCGCPCM->SetFillColor(0);
	legendRpPbCGCPCM->SetLineColor(0);
	legendRpPbCGCPCM->SetTextSize(0.04);
	legendRpPbCGCPCM->SetTextFont(42);
	legendRpPbCGCPCM->AddEntry(graphPi0CGC,"Color Glas Condensate","p");
    	legendRpPbESP09sPCM->Draw("same");
	legendRpPbCGCPCM->Draw("same");

  

    line->Draw("same");
 
    
     canvasCombinedWithModelsRpPb->Print(Form("%s/CombinedWithModels_RpPb.%s",outputDir.Data(),suffix.Data()));
 


    // **************************************************************************************
    // ************************* Plotting Combined R_pPb with Charged Particles ************************
    // **************************************************************************************
    TCanvas* canvasCombinedWithChargedParticlesRpPb = new TCanvas("canvasCombinedWithChargedParticlesRpPb","",200,10,1200,700);  // gives the page size
    DrawGammaCanvasSettings( canvasCombinedWithChargedParticlesRpPb, 0.08, 0.02, 0.02, 0.13);
    
    //  canvasCombinedWithChargedParticlesRpPb->SetLogx();
    // canvasCombinedWithChargedParticlesRpPb->SetLogy();
    TH2F * histo2DCombinedWithChargedParticles;
    histo2DCombinedWithChargedParticles = new TH2F("histo2DCombinedWithChargedParticles","histo2DCombinedWithChargedParticles",1000,0.,20.,1000,0.3,1.5);
     SetStyleHistoTH2ForGraphs(histo2DCombinedWithChargedParticles, "#it{p}_{T} (GeV/#it{c})","#it{R}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
    histo2DCombinedWithChargedParticles->DrawCopy();
    
    graphRpPbChargedParticlesSystErr->Draw("E2,same");  
    graphRpPbChargedParticlesStatErr->Draw("p,same");  
    graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
    graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");
  
    TLegend* legendRpPbChargedParticles = new TLegend(0.1,0.85,0.55,0.95);
    legendRpPbChargedParticles->SetFillColor(0);
    legendRpPbChargedParticles->SetLineColor(0);
    legendRpPbChargedParticles->SetTextSize(0.03);
    legendRpPbChargedParticles->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"#pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
    legendRpPbChargedParticles->AddEntry(graphRpPbChargedParticlesSystErr,"charged particles, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
   
    legendRpPbChargedParticles->Draw();

    line->Draw("same");
 
    
     canvasCombinedWithChargedParticlesRpPb->Print(Form("%s/CombinedWithChargedParticles_RpPb.%s",outputDir.Data(),suffix.Data()));


    // **************************************************************************************
    // ************************* Plotting Combined R_pPb with Charged Pions ************************
    // **************************************************************************************
    TCanvas* canvasCombinedWithChargedPionsRpPb = new TCanvas("canvasCombinedWithChargedPionsRpPb","",200,10,1200,700);  // gives the page size
    DrawGammaCanvasSettings( canvasCombinedWithChargedPionsRpPb, 0.08, 0.02, 0.02, 0.13);
    
    //  canvasCombinedWithChargedPionsRpPb->SetLogx();
    // canvasCombinedWithChargedPionsRpPb->SetLogy();
    TH2F * histo2DCombinedWithChargedPions;
    histo2DCombinedWithChargedPions = new TH2F("histo2DCombinedWithChargedPions","histo2DCombinedWithChargedPions",1000,0.,20.,1000,0.3,1.5);
     SetStyleHistoTH2ForGraphs(histo2DCombinedWithChargedPions, "#it{p}_{T} (GeV/#it{c})","#it{R}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
    histo2DCombinedWithChargedPions->DrawCopy();
    
    graphRpPbChargedPionsSystErr->Draw("E2,same");  
    graphRpPbChargedPionsStatErr->Draw("p,same");  
    graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
    graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");
  
    TLegend* legendRpPbChargedPions = new TLegend(0.1,0.85,0.55,0.95);
    legendRpPbChargedPions->SetFillColor(0);
    legendRpPbChargedPions->SetLineColor(0);
    legendRpPbChargedPions->SetTextSize(0.03);
    legendRpPbChargedPions->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"#pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
    legendRpPbChargedPions->AddEntry(graphRpPbChargedPionsSystErr,"#pi^{+/-}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
   
    legendRpPbChargedPions->Draw();

    line->Draw("same");
 
    
     canvasCombinedWithChargedPionsRpPb->Print(Form("%s/CombinedWithChargedPions_RpPb.%s",outputDir.Data(),suffix.Data()));


}
    
