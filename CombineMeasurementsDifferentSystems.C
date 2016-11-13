/****************************************************************************************************************************
******      provided by Gamma Conversion Group, PWGGA,                                                                  *****
******      Friederike Bock, friederike.bock@cern.ch                                                                    *****
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

extern TRandom* gRandom;
extern TBenchmark*  gBenchmark;
extern TSystem* gSystem;
extern TMinuit*     gMinuit;

//*********************************************************************************************************************************
//**************************************** Main function for combination macro ****************************************************
//*********************************************************************************************************************************
void CombineMeasurementsDifferentSystems(       TString fileNameCombMesonPbPb   = "", 
                                                TString fileNameCombGammaPbPb   = "", 
                                                TString fileNameCombMesonpPb    = "", 
                                                TString suffix                  = "eps"){  

    //******************************** general style settings ***********************************************
    gROOT->Reset(); 
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();  
    SetPlotStyle();

    //******************************** label declaration ****************************************************
    TString collisionSystemPbPb0010         = "0-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";        
    TString collisionSystemPbPb0020         = "0-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";        
    TString collisionSystempPb              = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";        
    TString dateForOutput                   = ReturnDateStringForOutput();
    TString date                            = ReturnDateString();
    
    //******************************** set common colors/markers ********************************************
    Color_t colorComb0010PbPb               = GetColorDefaultColor( "PbPb_2.76TeV", "", "0-10%", kFALSE ); 
    Color_t colorComb0010PbPbBox            = GetColorDefaultColor( "PbPb_2.76TeV", "", "0-10%", kTRUE ); 
    Color_t colorComb0005PbPb               = GetColorDefaultColor( "PbPb_2.76TeV", "", "0-5%", kFALSE ); 
    Color_t colorChPi0005PbPb               = GetColorDefaultColor( "PbPb_2.76TeV", "MC", "0-5%", kFALSE ); 
    Color_t colorComb0005PbPbBox            = GetColorDefaultColor( "PbPb_2.76TeV", "", "0-5%", kTRUE ); 
    Color_t colorComb0020PbPb               = GetColorDefaultColor( "PbPb_2.76TeV", "", "10-20%", kFALSE );
    Color_t colorComb0020PbPbBox            = GetColorDefaultColor( "PbPb_2.76TeV", "", "10-20%", kTRUE ); 
    Style_t markerStyleComb0010PbPb         = GetDefaultMarkerStyle( "PbPb_2.76TeV", "", "0-10%" ); 
    Style_t markerStyleComb0005PbPb         = GetDefaultMarkerStyle( "PbPb_2.76TeV", "", "0-5%" ); 
    Style_t markerStyleChPi0005PbPb         = GetDefaultMarkerStyle( "PbPb_2.76TeV", "MC", "0-5%" ); 
    Style_t markerStyleCombGamma0020PbPb    = GetDefaultMarkerStyle( "PbPb_2.76TeV", "", "5-10%" ); 
    Size_t markerSizeComb0010PbPb           = GetDefaultMarkerSize( "PbPb_2.76TeV", "", "0-10%" ); 
    Size_t markerSizeComb0005PbPb           = GetDefaultMarkerSize( "PbPb_2.76TeV", "", "0-5%" ); 
    Size_t markerSizeCombGamma0020PbPb      = GetDefaultMarkerSize( "PbPb_2.76TeV", "", "5-10%" ); 
    Color_t colorCombpPb                    = GetColorDefaultColor( "pPb_5.023TeV", "", "", kFALSE ); 
    Color_t colorCombpPbBox                 = GetColorDefaultColor( "pPb_5.023TeV", "", "", kTRUE ); 
    Style_t markerStyleCombpPb              = GetDefaultMarkerStyle( "pPb_5.023TeV", "", "" ); 
    Size_t markerSizeCombpPb                = GetDefaultMarkerSize( "pPb_5.023TeV", "", "" )*1.5; 
    
    Style_t markerStylePHENIX200GeV         = 25 ;
    Style_t markerStylePHENIX62GeV          = 27 ;
    Style_t markerStylePHENIX39GeV          = 24 ;
    Style_t markerStyleWA98                 = 28 ;
    Size_t markerSizePHENIX200GeV           = 1.95;
    Size_t markerSizePHENIX62GeV            = 3;
    Size_t markerSizePHENIX39GeV            = 1.95;
    Size_t markerSizeWA98                   = 1.95;
    Color_t colorVitevBas0005               = kRed-6;
    Color_t colorWHDG0005                   = kRed-4;
    Style_t fillStyleVitev                  = 3766;
    Style_t fillStyleWHDG                   = 3545;

    Width_t widthLinesBoxes                 = 2.3;
    if (suffix.CompareTo("eps")==0)
        widthLinesBoxes                     = 1.4;
    
    //********************************** xSections and errors ***************************************
    Double_t xSection2760GeVpp              = 55.416*1e-3;
    Double_t xSection2760GeVErrpp           = 3.9;
    Double_t xSection2769GeVppINEL          = 62.8*1e9;
    Double_t recalcBarn                     = 1e12; //NLO in pbarn!!!!
    Double_t commonCentralityErr0010        = 0.25;
    Double_t nCollErr0010                   = GetNCollErrFromName("0010");
    Double_t nCollErr0020                   = GetNCollErrFromName("0020");
    Double_t nColl0010                      = GetNCollFromName("0010");    
    Double_t nColl0020                      = GetNCollFromName("0020");    
    Double_t normErr0010                    = nCollErr0010/nColl0010;//pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0010/tAA0010),2)+pow((commonCentralityErr0010/100),2),0.5);
    Double_t normErr0020                    = nCollErr0020/nColl0020;
    Double_t normErrpPb                     = TMath::Sqrt(pow(0.031,2)+pow(0.036,2)+pow(0.036,2));//pPb normalization,TpPb and pp normalization errors

    //******************************** Declaration of files *****************************************
    TString fileNameDataOtherEnergyInput    = "ExternalInputPbPb/OtherExperiments/DataCompilationFromOtherEnergiesPbPb.root";
    TString fileNameTheoryInput             = "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";
    TString fileNameChargedPiRAA            = "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/RAA_Pion_08052014.root";
    
    TString outputDir                       = Form("%s/%s/CombineMeasurementsDifferentSystems",suffix.Data(),dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputCombMesonPbPb.root", fileNameCombMesonPbPb.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputCombGammaPbPb.root", fileNameCombGammaPbPb.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputCombMesonpPb.root", fileNameCombMesonpPb.Data(), outputDir.Data()));


    //******************************** Read PbPb meson file *****************************************
    TFile* fileFinalResultsMesonPbPb        = new TFile(fileNameCombMesonPbPb);
    TGraphAsymmErrors* graphPi0RAAPbPbStatComb0010      = (TGraphAsymmErrors*)fileFinalResultsMesonPbPb->Get("graphRAAStatErr_0010");
    TGraphAsymmErrors* graphPi0RAAPbPbSysComb0010       = (TGraphAsymmErrors*)fileFinalResultsMesonPbPb->Get("graphRAASysErr_0010");
    TGraphAsymmErrors* graphPi0RAAPbPbStatComb0010WX    = (TGraphAsymmErrors*)graphPi0RAAPbPbStatComb0010->Clone("graphRAAStatErr_0010_wXErr");
    ProduceGraphAsymmWithoutXErrors(graphPi0RAAPbPbStatComb0010WX);
    TGraphAsymmErrors* graphPi0RAAPbPbStatComb0005      = (TGraphAsymmErrors*)fileFinalResultsMesonPbPb->Get("graphRAAStatErr_0005");
    TGraphAsymmErrors* graphPi0RAAPbPbSysComb0005       = (TGraphAsymmErrors*)fileFinalResultsMesonPbPb->Get("graphRAASysErr_0005");
    TGraphAsymmErrors* graphPi0RAAPbPbStatComb0005WX    = (TGraphAsymmErrors*)graphPi0RAAPbPbStatComb0005->Clone("graphRAAStatErr_0005_wXErr");
    ProduceGraphAsymmWithoutXErrors(graphPi0RAAPbPbStatComb0005WX);
    
    //******************************** Read PbPb gamma file *****************************************
    TFile* fileFinalResultsGammaPbPb        = new TFile(fileNameCombGammaPbPb);
    TDirectory* dirGamma0020PbPb            = (TDirectory*)fileFinalResultsGammaPbPb->Get("Gamma_PbPb_2.76TeV_0-20%"); 
    TGraphAsymmErrors* graphGammaRAAPbPbStatComb0020    = (TGraphAsymmErrors*)dirGamma0020PbPb->Get("RAA_comb_StatErr");
    TGraphAsymmErrors* graphGammaRAAPbPbSysComb0020     = (TGraphAsymmErrors*)dirGamma0020PbPb->Get("RAA_comb_SysErr");
    TGraphAsymmErrors* graphGammaRAAPbPbStatComb0020WX  = (TGraphAsymmErrors*)graphGammaRAAPbPbStatComb0020->Clone("RAA_comb_StatErr_wXErr");
    ProduceGraphAsymmWithoutXErrors(graphGammaRAAPbPbStatComb0020WX);
        
    //******************************** Read pPb meson file ******************************************
    TFile* fileFinalResultsMesonpPb         = new TFile(fileNameCombMesonpPb);
    TGraphAsymmErrors* graphPi0RpApPbStatComb   = (TGraphAsymmErrors*)fileFinalResultsMesonpPb->Get("CombinedPi0RpPbStatErr");
    TGraphAsymmErrors* graphPi0RpApPbSysComb    = (TGraphAsymmErrors*)fileFinalResultsMesonpPb->Get("CombinedPi0RpPbSystErr");
    graphPi0RpApPbSysComb->RemovePoint(25);
    graphPi0RpApPbStatComb->RemovePoint(25);
    TGraphAsymmErrors* graphPi0RpApPbStatCombWX = (TGraphAsymmErrors*)graphPi0RpApPbStatComb->Clone("CombinedPi0RpPbStatErrr_wXErr");
    ProduceGraphAsymmWithoutXErrors(graphPi0RpApPbStatCombWX);
    
    //******************************** Read other energy file ******************************************
    TFile* fileDataOtherEnergies            = new TFile(fileNameDataOtherEnergyInput);
    TGraphErrors* graphWA98_17_3GeVRAA_0013 = (TGraphErrors*)fileDataOtherEnergies->Get("graphWA98RAA_0013");
    TGraphErrors* graphPHENIX200GeVRAA_0010 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_0010");
    TGraphErrors* graphPHENIX39GeVRAA_0010  = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_0010");
    TGraphErrors* graphPHENIX62GeVRAA_0010  = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_0010");


    //******************************** Theory input files **********************************************
    TFile* fileTheoryGraphs                 = new TFile(fileNameTheoryInput);
    TGraphErrors* Vitev_Bas_Raa_0005        = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA0005");
//     TGraphErrors* Xiao_Raa_0005             = (TGraphErrors*)Xiao_Raa_0020->Clone("Xiao_Raa_0005");
    TGraphAsymmErrors* gWHDG_Raa_0005       = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0005");

    TFile* fileChargedPiRAA                 = new TFile(fileNameChargedPiRAA);
    TH1F* histoChargedPiRAAStat0005         = (TH1F*)fileChargedPiRAA->Get("RAAPion_Stat_0_5");
    TH1F* histoChargedPiRAASys0005          = (TH1F*)fileChargedPiRAA->Get("RAAPion_Syst_0_5");
    TGraphAsymmErrors* graphChargedPiRAAStat0005    = new TGraphAsymmErrors(histoChargedPiRAAStat0005);
    TGraphAsymmErrors* graphChargedPiRAASys0005     = new TGraphAsymmErrors(histoChargedPiRAASys0005);
    graphChargedPiRAAStat0005->RemovePoint(0);;
    graphChargedPiRAASys0005->RemovePoint(0);;
    TGraphAsymmErrors* graphChargedPiRAAStat0005WX  = (TGraphAsymmErrors*)graphChargedPiRAAStat0005->Clone("RAAPion_Stat_0_5_wXErr");
    ProduceGraphAsymmWithoutXErrors(graphChargedPiRAAStat0005WX);
    
    
    //**************************************************************************************************
    //******************************** Plotting RAA with other energies ********************************
    //**************************************************************************************************
    TCanvas* canvasRAA_Comp             = new TCanvas("canvasRAA_Comp","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRAA_Comp,  0.1, 0.01, 0.015, 0.1);

    Int_t textSizeLabelsPixelRAA        = 50;
    Double_t marginRAA                  = 0.14*1200;
    Double_t textsizeLabelsRAA          = 0;
    Double_t textsizeFacRAA             = 0;

    if (canvasRAA_Comp->XtoPixel(canvasRAA_Comp->GetX2()) < canvasRAA_Comp->YtoPixel(canvasRAA_Comp->GetY1())){
        textsizeLabelsRAA               = (Double_t)textSizeLabelsPixelRAA/canvasRAA_Comp->XtoPixel(canvasRAA_Comp->GetX2()) ;
        textsizeFacRAA                  = (Double_t)1./canvasRAA_Comp->XtoPixel(canvasRAA_Comp->GetX2()) ;
    } else {
        textsizeLabelsRAA               = (Double_t)textSizeLabelsPixelRAA/canvasRAA_Comp->YtoPixel(canvasRAA_Comp->GetY1());
        textsizeFacRAA                  = (Double_t)1./canvasRAA_Comp->YtoPixel(canvasRAA_Comp->GetY1());
    }

    TH2F * histo2DRAAAll3Up             = new TH2F("histo2DRAAAll3Up","histo2DRAAAll3Up",1000,0.,20.,1000,-0.05,10.    );
    SetStyleHistoTH2ForGraphs(histo2DRAAAll3Up, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 512, 505); //#frac{#frac{1
    histo2DRAAAll3Up->GetYaxis()->SetLabelOffset(0.005);
    histo2DRAAAll3Up->GetXaxis()->SetRangeUser(0.,14.2);
    histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.00,2.1);
    histo2DRAAAll3Up->DrawCopy("");

    DrawGammaSetMarkerTGraphAsym(graphPi0RAAPbPbSysComb0010, markerStyleComb0010PbPb,markerSizeComb0010PbPb, colorComb0010PbPb , colorComb0010PbPb, widthLinesBoxes, kTRUE, 0);
    graphPi0RAAPbPbSysComb0010->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphPi0RAAPbPbStatComb0010WX, markerStyleComb0010PbPb,markerSizeComb0010PbPb, colorComb0010PbPb , colorComb0010PbPb);
    graphPi0RAAPbPbStatComb0010WX->Draw("p,same,e1,z");
    DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAA_0010, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
    DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAA_0010, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
    DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAA_0010, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);
    DrawGammaSetMarkerTGraphErr(graphWA98_17_3GeVRAA_0013, markerStyleWA98,markerSizeWA98, kGray+2 , kGray+2);

    graphPHENIX200GeVRAA_0010->Draw("p,same,e1");   
    graphPHENIX39GeVRAA_0010->Draw("p,same,e1");    
    graphPHENIX62GeVRAA_0010->Draw("p,same,e1");    
    graphWA98_17_3GeVRAA_0013->Draw("p,same,e1");   


    histo2DRAAAll3Up->Draw("axis,same");

    TLatex *labelRAAALICEPbPb0010   = new TLatex(0.35,0.93,"#pi^{0} ALICE    0-10% Pb-Pb");
    SetStyleTLatex( labelRAAALICEPbPb0010, 0.85*textsizeLabelsRAA,4);
    labelRAAALICEPbPb0010->Draw();
    TLegend* legendRAASinglePbPb0010 = GetAndSetLegend2(0.35, 0.88, 0.35+0.3, 0.92,0.85*textsizeLabelsRAA, 1, "", 42);new TLegend(0.35,0.88,0.65,0.92);
    legendRAASinglePbPb0010->AddEntry(graphPi0RAAPbPbSysComb0010,"#sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
    legendRAASinglePbPb0010->Draw();

    TLatex *labelRAAPHENIXPbPb0010  = new TLatex(0.35,0.83,"#pi^{0} PHENIX 0-10% Au-Au");
    SetStyleTLatex( labelRAAPHENIXPbPb0010, 0.85*textsizeLabelsRAA,4);
    labelRAAPHENIXPbPb0010->Draw();

    TLegend* legendRAARHICPbPb0010  = GetAndSetLegend2(0.35, 0.73, 0.35+0.6, 0.82,0.85*textsizeLabelsRAA, 2, "", 42);
    legendRAARHICPbPb0010->AddEntry(graphPHENIX200GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
    legendRAARHICPbPb0010->AddEntry(graphPHENIX62GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
    legendRAARHICPbPb0010->AddEntry(graphPHENIX39GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
    legendRAARHICPbPb0010->Draw();

    TLatex *labelRAAWA98PbPb0010    = new TLatex(0.35,0.68,"#pi^{0} WA98     0-13% Pb-Pb");
    SetStyleTLatex( labelRAAWA98PbPb0010, 0.85*textsizeLabelsRAA,4);
    labelRAAWA98PbPb0010->Draw();

    TLegend* legendRAASPSPbPb0010   = GetAndSetLegend2(0.35, 0.67-0.04, 0.35+0.6, 0.67,0.85*textsizeLabelsRAA, 2, "", 42);
    legendRAASPSPbPb0010->AddEntry(graphWA98_17_3GeVRAA_0013,"#sqrt{#it{s}_{_{NN}}} = 17.3 GeV","p");
    legendRAASPSPbPb0010->Draw();

    TBox* boxErrorNorm0010_Single   = CreateBoxConv(colorComb0010PbPbBox, 0.2, 1.-normErr0010 , 0.5, 1.+normErr0010);
    boxErrorNorm0010_Single->Draw();

    DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);


    canvasRAA_Comp->Update();
    canvasRAA_Comp->Print(Form("%s/RAA_0010_OtherSystems.%s",outputDir.Data(),suffix.Data()));

    //**************************************************************************************************
    //******************************** Plotting RAA 0-5% with control probes ***************************
    //**************************************************************************************************
    canvasRAA_Comp->SetLogy(1);
    
    TH1F * histo2DRAAAllUp             = new TH1F("histo2DRAAAllUp","histo2DRAAAllUp",1000,0.,21.5);
    SetStyleHistoTH1ForGraphs(histo2DRAAAllUp, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 512, 505); //#frac{#frac{1
    histo2DRAAAllUp->GetYaxis()->SetLabelOffset(0.005);
    histo2DRAAAllUp->GetYaxis()->SetRangeUser(0.03,35);
    histo2DRAAAllUp->DrawCopy("");

    DrawGammaSetMarkerTGraphAsym(graphPi0RAAPbPbSysComb0005, markerStyleComb0005PbPb,markerSizeComb0005PbPb, colorComb0005PbPb , colorComb0005PbPb, widthLinesBoxes, kTRUE, 0);
    graphPi0RAAPbPbSysComb0005->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphPi0RAAPbPbStatComb0005WX, markerStyleComb0005PbPb,markerSizeComb0005PbPb, colorComb0005PbPb , colorComb0005PbPb);
    graphPi0RAAPbPbStatComb0005WX->Draw("p,same,e1,z");

    DrawGammaSetMarkerTGraphAsym(graphPi0RpApPbSysComb, markerStyleCombpPb,markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE, 0);
    graphPi0RpApPbSysComb->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphPi0RpApPbStatCombWX, markerStyleCombpPb,markerSizeCombpPb, colorCombpPb , colorCombpPb);
    graphPi0RpApPbStatCombWX->Draw("p,same,e1,z");
    
    DrawGammaSetMarkerTGraphAsym(graphGammaRAAPbPbSysComb0020, markerStyleCombGamma0020PbPb,markerSizeCombGamma0020PbPb, colorComb0020PbPb , colorComb0020PbPb, widthLinesBoxes, kTRUE, 0);
    graphGammaRAAPbPbSysComb0020->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphGammaRAAPbPbStatComb0020WX, markerStyleCombGamma0020PbPb,markerSizeCombGamma0020PbPb, colorComb0020PbPb , colorComb0020PbPb);
    graphGammaRAAPbPbStatComb0020WX->Draw("p,same,e1,z");

    TLatex *labelRAAALICEPbPb       = new TLatex(0.52,0.92,"Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV ");
    SetStyleTLatex( labelRAAALICEPbPb, 0.85*textsizeLabelsRAA,4);
    labelRAAALICEPbPb->Draw();
    TLegend* legendRAASinglePbPb    = GetAndSetLegend2(0.52, 0.91-4*0.04, 0.52+0.3, 0.91,0.75*textsizeLabelsRAA, 1, "", 42);
    legendRAASinglePbPb->AddEntry(graphPi0RAAPbPbSysComb0005,"0-5% #pi^{0} ALICE","pf");
    legendRAASinglePbPb->AddEntry((TObject*)0,"EPJC 74 (2014) 3108","");
    legendRAASinglePbPb->AddEntry(graphGammaRAAPbPbSysComb0020,"0-20% #gamma_{dir} ALICE","pf");
    legendRAASinglePbPb->AddEntry((TObject*)0,"PLB 754 (2016) 235-248","");
    legendRAASinglePbPb->Draw();

    TLatex *labelRAAALICEpPb        = new TLatex(0.52,0.71,"p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.023 TeV");
    SetStyleTLatex( labelRAAALICEpPb, 0.85*textsizeLabelsRAA,4);
    labelRAAALICEpPb->Draw();
    TLegend* legendRAASinglePpPb    = GetAndSetLegend2(0.52, 0.70-0.04, 0.52+0.3, 0.70,0.75*textsizeLabelsRAA, 1, "", 42);
    legendRAASinglePpPb->AddEntry(graphPi0RpApPbSysComb,"0-100% #pi^{0} ALICE","pf");
    legendRAASinglePpPb->Draw();

    boxErrorNorm0010_Single->Draw();
    
    TBox* boxErrorNorm0020_Single   = CreateBoxConv(colorComb0020PbPbBox, 0.55, 1.-normErr0020 , 0.85, 1.+normErr0020);
    boxErrorNorm0020_Single->Draw();
    TBox* boxErrorNormpPb_Single    = CreateBoxConv(colorCombpPbBox, 20.4, 1.-normErrpPb , 20.7, 1.+normErrpPb);
    boxErrorNormpPb_Single->Draw();

    DrawGammaLines(0., 21.5 , 1, 1 ,1,kGray+1);


    TLatex *labelALICEPrelim        = new TLatex(0.13,0.92,"ALICE preliminary");
    SetStyleTLatex( labelALICEPrelim, 0.85*textsizeLabelsRAA,4);
    labelALICEPrelim->Draw();
    histo2DRAAAllUp->Draw("axis,same");
    
    canvasRAA_Comp->Update();
    canvasRAA_Comp->Print(Form("%s/RAA0005_WithControlProbes.%s",outputDir.Data(),suffix.Data()));

    //**************************************************************************************************
    //************************* Plotting RAA 0-5% with control probes with theory **********************
    //**************************************************************************************************
    canvasRAA_Comp->SetLogy(1);
    histo2DRAAAllUp->DrawCopy("");
    
    graphPi0RAAPbPbSysComb0005->Draw("E2same");
    graphPi0RAAPbPbStatComb0005WX->Draw("p,same,e1,z");
    graphPi0RpApPbSysComb->Draw("E2same");
    graphPi0RpApPbStatCombWX->Draw("p,same,e1,z");
    graphGammaRAAPbPbSysComb0020->Draw("E2same");
    graphGammaRAAPbPbStatComb0020WX->Draw("p,same,e1,z");

    DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_0005, 1,2, colorVitevBas0005 ,colorVitevBas0005,0.8);
    Vitev_Bas_Raa_0005->SetFillStyle(fillStyleVitev);
    Vitev_Bas_Raa_0005->SetFillColor(colorVitevBas0005);
    Vitev_Bas_Raa_0005->Draw("3 same");
    DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_0005, 1,2, colorWHDG0005, colorWHDG0005,0.8);
    gWHDG_Raa_0005->SetFillStyle(fillStyleWHDG);
    gWHDG_Raa_0005->SetFillColor(colorWHDG0005);
    gWHDG_Raa_0005->Draw("3 same");    
    

    labelRAAALICEPbPb->Draw();
    legendRAASinglePbPb->Draw();

    labelRAAALICEpPb->Draw();
    legendRAASinglePpPb->Draw();

    boxErrorNorm0010_Single->Draw();    
    boxErrorNorm0020_Single->Draw();
    boxErrorNormpPb_Single->Draw();

    DrawGammaLines(0., 21.5 , 1, 1 ,1,kGray+1);

    TLegend* legendRAAPi0TheoryPbPb     = GetAndSetLegend2(0.62, 0.13, 0.62+0.3, 0.13+0.04,0.75*textsizeLabelsRAA, 2, "", 42, 0.45);
    legendRAAPi0TheoryPbPb->AddEntry(Vitev_Bas_Raa_0005,"GLV","f");
    legendRAAPi0TheoryPbPb->AddEntry(gWHDG_Raa_0005,"WHDG","f");
    legendRAAPi0TheoryPbPb->Draw();    

    labelALICEPrelim->Draw();
    histo2DRAAAllUp->Draw("axis,same");
    
    canvasRAA_Comp->Update();
    canvasRAA_Comp->Print(Form("%s/RAA0005_WithControlProbes_Theory.%s",outputDir.Data(),suffix.Data()));

    //**************************************************************************************************
    //************************** Plotting RAA 0-5% with control probes, pi+-, theory *******************
    //**************************************************************************************************
    canvasRAA_Comp->SetLogy(1);
    
    histo2DRAAAllUp->DrawCopy("");

    graphPi0RAAPbPbSysComb0005->Draw("E2same");
    graphPi0RAAPbPbStatComb0005WX->Draw("p,same,e1,z");

    DrawGammaSetMarkerTGraphAsym(graphChargedPiRAASys0005, markerStyleChPi0005PbPb,markerSizeComb0005PbPb, colorChPi0005PbPb , colorChPi0005PbPb, widthLinesBoxes, kTRUE, 0);
    graphChargedPiRAASys0005->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphChargedPiRAAStat0005WX, markerStyleChPi0005PbPb,markerSizeComb0005PbPb, colorChPi0005PbPb , colorChPi0005PbPb);
    graphChargedPiRAAStat0005WX->Draw("p,same,e1,z");
    
    graphPi0RpApPbSysComb->Draw("E2same");
    graphPi0RpApPbStatCombWX->Draw("p,same,e1,z");
    
    graphGammaRAAPbPbSysComb0020->Draw("E2same");
    graphGammaRAAPbPbStatComb0020WX->Draw("p,same,e1,z");

    Vitev_Bas_Raa_0005->Draw("3 same");
    gWHDG_Raa_0005->Draw("3 same");    
    
    boxErrorNorm0010_Single->Draw();
    boxErrorNorm0020_Single->Draw();
    boxErrorNormpPb_Single->Draw();
    
    TLatex *labelRAAALICEPbPb2      = new TLatex(0.52,0.92,"Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV ");
    SetStyleTLatex( labelRAAALICEPbPb2, 0.85*textsizeLabelsRAA,4);
    labelRAAALICEPbPb2->Draw();
    TLegend* legendRAASinglePbPb2   = GetAndSetLegend2(0.52, 0.91-6*0.04, 0.52+0.3, 0.91,0.75*textsizeLabelsRAA, 1, "", 42);
    legendRAASinglePbPb2->AddEntry(graphPi0RAAPbPbSysComb0005,"0-5% #pi^{0} ALICE","pf");
    legendRAASinglePbPb2->AddEntry((TObject*)0,"EPJC 74 (2014) 3108","");
    legendRAASinglePbPb2->AddEntry(graphChargedPiRAASys0005,"0-5% #pi^{#pm} ALICE","pf");
    legendRAASinglePbPb2->AddEntry((TObject*)0,"PLB 736 (2014) 196-207","");
    legendRAASinglePbPb2->AddEntry(graphGammaRAAPbPbSysComb0020,"0-20% #gamma_{dir} ALICE","pf");
    legendRAASinglePbPb2->AddEntry((TObject*)0,"PLB 754 (2016) 235-248","");
    legendRAASinglePbPb2->Draw();

    TLatex *labelRAAALICEpPb2       = new TLatex(0.13,0.92,"p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.023 TeV");
    SetStyleTLatex( labelRAAALICEpPb2, 0.85*textsizeLabelsRAA,4);
    labelRAAALICEpPb2->Draw();
    TLegend* legendRAASinglePpPb2   = GetAndSetLegend2(0.13, 0.91-1*0.04, 0.13+0.3, 0.91,0.75*textsizeLabelsRAA, 1, "", 42);
    legendRAASinglePpPb2->AddEntry(graphPi0RpApPbSysComb,"0-100% #pi^{0} ALICE","pf");
    legendRAASinglePpPb2->Draw();


    TLatex *labelALICEPrelim2 = new TLatex(0.13,0.135,"ALICE preliminary");
    SetStyleTLatex( labelALICEPrelim2, 0.85*textsizeLabelsRAA,4);
    labelALICEPrelim2->Draw();
    
    boxErrorNorm0010_Single->Draw();    
    boxErrorNorm0020_Single->Draw();
    boxErrorNormpPb_Single->Draw();

    DrawGammaLines(0., 21.5 , 1, 1 ,1,kGray+1);

    legendRAAPi0TheoryPbPb->Draw();    
    histo2DRAAAllUp->Draw("axis,same");
    
    canvasRAA_Comp->Update();
    canvasRAA_Comp->Print(Form("%s/RAA0005_WithControlProbes_Theory_ChargedPi.%s",outputDir.Data(),suffix.Data()));
    
    //**************************************************************************************************
    //******************************** Plotting RAA 0-5% with control probes ***************************
    //**************************************************************************************************
    canvasRAA_Comp->SetLogy(0);

    TH1F * histo1DRAAAllUp             = new TH1F("histo1DRAAAllUp","histo1DRAAAllUp",1000,0.,21.5);
    SetStyleHistoTH1ForGraphs(histo1DRAAAllUp, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 512, 505); //#frac{#frac{1
    histo1DRAAAllUp->GetYaxis()->SetLabelOffset(0.005);    
    histo1DRAAAllUp->GetYaxis()->SetRangeUser(0.0,2.49);
    histo1DRAAAllUp->DrawCopy("");

    graphPi0RAAPbPbSysComb0005->Draw("E2same");
    graphPi0RAAPbPbStatComb0005WX->Draw("p,same,e1,z");

    graphPi0RpApPbSysComb->Draw("E2same");
    graphPi0RpApPbStatCombWX->Draw("p,same,e1,z");

    while(graphGammaRAAPbPbSysComb0020->GetX()[0] < 5) graphGammaRAAPbPbSysComb0020->RemovePoint(0);
    while(graphGammaRAAPbPbStatComb0020WX->GetX()[0] < 5) graphGammaRAAPbPbStatComb0020WX->RemovePoint(0);    
    graphGammaRAAPbPbSysComb0020->Draw("E2same");
    graphGammaRAAPbPbStatComb0020WX->Draw("p,same,e1,z");


    labelRAAALICEPbPb->Draw();
    legendRAASinglePbPb->Draw();

    labelRAAALICEpPb->Draw();
    legendRAASinglePpPb->Draw();

    boxErrorNorm0010_Single->Draw();    
    boxErrorNorm0020_Single->Draw();
    boxErrorNormpPb_Single->Draw();

    DrawGammaLines(0., 21.5 , 1, 1 ,1,kGray+1);
    labelALICEPrelim->Draw();

    histo1DRAAAllUp->Draw("axis,same");
    
    canvasRAA_Comp->Update();
    canvasRAA_Comp->Print(Form("%s/RAA0005_WithControlProbes_GammaCut.%s",outputDir.Data(),suffix.Data()));

    //**************************************************************************************************
    //******************************** Plotting RAA 0-5% with control probes gamma cut *****************
    //**************************************************************************************************
    canvasRAA_Comp->SetLogy(0);
    
    histo1DRAAAllUp->DrawCopy("");

    graphPi0RAAPbPbSysComb0005->Draw("E2same");
    graphPi0RAAPbPbStatComb0005WX->Draw("p,same,e1,z");

    graphPi0RpApPbSysComb->Draw("E2same");
    graphPi0RpApPbStatCombWX->Draw("p,same,e1,z");
    
    graphGammaRAAPbPbSysComb0020->Draw("E2same");
    graphGammaRAAPbPbStatComb0020WX->Draw("p,same,e1,z");

    Vitev_Bas_Raa_0005->Draw("3 same");
    gWHDG_Raa_0005->Draw("3 same");    
    
    labelRAAALICEPbPb->Draw();
    legendRAASinglePbPb->Draw();

    TLatex *labelRAAALICEpPb3       = new TLatex(0.13,0.87,"p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.023 TeV");
    SetStyleTLatex( labelRAAALICEpPb3, 0.85*textsizeLabelsRAA,4);
    labelRAAALICEpPb3->Draw();
    TLegend* legendRAASinglePpPb3   = GetAndSetLegend2(0.13, 0.86-1*0.04, 0.13+0.3, 0.86,0.75*textsizeLabelsRAA, 1, "", 42);
    legendRAASinglePpPb3->AddEntry(graphPi0RpApPbSysComb,"0-100% #pi^{0} ALICE","pf");
    legendRAASinglePpPb3->Draw();

    boxErrorNorm0010_Single->Draw();    
    boxErrorNorm0020_Single->Draw();
    boxErrorNormpPb_Single->Draw();

    DrawGammaLines(0., 21.5 , 1, 1 ,1,kGray+1);

    TLegend* legendRAAPi0TheoryPbPb2 = GetAndSetLegend2(0.52, 0.74-1*0.04, 0.52+0.32, 0.74,0.75*textsizeLabelsRAA, 2, "", 42, 0.45 );
    legendRAAPi0TheoryPbPb2->AddEntry(Vitev_Bas_Raa_0005,"GLV","f");
    legendRAAPi0TheoryPbPb2->AddEntry(gWHDG_Raa_0005,"WHDG","f");
    legendRAAPi0TheoryPbPb2->Draw();    

    labelALICEPrelim->Draw();
    histo1DRAAAllUp->Draw("axis,same");
    
    canvasRAA_Comp->Update();
    canvasRAA_Comp->Print(Form("%s/RAA0005_WithControlProbes_Theory_GammaCut.%s",outputDir.Data(),suffix.Data()));

    //**************************************************************************************************
    //************************** Plotting RAA 0-5% with control probes, pi+-, theory *******************
    //**************************************************************************************************
    canvasRAA_Comp->SetLogy(0);
    
    histo1DRAAAllUp->DrawCopy("");

    graphPi0RAAPbPbSysComb0005->Draw("E2same");
    graphPi0RAAPbPbStatComb0005WX->Draw("p,same,e1,z");

    graphChargedPiRAASys0005->Draw("E2same");
    graphChargedPiRAAStat0005WX->Draw("p,same,e1,z");
    
    graphPi0RpApPbSysComb->Draw("E2same");
    graphPi0RpApPbStatCombWX->Draw("p,same,e1,z");
    
    graphGammaRAAPbPbSysComb0020->Draw("E2same");
    graphGammaRAAPbPbStatComb0020WX->Draw("p,same,e1,z");

    Vitev_Bas_Raa_0005->Draw("3 same");
    gWHDG_Raa_0005->Draw("3 same");    
    
    boxErrorNorm0010_Single->Draw();
    boxErrorNorm0020_Single->Draw();
    boxErrorNormpPb_Single->Draw();
    
    labelRAAALICEPbPb2->Draw();
    legendRAASinglePbPb2->Draw();

    labelRAAALICEpPb3->Draw();
    legendRAASinglePpPb3->Draw();

    DrawGammaLines(0., 21.5 , 1, 1 ,1,kGray+1);

    labelALICEPrelim->Draw();
    
    TLegend* legendRAAPi0TheoryPbPb3    = GetAndSetLegend2(0.52, 0.66-1*0.04, 0.52+0.32, 0.66,0.75*textsizeLabelsRAA, 2, "", 42, 0.45 );
    legendRAAPi0TheoryPbPb3->AddEntry(Vitev_Bas_Raa_0005,"GLV","f");
    legendRAAPi0TheoryPbPb3->AddEntry(gWHDG_Raa_0005,"WHDG","f");
    legendRAAPi0TheoryPbPb3->Draw();    

    histo1DRAAAllUp->Draw("axis,same");
    
    canvasRAA_Comp->Update();
    canvasRAA_Comp->Print(Form("%s/RAA0005_WithControlProbes_Theory_ChargedPi_GammaCut.%s",outputDir.Data(),suffix.Data()));
    
    return;    
}
    
