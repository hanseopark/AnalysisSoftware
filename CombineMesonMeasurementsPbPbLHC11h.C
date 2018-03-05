/****************************************************************************************************************************
********    provided by Gamma Conversion Group, PWGGA,                                                                  *****
********    Lucia Leardini, leardini@cern.ch                                                                            *****
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
#include "CombineMesonMeasurementsPbPbLHC11h.h"


void CombineMesonMeasurementsPbPbLHC11h(TString meson = "Pi0",
                                        Bool_t noXerrorBars = kFALSE,
                                        TString thisthesis="", //"This thesis",//"ALICE work in progress", //"This thesis",
                                        Bool_t PaperPi0 = kTRUE,
//                                         Bool_t thesisPlotting = kFALSE, //kTRUE,
                                        TString suffix = "pdf",
                                        TString bWCorrection="X",
                                        Bool_t quietMode = kTRUE){

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetEndErrorSize(0);
    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput           = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;

    //___________________________________ Labels definition _____________________________________________
//     TLatex *labelFactorPi0104 = new TLatex(0.85,0.405,"#times 4#upoint10^{2}");
//     SetStyleTLatex( labelFactorPi0104,FontSize,4,colorCombo0010);
//     TLatex *labelFactorPi0100 = new TLatex(0.85,0.34,"#times 10^{2}");
//     SetStyleTLatex( labelFactorPi0100,FontSize,4,colorCombo2050);
//     TLatex *labelFactorEta4 = new TLatex(0.85,0.247,"#times 4");
//     SetStyleTLatex( labelFactorEta4,FontSize,4,colorCombo0010);
//     TLatex *labelFactorEta = new TLatex(0.85,0.2,"#times 10^{-1}");
//     SetStyleTLatex( labelFactorEta,FontSize,4,colorCombo2050);
    TLatex *labelFactorPi0104 = new TLatex(0.83,0.41,"#times 4#upoint10^{2}");
    SetStyleTLatex( labelFactorPi0104,FontSize,4,colorCombo0010);
    TLatex *labelFactorPi0100 = new TLatex(0.83,0.35,"#times 10^{2}");
    SetStyleTLatex( labelFactorPi0100,FontSize,4,colorCombo2050);
    TLatex *labelFactorEta = new TLatex(0.83,0.16,"#times 2#upoint10^{-1}");
    SetStyleTLatex( labelFactorEta,FontSize,4,colorCombo2050);

    Bool_t thesisPlotting;
    if(thisthesis.CompareTo("")==0) thesisPlotting = kFALSE;
    else thesisPlotting = kTRUE;

    Double_t maxX, minX;
    if(thesisPlotting || PaperPi0){
        maxX = 35;
        minX = 0.5;
    } else {
        maxX = 35;
        minX = 0.25;
    }

    TLatex *thesisLabelHigh = new TLatex(0.8,0.9,thisthesis.Data());
    SetStyleTLatex( thesisLabelHigh,FontSize,4);
    TLatex *thesisLabelHighRight = new TLatex(0.73,0.927,thisthesis.Data());
    SetStyleTLatex( thesisLabelHighRight,FontSize,4);
    TLatex *thesisLabelHighLeft = new TLatex(0.17,0.88,thisthesis.Data());
    SetStyleTLatex( thesisLabelHighLeft,FontSize,4);
    TLatex *thesisLabelHighLeft2 = new TLatex(0.132,0.91,thisthesis.Data());
    SetStyleTLatex( thesisLabelHighLeft2,FontSize,4);
    TLatex *thesisLabelLowRight = new TLatex(0.8,0.12,thisthesis.Data());
    SetStyleTLatex( thesisLabelLowRight,FontSize,4);
    TLatex *thesisLabelLowLeft = new TLatex(0.205,0.11,thisthesis.Data());
    SetStyleTLatex( thesisLabelLowLeft,FontSize,4);

    if(meson.CompareTo("Pi0")==0){     labelFactorLower = new TLatex(0.81,0.4615,"#times 4");}
    else if(meson.CompareTo("Eta")==0){labelFactorLower = new TLatex(0.81,0.382,"#times 4");}
    SetStyleTLatex( labelFactorLower, FontSize,4,colorCombo0010);

    if(meson.CompareTo("Pi0")==0){     labelFactorLowerOnlyPbPb = new TLatex(0.81,0.268,"#times 4");}
    else if(meson.CompareTo("Eta")==0){labelFactorLowerOnlyPbPb = new TLatex(0.81,0.32,"#times 4");}
    SetStyleTLatex( labelFactorLowerOnlyPbPb,FontSize,4,colorCombo0010);

    if(meson.CompareTo("Pi0")==0){     labelFactorUpper = new TLatex(0.415,0.89,"#times 4");}
    else if(meson.CompareTo("Eta")==0){labelFactorUpper = new TLatex(0.43,0.835,"#times 4");}
    SetStyleTLatex( labelFactorUpper,FontSize,4,colorCombo0010);

    if(meson.CompareTo("Pi0")==0){     labelSystOnlyPbPb= new TLatex(0.75,0.93,"#pi^{0} #rightarrow #gamma#gamma");}
    else if(meson.CompareTo("Eta")==0){labelSystOnlyPbPb= new TLatex(0.75,0.93,"#eta #rightarrow #gamma#gamma");}
    SetStyleTLatex(labelSystOnlyPbPb,FontSize,4);

    if(meson.CompareTo("Pi0")==0){     labelSyst= new TLatex(0.6,0.93,"#pi^{0} #rightarrow #gamma#gamma");}
    else if(meson.CompareTo("Eta")==0){labelSyst= new TLatex(0.6,0.93,"#eta #rightarrow #gamma#gamma");}
    SetStyleTLatex(labelSyst,FontSize,4);

    if(meson.CompareTo("Pi0")==0){     labelSystRaa= new TLatex(0.6,0.89,"#pi^{0} #rightarrow #gamma#gamma");}
    else if(meson.CompareTo("Eta")==0){labelSystRaa= new TLatex(0.6,0.89,"#eta #rightarrow #gamma#gamma ");}
    SetStyleTLatex( labelSystRaa,FontSize,4);

    TLatex *labelEnergyIndMeasRAA0010 = new TLatex(0.15,0.88,collisionSystemPbPb0010.Data());
    SetStyleTLatex( labelEnergyIndMeasRAA0010,FontSize,4);
    TLatex *labelEnergyIndMeasRAA2050 = new TLatex(0.15,0.88,collisionSystemPbPb2050.Data());
    SetStyleTLatex( labelEnergyIndMeasRAA2050,FontSize,4);
    TLatex *labelDetSysIndMeasRAA = new TLatex(0.15,0.84,"#pi^{0} #rightarrow #gamma#gamma");
    if(meson.CompareTo("Eta")==0) labelDetSysIndMeasRAA= new TLatex(0.15,0.84,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysIndMeasRAA, FontSize,4);

    TLatex *labelEnergyInvYieldSectionPi0LHC11h = new TLatex(0.6,0.88,collisionSystem2760GeV.Data());
    SetStyleTLatex( labelEnergyInvYieldSectionPi0LHC11h,FontSize,4);
    if(meson.CompareTo("Pi0")==0){
        labelDetSysInvYieldSectionPi0LHC11h= new TLatex(0.6,0.84,"#pi^{0} #rightarrow #gamma#gamma");
    } else if(meson.CompareTo("Eta")==0){
        labelDetSysInvYieldSectionPi0LHC11h= new TLatex(0.6,0.84,"#eta #rightarrow #gamma#gamma");
    }
    SetStyleTLatex( labelDetSysInvYieldSectionPi0LHC11h, FontSize,4);

    TLatex *labelEnergyInvYieldSectionPi0LHC11hnoPrelim = new TLatex(0.6,0.92,collisionSystem2760GeV.Data());
    SetStyleTLatex( labelEnergyInvYieldSectionPi0LHC11hnoPrelim, 0.035,4);
    if(meson.CompareTo("Pi0")==0){
        labelDetSysInvYieldSectionPi0LHC11hnoPrelim= new TLatex(0.6,0.88,"#pi^{0} #rightarrow #gamma#gamma");
    } else if(meson.CompareTo("Eta")==0){
        labelDetSysInvYieldSectionPi0LHC11hnoPrelim= new TLatex(0.6,0.88,"#eta #rightarrow #gamma#gamma");
    }
    SetStyleTLatex( labelDetSysInvYieldSectionPi0LHC11hnoPrelim, 0.035,4);

    if(meson.CompareTo("Pi0")==0){
        labelDetSysInvYieldSectionPi0LHC11hwithPP= new TLatex(0.62,0.88,"#pi^{0} #rightarrow #gamma#gamma");
    } else if(meson.CompareTo("Eta")==0){
        labelDetSysInvYieldSectionPi0LHC11hwithPP= new TLatex(0.62,0.88,"#eta #rightarrow #gamma#gamma");
    }
    SetStyleTLatex( labelDetSysInvYieldSectionPi0LHC11hwithPP, 0.035,4);

    TLatex *labelEtaToPi0Energy = new TLatex(0.12,0.89,collisionSystem2760GeV.Data());
    SetStyleTLatex( labelEtaToPi0Energy,0.04,4);

    for (Int_t i = 0; i < 11; i++){
        colorDet[i]                 = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
        colorDetMC[i]               = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
        markerStyleDet[i]           = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerStyleDetMC[i]         = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
        markerSizeDet[i]            = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
        markerSizeDetMC[i]          = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
    }

    Double_t mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    Double_t mass = 0.;
    if (meson.CompareTo("Pi0")==0){
        mass = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    } else if (meson.CompareTo("Eta")==0){
        mass = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    }
    Double_t normErr0010 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0010/tAA0010),2)+pow((commonCentralityErr0010/100),2),0.5);
    Double_t normErr2040 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2040/tAA2040),2)+pow((commonCentralityErr2040/100),2),0.5);
    Double_t normErr2050 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2050/tAA2050),2)+pow((commonCentralityErr2040/100),2),0.5);

    TBox* boxErrorNorm0010          = CreateBoxConv(kRed-7, 0.25, 1.-normErr0010 , 0.5, 1.+normErr0010);
    TBox* boxErrorNorm2050          = CreateBoxConv(kAzure-4, 0.5, 1.-normErr2050 , 0.75, 1.+normErr2050);
    TBox* boxErrorNorm0010_Single   = CreateBoxConv(colorComb0005Box, 0.2, 1.-normErr0010 , 0.5, 1.+normErr0010);
    TBox* boxErrorNorm2040_Single   = CreateBoxConv(colorComb2040Box, 0.2, 1.-normErr2040 , 0.5, 1.+normErr2040);
    TBox* boxErrorNorm2050_Single   = CreateBoxConv(colorCombo2050, 0.2, 1.-normErr2050 , 0.5, 1.+normErr2050);
    TBox* boxErrorNorm0010Only      = CreateBoxConv(kRed-7, 0.25, 1.-normErr0010 , 0.5, 1.+normErr0010);
    TBox* boxErrorNorm2050Only      = CreateBoxConv(kAzure-4, 0.25, 1.-normErr2050 , 0.5, 1.+normErr2050);
    TBox* boxErrorNorm2050EPOnly    = CreateBoxConv(kBlue-7, 0.25, 1.-normErr2050 , 0.5, 1.+normErr2050);

    //___________________________________ Declaration of files _____________________________________________

    TString fileNameTheory                  = "ExternalInputPbPb/Theory/TheoryCompilationPbPbforLHC11h.root";
    TString fileNameALICEData               = Form("%s/%s/CombineMesonMeasurementsPbPb2760GeVX/InputALICEResultsPbPb2760GeV_%s.root",suffix.Data(),dateForOutput.Data(),dateForOutput.Data());
    TString fileNameDataOtherEnergyInput    = "ExternalInputPbPb/OtherExperiments/DataCompilationFromOtherEnergiesPbPbPi0andEta.root";

    TString outputDir                       = Form("%s/%s/CombineMesonMeasurementsPbPb2760GeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
    TString paperPlots                      = Form("%s/PaperPlots_%s",outputDir.Data(),dateForOutput.Data());
    TString PubNotePlots                    = Form("%s/PubNotePlots_%s",outputDir.Data(),dateForOutput.Data());
    TString rootFiles                       = Form("%s/rootFiles",outputDir.Data());
    TString nameFinalResDat                 = Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+paperPlots);
    gSystem->Exec("mkdir -p "+PubNotePlots);
    gSystem->Exec("mkdir -p "+rootFiles);
    gSystem->Exec(Form("cp %s %s/Theory.root",                fileNameTheory.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/InputALICEData.root",        fileNameALICEData.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/OtherExperimentsInput.root", fileNameDataOtherEnergyInput.Data(), rootFiles.Data()));

    //*********************************************************************************************************//
    //*************************************        Theory inputs      *****************************************//
    //*********************************************************************************************************//
    TFile* fileTheoryGraphs = new TFile(fileNameTheory);

    //**************************** pQCD calculation ****************************//
    //inv yield = inv cross sec/xSection2760GeVppINEL
    //****************************     Pi0 DSS14    ****************************//
    TGraphAsymmErrors *graphNLOCalcmuHalfDSS14InvSecPi02760GeV  = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuHalfDSS14InvSecPi02760GeV");
    TGraphAsymmErrors *graphNLOCalcmuTwoDSS14InvSecPi02760GeV   = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuTwoDSS14InvSecPi02760GeV");
    TGraphAsymmErrors *graphNLOCalcDSS14InvSecPi02760GeV        = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS14InvCrossSecPi02760GeV");

    TGraphAsymmErrors* graphNLOCalcDSS14InvYieldPi0Band = (TGraphAsymmErrors*)graphNLOCalcDSS14InvSecPi02760GeV->Clone("graphNLOCalcDSS14InvYieldPi0Band");
    Double_t *yPi0PointDSS14muHalf = graphNLOCalcmuHalfDSS14InvSecPi02760GeV->GetY();
    Double_t *yPi0PointDSS14muOne = graphNLOCalcDSS14InvSecPi02760GeV->GetY();
    Double_t *yPi0PointDSS14muTwo  = graphNLOCalcmuTwoDSS14InvSecPi02760GeV->GetY();
    for(Int_t i =0; i<graphNLOCalcDSS14InvYieldPi0Band->GetN(); i++){
        graphNLOCalcDSS14InvYieldPi0Band->SetPointEYhigh(i,yPi0PointDSS14muTwo[i]-yPi0PointDSS14muOne[i]);
        graphNLOCalcDSS14InvYieldPi0Band->SetPointEYlow(i,yPi0PointDSS14muOne[i]-yPi0PointDSS14muHalf[i]);
    }
    graphNLOCalcDSS14InvYieldPi0Band->RemovePoint(0);
    graphNLOCalcDSS14InvYieldPi0Band = ScaleGraph(graphNLOCalcDSS14InvYieldPi0Band,1./xSection2760GeVppINEL);

    //****************************     Eta DSS07    ****************************//
    TGraphAsymmErrors* graphNLOCalcmuHalfDSS07InvSecEta2760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuHalfDSS07InvSecEta2760GeV");
    TGraphAsymmErrors* graphNLOCalcmuTwoDSS07InvSecEta2760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuTwoDSS07InvSecEta2760GeV");
    TGraphAsymmErrors* graphNLOCalcDSS07InvSecEta2760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS07InvCrossSecEta2760GeV");

    TGraphAsymmErrors* graphNLOCalcDSS07InvYieldEtaBand = (TGraphAsymmErrors*)graphNLOCalcDSS07InvSecEta2760GeV->Clone("graphNLOCalcDSS07InvYieldEtaBand");
    Double_t *yEtaPointDSS07muHalf = graphNLOCalcmuHalfDSS07InvSecEta2760GeV->GetY();
    Double_t *yEtaPointDSS07muOne = graphNLOCalcDSS07InvSecEta2760GeV->GetY();
    Double_t *yEtaPointDSS07muTwo  = graphNLOCalcmuTwoDSS07InvSecEta2760GeV->GetY();
    for(Int_t i =0; i<graphNLOCalcDSS07InvYieldEtaBand->GetN(); i++){
        graphNLOCalcDSS07InvYieldEtaBand->SetPointEYhigh(i,yEtaPointDSS07muTwo[i]-yEtaPointDSS07muOne[i]);
        graphNLOCalcDSS07InvYieldEtaBand->SetPointEYlow(i,yEtaPointDSS07muOne[i]-yEtaPointDSS07muHalf[i]);
    }
    graphNLOCalcDSS07InvYieldEtaBand->RemovePoint(0);
    graphNLOCalcDSS07InvYieldEtaBand = ScaleGraph(graphNLOCalcDSS07InvYieldEtaBand,1./xSection2760GeVppINEL);


    //**************************************** EPOS **********************************************
    //New EPOS (2?) (file from Anders email of the 19th Feb. 2016)
    TGraphErrors *graphEPOSPi0_0010 = (TGraphErrors*)fileTheoryGraphs->Get("graphEPOS_Pi0_0010");
    TGraphErrors *graphEPOSEta_0010 = (TGraphErrors*)fileTheoryGraphs->Get("graphEPOS_Eta_0010");
    TGraphErrors *graphEPOSPi0_2050 = (TGraphErrors*)fileTheoryGraphs->Get("graphEPOS_Pi0_2050");
    TGraphErrors *graphEPOSEta_2050 = (TGraphErrors*)fileTheoryGraphs->Get("graphEPOS_Eta_2050");
    TGraphErrors *graphEPOSEtaToPi0_0010 = (TGraphErrors*)fileTheoryGraphs->Get("graphEPOS_EtaToPi0_0010");
    TGraphErrors *graphEPOSEtaToPi0_2050 = (TGraphErrors*)fileTheoryGraphs->Get("graphEPOS_EtaToPi0_2050");

    DrawGammaSetMarkerTGraphErr(graphEPOSPi0_0010, 0, 0, colorEPOS0010,colorEPOS0010, 2, kFALSE, colorEPOS0010);
    graphEPOSPi0_0010->SetLineStyle(5);
    DrawGammaSetMarkerTGraphErr(graphEPOSPi0_2050, 0, 0, colorEPOS2050,colorEPOS2050, 2, kFALSE, colorEPOS2050);
    graphEPOSPi0_2050->SetLineStyle(5);
    DrawGammaSetMarkerTGraphErr(graphEPOSEta_0010, 0, 0, colorEPOS0010,colorEPOS0010, 2, kFALSE, colorEPOS0010);
    graphEPOSEta_0010->SetLineStyle(5);
    DrawGammaSetMarkerTGraphErr(graphEPOSEta_2050, 0, 0, colorEPOS2050,colorEPOS2050, 2, kFALSE, colorEPOS2050);
    graphEPOSEta_2050->SetLineStyle(5);
    if(meson.CompareTo("Pi0")==0){
        graphEPOSpred_0010 = (TGraphErrors*)graphEPOSPi0_0010->Clone("graphEPOSpred_0010");
        graphEPOSpred_2050 = (TGraphErrors*)graphEPOSPi0_2050->Clone("graphEPOSpred_2050");
    } else if(meson.CompareTo("Eta")==0){
        graphEPOSpred_0010 = (TGraphErrors*)graphEPOSEta_0010->Clone("graphEPOSpred_0010");
        graphEPOSpred_2050 = (TGraphErrors*)graphEPOSEta_2050->Clone("graphEPOSpred_2050");
    }
    DrawGammaSetMarkerTGraphErr(graphEPOSpred_0010, 0, 0, colorEPOS0010,colorEPOS0010, 2, kTRUE, colorEPOS0010);
    ProduceGraphErrWithoutXErrors(graphEPOSpred_0010);
    graphEPOSpred_0010->SetLineStyle(5);
    TGraphErrors *graphEPOSpredScaledforPlot_0010 = (TGraphErrors*)graphEPOSpred_0010->Clone("graphEPOSpredScaledforPlot_0010");
    graphEPOSpredScaledforPlot_0010 = ScaleGraph(graphEPOSpredScaledforPlot_0010,4);
    DrawGammaSetMarkerTGraphErr(graphEPOSpredScaledforPlot_0010, 0, 0, colorEPOS0010,colorEPOS0010, 2, kTRUE, colorEPOS0010);
    ProduceGraphErrWithoutXErrors(graphEPOSpredScaledforPlot_0010);
    graphEPOSpredScaledforPlot_0010->SetLineStyle(5);
    DrawGammaSetMarkerTGraphErr(graphEPOSpred_2050, 0, 0, colorEPOS2050,colorEPOS2050, 2, kTRUE, colorEPOS2050);
    ProduceGraphErrWithoutXErrors(graphEPOSpred_2050);
    graphEPOSpred_2050->SetLineStyle(5);
    DrawGammaSetMarkerTGraphErr(graphEPOSEtaToPi0_0010, 0, 0, colorEPOS0010,colorEPOS0010, 2, kTRUE, colorEPOS0010);
    ProduceGraphErrWithoutXErrors(graphEPOSEtaToPi0_0010);
    graphEPOSEtaToPi0_0010->SetLineStyle(5);
    DrawGammaSetMarkerTGraphErr(graphEPOSEtaToPi0_2050, 0, 0, colorEPOS2050,colorEPOS2050, 2, kTRUE, colorEPOS2050);
    ProduceGraphErrWithoutXErrors(graphEPOSEtaToPi0_2050);
    graphEPOSEtaToPi0_2050->SetLineStyle(5);

    //**************************************** Cracow **********************************************
    TGraphAsymmErrors *TheoryCracowPi0LowPt_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowPi0LowPt_0010");
    TGraphAsymmErrors *TheoryCracowPi0LowPt_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowPi0LowPt_2050");
    TGraphAsymmErrors *TheoryCracowEtaLowPt_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowEtaLowPt_0010");
    TGraphAsymmErrors *TheoryCracowEtaLowPt_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowEtaLowPt_2050");
    TGraphAsymmErrors* TheoryCracowEtaToPi0LowPt_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowEtaToPi0LowPt_0010");
    TGraphAsymmErrors* TheoryCracowEtaToPi0LowPt_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowEtaToPi0LowPt_2050");
    TGraphAsymmErrors *TheoryCracowChargedPionLowPt_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowChargedPionLowPt_0010");
    TGraphAsymmErrors *TheoryCracowChargedKaonLowPt_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowChargedKaonLowPt_0010");
    while(TheoryCracowPi0LowPt_0010->GetX()[TheoryCracowPi0LowPt_0010->GetN()-1]>3.) TheoryCracowPi0LowPt_0010->RemovePoint(TheoryCracowPi0LowPt_0010->GetN()-1);
    while(TheoryCracowPi0LowPt_2050->GetX()[TheoryCracowPi0LowPt_2050->GetN()-1]>3.) TheoryCracowPi0LowPt_2050->RemovePoint(TheoryCracowPi0LowPt_2050->GetN()-1);
    while(TheoryCracowEtaLowPt_0010->GetX()[TheoryCracowEtaLowPt_0010->GetN()-1]>3.) TheoryCracowEtaLowPt_0010->RemovePoint(TheoryCracowEtaLowPt_0010->GetN()-1);
    while(TheoryCracowEtaLowPt_2050->GetX()[TheoryCracowEtaLowPt_2050->GetN()-1]>3.) TheoryCracowEtaLowPt_2050->RemovePoint(TheoryCracowEtaLowPt_2050->GetN()-1);
    while(TheoryCracowEtaToPi0LowPt_0010->GetX()[TheoryCracowEtaToPi0LowPt_0010->GetN()-1]>3.) TheoryCracowEtaToPi0LowPt_0010->RemovePoint(TheoryCracowEtaToPi0LowPt_0010->GetN()-1);
    while(TheoryCracowEtaToPi0LowPt_2050->GetX()[TheoryCracowEtaToPi0LowPt_2050->GetN()-1]>3.) TheoryCracowEtaToPi0LowPt_2050->RemovePoint(TheoryCracowEtaToPi0LowPt_2050->GetN()-1);

    TGraphAsymmErrors *TheoryBegunEQPi0_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryBegunEQPi0_0010");
    TGraphAsymmErrors *TheoryBegunEQPi0_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryBegunEQPi0_2050");
    TGraphAsymmErrors *TheoryBegunEQEta_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryBegunEQEta_0010");
    TGraphAsymmErrors *TheoryBegunEQEta_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryBegunEQEta_2050");
    while(TheoryBegunEQPi0_0010->GetX()[TheoryBegunEQPi0_0010->GetN()-1]>3.) TheoryBegunEQPi0_0010->RemovePoint(TheoryBegunEQPi0_0010->GetN()-1);
    while(TheoryBegunEQPi0_2050->GetX()[TheoryBegunEQPi0_2050->GetN()-1]>3.) TheoryBegunEQPi0_2050->RemovePoint(TheoryBegunEQPi0_2050->GetN()-1);
    while(TheoryBegunEQEta_0010->GetX()[TheoryBegunEQEta_0010->GetN()-1]>3.) TheoryBegunEQEta_0010->RemovePoint(TheoryBegunEQEta_0010->GetN()-1);
    while(TheoryBegunEQEta_2050->GetX()[TheoryBegunEQEta_2050->GetN()-1]>3.) TheoryBegunEQEta_2050->RemovePoint(TheoryBegunEQEta_2050->GetN()-1);
    DrawGammaSetMarkerTGraph(TheoryBegunEQPi0_0010, 0, 0, kViolet+1,kViolet+1, 5, kTRUE, kViolet+1);
    DrawGammaSetMarkerTGraph(TheoryBegunEQPi0_2050, 0, 0, kCyan+3,kCyan+3, 5, kTRUE, kCyan+3);
    DrawGammaSetMarkerTGraph(TheoryBegunEQEta_0010, 0, 0, kPink+10,kPink+10, 5, kTRUE, kPink+10);
    DrawGammaSetMarkerTGraph(TheoryBegunEQEta_2050, 0, 0, kTeal-7,kTeal-7, 5, kTRUE,kTeal-7);
    TGraphAsymmErrors* TheoryCracowLowPt_0010 = NULL;
    TGraphAsymmErrors* TheoryCracowLowPt_2050 = NULL;
    TGraphAsymmErrors* TheoryBegunEQ_0010 = NULL;
    TGraphAsymmErrors* TheoryBegunEQ_2050 = NULL;
    if(meson.CompareTo("Pi0")==0){
        TheoryCracowLowPt_0010 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_0010->Clone("TheoryCracowLowPt_0010");
        TheoryCracowLowPt_2050 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_2050->Clone("TheoryCracowLowPt_2050");
        TheoryBegunEQ_0010 = (TGraphAsymmErrors*)TheoryBegunEQPi0_0010->Clone("TheoryBegunEQ_0010");
        TheoryBegunEQ_2050 = (TGraphAsymmErrors*)TheoryBegunEQPi0_2050->Clone("TheoryBegunEQ_2050");
    } else if(meson.CompareTo("Eta")==0){
        TheoryCracowLowPt_0010 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_0010->Clone("TheoryCracowLowPt_0010");
        TheoryCracowLowPt_2050 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_2050->Clone("TheoryCracowLowPt_2050");
        TheoryBegunEQ_0010 = (TGraphAsymmErrors*)TheoryBegunEQEta_0010->Clone("TheoryBegunEQ_0010");
        TheoryBegunEQ_2050 = (TGraphAsymmErrors*)TheoryBegunEQEta_2050->Clone("TheoryBegunEQ_2050");
    }
    TGraphAsymmErrors *TheoryCracowLowPtScaledforPlot_0010 = (TGraphAsymmErrors*)TheoryCracowLowPt_0010->Clone("TheoryCracowLowPtScaledforPlot_0010");
    TheoryCracowLowPtScaledforPlot_0010 = ScaleGraph(TheoryCracowLowPtScaledforPlot_0010,4);
    DrawGammaSetMarkerTGraphAsym(TheoryCracowLowPtScaledforPlot_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE, colorCracow0010);
    DrawGammaSetMarkerTGraphAsym(TheoryCracowLowPt_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE,colorCracow0010);
    DrawGammaSetMarkerTGraphAsym(TheoryCracowLowPt_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
    DrawGammaSetMarkerTGraphAsym(TheoryCracowPi0LowPt_0010, 0, 0, kBlue+1,kBlue+1, 5, kTRUE, kBlue+1);
    DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaLowPt_0010, 0, 0, kRed+1,kRed+1, 5, kTRUE, kRed+1);
    DrawGammaSetMarkerTGraphAsym(TheoryCracowChargedPionLowPt_0010, 0, 0, kCyan+2,kCyan+2, 5, kTRUE, kCyan+2);
    DrawGammaSetMarkerTGraphAsym(TheoryCracowChargedKaonLowPt_0010, 0, 0, kMagenta+1,kMagenta+1, 5, kTRUE, kMagenta+1);

    Int_t nBins = TheoryCracowPi0LowPt_0010->GetN();
    Double_t *xBins = TheoryCracowPi0LowPt_0010->GetX();
    Double_t *yKaons = TheoryCracowChargedKaonLowPt_0010->GetY();
    Double_t *yPions = TheoryCracowChargedPionLowPt_0010->GetY();
    Double_t *yPi0 = TheoryCracowPi0LowPt_0010->GetY();
    Double_t *yEta = TheoryCracowEtaLowPt_0010->GetY();
    Double_t *yPi0_2050 = TheoryCracowPi0LowPt_2050->GetY();
    Double_t *yEta_2050 = TheoryCracowEtaLowPt_2050->GetY();
    Double_t yKaonsToPions[nBins];
    Double_t yEtaToPi0[nBins];
    Double_t yEtaToPi0_2050[nBins];
    Double_t yEtaToKaon[nBins];
    Double_t yPi0ToPions[nBins];
    for(Int_t n=0; n<TheoryCracowPi0LowPt_0010->GetN(); n++){
      yKaonsToPions[n] = yKaons[n]/yPions[n];
      yEtaToPi0[n] = yEta[n]/yPi0[n];
      yEtaToPi0_2050[n] = yEta_2050[n]/yPi0_2050[n];
      yEtaToKaon[n] = yEta[n]/yKaons[n];
      yPi0ToPions[n] = yPi0[n]/yPions[n];
    }
    TGraphAsymmErrors *TheoryCracowChargedKaonToPionLowPt_0010 = new TGraphAsymmErrors(nBins,xBins,yKaonsToPions,0,0,0,0);
    TGraphAsymmErrors *TheoryCracowEtaToPi0_0010 = new TGraphAsymmErrors(nBins,xBins,yEtaToPi0,0,0,0,0);
    TGraphAsymmErrors *TheoryCracowEtaToPi0_2050 = new TGraphAsymmErrors(nBins,xBins,yEtaToPi0_2050,0,0,0,0);
    TGraphAsymmErrors *TheoryCracowChargedEtaToKaonLowPt_0010 = new TGraphAsymmErrors(nBins,xBins,yEtaToKaon,0,0,0,0);
    TGraphAsymmErrors *TheoryCracowPi0ToChargedPionLowPt_0010 = new TGraphAsymmErrors(nBins,xBins,yPi0ToPions,0,0,0,0);

    while(TheoryBegunEQ_0010->GetX()[TheoryBegunEQ_0010->GetN()-1]>3.) TheoryBegunEQ_0010->RemovePoint(TheoryBegunEQ_0010->GetN()-1);
    while(TheoryBegunEQ_2050->GetX()[TheoryBegunEQ_2050->GetN()-1]>3.) TheoryBegunEQ_2050->RemovePoint(TheoryBegunEQ_2050->GetN()-1);
    TGraphAsymmErrors *TheoryBegunEQ_0010ScaledforPlot = (TGraphAsymmErrors*)TheoryBegunEQ_0010->Clone("TheoryBegunEQ_0010ScaledforPlot");
    TheoryBegunEQ_0010ScaledforPlot = ScaleGraph(TheoryBegunEQ_0010ScaledforPlot,4);
    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQ_0010ScaledforPlot, 0, 0, colorBegun0010,colorBegun0010, 5, kTRUE, colorBegun0010);
    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQ_0010, 0, 0, colorBegun0010,colorBegun0010, 5, kTRUE, colorBegun0010);
    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQ_2050, 0, 0, colorBegun2050,colorBegun2050, 5, kTRUE, colorBegun2050);
    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQPi0_0010, 0, 0, colorBegun0010,colorBegun0010, 5, kTRUE, colorBegun0010);
    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQPi0_2050, 0, 0, colorBegun2050,colorBegun2050, 5, kTRUE, colorBegun2050);
    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQEta_0010, 0, 0, kOrange-4,kOrange-4, 5, kTRUE, kOrange-4);
    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQEta_2050, 0, 0, kAzure+6,kAzure+6, 5, kTRUE,kAzure+6);

    TGraphAsymmErrors* TheoryBegunEQEtaToPi0_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryBegunEQEtaToPi0_0010");
    TGraphAsymmErrors* TheoryBegunEQEtaToPi0_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryBegunEQEtaToPi0_2050");
    while(TheoryBegunEQEtaToPi0_0010->GetX()[TheoryBegunEQEtaToPi0_0010->GetN()-1]>3.) TheoryBegunEQEtaToPi0_0010->RemovePoint(TheoryBegunEQEtaToPi0_0010->GetN()-1);
    while(TheoryBegunEQEtaToPi0_2050->GetX()[TheoryBegunEQEtaToPi0_2050->GetN()-1]>3.) TheoryBegunEQEtaToPi0_2050->RemovePoint(TheoryBegunEQEtaToPi0_2050->GetN()-1);
    DrawGammaSetMarkerTGraph(TheoryBegunEQEtaToPi0_0010, 0, 0, kViolet+1,kViolet+1, 5, kTRUE, kViolet+1);
    DrawGammaSetMarkerTGraph(TheoryBegunEQEtaToPi0_2050, 0, 0, kCyan+1,kCyan+1, 5, kTRUE, kCyan+1);


    //**************************** Jet Quenching NLO arXiv:1506.00838 ****************************//
    TGraphAsymmErrors* graphPi0RAAJetQuenching_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphPi0RAAJetQuenching_0010");
    TGraphAsymmErrors* graphEtaRAAJetQuenching_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphEtaRAAJetQuenching_0010");
    TGraphAsymmErrors* graphEtaToPi0JetQuenching_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphEtaToPi0JetQuenching_0010");


    //**************************** WHDG ****************************//
    gWHDG_Pi0_Raa_0510 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0510");
    gWHDG_Pi0_Raa_0005 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0005");
    gWHDG_Pi0_Raa_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0010");
    gWHDG_Pi0_Raa_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA2050");
    gWHDG_Eta_Raa_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGetaRAA0010");
    gWHDG_Eta_Raa_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGetaRAA2050");


    //**************************** Djordjevic ****************************//
    TGraphAsymmErrors* graphPi0Djordjevic_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphPi0Djordjevic_0010");
    TGraphAsymmErrors* graphPi0Djordjevic_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphPi0Djordjevic_2050");



    //*********************************************************************************************************//
    //********************************        Other experiment inputs     *************************************//
    //*********************************************************************************************************//
    TFile* fileDataOtherEnergies = new TFile(fileNameDataOtherEnergyInput);
	TGraphErrors* graphPHENIX200GeVInvYield_0010 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVInvYield_0010");
	TGraphErrors* graphPHENIX200GeVInvYield_2040 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVInvYield_2040");
	TGraphErrors* graphPHENIX200GeVInvYield_2050 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVInvYield_2050");
	TGraphErrors* graphPHENIX200GeVEtaInvYield_0020 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaInvYield_0020");
	TGraphErrors* graphPHENIX200GeVEtaInvYield_2060 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaInvYield_2060");

    TGraphErrors* graphPHENIX39GeVInvYield_0010 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVInvYield_0010");
	TGraphErrors* graphPHENIX39GeVInvYield_2040 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVInvYield_2040");
	TGraphErrors* graphPHENIX62GeVInvYield_0010 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVInvYield_0010");
	TGraphErrors* graphPHENIX62GeVInvYield_2040 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVInvYield_2040");

    graphWA98_17_3GeVPi0RAA_0013 = (TGraphErrors*)fileDataOtherEnergies->Get("graphWA98RAA_0013");
    graphPHENIX200GeVPi0RAA_0010 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_0010");
    graphPHENIX200GeVPi0RAA_2040 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_2040");
    graphPHENIX39GeVPi0RAA_0010 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_0010");
    graphPHENIX39GeVPi0RAA_2040 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_2040");
    graphPHENIX62GeVPi0RAA_0010 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_0010");
    graphPHENIX62GeVPi0RAA_2040 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_2040");

    graphPHENIX200GeVEtaToPi0Ratio_0020 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaToPi0Ratio_0020");
    graphPHENIX200GeVEtaToPi0Ratio_2060 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaToPi0Ratio_2060");

    graphPHENIX200GeVEtaRAA_0020 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_0020");
    graphPHENIX200GeVEtaRAA_2060 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_2060");
    graphPHENIX200GeVEtaRAA_0010 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_0010");
    graphPHENIX200GeVEtaRAA_2040 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_2040");
    graphPHENIX200GeVEtaHighPtRAA_2060 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAAHighPt_2060");


    TFile* fileDataALICE = new TFile(fileNameALICEData);
    //*********************************************************************************************************//
    //**********************************     Charged pions and kaons     **************************************//
    //*********************************************************************************************************//
    TDirectoryFile* directoryChargedPbPb = (TDirectoryFile*)fileDataALICE->Get("ChargedParticles_PbPb_2.76TeV");

    histoChargedPionSpectraStat0010   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraStat0010");
    histoChargedPionSpectraSyst0010   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraSyst0010");
    histoChargedPionSpectraStat2040   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraStat2040");
    histoChargedPionSpectraSyst2040   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraSyst2040");

    histoChargedKaonSpectraStat0010   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraStat0010");
    histoChargedKaonSpectraSyst0010   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraSyst0010");
    histoChargedKaonSpectraStat2040   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraStat2040");
    histoChargedKaonSpectraSyst2040   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraSyst2040");

    graphChargedRatioKaonToPion0010      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPion0010");
    graphChargedRatioKaonToPionSys0010   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPionSys0010");
    DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPion0010, 24,markerSizeComb, colorCharged, colorCharged);
    DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPionSys0010,24,markerSizeComb, colorCharged, colorCharged, 1, kTRUE);
    graphChargedRatioKaonToPion2040      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPion2040");
    graphChargedRatioKaonToPionSys2040   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPionSys2040");

    graphChargedPionRAA0010      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAA0010");
    graphChargedPionRAASys0010   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAASys0010");
    graphChargedPionRAA2040      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAA2040");
    graphChargedPionRAASys2040   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAASys2040");
    graphChargedKaonRAA0010      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAA0010");
    graphChargedKaonRAASys0010   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAASys0010");
    graphChargedKaonRAA2040      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAA2040");
    graphChargedKaonRAASys2040   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAASys2040");

    while(graphChargedKaonRAA0010->GetX()[0] < .3){
      graphChargedKaonRAA0010->RemovePoint(0);
      graphChargedKaonRAASys0010->RemovePoint(0);
      graphChargedKaonRAA2040->RemovePoint(0);
      graphChargedKaonRAASys2040->RemovePoint(0);
      graphChargedPionRAA0010->RemovePoint(0);
      graphChargedPionRAASys0010->RemovePoint(0);
      graphChargedPionRAA2040->RemovePoint(0);
      graphChargedPionRAASys2040->RemovePoint(0);
    }

    TGraphAsymmErrors* graphChargedPionRAATotErr0010 = AddErrorsOfGraphsQuadratically(graphChargedPionRAA0010,graphChargedPionRAASys0010);
        DrawGammaSetMarkerTGraphAsym(graphChargedPionRAATotErr0010, 27,3, colorCharged,colorCharged, 0.1, kFALSE);
    TGraphAsymmErrors* graphChargedKaonRAATotErr0010 = AddErrorsOfGraphsQuadratically(graphChargedKaonRAA0010,graphChargedKaonRAASys0010);
         DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAATotErr0010, 33,3, colorCharged,colorCharged, 0.1, kFALSE);
   TGraphAsymmErrors* graphChargedPionRAATotErr2040 = AddErrorsOfGraphsQuadratically(graphChargedPionRAA2040,graphChargedPionRAASys2040);
        DrawGammaSetMarkerTGraphAsym(graphChargedPionRAATotErr2040, 27,3, colorCharged,colorCharged, 0.1, kFALSE);
    TGraphAsymmErrors* graphChargedKaonRAATotErr2040 = AddErrorsOfGraphsQuadratically(graphChargedKaonRAA2040,graphChargedKaonRAASys2040);
        DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAATotErr2040, 33,3, colorCharged,colorCharged, 0.1, kFALSE);

    //*********************************************************************************************************//
    //*****************************************  Neutral mesons   *********************************************//
    //*********************************************************************************************************//

    TDirectoryFile* directoryNeutralMesonPP7TeV = (TDirectoryFile*)fileDataALICE->Get("NeutralMesons_PP_7TeV");
        graphCombEtaToPi0Ratiopp7TeVNoXErrors = (TGraphAsymmErrors*)directoryNeutralMesonPP7TeV->Get("graphCombEtaToPi0Ratiopp7TeVNoXErrors");
        graphCombEtaToPi0RatioSysErrpp7TeV = (TGraphAsymmErrors*)directoryNeutralMesonPP7TeV->Get("graphCombEtaToPi0RatioSysErrpp7TeV");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RatioSysErrpp7TeV, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Ratiopp7TeVNoXErrors, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);

    TDirectoryFile* directoryNeutralMesonPP = (TDirectoryFile*)fileDataALICE->Get("NeutralMesons_PP_2.76TeV");
        //already scaled by 1./xSection2760GeVppINEL
        graphInvSectionCombStatPi02760GeV           = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0Comb2760GeVAStatErr");
        graphInvSectionCombSysPi02760GeV            = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0Comb2760GeVASysErr");
        graphInvSectionCombStatPi02760GeVPlot       = (TGraphAsymmErrors*)graphInvSectionCombStatPi02760GeV->Clone("graphInvSectionCombStatPi02760GeVPlot");
        graphInvSectionCombSysPi02760GeVPlot        = (TGraphAsymmErrors*)graphInvSectionCombSysPi02760GeV->Clone("graphInvSectionCombSysPi02760GeVPlot");
        DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatPi02760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvSectionCombStatPi02760GeVPlot);
        DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysPi02760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);
        TF1 *fitInvCrossSectionTsallisPi0Comb2760GeV = (TF1*)directoryNeutralMesonPP->Get("TsallisFitPi0");// ---> Tsallis fit is used for both bin-shifts
        TF1 *fitInvCrossSectionTCMPi0Comb2760GeV     = (TF1*)directoryNeutralMesonPP->Get("TwoComponentModelFitPi0");
        graphInvSectionCombStatPi02760GeVforRAA     = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0Comb2760GeVAStatErr_yShifted");
        graphInvSectionCombSysPi02760GeVforRAA      = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0Comb2760GeVASysErr_yShifted");

        graphInvSectionPCMStatPi02760GeV            = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
        graphInvSectionPCMSysPi02760GeV             = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
        graphInvSectionPCMStatPi02760GeVforRAA      = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PCM2760GeVStatErr_yShifted");
        graphInvSectionPCMSysPi02760GeV_yShifted    = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PCM2760GeVSysErr_yShifted");//WITH material budget error
        graphInvSectionPCMSysPi02760GeVforRAA       = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PCM2760GeVSysErr_forRAA");//WITHOUT common errors

        graphInvSectionPHOSStatPi02760GeV           = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PHOS2760GeVStatErr");
        graphInvSectionPHOSSysPi02760GeV            = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr");
        graphInvSectionPHOSStatPi02760GeVforRAA     = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PHOS2760GeVStatErr_yShifted");
        graphInvSectionPHOSSysPi02760GeV_yShifted   = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr_yShifted");//WITH common errors
        graphInvSectionPHOSSysPi02760GeVforRAA      = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr_forRAA");//WITHOUT common errors

        graphInvSectionEMCalStatPi02760GeV          = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0EMCAL2760GeVStatErr");
        graphInvSectionEMCalSysPi02760GeV           = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0EMCAL2760GeVSysErr");
        graphInvSectionEMCalStatPi02760GeVforRAA    = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0EMCAL2760GeVStatErr_yShifted");
        graphInvSectionEMCalSysPi02760GeV_yShifted  = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0EMCAL2760GeVSysErr_yShifted");//WITH common errors
        graphInvSectionEMCalSysPi02760GeVforRAA     = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionPi0EMCAL2760GeVSysErr_forRAA");//WITHOUT common errors

        graphInvSectionCombStatEta2760GeV           = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaComb2760GeVAStatErr");
        graphInvSectionCombSysEta2760GeV            = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaComb2760GeVASysErr");
        graphInvSectionCombStatEta2760GeVPlot       = (TGraphAsymmErrors*)graphInvSectionCombStatEta2760GeV->Clone("graphInvSectionCombStatEta2760GeVPlot");
        graphInvSectionCombSysEta2760GeVPlot        = (TGraphAsymmErrors*)graphInvSectionCombSysEta2760GeV->Clone("graphInvSectionCombSysEta2760GeVPlot");
        DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatEta2760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvSectionCombStatEta2760GeVPlot);
        DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysEta2760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);
        TF1 *fitInvCrossSectionTsallisEtaComb2760GeV = (TF1*)directoryNeutralMesonPP->Get("TsallisFitEta");// ---> Tsallis fit is used for both bin-shifts
        TF1 *fitInvCrossSectionTCMEtaComb2760GeV     = (TF1*)directoryNeutralMesonPP->Get("TwoComponentModelFitEta");
        graphInvSectionCombStatEta2760GeVforRAA     = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaComb2760GeVAStatErr_yShifted");
        graphInvSectionCombSysEta2760GeVforRAA      = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaComb2760GeVASysErr_yShifted");

        graphInvSectionPCMStatEta2760GeV            = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaPCM2760GeVStatErr");
        graphInvSectionPCMSysEta2760GeV             = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaPCM2760GeVSysErr");
        graphInvSectionPCMStatEta2760GeVforRAA      = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaPCM2760GeVStatErr_yShifted");
        graphInvSectionPCMSysEta2760GeV_yShifted    = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaPCM2760GeVSysErr_yShifted");//WITH material budget error
        graphInvSectionPCMSysEta2760GeVforRAA       = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaPCM2760GeVSysErr_forRAA");//WITHOUT material budget error

        graphInvSectionEMCalStatEta2760GeV          = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaEMCAL2760GeVStatErr");
        graphInvSectionEMCalSysEta2760GeV           = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaEMCAL2760GeVSysErr");
        graphInvSectionEMCalStatEta2760GeVforRAA    = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaEMCAL2760GeVStatErr_yShifted");
        graphInvSectionEMCalSysEta2760GeV_yShifted  = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaEMCAL2760GeVSysErr_yShifted");//WITH common errors
        graphInvSectionEMCalSysEta2760GeVforRAA     = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphInvCrossSectionEtaEMCAL2760GeVSysErr_forRAA");//WITHOUT common errors

        graphRatioEtaToPi0Comb2760GeVStatErr        = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphRatioEtaToPi0Comb2760GeVStatErr");
        graphRatioEtaToPi0Comb2760GeVSysErr         = (TGraphAsymmErrors*)directoryNeutralMesonPP->Get("graphRatioEtaToPi0Comb2760GeVSysErr");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaToPi0Comb2760GeVSysErr, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaToPi0Comb2760GeVStatErr, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioEtaToPi0Comb2760GeVStatErr);

    TDirectoryFile* directoryNeutralMesonPbPb   = (TDirectoryFile*)fileDataALICE->Get("NeutralMesons_PbPb_2.76TeV");

        TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_2040          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVStatErr_2040");
        TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_2040           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVSysErr_2040");
        TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_2040          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVStatErr_2040");
        TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_2040           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVSysErr_2040");

        TF1 *etapi0RatioFromMtScalingPCMPbPb2760GeV_0010 = (TF1*)directoryNeutralMesonPbPb->Get("etapi0RatioFromMtScalingPCMPbPb2760GeV_0010");
        TF1 *etapi0RatioFromMtScalingCombPbPb2760GeV_0010 = (TF1*)directoryNeutralMesonPbPb->Get("etapi0RatioFromMtScalingCombPbPb2760GeV_0010");
        TF1 *etapi0RatioFromMtScalingPCMPbPb2760GeV_2050 = (TF1*)directoryNeutralMesonPbPb->Get("etapi0RatioFromMtScalingPCMPbPb2760GeV_2050");
        TF1 *etapi0RatioFromMtScalingCombPbPb2760GeV_2050 = (TF1*)directoryNeutralMesonPbPb->Get("etapi0RatioFromMtScalingCombPbPb2760GeV_2050");

        graphPCMPubPi0InvYieldStatPbPb2760GeV_0010    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPubPbPb2760GeVStatErr_0010");
        graphPCMPubPi0InvYieldSysPbPb2760GeV_0010     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPubPbPb2760GeVSysErr_0010");
        graphPCMPubPi0InvYieldStatPbPb2760GeV_2040    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPubPbPb2760GeVStatErr_2040");
        graphPCMPubPi0InvYieldSysPbPb2760GeV_2040     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPubPbPb2760GeVSysErr_2040");

        histoPi0PHOSPbPb0010                        = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0PHOSPbPb2760GeVStatErr_0010");
        graphPHOSPi0InvYieldStatPbPb2760GeV_0010    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PHOSPbPb2760GeVStatErr_0010");
        graphPHOSPi0InvYieldSysPbPb2760GeV_0010     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PHOSPbPb2760GeVSysErr_0010");
        graphSysErrRAAYieldPi0PHOSPbPb0010          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PHOSPbPb2760GeVSysErr_forRAA_0010");

//         cout << __LINE__ << endl;
//         graphPHOSPi0InvYieldStatPbPb2760GeV_0010->Print();
//         graphPHOSPi0InvYieldSysPbPb2760GeV_0010->Print();
//         graphSysErrRAAYieldPi0PHOSPbPb0010->Print();
//         return;

        histoPCMPi0InvYieldPbPb2760GeV_0010         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0PCMPbPb2760GeVStatErr_0010");
        graphPCMPi0InvYieldStatPbPb2760GeV_0010     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_0010");
        graphPCMPi0InvYieldSysPbPb2760GeV_0010      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_0010");
        graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_yShifted_0010");
        graphPCMPi0InvYieldSysWOMat2760GeV_0010     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_forRAA_0010");
        histoPCMPi0RAAStatPbPb2760GeV_0010          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAPi0PCMPbPb2760GeVStatErr_0010");
        graphPCMPi0RAAStatPbPb2760GeV_0010          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVStatErr_0010");
        graphPCMPi0RAASysPbPb2760GeV_0010           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVSysErr_0010");

//         cout << __LINE__ << endl;
//         graphPCMPi0InvYieldStatPbPb2760GeV_0010->Print();
//         graphPCMPi0InvYieldSysPbPb2760GeV_0010->Print();
//         graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010->Print();
//         graphPCMPi0InvYieldSysWOMat2760GeV_0010->Print();
//         graphPCMPi0RAAStatPbPb2760GeV_0010->Print();
//         graphPCMPi0RAASysPbPb2760GeV_0010->Print();
//         return;

        histoPCMPi0InvYieldPbPb2760GeV_2050         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0PCMPbPb2760GeVStatErr_2050");
        graphPCMPi0InvYieldStatPbPb2760GeV_2050     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_2050");
        graphPCMPi0InvYieldSysPbPb2760GeV_2050      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_2050");
        graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_yShifted_2050");
        graphPCMPi0InvYieldSysWOMat2760GeV_2050     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_forRAA_2050");
        histoPCMPi0RAAStatPbPb2760GeV_2050          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAPi0PCMPbPb2760GeVStatErr_2050");
        graphPCMPi0RAAStatPbPb2760GeV_2050          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVStatErr_2050");
        graphPCMPi0RAASysPbPb2760GeV_2050           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVSysErr_2050");

//         cout << __LINE__ << endl;
//         graphPCMPi0InvYieldStatPbPb2760GeV_2050->Print();
//         graphPCMPi0InvYieldSysPbPb2760GeV_2050->Print();
//         graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050->Print();
//         graphPCMPi0InvYieldSysWOMat2760GeV_2050->Print();
//         graphPCMPi0RAAStatPbPb2760GeV_2050->Print();
//         graphPCMPi0RAASysPbPb2760GeV_2050->Print();
//         return;

        histoPCMPi0InvYieldPbPb2760GeV_2040         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0PCMPbPb2760GeVStatErr_2040");
        graphPCMPi0InvYieldStatPbPb2760GeV_2040     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_2040");
        graphPCMPi0InvYieldSysPbPb2760GeV_2040      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_2040");

        histoPCMEtaInvYieldPbPb2760GeV_0010         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaPCMPbPb2760GeVStatErr_0010");
        graphPCMEtaInvYieldStatPbPb2760GeV_0010     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_0010");
        graphPCMEtaInvYieldSysPbPb2760GeV_0010      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_0010");
        graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_yShifted_0010");
        graphPCMEtaInvYieldSysWOMat2760GeV_0010     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_forRAA_0010");
        histoPCMEtaRAAStatPbPb2760GeV_0010          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAEtaPCMPbPb2760GeVStatErr_0010");
        graphPCMEtaRAAStatPbPb2760GeV_0010          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVStatErr_0010");
        graphPCMEtaRAASysPbPb2760GeV_0010           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVSysErr_0010");
        histoPCMEtatoPi0Stat2760GeV_0010            = (TH1D*)directoryNeutralMesonPbPb->Get("histoEtaToPi0RatioPCMPbPb2760GeVStatErr_0010");
        graphPCMEtatoPi0Stat2760GeV_0010            = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVStatErr_0010");
        graphPCMEtatoPi0Sys2760GeV_0010             = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVSysErr_0010");

        histoPCMEtaInvYieldPbPb2760GeV_2050         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaPCMPbPb2760GeVStatErr_2050");
        graphPCMEtaInvYieldStatPbPb2760GeV_2050     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_2050");
        graphPCMEtaInvYieldSysPbPb2760GeV_2050      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_2050");
        graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_yShifted_2050");
        graphPCMEtaInvYieldSysWOMat2760GeV_2050     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_forRAA_2050");
        histoPCMEtaRAAStatPbPb2760GeV_2050          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAEtaPCMPbPb2760GeVStatErr_2050");
        graphPCMEtaRAAStatPbPb2760GeV_2050          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVStatErr_2050");
        graphPCMEtaRAASysPbPb2760GeV_2050           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVSysErr_2050");
        histoPCMEtatoPi0Stat2760GeV_2050            = (TH1D*)directoryNeutralMesonPbPb->Get("histoEtaToPi0RatioPCMPbPb2760GeVStatErr_2050");
        graphPCMEtatoPi0Stat2760GeV_2050            = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVStatErr_2050");
        graphPCMEtatoPi0Sys2760GeV_2050             = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVSysErr_2050");

        histoPCMEtaInvYieldPbPb2760GeV_2040         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaPCMPbPb2760GeVStatErr_2040");
        graphPCMEtaInvYieldStatPbPb2760GeV_2040     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_2040");
        graphPCMEtaInvYieldSysPbPb2760GeV_2040      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_2040");

        histoEMCalPi0InvYieldStatPbPb2760GeV_0010   = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0EMCalPbPb2760GeVStatErr_0010");
        graphEMCalPi0InvYieldStatPbPb2760GeV_0010   = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVStatErr_0010");
        graphEMCalPi0InvYieldSysPbPb2760GeV_0010    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErr_0010");
        graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_0010 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErr_forRAA_0010");
        histoEMCalPi0RAAStatPbPb2760GeV_0010        = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAPi0EMCalPbPb2760GeVStatErr_0010");
//         graphEMCalPi0RAAStatPbPb2760GeV_0010        = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0EMCalPbPb2760GeVStatErr_0010");
//         graphEMCalPi0RAASysPbPb2760GeV_0010         = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0EMCalPbPb2760GeVSysErr_0010");

        histoEMCalPi0InvYieldStatPbPb2760GeV_2050   = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0EMCalPbPb2760GeVStatErr_2050");
        graphEMCalPi0InvYieldStatPbPb2760GeV_2050   = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVStatErr_2050");
        graphEMCalPi0InvYieldSysPbPb2760GeV_2050    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErr_2050");
        graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_2050 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErr_forRAA_2050");
        histoEMCalPi0RAAStatPbPb2760GeV_2050        = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAPi0EMCalPbPb2760GeVStatErr_2050");
//         graphEMCalPi0RAAStatPbPb2760GeV_2050        = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0EMCalPbPb2760GeVStatErr_2050");
//         graphEMCalPi0RAASysPbPb2760GeV_2050         = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0EMCalPbPb2760GeVSysErr_2050");

        histoEMCalEtaInvYieldStatPbPb2760GeV_0010   = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtaInvYieldStatPbPb2760GeV_0010   = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtaInvYieldSysPbPb2760GeV_0010    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErr_0010");
        graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErr_forRAA_0010");
        histoEMCalEtaRAAStatPbPb2760GeV_0010        = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAEtaEMCalPbPb2760GeVStatErr_0010");
//         graphEMCalEtaRAAStatPbPb2760GeV_0010        = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaEMCalPbPb2760GeVStatErr_0010");
//         graphEMCalEtaRAASysPbPb2760GeV_0010         = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaEMCalPbPb2760GeVSysErr_0010");
        histoEMCalEtatoPi0Stat2760GeV_0010          = (TH1D*)directoryNeutralMesonPbPb->Get("histoEtaToPi0RatioEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtatoPi0Stat2760GeV_0010          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtatoPi0Sys2760GeV_0010           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioEMCalPbPb2760GeVSysErr_0010");

        histoEMCalEtaInvYieldStatPbPb2760GeV_2050   = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtaInvYieldStatPbPb2760GeV_2050   = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtaInvYieldSysPbPb2760GeV_2050    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErr_2050");
        graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_2050 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErr_forRAA_2050");
        histoEMCalEtaRAAStatPbPb2760GeV_2050        = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAEtaEMCalPbPb2760GeVStatErr_2050");
//         graphEMCalEtaRAAStatPbPb2760GeV_2050        = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaEMCalPbPb2760GeVStatErr_2050");
//         graphEMCalEtaRAASysPbPb2760GeV_2050         = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaEMCalPbPb2760GeVSysErr_2050");
        histoEMCalEtatoPi0Stat2760GeV_2050          = (TH1D*)directoryNeutralMesonPbPb->Get("histoEtaToPi0RatioEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtatoPi0Stat2760GeV_2050          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtatoPi0Sys2760GeV_2050           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioEMCalPbPb2760GeVSysErr_2050");


        TH1D * ratioRCSC = (TH1D*)histoEMCalPi0InvYieldStatPbPb2760GeV_0010->Clone();
        ratioRCSC->Divide(histoEMCalPi0InvYieldStatPbPb2760GeV_0010,histoEMCalPi0InvYieldStatPbPb2760GeV_2050,nColl2050,nColl0010,"");

    // *******************************************************************************************************
    // ************************** Combination of different measurements **************************************
    // *******************************************************************************************************
    // REMARKS:
    //      - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //      - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //      - currently only PCM - EMCal vs others fully implemeted energy independent
    //      - extendable to other energies
    //      - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented


    // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
    // for statistical error and final value from respective method
    for (Int_t i = 0; i< 11; i++){
        statErrorCollectionLHC11h_0010[i] = NULL;
        statErrorCollectionLHC11h_2050[i] = NULL;
        statErrorCollectionEtatoPi0LHC11h_0010[i] = NULL;
        statErrorCollectionEtatoPi0LHC11h_2050[i] = NULL;
        statErrorCollectionRaaLHC11h_0010[i] = NULL;
        statErrorCollectionRaaLHC11h_2050[i] = NULL;
    }
    if(meson.CompareTo("Pi0")==0){

        statErrorCollectionLHC11h_0010[0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_0010->Clone("statErrPCMPi0_0010");
        statErrorCollectionLHC11h_0010[1] = (TH1D*)histoPi0PHOSPbPb0010->Clone("statErrPHOSPi0_0010");
        statErrorCollectionLHC11h_0010[2] = (TH1D*)histoEMCalPi0InvYieldStatPbPb2760GeV_0010->Clone("statErrEMCalPi0_0010");
        statErrorCollectionLHC11h_2050[0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_2050->Clone("statErrPCMPi0_2050");
        statErrorCollectionLHC11h_2050[2] = (TH1D*)histoEMCalPi0InvYieldStatPbPb2760GeV_2050->Clone("statErrEMCalPi0_2050");

        statErrorCollectionRaaLHC11h_0010[0] = (TH1D*)histoPCMPi0RAAStatPbPb2760GeV_0010->Clone("statErrPCMPi0RAA_0010");
        statErrorCollectionRaaLHC11h_0010[2] = (TH1D*)histoEMCalPi0RAAStatPbPb2760GeV_0010->Clone("statErrEMCalPi0RAA_0010");
        statErrorCollectionRaaLHC11h_2050[0] = (TH1D*)histoPCMPi0RAAStatPbPb2760GeV_2050->Clone("statErrPCMPi0RAA_2050");
        statErrorCollectionRaaLHC11h_2050[2] = (TH1D*)histoEMCalPi0RAAStatPbPb2760GeV_2050->Clone("statErrEMCalPi0RAA_2050");

    } else if(meson.CompareTo("Eta")==0) {

        statErrorCollectionLHC11h_0010[0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_0010->Clone("statErrPCMEta_0010");
        statErrorCollectionLHC11h_0010[2] = (TH1D*)histoEMCalEtaInvYieldStatPbPb2760GeV_0010->Clone("statErrEMCalEta_0010");
        statErrorCollectionLHC11h_2050[0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_2050->Clone("statErrPCMEta_2050");
        statErrorCollectionLHC11h_2050[2] = (TH1D*)histoEMCalEtaInvYieldStatPbPb2760GeV_2050->Clone("statErrEMCalEta_2050");

        statErrorCollectionEtatoPi0LHC11h_0010[0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_0010->Clone("statErrPCMEtatoPi0_0010");
        statErrorCollectionEtatoPi0LHC11h_0010[2] = (TH1D*)histoEMCalEtatoPi0Stat2760GeV_0010->Clone("statErrEMCalEtatoPi0_0010");
        statErrorCollectionEtatoPi0LHC11h_2050[0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_2050->Clone("statErrPCMEtatoPi0_2050");
        statErrorCollectionEtatoPi0LHC11h_2050[2] = (TH1D*)histoEMCalEtatoPi0Stat2760GeV_2050->Clone("statErrEMCalEtatoPi0_2050");

        statErrorCollectionRaaLHC11h_0010[0] = (TH1D*)histoPCMEtaRAAStatPbPb2760GeV_0010->Clone("statErrPCMEtaRAA_0010");
        statErrorCollectionRaaLHC11h_0010[2] = (TH1D*)histoEMCalEtaRAAStatPbPb2760GeV_0010->Clone("statErrEMCalEtaRAA_0010");
        statErrorCollectionRaaLHC11h_2050[0] = (TH1D*)histoPCMEtaRAAStatPbPb2760GeV_2050->Clone("statErrPCMEtaRAA_2050");
        statErrorCollectionRaaLHC11h_2050[2] = (TH1D*)histoEMCalEtaRAAStatPbPb2760GeV_2050->Clone("statErrEMCalEtaRAA_2050");

    }
    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
    // for systematic error from respective method
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollectionLHC11h_0010[i] = NULL;
        sysErrorCollectionLHC11h_2050[i] = NULL;

        sysErrorCollectionforRatiosLHC11h_0010[i] = NULL;
        sysErrorCollectionforRatiosLHC11h_2050[i] = NULL;

        sysErrorCollectionEtatoPi0LHC11h_0010[i] = NULL;
        sysErrorCollectionEtatoPi0LHC11h_2050[i] = NULL;

        sysErrorCollectionRaaLHC11h_0010[i] = NULL;
        sysErrorCollectionRaaLHC11h_2050[i] = NULL;

    }
    if(meson.CompareTo("Pi0")==0){

        sysErrorCollectionLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0010->Clone("sysErrPCMPi0_0010");
        sysErrorCollectionLHC11h_0010[1] = (TGraphAsymmErrors*)graphPHOSPi0InvYieldSysPbPb2760GeV_0010->Clone("sysErrPHOSPi0_0010");
        sysErrorCollectionLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Clone("sysErrEMCalPi0_0010");
        sysErrorCollectionLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2050->Clone("sysErrPCMPi0_2050");
        sysErrorCollectionLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Clone("sysErrEMCalPi0_2050");

        sysErrorCollectionforRatiosLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysWOMat2760GeV_0010->Clone("sysErrPCMWOMatPi0_0010");
        sysErrorCollectionforRatiosLHC11h_0010[1] = (TGraphAsymmErrors*)graphSysErrRAAYieldPi0PHOSPbPb0010->Clone("sysErrPHOSTypeBPi0_0010");
        sysErrorCollectionforRatiosLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_0010->Clone("sysErrEMCalforRAAPi0_0010");
        sysErrorCollectionforRatiosLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysWOMat2760GeV_2050->Clone("sysErrPCMWOMatPi0_2050");
        sysErrorCollectionforRatiosLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_2050->Clone("sysErrEMCalforRAAPi0_2050");

//         sysErrorCollectionRaaLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_0010->Clone("sysErrPCMPi0Raa_0010");
//         sysErrorCollectionRaaLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalPi0RAASysPbPb2760GeV_0010->Clone("sysErrEMCalPi0Raa_0010");
//         sysErrorCollectionRaaLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_2050->Clone("sysErrPCMPi0Raa_2050");
//         sysErrorCollectionRaaLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalPi0RAASysPbPb2760GeV_2050->Clone("sysErrEMCalPi0Raa_2050");

    } else if(meson.CompareTo("Eta")==0) {

        sysErrorCollectionLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0010->Clone("sysErrPCMEta_0010");
        sysErrorCollectionLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Clone("sysErrEMCalEta_0010");
        sysErrorCollectionLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2050->Clone("sysErrPCMEta_2050");
        sysErrorCollectionLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Clone("sysErrEMCalEta_2050");

        sysErrorCollectionforRatiosLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysWOMat2760GeV_0010->Clone("sysErrPCMWOMatEta_0010");
        sysErrorCollectionforRatiosLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010->Clone("sysErrEMCalforRAAEta_0010");
        sysErrorCollectionforRatiosLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysWOMat2760GeV_2050->Clone("sysErrPCMWOMatEta_2050");
        sysErrorCollectionforRatiosLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_2050->Clone("sysErrEMCalforRAAEta_2050");

        sysErrorCollectionEtatoPi0LHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_0010->Clone("sysErrPCMEtatoPi0_0010");
        sysErrorCollectionEtatoPi0LHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_0010->Clone("sysErrEMCalEtatoPi0_0010");
        sysErrorCollectionEtatoPi0LHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_2050->Clone("sysErrPCMEtatoPi0_2050");
        sysErrorCollectionEtatoPi0LHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_2050->Clone("sysErrEMCalEtatoPi0_2050");

//         sysErrorCollectionRaaLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_0010->Clone("sysErrPCMEtaRaa_0010");
//         sysErrorCollectionRaaLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalEtaRAASysPbPb2760GeV_0010->Clone("sysErrEMCalEtaRaa_0010");
//         sysErrorCollectionRaaLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_2050->Clone("sysErrPCMEtaRaa_2050");
//         sysErrorCollectionRaaLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalEtaRAASysPbPb2760GeV_2050->Clone("sysErrEMCalEtaRaa_2050");

    }
    cout << __LINE__ << endl;

    Int_t textSizeLabelsPixel = 900*0.04;
    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();

    TH2F * histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,0.23,70.,1000,-0.6,1.1);
    SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DWeights->GetXaxis()->SetMoreLogLabels();
    histo2DWeights->GetXaxis()->SetLabelOffset(-0.01);
    histo2DWeights->Draw("copy");

    TLatex *labelWeightsEnergy = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
    SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
    labelWeightsEnergy->SetTextFont(43);

    if(meson.CompareTo("Pi0")==0){
        labelWeightsPi0 = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
    } else  if(meson.CompareTo("Eta")==0){
        labelWeightsPi0 = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
    }
    SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
    labelWeightsPi0->SetTextFont(43);

    //**********************************************************************************************************************//
    //**************************************** Combination WITHOUT common errors ********************************************//
    //**********************************************************************************************************************//
    for (Int_t i = 0; i< 11; i++){
            graphWeightsALHC11h_0010[i] = NULL;
            graphWeightsALHC11h_2050[i] = NULL;
    }

    TString fileNameOutputWeightingALHC11h_0010                                           = Form("%s/0010LHC11h_WeightingMethodA%s.dat",outputDir.Data(),meson.Data());
    TString fileNameOutputWeightingALHC11h_2050                                           = Form("%s/2050LHC11h_WeightingMethodA%s.dat",outputDir.Data(),meson.Data());

    cout << __LINE__ << endl;
    if(meson.CompareTo("Pi0")==0){
        // Declaration & calculation of combined spectrum

        graphCombInvYieldTotPbPb2760GeVA_0010       = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h_0010, sysErrorCollectionforRatiosLHC11h_0010,
                                                                                           xPtLimitsPi0, 23, offSetsPi0, offSetsPi0Sys,
                                                                                           graphCombInvYieldStatPbPb2760GeVA_0010, graphCombInvYieldSysPbPb2760GeVA_0010,
                                                                                           fileNameOutputWeightingALHC11h_0010);
        graphCombInvYieldTotPbPb2760GeVA_2050       = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h_2050, sysErrorCollectionforRatiosLHC11h_2050,
                                                                                           xPtLimitsPi0, 23, offSetsPi0, offSetsPi0Sys,
                                                                                           graphCombInvYieldStatPbPb2760GeVA_2050, graphCombInvYieldSysPbPb2760GeVA_2050,
                                                                                           fileNameOutputWeightingALHC11h_2050);

        while(graphCombInvYieldTotPbPb2760GeVA_0010->GetY()[0] < 1e-14){
          graphCombInvYieldTotPbPb2760GeVA_0010->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeVA_0010->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeVA_0010->RemovePoint(0);
          graphCombInvYieldTotPbPb2760GeVA_2050->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeVA_2050->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeVA_2050->RemovePoint(0);
        }

//         graphCombInvYieldSysPbPb2760GeVA_0010->Print();
//         graphCombInvYieldSysPbPb2760GeVA_2050->Print();
//         return;


    } else  if(meson.CompareTo("Eta")==0){
        cout << __LINE__ << endl;
        // Declaration & calculation of combined spectrum

        graphCombInvYieldTotPbPb2760GeVA_0010       = CombinePtPointsSpectraFullCorrMat(  statErrorCollectionLHC11h_0010,       sysErrorCollectionforRatiosLHC11h_0010,
                                                                                            xPtLimitsEta, 13, offSetsEta, offSetsEtaSys,
                                                                                            graphCombInvYieldStatPbPb2760GeVA_0010, graphCombInvYieldSysPbPb2760GeVA_0010,
                                                                                            fileNameOutputWeightingALHC11h_0010);
        graphCombInvYieldTotPbPb2760GeVA_2050       = CombinePtPointsSpectraFullCorrMat(  statErrorCollectionLHC11h_2050,       sysErrorCollectionforRatiosLHC11h_2050,
                                                                                            xPtLimitsEta, 13, offSetsEta, offSetsEtaSys,
                                                                                            graphCombInvYieldStatPbPb2760GeVA_2050, graphCombInvYieldSysPbPb2760GeVA_2050,
                                                                                            fileNameOutputWeightingALHC11h_2050);

        while(graphCombInvYieldTotPbPb2760GeVA_0010->GetY()[0]  < 1e-14){
          graphCombInvYieldTotPbPb2760GeVA_0010->RemovePoint(0);
          graphCombInvYieldTotPbPb2760GeVA_2050->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeVA_0010->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeVA_2050->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeVA_0010->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeVA_2050->RemovePoint(0);
        }

//         graphCombInvYieldTotPbPb2760GeVA_0010->Print();
//         graphCombInvYieldTotPbPb2760GeVA_2050->Print();
//         return;


    }

    // Reading weights from output file for plotting
    ifstream fileWeightsReadALHC11h_0010;
    fileWeightsReadALHC11h_0010.open(fileNameOutputWeightingALHC11h_0010,ios_base::in);
    ifstream fileWeightsReadALHC11h_2050;
    fileWeightsReadALHC11h_2050.open(fileNameOutputWeightingALHC11h_2050,ios_base::in);
    cout << "reading " << fileNameOutputWeightingALHC11h_0010 << " and " << fileNameOutputWeightingALHC11h_2050 << endl;
    Double_t xValuesReadALHC11h_0010[50];
    Double_t weightsReadALHC11h_0010[11][50];
    Int_t availableMeasALHC11h_0010[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    Double_t xValuesReadALHC11h_2050[50];
    Double_t weightsReadALHC11h_2050[11][50];
    Int_t availableMeasALHC11h_2050[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    Int_t nMeasSetALHC11h0010;
    Int_t nMeasSetALHC11h2050 = 2;
    if(meson.CompareTo("Eta")==0) nMeasSetALHC11h0010 = 2;
    else nMeasSetALHC11h0010 = 3;
    Int_t nPtBinsReadALHC11h              = 0;
    while(!fileWeightsReadALHC11h_0010.eof() && nPtBinsReadALHC11h < 50){
            TString garbage = "";
            if (nPtBinsReadALHC11h == 0){
                    fileWeightsReadALHC11h_0010 >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
                    fileWeightsReadALHC11h_2050 >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
                    for (Int_t i = 0; i < nMeasSetALHC11h0010; i++){
                      fileWeightsReadALHC11h_0010 >> availableMeasALHC11h_0010[i] ;
                    }
                    for (Int_t i = 0; i < nMeasSetALHC11h2050; i++){
                      fileWeightsReadALHC11h_2050 >> availableMeasALHC11h_2050[i] ;
                    }
                    cout << "read following measurements: ";
                    for (Int_t i = 0; i < 11; i++){
                            cout << availableMeasALHC11h_0010[i] << "\t" ;
                            cout << availableMeasALHC11h_2050[i] << "\t" ;
                    }
                    cout << endl;
            } else {
                    fileWeightsReadALHC11h_0010 >> xValuesReadALHC11h_0010[nPtBinsReadALHC11h-1];
                    fileWeightsReadALHC11h_2050 >> xValuesReadALHC11h_2050[nPtBinsReadALHC11h-1];
                    for (Int_t i = 0; i < nMeasSetALHC11h0010; i++){
                      fileWeightsReadALHC11h_0010 >> weightsReadALHC11h_0010[availableMeasALHC11h_0010[i]][nPtBinsReadALHC11h-1] ;
                    }
                    for (Int_t i = 0; i < nMeasSetALHC11h2050; i++){
                      fileWeightsReadALHC11h_2050 >> weightsReadALHC11h_2050[availableMeasALHC11h_2050[i]][nPtBinsReadALHC11h-1] ;
                    }
                    cout << "read: "<<  nPtBinsReadALHC11h << " xValuesReadALHC11h_0010 "<< xValuesReadALHC11h_0010[nPtBinsReadALHC11h-1] << "\t" ;
                    cout << "read: "<<  nPtBinsReadALHC11h << " xValuesReadALHC11h_2050 "<< xValuesReadALHC11h_2050[nPtBinsReadALHC11h-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetALHC11h0010; i++){
                      cout << weightsReadALHC11h_0010[availableMeasALHC11h_0010[i]][nPtBinsReadALHC11h-1] << "\t";
                    }
                    for (Int_t i = 0; i < nMeasSetALHC11h0010; i++){
                      cout << weightsReadALHC11h_2050[availableMeasALHC11h_2050[i]][nPtBinsReadALHC11h-1] << "\t";
                    }
                    cout << endl;
            }
            nPtBinsReadALHC11h++;
    }
    nPtBinsReadALHC11h = nPtBinsReadALHC11h-2 ;

    fileWeightsReadALHC11h_0010.close();
    fileWeightsReadALHC11h_2050.close();

    for (Int_t i = 0; i < nMeasSetALHC11h0010; i++){
            graphWeightsALHC11h_0010[availableMeasALHC11h_0010[i]] = new TGraph(nPtBinsReadALHC11h,xValuesReadALHC11h_0010,weightsReadALHC11h_0010[availableMeasALHC11h_0010[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsReadALHC11h; n++){
                    if (graphWeightsALHC11h_0010[availableMeasALHC11h_0010[i]]->GetY()[bin] == 0){
                            graphWeightsALHC11h_0010[availableMeasALHC11h_0010[i]]->RemovePoint(bin);
                    } else bin++;
            }
    }
    for (Int_t i = 0; i < nMeasSetALHC11h2050; i++){
            graphWeightsALHC11h_2050[availableMeasALHC11h_2050[i]] = new TGraph(nPtBinsReadALHC11h,xValuesReadALHC11h_2050,weightsReadALHC11h_2050[availableMeasALHC11h_2050[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsReadALHC11h; n++){
                    if (graphWeightsALHC11h_2050[availableMeasALHC11h_2050[i]]->GetY()[bin] == 0){
                            graphWeightsALHC11h_2050[availableMeasALHC11h_2050[i]]->RemovePoint(bin);
                    } else bin++;
            }
    }


    //**********************************************************************************************************************//
    //******************************************* Calculation of combined spectrum LHC11h **********************************//
    //**********************************************************************************************************************//

    for (Int_t i = 0; i< 11; i++){
        graphWeightsLHC11h_0010[i] = NULL;
        graphWeightsLHC11h_2050[i] = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNameOutputWeightingLHC11h_0010                  = Form("%s/0010LHC11h_WeightingComb%s.dat",outputDir.Data(),meson.Data());
    TString fileNameOutputWeightingLHC11h_2050                  = Form("%s/2050LHC11h_WeightingComb%s.dat",outputDir.Data(),meson.Data());

    cout << __LINE__ << endl;
    if(meson.CompareTo("Pi0")==0){

        graphCombInvYieldTotPbPb2760GeV_0010  = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h_0010,    sysErrorCollectionLHC11h_0010,
                                                                                        xPtLimitsPi0, 23,
                                                                                        offSetsPi0, offSetsPi0Sys,
                                                                                        graphCombInvYieldStatPbPb2760GeV_0010, graphCombInvYieldSysPbPb2760GeV_0010,
                                                                                        fileNameOutputWeightingLHC11h_0010);

        graphCombInvYieldTotPbPb2760GeV_2050 = CombinePtPointsSpectraFullCorrMat(  statErrorCollectionLHC11h_2050,    sysErrorCollectionLHC11h_2050,
                                                                                        xPtLimitsPi0, 23,
                                                                                        offSetsPi0, offSetsPi0Sys,
                                                                                        graphCombInvYieldStatPbPb2760GeV_2050, graphCombInvYieldSysPbPb2760GeV_2050,
                                                                                        fileNameOutputWeightingLHC11h_2050);
        while(graphCombInvYieldTotPbPb2760GeV_0010->GetY()[0]  < 1e-14){
          graphCombInvYieldTotPbPb2760GeV_0010->RemovePoint(0);
          graphCombInvYieldTotPbPb2760GeV_2050->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeV_0010->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeV_2050->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeV_0010->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeV_2050->RemovePoint(0);
        }

    } else  if(meson.CompareTo("Eta")==0){
        cout << __LINE__ << endl;

        graphCombInvYieldTotPbPb2760GeV_0010  = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h_0010,    sysErrorCollectionLHC11h_0010,
                                                                                        xPtLimitsEta, 13,
                                                                                        offSetsEta, offSetsEtaSys,
                                                                                        graphCombInvYieldStatPbPb2760GeV_0010, graphCombInvYieldSysPbPb2760GeV_0010,
                                                                                        fileNameOutputWeightingLHC11h_0010);

        graphCombInvYieldTotPbPb2760GeV_2050 = CombinePtPointsSpectraFullCorrMat(  statErrorCollectionLHC11h_2050,    sysErrorCollectionLHC11h_2050,
                                                                                        xPtLimitsEta, 13,
                                                                                        offSetsEta, offSetsEtaSys,
                                                                                        graphCombInvYieldStatPbPb2760GeV_2050, graphCombInvYieldSysPbPb2760GeV_2050,
                                                                                        fileNameOutputWeightingLHC11h_2050);

        while(graphCombInvYieldTotPbPb2760GeV_0010->GetY()[0]  < 1e-14){
          graphCombInvYieldTotPbPb2760GeV_0010->RemovePoint(0);
          graphCombInvYieldTotPbPb2760GeV_2050->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeV_0010->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeV_2050->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeV_0010->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeV_2050->RemovePoint(0);
        }

    }

    // Reading weights from output file for plotting
    ifstream fileWeightsReadLHC11h_0010;
    fileWeightsReadLHC11h_0010.open(fileNameOutputWeightingLHC11h_0010,ios_base::in);
    ifstream fileWeightsReadLHC11h_2050;
    fileWeightsReadLHC11h_2050.open(fileNameOutputWeightingLHC11h_2050,ios_base::in);

    cout << "reading " << fileNameOutputWeightingLHC11h_0010 << " and " << fileNameOutputWeightingLHC11h_2050 << endl;
    Double_t xValuesReadLHC11h_0010[50];
    Double_t xValuesReadLHC11h_2050[50];
    Double_t weightsReadLHC11h_0010[11][50];
    Double_t weightsReadLHC11h_2050[11][50];
    Int_t availableMeasLHC11h_0010[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    Int_t availableMeasLHC11h_2050[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    Int_t nMeasSetLHC11h0010;
    Int_t nMeasSetLHC11h2050 = 2;
    if(meson.CompareTo("Eta")==0) nMeasSetLHC11h0010 = 2;
    else nMeasSetLHC11h0010 = 3;
    Int_t nPtBinsReadLHC11h         = 0;
    while(!fileWeightsReadLHC11h_0010.eof() && nPtBinsReadLHC11h < 50){
        TString garbage = "";
        if (nPtBinsReadLHC11h == 0){
            fileWeightsReadLHC11h_0010 >> garbage;
            fileWeightsReadLHC11h_2050 >> garbage;
            for (Int_t i = 0; i < nMeasSetLHC11h0010; i++){
                fileWeightsReadLHC11h_0010 >> availableMeasLHC11h_0010[i];
            }
            for (Int_t i = 0; i < nMeasSetLHC11h2050; i++){
                fileWeightsReadLHC11h_2050 >> availableMeasLHC11h_2050[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableMeasLHC11h_0010[i] << "\t" ;
                cout << availableMeasLHC11h_2050[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsReadLHC11h_0010 >> xValuesReadLHC11h_0010[nPtBinsReadLHC11h-1];
            fileWeightsReadLHC11h_2050 >> xValuesReadLHC11h_2050[nPtBinsReadLHC11h-1];
            for (Int_t i = 0; i < nMeasSetLHC11h0010; i++){
                fileWeightsReadLHC11h_0010 >> weightsReadLHC11h_0010[availableMeasLHC11h_0010[i]][nPtBinsReadLHC11h-1];
            }
            for (Int_t i = 0; i < nMeasSetLHC11h2050; i++){
                fileWeightsReadLHC11h_2050 >> weightsReadLHC11h_2050[availableMeasLHC11h_2050[i]][nPtBinsReadLHC11h-1] ;
            }
            cout << "read: "<<  nPtBinsReadLHC11h << "\t"<< xValuesReadLHC11h_0010[nPtBinsReadLHC11h-1] << "\t" ;
            cout << "read: "<<  nPtBinsReadLHC11h << "\t"<< xValuesReadLHC11h_2050[nPtBinsReadLHC11h-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetLHC11h0010; i++){
                cout << weightsReadLHC11h_0010[availableMeasLHC11h_0010[i]][nPtBinsReadLHC11h-1] << "\t";
            }
            for (Int_t i = 0; i < nMeasSetLHC11h2050; i++){
                cout << weightsReadLHC11h_2050[availableMeasLHC11h_2050[i]][nPtBinsReadLHC11h-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsReadLHC11h++;
    }
    nPtBinsReadLHC11h = nPtBinsReadLHC11h-2 ;
    fileWeightsReadLHC11h_0010.close();
    fileWeightsReadLHC11h_2050.close();

    for (Int_t i = 0; i < nMeasSetLHC11h0010; i++){
        graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]] = new TGraph(nPtBinsReadLHC11h,xValuesReadLHC11h_0010,weightsReadLHC11h_0010[availableMeasLHC11h_0010[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsReadLHC11h; n++){
            if (graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]]->GetY()[bin] == 0 ){
                    graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]]->RemovePoint(bin);
            } else bin++;

        }
    }
    for (Int_t i = 0; i < nMeasSetLHC11h2050; i++){
        graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]] = new TGraph(nPtBinsReadLHC11h,xValuesReadLHC11h_2050,weightsReadLHC11h_2050[availableMeasLHC11h_2050[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsReadLHC11h; n++){
            if (graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]]->GetY()[bin] == 0){
                    graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]]->RemovePoint(bin);
            } else bin++;

        }
    }

    //**********************************************************************************************************************//
    //*********************************************** Plotting weights method **********************************************//
    //**********************************************************************************************************************//

    canvasWeights->cd();
    histo2DWeights->Draw("copy");

    TLegend* legendAccWeightsLHC11h = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*nMeasSetLHC11h0010*1.35), 32);
    for (Int_t i = 0; i < nMeasSetLHC11h0010; i++){
        DrawGammaSetMarkerTGraph(graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]],
                                markerStyleDet[availableMeasLHC11h_0010[i]],
                                markerSizeDet[availableMeasLHC11h_0010[i]]*0.5,
                                colorDet[availableMeasLHC11h_0010[i]] ,
                                colorDet[availableMeasLHC11h_0010[i]]);

        while(graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]]->GetX()[0] < 1.){
            graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]]->RemovePoint(0);
        }

        graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same");
        legendAccWeightsLHC11h->AddEntry(graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]],nameMeasGlobal[availableMeasLHC11h_0010[i]],"p");

    }
    legendAccWeightsLHC11h->Draw();
    labelWeightsEnergy->Draw();
    labelWeightsPi0->Draw();

    DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/%s_WeightsPCM-PHOS-EMCal_LHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    canvasWeights->cd();
    histo2DWeights->Draw("copy");

    for (Int_t i = 0; i < nMeasSetLHC11h2050; i++){

        DrawGammaSetMarkerTGraph(graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]],
                                markerStyleDet[availableMeasLHC11h_2050[i]],
                                markerSizeDet[availableMeasLHC11h_2050[i]]*0.5,
                                colorDet[availableMeasLHC11h_2050[i]],
                                colorDet[availableMeasLHC11h_2050[i]]);
        graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same");

        while(graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]]->GetX()[0] < 1.){
            graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]]->RemovePoint(0);
        }

    }
    legendAccWeightsLHC11h->Draw();
    labelWeightsEnergy->Draw();
    labelWeightsPi0->Draw();

    DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/%s_WeightsPCM-EMCal_LHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    //*********************************************************************************************************************//
    //************************************* Visualize relative systematic errors ******************************************//
    //*********************************************************************************************************************//
    //relative errors
    TGraphAsymmErrors* sysErrorRelCollectionLHC11h_0010[11];
    TGraphAsymmErrors* sysErrorRelCollectionLHC11h_2050[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionLHC11h_0010[i] = NULL;
        sysErrorRelCollectionLHC11h_2050[i] = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (sysErrorCollectionLHC11h_0010[i]) sysErrorRelCollectionLHC11h_0010[i] = CalculateRelErrUpAsymmGraph( sysErrorCollectionLHC11h_0010[i], Form("relativeSysErrorLHC11h_%s_0010", nameMeasGlobal[i].Data()));

        if (sysErrorCollectionLHC11h_2050[i]) sysErrorRelCollectionLHC11h_2050[i] = CalculateRelErrUpAsymmGraph( sysErrorCollectionLHC11h_2050[i], Form("relativeSysErrorLHC11h_%s_2050", nameMeasGlobal[i].Data()));
    }


    TCanvas* canvasRelSysErrLHC11h = new TCanvas("canvasRelSysErrLHC11h","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErrLHC11h, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErrLHC11h->SetLogx();

        TH2F * histo2DRelSysErrLHC11h = new TH2F("histo2DRelSysErrLHC11h","histo2DRelSysErrLHC11h",11000,0.23,70.,1000,0,80.5);
        SetStyleHistoTH2ForGraphs(histo2DRelSysErrLHC11h, "#it{p}_{T} (GeV/#it{c})","systematic error (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelSysErrLHC11h->GetXaxis()->SetMoreLogLabels();
        histo2DRelSysErrLHC11h->GetXaxis()->SetLabelOffset(-0.005);
        if(thesisPlotting || PaperPi0) histo2DRelSysErrLHC11h->GetYaxis()->SetRangeUser(0.,60.);
        if(thesisPlotting || PaperPi0) histo2DRelSysErrLHC11h->GetXaxis()->SetRangeUser(0.5,30.);
        histo2DRelSysErrLHC11h->Draw("copy");

        if(meson.CompareTo("Pi0")==0 && (thesisPlotting || PaperPi0)){
          while(sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[0]]->GetX()[0] < 1.){
            sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[0]]->RemovePoint(0);
          }
        }
        TLegend* legendRelSysErrLHC11h = GetAndSetLegend2(0.12, 0.92-(0.04*nMeasSetLHC11h0010*1.3), 0.3, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetLHC11h0010; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]], markerStyleDet[availableMeasLHC11h_0010[i]],
                                        markerSizeDet[availableMeasLHC11h_0010[i]]*0.5,
                                        colorDet[availableMeasLHC11h_0010[i]], colorDet[availableMeasLHC11h_0010[i]]);
            sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same");
            if(i == 0)
                legendRelSysErrLHC11h->AddEntry(sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]],Form("%s", /*(this thesis)",*/
(nameMeasGlobal[availableMeasLHC11h_0010[i]]).Data()),"p");
            else if(i == 2 || (meson.CompareTo("Eta")==0 && i == 1))
                legendRelSysErrLHC11h->AddEntry(sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]],Form("%s", /*(*)",*/(nameMeasGlobal[availableMeasLHC11h_0010[i]]).Data()),"p");
            else legendRelSysErrLHC11h->AddEntry(sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]],nameMeasGlobal[availableMeasLHC11h_0010[i]],"p");
        }
        legendRelSysErrLHC11h->Draw();

        TLatex *labelRelSysErrEnergy0010 = new TLatex(0.6,0.89,collisionSystemPbPb0010.Data());
        SetStyleTLatex( labelRelSysErrEnergy0010, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrEnergy0010->SetTextFont(43);
        labelRelSysErrEnergy0010->Draw();

        TLatex *labelthisthesisErr = new TLatex(0.28,0.88,"(this thesis)");
        SetStyleTLatex( labelthisthesisErr, 0.85*textSizeLabelsPixel,4);
        labelthisthesisErr->SetTextFont(43);
//         if(thesisPlotting) labelthisthesisErr->Draw();

        if(meson.CompareTo("Pi0")==0){
            labelRelSysErrPi0= new TLatex(0.6,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        } else if(meson.CompareTo("Eta")==0){
            labelRelSysErrPi0= new TLatex(0.6,0.85,"#eta #rightarrow #gamma#gamma");
        }
        SetStyleTLatex( labelRelSysErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrPi0->SetTextFont(43);
        labelRelSysErrPi0->Draw();

        TLatex *textEMCalNote = new TLatex(0.15,0.73,"(*) ALICE preliminary input");
        SetStyleTLatex( textEMCalNote,  0.8*textSizeLabelsPixel,4);
        textEMCalNote->SetTextFont(43);
        TLatex *textEMCalNote1 = new TLatex(0.17,0.69,"arXiv:1609.06106"); //NPA (2016) 956
        SetStyleTLatex( textEMCalNote1,  0.8*textSizeLabelsPixel,4);
        textEMCalNote1->SetTextFont(43);
        TLatex *textEMCalNote2 = new TLatex(0.17,0.655,"analysis by A. Morreale");
        SetStyleTLatex( textEMCalNote2,  0.8*textSizeLabelsPixel,4);
        textEMCalNote2->SetTextFont(43);
//         textEMCalNote->Draw();
//         textEMCalNote1->Draw();
//         textEMCalNote2->Draw();

    canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrLHC11h_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrLHC11h_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


        canvasRelSysErrLHC11h->cd();
        histo2DRelSysErrLHC11h->Draw("copy");

        if(meson.CompareTo("Pi0")==0 && (thesisPlotting || PaperPi0)){
          while(sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[0]]->GetX()[0] < 1.){
            sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[0]]->RemovePoint(0);
          }
        }

        TLegend* legendRelSysErrLHC11h2050 = GetAndSetLegend2(0.12, 0.92-(0.04*nMeasSetLHC11h2050*1.35), 0.5, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetLHC11h2050; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]], markerStyleDet[availableMeasLHC11h_2050[i]],
                                        markerSizeDet[availableMeasLHC11h_2050[i]]*0.5,
                                        colorDet[availableMeasLHC11h_2050[i]], colorDet[availableMeasLHC11h_2050[i]]);
            sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same");
            if(i == 0)
                legendRelSysErrLHC11h2050->AddEntry(sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]],Form("%s",/* (this thesis)",*/(nameMeasGlobal[availableMeasLHC11h_2050[i]]).Data()),"p");
            else if(i == 1)
                legendRelSysErrLHC11h2050->AddEntry(sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]],Form("%s", /*(*)",*/(nameMeasGlobal[availableMeasLHC11h_2050[i]]).Data()),"p");
            else legendRelSysErrLHC11h2050->AddEntry(sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]],nameMeasGlobal[availableMeasLHC11h_2050[i]],"p");
        }
        TLatex *labelRelSysErrEnergy2050 = new TLatex(0.6,0.89,collisionSystemPbPb2050.Data());
        SetStyleTLatex( labelRelSysErrEnergy2050, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrEnergy2050->SetTextFont(43);
        labelRelSysErrEnergy2050->Draw();
//         if(thesisPlotting) labelthisthesisErr->Draw();

        legendRelSysErrLHC11h2050->Draw();
        labelRelSysErrPi0->Draw();
//         textEMCalNote->Draw();
//         textEMCalNote1->Draw();
//         textEMCalNote2->Draw();

    canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrLHC11h_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrLHC11h_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasRelSysErrLHC11h->cd();
        histo2DRelSysErrLHC11h->GetYaxis()->SetRangeUser(0,30.5);
        histo2DRelSysErrLHC11h->Draw("copy");

        for (Int_t i = 0; i < nMeasSetLHC11h0010; i++){
            sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same");
        }

        legendRelSysErrLHC11h->Draw();
        labelRelSysErrEnergy0010->Draw();
        labelRelSysErrPi0->Draw();

//     canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrZoomedLHC11h_0010.%s",outputDir.Data(),meson.Data(), suffix.Data()));
    canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrZoomedLHC11h_0010.%s",PubNotePlots.Data(),meson.Data(), suffix.Data()));

    canvasRelSysErrLHC11h->cd();
        histo2DRelSysErrLHC11h->GetYaxis()->SetRangeUser(0,30.5);
        histo2DRelSysErrLHC11h->Draw("copy");

        for (Int_t i = 0; i < nMeasSetLHC11h2050; i++){
            sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same");

        }
        legendRelSysErrLHC11h->Draw();
        labelRelSysErrEnergy2050->Draw();
        labelRelSysErrPi0->Draw();

//     canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrZoomedLHC11h_2050.%s",outputDir.Data(),meson.Data(), suffix.Data()));
    canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrZoomedLHC11h_2050.%s",PubNotePlots.Data(),meson.Data(), suffix.Data()));

      histo2DRelSysErrLHC11h->GetYaxis()->SetRangeUser(0.,60.5);
      histo2DRelSysErrLHC11h->GetXaxis()->SetRangeUser(0.3,25.);
      if(meson.CompareTo("Eta")==0){
        histo2DRelSysErrLHC11h->GetYaxis()->SetRangeUser(0.,60.5);
        histo2DRelSysErrLHC11h->GetXaxis()->SetRangeUser(0.8,15.);
      }
      histo2DRelSysErrLHC11h->Draw("copy");

          DrawGammaSetMarkerTGraph(sysErrorRelCollectionLHC11h_0010[0], markerStyleDet[0], markerSizeDet[0]*0.5, colorCombo0010 ,colorCombo0010);
          DrawGammaSetMarkerTGraph(sysErrorRelCollectionLHC11h_2050[0], markerStyleDet[0], markerSizeDet[0]*0.5, colorCombo2050, colorCombo2050);

          sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[0]]->Draw("p,same");
          sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[0]]->Draw("p,same");

          TLegend* legendRelSysErrLHC11hforAN = GetAndSetLegend2(0.12, 0.92-(0.06*2*1.35), 0.45, 0.92, 32);
          legendRelSysErrLHC11hforAN->AddEntry(sysErrorRelCollectionLHC11h_0010[0],"PCM -   0#font[122]{-}10%","p");
          legendRelSysErrLHC11hforAN->AddEntry(sysErrorRelCollectionLHC11h_2050[0],"PCM - 20#font[122]{-}50%","p");
          legendRelSysErrLHC11hforAN->Draw();

          TLatex *labelRelSysErrEnergy = new TLatex(0.7,0.89,collisionSystem2760GeV.Data());
          SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel,4);
          labelRelSysErrEnergy->SetTextFont(43);
          labelRelSysErrEnergy->Draw();
          labelRelSysErrPi0->Draw();

    canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrZoomedLHC11hforAN.%s",outputDir.Data(),meson.Data(),suffix.Data()));


//     if(meson.CompareTo("Pi0")==0){
//       for(Int_t i = 0; i< 2; i++){
//         sysErrorRel2010and2011_0010[i] = NULL;
//         sysErrorRel2010and2011_2040[i] = NULL;
//         sysErrorRelRAA2010and2011_0010[i] = NULL;
//         sysErrorRelRAA2010and2011_2040[i] = NULL;
//       }
//
//       sysErrorRel2010and2011_0010[0] = CalculateRelErrUpAsymmGraph( graphPCMPi0InvYieldSysPbPb2760GeV_0010,"relativeSysError2011_0010");
// //       sysErrorRelRAA2010and2011_0010[0] = CalculateRelErrUpAsymmGraph( graphPCMPi0RAASysPbPb2760GeV_0010,"relativeSysError2011_0010");
//       sysErrorRel2010and2011_2040[0] = CalculateRelErrUpAsymmGraph( graphPCMPi0InvYieldSysPbPb2760GeV_2040,"relativeSysError2011_2040");
// //       sysErrorRelRAA2010and2011_2040[0] = CalculateRelErrUpAsymmGraph( graphPCMPi0RAASysPbPb2760GeV_2040,"relativeSysError2011_2040");
//
//       sysErrorRel2010and2011_0010[1] = CalculateRelErrUpAsymmGraph( graphPCMPubPi0InvYieldSysPbPb2760GeV_0010,"relativeSysError2010_0010");
// //       sysErrorRelRAA2010and2011_0010[1] = CalculateRelErrUpAsymmGraph( graphPCMPi0RAASys2010PbPb2760GeV_0010,"relativeSysError2010_0010");
//       sysErrorRel2010and2011_2040[1] = CalculateRelErrUpAsymmGraph( graphPCMPubPi0InvYieldSysPbPb2760GeV_2040,"relativeSysError2010_2040");
// //       sysErrorRelRAA2010and2011_2040[1] = CalculateRelErrUpAsymmGraph( graphPCMPi0RAASys2010PbPb2760GeV_2040,"relativeSysError2010_2040");
//
//       for(Int_t i = 0; i< 2; i++){
//         sysErrorRel2010and2011_0010[i]->Print();
//         sysErrorRel2010and2011_2040[i]->Print();
// //         sysErrorRelRAA2010and2011_0010[i]->Print();
// //         sysErrorRelRAA2010and2011_2040[i]->Print();
//       }
//
//       canvasRelSysErrLHC11h->cd();
//       histo2DRelSysErrLHC11h->Draw("copy");
//
//           DrawGammaSetMarkerTGraph(sysErrorRel2010and2011_0010[0], 20, 2, kRed, kRed);
//           sysErrorRel2010and2011_0010[0]->Draw("p,same");
//           DrawGammaSetMarkerTGraph(sysErrorRel2010and2011_0010[1], 24, 2, colorCombo0010 , colorCombo0010);
//           sysErrorRel2010and2011_0010[1]->Draw("p,same");
//           DrawGammaSetMarkerTGraph(sysErrorRel2010and2011_2040[0], 21, 2, kAzure-2, kAzure-2);
//           sysErrorRel2010and2011_2040[0]->Draw("p,same");
//           DrawGammaSetMarkerTGraph(sysErrorRel2010and2011_2040[1], 25, 2, colorEMCal2050 , colorEMCal2050);
//           sysErrorRel2010and2011_2040[1]->Draw("p,same");
//
//           TLegend* legendRelSysErr20102011 = GetAndSetLegend2(0.12, 0.96-(0.06*2*1.35), 0.45, 0.96, 32);
//           legendRelSysErr20102011->AddEntry(sysErrorRel2010and2011_0010[0],"PCM 2011 - 0#font[122]{-}10%","p");
//           legendRelSysErr20102011->AddEntry(sysErrorRel2010and2011_0010[1],"PCM 2010 - 0#font[122]{-}10%","p");
//           legendRelSysErr20102011->AddEntry(sysErrorRel2010and2011_2040[0],"PCM 2011 - 20-40%","p");
//           legendRelSysErr20102011->AddEntry(sysErrorRel2010and2011_2040[1],"PCM 2010 - 20-40%","p");
//           legendRelSysErr20102011->Draw();
//           labelRelSysErrEnergy->Draw();
//           labelRelSysErrPi0->Draw();
//
//       canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErr2010and2011.%s",outputDir.Data(),meson.Data(),suffix.Data()));

//         TH2F * histo2DRelSysErrAllRaa;
//         histo2DRelSysErrAllRaa = new TH2F("histo2DRelSysErrAllRaa","histo2DRelSysErrAllRaa",11000,0.23,70.,1000,0,80.5);
//         SetStyleHistoTH2ForGraphs(histo2DRelSysErrAllRaa, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA} sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
//         histo2DRelSysErrAllRaa->GetXaxis()->SetMoreLogLabels();
//         histo2DRelSysErrAllRaa->GetXaxis()->SetLabelOffset(-0.01);
//         histo2DRelSysErrAllRaa->GetYaxis()->SetRangeUser(0,80.5);
//         histo2DRelSysErrAllRaa->GetXaxis()->SetRangeUser(0.2,20);
//         histo2DRelSysErrAllRaa->Draw("copy");
//
//            DrawGammaSetMarkerTGraph(sysErrorRelRAA2010and2011_0010[0], 20, 2, kRed, kRed);
//            sysErrorRelRAA2010and2011_0010[0]->Draw("p,same");
//            DrawGammaSetMarkerTGraph(sysErrorRelRAA2010and2011_0010[1], 24, 2, colorCombo0010 , colorCombo0010);
//            sysErrorRelRAA2010and2011_0010[1]->Draw("p,same");
//            DrawGammaSetMarkerTGraph(sysErrorRelRAA2010and2011_2040[0], 21, 2, kAzure-2, kAzure-2);
//            sysErrorRelRAA2010and2011_2040[0]->Draw("p,same");
//            DrawGammaSetMarkerTGraph(sysErrorRelRAA2010and2011_2040[1], 25, 2, colorEMCal2050 , colorEMCal2050);
//            sysErrorRelRAA2010and2011_2040[1]->Draw("p,same");
//
//            legendRelSysErr20102011->Draw();
//            labelRelSysErrEnergy->Draw();
//            labelRelSysErrPi0->Draw();
//
//       canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RaaRelSysErr2010and2011.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     }


    //*********************************************************************************************************************//
    //************************************* Visualize relative statistical errors *****************************************//
    //*********************************************************************************************************************//

    TH1D* statErrorRelCollectionLHC11h_0010[11];
    TH1D* statErrorRelCollectionLHC11h_2050[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionLHC11h_0010[i] = NULL;
        statErrorRelCollectionLHC11h_2050[i] = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (statErrorCollectionLHC11h_0010[i]) statErrorRelCollectionLHC11h_0010[i] = CalculateRelErrUpTH1D( statErrorCollectionLHC11h_0010[i], Form("relativeStatErrorLHC11h_%s_0010", nameMeasGlobal[i].Data()));

        if (statErrorCollectionLHC11h_2050[i]) statErrorRelCollectionLHC11h_2050[i] = CalculateRelErrUpTH1D( statErrorCollectionLHC11h_2050[i], Form("relativeStatErrorLHC11h_%s_2050", nameMeasGlobal[i].Data()));
    }


    TCanvas* canvasRelStatErrLHC11h = new TCanvas("canvasRelStatErrLHC11h","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelStatErrLHC11h, 0.08, 0.02, 0.035, 0.09);
        canvasRelStatErrLHC11h->SetLogx();

        TH2F * histo2DRelStatErrLHC11h = new TH2F("histo2DRelStatErrLHC11h","histo2DRelStatErrLHC11h",11000,0.23,70.,1000,0,80.5);
        SetStyleHistoTH2ForGraphs(histo2DRelStatErrLHC11h, "#it{p}_{T} (GeV/#it{c})","statistical error (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelStatErrLHC11h->GetXaxis()->SetMoreLogLabels();
        histo2DRelStatErrLHC11h->GetXaxis()->SetLabelOffset(-0.005);
        if(thesisPlotting || PaperPi0) histo2DRelStatErrLHC11h->GetYaxis()->SetRangeUser(0.,60.);
        if(thesisPlotting || PaperPi0) histo2DRelStatErrLHC11h->GetXaxis()->SetRangeUser(0.5,30.);
        histo2DRelStatErrLHC11h->Draw("copy");

        if(meson.CompareTo("Pi0")==0 && thesisPlotting){
            statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[0]]->GetXaxis()->SetRangeUser(1,14.);
        } else if(meson.CompareTo("Eta")==0){
            statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[0]]->GetXaxis()->SetRangeUser(1,10.);
        }

        TLegend* legendRelStatErrLHC11h = GetAndSetLegend2(0.12, 0.92-(0.04*nMeasSetLHC11h0010*1.3), 0.3, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetLHC11h0010; i++){
            DrawGammaSetMarker(statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]],
                            markerStyleDet[availableMeasLHC11h_0010[i]], markerSizeDet[availableMeasLHC11h_0010[i]]*0.5,
                            colorDet[availableMeasLHC11h_0010[i]] , colorDet[availableMeasLHC11h_0010[i]]);
            statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same");
            if(i == 0)
                legendRelStatErrLHC11h->AddEntry(statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]],Form("%s",/* (this thesis)",*/(nameMeasGlobal[availableMeasLHC11h_0010[i]]).Data()),"p");
            else if(i == 2 || (meson.CompareTo("Eta")==0 && i == 1))
                legendRelStatErrLHC11h->AddEntry(statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]],Form("%s",/* (*)",*/(nameMeasGlobal[availableMeasLHC11h_0010[i]]).Data()),"p");
            else legendRelStatErrLHC11h->AddEntry(statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]],nameMeasGlobal[availableMeasLHC11h_0010[i]],"p");
        }
        legendRelStatErrLHC11h->Draw();

        TLatex *labelRelStatErrEnergy0010 = new TLatex(0.6,0.89,collisionSystemPbPb0010.Data());
        SetStyleTLatex( labelRelStatErrEnergy0010, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrEnergy0010->SetTextFont(43);
        labelRelStatErrEnergy0010->Draw();
        if(meson.CompareTo("Pi0")==0){
            labelRelStatErrPi0= new TLatex(0.6,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        } else if(meson.CompareTo("Eta")==0){
            labelRelStatErrPi0= new TLatex(0.6,0.85,"#eta #rightarrow #gamma#gamma");
        }
        SetStyleTLatex( labelRelStatErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrPi0->SetTextFont(43);
        labelRelStatErrPi0->Draw();
//         if(thesisPlotting) labelthisthesisErr->Draw();

//         TLatex *textEMCalNote = new TLatex(0.17,0.75,"(*) ALICE preliminary input");
//         SetStyleTLatex( textEMCalNote,  0.8*textSizeLabelsPixel,4);
//         textEMCalNote->SetTextFont(43);
//         TLatex *textEMCalNote1 = new TLatex(0.17,0.7,"arXiv:1609.06106");
//         SetStyleTLatex( textEMCalNote1,  0.8*textSizeLabelsPixel,4);
//         textEMCalNote1->SetTextFont(43);
//         TLatex *textEMCalNote2 = new TLatex(0.17,0.65,"analysis by A. Morreale");
//         SetStyleTLatex( textEMCalNote2,  0.8*textSizeLabelsPixel,4);
//         textEMCalNote2->SetTextFont(43);

//         textEMCalNote->Draw();
//         textEMCalNote1->Draw();
//         textEMCalNote2->Draw();

    canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrLHC11h_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrLHC11h_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasRelStatErrLHC11h->cd();
        histo2DRelStatErrLHC11h->Draw("copy");

        if(meson.CompareTo("Pi0")==0  && (thesisPlotting || PaperPi0)){
            statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_0010[0]]->GetXaxis()->SetRangeUser(1,14.);
        } else if(meson.CompareTo("Eta")==0){
            statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_0010[0]]->GetXaxis()->SetRangeUser(1,10.);
        }

        TLegend* legendRelStatErrLHC11h2050 = GetAndSetLegend2(0.12, 0.92-(0.04*nMeasSetLHC11h2050*1.35), 0.5, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetLHC11h2050; i++){
            DrawGammaSetMarker(statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]],
                            markerStyleDet[availableMeasLHC11h_2050[i]], markerSizeDet[availableMeasLHC11h_2050[i]]*0.5,
                            colorDet[availableMeasLHC11h_2050[i]], colorDet[availableMeasLHC11h_2050[i]]);
            statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same");
            if(i == 0)
                legendRelStatErrLHC11h2050->AddEntry(statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]],Form("%s",/* (this thesis)",*/(nameMeasGlobal[availableMeasLHC11h_2050[i]]).Data()),"p");
            else if(i == 2 || (meson.CompareTo("Eta")==0 && i == 1))
                legendRelStatErrLHC11h2050->AddEntry(statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]],Form("%s",/* (*)",*/(nameMeasGlobal[availableMeasLHC11h_2050[i]]).Data()),"p");
            else legendRelStatErrLHC11h2050->AddEntry(statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]],nameMeasGlobal[availableMeasLHC11h_2050[i]],"p");
        }
        legendRelStatErrLHC11h2050->Draw();
        TLatex *labelRelStatErrEnergy2050 = new TLatex(0.6,0.89,collisionSystemPbPb2050.Data());
        SetStyleTLatex( labelRelStatErrEnergy2050, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrEnergy2050->SetTextFont(43);
        labelRelStatErrEnergy2050->Draw();

        labelRelStatErrPi0->Draw();
//         textEMCalNote->Draw();
//         textEMCalNote1->Draw();
//         textEMCalNote2->Draw();
//         if(thesisPlotting) labelthisthesisErr->Draw();

    canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrLHC11h_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrLHC11h_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    canvasRelStatErrLHC11h->cd();
        histo2DRelStatErrLHC11h->GetYaxis()->SetRangeUser(0,30.5);
        histo2DRelStatErrLHC11h->Draw("copy");
        for (Int_t i = 0; i < nMeasSetLHC11h0010; i++){
            statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same");
        }
        legendRelStatErrLHC11h->Draw();
        labelRelStatErrEnergy0010->Draw();
        labelRelStatErrPi0->Draw();

//     canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrZoomedLHC11h_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrZoomedLHC11h_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasRelStatErrLHC11h->cd();
        histo2DRelStatErrLHC11h->GetYaxis()->SetRangeUser(0,30.5);
        histo2DRelStatErrLHC11h->Draw("copy");
        for (Int_t i = 0; i < nMeasSetLHC11h2050; i++){
                statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same");
        }
        legendRelStatErrLHC11h->Draw();
        labelRelStatErrEnergy2050->Draw();
        labelRelStatErrPi0->Draw();

//     canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrZoomedLHC11h_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrZoomedLHC11h_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

        histo2DRelStatErrLHC11h->GetYaxis()->SetRangeUser(0.,60.5);
        histo2DRelStatErrLHC11h->GetXaxis()->SetRangeUser(0.3,25.);
        if(meson.CompareTo("Eta")==0){
          histo2DRelStatErrLHC11h->GetYaxis()->SetRangeUser(0.,60.5);
          histo2DRelStatErrLHC11h->GetXaxis()->SetRangeUser(0.8,15.);
        }
        histo2DRelStatErrLHC11h->Draw("copy");

          DrawGammaSetMarker(statErrorRelCollectionLHC11h_0010[0], markerStyleDet[0], markerSizeDet[0]*0.5, colorCombo0010 ,colorCombo0010);
          DrawGammaSetMarker(statErrorRelCollectionLHC11h_2050[0], markerStyleDet[0], markerSizeDet[0]*0.5, colorCombo2050, colorCombo2050);

          statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[0]]->Draw("p,same");
          statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[0]]->Draw("p,same");

          TLegend* legendRelStatErrLHC11hforAN = GetAndSetLegend2(0.14, 0.92-(0.06*2*1.35), 0.45, 0.92, 32);
          legendRelStatErrLHC11hforAN->AddEntry(statErrorRelCollectionLHC11h_0010[0],"PCM -   0#font[122]{-}10%","p");
          legendRelStatErrLHC11hforAN->AddEntry(statErrorRelCollectionLHC11h_2050[0],"PCM - 20#font[122]{-}50%","p");
          legendRelStatErrLHC11hforAN->Draw();
          TLatex *labelRelStatErrEnergy = new TLatex(0.7,0.89,collisionSystem2760GeV.Data());
          SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel,4);
          labelRelStatErrEnergy->SetTextFont(43);
          labelRelStatErrEnergy->Draw();
          labelRelStatErrPi0->Draw();

    canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrZoomedLHC11hforAN.%s",outputDir.Data(),meson.Data(),suffix.Data()));



    if(meson.CompareTo("Pi0")==0){
      for(Int_t i = 0; i< 2; i++){
        statErrorRel2010and2011_0010[i] = NULL;
        statErrorRel2010and2011_2040[i] = NULL;
        statErrorRelRAA2010and2011_0010[i] = NULL;
        statErrorRelRAA2010and2011_2040[i] = NULL;
      }

      statErrorRel2010and2011_0010[0] = CalculateRelErrUpAsymmGraph( graphPCMPi0InvYieldStatPbPb2760GeV_0010,"relativeStatError2011_0010");
//       statErrorRelRAA2010and2011_0010[0] = CalculateRelErrUpAsymmGraph( graphPCMPi0RAAStatPbPb2760GeV_0010,"relativeStatError2011_0010");
      statErrorRel2010and2011_2040[0] = CalculateRelErrUpAsymmGraph( graphPCMPi0InvYieldStatPbPb2760GeV_2040,"relativeStatError2011_2040");
//       statErrorRelRAA2010and2011_2040[0] = CalculateRelErrUpAsymmGraph( graphPCMPi0RAAStatPbPb2760GeV_2040,"relativeStatError2011_2040");

      statErrorRel2010and2011_0010[1] = CalculateRelErrUpAsymmGraph( graphPCMPubPi0InvYieldStatPbPb2760GeV_0010,"relativeStatError2010_0010");
//       statErrorRelRAA2010and2011_0010[1] = CalculateRelErrUpAsymmGraph( graphPCMPi0RAAStat2010PbPb2760GeV_0010,"relativeStatError2010_0010");
      statErrorRel2010and2011_2040[1] = CalculateRelErrUpAsymmGraph( graphPCMPubPi0InvYieldStatPbPb2760GeV_2040,"relativeStatError2010_2040");
//       statErrorRelRAA2010and2011_2040[1] = CalculateRelErrUpAsymmGraph( graphPCMPi0RAAStat2010PbPb2760GeV_2040,"relativeStatError2010_2040");

//       for(Int_t i = 0; i< 2; i++){
//         statErrorRel2010and2011_0010[i]->Print();
//         statErrorRel2010and2011_2040[i]->Print();
//         statErrorRelRAA2010and2011_0010[i]->Print();
//         statErrorRelRAA2010and2011_2040[i]->Print();
//       }

      canvasRelStatErrLHC11h->cd();

      if(thesisPlotting || PaperPi0){
          while(statErrorRel2010and2011_0010[0]->GetX()[0] < 1.) statErrorRel2010and2011_0010[0]->RemovePoint(0);
          while(statErrorRel2010and2011_2040[0]->GetX()[0] < 1.) statErrorRel2010and2011_2040[0]->RemovePoint(0);
          while(statErrorRel2010and2011_0010[1]->GetX()[0] < 1.) statErrorRel2010and2011_0010[1]->RemovePoint(0);
          while(statErrorRel2010and2011_2040[1]->GetX()[0] < 1.) statErrorRel2010and2011_2040[1]->RemovePoint(0);

        histo2DRelStatErrLHC11h->GetYaxis()->SetRangeUser(0,70.5);
        histo2DRelStatErrLHC11h->GetXaxis()->SetRangeUser(0.7,20);
      }
      histo2DRelStatErrLHC11h->Draw("copy");

          DrawGammaSetMarkerTGraph(statErrorRel2010and2011_0010[0], 20, 2, colorCombo0010, colorCombo0010);
          statErrorRel2010and2011_0010[0]->Draw("p,same");
          DrawGammaSetMarkerTGraph(statErrorRel2010and2011_2040[0], 21, 2, colorCombo2050, colorCombo2050);
          statErrorRel2010and2011_2040[0]->Draw("p,same");

          DrawGammaSetMarkerTGraph(statErrorRel2010and2011_0010[1], 24, 2, colorCombo0010+1, colorCombo0010+1);
          statErrorRel2010and2011_0010[1]->Draw("p,same");
          DrawGammaSetMarkerTGraph(statErrorRel2010and2011_2040[1], 25, 2, kBlue+1 ,  kBlue+1);
          statErrorRel2010and2011_2040[1]->Draw("p,same");

          TLegend* legendlabel = GetAndSetLegend2(0.11, 0.86-(0.037*2), 0.5, 0.86, 32);
          legendlabel->SetNColumns(2);
          legendlabel->SetMargin(0.1);
//           legendlabel->SetHeader(collisionSystem2760GeV.Data());
          legendlabel->AddEntry((TObject*)0,"#pi^{0} 2010","");
          legendlabel->AddEntry((TObject*)0,"#pi^{0} 2011","");
          legendlabel->AddEntry((TObject*)0,"EPJC 74 (2014) 3108","");
          legendlabel->AddEntry((TObject*)0,"This thesis","");
          legendlabel->Draw();

          TLegend* legendRelStatErr20102011 = GetAndSetLegend2(0.12, 0.92-(0.06*3*1.35), 0.6, 0.92, 32);
          legendRelStatErr20102011->SetNColumns(2);
          legendRelStatErr20102011->SetHeader(collisionSystem2760GeV.Data());
          legendRelStatErr20102011->AddEntry((TObject*)0,"","");
//           legendRelStatErr20102011->AddEntry((TObject*)0,"","");
//           legendRelStatErr20102011->AddEntry((TObject*)0,"","");
          legendRelStatErr20102011->AddEntry((TObject*)0,"","");
          legendRelStatErr20102011->AddEntry(statErrorRel2010and2011_0010[1],"0#font[122]{-}10%","p");
          legendRelStatErr20102011->AddEntry(statErrorRel2010and2011_0010[0],"0#font[122]{-}10%","p");
          legendRelStatErr20102011->AddEntry(statErrorRel2010and2011_2040[1],"20#font[122]{-}40%","p");
          legendRelStatErr20102011->AddEntry(statErrorRel2010and2011_2040[0],"20#font[122]{-}40%","p");
          legendRelStatErr20102011->Draw();
//           labelRelStatErrEnergy->Draw();
//           labelRelStatErrPi0->Draw();
//            thesisLabelLowRight->Draw();

      canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErr2010and2011.%s",outputDir.Data(),meson.Data(),suffix.Data()));

//         TH2F * histo2DRelStatErrAllRaa;
//         histo2DRelStatErrAllRaa = new TH2F("histo2DRelStatErrAllRaa","histo2DRelStatErrAllRaa",11000,0.23,70.,1000,0,80.5);
//         SetStyleHistoTH2ForGraphs(histo2DRelStatErrAllRaa, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA} stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
//         histo2DRelStatErrAllRaa->GetXaxis()->SetMoreLogLabels();
//         histo2DRelStatErrAllRaa->GetXaxis()->SetLabelOffset(-0.01);
//         histo2DRelStatErrAllRaa->GetYaxis()->SetRangeUser(0,80.5);
//         histo2DRelStatErrAllRaa->GetXaxis()->SetRangeUser(0.2,20);
//         histo2DRelStatErrAllRaa->Draw("copy");
//
//            DrawGammaSetMarkerTGraph(statErrorRelRAA2010and2011_0010[0], 20, 2, kRed, kRed);
//            statErrorRelRAA2010and2011_0010[0]->Draw("p,same");
//            DrawGammaSetMarkerTGraph(statErrorRelRAA2010and2011_0010[1], 24, 2, colorCombo0010 , colorCombo0010);
//            statErrorRelRAA2010and2011_0010[1]->Draw("p,same");
//            DrawGammaSetMarkerTGraph(statErrorRelRAA2010and2011_2040[0], 21, 2, kAzure-2, kAzure-2);
//            statErrorRelRAA2010and2011_2040[0]->Draw("p,same");
//            DrawGammaSetMarkerTGraph(statErrorRelRAA2010and2011_2040[1], 25, 2, colorEMCal2050 , colorEMCal2050);
//            statErrorRelRAA2010and2011_2040[1]->Draw("p,same");
//
//            legendRelStatErr20102011->Draw();
//            labelRelStatErrEnergy->Draw();
//            labelRelStatErrPi0->Draw();
//
//       canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RaaRelStatErr2010and2011.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    }



    //  *********************************************************************************************************************
    //  ************************ Visualize relative total errors of different combination methods Pi0 ***********************
    //  *********************************************************************************************************************
    TGraphAsymmErrors* graphCombInvYieldsRelStat_0010  = CalculateRelErrUpAsymmGraph( graphCombInvYieldStatPbPb2760GeV_0010, "relativeStatError_0010");
    TGraphAsymmErrors* graphCombInvYieldsRelSys_0010   = CalculateRelErrUpAsymmGraph( graphCombInvYieldSysPbPb2760GeV_0010, "relativeSysError_0010");
    TGraphAsymmErrors* graphCombInvYieldsRelTot_0010 = CalculateRelErrUpAsymmGraph( graphCombInvYieldTotPbPb2760GeV_0010, "relativeTotalError_0010");

    TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();

    TH2F * histo2DRelTotErr = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,0.3,50.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelTotErr->GetYaxis()->SetRangeUser(0,50.5);
    histo2DRelTotErr->Draw("copy");

        if(meson.CompareTo("Pi0")==0 && (thesisPlotting || PaperPi0)){
          while(graphCombInvYieldsRelTot_0010->GetX()[0] < /* 0.8 */ 1.){
            graphCombInvYieldsRelTot_0010->RemovePoint(0);
            graphCombInvYieldsRelStat_0010->RemovePoint(0);
            graphCombInvYieldsRelSys_0010->RemovePoint(0);
          }
        }

        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldsRelTot_0010, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombInvYieldsRelTot_0010->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldsRelStat_0010, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombInvYieldsRelStat_0010->Draw("l,x0,same");
        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldsRelSys_0010, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombInvYieldsRelSys_0010->SetLineStyle(7);
        graphCombInvYieldsRelSys_0010->Draw("l,x0,same");

        TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3/*4*/), 0.45, 0.92, 32);
//         legendRelTotErr3->SetHeader("This thesis");
        legendRelTotErr3->AddEntry(graphCombInvYieldsRelTot_0010,"tot","p");
        legendRelTotErr3->AddEntry(graphCombInvYieldsRelStat_0010,"stat","l");
        legendRelTotErr3->AddEntry(graphCombInvYieldsRelSys_0010,"sys","l");
        legendRelTotErr3->Draw();
        labelRelStatErrEnergy0010->Draw();
        labelRelStatErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/%s_RelErrdecomp_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRelTotErr->SaveAs(Form("%s/%s_RelErrdecomp_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    TGraphAsymmErrors* graphCombInvYieldsRelStat_2050  = CalculateRelErrUpAsymmGraph( graphCombInvYieldStatPbPb2760GeV_2050, "relativeStatError_2050");
    TGraphAsymmErrors* graphCombInvYieldsRelSys_2050   = CalculateRelErrUpAsymmGraph( graphCombInvYieldSysPbPb2760GeV_2050, "relativeSysError_2050");
    TGraphAsymmErrors* graphCombInvYieldsRelTot_2050 = CalculateRelErrUpAsymmGraph( graphCombInvYieldTotPbPb2760GeV_2050, "relativeTotalError_2050");

    canvasRelTotErr->cd();
      histo2DRelTotErr->Draw("copy");

        if(meson.CompareTo("Pi0")==0 && (thesisPlotting || PaperPi0)){
          while(graphCombInvYieldsRelTot_2050->GetX()[0] < /* 0.8 */ 1.){
            graphCombInvYieldsRelTot_2050->RemovePoint(0);
            graphCombInvYieldsRelStat_2050->RemovePoint(0);
            graphCombInvYieldsRelSys_2050->RemovePoint(0);
          }
        }

        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldsRelTot_2050, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombInvYieldsRelTot_2050->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldsRelStat_2050, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombInvYieldsRelStat_2050->Draw("l,x0,same");
        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldsRelSys_2050, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombInvYieldsRelSys_2050->SetLineStyle(7);
        graphCombInvYieldsRelSys_2050->Draw("l,x0,same");

        legendRelTotErr3->Draw();
        labelRelStatErrEnergy2050->Draw();
        labelRelStatErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/%s_RelErrdecomp_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRelTotErr->SaveAs(Form("%s/%s_RelErrdecomp_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    //**********************************************************************************************************************//
    //************************************* Calculating bin shifted spectra & fitting **************************************//
    //**********************************************************************************************************************//
    if(meson.CompareTo("Pi0")==0 && (thesisPlotting || PaperPi0)){
        while(graphCombInvYieldTotPbPb2760GeV_0010->GetX()[0] < 1.){
          graphCombInvYieldTotPbPb2760GeV_0010->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeV_0010->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeV_0010->RemovePoint(0);
          graphCombInvYieldTotPbPb2760GeV_2050->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeV_2050->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeV_2050->RemovePoint(0);

          graphCombInvYieldTotPbPb2760GeVA_0010->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeVA_0010->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeVA_0010->RemovePoint(0);
          graphCombInvYieldTotPbPb2760GeVA_2050->RemovePoint(0);
          graphCombInvYieldStatPbPb2760GeVA_2050->RemovePoint(0);
          graphCombInvYieldSysPbPb2760GeVA_2050->RemovePoint(0);
        }
    }

    // Cloning spectra for x shift (PCM with MatBudget error)
    TGraphAsymmErrors* graphCombInvYieldTotPbPb2760GeVUnShifted_0010         = (TGraphAsymmErrors*)graphCombInvYieldTotPbPb2760GeV_0010->Clone("Unshifted_0010");
    TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeVUnShifted_0010        = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStat_0010");
    TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeVUnShifted_0010         = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSys_0010");
    TGraphAsymmErrors* graphCombInvYieldTotPbPb2760GeVUnShifted_2050         = (TGraphAsymmErrors*)graphCombInvYieldTotPbPb2760GeV_2050->Clone("Unshifted_2050");
    TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeVUnShifted_2050        = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeV_2050->Clone("UnshiftedStat_2050");
    TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeVUnShifted_2050         = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeV_2050->Clone("UnshiftedSys_2050");

    // Cloning spectra for Y shift (PCM without MatBudget error)
    TGraphAsymmErrors* graphCombInvYieldTotPbPb2760GeVAUnShifted_0010         = (TGraphAsymmErrors*)graphCombInvYieldTotPbPb2760GeVA_0010->Clone("UnshiftedA_0010");
    TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeVAUnShifted_0010        = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeVA_0010->Clone("UnshiftedAStat_0010");
    TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeVAUnShifted_0010         = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeVA_0010->Clone("UnshiftedASys_0010");
    TGraphAsymmErrors* graphCombInvYieldTotPbPb2760GeVAUnShifted_2050         = (TGraphAsymmErrors*)graphCombInvYieldTotPbPb2760GeVA_2050->Clone("UnshiftedA_2050");
    TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeVAUnShifted_2050        = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeVA_2050->Clone("UnshiftedAStat_2050");
    TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeVAUnShifted_2050         = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeVA_2050->Clone("UnshiftedASys_2050");

    if(meson.CompareTo("Pi0")==0){

        graphPCMInvYieldStatPbPb2760GeVUnShifted_0010       = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStatPCMPi0_0010");
        graphPCMInvYieldSysPbPb2760GeVUnShifted_0010        = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSysPCMPi0_0010");
        graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_0010  = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysWOMat2760GeV_0010->Clone("UnshiftedSysforRAAPCMPi0_0010");

        graphPHOSInvYieldStatPbPb2760GeVUnshifted_0010      = (TGraphAsymmErrors*)graphPHOSPi0InvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStatPHOSPi0_0010");
        graphPHOSInvYieldSysPbPb2760GeVUnshifted_0010       = (TGraphAsymmErrors*)graphPHOSPi0InvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSysPHOSPi0_0010");
        graphPHOSInvYieldSysPbPb2760GeVforRAAUnshifted_0010 = (TGraphAsymmErrors*)graphSysErrRAAYieldPi0PHOSPbPb0010->Clone("UnshiftedSysforRAAPHOSPi0_0010");

        graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010     = (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStatEMCalPi0_0010");
        graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010      = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSysEMCalPi0_0010");
        graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_0010 = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_0010->Clone("UnshiftedSysforRAAEMCalPi0_0010");

        graphPCMInvYieldStatPbPb2760GeVUnShifted_2050       = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_2050->Clone("UnshiftedStatPCMPi0_2050");
        graphPCMInvYieldSysPbPb2760GeVUnShifted_2050        = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2050->Clone("UnshiftedSysPCMPi0_2050");
        graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_2050  = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysWOMat2760GeV_2050->Clone("UnshiftedSysWOMatPCMPi0_2050");

        graphPCMInvYieldStatPbPb2760GeVUnShifted_2040       = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_2040->Clone("UnshiftedStatPCMPi0_2040");
        graphPCMInvYieldSysPbPb2760GeVUnShifted_2040        = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2040->Clone("UnshiftedSysPCMPi0_2040");

        graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050     = (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Clone("UnshiftedStatEMCalPi0_2050");
        graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050      = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Clone("UnshiftedSysEMCalPi0_2050");
        graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_2050 = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_2050->Clone("UnshiftedSysforRAAEMCalPi0_2050");


    } else if(meson.CompareTo("Eta")==0){

        graphPCMInvYieldStatPbPb2760GeVUnShifted_0010       = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStatPCMEta_0010");
        graphPCMInvYieldSysPbPb2760GeVUnShifted_0010        = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSysPCMEta_0010");
        graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_0010  = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysWOMat2760GeV_0010->Clone("UnshiftedSysWOMatPCMEta_0010");

        graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010     = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStatEMCalEta_0010");
        graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010      = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSysEMCalEta_0010");
        graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_0010 = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010->Clone("UnshiftedSysforRAAEMCalPi0_0010");

        graphPCMInvYieldStatPbPb2760GeVUnShifted_2040       = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_2040->Clone("UnshiftedStatPCMEta_2040");
        graphPCMInvYieldSysPbPb2760GeVUnShifted_2040        = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2040->Clone("UnshiftedSysPCMEta_2040");

        graphPCMInvYieldStatPbPb2760GeVUnShifted_2050       = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_2050->Clone("UnshiftedStatPCMEta_2050");
        graphPCMInvYieldSysPbPb2760GeVUnShifted_2050        = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2050->Clone("UnshiftedSysPCMEta_2050");
        graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_2050  = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysWOMat2760GeV_2050->Clone("UnshiftedSysWOMatPCMEta_2050");

        graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050     = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Clone("UnshiftedStatEMCalEta_2050");
        graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050      = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Clone("UnshiftedSysEMCalEta_2050");
        graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_2050 = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010->Clone("UnshiftedSysforRAAEMCalPi0_2050");

    }

    if(meson.CompareTo("Pi0")==0 && (thesisPlotting || PaperPi0)){

        while(graphPCMPi0InvYieldStatPbPb2760GeV_0010->GetX()[0] < /* 0.8 */ 1.)
          graphPCMPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(0);
        while(graphPCMPi0InvYieldSysPbPb2760GeV_0010->GetX()[0] < /* 0.8 */ 1.)
          graphPCMPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(0);
        while(graphPCMPi0InvYieldSysWOMat2760GeV_0010->GetX()[0] < /* 0.8 */ 1.)
          graphPCMPi0InvYieldSysWOMat2760GeV_0010->RemovePoint(0);

        while(graphPCMPi0InvYieldStatPbPb2760GeV_2050->GetX()[0] < /* 0.8 */ 1.)
          graphPCMPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(0);
        while(graphPCMPi0InvYieldSysPbPb2760GeV_2050->GetX()[0] < /* 0.8 */ 1.)
          graphPCMPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(0);
        while(graphPCMPi0InvYieldSysWOMat2760GeV_2050->GetX()[0] < /* 0.8 */ 1.)
          graphPCMPi0InvYieldSysWOMat2760GeV_2050->RemovePoint(0);

        while(graphPCMPi0InvYieldStatPbPb2760GeV_2040->GetX()[0] < /* 0.8 */ 1.)
          graphPCMPi0InvYieldStatPbPb2760GeV_2040->RemovePoint(0);
        while(graphPCMPi0InvYieldSysPbPb2760GeV_2040->GetX()[0] < /* 0.8 */ 1.)
          graphPCMPi0InvYieldSysPbPb2760GeV_2040->RemovePoint(0);

    }

    TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeVUnShiftedCopy_0010 = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeVUnShifted_0010->Clone("graphCombInvYieldStatPbPb2760GeVUnShiftedCopy_0010");
    TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeVUnShifted_0010->Clone("graphCombInvYieldSysPbPb2760GeVUnShiftedCopy_0010");
    TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeVUnShiftedCopy_2050 = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeVUnShifted_2050->Clone("graphCombInvYieldStatPbPb2760GeVUnShiftedCopy_2050");
    TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeVUnShiftedCopy_2050  = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeVUnShifted_2050->Clone("graphCombInvYieldSysPbPb2760GeVUnShiftedCopy_2050");

    if(meson.CompareTo("Pi0")==0){

        graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_0010->Clone("graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_0010");
        graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0010->Clone("graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_0010");
        graphPHOSInvYieldStatPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphPHOSPi0InvYieldStatPbPb2760GeV_0010->Clone("graphPHOSInvYieldStatPbPb2760GeVUnShiftedCopy_0010");
        graphPHOSInvYieldSysPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphPHOSPi0InvYieldSysPbPb2760GeV_0010->Clone("graphPHOSInvYieldSysPbPb2760GeVUnShiftedCopy_0010");
        graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Clone("graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_0010");
        graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Clone("graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_0010");

        graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_2050  = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_2040->Clone("graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_2050");
        graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_2050  = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2040->Clone("graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_2050");
        graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_2050  = (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Clone("graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_2050");
        graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_2050  = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Clone("graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_2050");

    } else if(meson.CompareTo("Eta")==0){

        graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_0010->Clone("graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_0010");
        graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0010->Clone("graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_0010");
        graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_0010->Clone("graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_0010");
        graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_0010  = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Clone("graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_0010");

        graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_2050  = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_2040->Clone("graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_2050");
        graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_2050  = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2040->Clone("graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_2050");
        graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_2050  = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Clone("graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_2050");
        graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_2050  = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Clone("graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_2050");

    }



    if(meson.CompareTo("Pi0")==0){

        graphRatioNeutralToChargedComb_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphCombInvYieldStatPbPb2760GeVUnShiftedCopy_0010, graphCombInvYieldSysPbPb2760GeVUnShiftedCopy_0010,
                                                                                               histoChargedPionSpectraStat0010, histoChargedPionSpectraSyst0010,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldPi0CombDummyStatForRatioToCharged_0010, &graphYieldPi0CombDummySystForRatioToCharged_0010,
                                                                                               &graphYieldPionDummyStatForRatioToCharged_0010, &graphYieldPionDummySystForRatioToCharged_0010);
        cout << "ratio comb Pi0 to charged pion" << endl;
        graphRatioNeutralToChargedComb_0010->Print();

        graphRatioPi0PCMToChargedPion_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_0010, graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_0010,
                                                                                               histoChargedPionSpectraStat0010, histoChargedPionSpectraSyst0010,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldPi0PCMDummyStatForRatioToCharged_0010, &graphYieldPi0PCMDummySystForRatioToCharged_0010,
                                                                                               &graphYieldPionDummyStatForRatioToCharged_0010, &graphYieldPionDummySystForRatioToCharged_0010);
        cout << "ratio PCM Pi0 to charged pion" << endl;
        graphRatioPi0PCMToChargedPion_0010->Print();

        graphRatioPi0PHOSToChargedPion_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphPHOSInvYieldStatPbPb2760GeVUnShiftedCopy_0010, graphPHOSInvYieldSysPbPb2760GeVUnShiftedCopy_0010,
                                                                                               histoChargedPionSpectraStat0010, histoChargedPionSpectraSyst0010,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldPi0PHOSDummyStatForRatioToCharged_0010, &graphYieldPi0PHOSDummySystForRatioToCharged_0010,
                                                                                               &graphYieldPionDummyStatForRatioToCharged_0010, &graphYieldPionDummySystForRatioToCharged_0010);
        cout << "ratio PHOS Pi0 to charged pion" << endl;
        graphRatioPi0PHOSToChargedPion_0010->Print();

        graphRatioPi0EMCalToChargedPion_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_0010, graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_0010,
                                                                                               histoChargedPionSpectraStat0010, histoChargedPionSpectraSyst0010,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldPi0EMCalDummyStatForRatioToCharged_0010, &graphYieldPi0EMCalDummySystForRatioToCharged_0010,
                                                                                               &graphYieldPionDummyStatForRatioToCharged_0010, &graphYieldPionDummySystForRatioToCharged_0010);
        cout << "ratio EMCal Pi0 to charged pion" << endl;
        graphRatioPi0EMCalToChargedPion_0010->Print();

        graphRatioNeutralToChargedComb_2050 = CalculateRatioBetweenSpectraWithDifferentBinning(graphCombInvYieldStatPbPb2760GeVUnShiftedCopy_2050, graphCombInvYieldSysPbPb2760GeVUnShiftedCopy_2050,
                                                                                               histoChargedPionSpectraStat2040, histoChargedPionSpectraSyst2040,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldPi0CombDummyStatForRatioToCharged_2050, &graphYieldPi0CombDummySystForRatioToCharged_2050,
                                                                                               &graphYieldPionDummyStatForRatioToCharged_2040, &graphYieldPionDummySystForRatioToCharged_2040);
        cout << "ratio comb Pi0 to charged pion" << endl;
        graphRatioNeutralToChargedComb_2050->Print();

        graphRatioPi0PCMToChargedPion_2050 = CalculateRatioBetweenSpectraWithDifferentBinning(graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_2050, graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_2050,
                                                                                               histoChargedPionSpectraStat2040, histoChargedPionSpectraSyst2040,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldPi0PCMDummyStatForRatioToCharged_2050, &graphYieldPi0PCMDummySystForRatioToCharged_2050,
                                                                                               &graphYieldPionDummyStatForRatioToCharged_2040, &graphYieldPionDummySystForRatioToCharged_2040);
        cout << "ratio PCM Pi0 to charged pion" << endl;
        graphRatioPi0PCMToChargedPion_2050->Print();

        graphRatioPi0EMCalToChargedPion_2050 = CalculateRatioBetweenSpectraWithDifferentBinning(graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_2050, graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_2050,
                                                                                               histoChargedPionSpectraStat2040, histoChargedPionSpectraSyst2040,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldPi0EMCalDummyStatForRatioToCharged_2050, &graphYieldPi0EMCalDummySystForRatioToCharged_2050,
                                                                                               &graphYieldPionDummyStatForRatioToCharged_2040, &graphYieldPionDummySystForRatioToCharged_2040);
        cout << "ratio EMCal Pi0 to charged pion" << endl;
        graphRatioPi0EMCalToChargedPion_2050->Print();


    }
    else if(meson.CompareTo("Eta")==0){

        graphRatioNeutralToChargedComb_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphCombInvYieldStatPbPb2760GeVUnShiftedCopy_0010, graphCombInvYieldSysPbPb2760GeVUnShiftedCopy_0010,
                                                                                               histoChargedKaonSpectraStat0010, histoChargedKaonSpectraSyst0010,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldEtaCombDummyStatForRatioToCharged_0010, &graphYieldEtaCombDummySystForRatioToCharged_0010,
                                                                                               &graphYieldKaonDummyStatForRatioToCharged_0010, &graphYieldKaonDummySystForRatioToCharged_0010);
        cout << "ratio comb Eta to charged pion" << endl;
        graphRatioNeutralToChargedComb_0010->Print();

        graphRatioEtaPCMToChargedKaon_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_0010, graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_0010,
                                                                                               histoChargedKaonSpectraStat0010, histoChargedKaonSpectraSyst0010,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldEtaPCMDummyStatForRatioToCharged_0010, &graphYieldEtaPCMDummySystForRatioToCharged_0010,
                                                                                               &graphYieldKaonDummyStatForRatioToCharged_0010, &graphYieldKaonDummySystForRatioToCharged_0010);
        cout << "ratio PCM Eta to charged pion" << endl;
        graphRatioEtaPCMToChargedKaon_0010->Print();

        graphRatioEtaEMCalToChargedKaon_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_0010, graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_0010,
                                                                                               histoChargedKaonSpectraStat0010, histoChargedKaonSpectraSyst0010,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldEtaEMCalDummyStatForRatioToCharged_0010, &graphYieldEtaEMCalDummySystForRatioToCharged_0010,
                                                                                               &graphYieldKaonDummyStatForRatioToCharged_0010, &graphYieldKaonDummySystForRatioToCharged_0010);
        cout << "ratio EMCal Eta to charged pion" << endl;
        graphRatioEtaEMCalToChargedKaon_0010->Print();

        graphRatioNeutralToChargedComb_2050 = CalculateRatioBetweenSpectraWithDifferentBinning(graphCombInvYieldStatPbPb2760GeVUnShiftedCopy_2050, graphCombInvYieldSysPbPb2760GeVUnShiftedCopy_2050,
                                                                                               histoChargedKaonSpectraStat2040, histoChargedKaonSpectraSyst2040,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldEtaCombDummyStatForRatioToCharged_2050, &graphYieldEtaCombDummySystForRatioToCharged_2050,
                                                                                               &graphYieldKaonDummyStatForRatioToCharged_2040, &graphYieldKaonDummySystForRatioToCharged_2040);
        cout << "ratio comb Eta to charged pion" << endl;
        graphRatioNeutralToChargedComb_2050->Print();

        graphRatioEtaPCMToChargedKaon_2050 = CalculateRatioBetweenSpectraWithDifferentBinning(graphPCMInvYieldStatPbPb2760GeVUnShiftedCopy_2050, graphPCMInvYieldSysPbPb2760GeVUnShiftedCopy_2050,
                                                                                               histoChargedKaonSpectraStat2040, histoChargedKaonSpectraSyst2040,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldEtaPCMDummyStatForRatioToCharged_2050, &graphYieldEtaPCMDummySystForRatioToCharged_2050,
                                                                                               &graphYieldKaonDummyStatForRatioToCharged_2040, &graphYieldKaonDummySystForRatioToCharged_2040);
        cout << "ratio PCM Eta to charged pion" << endl;
        graphRatioEtaPCMToChargedKaon_2050->Print();

        graphRatioEtaEMCalToChargedKaon_2050 = CalculateRatioBetweenSpectraWithDifferentBinning(graphEMCalInvYieldStatPbPb2760GeVUnShiftedCopy_2050, graphEMCalInvYieldSysPbPb2760GeVUnShiftedCopy_2050,
                                                                                               histoChargedKaonSpectraStat2040, histoChargedKaonSpectraSyst2040,
                                                                                               kTRUE,  kTRUE,
                                                                                               &graphYieldEtaEMCalDummyStatForRatioToCharged_2050, &graphYieldEtaEMCalDummySystForRatioToCharged_2050,
                                                                                               &graphYieldKaonDummyStatForRatioToCharged_2040, &graphYieldKaonDummySystForRatioToCharged_2040);
        cout << "ratio EMCal Eta to charged pion" << endl;
        graphRatioEtaEMCalToChargedKaon_2050->Print();

    }

    if(thesisPlotting || PaperPi0) minfitPt = /*0.8*/1.0;
    else minfitPt = 0.4;

    // Calculating binshifts
    if(meson.CompareTo("Pi0")==0){

        cout << __LINE__ << endl;
        //Bylinkin
        Double_t paramTCM_0010[5] = {graphCombInvYieldTotPbPb2760GeVUnShifted_0010->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeVUnShifted_0010->GetY()[0],0.3,8};
        fitBylinkinPbPb2760GeVPtLHC11h_0010 = FitObject("tcm","BylinkinFitPi00010","Pi0",graphCombInvYieldTotPbPb2760GeVUnShifted_0010,minfitPt,20.,paramTCM_0010);
        graphCombInvYieldTotPbPb2760GeVUnShifted_0010->Fit(fitBylinkinPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",minfitPt,20.);
        Double_t paramTCM_2050[5] = {graphCombInvYieldTotPbPb2760GeVUnShifted_2050->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeVUnShifted_2050->GetY()[0],0.3,8};
        fitBylinkinPbPb2760GeVPtLHC11h_2050 = FitObject("tcm","BylinkinFitPi02050","Pi0",graphCombInvYieldTotPbPb2760GeVUnShifted_2050,minfitPt,20.,paramTCM_2050);
        graphCombInvYieldTotPbPb2760GeVUnShifted_2050->Fit(fitBylinkinPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",minfitPt,20.);
//         cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_0010)<< endl << endl;
//         cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_2050)<< endl << endl;

        //QCD
        Double_t paramQCD_0010[5] = {24,5.,-20.,2.,20};
        Double_t paramQCD_2050[5] = {24,5.,-20.,2.,20};
        fitQCDPbPb2760GeVPtLHC11h_0010 = FitObject("qcd","fitQCD","Pi0",graphCombInvYieldTotPbPb2760GeVUnShifted_0010,minfitPt,20.,paramQCD_0010);
        graphCombInvYieldTotPbPb2760GeVUnShifted_0010->Fit(fitQCDPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",minfitPt,20.);
        fitQCDPbPb2760GeVPtLHC11h_2050 = FitObject("qcd","fitQCD","Pi0",graphCombInvYieldTotPbPb2760GeVUnShifted_2050,minfitPt,20.,paramQCD_2050);
        graphCombInvYieldTotPbPb2760GeVUnShifted_2050->Fit(fitQCDPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",minfitPt,20.);
//         cout << WriteParameterToFile(fitQCDPbPb2760GeVPtLHC11h_0010)<< endl << endl;
//         cout << WriteParameterToFile(fitQCDPbPb2760GeVPtLHC11h_2050)<< endl << endl;

        //Tsallis
//         fitTsallisPbPb2760GeVPtLHC11h_0010 = FitObject("l","fitTsallis","Pi0",graphCombInvYieldTotPbPb2760GeVUnShifted_0010,1.0,20.);
//         fitTsallisPbPb2760GeVPtLHC11h_2050 = FitObject("l","fitTsallis","Pi0",graphCombInvYieldTotPbPb2760GeVUnShifted_2050,1.0,20.);
//         graphCombInvYieldTotPbPb2760GeVUnShifted_0010->Fit(fitTsallisPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",1.0,20.);
//         graphCombInvYieldTotPbPb2760GeVUnShifted_2050->Fit(fitTsallisPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",1.0,20.);
//         cout << WriteParameterToFile(fitTsallisPbPb2760GeVPtLHC11h_0010)<< endl << endl;
//         cout << WriteParameterToFile(fitTsallisPbPb2760GeVPtLHC11h_2050)<< endl << endl;

        Int_t binsPi0PCM,binsPi0PHOS,binsPi0EMCal;
        if(thesisPlotting || PaperPi0){
          binsPi0PCM = /*18*/17;
          binsPi0PHOS = /*1*/0;
          binsPi0EMCal = /*12*/11;
        } else {
          binsPi0PCM = 19;
          binsPi0PHOS = 3;
          binsPi0EMCal = 14;
        }

        if(bWCorrection.CompareTo("X")==0 ){
//          cout << __LINE__ << endl;
            //shifting in X
            //Bylinkin
            TF1* fitFunctionShiftingX_0010 = FitObject("tcmpt","tcmptPi00010","Pi0",graphCombInvYieldTotPbPb2760GeV_0010);

            //combined
            graphCombInvYieldTotPbPb2760GeV_0010         = ApplyXshift(graphCombInvYieldTotPbPb2760GeV_0010, fitFunctionShiftingX_0010,"Pi0");

            graphCombInvYieldStatPbPb2760GeV_0010        = ApplyXshiftIndividualSpectra (graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphCombInvYieldStatPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            0, graphCombInvYieldTotPbPb2760GeV_0010->GetN(),"Pi0");

            graphCombInvYieldSysPbPb2760GeV_0010         = ApplyXshiftIndividualSpectra (graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphCombInvYieldSysPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            0, graphCombInvYieldTotPbPb2760GeV_0010->GetN(),"Pi0");

            //  PCM
            graphPCMPi0InvYieldStatPbPb2760GeV_0010         = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphPCMPi0InvYieldStatPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            0, binsPi0PCM,"Pi0");

            graphPCMPi0InvYieldSysPbPb2760GeV_0010          = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphPCMPi0InvYieldSysPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            0, binsPi0PCM,"Pi0");

            //  PHOS
            graphPHOSPi0InvYieldStatPbPb2760GeV_0010        = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphPHOSPi0InvYieldStatPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            binsPi0PHOS, binsPi0PHOS+16,"Pi0");

            graphPHOSPi0InvYieldSysPbPb2760GeV_0010         = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphPHOSPi0InvYieldSysPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            binsPi0PHOS, binsPi0PHOS+16,"Pi0");

            //  EMCal
            graphEMCalPi0InvYieldStatPbPb2760GeV_0010       = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphEMCalPi0InvYieldStatPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            binsPi0EMCal, binsPi0EMCal+8,"Pi0");

            graphEMCalPi0InvYieldSysPbPb2760GeV_0010        = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphEMCalPi0InvYieldSysPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            binsPi0EMCal, binsPi0EMCal+8,"Pi0");


            TF1* fitFunctionShiftingX_2050 = FitObject("tcmpt","tcmptPi02050","Pi0",graphCombInvYieldStatPbPb2760GeV_2050);

            //combined
            graphCombInvYieldTotPbPb2760GeV_2050         = ApplyXshift(graphCombInvYieldTotPbPb2760GeV_2050, fitFunctionShiftingX_2050,"Pi0");

            graphCombInvYieldStatPbPb2760GeV_2050        = ApplyXshiftIndividualSpectra (graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphCombInvYieldStatPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            0, graphCombInvYieldTotPbPb2760GeV_2050->GetN(),"Pi0");

            graphCombInvYieldSysPbPb2760GeV_2050         = ApplyXshiftIndividualSpectra (graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphCombInvYieldSysPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            0, graphCombInvYieldTotPbPb2760GeV_2050->GetN(),"Pi0");

            //  PCM
            graphPCMPi0InvYieldStatPbPb2760GeV_2050         = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphPCMPi0InvYieldStatPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            0, binsPi0PCM,"Pi0");

            graphPCMPi0InvYieldSysPbPb2760GeV_2050          = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphPCMPi0InvYieldSysPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            0, binsPi0PCM,"Pi0");

            //  EMCal
            graphEMCalPi0InvYieldStatPbPb2760GeV_2050       = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphEMCalPi0InvYieldStatPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            binsPi0EMCal, binsPi0EMCal+8,"Pi0");

            graphEMCalPi0InvYieldSysPbPb2760GeV_2050        = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphEMCalPi0InvYieldSysPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            binsPi0EMCal, binsPi0EMCal+8,"Pi0");


            TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1350,1350*1.15);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy2, 0.16, 0.02, 0.02, 0.09);
            canvasDummy2->SetLogx();
            canvasDummy2->SetLogy();

            TH2F * histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",11000,minPtRange-0.2,maxPtRange,1000,2e-8,1e4);
                SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
            histo2DDummy2->GetXaxis()->SetMoreLogLabels();
            histo2DDummy2->GetXaxis()->SetLabelOffset(-0.01);
            histo2DDummy2->Draw("copy");

            fitBylinkinPbPb2760GeVPtLHC11h_0010->SetLineColor(kBlue);
            fitBylinkinPbPb2760GeVPtLHC11h_2050->SetLineColor(kBlue);
            fitQCDPbPb2760GeVPtLHC11h_0010->SetLineColor(kGreen+2);
            fitQCDPbPb2760GeVPtLHC11h_2050->SetLineColor(kGreen+2);
//             fitTsallisPbPb2760GeVPtLHC11h_0010->SetLineColor(kMagenta-7);
//             fitTsallisPbPb2760GeVPtLHC11h_2050->SetLineColor(kMagenta-7);
            fitBylinkinPbPb2760GeVPtLHC11h_0010->Draw("same");
            fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");
//             fitQCDPbPb2760GeVPtLHC11h_0010->Draw("same");
//             fitQCDPbPb2760GeVPtLHC11h_2050->Draw("same");
//             fitTsallisPbPb2760GeVPtLHC11h_0010->Draw("same");
//             fitTsallisPbPb2760GeVPtLHC11h_2050->Draw("same");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVUnShifted_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVUnShifted_0010->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_0010, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeV_0010->Draw("pEsame");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVUnShifted_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVUnShifted_2050->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_2050, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeV_2050->Draw("pEsame");

            TLegend* legendXdummy = new TLegend(0.55,0.85,0.85,0.95);
            legendXdummy->SetFillColor(0);
            legendXdummy->SetLineColor(0);
            legendXdummy->SetTextFont(42);
            legendXdummy->SetTextSize(FontSize);
            legendXdummy->AddEntry(graphCombInvYieldStatPbPb2760GeVUnShifted_0010,"combined unshifted","lp");
            legendXdummy->AddEntry(graphCombInvYieldStatPbPb2760GeV_0010,"combined shifted","lp");
            legendXdummy->AddEntry(fitBylinkinPbPb2760GeVPtLHC11h_0010,"Bylinkin","l");
            legendXdummy->AddEntry(fitQCDPbPb2760GeVPtLHC11h_0010,"QCD","l");
//             legendXdummy->AddEntry(fitTsallisPbPb2760GeVPtLHC11h_0010,"Tsallis","l");
            legendXdummy->Draw();

            canvasDummy2->Update();
            canvasDummy2->Print(Form("%s/%s_ComparisonXShifted.%s",outputDir.Data(),meson.Data(),suffix.Data()));
            delete canvasDummy2;

            //shifting in Y
            TF1* fitFunctionShiftingY_0010 = (TF1*)fitFunctionShiftingX_0010->Clone("fitFunctionShiftingY_0010");
            TF1* fitFunctionShiftingY_2050 = (TF1*)fitFunctionShiftingX_2050->Clone("fitFunctionShiftingY_2050");

            cout << "combined binshift Y" << endl;
            graphCombInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeVAUnShifted_0010->Clone("YShiftedCombSys0010");
            graphCombInvYieldSysPbPb2760GeVYShifted_0010= ApplyYshiftIndividualSpectra( graphCombInvYieldSysPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);
            graphCombInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeVAUnShifted_0010->Clone("YShiftedCombStat0010");
            graphCombInvYieldStatPbPb2760GeVYShifted_0010= ApplyYshiftIndividualSpectra( graphCombInvYieldStatPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);

            cout << "PCM binshift Y" << endl;
            graphPCMInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_0010->Clone("YShiftedPCMSys0010");
            graphPCMInvYieldSysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMInvYieldSysPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);
            graphPCMInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMInvYieldStatPbPb2760GeVUnShifted_0010->Clone("YShiftedPCMStat0010");
            graphPCMInvYieldStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMInvYieldStatPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);

            cout << "PHOS binshift Y" << endl;
            graphPHOSInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPHOSInvYieldSysPbPb2760GeVforRAAUnshifted_0010->Clone("YShiftedPHOSSys0010");
            graphPHOSInvYieldSysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPHOSInvYieldSysPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);
            graphPHOSInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPHOSInvYieldStatPbPb2760GeVUnshifted_0010->Clone("YShiftedPHOSStat0010");
            graphPHOSInvYieldStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPHOSInvYieldStatPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);

            cout << "EMCal binshift Y" << endl;
            graphEMCalInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_0010->Clone("YShiftedEMCalSys0010");
            graphEMCalInvYieldSysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldSysPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);
            graphEMCalInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010->Clone("YShiftedEMCalStat0010");
            graphEMCalInvYieldStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldStatPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);

            cout << "combined binshift Y" << endl;
            graphCombInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeVAUnShifted_2050->Clone("YShiftedCombSys2050");
            graphCombInvYieldSysPbPb2760GeVYShifted_2050= ApplyYshiftIndividualSpectra( graphCombInvYieldSysPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);
            graphCombInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeVAUnShifted_2050->Clone("YShiftedCombStat2050");
            graphCombInvYieldStatPbPb2760GeVYShifted_2050= ApplyYshiftIndividualSpectra( graphCombInvYieldStatPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);

            cout << "PCM binshift Y" << endl;
            graphPCMInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_2050->Clone("YShiftedPCMSys2050");
            graphPCMInvYieldSysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMInvYieldSysPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);
            graphPCMInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMInvYieldStatPbPb2760GeVUnShifted_2050->Clone("YShiftedPCMStat2050");
            graphPCMInvYieldStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMInvYieldStatPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);

            cout << "EMCal binshift Y" << endl;
            graphEMCalInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_2050->Clone("YShiftedEMCalSys2050");
            graphEMCalInvYieldSysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldSysPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);
            graphEMCalInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050->Clone("YShiftedEMCalStat2050");
            graphEMCalInvYieldStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);

            TCanvas* canvasDummy3 = new TCanvas("canvasDummy3","",200,10,1350,1350*1.15);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy3, 0.16, 0.02, 0.02, 0.09);
            canvasDummy3->SetLogx();
            canvasDummy3->SetLogy();

            TH2F * histo2DDummy3 = new TH2F("histo2DDummy3","histo2DDummy3",11000,minPtRange-0.2,maxPtRange,1000,2e-8,1e4);
                SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
            histo2DDummy3->GetXaxis()->SetMoreLogLabels();
            histo2DDummy3->GetXaxis()->SetLabelOffset(-0.01);
            histo2DDummy3->Draw("copy");

            fitBylinkinPbPb2760GeVPtLHC11h_0010->SetLineColor(kBlue);
            fitBylinkinPbPb2760GeVPtLHC11h_0010->Draw("same");
            fitBylinkinPbPb2760GeVPtLHC11h_2050->SetLineColor(kBlue);
            fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVAUnShifted_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVAUnShifted_0010->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVAUnShifted_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVAUnShifted_2050->Draw("pEsame");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVYShifted_0010, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVYShifted_0010->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVYShifted_2050, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVYShifted_2050->Draw("pEsame");

            TLegend* legendYdummy = new TLegend(0.55,0.8,0.85,0.95);
            legendYdummy->SetFillColor(0);
            legendYdummy->SetLineColor(0);
            legendYdummy->SetTextFont(42);
            legendYdummy->SetTextSize(FontSize);
            legendYdummy->AddEntry(graphCombInvYieldStatPbPb2760GeVAUnShifted_0010,"combined unshifted","lp");
            legendYdummy->AddEntry(graphCombInvYieldStatPbPb2760GeVYShifted_0010,"combined shifted","lp");
            legendYdummy->Draw();

            canvasDummy3->Update();
            canvasDummy3->Print(Form("%s/%s_ComparisonYShifted.%s",outputDir.Data(),meson.Data(),suffix.Data()));


            canvasDummy3->cd();
            histo2DDummy3->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeV_0010->Draw("pEsame");
            graphCombInvYieldStatPbPb2760GeV_2050->Draw("pEsame");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVYShifted_0010, 21, 1.5, kBlue, kBlue, widthLinesBoxes, kTRUE);
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVYShifted_2050, 21, 1.5, kBlue, kBlue, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVYShifted_0010->Draw("pEsame");
            graphCombInvYieldStatPbPb2760GeVYShifted_2050->Draw("pEsame");

            TLegend* legendXYdummy = new TLegend(0.6,0.8,0.9,0.95);
            legendXYdummy->SetFillColor(0);
            legendXYdummy->SetLineColor(0);
            legendXYdummy->SetTextFont(42);
            legendXYdummy->SetTextSize(FontSize);
            legendXYdummy->AddEntry(graphCombInvYieldStatPbPb2760GeV_0010,"x-shifted","lp");
            legendXYdummy->AddEntry(graphCombInvYieldStatPbPb2760GeVYShifted_0010,"y-shifted","lp");
            legendXYdummy->Draw();

            canvasDummy3->Update();
            canvasDummy3->Print(Form("%s/%s_ComparisonXYShifted.%s",outputDir.Data(),meson.Data(),suffix.Data()));

            delete canvasDummy3;

        }


    }
    else if(meson.CompareTo("Eta")==0){

        //Bylinkin
        Double_t paramTCM_0010[5] = {graphCombInvYieldTotPbPb2760GeVUnShifted_0010->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeVUnShifted_0010->GetY()[0],0.3,8};
        fitBylinkinPbPb2760GeVPtLHC11h_0010 = FitObject("tcm","BylinkinFitPi00010","Eta",graphCombInvYieldTotPbPb2760GeVUnShifted_0010,1.0,20.,paramTCM_0010);
        graphCombInvYieldTotPbPb2760GeVUnShifted_0010->Fit(fitBylinkinPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",1.0,20.);
        Double_t paramTCM_2050[5] = {graphCombInvYieldTotPbPb2760GeVUnShifted_2050->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeVUnShifted_2050->GetY()[0],0.3,8};
        fitBylinkinPbPb2760GeVPtLHC11h_2050 = FitObject("tcm","BylinkinFitPi02050","Eta",graphCombInvYieldTotPbPb2760GeVUnShifted_2050,1.0,20.,paramTCM_2050);
        graphCombInvYieldTotPbPb2760GeVUnShifted_2050->Fit(fitBylinkinPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",1.0,20.);
//         cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_0010)<< endl << endl;
//         cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_2050)<< endl << endl;

        //QCD
        Double_t paramQCD_0010[5] = {24,5.,-20.,2.,20};
        Double_t paramQCD_2050[5] = {24,5.,-20.,2.,20};
        fitQCDPbPb2760GeVPtLHC11h_0010 = FitObject("qcd","fitQCD","Eta",graphCombInvYieldTotPbPb2760GeVUnShifted_0010,1.0,20.,paramQCD_0010);
        graphCombInvYieldTotPbPb2760GeVUnShifted_0010->Fit(fitQCDPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",1.0,20.);
        fitQCDPbPb2760GeVPtLHC11h_2050 = FitObject("qcd","fitQCD","Eta",graphCombInvYieldTotPbPb2760GeVUnShifted_2050,1.0,20.,paramQCD_2050);
        graphCombInvYieldTotPbPb2760GeVUnShifted_2050->Fit(fitQCDPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",1.0,20.);
//         cout << WriteParameterToFile(fitQCDPbPb2760GeVPtLHC11h_0010)<< endl << endl;
//         cout << WriteParameterToFile(fitQCDPbPb2760GeVPtLHC11h_2050)<< endl << endl;

        //Tsallis
//         fitTsallisPbPb2760GeVPtLHC11h_0010 = FitObject("l","fitTsallis","Eta",graphCombInvYieldTotPbPb2760GeVUnShifted_0010,1.0,20.);
//         fitTsallisPbPb2760GeVPtLHC11h_2050 = FitObject("l","fitTsallis","Eta",graphCombInvYieldTotPbPb2760GeVUnShifted_2050,1.0,20.);
//         graphCombInvYieldTotPbPb2760GeVUnShifted_0010->Fit(fitTsallisPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",1.0,20.);
//         graphCombInvYieldTotPbPb2760GeVUnShifted_2050->Fit(fitTsallisPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",1.0,20.);
//         cout << WriteParameterToFile(fitTsallisPbPb2760GeVPtLHC11h_0010)<< endl << endl;
//         cout << WriteParameterToFile(fitTsallisPbPb2760GeVPtLHC11h_2050)<< endl << endl;

        if(bWCorrection.CompareTo("X")==0 ){

            // binshift in x
            TF1* fitFunctionShiftingX_0010 = FitObject("tcmpt","tcmptEta0010","Eta",graphCombInvYieldStatPbPb2760GeV_0010);

            // combined
            graphCombInvYieldTotPbPb2760GeV_0010         = ApplyXshift(graphCombInvYieldTotPbPb2760GeV_0010, fitFunctionShiftingX_0010,"Eta");

            graphCombInvYieldStatPbPb2760GeV_0010        = ApplyXshiftIndividualSpectra (graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphCombInvYieldStatPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            0, graphCombInvYieldTotPbPb2760GeV_0010->GetN(),"Eta");
            graphCombInvYieldSysPbPb2760GeV_0010         = ApplyXshiftIndividualSpectra (graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphCombInvYieldSysPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            0, graphCombInvYieldTotPbPb2760GeV_0010->GetN(),"Eta");

            // PCM
            graphPCMEtaInvYieldStatPbPb2760GeV_0010         = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphPCMEtaInvYieldStatPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            0, 7,"Eta");
            graphPCMEtaInvYieldSysPbPb2760GeV_0010          = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphPCMEtaInvYieldSysPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            0, 7,"Eta");

            // EMCal
            graphEMCalEtaInvYieldStatPbPb2760GeV_0010       = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphEMCalEtaInvYieldStatPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            4, 11,"Eta");

            graphEMCalEtaInvYieldSysPbPb2760GeV_0010        = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_0010,
                                                                                            graphEMCalEtaInvYieldSysPbPb2760GeV_0010,
                                                                                            fitFunctionShiftingX_0010,
                                                                                            4, 11,"Eta");

            TF1* fitFunctionShiftingX_2050 = FitObject("tcmpt","tcmptEta2050","Eta",graphCombInvYieldStatPbPb2760GeV_2050);

            // combined
            graphCombInvYieldTotPbPb2760GeV_2050         = ApplyXshift(graphCombInvYieldTotPbPb2760GeV_2050, fitFunctionShiftingX_2050,"Eta");

            graphCombInvYieldStatPbPb2760GeV_2050        = ApplyXshiftIndividualSpectra (graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphCombInvYieldStatPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            0, graphCombInvYieldTotPbPb2760GeV_2050->GetN(),"Eta");

            graphCombInvYieldSysPbPb2760GeV_2050         = ApplyXshiftIndividualSpectra (graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphCombInvYieldSysPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            0, graphCombInvYieldTotPbPb2760GeV_2050->GetN(),"Eta");

            // PCM
            graphPCMEtaInvYieldStatPbPb2760GeV_2050         = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphPCMEtaInvYieldStatPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            0, 7,"Eta");

            graphPCMEtaInvYieldSysPbPb2760GeV_2050          = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphPCMEtaInvYieldSysPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            0, 7,"Eta");

            // EMCal
            graphEMCalEtaInvYieldStatPbPb2760GeV_2050       = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphEMCalEtaInvYieldStatPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            4, 11,"Eta");

            graphEMCalEtaInvYieldSysPbPb2760GeV_2050        = ApplyXshiftIndividualSpectra( graphCombInvYieldTotPbPb2760GeV_2050,
                                                                                            graphEMCalEtaInvYieldSysPbPb2760GeV_2050,
                                                                                            fitFunctionShiftingX_2050,
                                                                                            4, 11,"Eta");


            TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1350,1350*1.15);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy2, 0.16, 0.02, 0.02, 0.09);
            canvasDummy2->SetLogx();
            canvasDummy2->SetLogy();

            TH2F * histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",11000,minPtRange-0.2,maxPtRange,1000,2e-8,1e4);
                SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
            histo2DDummy2->GetXaxis()->SetMoreLogLabels();
            histo2DDummy2->GetXaxis()->SetLabelOffset(-0.01);
            histo2DDummy2->Draw("copy");

            fitBylinkinPbPb2760GeVPtLHC11h_0010->SetLineColor(kBlue);
            fitBylinkinPbPb2760GeVPtLHC11h_2050->SetLineColor(kBlue);
            fitQCDPbPb2760GeVPtLHC11h_0010->SetLineColor(kGreen+2);
            fitQCDPbPb2760GeVPtLHC11h_2050->SetLineColor(kGreen+2);
//             fitTsallisPbPb2760GeVPtLHC11h_0010->SetLineColor(kMagenta-7);
//             fitTsallisPbPb2760GeVPtLHC11h_2050->SetLineColor(kMagenta-7);
            fitBylinkinPbPb2760GeVPtLHC11h_0010->Draw("same");
            fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");
            fitQCDPbPb2760GeVPtLHC11h_0010->Draw("same");
            fitQCDPbPb2760GeVPtLHC11h_2050->Draw("same");
//             fitTsallisPbPb2760GeVPtLHC11h_0010->Draw("same");
//             fitTsallisPbPb2760GeVPtLHC11h_2050->Draw("same");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVUnShifted_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVUnShifted_0010->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_0010, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeV_0010->Draw("pEsame");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVUnShifted_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVUnShifted_2050->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_2050, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeV_2050->Draw("pEsame");

            TLegend* legendXdummy = new TLegend(0.55,0.8,0.85,0.95);
            legendXdummy->SetFillColor(0);
            legendXdummy->SetLineColor(0);
            legendXdummy->SetTextFont(42);
            legendXdummy->SetTextSize(FontSize);
            legendXdummy->AddEntry(graphCombInvYieldStatPbPb2760GeVUnShifted_0010,"combined unshifted","lp");
            legendXdummy->AddEntry(graphCombInvYieldStatPbPb2760GeV_0010,"combined shifted","lp");
            legendXdummy->AddEntry(fitBylinkinPbPb2760GeVPtLHC11h_0010,"Bylinkin","l");
            legendXdummy->AddEntry(fitQCDPbPb2760GeVPtLHC11h_0010,"QCD","l");
//             legendXdummy->AddEntry(fitTsallisPbPb2760GeVPtLHC11h_0010,"Tsallis","l");
            legendXdummy->Draw();

            canvasDummy2->Update();
            canvasDummy2->Print(Form("%s/%s_ComparisonShiftedEta_PbPbPbPb2760GeV.%s",outputDir.Data(),meson.Data(),suffix.Data()));
            delete canvasDummy2;

            //shifting in Y
            TF1* fitFunctionShiftingY_0010 = (TF1*)fitFunctionShiftingX_0010->Clone("fitFunctionShiftingY_0010");
            TF1* fitFunctionShiftingY_2050 = (TF1*)fitFunctionShiftingX_2050->Clone("fitFunctionShiftingY_2050");

            cout << "combined binshift Y" << endl;
            graphCombInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeVAUnShifted_0010->Clone("YShiftedCombSys0010");
            graphCombInvYieldSysPbPb2760GeVYShifted_0010= ApplyYshiftIndividualSpectra( graphCombInvYieldSysPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);
            graphCombInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeVAUnShifted_0010->Clone("YShiftedCombStat0010");
            graphCombInvYieldStatPbPb2760GeVYShifted_0010= ApplyYshiftIndividualSpectra( graphCombInvYieldStatPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);

            cout << "PCM binshift Y" << endl;
            graphPCMInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_0010->Clone("YShiftedPCMSys0010");
            graphPCMInvYieldSysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMInvYieldSysPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);
            graphPCMInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMInvYieldStatPbPb2760GeVUnShifted_0010->Clone("YShiftedPCMStat0010");
            graphPCMInvYieldStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMInvYieldStatPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);

            cout << "EMCal binshift Y" << endl;
            graphEMCalInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_0010->Clone("YShiftedEMCalSys0010");
            graphEMCalInvYieldSysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldSysPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);
            graphEMCalInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010->Clone("YShiftedEMCalStat0010");
            graphEMCalInvYieldStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldStatPbPb2760GeVYShifted_0010, fitFunctionShiftingY_0010);

            cout << "combined binshift Y" << endl;
            graphCombInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeVAUnShifted_2050->Clone("YShiftedCombSys2050");
            graphCombInvYieldSysPbPb2760GeVYShifted_2050= ApplyYshiftIndividualSpectra( graphCombInvYieldSysPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);
            graphCombInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeVAUnShifted_2050->Clone("YShiftedCombStat2050");
            graphCombInvYieldStatPbPb2760GeVYShifted_2050= ApplyYshiftIndividualSpectra( graphCombInvYieldStatPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);

            cout << "PCM binshift Y" << endl;
            graphPCMInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_2050->Clone("YShiftedPCMSys2050");
            graphPCMInvYieldSysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMInvYieldSysPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);
            graphPCMInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMInvYieldStatPbPb2760GeVUnShifted_2050->Clone("YShiftedPCMStat2050");
            graphPCMInvYieldStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMInvYieldStatPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);

            cout << "EMCal binshift Y" << endl;
            graphEMCalInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_2050->Clone("YShiftedEMCalSys2050");
            graphEMCalInvYieldSysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldSysPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);
            graphEMCalInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050->Clone("YShiftedEMCalStat2050");
            graphEMCalInvYieldStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, fitFunctionShiftingY_2050);


            TCanvas* canvasDummy3 = new TCanvas("canvasDummy3","",200,10,1350,1350*1.15);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy3, 0.16, 0.02, 0.02, 0.09);
            canvasDummy3->SetLogx();
            canvasDummy3->SetLogy();

            TH2F * histo2DDummy3 = new TH2F("histo2DDummy3","histo2DDummy3",11000,minPtRange-0.2,maxPtRange,1000,2e-8,1e4);
                SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
            histo2DDummy3->GetXaxis()->SetMoreLogLabels();
            histo2DDummy3->GetXaxis()->SetLabelOffset(-0.01);
            histo2DDummy3->Draw("copy");

            fitBylinkinPbPb2760GeVPtLHC11h_0010->SetLineColor(kBlue);
            fitBylinkinPbPb2760GeVPtLHC11h_0010->Draw("same");
            fitBylinkinPbPb2760GeVPtLHC11h_2050->SetLineColor(kBlue);
            fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVAUnShifted_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVAUnShifted_0010->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVAUnShifted_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVAUnShifted_2050->Draw("pEsame");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVYShifted_0010, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVYShifted_0010->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVYShifted_2050, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVYShifted_2050->Draw("pEsame");

            TLegend* legendYdummy = new TLegend(0.55,0.8,0.85,0.95);
            legendYdummy->SetFillColor(0);
            legendYdummy->SetLineColor(0);
            legendYdummy->SetTextFont(42);
            legendYdummy->SetTextSize(FontSize);
            legendYdummy->AddEntry(graphCombInvYieldStatPbPb2760GeVAUnShifted_0010,"combined unshifted","lp");
            legendYdummy->AddEntry(graphCombInvYieldStatPbPb2760GeVYShifted_0010,"combined shifted","lp");
            legendYdummy->Draw();

            canvasDummy3->Update();
            canvasDummy3->Print(Form("%s/%s_ComparisonYShifted.%s",outputDir.Data(),meson.Data(),suffix.Data()));


            canvasDummy3->cd();
            histo2DDummy3->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeV_0010->Draw("pEsame");
            graphCombInvYieldStatPbPb2760GeV_2050->Draw("pEsame");

            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVYShifted_0010, 21, 1.5, kBlue, kBlue, widthLinesBoxes, kTRUE);
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVYShifted_2050, 21, 1.5, kBlue, kBlue, widthLinesBoxes, kTRUE);
            graphCombInvYieldStatPbPb2760GeVYShifted_0010->Draw("pEsame");
            graphCombInvYieldStatPbPb2760GeVYShifted_2050->Draw("pEsame");

            TLegend* legendXYdummy = new TLegend(0.6,0.8,0.9,0.95);
            legendXYdummy->SetFillColor(0);
            legendXYdummy->SetLineColor(0);
            legendXYdummy->SetTextFont(42);
            legendXYdummy->SetTextSize(FontSize);
            legendXYdummy->AddEntry(graphCombInvYieldStatPbPb2760GeV_0010,"x-shifted","lp");
            legendXYdummy->AddEntry(graphCombInvYieldStatPbPb2760GeVYShifted_0010,"y-shifted","lp");
            legendXYdummy->Draw();

            canvasDummy3->Update();
            canvasDummy3->Print(Form("%s/%s_ComparisonXYShifted.%s",outputDir.Data(),meson.Data(),suffix.Data()));

            delete canvasDummy3;
        }
    }

    Double_t limitPar = kTRUE;
    TString fitOpt = "QNRMEX0+";
    if(meson.CompareTo("Pi0")==0){

        //Bylinkin fit ---------------------------------------
        Double_t parampart1TCM_0010[2] = {graphCombInvYieldTotPbPb2760GeV_0010->GetY()[0],0.3};
        fitLowPtBylinkin_0010 = FitObject("Lowtcm","Lowtcm0010","Pi0",graphCombInvYieldTotPbPb2760GeV_0010,0.,6.,parampart1TCM_0010,fitOpt,"",limitPar);
        graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitLowPtBylinkin_0010,fitOpt,"",0.8,2.);
        Double_t parampart2TCM_0010[3] = {graphCombInvYieldTotPbPb2760GeV_0010->GetY()[3],0.3,8};
        cout << WriteParameterToFile(fitLowPtBylinkin_0010)<< endl << endl;
        fitHighPtBylinkin_0010 = FitObject("Hightcm","Hightcm0010","Pi0",graphCombInvYieldTotPbPb2760GeV_0010,0.,30.,parampart2TCM_0010,fitOpt,"",limitPar);
        graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitHighPtBylinkin_0010,fitOpt,"",1.,20.);
        cout << WriteParameterToFile(fitHighPtBylinkin_0010)<< endl << endl;

        Double_t paramTCM_0010[5] = {fitLowPtBylinkin_0010->GetParameter(0),fitLowPtBylinkin_0010->GetParameter(1),fitHighPtBylinkin_0010->GetParameter(0),fitHighPtBylinkin_0010->GetParameter(1),fitHighPtBylinkin_0010->GetParameter(2)};
//         Double_t paramTCM_0010[5] = {graphCombInvYieldTotPbPb2760GeV_0010->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeV_0010->GetY()[3],0.3,8};
        fitBylinkinPbPb2760GeVPtLHC11h_0010 = FitObject("tcm","BylinkinFitToYieldsPi00010","Pi0",graphCombInvYieldTotPbPb2760GeV_0010,0.9,21.,paramTCM_0010,fitOpt,"",limitPar);
        graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitBylinkinPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",minfitPt,20.);
        cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_0010)<< endl << endl;

        cout << __LINE__ << endl;
        mTScaledEtaFromPi0 = (TF1*)MtScaledParam(fitBylinkinPbPb2760GeVPtLHC11h_0010, 221, 0.48);

//         Double_t paramNormTCM_0010[5] = {graphCombInvYieldTotPbPb2760GeV_0010->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeV_0010->GetY()[3],0.3,8};
//         fitNormBylinkinPbPb2760GeVPtLHC11h_0010 = FitObject("ntcm","NormBylinkinFitToYieldsPi00010","Pi0",graphCombInvYieldTotPbPb2760GeV_0010,0.9,21.,paramTCM_0010,fitOpt,"",limitPar);
//         graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitNormBylinkinPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",minfitPt,20.);
//         cout << WriteParameterToFile(fitNormBylinkinPbPb2760GeVPtLHC11h_0010)<< endl << endl;

        Double_t parampart1TCM_2050[2] = {graphCombInvYieldTotPbPb2760GeV_2050->GetY()[0],0.3};
        fitLowPtBylinkin_2050 = FitObject("Lowtcm","Lowtcm2050","Pi0",graphCombInvYieldTotPbPb2760GeV_2050,0.,6.,parampart1TCM_2050,fitOpt,"",limitPar);
        graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitLowPtBylinkin_2050,fitOpt,"",0.8,2.);
//         cout << WriteParameterToFile(fitLowPtBylinkin_2050)<< endl << endl;
        Double_t parampart2TCM_2050[3] = {graphCombInvYieldTotPbPb2760GeV_2050->GetY()[3],0.3,8};
        fitHighPtBylinkin_2050 = FitObject("Hightcm","Hightcm2050","Pi0",graphCombInvYieldTotPbPb2760GeV_2050,0.,30.,parampart2TCM_2050,fitOpt,"",limitPar);
        graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitHighPtBylinkin_2050,fitOpt,"",1.,20.);
//         cout << WriteParameterToFile(fitHighPtBylinkin_2050)<< endl << endl;

        Double_t paramTCM_2050[5] = {fitLowPtBylinkin_2050->GetParameter(0),fitLowPtBylinkin_2050->GetParameter(1),fitHighPtBylinkin_2050->GetParameter(0),fitHighPtBylinkin_2050->GetParameter(1),fitHighPtBylinkin_2050->GetParameter(2)};
//         Double_t paramTCM_2050[5] = {graphCombInvYieldTotPbPb2760GeV_2050->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeV_2050->GetY()[3],0.3,8};
        fitBylinkinPbPb2760GeVPtLHC11h_2050 = FitObject("tcm","BylinkinFitToYieldsPi02050","Pi0",graphCombInvYieldTotPbPb2760GeV_2050,0.9,21.,paramTCM_2050,fitOpt,"",limitPar);
        graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitBylinkinPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",minfitPt,20.);
        cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_2050)<< endl << endl;

//         Double_t paramNormTCM_2050[5] = {graphCombInvYieldTotPbPb2760GeV_2050->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeV_2050->GetY()[3],0.3,8};
//         fitNormBylinkinPbPb2760GeVPtLHC11h_2050 = FitObject("ntcm","NormBylinkinFitToYieldsPi02050","Pi0",graphCombInvYieldTotPbPb2760GeV_2050,0.9,21.,paramTCM_2050,fitOpt,"",limitPar);
//         graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitNormBylinkinPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",minfitPt,20.);
//         cout << WriteParameterToFile(fitNormBylinkinPbPb2760GeVPtLHC11h_2050)<< endl << endl;

        //QCD fit ---------------------------------------
        fitQCDInvYieldPbPb2760GeV_0010 = FitObject("qcd","fitQCD","Pi0",graphCombInvYieldTotPbPb2760GeV_0010,0.9,21.,NULL,fitOpt);
        graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitQCDInvYieldPbPb2760GeV_0010,fitOpt,"",minfitPt,20.);
//         cout << WriteParameterToFile(fitQCDInvYieldPbPb2760GeV_0010)<< endl << endl;
        fitQCDInvYieldPbPb2760GeV_2050 = FitObject("qcd","fitQCD","Pi0",graphCombInvYieldTotPbPb2760GeV_2050,0.9,21.,NULL,fitOpt);
        graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitQCDInvYieldPbPb2760GeV_2050,fitOpt,"",minfitPt,20.);
//         cout << WriteParameterToFile(fitQCDInvYieldPbPb2760GeV_2050)<< endl << endl;


    }
    else if(meson.CompareTo("Eta")==0){

        //Bylinkin fit ---------------------------------------
        Double_t parampart1TCM_0010[2] = {graphCombInvYieldTotPbPb2760GeV_0010->GetY()[0],0.3};
        fitLowPtBylinkin_0010 = FitObject("Lowtcm","Lowtcm","Eta",graphCombInvYieldTotPbPb2760GeV_0010,0.,6.,parampart1TCM_0010,"NRMEX0+","",limitPar);
        graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitLowPtBylinkin_0010,fitOpt,"",0.8,2.);
        cout << WriteParameterToFile(fitLowPtBylinkin_0010)<< endl << endl;
        Double_t parampart2TCM_0010[3] = {graphCombInvYieldTotPbPb2760GeV_0010->GetY()[3],0.3,8};
        fitHighPtBylinkin_0010 = FitObject("Hightcm","Hightcm","Eta",graphCombInvYieldTotPbPb2760GeV_0010,0.,30.,parampart2TCM_0010,"NRMEX0+","",limitPar);
        graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitHighPtBylinkin_0010,fitOpt,"",1.,20.);
        cout << WriteParameterToFile(fitHighPtBylinkin_0010)<< endl << endl;

        Double_t paramTCM_0010[5] = {fitLowPtBylinkin_0010->GetParameter(0),fitLowPtBylinkin_0010->GetParameter(1),fitHighPtBylinkin_0010->GetParameter(0),fitHighPtBylinkin_0010->GetParameter(1),fitHighPtBylinkin_0010->GetParameter(2)};
//         Double_t paramTCM_0010[5] = {graphCombInvYieldTotPbPb2760GeV_0010->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeV_0010->GetY()[3],0.3,8};
        fitBylinkinPbPb2760GeVPtLHC11h_0010 = FitObject("tcm","BylinkinFitToYieldsEta0010","Eta",graphCombInvYieldTotPbPb2760GeV_0010,0.9,21.,paramTCM_0010,"NRMEX0+","",limitPar);
        graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitBylinkinPbPb2760GeVPtLHC11h_0010,fitOpt,"",1.,20.);
        cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_0010)<< endl << endl;

//         Double_t paramNormTCM_0010[5] = {graphCombInvYieldTotPbPb2760GeV_0010->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeV_0010->GetY()[3],0.3,8};
//         fitNormBylinkinPbPb2760GeVPtLHC11h_0010 = FitObject("ntcm","NormBylinkinFitToYieldsEta0010","Eta",graphCombInvYieldTotPbPb2760GeV_0010,0.9,21.,paramNormTCM_0010,fitOpt,"",limitPar);
//         graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitNormBylinkinPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",minfitPt,20.);
//         cout << WriteParameterToFile(fitNormBylinkinPbPb2760GeVPtLHC11h_0010)<< endl << endl;

        Double_t parampart1TCM_2050[2] = {graphCombInvYieldTotPbPb2760GeV_2050->GetY()[0],0.3};
        fitLowPtBylinkin_2050 = FitObject("Lowtcm","Lowtcm","Eta",graphCombInvYieldTotPbPb2760GeV_2050,0.,6.,parampart1TCM_2050,"NRMEX0+","",limitPar);
        graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitLowPtBylinkin_2050,fitOpt,"",0.8,2.);
        cout << WriteParameterToFile(fitLowPtBylinkin_2050)<< endl << endl;
        Double_t parampart2TCM_2050[3] = {graphCombInvYieldTotPbPb2760GeV_2050->GetY()[3],0.3,8};
        fitHighPtBylinkin_2050 = FitObject("Hightcm","Hightcm","Eta",graphCombInvYieldTotPbPb2760GeV_2050,0.,30.,parampart2TCM_2050,"NRMEX0+","",limitPar);
        graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitHighPtBylinkin_2050,fitOpt,"",1.,20.);
        cout << WriteParameterToFile(fitHighPtBylinkin_2050)<< endl << endl;

        Double_t paramTCM_2050[5] = {fitLowPtBylinkin_2050->GetParameter(0),fitLowPtBylinkin_2050->GetParameter(1),fitHighPtBylinkin_2050->GetParameter(0),fitHighPtBylinkin_2050->GetParameter(1),fitHighPtBylinkin_2050->GetParameter(2)};
//         Double_t paramTCM_2050[5] = {graphCombInvYieldTotPbPb2760GeV_2050->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeV_2050->GetY()[3],0.3,8};
        fitBylinkinPbPb2760GeVPtLHC11h_2050 = FitObject("tcm","BylinkinFitToYieldsEta2050","Eta",graphCombInvYieldTotPbPb2760GeV_2050,0.9,21.,paramTCM_2050,"NRMEX0+","",limitPar);
        graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitBylinkinPbPb2760GeVPtLHC11h_2050,fitOpt,"",1.,20.);
        cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_2050)<< endl << endl;

//         Double_t paramNormTCM_2050[5] = {graphCombInvYieldTotPbPb2760GeV_2050->GetY()[0],0.3,graphCombInvYieldTotPbPb2760GeV_2050->GetY()[3],0.3,8};
//         fitNormBylinkinPbPb2760GeVPtLHC11h_2050 = FitObject("ntcm","NormBylinkinFitToYieldsEta2050","Eta",graphCombInvYieldTotPbPb2760GeV_2050,0.9,21.,paramNormTCM_2050,fitOpt,"",limitPar);
//         graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitNormBylinkinPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",minfitPt,20.);
//         cout << WriteParameterToFile(fitNormBylinkinPbPb2760GeVPtLHC11h_2050)<< endl << endl;

        //QCD fit ---------------------------------------
        fitQCDInvYieldPbPb2760GeV_0010 = FitObject("qcd","fitQCD","Eta",graphCombInvYieldTotPbPb2760GeV_0010,0.9,21.,NULL,fitOpt);
        graphCombInvYieldTotPbPb2760GeV_0010->Fit(fitQCDInvYieldPbPb2760GeV_0010,fitOpt,"",1,20.);
//         cout << WriteParameterToFile(fitQCDInvYieldPbPb2760GeV_0010)<< endl << endl;
        fitQCDInvYieldPbPb2760GeV_2050 = FitObject("qcd","fitQCD","Eta",graphCombInvYieldTotPbPb2760GeV_2050,0.9,21.,NULL,fitOpt);
        graphCombInvYieldTotPbPb2760GeV_2050->Fit(fitQCDInvYieldPbPb2760GeV_2050,fitOpt,"",1,20.);
//         cout << WriteParameterToFile(fitQCDInvYieldPbPb2760GeV_2050)<< endl << endl;

    }

    // Ratio fit to data
    TGraphAsymmErrors* graphRatioCombCombFitTotPbPb2760GeV_0010   = (TGraphAsymmErrors*)graphCombInvYieldTotPbPb2760GeV_0010->Clone();
    graphRatioCombCombFitTotPbPb2760GeV_0010                      = CalculateGraphErrRatioToFit(graphRatioCombCombFitTotPbPb2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
    TGraphAsymmErrors* graphRatioCombCombFitStatPbPb2760GeV_0010  = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeV_0010->Clone();
    graphRatioCombCombFitStatPbPb2760GeV_0010                         = CalculateGraphErrRatioToFit(graphRatioCombCombFitStatPbPb2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
    TGraphAsymmErrors* graphRatioCombCombFitSysPbPb2760GeV_0010   = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeV_0010->Clone();
    graphRatioCombCombFitSysPbPb2760GeV_0010                      = CalculateGraphErrRatioToFit(graphRatioCombCombFitSysPbPb2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);

    TGraphAsymmErrors* graphRatioCombCombFitTotPbPb2760GeV_2050   = (TGraphAsymmErrors*)graphCombInvYieldTotPbPb2760GeV_2050->Clone();
    graphRatioCombCombFitTotPbPb2760GeV_2050                      = CalculateGraphErrRatioToFit(graphRatioCombCombFitTotPbPb2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);
    TGraphAsymmErrors* graphRatioCombCombFitStatPbPb2760GeV_2050  = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeV_2050->Clone();
    graphRatioCombCombFitStatPbPb2760GeV_2050                         = CalculateGraphErrRatioToFit(graphRatioCombCombFitStatPbPb2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);
    TGraphAsymmErrors* graphRatioCombCombFitSysPbPb2760GeV_2050   = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeV_2050->Clone();
    graphRatioCombCombFitSysPbPb2760GeV_2050                      = CalculateGraphErrRatioToFit(graphRatioCombCombFitSysPbPb2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);

    TGraphAsymmErrors* graphRatioCombQCDFitTotPbPb2760GeV_0010   = (TGraphAsymmErrors*)graphCombInvYieldTotPbPb2760GeV_0010->Clone();
    graphRatioCombQCDFitTotPbPb2760GeV_0010                      = CalculateGraphErrRatioToFit(graphRatioCombQCDFitTotPbPb2760GeV_0010, fitQCDInvYieldPbPb2760GeV_0010);
    TGraphAsymmErrors* graphRatioCombQCDFitStatPbPb2760GeV_0010  = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeV_0010->Clone();
    graphRatioCombQCDFitStatPbPb2760GeV_0010                         = CalculateGraphErrRatioToFit(graphRatioCombQCDFitStatPbPb2760GeV_0010, fitQCDInvYieldPbPb2760GeV_0010);
    TGraphAsymmErrors* graphRatioCombQCDFitSysPbPb2760GeV_0010   = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeV_0010->Clone();
    graphRatioCombQCDFitSysPbPb2760GeV_0010                      = CalculateGraphErrRatioToFit(graphRatioCombQCDFitSysPbPb2760GeV_0010, fitQCDInvYieldPbPb2760GeV_0010);

    TGraphAsymmErrors* graphRatioCombQCDFitTotPbPb2760GeV_2050   = (TGraphAsymmErrors*)graphCombInvYieldTotPbPb2760GeV_2050->Clone();
    graphRatioCombQCDFitTotPbPb2760GeV_2050                      = CalculateGraphErrRatioToFit(graphRatioCombQCDFitTotPbPb2760GeV_2050, fitQCDInvYieldPbPb2760GeV_2050);
    TGraphAsymmErrors* graphRatioCombQCDFitStatPbPb2760GeV_2050  = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeV_2050->Clone();
    graphRatioCombQCDFitStatPbPb2760GeV_2050                         = CalculateGraphErrRatioToFit(graphRatioCombQCDFitStatPbPb2760GeV_2050, fitQCDInvYieldPbPb2760GeV_2050);
    TGraphAsymmErrors* graphRatioCombQCDFitSysPbPb2760GeV_2050   = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeV_2050->Clone();
    graphRatioCombQCDFitSysPbPb2760GeV_2050                      = CalculateGraphErrRatioToFit(graphRatioCombQCDFitSysPbPb2760GeV_2050, fitQCDInvYieldPbPb2760GeV_2050);

    if(meson.CompareTo("Pi0")==0){

        graphRatioPCMCombFitStat2760GeV_0010        = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_0010->Clone();
        graphRatioPCMCombFitStat2760GeV_0010                        = CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
        graphRatioPCMCombFitSys2760GeV_0010         = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0010->Clone();
        graphRatioPCMCombFitSys2760GeV_0010                         = CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);

        graphRatioPHOSCombFitStat2760GeV_0010   = (TGraphAsymmErrors*)graphPHOSPi0InvYieldStatPbPb2760GeV_0010->Clone();
        graphRatioPHOSCombFitStat2760GeV_0010                       = CalculateGraphErrRatioToFit(graphRatioPHOSCombFitStat2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
        graphRatioPHOSCombFitSys2760GeV_0010    = (TGraphAsymmErrors*)graphPHOSPi0InvYieldSysPbPb2760GeV_0010->Clone();
        graphRatioPHOSCombFitSys2760GeV_0010                        = CalculateGraphErrRatioToFit(graphRatioPHOSCombFitSys2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);

        graphRatioEMCalCombFitStat2760GeV_0010  = (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Clone();
        graphRatioEMCalCombFitStat2760GeV_0010                      = CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
        graphRatioEMCalCombFitSys2760GeV_0010   = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Clone();
        graphRatioEMCalCombFitSys2760GeV_0010                       = CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);

        graphRatioPCMCombFitStat2760GeV_2050        = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_2050->Clone();
        graphRatioPCMCombFitStat2760GeV_2050                        = CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);
        graphRatioPCMCombFitSys2760GeV_2050         = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2050->Clone();
        graphRatioPCMCombFitSys2760GeV_2050                         = CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);

        graphRatioEMCalCombFitStat2760GeV_2050  = (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Clone();
        graphRatioEMCalCombFitStat2760GeV_2050                      = CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);
        graphRatioEMCalCombFitSys2760GeV_2050   = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Clone();
        graphRatioEMCalCombFitSys2760GeV_2050                       = CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);

    }
    else if(meson.CompareTo("Eta")==0){

        graphRatioPCMCombFitStat2760GeV_0010        = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_0010->Clone();
        graphRatioPCMCombFitStat2760GeV_0010                        = CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
        graphRatioPCMCombFitSys2760GeV_0010         = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0010->Clone();
        graphRatioPCMCombFitSys2760GeV_0010                         = CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);

        graphRatioEMCalCombFitStat2760GeV_0010  = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_0010->Clone();
        graphRatioEMCalCombFitStat2760GeV_0010                      = CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
        graphRatioEMCalCombFitSys2760GeV_0010   = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Clone();
        graphRatioEMCalCombFitSys2760GeV_0010                       = CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);

        graphRatioPCMCombFitStat2760GeV_2050        = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_2050->Clone();
        graphRatioPCMCombFitStat2760GeV_2050                        = CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);
        graphRatioPCMCombFitSys2760GeV_2050         = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2050->Clone();
        graphRatioPCMCombFitSys2760GeV_2050                         = CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);

        graphRatioEMCalCombFitStat2760GeV_2050  = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Clone();
        graphRatioEMCalCombFitStat2760GeV_2050                      = CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);
        graphRatioEMCalCombFitSys2760GeV_2050   = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Clone();
        graphRatioEMCalCombFitSys2760GeV_2050                       = CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);

    }

    textSizeLabelsPixel = 48;
    TCanvas* canvasRatioToCombFit = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFit, 0.12, 0.03, 0.01, 0.11);
    canvasRatioToCombFit->SetLogx();

        Double_t textsizeLabelsPP = 0;
        Double_t textsizeFacPP= 0;
        if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
            textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
            textsizeFacPP = (Double_t)1./canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
        } else {
            textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
            textsizeFacPP = (Double_t)1./canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
        }
        cout << textsizeLabelsPP << endl;

        Double_t push;
        if(thesisPlotting || PaperPi0) push = 0.;
        else push = 0.25;

        TH2F * histo2DRatioToCombFit = new TH2F("histo2DRatioToCombFit","histo2DRatioToCombFit",1000,0.25,30.,1000,0.2,4.  );
        SetStyleHistoTH2ForGraphs(histo2DRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb fit} ", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                                0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
        histo2DRatioToCombFit->GetXaxis()->SetMoreLogLabels();
        histo2DRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);

        histo2DRatioToCombFit->GetXaxis()->SetRangeUser(0.5-push,30.);
        histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
        histo2DRatioToCombFit->Draw("copy");

        DrawGammaLines(0.5-push, 30. , 1., 1.,0.5,   kGray+2);
        DrawGammaLines(0.5-push, 30. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.5-push, 30. , 0.9, 0.9,0.5, kGray, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysPbPb2760GeV_0010, markerStyleComb, markerSizeComb, colorCombo0010 , colorCombo0010, widthLinesBoxes, kTRUE);
        graphRatioCombCombFitSysPbPb2760GeV_0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatPbPb2760GeV_0010, markerStyleComb, markerSizeComb, colorCombo0010 , colorCombo0010);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitStatPbPb2760GeV_0010);
        graphRatioCombCombFitStatPbPb2760GeV_0010->Draw("p,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysPbPb2760GeV_2050, markerStyleComb, markerSizeComb, colorCombo2050, colorCombo2050, widthLinesBoxes, kTRUE);
        graphRatioCombCombFitSysPbPb2760GeV_2050->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatPbPb2760GeV_2050, markerStyleComb, markerSizeComb, colorCombo2050, colorCombo2050);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitStatPbPb2760GeV_2050);
        graphRatioCombCombFitStatPbPb2760GeV_2050->Draw("p,same");

        TLegend* legendCombToCombFit = new TLegend(0.2,0.75,0.5,0.94);
        legendCombToCombFit->SetFillColor(0);
        legendCombToCombFit->SetLineColor(0);
        legendCombToCombFit->SetTextFont(42);
        legendCombToCombFit->SetTextSize(0.04);
        legendCombToCombFit->AddEntry(graphRatioCombCombFitSysPbPb2760GeV_0010,"0#font[122]{-}10%","fp");
        legendCombToCombFit->AddEntry(graphRatioCombCombFitSysPbPb2760GeV_2050,"20#font[122]{-}50%","fp");
        legendCombToCombFit->Draw();

        TLatex *labelRatioToFitEnergy = new TLatex(0.65,0.92,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioToFitEnergy, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitEnergy->SetTextFont(43);
        if(meson.CompareTo("Pi0")==0){
            labelRatioToFitPi0= new TLatex(0.73,0.87,"#pi^{0} #rightarrow #gamma#gamma");
        } else if(meson.CompareTo("Eta")==0){
            labelRatioToFitPi0= new TLatex(0.73,0.87,"#eta #rightarrow #gamma#gamma");
        }
        SetStyleTLatex( labelRatioToFitPi0, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitPi0->SetTextFont(43);
        labelRatioToFitEnergy->Draw();
        labelRatioToFitPi0->Draw();
        histo2DRatioToCombFit->Draw("axis,same");

    canvasRatioToCombFit->SaveAs(Form("%s/%s_RatioOfCombToCombFit_PbPbPbPb2760GeV.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/%s_RatioOfCombToCombFit_PbPbPbPb2760GeV.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

        histo2DRatioToCombFit->Draw("copy");

        DrawGammaLines(0.5-push, 30. , 1., 1.,0.5,   kGray+2);
        DrawGammaLines(0.5-push, 30. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.5-push, 30. , 0.9, 0.9,0.5, kGray, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombQCDFitSysPbPb2760GeV_0010, markerStyleComb, markerSizeComb, colorCombo0010 , colorCombo0010, widthLinesBoxes, kTRUE);
        graphRatioCombQCDFitSysPbPb2760GeV_0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioCombQCDFitStatPbPb2760GeV_0010, markerStyleComb, markerSizeComb, colorCombo0010 , colorCombo0010);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioCombQCDFitStatPbPb2760GeV_0010);
        graphRatioCombQCDFitStatPbPb2760GeV_0010->Draw("p,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombQCDFitSysPbPb2760GeV_2050, markerStyleComb, markerSizeComb, colorCombo2050, colorCombo2050, widthLinesBoxes, kTRUE);
        graphRatioCombQCDFitSysPbPb2760GeV_2050->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioCombQCDFitStatPbPb2760GeV_2050, markerStyleComb, markerSizeComb, colorCombo2050, colorCombo2050);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioCombQCDFitStatPbPb2760GeV_2050);
        graphRatioCombQCDFitStatPbPb2760GeV_2050->Draw("p,same");

        legendCombToCombFit->Draw();
        labelRatioToFitPi0->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/%s_RatioOfCombToQCDFit_PbPbPbPb2760GeV.%s",outputDir.Data(),meson.Data(),suffix.Data()));

    canvasRatioToCombFit->cd();
    histo2DRatioToCombFit->Draw("copy");

        DrawGammaLines(0.5-push, 30. , 1., 1.,0.5,   kGray+2);
        DrawGammaLines(0.5-push, 30. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.5-push, 30. , 0.9, 0.9,0.5, kGray, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitSys2760GeV_0010, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);//colorPCM0010, colorPCM0010, widthLinesBoxes, kTRUE);
        graphRatioPCMCombFitSys2760GeV_0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitStat2760GeV_0010, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);//colorPCM0010, colorPCM0010);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioPCMCombFitStat2760GeV_0010);
        graphRatioPCMCombFitStat2760GeV_0010->Draw("p,same");

        if(meson.CompareTo("Pi0")==0){
            DrawGammaSetMarkerTGraphAsym(graphRatioPHOSCombFitSys2760GeV_0010, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
            graphRatioPHOSCombFitSys2760GeV_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphRatioPHOSCombFitStat2760GeV_0010, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioPHOSCombFitStat2760GeV_0010);
            graphRatioPHOSCombFitStat2760GeV_0010->Draw("p,same");
        }

        DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitSys2760GeV_0010, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);//colorEMCal0010, colorEMCal0010, widthLinesBoxes, kTRUE);
        graphRatioEMCalCombFitSys2760GeV_0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitStat2760GeV_0010, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);//colorEMCal0010, colorEMCal0010);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioEMCalCombFitStat2760GeV_0010);
        graphRatioEMCalCombFitStat2760GeV_0010->Draw("p,same");

        // ****************************** Definition of the Legend ******************************************
        // **************** Row def ************************
//         Double_t rowsLegendOnlyPi0Ratio[5]      = {0.92,0.88,0.84,0.80,0.76};
//         Double_t rowsLegendOnlyPi0RatioAbs[5]   = {0.91,2.2,2.1,2.0,1.9};
//         Double_t columnsLegendOnlyPi0Ratio[3]   = {0.15,0.25, 0.305};
//         Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.15,1.04, 1.37};
//         Double_t lengthBox                      = 0.2/2;
//         Double_t heightBox                      = 0.08/2;
        //for thesis:
        Double_t rowsLegendOnlyPi0Ratio[5]      = {0.92,0.88,0.84,0.80,0.76};
        Double_t rowsLegendOnlyPi0RatioAbs[5]   = {0.91,2.2,2.1,2.0,1.9};
        Double_t columnsLegendOnlyPi0Ratio[3]   = {0.17,0.38, 0.43};
        Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.17,1.3+0.17+0.38, 1.7+0.17+0.43};
        Double_t lengthBox                      = 0.2/2;
        Double_t heightBox                      = 0.08/2;
        if(meson.CompareTo("Eta")==0){
          lengthBox                      = 0.3/2;
          columnsLegendOnlyPi0RatioAbs[0]= 0.15;
          columnsLegendOnlyPi0RatioAbs[1] = 2.;
          columnsLegendOnlyPi0RatioAbs[2] = 2.5;
        }
        // ****************** first Column **************************************************
        TLatex *textPCMOnlyRatioPi0LHC11h = NULL;
        if(thesisPlotting) textPCMOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],"PCM (this thesis)");
        else textPCMOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],"PCM");
        SetStyleTLatex( textPCMOnlyRatioPi0LHC11h, 0.85*textSizeLabelsPixel,4);
        textPCMOnlyRatioPi0LHC11h->SetTextFont(43);
        textPCMOnlyRatioPi0LHC11h->Draw();
        TLatex *textEMCalOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[2],"EMCal");//" (*)");
        SetStyleTLatex( textEMCalOnlyRatioPi0LHC11h,  0.85*textSizeLabelsPixel,4);
        textEMCalOnlyRatioPi0LHC11h->SetTextFont(43);
        textEMCalOnlyRatioPi0LHC11h->Draw();
        TLatex *textPHOSOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],"PHOS");
        SetStyleTLatex( textPHOSOnlyRatioPi0LHC11h,  0.85*textSizeLabelsPixel,4);
        textPHOSOnlyRatioPi0LHC11h->SetTextFont(43);
        if(meson.CompareTo("Pi0")==0){
          textPHOSOnlyRatioPi0LHC11h->Draw();
        }
        // ****************** second Column *************************************************
        TLatex *textStatOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioPi0LHC11h, 0.85*textSizeLabelsPixel,4);
        textStatOnlyRatioPi0LHC11h->SetTextFont(43);
        textStatOnlyRatioPi0LHC11h->Draw();
        TLatex *textSysOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioPi0LHC11h, 0.85*textSizeLabelsPixel,4);
        textSysOnlyRatioPi0LHC11h->SetTextFont(43);
        textSysOnlyRatioPi0LHC11h->Draw();


        TMarker* markerPCMPi0OnlyRatioPi0LHC11h = CreateMarkerFromGraph(graphRatioPCMCombFitSys2760GeV_0010,columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
        markerPCMPi0OnlyRatioPi0LHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
        TBox* boxPCMPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioPCMCombFitSys2760GeV_0010, columnsLegendOnlyPi0RatioAbs[2]-1.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 2*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
        boxPCMPi0OnlyRatioPi0->Draw("l");
        TMarker* markerEMCalPi0OnlyRatioPi0LHC11h = CreateMarkerFromGraph(graphRatioEMCalCombFitSys2760GeV_0010, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[2],1);
        markerEMCalPi0OnlyRatioPi0LHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[2]);
        TBox* boxEMCalPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioEMCalCombFitSys2760GeV_0010, columnsLegendOnlyPi0RatioAbs[2]-1.5*lengthBox , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 2*lengthBox, rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
        boxEMCalPi0OnlyRatioPi0->Draw("l");
        if(meson.CompareTo("Pi0")==0){
            TMarker* markerPHOSPi0OnlyRatioPi0LHC11h = CreateMarkerFromGraph(graphRatioPHOSCombFitSys2760GeV_0010, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
            markerPHOSPi0OnlyRatioPi0LHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);
            TBox* boxPHOSPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioPHOSCombFitSys2760GeV_0010, columnsLegendOnlyPi0RatioAbs[2]-1.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
                                                            columnsLegendOnlyPi0RatioAbs[2]+ 2*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
            boxPHOSPi0OnlyRatioPi0->Draw("l");
        }

        TLatex *labelRatioToFitEnergy0010 = new TLatex(0.55,0.92,collisionSystemPbPb0010.Data());
        SetStyleTLatex( labelRatioToFitEnergy0010, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitEnergy0010->SetTextFont(43);
        labelRatioToFitEnergy0010->Draw();
        TLatex *labelRatioToFitPi0sepMeas = new TLatex(0.5,0.89,"");
        if(meson.CompareTo("Pi0")==0){
            labelRatioToFitPi0sepMeas= new TLatex(0.5,0.89,"#pi^{0} #rightarrow #gamma#gamma");
        } else if(meson.CompareTo("Eta")==0){
            labelRatioToFitPi0sepMeas= new TLatex(0.5,0.89,"#eta #rightarrow #gamma#gamma");
        }
        SetStyleTLatex( labelRatioToFitPi0sepMeas, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitPi0sepMeas->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/%s_RatioOfIndividualMeasToCombFitLHC11h_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/%s_RatioOfIndividualMeasToCombFitLHC11h_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    canvasRatioToCombFit->cd();
    histo2DRatioToCombFit->Draw("copy");

        DrawGammaLines(0.5-push, 30. , 1., 1.,0.5,   kGray+2);
        DrawGammaLines(0.5-push, 30. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.5-push, 30. , 0.9, 0.9,0.5, kGray, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitSys2760GeV_2050, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);//colorPCM2050, colorPCM2050, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitSys2760GeV_2050, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);//colorEMCal2050, colorEMCal2050, widthLinesBoxes, kTRUE);

        DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitStat2760GeV_2050, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);//colorPCM2050, colorPCM2050);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioPCMCombFitStat2760GeV_2050);
        DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitStat2760GeV_2050, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);//colorEMCal2050, colorEMCal2050);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioEMCalCombFitStat2760GeV_2050);

        graphRatioPCMCombFitSys2760GeV_2050->Draw("E2same");
        graphRatioEMCalCombFitSys2760GeV_2050->Draw("E2same");
        graphRatioPCMCombFitStat2760GeV_2050->Draw("p,same");
        graphRatioEMCalCombFitStat2760GeV_2050->Draw("p,same");

        textPCMOnlyRatioPi0LHC11h->Draw();
        textEMCalOnlyRatioPi0LHC11h->Draw();

        textStatOnlyRatioPi0LHC11h->Draw();
        textSysOnlyRatioPi0LHC11h->Draw();

        TMarker* markerPCMPi0OnlyRatioPi0LHC11h_2050 = CreateMarkerFromGraph(graphRatioPCMCombFitSys2760GeV_2050,columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
        markerPCMPi0OnlyRatioPi0LHC11h_2050->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
        TBox* boxPCMPi0OnlyRatioPi0_2050 = CreateBoxFromGraph(graphRatioPCMCombFitSys2760GeV_2050, columnsLegendOnlyPi0RatioAbs[2]-1.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 2*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
        boxPCMPi0OnlyRatioPi0_2050->Draw("l");

        TMarker* markerEMCalPi0OnlyRatioPi0LHC11h_2050 = CreateMarkerFromGraph(graphRatioEMCalCombFitSys2760GeV_2050, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
        markerEMCalPi0OnlyRatioPi0LHC11h_2050->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);
        TBox* boxEMCalPi0OnlyRatioPi0_2050 = CreateBoxFromGraph(graphRatioEMCalCombFitSys2760GeV_2050, columnsLegendOnlyPi0RatioAbs[2]-1.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 2*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
        boxEMCalPi0OnlyRatioPi0_2050->Draw("l");

        TLatex *labelRatioToFitEnergy2050 = new TLatex(0.5,0.92,collisionSystemPbPb2050.Data());
        SetStyleTLatex( labelRatioToFitEnergy2050, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitEnergy2050->SetTextFont(43);
        labelRatioToFitEnergy2050->Draw();
        labelRatioToFitPi0sepMeas->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/%s_RatioOfIndividualMeasToCombFitLHC11h_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/%s_RatioOfIndividualMeasToCombFitLHC11h_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    //4pads yields ratio to model
    Double_t arrayXMeas4padsratios[2];
    Double_t arrayYMeas4padsratios[3];
    Double_t relXMeas4padsratios[3];
    Double_t relYMeas4padsratios[3];
    Int_t textSizeLabelsPixel4pads = 70;
    Double_t textsizeLabelsXSecDown = 0;
    Double_t textsizeFacXSecDown = 0;
    Double_t textsizeLabelsXSecMiddle = 0;
    Double_t textsizeFacXSecMiddle = 0;
    ReturnCorrectValuesForCanvasScaling(2500, 2000, 1, 2, 0.1, 0.025, 0.025, 0.08, arrayXMeas4padsratios, arrayYMeas4padsratios, relXMeas4padsratios, relYMeas4padsratios);
    TCanvas* canvasInvYieldMeasRatios4pads = new TCanvas("canvasInvYieldMeasRatios4pads","",0,0,2500,2000);
        TPad* padInvYieldMeasLowerRatio = new TPad("","",arrayXMeas4padsratios[0], arrayYMeas4padsratios[2], arrayXMeas4padsratios[1], arrayYMeas4padsratios[1],-1, -1, -2);
        DrawGammaPadSettings( padInvYieldMeasLowerRatio, relXMeas4padsratios[0], relXMeas4padsratios[2], relYMeas4padsratios[1], relYMeas4padsratios[2]);
        TPad* padInvYieldMeasUpperRatio = new TPad("","",arrayXMeas4padsratios[0], arrayYMeas4padsratios[1], arrayXMeas4padsratios[1], arrayYMeas4padsratios[0],-1, -1, -2);
        DrawGammaPadSettings( padInvYieldMeasUpperRatio, relXMeas4padsratios[0], relXMeas4padsratios[2], relYMeas4padsratios[0], relYMeas4padsratios[1]);
        Double_t marginXRatio = relXMeas4padsratios[0]*2500;
        padInvYieldMeasLowerRatio->Draw();
        if (padInvYieldMeasLowerRatio->XtoPixel(padInvYieldMeasLowerRatio->GetX2()) < padInvYieldMeasLowerRatio->YtoPixel(padInvYieldMeasLowerRatio->GetY1())){
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel4pads/padInvYieldMeasLowerRatio->XtoPixel(padInvYieldMeasLowerRatio->GetX2()) ;
            textsizeFacXSecDown = (Double_t)1./padInvYieldMeasLowerRatio->XtoPixel(padInvYieldMeasLowerRatio->GetX2()) ;
        } else {
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel4pads/padInvYieldMeasLowerRatio->YtoPixel(padInvYieldMeasLowerRatio->GetY1());
            textsizeFacXSecDown = (Double_t)1./padInvYieldMeasLowerRatio->YtoPixel(padInvYieldMeasLowerRatio->GetY1());
        }
        padInvYieldMeasUpperRatio->Draw();
        if (padInvYieldMeasUpperRatio->XtoPixel(padInvYieldMeasUpperRatio->GetX2()) < padInvYieldMeasUpperRatio->YtoPixel(padInvYieldMeasUpperRatio->GetY1())){
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel4pads/padInvYieldMeasUpperRatio->XtoPixel(padInvYieldMeasUpperRatio->GetX2()) ;
            textsizeFacXSecMiddle = (Double_t)1./padInvYieldMeasUpperRatio->XtoPixel(padInvYieldMeasUpperRatio->GetX2()) ;
        } else {
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel4pads/padInvYieldMeasUpperRatio->YtoPixel(padInvYieldMeasUpperRatio->GetY1());
            textsizeFacXSecMiddle = (Double_t)1./padInvYieldMeasUpperRatio->YtoPixel(padInvYieldMeasUpperRatio->GetY1());
        }

        TH2F * ratio2DUpperPadMeasRatio = new TH2F("ratio2DUpperPadMeasRatio","ratio2DUpperPadMeasRatio",1000,0.25,70.,1000,0.01,2.25);
        SetStyleHistoTH2ForGraphs(ratio2DUpperPadMeasRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb. fit} ", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                    0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.15/(textsizeFacXSecMiddle*marginXRatio), 510, 505);
        if(thesisPlotting) ratio2DUpperPadMeasRatio->GetXaxis()->SetRangeUser(0.6,30.);
        else ratio2DUpperPadMeasRatio->GetXaxis()->SetRangeUser(0.5,30);
        ratio2DUpperPadMeasRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DUpperPadMeasRatio->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DUpperPadMeasRatio->GetYaxis()->SetLabelOffset(+0.01);
        TH2F * ratio2DLowerPadMeasRatio =  new TH2F("ratio2DLowerPadMeasRatio","ratio2DLowerPadMeasRatio",1000,0.25,70.,1000,0.01,2.25);
        SetStyleHistoTH2ForGraphs(ratio2DLowerPadMeasRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb. fit} ", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                    0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 1,0.15/(textsizeFacXSecDown*marginXRatio), 510, 505);
        if(thesisPlotting) ratio2DLowerPadMeasRatio->GetXaxis()->SetRangeUser(0.6,30);
        else ratio2DLowerPadMeasRatio->GetXaxis()->SetRangeUser(0.5,30);
        ratio2DLowerPadMeasRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DLowerPadMeasRatio->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DLowerPadMeasRatio->GetYaxis()->SetLabelOffset(+0.01);

        Double_t rowsLegendMeasRatio[5]      = {0.87,0.82,0.76,0.7,0.68};
        Double_t rowsLegendMeasRatioAbs[5]   = {0.88,1.95,1.83,1.72,1.9};
        Double_t columnsLegendMeasRatio[3]   = {0.15,0.33, 0.38};
        Double_t columnsLegendMeasRatioAbs[3]= {0.15,1.77, 2.3};
        Double_t lengthMeasBox                      = 0.2/2;
        Double_t heightMeasBox                      = 0.08/2;
        if(meson.CompareTo("Eta")==0){
          lengthMeasBox                      = 0.2/2;
          columnsLegendMeasRatioAbs[0]= 0.15;
          columnsLegendMeasRatioAbs[1] = 1.77;
          columnsLegendMeasRatioAbs[2] = 2.3;
        }
        // ****************** first Column **************************************************
        TLatex *textPCMMeasRatioLHC11h = NULL;
        if(thesisPlotting) textPCMMeasRatioLHC11h = new TLatex(columnsLegendMeasRatio[0],rowsLegendMeasRatio[1],"PCM (this thesis)");
        else textPCMMeasRatioLHC11h = new TLatex(columnsLegendMeasRatio[0],rowsLegendMeasRatio[1],"PCM");
        SetStyleTLatex( textPCMMeasRatioLHC11h, 0.85*textSizeLabelsPixel4pads,4);
        textPCMMeasRatioLHC11h->SetTextFont(43);
        TLatex *textEMCalMeasRatioLHC11h = new TLatex(columnsLegendMeasRatio[0],rowsLegendMeasRatio[2],"EMCal"); //" (*)");
        SetStyleTLatex( textEMCalMeasRatioLHC11h,  0.85*textSizeLabelsPixel4pads,4);
        textEMCalMeasRatioLHC11h->SetTextFont(43);
        TLatex *textPHOSMeasRatioLHC11h = new TLatex(columnsLegendMeasRatio[0],rowsLegendMeasRatio[3],"PHOS");
        SetStyleTLatex( textPHOSMeasRatioLHC11h,  0.85*textSizeLabelsPixel4pads,4);
        textPHOSMeasRatioLHC11h->SetTextFont(43);
        // ****************** second Column *************************************************
        TLatex *textStatMeasRatioLHC11h = new TLatex(columnsLegendMeasRatio[1],rowsLegendMeasRatio[0] ,"stat");
        SetStyleTLatex( textStatMeasRatioLHC11h, 0.85*textSizeLabelsPixel4pads,4);
        textStatMeasRatioLHC11h->SetTextFont(43);
        TLatex *textSysMeasRatioLHC11h = new TLatex(columnsLegendMeasRatio[2] ,rowsLegendMeasRatio[0],"syst");
        SetStyleTLatex( textSysMeasRatioLHC11h, 0.85*textSizeLabelsPixel4pads,4);
        textSysMeasRatioLHC11h->SetTextFont(43);


        padInvYieldMeasUpperRatio->cd();
        padInvYieldMeasUpperRatio->SetLogx(1);
        ratio2DUpperPadMeasRatio->DrawCopy();
            if(thesisPlotting)  DrawGammaLines(0.6,30,1., 1.,0.1,kGray);
            else DrawGammaLines(0.4,30,1., 1.,0.1,kGray);

            graphRatioPCMCombFitSys2760GeV_0010->Draw("E2same");
            graphRatioPCMCombFitStat2760GeV_0010->Draw("p,same");
            if(meson.CompareTo("Pi0")==0){
                graphRatioPHOSCombFitSys2760GeV_0010->Draw("E2same");
                graphRatioPHOSCombFitStat2760GeV_0010->Draw("p,same");
            }
            graphRatioEMCalCombFitSys2760GeV_0010->Draw("E2same");
            graphRatioEMCalCombFitStat2760GeV_0010->Draw("p,same");

            textStatMeasRatioLHC11h->Draw();
            textSysMeasRatioLHC11h->Draw();

            TMarker* markerPCMPi0MeasRatioLHC11h = CreateMarkerFromGraph(graphRatioPCMCombFitSys2760GeV_0010,columnsLegendMeasRatio[1] ,rowsLegendMeasRatio[1],1);
            markerPCMPi0MeasRatioLHC11h->DrawMarker(columnsLegendMeasRatioAbs[1] ,rowsLegendMeasRatioAbs[1]);
            TBox* boxPCMPi0MeasRatio = CreateBoxFromGraph(graphRatioPCMCombFitSys2760GeV_0010, columnsLegendMeasRatioAbs[2]-1.5*lengthMeasBox , rowsLegendMeasRatioAbs[1]- heightMeasBox,
                                                            columnsLegendMeasRatioAbs[2]+ 2*lengthMeasBox, rowsLegendMeasRatioAbs[1]+ heightMeasBox);
            TMarker* markerEMCalPi0MeasRatioLHC11h = CreateMarkerFromGraph(graphRatioEMCalCombFitSys2760GeV_0010, columnsLegendMeasRatio[1] ,rowsLegendMeasRatio[2],1);
            markerEMCalPi0MeasRatioLHC11h->DrawMarker(columnsLegendMeasRatioAbs[1] ,rowsLegendMeasRatioAbs[2]);
            TBox* boxEMCalPi0MeasRatio = CreateBoxFromGraph(graphRatioEMCalCombFitSys2760GeV_0010, columnsLegendMeasRatioAbs[2]-1.5*lengthMeasBox , rowsLegendMeasRatioAbs[2]- heightMeasBox,
                                                            columnsLegendMeasRatioAbs[2]+ 2*lengthMeasBox, rowsLegendMeasRatioAbs[2]+ heightMeasBox);
            textPCMMeasRatioLHC11h->Draw();
            textEMCalMeasRatioLHC11h->Draw();
            boxPCMPi0MeasRatio->Draw("l");
            boxEMCalPi0MeasRatio->Draw("l");
            if(meson.CompareTo("Pi0")==0){
                textPHOSMeasRatioLHC11h->Draw();

                TMarker* markerPHOSPi0MeasRatioLHC11h = CreateMarkerFromGraph(graphRatioPHOSCombFitSys2760GeV_0010, columnsLegendMeasRatio[1] ,rowsLegendMeasRatio[3],1);
                markerPHOSPi0MeasRatioLHC11h->DrawMarker(columnsLegendMeasRatioAbs[1] ,rowsLegendMeasRatioAbs[3]);
                TBox* boxPHOSPi0MeasRatio = CreateBoxFromGraph(graphRatioPHOSCombFitSys2760GeV_0010, columnsLegendMeasRatioAbs[2]-1.5*lengthMeasBox , rowsLegendMeasRatioAbs[3]- heightMeasBox,
                                                                columnsLegendMeasRatioAbs[2]+ 2*lengthMeasBox, rowsLegendMeasRatioAbs[3]+ heightMeasBox);
                boxPHOSPi0MeasRatio->Draw("l");
            }
//             labelRatioToFitPi0sepMeas->Draw();

            TLatex *labelMeasEnergy0010 = new TLatex(0.55,0.85,collisionSystemPbPb0010.Data());
            SetStyleTLatex( labelMeasEnergy0010, 0.85*textSizeLabelsPixel4pads,4);
            labelMeasEnergy0010->SetTextFont(43);
            labelMeasEnergy0010->Draw();

        ratio2DUpperPadMeasRatio->Draw("axis,same");
        padInvYieldMeasLowerRatio->cd();
        padInvYieldMeasLowerRatio->SetLogx(1);
        ratio2DLowerPadMeasRatio->DrawCopy();
            if(thesisPlotting)  DrawGammaLines(0.6,30,1., 1.,0.1,kGray);
            else DrawGammaLines(0.4,30,1., 1.,0.1,kGray);

            TLatex *textEMCalNoteRatioToFit = new TLatex(0.15,0.88,"(*) ALICE preliminary input");
            SetStyleTLatex( textEMCalNoteRatioToFit,  0.8*textSizeLabelsPixel4pads,4);
            textEMCalNoteRatioToFit->SetTextFont(43);
            TLatex *textEMCalNoteRatioToFit1 = new TLatex(0.15,0.83,"arXiv:1609.06106");
            SetStyleTLatex( textEMCalNoteRatioToFit1,  0.8*textSizeLabelsPixel4pads,4);
            textEMCalNoteRatioToFit1->SetTextFont(43);
            TLatex *textEMCalNoteRatioToFit2 = new TLatex(0.15,0.78,"analysis by A. Morreale");
            SetStyleTLatex( textEMCalNoteRatioToFit2,  0.8*textSizeLabelsPixel4pads,4);
            textEMCalNoteRatioToFit2->SetTextFont(43);
//             textEMCalNoteRatioToFit->Draw();
//             textEMCalNoteRatioToFit1->Draw();
//             textEMCalNoteRatioToFit2->Draw();

            graphRatioPCMCombFitSys2760GeV_2050->Draw("E2same");
            graphRatioEMCalCombFitSys2760GeV_2050->Draw("E2same");
            graphRatioPCMCombFitStat2760GeV_2050->Draw("p,same");
            graphRatioEMCalCombFitStat2760GeV_2050->Draw("p,same");

            TLatex *labelMeasEnergy2050 = new TLatex(0.55,0.9,collisionSystemPbPb2050.Data());
            SetStyleTLatex( labelMeasEnergy2050, 0.85*textSizeLabelsPixel4pads,4);
            labelMeasEnergy2050->SetTextFont(43);
            labelMeasEnergy2050->Draw();

        ratio2DLowerPadMeasRatio->Draw("axis,same");
    canvasInvYieldMeasRatios4pads->Update();
    canvasInvYieldMeasRatios4pads->SaveAs(Form("%s/%s_RatioOfIndividualMeasToCombFitLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasInvYieldMeasRatios4pads->SaveAs(Form("%s/%s_RatioOfIndividualMeasToCombFitLHC11h.%s",paperPlots.Data(),meson.Data(),suffix.Data()));

    //**********************************************************************************************************************//
    //***************************************** Plotting Combined Invariant Yields *****************************************//
    //**********************************************************************************************************************//
    // yields plotting

    TCanvas* canvasInvariantYield = new TCanvas("canvasInvariantYield","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasInvariantYield, 0.16, 0.02, 0.02, 0.09);
    canvasInvariantYield->SetLogx();
    canvasInvariantYield->SetLogy();

        if(meson.CompareTo("Pi0")==0){
            histo2DInvariantYield = new TH2F("histo2DInvariantYield","histo2DInvariantYield",11000,0.25,70.,1000,2e-8,1e4);
            SetStyleHistoTH2ForGraphs(histo2DInvariantYield, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.5);
        } else if(meson.CompareTo("Eta")==0){
            histo2DInvariantYield = new TH2F("histo2DInvariantYield","histo2DInvariantYield",11000,0.25,70.,1000,2e-9,1e3);
            SetStyleHistoTH2ForGraphs(histo2DInvariantYield, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
        }
        histo2DInvariantYield->GetXaxis()->SetMoreLogLabels();
        histo2DInvariantYield->GetXaxis()->SetLabelOffset(-0.01);
        histo2DInvariantYield->Draw("copy");

        if(meson.CompareTo("Pi0")==0){
            DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldSysPbPb2760GeV_0010, markerStyleDet[0], markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
            graphPCMPi0InvYieldSysPbPb2760GeV_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldStatPbPb2760GeV_0010,markerStyleDet[0], markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphPCMPi0InvYieldStatPbPb2760GeV_0010);
            graphPCMPi0InvYieldStatPbPb2760GeV_0010->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphPHOSPi0InvYieldSysPbPb2760GeV_0010, markerStyleDet[1], markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
            graphPHOSPi0InvYieldSysPbPb2760GeV_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphPHOSPi0InvYieldStatPbPb2760GeV_0010, markerStyleDet[1], markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphPHOSPi0InvYieldStatPbPb2760GeV_0010);
            graphPHOSPi0InvYieldStatPbPb2760GeV_0010->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphEMCalPi0InvYieldSysPbPb2760GeV_0010, markerStyleDet[2], markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
            graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphEMCalPi0InvYieldStatPbPb2760GeV_0010, markerStyleDet[2], markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphEMCalPi0InvYieldStatPbPb2760GeV_0010);
            graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Draw("p,same");

        }
        else if(meson.CompareTo("Eta")==0){

            DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldSysPbPb2760GeV_0010, markerStyleDet[0], markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);            graphPCMEtaInvYieldSysPbPb2760GeV_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldStatPbPb2760GeV_0010,markerStyleDet[0], markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphPCMEtaInvYieldStatPbPb2760GeV_0010);
            graphPCMEtaInvYieldStatPbPb2760GeV_0010->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphEMCalEtaInvYieldSysPbPb2760GeV_0010, markerStyleDet[2], markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
            graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphEMCalEtaInvYieldStatPbPb2760GeV_0010, markerStyleDet[2], markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphEMCalEtaInvYieldStatPbPb2760GeV_0010);
            graphEMCalEtaInvYieldStatPbPb2760GeV_0010->Draw("p,same");

        }
        Float_t xspace;
        if(meson.CompareTo("Pi0")==0) xspace = 0.1;
        if(meson.CompareTo("Eta")==0) xspace = 0.13;
        TLegend* legendInvYieldSectionPi0LHC11h_separate = new TLegend(0.2,0.15,0.5,xspace*nMeasSetALHC11h0010);
        legendInvYieldSectionPi0LHC11h_separate->SetFillColor(0);
        legendInvYieldSectionPi0LHC11h_separate->SetLineColor(0);
        legendInvYieldSectionPi0LHC11h_separate->SetTextFont(42);
        legendInvYieldSectionPi0LHC11h_separate->SetTextSize(FontSize);
        if(meson.CompareTo("Pi0")==0){
            legendInvYieldSectionPi0LHC11h_separate->AddEntry(graphPCMPi0InvYieldSysPbPb2760GeV_0010,"PCM","fp");
            legendInvYieldSectionPi0LHC11h_separate->AddEntry(graphPHOSPi0InvYieldSysPbPb2760GeV_0010,"PHOS","fp");
            legendInvYieldSectionPi0LHC11h_separate->AddEntry(graphEMCalPi0InvYieldSysPbPb2760GeV_0010,"EMCal","fp");
        } else if(meson.CompareTo("Eta")==0){
            legendInvYieldSectionPi0LHC11h_separate->AddEntry(graphPCMEtaInvYieldSysPbPb2760GeV_0010,"PCM","fp");
            legendInvYieldSectionPi0LHC11h_separate->AddEntry(graphEMCalEtaInvYieldSysPbPb2760GeV_0010,"EMCal","fp");
        }
        legendInvYieldSectionPi0LHC11h_separate->Draw();
        TLatex *labelYieldsEnergysepMeas0010 = new TLatex(0.45,0.92,collisionSystemPbPb0010.Data());
        SetStyleTLatex( labelYieldsEnergysepMeas0010,0.035,4);
        labelYieldsEnergysepMeas0010->Draw();
        TLatex *labelyieldsMesonsepMeas = NULL;
        if(meson.CompareTo("Pi0")==0){
            labelyieldsMesonsepMeas= new TLatex(0.45,0.89,"#pi^{0} #rightarrow #gamma#gamma");
        } else if(meson.CompareTo("Eta")==0){
            labelyieldsMesonsepMeas= new TLatex(0.45,0.89,"#eta #rightarrow #gamma#gamma");
        }
        SetStyleTLatex( labelyieldsMesonsepMeas, 0.035,4);
        labelyieldsMesonsepMeas->Draw();

    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCompareAllSystemsLHC11h_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCompareAllSystemsLHC11h_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasInvariantYield->cd();
        histo2DInvariantYield->Draw("copy");

        if(meson.CompareTo("Pi0")==0){
            DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldSysPbPb2760GeV_2050, markerStyleDet[0], markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
            graphPCMPi0InvYieldSysPbPb2760GeV_2050->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldStatPbPb2760GeV_2050,markerStyleDet[0], markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphPCMPi0InvYieldStatPbPb2760GeV_2050);
            graphPCMPi0InvYieldStatPbPb2760GeV_2050->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphEMCalPi0InvYieldSysPbPb2760GeV_2050, markerStyleDet[2], markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
            graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphEMCalPi0InvYieldStatPbPb2760GeV_2050, markerStyleDet[2], markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphEMCalPi0InvYieldStatPbPb2760GeV_2050);
            graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Draw("p,same");
        }
        else if(meson.CompareTo("Eta")==0){
            DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldSysPbPb2760GeV_2050, markerStyleDet[0], markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
            graphPCMEtaInvYieldSysPbPb2760GeV_2050->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldStatPbPb2760GeV_2050,markerStyleDet[0], markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphPCMEtaInvYieldStatPbPb2760GeV_2050);
            graphPCMEtaInvYieldStatPbPb2760GeV_2050->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphEMCalEtaInvYieldSysPbPb2760GeV_2050, markerStyleDet[2], markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
            graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphEMCalEtaInvYieldStatPbPb2760GeV_2050, markerStyleDet[2], markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphEMCalEtaInvYieldStatPbPb2760GeV_2050);
            graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Draw("p,same");
        }

        TLegend* legendInvYieldSectionPi0LHC11h_separate2 = new TLegend(0.2,0.15,0.5,0.13*nMeasSetALHC11h2050);
        legendInvYieldSectionPi0LHC11h_separate2->SetFillColor(0);
        legendInvYieldSectionPi0LHC11h_separate2->SetLineColor(0);
        legendInvYieldSectionPi0LHC11h_separate2->SetTextFont(42);
        legendInvYieldSectionPi0LHC11h_separate2->SetTextSize(FontSize);
        if(meson.CompareTo("Pi0")==0){
          legendInvYieldSectionPi0LHC11h_separate2->AddEntry(graphPCMPi0InvYieldSysPbPb2760GeV_2050,"PCM","fp");
          legendInvYieldSectionPi0LHC11h_separate2->AddEntry(graphEMCalPi0InvYieldSysPbPb2760GeV_2050,"EMCal","fp");
        } else if(meson.CompareTo("Eta")==0){
          legendInvYieldSectionPi0LHC11h_separate2->AddEntry(graphPCMEtaInvYieldSysPbPb2760GeV_2050,"PCM","fp");
          legendInvYieldSectionPi0LHC11h_separate2->AddEntry(graphEMCalEtaInvYieldSysPbPb2760GeV_2050,"EMCal","fp");
        }
        legendInvYieldSectionPi0LHC11h_separate2->Draw();
        TLatex *labelYieldsEnergysepMeas2050 = new TLatex(0.45,0.92,collisionSystemPbPb2050.Data());
        SetStyleTLatex( labelYieldsEnergysepMeas2050,0.035,4);
        labelYieldsEnergysepMeas2050->Draw();
        labelyieldsMesonsepMeas->Draw();

    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCompareAllSystemsLHC11h_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCompareAllSystemsLHC11h_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasInvariantYield->cd();
        histo2DInvariantYield->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldSysPbPb2760GeV_0010, markerStyleComb, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_0010, markerStyleComb, markerSizeComb, colorCombo0010, colorCombo0010);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStatPbPb2760GeV_0010);

        //scaling of yields for plotting sake: ////////////////////////////
        TGraphAsymmErrors *graphCombInvYieldSysPbPb2760GeVScaled_0010 = (TGraphAsymmErrors*)graphCombInvYieldSysPbPb2760GeV_0010->Clone("graphCombInvYieldSysPbPb2760GeVScaled_0010");
        graphCombInvYieldSysPbPb2760GeVScaled_0010 = ScaleGraph(graphCombInvYieldSysPbPb2760GeVScaled_0010, 4);
        TGraphAsymmErrors *graphCombInvYieldStatPbPb2760GeVScaled_0010 = (TGraphAsymmErrors*)graphCombInvYieldStatPbPb2760GeV_0010->Clone("graphCombInvYieldStatPbPb2760GeVScaled_0010");
        graphCombInvYieldStatPbPb2760GeVScaled_0010 = ScaleGraph(graphCombInvYieldStatPbPb2760GeVScaled_0010,4);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStatPbPb2760GeVScaled_0010);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStatPbPb2760GeV_2050);
        //////////////////////////////////////////////////////////////////

        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldSysPbPb2760GeVScaled_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeVScaled_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010);
        graphCombInvYieldSysPbPb2760GeVScaled_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeVScaled_0010->Draw("p,same");

        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldSysPbPb2760GeV_2050, markerStyle2050, markerSizeComb, colorCombo2050 ,colorCombo2050, 2, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStatPbPb2760GeV_2050, markerStyle2050, markerSizeComb,  colorCombo2050 ,colorCombo2050);
        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

        TLegend* legendInvYieldSectionPi0LHC11h_onlyPbPb = new TLegend(0.2,0.13,0.53,0.23);
        legendInvYieldSectionPi0LHC11h_onlyPbPb->SetMargin(0.17);
        legendInvYieldSectionPi0LHC11h_onlyPbPb->SetFillColor(0);
        legendInvYieldSectionPi0LHC11h_onlyPbPb->SetLineColor(0);
        legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextFont(42);
        legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextSize(FontSize);
        legendInvYieldSectionPi0LHC11h_onlyPbPb->SetHeader(collisionSystem2760GeV.Data());
        legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphCombInvYieldSysPbPb2760GeVScaled_0010,Form("  %s",cent0010.Data()),"fp");
        legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphCombInvYieldSysPbPb2760GeV_2050,Form("%s",cent2050.Data()),"fp");
        legendInvYieldSectionPi0LHC11h_onlyPbPb->Draw();

        labelSystOnlyPbPb->Draw();
        labelFactorLowerOnlyPbPb->Draw();

    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    canvasInvariantYield->cd();
        histo2DInvariantYield->Draw("copy");

        graphCombInvYieldSysPbPb2760GeVScaled_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeVScaled_0010->Draw("p,same");

        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

        TF1* fitTCMInvYieldPbPb2760GeVScaledforPlot_0010 = (TF1*)fitBylinkinPbPb2760GeVPtLHC11h_0010->Clone();
        fitTCMInvYieldPbPb2760GeVScaledforPlot_0010->SetLineStyle(4);
        fitTCMInvYieldPbPb2760GeVScaledforPlot_0010->SetParameter(0,fitTCMInvYieldPbPb2760GeVScaledforPlot_0010->GetParameter(0)*4);
        fitTCMInvYieldPbPb2760GeVScaledforPlot_0010->SetParameter(2,fitTCMInvYieldPbPb2760GeVScaledforPlot_0010->GetParameter(2)*4);
        fitBylinkinPbPb2760GeVPtLHC11h_2050->SetLineStyle(4);
        fitTCMInvYieldPbPb2760GeVScaledforPlot_0010->Draw("same");
        fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");

        TF1* fitQCDInvYieldPbPb2760GeVScaledforPlot_0010 = (TF1*)fitQCDInvYieldPbPb2760GeV_0010->Clone();
        fitQCDInvYieldPbPb2760GeVScaledforPlot_0010->SetLineStyle(2);
        fitQCDInvYieldPbPb2760GeVScaledforPlot_0010->SetLineColor(2);
        fitQCDInvYieldPbPb2760GeVScaledforPlot_0010->SetParameter(0,fitQCDInvYieldPbPb2760GeVScaledforPlot_0010->GetParameter(0)*4);
//         fitQCDInvYieldPbPb2760GeVScaledforPlot_0010->Draw("same");
        fitQCDInvYieldPbPb2760GeV_2050->SetLineStyle(2);
        fitQCDInvYieldPbPb2760GeV_2050->SetLineColor(2);
//         fitQCDInvYieldPbPb2760GeV_2050->Draw("same");

        TF1* fitTCMpart1ScaledforPlot_0010 = (TF1*)fitLowPtBylinkin_0010->Clone();
        TF1* fitTCMpart2ScaledforPlot_0010 = (TF1*)fitHighPtBylinkin_0010->Clone();
        fitTCMpart1ScaledforPlot_0010->SetLineStyle(2);
        fitTCMpart1ScaledforPlot_0010->SetLineColor(2);
        fitTCMpart1ScaledforPlot_0010->SetParameter(0,fitTCMpart1ScaledforPlot_0010->GetParameter(0)*4);
        fitTCMpart2ScaledforPlot_0010->SetLineStyle(4);
        fitTCMpart2ScaledforPlot_0010->SetLineColor(4);
        fitTCMpart2ScaledforPlot_0010->SetParameter(0,fitTCMpart2ScaledforPlot_0010->GetParameter(0)*4);
//         fitTCMpart1ScaledforPlot_0010->Draw("same");
//         fitTCMpart2ScaledforPlot_0010->Draw("same");

        TF1* fitTCMpart1ScaledforPlot_2050 = (TF1*)fitLowPtBylinkin_2050->Clone();
        TF1* fitTCMpart2ScaledforPlot_2050 = (TF1*)fitHighPtBylinkin_2050->Clone();
        fitTCMpart1ScaledforPlot_2050->SetLineStyle(2);
        fitTCMpart1ScaledforPlot_2050->SetLineColor(2);
        fitTCMpart2ScaledforPlot_2050->SetLineStyle(4);
        fitTCMpart2ScaledforPlot_2050->SetLineColor(4);
//         fitTCMpart1ScaledforPlot_2050->Draw("same");
//         fitTCMpart2ScaledforPlot_2050->Draw("same");

        TLegend* legendInvYieldSectionPi0LHC11h_WithFit= new TLegend(0.595,0.77,0.91,0.915); //0.17,0.13,0.5,0.24);
        legendInvYieldSectionPi0LHC11h_WithFit->SetFillColor(0);
        legendInvYieldSectionPi0LHC11h_WithFit->SetMargin(0.17);
        legendInvYieldSectionPi0LHC11h_WithFit->SetLineColor(0);
        legendInvYieldSectionPi0LHC11h_WithFit->SetTextFont(42);
        legendInvYieldSectionPi0LHC11h_WithFit->SetTextSize(FontSize);
        legendInvYieldSectionPi0LHC11h_WithFit->SetHeader(collisionSystem2760GeV.Data());
        legendInvYieldSectionPi0LHC11h_WithFit->AddEntry(graphCombInvYieldSysPbPb2760GeVScaled_0010,Form("  %s",cent0010.Data()),"fp");
        legendInvYieldSectionPi0LHC11h_WithFit->AddEntry(graphCombInvYieldSysPbPb2760GeV_2050,Form("%s",cent2050.Data()),"fp");
        legendInvYieldSectionPi0LHC11h_WithFit->AddEntry(fitTCMInvYieldPbPb2760GeVScaledforPlot_0010,"fits to Pb#font[122]{-}Pb","l");
        legendInvYieldSectionPi0LHC11h_WithFit->Draw();
        labelSyst->Draw();
        labelFactorLowerOnlyPbPb->Draw();
        if(thesisPlotting){
          thesisLabelHighRight->Draw();
        }

    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnlyWithFit.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnlyWithFit.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    canvasInvariantYield->cd();
      if(meson.CompareTo("Pi0")==0){
          histo2DInvariantYieldwithPP = new TH2F("histo2DInvariantYieldwithPP","histo2DInvariantYieldwithPP",11000,0.25,70.,1000,5e-12,3e3);
          SetStyleHistoTH2ForGraphs(histo2DInvariantYieldwithPP, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2} ",0.035,0.04, 0.035,0.04, 1.,1.65);

      } else if(meson.CompareTo("Eta")==0){
          histo2DInvariantYieldwithPP = new TH2F("histo2DInvariantYieldwithPP","histo2DInvariantYieldwithPP",11000,0.25,70.,1000,5e-10,3e2);
          SetStyleHistoTH2ForGraphs(histo2DInvariantYieldwithPP, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2} ",0.035,0.04, 0.035,0.04, 1.,1.6);
      }
      histo2DInvariantYieldwithPP->GetXaxis()->SetMoreLogLabels();
      histo2DInvariantYieldwithPP->GetXaxis()->SetLabelOffset(-0.01);
      histo2DInvariantYieldwithPP->Draw("copy");

        TLegend* legendSpectraPP = new TLegend(0.22,0.15,0.56,0.22);
        legendSpectraPP->SetFillColor(0);
        legendSpectraPP->SetLineColor(0);
        legendSpectraPP->SetTextFont(42);
        legendSpectraPP->SetMargin(0.17);
        legendSpectraPP->SetTextSize(FontSize);
        if(meson.CompareTo("Pi0")==0){
            graphInvSectionCombStatPi02760GeVPlot->SetMarkerSize(3.5);
            graphInvSectionCombStatPi02760GeVPlot->Draw("p,same");
            graphInvSectionCombSysPi02760GeVPlot->Draw("E2same");
            legendSpectraPP->SetHeader(collisionSystemPP2760GeV.Data());
            legendSpectraPP->AddEntry(graphInvSectionCombSysPi02760GeVPlot,"EPJC 77 (2017) 339"/*"arXiv: 1702.00917"*/,"pf");
        } else if(meson.CompareTo("Eta")==0){
            graphInvSectionCombStatEta2760GeVPlot->SetMarkerSize(3.5);
            graphInvSectionCombStatEta2760GeVPlot->Draw("p,same");
            graphInvSectionCombSysEta2760GeVPlot->Draw("E2same");
            legendSpectraPP->SetHeader(collisionSystemPP2760GeV.Data());
            legendSpectraPP->AddEntry(graphInvSectionCombSysEta2760GeVPlot,"EPJC 77 (2017) 339"/*"arXiv: 1702.00917"*/,"pf");
        }
        legendSpectraPP->Draw();

        fitTCMInvYieldPbPb2760GeVScaledforPlot_0010->Draw("same");
        fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");

        TLegend* legendInvYieldSectionPi0LHC11h_WitPP = new TLegend(0.595,0.92-(0.035*5-0.03),0.595+0.33,0.92); //0.17,0.13,0.5,0.24);
        legendInvYieldSectionPi0LHC11h_WitPP->SetFillColor(0);
        legendInvYieldSectionPi0LHC11h_WitPP->SetMargin(0.17);
        legendInvYieldSectionPi0LHC11h_WitPP->SetLineColor(0);
        legendInvYieldSectionPi0LHC11h_WitPP->SetTextFont(42);
        legendInvYieldSectionPi0LHC11h_WitPP->SetTextSize(FontSize);
        legendInvYieldSectionPi0LHC11h_WitPP->SetHeader(collisionSystem2760GeV.Data());
        legendInvYieldSectionPi0LHC11h_WitPP->AddEntry(graphCombInvYieldSysPbPb2760GeVScaled_0010,Form("  %s",cent0010.Data()),"fp");
        legendInvYieldSectionPi0LHC11h_WitPP->AddEntry(graphCombInvYieldSysPbPb2760GeV_2050,Form("%s",cent2050.Data()),"fp");
        legendInvYieldSectionPi0LHC11h_WitPP->AddEntry(fitTCMInvYieldPbPb2760GeVScaledforPlot_0010,"fits to Pb#font[122]{-}Pb","l");
        legendInvYieldSectionPi0LHC11h_WitPP->Draw();
        labelSyst->Draw();
        labelFactorLower->Draw();

        graphCombInvYieldStatPbPb2760GeV_2050->SetMarkerSize(3.);
        graphCombInvYieldStatPbPb2760GeVScaled_0010->SetMarkerSize(3.);

        graphCombInvYieldSysPbPb2760GeVScaled_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeVScaled_0010->Draw("p,same");

        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnlyWithPP.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnlyWithPP.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnlyWithPP.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasInvariantYield->cd();
      histo2DInvariantYield->Draw("copy");
//       histo2DInvariantYield->GetYaxis()->SetRangeUser(1e-9,);

        fitBylinkinPbPb2760GeVPtLHC11h_0010->SetLineStyle(1);
        fitBylinkinPbPb2760GeVPtLHC11h_0010->Draw("same");

        graphCombInvYieldSysPbPb2760GeV_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_0010->Draw("p,same");

        TF1 *fitLowPtpart_0010 = new TF1("fitLowPtpart_0010",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass));
        fitLowPtpart_0010->SetRange(0., 30.);
        fitLowPtpart_0010->FixParameter(0,fitBylinkinPbPb2760GeVPtLHC11h_0010->GetParameter(0));
        fitLowPtpart_0010->FixParameter(1,fitBylinkinPbPb2760GeVPtLHC11h_0010->GetParameter(1));
        fitLowPtpart_0010->SetLineStyle(2);
        fitLowPtpart_0010->SetLineColor(kRed);
//         fitLowPtpart_0010->Draw("same");

        TF1 *fitHighPtpart_0010 = new TF1("fitHighPtpart_0010","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]) )");
        fitHighPtpart_0010->SetRange(0., 30.);
        fitHighPtpart_0010->FixParameter(0,fitBylinkinPbPb2760GeVPtLHC11h_0010->GetParameter(2));
        fitHighPtpart_0010->FixParameter(1,fitBylinkinPbPb2760GeVPtLHC11h_0010->GetParameter(3));
        fitHighPtpart_0010->FixParameter(2,fitBylinkinPbPb2760GeVPtLHC11h_0010->GetParameter(4));
        fitHighPtpart_0010->SetLineStyle(4);
        fitHighPtpart_0010->SetLineColor(kBlue);
//         fitHighPtpart_0010->Draw("same");

        fitLowPtBylinkin_0010->SetLineStyle(2);
        fitLowPtBylinkin_0010->SetLineColor(kMagenta+1);
        fitLowPtBylinkin_0010->Draw("same");
        fitHighPtBylinkin_0010->SetLineStyle(4);
        fitHighPtBylinkin_0010->SetLineColor(kGreen+2);
        fitHighPtBylinkin_0010->Draw("same");

        TLegend* legendInvYieldSingle = new TLegend(0.2,0.2,0.5,0.35);
        legendInvYieldSingle->SetFillColor(0);
        legendInvYieldSingle->SetMargin(0.15);
        legendInvYieldSingle->SetLineColor(0);
        legendInvYieldSingle->SetTextFont(42);
        legendInvYieldSingle->SetTextSize(FontSize);
        legendInvYieldSingle->SetHeader(collisionSystem2760GeV.Data());
        legendInvYieldSingle->AddEntry(fitBylinkinPbPb2760GeVPtLHC11h_0010,"Bylinkin fit to Pb#font[122]{-}Pb","l");
//         legendInvYieldSingle->AddEntry(fitLowPtpart_0010,"Low p_{T} Bylinkin (drawn)","l");
//         legendInvYieldSingle->AddEntry(fitHighPtpart_0010,"High p_{T} Bylinkin (drawn)","l");
        legendInvYieldSingle->AddEntry(fitLowPtBylinkin_0010,"Low p_{T} Bylinkin (fitted)","l");
        legendInvYieldSingle->AddEntry(fitHighPtBylinkin_0010,"High p_{T} Bylinkin (fitted)","l");
        legendInvYieldSingle->Draw();

        TLatex *labelBylinkin = new TLatex(0.22,0.14,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}");
        SetStyleTLatex( labelBylinkin, 0.03,4);
        labelBylinkin->Draw();
        labelSystOnlyPbPb->Draw();

//         fitNormBylinkinPbPb2760GeVPtLHC11h_0010->SetLineStyle(4);
//         fitNormBylinkinPbPb2760GeVPtLHC11h_0010->SetLineColor(6);
//         fitNormBylinkinPbPb2760GeVPtLHC11h_0010->Draw("same");

    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    canvasInvariantYield->cd();
        histo2DInvariantYieldwithPP->DrawCopy();

        graphCombInvYieldSysPbPb2760GeV_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_0010->Draw("p,same");

        TLegend* legendInvYieldPi0= new TLegend(0.2,0.3,0.5,0.37);
        legendInvYieldPi0->SetFillColor(0);
        legendInvYieldPi0->SetMargin(0.17);
        legendInvYieldPi0->SetLineColor(0);
        legendInvYieldPi0->SetTextFont(42);
        legendInvYieldPi0->SetTextSize(FontSize);
        legendInvYieldPi0->SetHeader("ALICE 0#font[122]{-}10% Pb#font[122]{-}Pb");
        legendInvYieldPi0->AddEntry(graphCombInvYieldSysPbPb2760GeV_0010,Form("  %s",cent0010.Data()),"fp");
        legendInvYieldPi0->Draw();

        TLegend* legendInvYieldwithPHENIX = new TLegend(0.2,0.12,0.5,0.3);
        legendInvYieldwithPHENIX->SetFillColor(0);
        legendInvYieldwithPHENIX->SetMargin(0.15);
        legendInvYieldwithPHENIX->SetLineColor(0);
        legendInvYieldwithPHENIX->SetTextFont(42);
        legendInvYieldwithPHENIX->SetTextSize(FontSize);

        if(meson.CompareTo("Pi0")==0){

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVInvYield_0010, markerStylePHENIX200GeV, markerSizePHENIX200GeV, kBlack , kBlack);
            DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVInvYield_0010, markerStylePHENIX62GeV, markerSizePHENIX62GeV, kBlack , kBlack);
            DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVInvYield_0010, markerStylePHENIX39GeV, markerSizePHENIX39GeV, kBlack, kBlack);
            graphPHENIX200GeVInvYield_0010->Draw("p,same");
            graphPHENIX62GeVInvYield_0010->Draw("p,same");
            graphPHENIX39GeVInvYield_0010->Draw("p,same");

            legendInvYieldwithPHENIX->SetHeader("PHENIX 0#font[122]{-}10% Au#font[122]{-}Au");
            legendInvYieldwithPHENIX->AddEntry(graphPHENIX200GeVInvYield_0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendInvYieldwithPHENIX->AddEntry(graphPHENIX62GeVInvYield_0010,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
            legendInvYieldwithPHENIX->AddEntry(graphPHENIX39GeVInvYield_0010,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");

        } else {

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaInvYield_0020, markerStylePHENIX39GeV, markerSizePHENIX39GeV, kBlack, kBlack);
            graphPHENIX200GeVEtaInvYield_0020->Draw("p,same");

            legendInvYieldwithPHENIX->SetHeader("PHENIX 0#font[122]{-}20% Au#font[122]{-}Au");
            legendInvYieldwithPHENIX->AddEntry(graphPHENIX200GeVEtaInvYield_0020,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");

        }
        legendInvYieldwithPHENIX->Draw();

        labelSystOnlyPbPb->Draw();
        if(thesisPlotting)thesisLabelHighRight->Draw();

    canvasInvariantYield->Update();
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithPhenix_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithPhenix_0010.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithPhenix_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasInvariantYield->cd();
        histo2DInvariantYieldwithPP->DrawCopy();

        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

        TLegend* legendInvYieldPi02 = new TLegend(0.2,0.3,0.5,0.37);
        legendInvYieldPi02->SetFillColor(0);
        legendInvYieldPi02->SetMargin(0.17);
        legendInvYieldPi02->SetLineColor(0);
        legendInvYieldPi02->SetTextFont(42);
        legendInvYieldPi02->SetTextSize(FontSize);
        legendInvYieldPi02->SetHeader("ALICE 20#font[122]{-}50% Pb#font[122]{-}Pb");
        legendInvYieldPi02->AddEntry(graphCombInvYieldSysPbPb2760GeV_2050,Form("  %s",cent0010.Data()),"fp");
        legendInvYieldPi02->Draw();

        TLegend* legendInvYieldwithPHENIX2 = new TLegend(0.2,0.12,0.5,0.3);
        legendInvYieldwithPHENIX2->SetFillColor(0);
        legendInvYieldwithPHENIX2->SetMargin(0.15);
        legendInvYieldwithPHENIX2->SetLineColor(0);
        legendInvYieldwithPHENIX2->SetTextFont(42);
        legendInvYieldwithPHENIX2->SetTextSize(FontSize);

        if(meson.CompareTo("Pi0")==0){

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVInvYield_2040, markerStylePHENIX200GeV, markerSizePHENIX200GeV, kBlack , kBlack);
            DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVInvYield_2040, markerStylePHENIX62GeV, markerSizePHENIX62GeV, kBlack , kBlack);
            DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVInvYield_2040, markerStylePHENIX39GeV, markerSizePHENIX39GeV, kBlack, kBlack);
            graphPHENIX200GeVInvYield_2040->Draw("p,same");
            graphPHENIX62GeVInvYield_2040->Draw("p,same");
            graphPHENIX39GeVInvYield_2040->Draw("p,same");

            legendInvYieldwithPHENIX2->SetHeader("PHENIX 20#font[122]{-}40% Au#font[122]{-}Au");
            legendInvYieldwithPHENIX2->AddEntry(graphPHENIX200GeVInvYield_2040,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendInvYieldwithPHENIX2->AddEntry(graphPHENIX62GeVInvYield_2040,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
            legendInvYieldwithPHENIX2->AddEntry(graphPHENIX39GeVInvYield_2040,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");

        }
        else {

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaInvYield_2060, markerStylePHENIX39GeV, markerSizePHENIX39GeV, kBlack, kBlack);
            graphPHENIX200GeVEtaInvYield_2060->Draw("p,same");

            legendInvYieldwithPHENIX2->SetHeader("PHENIX 20#font[122]{-}60% Au#font[122]{-}Au");
            legendInvYieldwithPHENIX2->AddEntry(graphPHENIX200GeVEtaInvYield_2060,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");

        }
        legendInvYieldwithPHENIX2->Draw();
        labelSystOnlyPbPb->Draw();
        if(thesisPlotting)thesisLabelHighRight->Draw();

    canvasInvariantYield->Update();
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithPhenix_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithPhenix_2050.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithPhenix_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasInvariantYield->cd();
      histo2DInvariantYield->Draw("copy");

        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

        fitBylinkinPbPb2760GeVPtLHC11h_2050->SetLineStyle(1);
        fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");


        TF1 *fitLowPtpart_2050 = new TF1("fitLowPtpart_2050",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass));
        fitLowPtpart_2050->SetRange(0., 30.);
        fitLowPtpart_2050->FixParameter(0,fitBylinkinPbPb2760GeVPtLHC11h_2050->GetParameter(0));
        fitLowPtpart_2050->FixParameter(1,fitBylinkinPbPb2760GeVPtLHC11h_2050->GetParameter(1));
        fitLowPtpart_2050->SetLineStyle(2);
        fitLowPtpart_2050->SetLineColor(kRed);
//         fitLowPtpart_2050->Draw("same");

        TF1 *fitHighPtpart_2050 = new TF1("fitHighPtpart_2050","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]) )");
        fitHighPtpart_2050->SetRange(0., 30.);
        fitHighPtpart_2050->FixParameter(0,fitBylinkinPbPb2760GeVPtLHC11h_2050->GetParameter(2));
        fitHighPtpart_2050->FixParameter(1,fitBylinkinPbPb2760GeVPtLHC11h_2050->GetParameter(3));
        fitHighPtpart_2050->FixParameter(2,fitBylinkinPbPb2760GeVPtLHC11h_2050->GetParameter(4));
        fitHighPtpart_2050->SetLineStyle(4);
        fitHighPtpart_2050->SetLineColor(kBlue);
//         fitHighPtpart_2050->Draw("same");

        fitLowPtBylinkin_2050->SetLineStyle(2);
        fitLowPtBylinkin_2050->SetLineColor(kMagenta+1);
        fitLowPtBylinkin_2050->Draw("same");
        fitHighPtBylinkin_2050->SetLineStyle(4);
        fitHighPtBylinkin_2050->SetLineColor(kGreen+2);
        fitHighPtBylinkin_2050->Draw("same");

        legendInvYieldSingle->Draw();
        labelBylinkin->Draw();
        labelSystOnlyPbPb->Draw();

//         fitNormBylinkinPbPb2760GeVPtLHC11h_2050->SetLineStyle(4);
//         fitNormBylinkinPbPb2760GeVPtLHC11h_2050->SetLineColor(6);
//         fitNormBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");

    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    canvasInvariantYield->cd();

        TH2F * histo2DInvariantYieldWithModels = new TH2F("histo2DInvariantYieldWithModels","histo2DInvariantYieldWithModels",11000,0.25,70.,1000,7e-9,1e6);
        SetStyleHistoTH2ForGraphs(histo2DInvariantYieldWithModels, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}, #eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2} ",0.035,0.04, 0.035,0.04, 1.,1.6);
        if(thesisPlotting){
          histo2DInvariantYieldWithModels->GetXaxis()->SetRangeUser(0.5,40.);
          histo2DInvariantYieldWithModels->GetYaxis()->SetRangeUser(7e-9,3e7);
        } else {
          histo2DInvariantYieldWithModels->GetXaxis()->SetRangeUser(0.5,50.);
        }
        histo2DInvariantYieldWithModels->DrawCopy();

        //all mesons together
        TFile *MesonInput = new TFile(Form("%s/CombinedResultsPaperPbPb2760GeV_%s.root", outputDir.Data(),dateForOutput.Data()));
        if(MesonInput){

            TGraphAsymmErrors *graphInvYieldPi0StatBothMeson_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0CombPbPb2760GeVStatErr_0010");
            TGraphAsymmErrors *graphInvYieldPi0SysBothMeson_0010= (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0CombPbPb2760GeVSysErr_0010");
            TGraphAsymmErrors *graphInvYieldPi0StatBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0CombPbPb2760GeVStatErr_2050");
            TGraphAsymmErrors *graphInvYieldPi0SysBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0CombPbPb2760GeVSysErr_2050");

            TGraphAsymmErrors *graphInvYieldEtaStatBothMeson_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaCombPbPb2760GeVStatErr_0010");
            TGraphAsymmErrors *graphInvYieldEtaSysBothMeson_0010= (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaCombPbPb2760GeVSysErr_0010");
            TGraphAsymmErrors *graphInvYieldEtaStatBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaCombPbPb2760GeVStatErr_2050");
            TGraphAsymmErrors *graphInvYieldEtaSysBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaCombPbPb2760GeVSysErr_2050");

            TGraphErrors *graphEPOSPi0forPlot_0010 = (TGraphErrors*)MesonInput->Get("EPOSpredictionPi0_0010");
            TGraphErrors *graphEPOSPi0forPlot_2050 = (TGraphErrors*)MesonInput->Get("EPOSpredictionPi0_2050");
            TGraphErrors *graphEPOSEtaforPlot_0010 = (TGraphErrors*)MesonInput->Get("EPOSpredictionEta_0010");
            TGraphErrors *graphEPOSEtaforPlot_2050 = (TGraphErrors*)MesonInput->Get("EPOSpredictionEta_2050");

            if(graphInvYieldPi0StatBothMeson_0010 && graphInvYieldEtaStatBothMeson_0010){

                TGraphAsymmErrors *graphNLOPi0ScaledandNColl_0010 = (TGraphAsymmErrors*)graphNLOCalcDSS14InvYieldPi0Band->Clone("graphNLOPi0ScaledandNColl_0010");
                graphNLOPi0ScaledandNColl_0010 = ScaleGraph(graphNLOPi0ScaledandNColl_0010,nColl0010*4e2);
                TGraphAsymmErrors *graphNLOPi0ScaledandNColl_2050 = (TGraphAsymmErrors*)graphNLOCalcDSS14InvYieldPi0Band->Clone("graphNLOPi0ScaledandNColl_2050");
                graphNLOPi0ScaledandNColl_2050 = ScaleGraph(graphNLOPi0ScaledandNColl_2050,nColl2050*1e2);

                TGraphAsymmErrors *graphNLOEtaonlyNColl_0010 = (TGraphAsymmErrors*)graphNLOCalcDSS07InvYieldEtaBand->Clone("graphNLOEtaonlyNColl_0010");
                graphNLOEtaonlyNColl_0010 = ScaleGraph(graphNLOEtaonlyNColl_0010,nColl0010);
                TGraphAsymmErrors *graphNLOEtaScaledandNColl_2050 = (TGraphAsymmErrors*)graphNLOCalcDSS07InvYieldEtaBand->Clone("graphNLOEtaScaledandNColl_2050");
                graphNLOEtaScaledandNColl_2050 = ScaleGraph(graphNLOEtaScaledandNColl_2050,nColl2050*2e-1);


                for(Int_t i=0;i<19;i++){
                    graphNLOPi0ScaledandNColl_0010->RemovePoint(graphNLOPi0ScaledandNColl_0010->GetN()-1);
                    graphNLOPi0ScaledandNColl_2050->RemovePoint(graphNLOPi0ScaledandNColl_2050->GetN()-1);
                }
                for(Int_t i=0;i<4;i++){
                    graphNLOEtaonlyNColl_0010->RemovePoint(graphNLOEtaonlyNColl_0010->GetN()-1);
                    graphNLOEtaScaledandNColl_2050->RemovePoint(graphNLOEtaScaledandNColl_2050->GetN()-1);
                }

                DrawGammaSetMarkerTGraphAsym(graphNLOPi0ScaledandNColl_0010, 0, 0, colorNLO+1, colorNLO+1, 1, kTRUE, kRed-7);
                graphNLOPi0ScaledandNColl_0010->SetFillStyle(3001);
                DrawGammaSetMarkerTGraphAsym(graphNLOPi0ScaledandNColl_2050, 0, 0, colorNLO+1, colorNLO+1, 1, kTRUE, kAzure-4);
                graphNLOPi0ScaledandNColl_2050->SetFillStyle(3001);
                DrawGammaSetMarkerTGraphAsym(graphNLOEtaonlyNColl_0010, 0, 0, colorNLO, colorNLO, 1, kTRUE, kRed-9);
                DrawGammaSetMarkerTGraphAsym(graphNLOEtaScaledandNColl_2050, 0, 0, colorNLO, colorNLO, 1, kTRUE, kAzure-9);

                TGraphAsymmErrors *TheoryCracowPi0LowPt_scaledx400_0010 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_0010->Clone();
                TheoryCracowPi0LowPt_scaledx400_0010 = ScaleGraph(TheoryCracowPi0LowPt_scaledx400_0010,400);
                DrawGammaSetMarkerTGraphAsym(TheoryCracowPi0LowPt_scaledx400_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE, colorCracow0010);
                TGraphAsymmErrors *TheoryCracowPi0LowPt_scaledx100_2050 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_2050->Clone();
                TheoryCracowPi0LowPt_scaledx100_2050 = ScaleGraph(TheoryCracowPi0LowPt_scaledx100_2050,100);
                DrawGammaSetMarkerTGraphAsym(TheoryCracowPi0LowPt_scaledx100_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
                TGraphAsymmErrors *TheoryCracowEtaLowPt_notscaled_0010 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_0010->Clone();
//                 TheoryCracowEtaLowPt_notscaled_0010 = ScaleGraph(TheoryCracowEtaLowPt_notscaled_0010,4);
                DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaLowPt_notscaled_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE, colorCracow0010);
                TGraphAsymmErrors *TheoryCracowEtaLowPt_scaleddown_2050 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_2050->Clone();
                TheoryCracowEtaLowPt_scaleddown_2050 = ScaleGraph(TheoryCracowEtaLowPt_scaleddown_2050,2e-1);
                DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaLowPt_scaleddown_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);

                TGraphAsymmErrors *TheoryBegunEQPi0_scaledx400_0010 = (TGraphAsymmErrors*)TheoryBegunEQPi0_0010->Clone();
                TheoryBegunEQPi0_scaledx400_0010 = ScaleGraph(TheoryBegunEQPi0_scaledx400_0010,400);
                TGraphAsymmErrors *TheoryBegunEQPi0_scaledx100_2050 = (TGraphAsymmErrors*)TheoryBegunEQPi0_2050->Clone();
                TheoryBegunEQPi0_scaledx100_2050 = ScaleGraph(TheoryBegunEQPi0_scaledx100_2050,100);
                TGraphAsymmErrors *TheoryBegunEQEta_notscaled_0010 = (TGraphAsymmErrors*)TheoryBegunEQEta_0010->Clone();
//                 TheoryBegunEQEta_notscaled_0010 = ScaleGraph(TheoryBegunEQEta_notscaled_0010,4);
                TGraphAsymmErrors *TheoryBegunEQEta_scaleddown_2050 = (TGraphAsymmErrors*)TheoryBegunEQEta_2050->Clone();
                TheoryBegunEQEta_scaleddown_2050 = ScaleGraph(TheoryBegunEQEta_scaleddown_2050,2e-1);

                TGraphAsymmErrors *graphInvYieldPi0SysBothMeson_scaledx400_0010 = (TGraphAsymmErrors*)graphInvYieldPi0SysBothMeson_0010->Clone("");
                graphInvYieldPi0SysBothMeson_scaledx400_0010 = ScaleGraph(graphInvYieldPi0SysBothMeson_scaledx400_0010,400);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0SysBothMeson_scaledx400_0010, 24, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
                graphInvYieldPi0SysBothMeson_scaledx400_0010->Draw("E2same");

                TGraphAsymmErrors *graphInvYieldPi0StatBothMeson_scaledx400_0010 = (TGraphAsymmErrors*)graphInvYieldPi0StatBothMeson_0010->Clone("");
                graphInvYieldPi0StatBothMeson_scaledx400_0010 = ScaleGraph(graphInvYieldPi0StatBothMeson_scaledx400_0010,400);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0StatBothMeson_scaledx400_0010, 24, markerSizeComb, colorCombo0010, colorCombo0010);
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvYieldPi0StatBothMeson_scaledx400_0010);
                graphInvYieldPi0StatBothMeson_scaledx400_0010->Draw("p,same");

                TGraphAsymmErrors *graphInvYieldPi0SysBothMeson_scaledx100_2050 = (TGraphAsymmErrors*)graphInvYieldPi0SysBothMeson_2050->Clone("");
                graphInvYieldPi0SysBothMeson_scaledx100_2050 = ScaleGraph(graphInvYieldPi0SysBothMeson_scaledx100_2050,100);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0SysBothMeson_scaledx100_2050, 25, markerSizeComb, colorCombo2050 ,colorCombo2050, 2, kTRUE);
                graphInvYieldPi0SysBothMeson_scaledx100_2050->Draw("E2same");
                TGraphAsymmErrors *graphInvYieldPi0StatBothMeson_scaledx100_2050 = (TGraphAsymmErrors*)graphInvYieldPi0StatBothMeson_2050->Clone("");
                graphInvYieldPi0StatBothMeson_scaledx100_2050 = ScaleGraph(graphInvYieldPi0StatBothMeson_scaledx100_2050,100);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0StatBothMeson_scaledx100_2050, 25, markerSizeComb,  colorCombo2050 ,colorCombo2050);
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvYieldPi0StatBothMeson_scaledx100_2050);
                graphInvYieldPi0StatBothMeson_scaledx100_2050->Draw("p,same");

                TGraphAsymmErrors *graphInvYieldPi0SysBothMeson_notscaled_2050 = (TGraphAsymmErrors*)graphInvYieldEtaSysBothMeson_0010->Clone("");
//              graphInvYieldPi0SysBothMeson_notscaled_2050 = ScaleGraph(graphInvYieldPi0SysBothMeson_notscaled_2050,4);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0SysBothMeson_notscaled_2050, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
                graphInvYieldPi0SysBothMeson_notscaled_2050->Draw("E2same");

                TGraphAsymmErrors *graphInvYieldEtaStatBothMeson_notscaled_0010 = (TGraphAsymmErrors*)graphInvYieldEtaStatBothMeson_0010->Clone("");
//              graphInvYieldEtaStatBothMeson_notscaled_0010 = ScaleGraph(graphInvYieldEtaStatBothMeson_notscaled_0010,4);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaStatBothMeson_notscaled_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010);
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvYieldEtaStatBothMeson_notscaled_0010);
                graphInvYieldEtaStatBothMeson_notscaled_0010->Draw("p,same");

                TGraphAsymmErrors *graphInvYieldEtaSysBothMeson_scaleddown_2050 = (TGraphAsymmErrors*)graphInvYieldEtaSysBothMeson_2050->Clone("");
                graphInvYieldEtaSysBothMeson_scaleddown_2050 = ScaleGraph(graphInvYieldEtaSysBothMeson_scaleddown_2050,2e-1);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaSysBothMeson_scaleddown_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050, 2, kTRUE);
                graphInvYieldEtaSysBothMeson_scaleddown_2050->Draw("E2same");

                TGraphAsymmErrors *graphInvYieldEtaStatBothMeson_scaleddown_2050 = (TGraphAsymmErrors*)graphInvYieldEtaStatBothMeson_2050->Clone("");
                graphInvYieldEtaStatBothMeson_scaleddown_2050 = ScaleGraph(graphInvYieldEtaStatBothMeson_scaleddown_2050,2e-1);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaStatBothMeson_scaleddown_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050);
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvYieldEtaStatBothMeson_scaleddown_2050);
                graphInvYieldEtaStatBothMeson_scaleddown_2050->Draw("p,same");

                TGraphAsymmErrors* copydummyPi0= (TGraphAsymmErrors*)graphInvYieldPi0StatBothMeson_0010->Clone();
                DrawGammaSetMarkerTGraphAsym(copydummyPi0, 0, 0, colorNLO+1, colorNLO+1, 1, kTRUE, colorNLO+1);
                copydummyPi0->SetFillStyle(3001);
                TGraphAsymmErrors* copydummyEta= (TGraphAsymmErrors*)graphInvYieldEtaStatBothMeson_0010->Clone();
                DrawGammaSetMarkerTGraphAsym(copydummyEta, 0, 0, colorNLO+1, colorNLO+1, 1, kTRUE, colorNLO+1);


                TLegend* legendInvYieldSectionPi0LHC11h_onlyPbPb = new TLegend(0.2,0.13,0.7,0.23);
                legendInvYieldSectionPi0LHC11h_onlyPbPb->SetMargin(0.17);
                legendInvYieldSectionPi0LHC11h_onlyPbPb->SetFillColor(0);
                legendInvYieldSectionPi0LHC11h_onlyPbPb->SetLineColor(0);
                legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextFont(42);
                legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextSize(FontSize);
                legendInvYieldSectionPi0LHC11h_onlyPbPb->SetNColumns(2);
                legendInvYieldSectionPi0LHC11h_onlyPbPb->SetHeader(collisionSystem2760GeV.Data());
                legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphInvYieldPi0SysBothMeson_notscaled_2050,Form("#eta,   %s",cent0010.Data()),"fp");
                legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphInvYieldPi0SysBothMeson_scaledx400_0010,Form("#pi^{0},   %s",cent0010.Data()),"fp");
                legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphInvYieldEtaSysBothMeson_scaleddown_2050,Form("#eta, %s",cent2050.Data()),"fp");
                legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphInvYieldPi0SysBothMeson_scaledx100_2050,Form("#pi^{0}, %s",cent2050.Data()),"fp");
                legendInvYieldSectionPi0LHC11h_onlyPbPb->Draw();

                labelFactorPi0104->Draw();
                labelFactorPi0100->Draw();
//                 labelFactorEta4->Draw();
                labelFactorEta->Draw();
                if(thesisPlotting) thesisLabelLowLeft->Draw();

                canvasInvariantYield->SaveAs(Form("%s/YieldCombinedLHC11h_BothMesonsDataOnly.%s",outputDir.Data(),suffix.Data()));
//                 canvasInvariantYield->SaveAs(Form("%s/YieldCombinedLHC11h_BothMesonsDataOnly.%s",paperPlots.Data(),suffix.Data()));


                canvasInvariantYield->cd();
                histo2DInvariantYieldWithModels->DrawCopy();

                    graphInvYieldPi0SysBothMeson_scaledx400_0010->Draw("E2same");
                    graphInvYieldPi0StatBothMeson_scaledx400_0010->Draw("p,same");

                    graphInvYieldPi0SysBothMeson_scaledx100_2050->Draw("E2same");
                    graphInvYieldPi0StatBothMeson_scaledx100_2050->Draw("p,same");

                    graphInvYieldPi0SysBothMeson_notscaled_2050->Draw("E2same");
                    graphInvYieldEtaStatBothMeson_notscaled_0010->Draw("p,same");

                    graphInvYieldEtaSysBothMeson_scaleddown_2050->Draw("E2same");
                    graphInvYieldEtaStatBothMeson_scaleddown_2050->Draw("p,same");

                    legendInvYieldSectionPi0LHC11h_onlyPbPb->Draw();

                    graphNLOPi0ScaledandNColl_0010->Draw("3,same");
                    graphNLOPi0ScaledandNColl_2050->Draw("3,same");
                    graphNLOEtaonlyNColl_0010->Draw("3,same");
                    graphNLOEtaScaledandNColl_2050->Draw("3,same");
                    TLegend* legNLO = new TLegend(0.55,0.75,0.9,0.95);
                    legNLO->SetMargin(0.17);
                    legNLO->SetFillColor(0);
                    legNLO->SetLineColor(0);
                    legNLO->SetTextFont(42);
                    legNLO->SetTextSize(FontSize);
                    legNLO->SetHeader("pQCD NLO #mu = #it{p}_{T}");
    //              legNLO->AddEntry(graphNLOEtaonlyNColl_0010,"pQCD NLO #mu = #it{p}_{T}","l");
                    legNLO->AddEntry(copydummyPi0,"for #pi^{0}, MSTW + DSS14","f");
                    legNLO->AddEntry((TObject*)0,"PRD 91 014035","");
                    legNLO->AddEntry(copydummyEta,"for #eta, MSTW + DSS07","f");
                    legNLO->AddEntry((TObject*)0,"PRD 91 014035","");
                    legNLO->AddEntry((TObject*)0,"all scaled by <N_{coll}>","");
                    legNLO->Draw();
                    legNLO->Draw();

                    labelFactorPi0104->Draw();
                    labelFactorPi0100->Draw();
                    labelFactorEta->Draw();

                canvasInvariantYield->SaveAs(Form("%s/YieldCombinedLHC11h_BothMesonswithNLO.%s",outputDir.Data(),suffix.Data()));
//                 canvasInvariantYield->SaveAs(Form("%s/YieldCombinedLHC11h_BothMesonswithNLO.%s",paperPlots.Data(),suffix.Data()));

                canvasInvariantYield->cd();
                histo2DInvariantYieldWithModels->DrawCopy();

                    ProduceGraphErrWithoutXErrors(graphEPOSPi0forPlot_0010);
                    ProduceGraphErrWithoutXErrors(graphEPOSPi0forPlot_2050);
                    ProduceGraphErrWithoutXErrors(graphEPOSEtaforPlot_0010);
                    ProduceGraphErrWithoutXErrors(graphEPOSEtaforPlot_2050);

                    graphEPOSPi0forPlot_0010 = ScaleGraph(graphEPOSPi0forPlot_0010,400);
                    DrawGammaSetMarkerTGraphErr(graphEPOSPi0forPlot_0010,0, 0, colorEPOS0010,colorEPOS0010, 2, kTRUE, colorEPOS0010);
                    graphEPOSPi0forPlot_0010->SetLineStyle(5);
//                     graphEPOSPi0forPlot_0010->GetXaxis()->SetRangeUser(0.,13.);
                    graphEPOSPi0forPlot_0010->Draw("same,l,x");

                    graphEPOSPi0forPlot_2050 = ScaleGraph(graphEPOSPi0forPlot_2050,100);
                    DrawGammaSetMarkerTGraphErr(graphEPOSPi0forPlot_2050,5,2,colorEPOS2050,colorEPOS2050);
                    graphEPOSPi0forPlot_2050->SetLineStyle(5);
//                     graphEPOSPi0forPlot_2050->GetXaxis()->SetRangeUser(0.,13.);
                    graphEPOSPi0forPlot_2050->Draw("same,l,x");

                    DrawGammaSetMarkerTGraphErr(graphEPOSEtaforPlot_0010,5,2,colorEPOS0010,colorEPOS0010);
                    graphEPOSEtaforPlot_0010->SetLineStyle(5);
//                     graphEPOSEtaforPlot_0010->GetXaxis()->SetRangeUser(0.,13.);
                    graphEPOSEtaforPlot_0010->Draw("same,l,x");

                    graphEPOSEtaforPlot_2050 = ScaleGraph(graphEPOSEtaforPlot_2050,2e-1);
                    DrawGammaSetMarkerTGraphErr(graphEPOSEtaforPlot_2050,5,2,colorEPOS2050,colorEPOS2050);
                    graphEPOSEtaforPlot_2050->SetLineStyle(5);
//                     graphEPOSEtaforPlot_2050->GetXaxis()->SetRangeUser(0.,13.);
                    graphEPOSEtaforPlot_2050->Draw("same,l,x");

                    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQPi0_scaledx400_0010, 0, 0, colorBegun0010,colorBegun0010, 5, kTRUE, colorBegun0010);
                    TheoryBegunEQPi0_scaledx400_0010->SetLineStyle(9);
                    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQEta_notscaled_0010, 0, 0, colorBegun0010,colorBegun0010, 5, kTRUE, colorBegun0010);
                    TheoryBegunEQEta_notscaled_0010->SetLineStyle(9);
                    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQPi0_scaledx100_2050, 0, 0, colorBegun2050,colorBegun2050, 5, kTRUE,colorBegun2050);
                    TheoryBegunEQPi0_scaledx100_2050->SetLineStyle(9);
                    DrawGammaSetMarkerTGraphAsym(TheoryBegunEQEta_scaleddown_2050, 0, 0, colorBegun2050,colorBegun2050, 5, kTRUE,colorBegun2050);
                    TheoryBegunEQEta_scaleddown_2050->SetLineStyle(9);

                    TheoryCracowPi0LowPt_scaledx400_0010->Draw("l,same");
                    TheoryCracowPi0LowPt_scaledx100_2050->Draw("l,same");
                    TheoryCracowEtaLowPt_notscaled_0010->Draw("l,same");
                    TheoryCracowEtaLowPt_scaleddown_2050->Draw("l,same");

                    TheoryBegunEQPi0_scaledx400_0010->Draw("l,same");
                    TheoryBegunEQEta_notscaled_0010->Draw("l,same");
                    TheoryBegunEQPi0_scaledx100_2050->Draw("l,same");
                    TheoryBegunEQEta_scaleddown_2050->Draw("l,same");

                    graphInvYieldPi0StatBothMeson_scaledx400_0010->SetMarkerSize(3);
                    graphInvYieldPi0StatBothMeson_scaledx100_2050->SetMarkerSize(3);
                    graphInvYieldEtaStatBothMeson_notscaled_0010->SetMarkerSize(3);
                    graphInvYieldEtaStatBothMeson_scaleddown_2050->SetMarkerSize(3);
                    graphInvYieldPi0SysBothMeson_scaledx400_0010->Draw("E2same");
                    graphInvYieldPi0StatBothMeson_scaledx400_0010->Draw("p,same");
                    graphInvYieldPi0SysBothMeson_scaledx100_2050->Draw("E2same");
                    graphInvYieldPi0StatBothMeson_scaledx100_2050->Draw("p,same");
                    graphInvYieldPi0SysBothMeson_notscaled_2050->Draw("E2same");
                    graphInvYieldEtaStatBothMeson_notscaled_0010->Draw("p,same");
                    graphInvYieldEtaSysBothMeson_scaleddown_2050->Draw("E2same");
                    graphInvYieldEtaStatBothMeson_scaleddown_2050->Draw("p,same");

                    legendInvYieldSectionPi0LHC11h_onlyPbPb->Draw();

                    TGraphAsymmErrors* dummyTheoryBegunEQPi0_scaledx400_0010 = (TGraphAsymmErrors*)TheoryBegunEQPi0_scaledx400_0010->Clone();
                    dummyTheoryBegunEQPi0_scaledx400_0010->SetLineStyle(7);
                    TGraphAsymmErrors* dummyTheoryBegunEQPi0_scaledx100_2050 = (TGraphAsymmErrors*)TheoryBegunEQPi0_scaledx100_2050->Clone();
                    dummyTheoryBegunEQPi0_scaledx100_2050->SetLineStyle(7);

                    TLegend* legtheoryNEQ = new TLegend(0.48,0.95-(4*0.035),0.95,0.95);
                    legtheoryNEQ->SetFillColor(0);
                    legtheoryNEQ->SetLineColor(0);
                    legtheoryNEQ->SetTextFont(42);
                    legtheoryNEQ->SetMargin(0.2);
                    legtheoryNEQ->SetNColumns(2);
                    legtheoryNEQ->SetTextSize(FontSize);
//                     legtheoryNEQ->SetTextAlign(31);
                    legtheoryNEQ->SetHeader("SHM - PRC 90, 014906 (2014)");
                    legtheoryNEQ->AddEntry((TObject*)0," 0#font[122]{-}10%","");
                    legtheoryNEQ->AddEntry((TObject*)0,"20#font[122]{-}50%","");
                    legtheoryNEQ->AddEntry(TheoryCracowPi0LowPt_scaledx400_0010,"NEQ ", "l");
                    legtheoryNEQ->AddEntry(dummyTheoryBegunEQPi0_scaledx400_0010,"EQ ", "l");
                    legtheoryNEQ->AddEntry(TheoryCracowPi0LowPt_scaledx100_2050,"NEQ ", "l");
                    legtheoryNEQ->AddEntry(dummyTheoryBegunEQPi0_scaledx100_2050,"EQ ", "l");
                    legtheoryNEQ->Draw("same");

                    TLegend* legtheoryEPOS = new TLegend(0.58,0.8-(0.035*2.7),0.95,0.8);
                    legtheoryEPOS->SetFillColor(0);
                    legtheoryEPOS->SetLineColor(0);
                    legtheoryEPOS->SetTextFont(42);
                    legtheoryEPOS->SetMargin(0.17);
                    legtheoryEPOS->SetTextSize(FontSize);
                    legtheoryEPOS->SetHeader("PRC 85, 064907 (2012)");
                    legtheoryEPOS->AddEntry(graphEPOSEtaforPlot_0010,"EPOS 0#font[122]{-}10%", "l");
                    legtheoryEPOS->AddEntry(graphEPOSEtaforPlot_2050,"EPOS 20#font[122]{-}50%", "l");
                    legtheoryEPOS->Draw("same");

                    labelFactorPi0104->Draw();
                    labelFactorPi0100->Draw();
                    labelFactorEta->Draw();
                    TLatex *thesisLabelLowLeft = new TLatex(0.2,0.24,thisthesis.Data());
                    SetStyleTLatex( thesisLabelLowLeft,FontSize,4);
//                     if(thesisPlotting) thesisLabelLowLeft->Draw();

                canvasInvariantYield->SaveAs(Form("%s/YieldCombinedLHC11h_BothMesonswithBothModels.%s",outputDir.Data(),suffix.Data()));
                canvasInvariantYield->SaveAs(Form("%s/YieldCombinedLHC11h_BothMesonswithBothModels.%s",paperPlots.Data(),suffix.Data()));

            }
        }


    canvasInvariantYield->cd();
    histo2DInvariantYieldwithPP->Draw("copy");

        TGraphAsymmErrors *graphNLOforPi0PP = (TGraphAsymmErrors*)graphNLOCalcDSS14InvYieldPi0Band->Clone("graphNLOforPi0PP");
        TGraphAsymmErrors *graphNLOforPi0PbPb0010 = (TGraphAsymmErrors*)graphNLOCalcDSS14InvYieldPi0Band->Clone("graphNLOforPi0PbPb0010");
        graphNLOforPi0PbPb0010 = ScaleGraph(graphNLOforPi0PbPb0010,nColl0010*4);
        TGraphAsymmErrors *graphNLOforPi0PbPb2050 = (TGraphAsymmErrors*)graphNLOCalcDSS14InvYieldPi0Band->Clone("graphNLOforPi0PbPb2050");
        graphNLOforPi0PbPb2050 = ScaleGraph(graphNLOforPi0PbPb2050,nColl2050);

        TGraphAsymmErrors *graphNLOforEtaPP = (TGraphAsymmErrors*)graphNLOCalcDSS07InvYieldEtaBand->Clone("graphNLOforEtaPP");
        TGraphAsymmErrors *graphNLOforEtaPbPb0010 = (TGraphAsymmErrors*)graphNLOCalcDSS07InvYieldEtaBand->Clone("graphNLOforEtaPbPb0010");
        graphNLOforEtaPbPb0010 = ScaleGraph(graphNLOforEtaPbPb0010,nColl0010*4);
        TGraphAsymmErrors *graphNLOforEtaPbPb2050 = (TGraphAsymmErrors*)graphNLOCalcDSS07InvYieldEtaBand->Clone("graphNLOforEtaPbPb2050");
        graphNLOforEtaPbPb2050 = ScaleGraph(graphNLOforEtaPbPb2050,nColl2050);

        TLatex *labelNLO;
        labelNLO= new TLatex(0.61,0.75,"PDF: MSTW, FF: DSS14");
        SetStyleTLatex( labelNLO, 0.035,4);

        TLegend* legendSpectraPPnlo = new TLegend(0.2,0.12,0.52,0.34);//0.6,0.7,0.96,0.87);
        legendSpectraPPnlo->SetFillColor(0);
        legendSpectraPPnlo->SetLineColor(0);
        legendSpectraPPnlo->SetMargin(0.17);
        legendSpectraPPnlo->SetTextFont(42);
        legendSpectraPPnlo->SetTextSize(FontSize);

        if(meson.CompareTo("Pi0")==0){

            legendSpectraPPnlo->SetHeader(collisionSystemPP2760GeV.Data());
            legendSpectraPPnlo->AddEntry(graphInvSectionCombSysPi02760GeVPlot,"EPJC 77 (2017) 339"/*"arXiv: 1702.00917"*/,"pf");
            legendSpectraPPnlo->AddEntry(graphNLOforPi0PP,"pQCD NLO #mu = #it{p}_{T}","l");
            legendSpectraPPnlo->AddEntry((TObject*)0,"PDF: MSTW, FF: DSS14","");
            legendSpectraPPnlo->AddEntry((TObject*)0,"PRD 91 014035","");
            legendSpectraPPnlo->AddEntry((TObject*)0,"scaled by <N_{coll}>","");
            legendSpectraPPnlo->Draw();

            DrawGammaSetMarkerTGraphAsym(graphNLOforPi0PP, 0, 0, colorNLO, colorNLO, 3, kTRUE, colorNLO);
            graphNLOforPi0PP->Draw("3,same");
            DrawGammaSetMarkerTGraphAsym(graphNLOforPi0PbPb0010, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
            graphNLOforPi0PbPb0010->Draw("3,same");
            DrawGammaSetMarkerTGraphAsym(graphNLOforPi0PbPb2050, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
            graphNLOforPi0PbPb2050->Draw("3,same");

            graphInvSectionCombStatPi02760GeVPlot->Draw("p,same");
            graphInvSectionCombSysPi02760GeVPlot->Draw("E2same");


        } else if(meson.CompareTo("Eta")==0){

            legendSpectraPPnlo->SetHeader(collisionSystemPP2760GeV.Data());
            legendSpectraPPnlo->AddEntry(graphInvSectionCombSysEta2760GeVPlot,"EPJC 77 (2017) 339"/*"arXiv: 1702.00917"*/,"pf");
            legendSpectraPPnlo->AddEntry(graphNLOforEtaPP,"pQCD NLO #mu = #it{p}_{T}","l");
            legendSpectraPPnlo->AddEntry((TObject*)0,"PDF: MSTW, FF: DSS07","");
            legendSpectraPPnlo->AddEntry((TObject*)0,"PRD 91 014035","");
            legendSpectraPPnlo->AddEntry((TObject*)0,"scaled by <N_{coll}>","");
            legendSpectraPPnlo->Draw();

            DrawGammaSetMarkerTGraphAsym(graphNLOforEtaPP, 0, 0, colorNLO, colorNLO, 3, kTRUE, colorNLO);
            graphNLOforEtaPP->Draw("3,same");
            DrawGammaSetMarkerTGraphAsym(graphNLOforEtaPbPb0010, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
            graphNLOforEtaPbPb0010->Draw("3,same");
            DrawGammaSetMarkerTGraphAsym(graphNLOforEtaPbPb2050, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
            graphNLOforEtaPbPb2050->Draw("3,same");

            graphInvSectionCombStatEta2760GeVPlot->Draw("p,same");
            graphInvSectionCombSysEta2760GeVPlot->Draw("E2same");
        }

        TLegend* legendInvYieldSectionPi0LHC11h_WitPPandNLO = new TLegend(0.595,0.79,0.91,0.92); //0.17,0.13,0.5,0.24);
        legendInvYieldSectionPi0LHC11h_WitPPandNLO->SetFillColor(0);
        legendInvYieldSectionPi0LHC11h_WitPPandNLO->SetMargin(0.17);
        legendInvYieldSectionPi0LHC11h_WitPPandNLO->SetLineColor(0);
        legendInvYieldSectionPi0LHC11h_WitPPandNLO->SetTextFont(42);
        legendInvYieldSectionPi0LHC11h_WitPPandNLO->SetTextSize(FontSize);
        legendInvYieldSectionPi0LHC11h_WitPPandNLO->SetHeader(collisionSystem2760GeV.Data());
        legendInvYieldSectionPi0LHC11h_WitPPandNLO->AddEntry(graphCombInvYieldSysPbPb2760GeVScaled_0010,Form("  %s scaled by 4",cent0010.Data()),"fp");
        legendInvYieldSectionPi0LHC11h_WitPPandNLO->AddEntry(graphCombInvYieldSysPbPb2760GeV_2050,Form("%s",cent2050.Data()),"fp");
        legendInvYieldSectionPi0LHC11h_WitPPandNLO->Draw();
        labelSyst->Draw();

        graphCombInvYieldSysPbPb2760GeVScaled_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeVScaled_0010->Draw("p,same");

        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

    canvasInvariantYield->Update();
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnlyWithPPandNLO.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnlyWithPPandNLO.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnlyWithPPandNLO.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    canvasInvariantYield->cd();
        histo2DInvariantYield->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(TheoryBegunEQ_0010ScaledforPlot, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE, colorCracow0010);
        DrawGammaSetMarkerTGraphAsym(TheoryBegunEQ_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
        TheoryBegunEQ_0010ScaledforPlot->Draw("l,same");
        TheoryBegunEQ_2050->Draw("l,same");

        graphEPOSpredScaledforPlot_0010->Draw("same,l,x");
        graphEPOSpred_2050->Draw("same,l,x");

        graphCombInvYieldSysPbPb2760GeVScaled_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeVScaled_0010->Draw("p,same");

        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

        TLegend* legtheory2 = new TLegend(0.2,0.12,0.52,0.12+(0.04*2.7));
        legtheory2->SetFillColor(0);
        legtheory2->SetLineColor(0);
        legtheory2->SetTextFont(42);
        legtheory2->SetTextSize(FontSize);
        legtheory2->SetHeader("PRC 85, 064907 (2012)");
        legtheory2->AddEntry(graphEPOSpred_0010,"EPOS 0#font[122]{-}10%", "l");
        legtheory2->AddEntry(graphEPOSpred_2050,"EPOS 20#font[122]{-}50%", "l");
        legtheory2->Draw("same");

        TLegend* legtheoryEQ = new TLegend(0.2,0.24,0.6,0.24+(0.04*3));
        legtheoryEQ->SetFillColor(0);
        legtheoryEQ->SetLineColor(0);
        legtheoryEQ->SetTextFont(42);
        legtheoryEQ->SetTextSize(FontSize);
        legtheoryEQ->SetHeader("PRC 90, 014906 (2014)");
        legtheoryEQ->AddEntry(TheoryBegunEQ_0010ScaledforPlot,"EQ SHM 0#font[122]{-}10%", "l");
        legtheoryEQ->AddEntry(TheoryBegunEQ_2050,"EQ SHM 20#font[122]{-}50%", "l");
        legtheoryEQ->Draw("same");

        legendInvYieldSectionPi0LHC11h_WitPP->Draw();
        labelFactorLowerOnlyPbPb->Draw();
        labelSyst->Draw();

        TLatex *thesisLabelModel = new TLatex(0.73,0.927,"(this thesis)");
        SetStyleTLatex( thesisLabelModel,FontSize,4);
        if(thesisPlotting) thesisLabelModel->Draw();

    canvasInvariantYield->Update();
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsandEQ.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsandEQ.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    canvasInvariantYield->cd();
        histo2DInvariantYield->DrawCopy();

        TheoryCracowLowPtScaledforPlot_0010->Draw("l,same");
        TheoryCracowLowPt_2050->Draw("l,same");

        graphEPOSpredScaledforPlot_0010->Draw("same,l,x");
        graphEPOSpred_2050->Draw("same,l,x");

        graphCombInvYieldSysPbPb2760GeVScaled_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeVScaled_0010->Draw("p,same");

        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

        legtheory2->Draw("same");

        TLegend* legtheoryNEQ = new TLegend(0.2,0.24,0.6,0.24+(0.04*3));
        legtheoryNEQ->SetFillColor(0);
        legtheoryNEQ->SetLineColor(0);
        legtheoryNEQ->SetTextFont(42);
        legtheoryNEQ->SetTextSize(FontSize);
        legtheoryNEQ->SetHeader("PRC 90, 014906 (2014)");
        legtheoryNEQ->AddEntry(TheoryCracowLowPtScaledforPlot_0010,"NEQ SHM 0#font[122]{-}10%", "l");
        legtheoryNEQ->AddEntry(TheoryCracowLowPt_2050,"NEQ SHM 20#font[122]{-}50%", "l");
        legtheoryNEQ->Draw("same");

        legendInvYieldSectionPi0LHC11h_WitPP->Draw();
        labelFactorLowerOnlyPbPb->Draw();
        labelSyst->Draw();

        if(thesisPlotting) thesisLabelModel->Draw();

    canvasInvariantYield->Update();
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsandNEQ.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsandNEQ.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    canvasInvariantYield->cd();
        histo2DInvariantYield->DrawCopy();

        TheoryCracowLowPtScaledforPlot_0010->Draw("l,same");
        TheoryCracowLowPt_2050->Draw("l,same");

//         SetStyleGammaNLOTGraphWithBand( TheoryBegunEQ_0010ScaledforPlot, 1.0, 1, colorBegun0010, 3001, colorBegun0010, 0);
        DrawGammaSetMarkerTGraphAsym(TheoryBegunEQ_0010ScaledforPlot, 0, 0, colorBegun0010,colorBegun0010, 5, kTRUE, colorBegun0010);
        TheoryBegunEQ_0010ScaledforPlot->SetLineStyle(9);
        DrawGammaSetMarkerTGraphAsym(TheoryBegunEQ_2050, 0, 0, colorBegun2050,colorBegun2050, 5, kTRUE,colorBegun2050);
        TheoryBegunEQ_2050->SetLineStyle(9);
        TheoryBegunEQ_0010ScaledforPlot->Draw("3pl,same");
        TheoryBegunEQ_2050->Draw("l,same");

        graphEPOSpredScaledforPlot_0010->Draw("same,l,x");
        graphEPOSpred_2050->Draw("same,l,x");

        graphCombInvYieldSysPbPb2760GeVScaled_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeVScaled_0010->Draw("p,same");

        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

        TLegend* legtheoryBothSHM = new TLegend(0.2,0.12,0.6,0.12+(0.04*3.7));
        legtheoryBothSHM->SetFillColor(0);
        legtheoryBothSHM->SetLineColor(0);
        legtheoryBothSHM->SetTextFont(42);
        legtheoryBothSHM->SetTextSize(FontSize);
        legtheoryBothSHM->SetHeader("SHM - PRC 90, 014906 (2014)");
        legtheoryBothSHM->SetNColumns(2);
//         legtheoryBothSHM->SetMargin(0.17);
        legtheoryBothSHM->AddEntry((TObject*)0,cent0010.Data(),"");
        legtheoryBothSHM->AddEntry((TObject*)0,cent2050.Data(),"");
        legtheoryBothSHM->AddEntry(TheoryCracowLowPtScaledforPlot_0010,"NEQ", "l");
        legtheoryBothSHM->AddEntry(TheoryBegunEQ_0010ScaledforPlot,"EQ", "l");
        legtheoryBothSHM->AddEntry(TheoryCracowLowPt_2050,"NEQ", "l");
        legtheoryBothSHM->AddEntry(TheoryBegunEQ_2050,"EQ", "l");
        legtheoryBothSHM->Draw("same");

        TLegend* legtheoryEPOS = new TLegend(0.2,0.27,0.52,0.27+(0.04*2.7));
        legtheoryEPOS->SetFillColor(0);
        legtheoryEPOS->SetLineColor(0);
        legtheoryEPOS->SetTextFont(42);
        legtheoryEPOS->SetTextSize(FontSize);
        legtheoryEPOS->SetHeader("PRC 85, 064907 (2012)");
        legtheoryEPOS->AddEntry(graphEPOSpred_0010,"EPOS 0#font[122]{-}10%", "l");
        legtheoryEPOS->AddEntry(graphEPOSpred_2050,"EPOS 20#font[122]{-}50%", "l");
        legtheoryEPOS->Draw("same");

        legendInvYieldSectionPi0LHC11h_WitPP->Draw();
        labelFactorLowerOnlyPbPb->Draw();
        labelSyst->Draw();
//      labelPreliminary->Draw();

        if(thesisPlotting)thesisLabelModel->Draw();


    canvasInvariantYield->Update();
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsV2.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsV2.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsV2.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    canvasInvariantYield->cd();
    histo2DInvariantYieldwithPP->Draw("copy");

        graphCombInvYieldSysPbPb2760GeVScaled_0010->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeVScaled_0010->Draw("p,same");

        graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
        graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

        legendInvYieldSectionPi0LHC11h_WitPP->Draw();

        TLegend* legendSpectraPP2 = new TLegend(0.595,0.72,0.91,0.79);
        legendSpectraPP2->SetFillColor(0);
        legendSpectraPP2->SetLineColor(0);
        legendSpectraPP2->SetTextFont(42);
        legendSpectraPP2->SetMargin(0.17);
        legendSpectraPP2->SetTextSize(FontSize);
        legendSpectraPP2->SetHeader(collisionSystemPP2760GeV.Data());
        if(meson.CompareTo("Pi0")==0){
            graphInvSectionCombStatPi02760GeVPlot->Draw("p,same");
            graphInvSectionCombSysPi02760GeVPlot->Draw("E2same");
            legendSpectraPP2->AddEntry(graphInvSectionCombSysPi02760GeVPlot,"EPJC 77 (2017) 339"/*"arXiv: 1702.00917"*/,"pf");
        } else if(meson.CompareTo("Eta")==0){
            graphInvSectionCombStatEta2760GeVPlot->Draw("p,same");
            graphInvSectionCombSysEta2760GeVPlot->Draw("E2same");
            legendSpectraPP2->AddEntry(graphInvSectionCombSysEta2760GeVPlot,"EPJC 77 (2017) 339"/*"arXiv: 1702.00917"*/,"pf");
        }
        legendSpectraPP2->Draw();
        labelSyst->Draw();
        labelFactorLower->Draw();

        TheoryCracowLowPtScaledforPlot_0010->Draw("l,same");
        TheoryCracowLowPt_2050->Draw("l,same");

        TheoryBegunEQ_0010ScaledforPlot->Draw("l,same");
        TheoryBegunEQ_2050->Draw("l,same");

        graphEPOSpredScaledforPlot_0010->Draw("same,l,x");
        graphEPOSpred_2050->Draw("same,l,x");

        legtheoryEQ->Draw("same");
        legtheory2->Draw("same");

    canvasInvariantYield->Update();
    canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithPPandModels.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasInvariantYield->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModels_withPP.%s",paperPlots.Data(),meson.Data(),suffix.Data()));


    //*************************************************************************************************************
    //***************************** Plot with data + models and ratios ***********************************************
    //*************************************************************************************************************
    //ratios to models
    Double_t arrayYieldandModelRatioX_4pads[2];
    Double_t arrayYieldandModelRatioY_4pads[6];
    Double_t relativeMarginsX_YieldModelRatio[3];
    Double_t relativeMarginsY_YieldModelRatio[3];
    textSizeLabelsPixel = 90;
    ReturnCorrectValuesForCanvasScaling(2500,-4000, 1, 5,0.13, 0.025, 0.003, 0.05, arrayYieldandModelRatioX_4pads, arrayYieldandModelRatioY_4pads, relativeMarginsX_YieldModelRatio, relativeMarginsY_YieldModelRatio);
    TCanvas* canvasInvYieldSectionRatios = new TCanvas("canvasInvYieldSectionRatios","",0,0,2500,4000);  // gives the page size

        TPad* padInvSectionSpec = new TPad("padInvSectionSpec", "", arrayYieldandModelRatioX_4pads[0], arrayYieldandModelRatioY_4pads[3], arrayYieldandModelRatioX_4pads[1], arrayYieldandModelRatioY_4pads[0],-1, -1, -2);
        DrawGammaPadSettings( padInvSectionSpec, relativeMarginsX_YieldModelRatio[0], relativeMarginsX_YieldModelRatio[2], relativeMarginsY_YieldModelRatio[0], relativeMarginsY_YieldModelRatio[1]);
        padInvSectionSpec->Draw();
        Double_t marginXSec = relativeMarginsX_YieldModelRatio[0]*2500;
        Double_t textsizeLabelsXSecUp = 0;
        Double_t textsizeFacXSecUp = 0;
        if (padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) < padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1())){
            textsizeLabelsXSecUp = (Double_t)textSizeLabelsPixel/padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
            textsizeFacXSecUp = (Double_t)1./padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
        } else {
            textsizeLabelsXSecUp = (Double_t)textSizeLabelsPixel/padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
            textsizeFacXSecUp = (Double_t)1./padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
        }

        TPad* padInvSectionLowPtRatio = new TPad("padInvSectionLowPtRatio", "", arrayYieldandModelRatioX_4pads[0], arrayYieldandModelRatioY_4pads[4], arrayYieldandModelRatioX_4pads[1], arrayYieldandModelRatioY_4pads[3],-1, -1, -2);
        DrawGammaPadSettings( padInvSectionLowPtRatio, relativeMarginsX_YieldModelRatio[0], relativeMarginsX_YieldModelRatio[2], relativeMarginsY_YieldModelRatio[1], relativeMarginsY_YieldModelRatio[1]);
        padInvSectionLowPtRatio->Draw();
        textsizeLabelsXSecMiddle = 0;
        textsizeFacXSecMiddle = 0;
        if (padInvSectionLowPtRatio->XtoPixel(padInvSectionLowPtRatio->GetX2()) < padInvSectionLowPtRatio->YtoPixel(padInvSectionLowPtRatio->GetY1())){
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionLowPtRatio->XtoPixel(padInvSectionLowPtRatio->GetX2()) ;
            textsizeFacXSecMiddle = (Double_t)1./padInvSectionLowPtRatio->XtoPixel(padInvSectionLowPtRatio->GetX2()) ;
        } else {
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionLowPtRatio->YtoPixel(padInvSectionLowPtRatio->GetY1());
            textsizeFacXSecMiddle = (Double_t)1./padInvSectionLowPtRatio->YtoPixel(padInvSectionLowPtRatio->GetY1());
        }

        TPad* padInvSectionEPOSRatio = new TPad("padInvSectionEPOSRatio", "", arrayYieldandModelRatioX_4pads[0], arrayYieldandModelRatioY_4pads[5], arrayYieldandModelRatioX_4pads[1], arrayYieldandModelRatioY_4pads[4],-1, -1, -2);
        DrawGammaPadSettings( padInvSectionEPOSRatio, relativeMarginsX_YieldModelRatio[0], relativeMarginsX_YieldModelRatio[2], relativeMarginsY_YieldModelRatio[1], relativeMarginsY_YieldModelRatio[2]);
        padInvSectionEPOSRatio->Draw();
        textsizeLabelsXSecDown = 0;
        textsizeFacXSecDown = 0;
        if (padInvSectionEPOSRatio->XtoPixel(padInvSectionEPOSRatio->GetX2()) < padInvSectionEPOSRatio->YtoPixel(padInvSectionEPOSRatio->GetY1())){
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionEPOSRatio->XtoPixel(padInvSectionEPOSRatio->GetX2()) ;
            textsizeFacXSecDown = (Double_t)1./padInvSectionEPOSRatio->XtoPixel(padInvSectionEPOSRatio->GetX2()) ;
        } else {
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionEPOSRatio->YtoPixel(padInvSectionEPOSRatio->GetY1());
            textsizeFacXSecDown = (Double_t)1./padInvSectionEPOSRatio->YtoPixel(padInvSectionEPOSRatio->GetY1());
        }

        if(meson.CompareTo("Pi0")==0){
            histo2DInvariantYieldAndModels = new TH2F("histo2DInvariantYieldAndModels","histo2DInvariantYieldAndModels",1000,0.25,70.,1000,2e-8,1e4);
            SetStyleHistoTH2ForGraphs(histo2DInvariantYieldAndModels, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
        } else {
            histo2DInvariantYieldAndModels = new TH2F("histo2DInvariantYieldAndModels","histo2DInvariantYieldAndModels",1000,0.25,70.,1000,2e-9,1e3);
            SetStyleHistoTH2ForGraphs(histo2DInvariantYieldAndModels, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
        }

        TH2F * ratio2DLowPt = new TH2F("ratio2DLowPt","ratio2DLowPt",1000,0.25,70.,1000,0.2,2.1);
        SetStyleHistoTH2ForGraphs(ratio2DLowPt, "#it{p}_{T} (GeV/#it{c})","#frac{NEQ SHM, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
        ratio2DLowPt->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DLowPt->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DLowPt->GetYaxis()->SetLabelOffset(+0.01);
        TH2F * ratio2DPythia =  new TH2F("ratio2DPythia","ratio2DPythia",1000,0.25,70.,1000,0.2,2.1);
        SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{EPOS, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                      0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
        ratio2DPythia->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DPythia->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DPythia->GetYaxis()->SetLabelOffset(+0.01);

        if(thesisPlotting || PaperPi0){
            histo2DInvariantYieldAndModels->GetXaxis()->SetRangeUser(0.5,30.);
            ratio2DLowPt->GetXaxis()->SetRangeUser(0.5,30.);
            ratio2DPythia->GetXaxis()->SetRangeUser(0.5,30.);
        } else {
            histo2DInvariantYieldAndModels->GetXaxis()->SetRangeUser(0.25,30.);
            ratio2DLowPt->GetXaxis()->SetRangeUser(0.25,30.);
            ratio2DPythia->GetXaxis()->SetRangeUser(0.25,30.);
        }

        padInvSectionSpec->cd();
        padInvSectionSpec->SetLogy();
        padInvSectionSpec->SetLogx();
        histo2DInvariantYieldAndModels->Draw("copy");

            TGraphAsymmErrors *TheoryCracowLowPtForPlot_0010 = (TGraphAsymmErrors*)TheoryCracowLowPt_0010->Clone("TheoryCracowLowPtForPlot_0010");
            DrawGammaSetMarkerTGraphAsym(TheoryCracowLowPtForPlot_0010, 0, 0, colorCracowRatio,colorCracowRatio, 5, kTRUE,colorCracowRatio);
//             TheoryCracowLowPtForPlot_0010->Draw("same");

            TheoryBegunEQ_0010ScaledforPlot->Draw("same");

            graphCombInvYieldSysPbPb2760GeV_0010->Draw("E2same");
            graphCombInvYieldStatPbPb2760GeV_0010->Draw("p,same");

            fitBylinkinPbPb2760GeVPtLHC11h_0010->Draw("same");
            fitBylinkinPbPb2760GeVPtLHC11h_0010->SetLineStyle(2);

            TGraphErrors *graphEPOSpredForPlot_0010 = (TGraphErrors*)graphEPOSpred_0010->Clone("graphEPOSpredForPlot_0010");
            DrawGammaSetMarkerTGraphErr(graphEPOSpredForPlot_0010,0,0,colorEPOSRatio,colorEPOSRatio,2,kTRUE,colorEPOSRatio);
            graphEPOSpredForPlot_0010->SetLineStyle(5);
            graphEPOSpredForPlot_0010->Draw("same,l,x");

            TLegend* legendYieldAndModelRatio = GetAndSetLegend2(0.17, 0.22-5*0.045, 0.5, 0.22, 0.85*textSizeLabelsPixel);
            legendYieldAndModelRatio->SetNColumns(1);
            legendYieldAndModelRatio->SetMargin(0.2);
            if(meson.CompareTo("Pi0")==0) legendYieldAndModelRatio->AddEntry(graphCombInvYieldSysPbPb2760GeV_0010,"#pi^{0} ALICE","pf");
            else if(meson.CompareTo("Eta")==0) legendYieldAndModelRatio->AddEntry(graphCombInvYieldSysPbPb2760GeV_0010,"#eta ALICE","pf");
            legendYieldAndModelRatio->AddEntry(TheoryCracowLowPtForPlot_0010,"NEQ SHM","l");
            legendYieldAndModelRatio->AddEntry(TheoryBegunEQ_0010ScaledforPlot,"EQ SHM","l");
            legendYieldAndModelRatio->AddEntry(graphEPOSpredForPlot_0010,"EPOS","l");
            legendYieldAndModelRatio->AddEntry(fitBylinkinPbPb2760GeVPtLHC11h_0010,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}","l");
            legendYieldAndModelRatio->Draw();

            if(meson.CompareTo("Pi0")==0){
                labelDetSysRatioToModelsLHC11h= new TLatex(0.5,0.93,"#pi^{0} #rightarrow #gamma#gamma");
            } else if(meson.CompareTo("Eta")==0){
                labelDetSysRatioToModelsLHC11h= new TLatex(0.5,0.93,"#eta #rightarrow #gamma#gamma");
            }
            SetStyleTLatex( labelDetSysRatioToModelsLHC11h, 0.035,4);
            labelDetSysRatioToModelsLHC11h->Draw();
            TLatex *labelEnergyRatioToModelsLHC11h = new TLatex(0.5,0.88,collisionSystemPbPb0010.Data());
            SetStyleTLatex( labelEnergyRatioToModelsLHC11h, 0.035,4);
            labelEnergyRatioToModelsLHC11h->Draw();

        padInvSectionLowPtRatio->cd();
        padInvSectionLowPtRatio->SetLogx(1);
        ratio2DLowPt->Draw("copy");
        if(thesisPlotting){
            DrawGammaLines(0.5, 30.,1., 1.,0.1,kGray);
        } else {
            DrawGammaLines(0.25, 30.,1., 1.,0.1,kGray);
        }

            TGraphAsymmErrors* graphRatioCombLowPtPbPb2760GeV_0010 = (TGraphAsymmErrors*)TheoryCracowLowPt_0010->Clone("graphRatioCombLowPtPbPb2760GeV_0010");
            graphRatioCombLowPtPbPb2760GeV_0010 = CalculateGraphErrRatioToFit(graphRatioCombLowPtPbPb2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioCombLowPtPbPb2760GeV_0010, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);
            graphRatioCombLowPtPbPb2760GeV_0010->Draw("l,same");

            TGraphAsymmErrors* graphRatioCombBegunEQPbPb2760GeV_0010 = (TGraphAsymmErrors*)TheoryBegunEQ_0010->Clone("graphRatioCombBegunEQPbPb2760GeV_0010");
            graphRatioCombBegunEQPbPb2760GeV_0010 = CalculateGraphErrRatioToFit(graphRatioCombBegunEQPbPb2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioCombBegunEQPbPb2760GeV_0010, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);
            graphRatioCombBegunEQPbPb2760GeV_0010->Draw("l,same");

            DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatPbPb2760GeV_0010, markerStyleComb, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kFALSE);
            graphRatioCombCombFitStatPbPb2760GeV_0010->SetLineWidth(widthLinesBoxes);
            graphRatioCombCombFitStatPbPb2760GeV_0010->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysPbPb2760GeV_0010, markerStyleComb, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
            graphRatioCombCombFitSysPbPb2760GeV_0010->SetLineWidth(0);
            graphRatioCombCombFitSysPbPb2760GeV_0010->Draw("2,same");

        ratio2DLowPt->Draw("axis,same");
        padInvSectionEPOSRatio->cd();
        padInvSectionEPOSRatio->SetLogx(1);
        ratio2DPythia->Draw("copy");
        if(thesisPlotting){
            DrawGammaLines(0.5, 30.,1., 1.,0.1,kGray);
        } else {
            DrawGammaLines(0.25, 30.,1., 1.,0.1,kGray);
        }


            TGraphErrors* graphRatioEPOSToFitPbPb2760GeV_0010 = (TGraphErrors*) graphEPOSpred_0010->Clone("graphRatioEPOSToFitPbPb2760GeV_0010");
            ProduceGraphErrWithoutXErrors(graphRatioEPOSToFitPbPb2760GeV_0010);
            graphRatioEPOSToFitPbPb2760GeV_0010 = CalculateGraphErrRatioToFit(graphRatioEPOSToFitPbPb2760GeV_0010, fitBylinkinPbPb2760GeVPtLHC11h_0010);
            DrawGammaSetMarkerTGraphErr(graphRatioEPOSToFitPbPb2760GeV_0010, 0, 0, colorEPOSRatio,colorEPOSRatio, 2, kTRUE, colorEPOSRatio);
            graphRatioEPOSToFitPbPb2760GeV_0010->SetLineStyle(5);
            graphRatioEPOSToFitPbPb2760GeV_0010->Draw("same,l,x");

            graphRatioCombCombFitStatPbPb2760GeV_0010->Draw("p,same");
            graphRatioCombCombFitSysPbPb2760GeV_0010->Draw("2,same");

        ratio2DPythia->Draw("axis,same");
    canvasInvYieldSectionRatios->Update();
    //  canvasInvYieldSectionRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    //  canvasInvYieldSectionRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_0010.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
    canvasInvYieldSectionRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasInvYieldSectionRatios->cd();
        padInvSectionSpec->cd();
        padInvSectionSpec->SetLogy();
        padInvSectionSpec->SetLogx();
        histo2DInvariantYieldAndModels->Draw("copy");

            TGraphAsymmErrors *TheoryCracowLowPtForPlot_2050 = (TGraphAsymmErrors*)TheoryCracowLowPt_2050->Clone("TheoryCracowLowPtForPlot_2050");
            DrawGammaSetMarkerTGraphAsym(TheoryCracowLowPtForPlot_2050, 0, 0, colorCracowRatio,colorCracowRatio, 5, kTRUE,colorCracowRatio);
            TheoryCracowLowPtForPlot_2050->Draw("same");

            TGraphAsymmErrors *TheoryBegunEQForPlot_2050 = (TGraphAsymmErrors*)TheoryBegunEQ_2050->Clone("TheoryBegunEQForPlot_2050");
            TheoryBegunEQForPlot_2050->Draw("same");

            graphCombInvYieldSysPbPb2760GeV_2050->Draw("E2same");
            graphCombInvYieldStatPbPb2760GeV_2050->Draw("p,same");

            fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");
            fitBylinkinPbPb2760GeVPtLHC11h_2050->SetLineStyle(2);

            TGraphErrors *graphEPOSpredForPlot_2050 = (TGraphErrors*)graphEPOSpred_2050->Clone("graphEPOSpredForPlot_2050");
            DrawGammaSetMarkerTGraphErr(graphEPOSpredForPlot_2050,0,0,colorEPOSRatio,colorEPOSRatio,2,kTRUE,colorEPOSRatio);
            graphEPOSpredForPlot_2050->SetLineStyle(5);
            graphEPOSpredForPlot_2050->Draw("same,l,x");

            TLegend* legendYieldAndModelRatio2 = GetAndSetLegend2(0.17, 0.22-5*0.045, 0.5, 0.22, 0.85*textSizeLabelsPixel);
            legendYieldAndModelRatio2->SetNColumns(1);
            legendYieldAndModelRatio2->SetMargin(0.2);
            if(meson.CompareTo("Pi0")==0) legendYieldAndModelRatio2->AddEntry(graphCombInvYieldSysPbPb2760GeV_2050,"#pi^{0} ALICE","pf");
            else if(meson.CompareTo("Eta")==0) legendYieldAndModelRatio2->AddEntry(graphCombInvYieldSysPbPb2760GeV_2050,"#eta ALICE","pf");
            legendYieldAndModelRatio2->AddEntry(TheoryCracowLowPtForPlot_2050,"NEQ SHM","l");
            legendYieldAndModelRatio2->AddEntry(TheoryBegunEQForPlot_2050,"EQ SHM","l");
            legendYieldAndModelRatio2->AddEntry(graphEPOSpredForPlot_2050,"EPOS","l");
            legendYieldAndModelRatio2->AddEntry(fitBylinkinPbPb2760GeVPtLHC11h_2050,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}","l");
            legendYieldAndModelRatio2->Draw();

            labelDetSysRatioToModelsLHC11h->Draw();
            labelEnergyRatioToModelsLHC11h->Draw();

        padInvSectionLowPtRatio->cd();
        padInvSectionLowPtRatio->SetLogx(1);
        ratio2DLowPt->Draw("copy");
        if(thesisPlotting){
            DrawGammaLines(0.5, 30.,1., 1.,0.1,kGray);
        } else {
            DrawGammaLines(0.25, 30.,1., 1.,0.1,kGray);
        }

            TGraphAsymmErrors* graphRatioCombLowPtPbPb2760GeV_2050 = (TGraphAsymmErrors*)TheoryCracowLowPt_2050->Clone("graphRatioCombLowPtPbPb2760GeV_2050");
            graphRatioCombLowPtPbPb2760GeV_2050 = CalculateGraphErrRatioToFit(graphRatioCombLowPtPbPb2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);
            DrawGammaSetMarkerTGraphAsym(graphRatioCombLowPtPbPb2760GeV_2050, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);
            graphRatioCombLowPtPbPb2760GeV_2050->Draw("l,same");

            TGraphAsymmErrors* graphRatioCombBegunEQPbPb2760GeV_2050 = (TGraphAsymmErrors*)TheoryBegunEQ_2050->Clone("graphRatioCombBegunEQPbPb2760GeV_2050");
            graphRatioCombBegunEQPbPb2760GeV_2050 = CalculateGraphErrRatioToFit(graphRatioCombBegunEQPbPb2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);
            graphRatioCombBegunEQPbPb2760GeV_2050->Draw("l,same");

            DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatPbPb2760GeV_2050, markerStyleComb, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kFALSE);
            graphRatioCombCombFitStatPbPb2760GeV_2050->SetLineWidth(widthLinesBoxes);
            graphRatioCombCombFitStatPbPb2760GeV_2050->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysPbPb2760GeV_2050, markerStyleComb, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
            graphRatioCombCombFitSysPbPb2760GeV_2050->SetLineWidth(0);
            graphRatioCombCombFitSysPbPb2760GeV_2050->Draw("2,same");

        ratio2DLowPt->Draw("axis,same");
        padInvSectionEPOSRatio->cd();
        padInvSectionEPOSRatio->SetLogx(1);
        ratio2DPythia->Draw("copy");
        if(thesisPlotting){
            DrawGammaLines(0.5, 30.,1., 1.,0.1,kGray);
        } else {
            DrawGammaLines(0.25, 30.,1., 1.,0.1,kGray);
        }

            TGraphErrors* graphRatioEPOSToFitPbPb2760GeV_2050 = (TGraphErrors*) graphEPOSpred_2050->Clone("graphRatioEPOSToFitPbPb2760GeV_2050");
            ProduceGraphErrWithoutXErrors(graphRatioEPOSToFitPbPb2760GeV_2050);
            graphRatioEPOSToFitPbPb2760GeV_2050 = CalculateGraphErrRatioToFit(graphRatioEPOSToFitPbPb2760GeV_2050, fitBylinkinPbPb2760GeVPtLHC11h_2050);
            DrawGammaSetMarkerTGraphErr(graphRatioEPOSToFitPbPb2760GeV_2050, 0, 0, colorEPOSRatio,colorEPOSRatio, 2, kTRUE, colorEPOSRatio);
            graphRatioEPOSToFitPbPb2760GeV_2050->SetLineStyle(5);
            graphRatioEPOSToFitPbPb2760GeV_2050->Draw("same,l,x");

            graphRatioCombCombFitStatPbPb2760GeV_2050->Draw("p,same");
            graphRatioCombCombFitSysPbPb2760GeV_2050->Draw("2,same");

        ratio2DPythia->Draw("axis,same");
    canvasInvYieldSectionRatios->Update();
    //  canvasInvYieldSectionRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    //  canvasInvYieldSectionRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_2050.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
    canvasInvYieldSectionRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));



    Double_t arrayYieldandModelOnlyRatioX[2];
    Double_t arrayYieldandModelOnlyRatioY[3];
    Double_t relativeMarginsX_YieldModelOnlyRatio[3];
    Double_t relativeMarginsY_YieldModelOnlyRatio[3];
    textSizeLabelsPixel = 90;
    ReturnCorrectValuesForCanvasScaling(2500, 2000, 1, 2, 0.13, 0.025, 0.025, 0.1, arrayYieldandModelOnlyRatioX, arrayYieldandModelOnlyRatioY, relativeMarginsX_YieldModelOnlyRatio, relativeMarginsY_YieldModelOnlyRatio);

    TCanvas* canvasInvYieldSectionOnlyRatios = new TCanvas("canvasInvYieldSectionOnlyRatios","",0,0,2500,2000);  //0,0,2500,4000);  // gives the page size

        TPad* padInvSectionLowerRatio = new TPad("","",arrayYieldandModelOnlyRatioX[0], arrayYieldandModelOnlyRatioY[2], arrayYieldandModelOnlyRatioX[1], arrayYieldandModelOnlyRatioY[1],-1, -1, -2);
        DrawGammaPadSettings( padInvSectionLowerRatio, relativeMarginsX_YieldModelOnlyRatio[0], relativeMarginsX_YieldModelOnlyRatio[2], relativeMarginsY_YieldModelOnlyRatio[1], relativeMarginsY_YieldModelOnlyRatio[2]);

        TPad* padInvSectionUpperRatio = new TPad("","",arrayYieldandModelOnlyRatioX[0], arrayYieldandModelOnlyRatioY[1], arrayYieldandModelOnlyRatioX[1], arrayYieldandModelOnlyRatioY[0],-1, -1, -2);
        DrawGammaPadSettings( padInvSectionUpperRatio, relativeMarginsX_YieldModelOnlyRatio[0], relativeMarginsX_YieldModelOnlyRatio[2], relativeMarginsY_YieldModelOnlyRatio[0], relativeMarginsY_YieldModelOnlyRatio[1]);

        padInvSectionLowerRatio->Draw();
        if (padInvSectionLowerRatio->XtoPixel(padInvSectionLowerRatio->GetX2()) < padInvSectionLowerRatio->YtoPixel(padInvSectionLowerRatio->GetY1())){
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionLowerRatio->XtoPixel(padInvSectionLowerRatio->GetX2()) ;
            textsizeFacXSecDown = (Double_t)1./padInvSectionLowerRatio->XtoPixel(padInvSectionLowerRatio->GetX2()) ;
        } else {
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionLowerRatio->YtoPixel(padInvSectionLowerRatio->GetY1());
            textsizeFacXSecDown = (Double_t)1./padInvSectionLowerRatio->YtoPixel(padInvSectionLowerRatio->GetY1());
        }
        padInvSectionUpperRatio->Draw();
        if (padInvSectionUpperRatio->XtoPixel(padInvSectionUpperRatio->GetX2()) < padInvSectionUpperRatio->YtoPixel(padInvSectionUpperRatio->GetY1())){
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionUpperRatio->XtoPixel(padInvSectionUpperRatio->GetX2()) ;
            textsizeFacXSecMiddle = (Double_t)1./padInvSectionUpperRatio->XtoPixel(padInvSectionUpperRatio->GetX2()) ;
        } else {
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionUpperRatio->YtoPixel(padInvSectionUpperRatio->GetY1());
            textsizeFacXSecMiddle = (Double_t)1./padInvSectionUpperRatio->YtoPixel(padInvSectionUpperRatio->GetY1());
        }

        TH2F * ratio2DUpperPadOnlyRatio = new TH2F("ratio2DUpperPadOnlyRatio","ratio2DUpperPadOnlyRatio",1000,0.25,70.,1000,0.01,2.1);
        SetStyleHistoTH2ForGraphs(ratio2DUpperPadOnlyRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                    0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
        if(thesisPlotting)
            ratio2DUpperPadOnlyRatio->GetXaxis()->SetRangeUser(0.5,30);
        else
            ratio2DUpperPadOnlyRatio->GetXaxis()->SetRangeUser(0.25,30.);
        ratio2DUpperPadOnlyRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DUpperPadOnlyRatio->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DUpperPadOnlyRatio->GetYaxis()->SetLabelOffset(+0.01);
//      ratio2DUpperPadOnlyRatio->GetXaxis()->SetTickLength(0.07);

        TH2F * ratio2DLowerPadOnlyRatio =  new TH2F("ratio2DLowerPadOnlyRatio","ratio2DLowerPadOnlyRatio",1000,0.25,70.,1000,0.01,2.1);
        SetStyleHistoTH2ForGraphs(ratio2DLowerPadOnlyRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                    0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
        if(thesisPlotting) ratio2DLowerPadOnlyRatio->GetXaxis()->SetRangeUser(0.25,30);
        else ratio2DLowerPadOnlyRatio->GetXaxis()->SetRangeUser(0.5,30);
        ratio2DLowerPadOnlyRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DLowerPadOnlyRatio->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DLowerPadOnlyRatio->GetYaxis()->SetLabelOffset(+0.01);
//      ratio2DLowerPadOnlyRatio->GetXaxis()->SetTickLength(0.07);


        padInvSectionUpperRatio->cd();
        padInvSectionUpperRatio->SetLogx(1);
        ratio2DUpperPadOnlyRatio->DrawCopy();
            if(thesisPlotting)
                DrawGammaLines(0.5,30,1., 1.,0.1,kGray);
            else
                DrawGammaLines(0.25,30,1., 1.,0.1,kGray);

            graphRatioCombLowPtPbPb2760GeV_0010->Draw("l,same");
            graphRatioCombCombFitStatPbPb2760GeV_0010->Draw("p,same");
            graphRatioCombCombFitSysPbPb2760GeV_0010->Draw("2,same");
            graphRatioEPOSToFitPbPb2760GeV_0010->Draw("same,l,x");
    //         histoRatioEPOS3ToFitPbPb2760GeV_0010->Draw("same,l");

            TLatex *labelcent0010 = new TLatex(0.8,0.82," 0#font[122]{-}10%");
            SetStyleTLatex( labelcent0010, textsizeLabelsXSecMiddle,4);
            labelcent0010->Draw();

                if(meson.CompareTo("Pi0")==0)legendXsectionPaperOnlyRatios = GetAndSetLegend2(0.55, 0.45-4*0.065, 0.95, 0.45, 0.85* textSizeLabelsPixel);
                else if(meson.CompareTo("Eta")==0)legendXsectionPaperOnlyRatios = GetAndSetLegend2(0.55, 0.45-4*0.065, 0.95, 0.45, 0.85* textSizeLabelsPixel);
                legendXsectionPaperOnlyRatios->SetNColumns(1);
                legendXsectionPaperOnlyRatios->SetMargin(0.2);
                legendXsectionPaperOnlyRatios->SetNColumns(2);
                legendXsectionPaperOnlyRatios->SetHeader(collisionSystem2760GeV.Data());
                if(meson.CompareTo("Pi0")==0) legendXsectionPaperOnlyRatios->AddEntry(graphRatioCombCombFitSysPbPb2760GeV_0010,"#pi^{0} ALICE","pf");
                else if(meson.CompareTo("Eta")==0) legendXsectionPaperOnlyRatios->AddEntry(graphRatioCombCombFitSysPbPb2760GeV_0010,"#eta ALICE","pf");
                legendXsectionPaperOnlyRatios->AddEntry(graphRatioCombLowPtPbPb2760GeV_0010,"NEQ SHM","l");
                legendXsectionPaperOnlyRatios->AddEntry((TObject*)0,"","");
                legendXsectionPaperOnlyRatios->AddEntry(graphRatioEPOSToFitPbPb2760GeV_0010,"EPOS","l");
        //         legendXsectionPaperOnlyRatios->AddEntry(histoRatioEPOS3ToFitPbPb2760GeV_0010,"EPOS3","l");

                if(thesisPlotting){
                    SetStyleTLatex( thesisLabelHighLeft,textsizeLabelsXSecMiddle,4);
                    thesisLabelHighLeft->Draw();
                }

            ratio2DUpperPadOnlyRatio->Draw("axis,same");
        padInvSectionLowerRatio->cd();
        padInvSectionLowerRatio->SetLogx(1);
        ratio2DLowerPadOnlyRatio->DrawCopy();
            if(thesisPlotting)
                DrawGammaLines(0.5,30,1., 1.,0.1,kGray);
            else
                DrawGammaLines(0.25,30,1., 1.,0.1,kGray);

            graphRatioCombLowPtPbPb2760GeV_2050->Draw("l,same");
            graphRatioCombCombFitStatPbPb2760GeV_2050->Draw("p,same");
            graphRatioCombCombFitSysPbPb2760GeV_2050->Draw("2,same");
            graphRatioEPOSToFitPbPb2760GeV_2050->Draw("same,l,x");

            TLatex *labelcent2050 = new TLatex(0.8,0.88,"20#font[122]{-}50%");
            SetStyleTLatex( labelcent2050, textsizeLabelsXSecDown,4);
            labelcent2050->Draw();
            legendXsectionPaperOnlyRatios->Draw();

        ratio2DLowerPadOnlyRatio->Draw("axis,same");
    canvasInvYieldSectionOnlyRatios->Update();
    canvasInvYieldSectionOnlyRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsOnlyRatios.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasInvYieldSectionOnlyRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsOnlyRatios.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
    canvasInvYieldSectionOnlyRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsOnlyRatios.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    if(MesonInput){
        TGraphAsymmErrors *graphInvYieldPi0Stat_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0CombPbPb2760GeVStatErr_0010");
        TGraphAsymmErrors *graphInvYieldPi0Syst_0010= (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0CombPbPb2760GeVSysErr_0010");

        TGraphAsymmErrors *graphInvYieldPi0Stat_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0CombPbPb2760GeVStatErr_2050");
        TGraphAsymmErrors *graphInvYieldPi0Syst_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0CombPbPb2760GeVSysErr_2050");

        TGraphAsymmErrors *graphInvYieldEtaStat_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaCombPbPb2760GeVStatErr_0010");
        TGraphAsymmErrors *graphInvYieldEtaSyst_0010= (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaCombPbPb2760GeVSysErr_0010");

        TGraphAsymmErrors *graphInvYieldEtaStat_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaCombPbPb2760GeVStatErr_2050");
        TGraphAsymmErrors *graphInvYieldEtaSyst_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaCombPbPb2760GeVSysErr_2050");

        TF1 *fitPi0InvYieldPbPb2760GeV_0010 = (TF1*)MesonInput->Get("FitToYieldPi0_0010");
        TF1 *fitPi0InvYieldPbPb2760GeV_2050 = (TF1*)MesonInput->Get("FitToYieldPi0_2050");
        TF1 *fitEtaInvYieldPbPb2760GeV_0010 = (TF1*)MesonInput->Get("FitToYieldEta_0010");
        TF1 *fitEtaInvYieldPbPb2760GeV_2050 = (TF1*)MesonInput->Get("FitToYieldEta_2050");

        if(graphInvYieldPi0Stat_0010 && graphInvYieldEtaStat_0010 && fitPi0InvYieldPbPb2760GeV_0010 && fitEtaInvYieldPbPb2760GeV_0010){

        ////////// data
            TGraphAsymmErrors* graphRatioPi0InvYieldFitStatPbPb2760GeV_0010 = (TGraphAsymmErrors*)graphInvYieldPi0Stat_0010->Clone();
            graphRatioPi0InvYieldFitStatPbPb2760GeV_0010 = CalculateGraphErrRatioToFit(graphRatioPi0InvYieldFitStatPbPb2760GeV_0010, fitPi0InvYieldPbPb2760GeV_0010);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioPi0InvYieldFitStatPbPb2760GeV_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0InvYieldFitStatPbPb2760GeV_0010, markerStyleComb, markerSizeComb, kBlack, kBlack);

            TGraphAsymmErrors* graphRatioPi0InvYieldFitSystPbPb2760GeV_0010 = (TGraphAsymmErrors*)graphInvYieldPi0Syst_0010->Clone();
            graphRatioPi0InvYieldFitSystPbPb2760GeV_0010 = CalculateGraphErrRatioToFit(graphRatioPi0InvYieldFitSystPbPb2760GeV_0010, fitPi0InvYieldPbPb2760GeV_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0InvYieldFitSystPbPb2760GeV_0010,markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);

            TGraphAsymmErrors* graphRatioPi0InvYieldFitStatPbPb2760GeV_2050 = (TGraphAsymmErrors*)graphInvYieldPi0Stat_2050->Clone();
            graphRatioPi0InvYieldFitStatPbPb2760GeV_2050 = CalculateGraphErrRatioToFit(graphRatioPi0InvYieldFitStatPbPb2760GeV_2050, fitPi0InvYieldPbPb2760GeV_2050);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioPi0InvYieldFitStatPbPb2760GeV_2050);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0InvYieldFitStatPbPb2760GeV_2050, markerStyleComb, markerSizeComb, kBlack, kBlack);

            TGraphAsymmErrors* graphRatioPi0InvYieldFitSystPbPb2760GeV_2050 = (TGraphAsymmErrors*)graphInvYieldPi0Syst_2050->Clone();
            graphRatioPi0InvYieldFitSystPbPb2760GeV_2050 = CalculateGraphErrRatioToFit(graphRatioPi0InvYieldFitSystPbPb2760GeV_2050, fitPi0InvYieldPbPb2760GeV_2050);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0InvYieldFitSystPbPb2760GeV_2050, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);

            TGraphAsymmErrors* graphRatioEtaInvYieldFitStatPbPb2760GeV_0010 = (TGraphAsymmErrors*)graphInvYieldEtaStat_0010->Clone();
            graphRatioEtaInvYieldFitStatPbPb2760GeV_0010 = CalculateGraphErrRatioToFit(graphRatioEtaInvYieldFitStatPbPb2760GeV_0010, fitEtaInvYieldPbPb2760GeV_0010);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioEtaInvYieldFitStatPbPb2760GeV_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaInvYieldFitStatPbPb2760GeV_0010, markerStyleComb, markerSizeComb, kBlack, kBlack);

            TGraphAsymmErrors* graphRatioEtaInvYieldFitSystPbPb2760GeV_0010 = (TGraphAsymmErrors*)graphInvYieldEtaSyst_0010->Clone();
            graphRatioEtaInvYieldFitSystPbPb2760GeV_0010 = CalculateGraphErrRatioToFit(graphRatioEtaInvYieldFitSystPbPb2760GeV_0010, fitEtaInvYieldPbPb2760GeV_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaInvYieldFitSystPbPb2760GeV_0010, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);

            TGraphAsymmErrors* graphRatioEtaInvYieldFitStatPbPb2760GeV_2050 = (TGraphAsymmErrors*)graphInvYieldEtaStat_2050->Clone();
            graphRatioEtaInvYieldFitStatPbPb2760GeV_2050 = CalculateGraphErrRatioToFit(graphRatioEtaInvYieldFitStatPbPb2760GeV_2050, fitEtaInvYieldPbPb2760GeV_2050);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioEtaInvYieldFitStatPbPb2760GeV_2050);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaInvYieldFitStatPbPb2760GeV_2050, markerStyleComb, markerSizeComb, kBlack, kBlack);

            TGraphAsymmErrors* graphRatioEtaInvYieldFitSystPbPb2760GeV_2050 = (TGraphAsymmErrors*)graphInvYieldEtaSyst_2050->Clone();
            graphRatioEtaInvYieldFitSystPbPb2760GeV_2050 = CalculateGraphErrRatioToFit(graphRatioEtaInvYieldFitSystPbPb2760GeV_2050, fitEtaInvYieldPbPb2760GeV_2050);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaInvYieldFitSystPbPb2760GeV_2050, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);


        ////////// Cracow
            TGraphAsymmErrors* graphRatioLowPtPi0_0010 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_0010->Clone();
            graphRatioLowPtPi0_0010 = CalculateGraphErrRatioToFit(graphRatioLowPtPi0_0010, fitPi0InvYieldPbPb2760GeV_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioLowPtPi0_0010, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);

            TGraphAsymmErrors* graphRatioLowPtPi0_2050 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_2050->Clone();
            graphRatioLowPtPi0_2050 = CalculateGraphErrRatioToFit(graphRatioLowPtPi0_2050, fitPi0InvYieldPbPb2760GeV_2050);
            DrawGammaSetMarkerTGraphAsym(graphRatioLowPtPi0_2050, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);

            TGraphAsymmErrors* graphRatioLowPtEta_0010 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_0010->Clone();
            graphRatioLowPtEta_0010 = CalculateGraphErrRatioToFit(graphRatioLowPtEta_0010, fitEtaInvYieldPbPb2760GeV_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioLowPtEta_0010, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);

            TGraphAsymmErrors* graphRatioLowPtEta_2050 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_2050->Clone();
            graphRatioLowPtEta_2050 = CalculateGraphErrRatioToFit(graphRatioLowPtEta_2050, fitEtaInvYieldPbPb2760GeV_2050);
            DrawGammaSetMarkerTGraphAsym(graphRatioLowPtEta_2050, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);

        ////////// Begun EQ
            TGraphAsymmErrors* graphRatioBegunEQPi0_0010 = (TGraphAsymmErrors*)TheoryBegunEQPi0_0010->Clone();
            graphRatioBegunEQPi0_0010 = CalculateGraphErrRatioToFit(graphRatioBegunEQPi0_0010, fitPi0InvYieldPbPb2760GeV_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioBegunEQPi0_0010, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
            graphRatioBegunEQPi0_0010->SetLineStyle(9);

            TGraphAsymmErrors* graphRatioBegunEQPi0_2050 = (TGraphAsymmErrors*)TheoryBegunEQPi0_2050->Clone();
            graphRatioBegunEQPi0_2050 = CalculateGraphErrRatioToFit(graphRatioBegunEQPi0_2050, fitPi0InvYieldPbPb2760GeV_2050);
            DrawGammaSetMarkerTGraphAsym(graphRatioBegunEQPi0_2050, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
            graphRatioBegunEQPi0_2050->SetLineStyle(9);

            TGraphAsymmErrors* graphRatioBegunEQEta_0010 = (TGraphAsymmErrors*)TheoryBegunEQEta_0010->Clone();
            graphRatioBegunEQEta_0010 = CalculateGraphErrRatioToFit(graphRatioBegunEQEta_0010, fitEtaInvYieldPbPb2760GeV_0010);
            DrawGammaSetMarkerTGraphAsym(graphRatioBegunEQEta_0010, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
            graphRatioBegunEQEta_0010->SetLineStyle(9);

            TGraphAsymmErrors* graphRatioBegunEQEta_2050 = (TGraphAsymmErrors*)TheoryBegunEQEta_2050->Clone();
            graphRatioBegunEQEta_2050 = CalculateGraphErrRatioToFit(graphRatioBegunEQEta_2050, fitEtaInvYieldPbPb2760GeV_2050);
            DrawGammaSetMarkerTGraphAsym(graphRatioBegunEQEta_2050, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
            graphRatioBegunEQEta_2050->SetLineStyle(9);

        ////////// EPOS
            TGraphErrors* graphRatioEPOSToFitPi0_0010 = (TGraphErrors*) graphEPOSPi0_0010->Clone();
            graphRatioEPOSToFitPi0_0010 = CalculateGraphErrRatioToFit(graphRatioEPOSToFitPi0_0010, fitPi0InvYieldPbPb2760GeV_0010);
            ProduceGraphErrWithoutXErrors(graphRatioEPOSToFitPi0_0010);
            graphRatioEPOSToFitPi0_0010->SetLineStyle(5);
            DrawGammaSetMarkerTGraphErr(graphRatioEPOSToFitPi0_0010, 0, 0,colorEPOSRatio,colorEPOSRatio, 2, kTRUE, colorEPOSRatio);
    //           graphRatioEPOSToFitPi0_0010->SetLineWidth(widthCommonFit);

            TGraphErrors* graphRatioEPOSToFitPi0_2050 = (TGraphErrors*) graphEPOSPi0_2050->Clone();
            graphRatioEPOSToFitPi0_2050 = CalculateGraphErrRatioToFit(graphRatioEPOSToFitPi0_2050, fitPi0InvYieldPbPb2760GeV_2050);
            ProduceGraphErrWithoutXErrors(graphRatioEPOSToFitPi0_2050);
            graphRatioEPOSToFitPi0_2050->SetLineStyle(5);
            DrawGammaSetMarkerTGraphErr(graphRatioEPOSToFitPi0_2050, 0, 0,colorEPOSRatio,colorEPOSRatio, 2, kTRUE, colorEPOSRatio);
    //           graphRatioEPOSToFitPi0_2050->SetLineWidth(widthCommonFit);

            TGraphErrors* graphRatioEPOSToFitEta_0010 = (TGraphErrors*) graphEPOSEta_0010->Clone();
            graphRatioEPOSToFitEta_0010 = CalculateGraphErrRatioToFit(graphRatioEPOSToFitEta_0010, fitEtaInvYieldPbPb2760GeV_0010);
            ProduceGraphErrWithoutXErrors(graphRatioEPOSToFitEta_0010);
            graphRatioEPOSToFitEta_0010->SetLineStyle(5);
            DrawGammaSetMarkerTGraphErr(graphRatioEPOSToFitEta_0010, 0, 0,colorEPOSRatio,colorEPOSRatio, 2, kTRUE, colorEPOSRatio);
    //           graphRatioEPOSToFitEta_0010->SetLineWidth(widthCommonFit);

            TGraphErrors* graphRatioEPOSToFitEta_2050 = (TGraphErrors*) graphEPOSEta_2050->Clone();
            graphRatioEPOSToFitEta_2050 = CalculateGraphErrRatioToFit(graphRatioEPOSToFitEta_2050, fitEtaInvYieldPbPb2760GeV_2050);
            ProduceGraphErrWithoutXErrors(graphRatioEPOSToFitEta_2050);
            graphRatioEPOSToFitEta_2050->SetLineStyle(5);
            DrawGammaSetMarkerTGraphErr(graphRatioEPOSToFitEta_2050, 0, 0,colorEPOSRatio,colorEPOSRatio, 2, kTRUE, colorEPOSRatio);
    //           graphRatioEPOSToFitEta_2050->SetLineWidth(widthCommonFit);


            //4pads yields ratio to model
            Double_t arrayX14padsratios[3];
            Double_t arrayY14padsratios[3];
            Double_t relX4padsratios[3];
            Double_t relY4padsratios[3];
            textSizeLabelsPixel4pads = 50;
            ReturnCorrectValuesForCanvasScaling(2500,1300, 2, 2,0.08, 0.025, 0.025,0.1,arrayX14padsratios,arrayY14padsratios,relX4padsratios,relY4padsratios);
            TCanvas* canvasInvYieldSectionOnlyRatios4pads = new TCanvas("canvasInvYieldSectionOnlyRatios4pads","",0,0,2500,1300);  //0,0,2500,4000);  // gives the page size

                textsizeLabelsXSecMiddle = 0;
                textsizeFacXSecMiddle = 0;
                textsizeLabelsXSecDown = 0;
                textsizeFacXSecDown = 0;

                TPad* padInvYieldLowerLeft = new TPad("padInvYieldLowerLeft", "", arrayX14padsratios[0], arrayY14padsratios[2], arrayX14padsratios[1], arrayY14padsratios[1],-1, -1, -2);
                DrawGammaPadSettings( padInvYieldLowerLeft, relX4padsratios[0], relX4padsratios[1], relY4padsratios[1], relY4padsratios[2]);
                padInvYieldLowerLeft->Draw();

                TPad* padInvYieldLowerRight = new TPad("padInvYieldLowerRight", "", arrayX14padsratios[1], arrayY14padsratios[2], arrayX14padsratios[2], arrayY14padsratios[1],-1, -1, -2);
                DrawGammaPadSettings( padInvYieldLowerRight, relX4padsratios[1], relX4padsratios[2], relY4padsratios[1], relY4padsratios[2]);
                padInvYieldLowerRight->Draw();

                TPad* padInvYieldUpperLeft = new TPad("padInvYieldUpperLeft", "", arrayX14padsratios[0], arrayY14padsratios[1], arrayX14padsratios[1], arrayY14padsratios[0],-1, -1, -2);
                DrawGammaPadSettings( padInvYieldUpperLeft, relX4padsratios[0], relX4padsratios[1], relY4padsratios[0], relY4padsratios[1]);
                padInvYieldUpperLeft->Draw();

                TPad* padInvYieldUpperRight = new TPad("padInvYieldUpperRight", "", arrayX14padsratios[1], arrayY14padsratios[1], arrayX14padsratios[2], arrayY14padsratios[0],-1, -1, -2);
                DrawGammaPadSettings( padInvYieldUpperRight, relX4padsratios[1], relX4padsratios[2], relY4padsratios[0], relY4padsratios[1]);
                padInvYieldUpperRight->Draw();

                if (padInvYieldUpperLeft->XtoPixel(padInvYieldUpperLeft->GetX2()) < padInvYieldUpperLeft->YtoPixel(padInvYieldUpperLeft->GetY1())){
                    textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel4pads/padInvYieldUpperLeft->XtoPixel(padInvYieldUpperLeft->GetX2()) ;
                    textsizeFacXSecMiddle = (Double_t)1./padInvYieldUpperLeft->XtoPixel(padInvYieldUpperLeft->GetX2()) ;
                } else {
                    textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel4pads/padInvYieldUpperLeft->YtoPixel(padInvYieldUpperLeft->GetY1());
                    textsizeFacXSecMiddle = (Double_t)1./padInvYieldUpperLeft->YtoPixel(padInvYieldUpperLeft->GetY1());
                }

                if (padInvYieldLowerLeft->XtoPixel(padInvYieldLowerLeft->GetX2()) < padInvYieldLowerLeft->YtoPixel(padInvYieldLowerLeft->GetY1())){
                    textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel4pads/padInvYieldLowerLeft->XtoPixel(padInvYieldLowerLeft->GetX2()) ;
                    textsizeFacXSecDown = (Double_t)1./padInvYieldLowerLeft->XtoPixel(padInvYieldLowerLeft->GetX2()) ;
                } else {
                    textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel4pads/padInvYieldLowerLeft->YtoPixel(padInvYieldLowerLeft->GetY1());
                    textsizeFacXSecDown = (Double_t)1./padInvYieldLowerLeft->YtoPixel(padInvYieldLowerLeft->GetY1());
                }

                TH2F *ratio2DUpPad1OnlyRatio = new TH2F("ratio2DUpPad1OnlyRatio","ratio2DUpPad1OnlyRatio",1000,0.25,70.,1000,0.01,2.1);
                SetStyleHistoTH2ForGraphs(ratio2DUpPad1OnlyRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                        0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.8, 510, 505);
                ratio2DUpPad1OnlyRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
                ratio2DUpPad1OnlyRatio->GetXaxis()->SetNoExponent(kTRUE);
                ratio2DUpPad1OnlyRatio->GetYaxis()->SetLabelOffset(+0.01);
                ratio2DUpPad1OnlyRatio->GetXaxis()->SetRangeUser(minX,maxX);

                TH2F * ratio2DLowerPad1OnlyRatio = new TH2F("ratio2DLowerPad1OnlyRatio","ratio2DLowerPad1OnlyRatio",1000,0.25,70.,1000,0.01,2.1);
                SetStyleHistoTH2ForGraphs(ratio2DLowerPad1OnlyRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                        0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 1,0.85, 510, 505);
                ratio2DLowerPad1OnlyRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
                ratio2DLowerPad1OnlyRatio->GetXaxis()->SetNoExponent(kTRUE);
                ratio2DLowerPad1OnlyRatio->GetYaxis()->SetLabelOffset(+0.01);
                ratio2DLowerPad1OnlyRatio->GetXaxis()->SetRangeUser(minX,maxX);

                padInvYieldUpperLeft->cd();
                padInvYieldUpperLeft->SetLogx();
                ratio2DUpPad1OnlyRatio->DrawCopy();
                DrawGammaLines(minX,maxX,1., 1.,0.1,kGray);

                    graphRatioLowPtPi0_0010->Draw("l,same");
                    graphRatioBegunEQPi0_0010->Draw("l,same");
                    graphRatioEPOSToFitPi0_0010->Draw("same,l,x");

                    graphRatioPi0InvYieldFitStatPbPb2760GeV_0010->Draw("p,same");
                    graphRatioPi0InvYieldFitSystPbPb2760GeV_0010->Draw("E2same");

                    TLegend* legen4PadOnlyRatios_pad1 = GetAndSetLegend2(0.62, 0.22-2*0.08, 0.9, 0.22, 0.85*textSizeLabelsPixel4pads);
                    legen4PadOnlyRatios_pad1->SetMargin(0.17);
                    legen4PadOnlyRatios_pad1->SetHeader(collisionSystem2760GeV.Data());
                    legen4PadOnlyRatios_pad1->AddEntry(graphRatioPi0InvYieldFitSystPbPb2760GeV_0010,"#pi^{0},  0#font[122]{-}10%","pf");
                    legen4PadOnlyRatios_pad1->Draw();

                    TLatex *labelThesis4Pad = new TLatex(0.19,0.83,"This thesis");
                    SetStyleTLatex( labelThesis4Pad, 1.1*textsizeLabelsXSecDown,4);
                    labelThesis4Pad->Draw();

                    TGraphAsymmErrors *dummygraphRatioBegunEQPi0_0010 = (TGraphAsymmErrors*)graphRatioBegunEQPi0_0010->Clone("dummygraphRatioBegunEQPi0_0010");
                    DrawGammaSetMarkerTGraphAsym(dummygraphRatioBegunEQPi0_0010, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
                    dummygraphRatioBegunEQPi0_0010->SetLineStyle(7);

                    TLegend* legen4PadOnlyRatiosTheory = GetAndSetLegend2(0.75, 0.67, 0.98, 0.67+(3*0.08), 0.85* textSizeLabelsPixel4pads);
                    legen4PadOnlyRatiosTheory->SetMargin(0.3);
                    legen4PadOnlyRatiosTheory->AddEntry(graphRatioCombLowPtPbPb2760GeV_0010,"NEQ SHM","l");
                    legen4PadOnlyRatiosTheory->AddEntry(/*graphRatioBegunEQPi0_0010*/dummygraphRatioBegunEQPi0_0010,"EQ SHM","l");
                    legen4PadOnlyRatiosTheory->AddEntry(graphRatioEPOSToFitPbPb2760GeV_0010,"EPOS","l");
                    legen4PadOnlyRatiosTheory->Draw();

                ratio2DUpPad1OnlyRatio->Draw("axis,same");
                padInvYieldLowerLeft->cd();
                padInvYieldLowerLeft->SetLogx();
                ratio2DLowerPad1OnlyRatio->DrawCopy();
                DrawGammaLines(minX,maxX,1., 1.,0.1,kGray);

                    graphRatioLowPtEta_0010->Draw("l,same");
                    graphRatioBegunEQEta_0010->Draw("l,same");
                    graphRatioEPOSToFitEta_0010->Draw("same,l,x");

                    graphRatioEtaInvYieldFitStatPbPb2760GeV_0010->Draw("p,same");
                    graphRatioEtaInvYieldFitSystPbPb2760GeV_0010->Draw("E2same");

                    TLatex *labelcent0010a = new TLatex(0.77,0.88," 0#font[122]{-}10%");
                    SetStyleTLatex( labelcent0010a, 0.85*textSizeLabelsPixel4pads,4);
        //                  labelcent0010a->Draw();

                    TLegend* legen4PadOnlyRatios_pad2 = GetAndSetLegend2(0.62, 0.37-2*0.07, 0.9, 0.37, 0.85*textSizeLabelsPixel4pads);
                    legen4PadOnlyRatios_pad2->SetMargin(0.17);
                    legen4PadOnlyRatios_pad2->SetHeader(collisionSystem2760GeV.Data());
                    legen4PadOnlyRatios_pad2->AddEntry(graphRatioEtaInvYieldFitSystPbPb2760GeV_0010,"#eta,  0#font[122]{-}10%","pf"); // \\sqrt{s_\\mathrm{nn}}
                    legen4PadOnlyRatios_pad2->Draw();

                ratio2DLowerPad1OnlyRatio->Draw("axis,same");
                padInvYieldUpperRight->cd();
                padInvYieldUpperRight->SetLogx();
                ratio2DUpPad1OnlyRatio->DrawCopy();
                DrawGammaLines(minX,maxX,1., 1.,0.1,kGray);

                    graphRatioLowPtPi0_2050->Draw("l,same");
                    graphRatioBegunEQPi0_2050->Draw("l,same");
                    graphRatioEPOSToFitPi0_2050->Draw("same,l,x");

                    graphRatioPi0InvYieldFitStatPbPb2760GeV_2050->Draw("p,same");
                    graphRatioPi0InvYieldFitSystPbPb2760GeV_2050->Draw("E2same");

                    TLatex *labelcent2050 = new TLatex(0.85,0.88,"20#font[122]{-}50%");
                    SetStyleTLatex( labelcent2050, 0.85*textsizeLabelsXSecDown,4);
        //                  labelcent2050->Draw();

                    TLegend* legen4PadOnlyRatios_pad3 = GetAndSetLegend2(0.53, 0.22-2*0.08, 0.85, 0.22, 0.85* textSizeLabelsPixel4pads);
                    legen4PadOnlyRatios_pad3->SetMargin(0.17);
                    legen4PadOnlyRatios_pad3->SetHeader(collisionSystem2760GeV.Data());
                    legen4PadOnlyRatios_pad3->AddEntry(graphRatioPi0InvYieldFitSystPbPb2760GeV_2050,"#pi^{0}, 20#font[122]{-}50%","pf");
                    legen4PadOnlyRatios_pad3->Draw();

                ratio2DUpPad1OnlyRatio->Draw("axis,same");
                padInvYieldLowerRight->cd();
                padInvYieldLowerRight->SetLogx(1);
                ratio2DLowerPad1OnlyRatio->DrawCopy();
                DrawGammaLines(minX,maxX,1., 1.,0.1,kGray);

                    graphRatioLowPtEta_2050->Draw("l,same");
                    graphRatioBegunEQEta_2050->Draw("l,same");
                    graphRatioEPOSToFitEta_2050->Draw("same,l,x");

                    graphRatioEtaInvYieldFitStatPbPb2760GeV_2050->Draw("p,same");
                    graphRatioEtaInvYieldFitSystPbPb2760GeV_2050->Draw("E2same");

                    TLatex *labelcent2050a = new TLatex(0.77,0.88,"20#font[122]{-}50%");
                    SetStyleTLatex( labelcent2050a, 0.85*textsizeLabelsXSecDown,4);
        //                  labelcent2050a->Draw();

                    TLegend* legen4PadOnlyRatios_pad4 = GetAndSetLegend2(0.53, 0.37-2*0.07, 0.85, 0.37, 0.85* textSizeLabelsPixel4pads);
                    legen4PadOnlyRatios_pad4->SetMargin(0.17);
                    legen4PadOnlyRatios_pad4->SetHeader(collisionSystem2760GeV.Data());
                    legen4PadOnlyRatios_pad4->AddEntry(graphRatioEtaInvYieldFitSystPbPb2760GeV_2050,"#eta, 20#font[122]{-}50%","pf");
                    legen4PadOnlyRatios_pad4->Draw();

                ratio2DLowerPad1OnlyRatio->Draw("axis,same");
            canvasInvYieldSectionOnlyRatios4pads->Update();
            canvasInvYieldSectionOnlyRatios4pads->SaveAs(Form("%s/YieldCombinedLHC11h_DataWithModelsOnlyRatios4Pad.%s",outputDir.Data(),suffix.Data()));
//             canvasInvYieldSectionOnlyRatios4pads->SaveAs(Form("%s/YieldCombinedLHC11h_DataWithModelsOnlyRatios4Pad.%s",paperPlots.Data(),suffix.Data()));

            //column ratios to models
            Double_t arrayBoundariesX_colRatio[2];
            Double_t arrayBoundariesY_colRatio[5];
            Double_t relativeMarginsX_colRatio[3];
            Double_t relativeMarginsY_colRatio[3];
            textSizeLabelsPixel = 40;
            ReturnCorrectValuesForCanvasScaling(1000,1400, 1, 4,0.14, 0.03, 0.02,0.07,arrayBoundariesX_colRatio,arrayBoundariesY_colRatio,relativeMarginsX_colRatio,relativeMarginsY_colRatio);
            TCanvas* canvasInvYieldOnlyRatios = new TCanvas("canvasInvYieldOnlyRatios","",0,0,1000,1400);  // gives the page size

            TPad* padPi0YieldModelRatio0010 = new TPad("padPi0YieldModelRatio0010", "", arrayBoundariesX_colRatio[0], arrayBoundariesY_colRatio[1], arrayBoundariesX_colRatio[1], arrayBoundariesY_colRatio[0],-1, -1, -2);
            DrawGammaPadSettings( padPi0YieldModelRatio0010, relativeMarginsX_colRatio[0], relativeMarginsX_colRatio[2], relativeMarginsY_colRatio[0], relativeMarginsY_colRatio[1]);

            TPad* padPi0YieldModelRatio2050 = new TPad("padPi0YieldModelRatio2050", "", arrayBoundariesX_colRatio[0], arrayBoundariesY_colRatio[2], arrayBoundariesX_colRatio[1], arrayBoundariesY_colRatio[1],-1, -1, -2);
            DrawGammaPadSettings( padPi0YieldModelRatio2050, relativeMarginsX_colRatio[0], relativeMarginsX_colRatio[2], relativeMarginsY_colRatio[1], relativeMarginsY_colRatio[1]);

            TPad* padEtaYieldModelRatio0010 = new TPad("padEtaYieldModelRatio0010", "", arrayBoundariesX_colRatio[0], arrayBoundariesY_colRatio[3], arrayBoundariesX_colRatio[1], arrayBoundariesY_colRatio[2],-1, -1, -2);
            DrawGammaPadSettings( padEtaYieldModelRatio0010, relativeMarginsX_colRatio[0], relativeMarginsX_colRatio[2], relativeMarginsY_colRatio[1], relativeMarginsY_colRatio[1]);

            TPad* padEtaYieldModelRatio2050 = new TPad("padEtaYieldModelRatio2050", "", arrayBoundariesX_colRatio[0], arrayBoundariesY_colRatio[4], arrayBoundariesX_colRatio[1], arrayBoundariesY_colRatio[3],-1, -1, -2);
            DrawGammaPadSettings( padEtaYieldModelRatio2050, relativeMarginsX_colRatio[0], relativeMarginsX_colRatio[2], relativeMarginsY_colRatio[1], relativeMarginsY_colRatio[2]);

            padPi0YieldModelRatio0010->Draw();
            Double_t margincolYR = relativeMarginsX_colRatio[0]*1000;
            Double_t textsizeLabelscolYPi00010 = 0;
            Double_t textsizeFaccolYPi00010 = 0;
            if (padPi0YieldModelRatio0010->XtoPixel(padPi0YieldModelRatio0010->GetX2()) < padPi0YieldModelRatio0010->YtoPixel(padPi0YieldModelRatio0010->GetY1())){
                textsizeLabelscolYPi00010 = (Double_t)textSizeLabelsPixel/padPi0YieldModelRatio0010->XtoPixel(padPi0YieldModelRatio0010->GetX2()) ;
                textsizeFaccolYPi00010 = (Double_t)1./padPi0YieldModelRatio0010->XtoPixel(padPi0YieldModelRatio0010->GetX2()) ;
            } else {
                textsizeLabelscolYPi00010 = (Double_t)textSizeLabelsPixel/padPi0YieldModelRatio0010->YtoPixel(padPi0YieldModelRatio0010->GetY1());
                textsizeFaccolYPi00010 = (Double_t)1./padPi0YieldModelRatio0010->YtoPixel(padPi0YieldModelRatio0010->GetY1());
            }

            padPi0YieldModelRatio2050->Draw();
            Double_t textsizeLabelscolYPi02050 = 0;
            Double_t textsizeFaccolYPi02050 = 0;
            if (padPi0YieldModelRatio2050->XtoPixel(padPi0YieldModelRatio2050->GetX2()) < padPi0YieldModelRatio2050->YtoPixel(padPi0YieldModelRatio2050->GetY1())){
                textsizeLabelscolYPi02050 = (Double_t)textSizeLabelsPixel/padPi0YieldModelRatio2050->XtoPixel(padPi0YieldModelRatio2050->GetX2()) ;
                textsizeFaccolYPi02050 = (Double_t)1./padPi0YieldModelRatio2050->XtoPixel(padPi0YieldModelRatio2050->GetX2()) ;
            } else {
                textsizeLabelscolYPi02050 = (Double_t)textSizeLabelsPixel/padPi0YieldModelRatio2050->YtoPixel(padPi0YieldModelRatio2050->GetY1());
                textsizeFaccolYPi02050 = (Double_t)1./padPi0YieldModelRatio2050->YtoPixel(padPi0YieldModelRatio2050->GetY1());
            }

            padEtaYieldModelRatio0010->Draw();
            Double_t textsizeLabelscolYEta0010 = 0;
            Double_t textsizeFaccolYEta0010 = 0;
            if (padEtaYieldModelRatio0010->XtoPixel(padEtaYieldModelRatio0010->GetX2()) < padEtaYieldModelRatio0010->YtoPixel(padEtaYieldModelRatio0010->GetY1())){
                textsizeLabelscolYEta0010 = (Double_t)textSizeLabelsPixel/padEtaYieldModelRatio0010->XtoPixel(padEtaYieldModelRatio0010->GetX2()) ;
                textsizeFaccolYEta0010 = (Double_t)1./padEtaYieldModelRatio0010->XtoPixel(padEtaYieldModelRatio0010->GetX2()) ;
            } else {
                textsizeLabelscolYEta0010 = (Double_t)textSizeLabelsPixel/padEtaYieldModelRatio0010->YtoPixel(padEtaYieldModelRatio0010->GetY1());
                textsizeFaccolYEta0010 = (Double_t)1./padEtaYieldModelRatio0010->YtoPixel(padEtaYieldModelRatio0010->GetY1());
            }

            padEtaYieldModelRatio2050->Draw();
            Double_t textsizeLabelscolYEta2050 = 0;
            Double_t textsizeFaccolYEta2050 = 0;
            if (padEtaYieldModelRatio2050->XtoPixel(padEtaYieldModelRatio2050->GetX2()) < padEtaYieldModelRatio2050->YtoPixel(padEtaYieldModelRatio2050->GetY1())){
                textsizeLabelscolYEta2050 = (Double_t)textSizeLabelsPixel/padEtaYieldModelRatio2050->XtoPixel(padEtaYieldModelRatio2050->GetX2()) ;
                textsizeFaccolYEta2050 = (Double_t)1./padEtaYieldModelRatio2050->XtoPixel(padEtaYieldModelRatio2050->GetX2()) ;
            } else {
                textsizeLabelscolYEta2050 = (Double_t)textSizeLabelsPixel/padEtaYieldModelRatio2050->YtoPixel(padEtaYieldModelRatio2050->GetY1());
                textsizeFaccolYEta2050 = (Double_t)1./padEtaYieldModelRatio2050->YtoPixel(padEtaYieldModelRatio2050->GetY1());
            }

            TH2F * ratio2DPi0010 = new TH2F("ratio2DPi0010","ratio2DPi0010",1000,0.25,70.,1000,0.2,2.1);
            SetStyleHistoTH2ForGraphs(ratio2DPi0010, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Fit}", 0.85*textsizeLabelscolYPi00010, textsizeLabelscolYPi00010, 0.85*textsizeLabelscolYPi00010,textsizeLabelscolYPi00010, 1.2,0.2/(textsizeFaccolYPi00010*margincolYR), 510, 505);
            ratio2DPi0010->GetXaxis()->SetMoreLogLabels(kTRUE);
            ratio2DPi0010->GetXaxis()->SetNoExponent(kTRUE);
            ratio2DPi0010->GetXaxis()->SetRangeUser(minX,maxX);

            TH2F * ratio2DPi2050 = new TH2F("ratio2DPi2050","ratio2DPi2050",1000,0.25,70.,1000,0.2,2.1);
            SetStyleHistoTH2ForGraphs(ratio2DPi2050, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Fit}", 0.85*textsizeLabelscolYPi02050, textsizeLabelscolYPi02050, 0.85*textsizeLabelscolYPi02050,textsizeLabelscolYPi02050, 1,0.2/(textsizeFaccolYPi02050*margincolYR), 510, 505);
            ratio2DPi2050->GetXaxis()->SetMoreLogLabels(kTRUE);
            ratio2DPi2050->GetXaxis()->SetNoExponent(kTRUE);
            ratio2DPi2050->GetXaxis()->SetRangeUser(minX,maxX);

            TH2F * ratio2DEta0010 = new TH2F("ratio2DEta0010","ratio2DEta0010",1000,0.25,70.,1000,0.,5);
            SetStyleHistoTH2ForGraphs(ratio2DEta0010, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Fit}", 0.85*textsizeLabelscolYEta0010, textsizeLabelscolYEta0010, 0.85*textsizeLabelscolYEta0010,textsizeLabelscolYEta0010, 1,0.2/(textsizeFaccolYEta0010*margincolYR), 510, 505);
            ratio2DEta0010->GetXaxis()->SetMoreLogLabels(kTRUE);
            ratio2DEta0010->GetXaxis()->SetNoExponent(kTRUE);
            ratio2DEta0010->GetYaxis()->SetRangeUser(0.08,2.8);
            ratio2DEta0010->GetXaxis()->SetRangeUser(minX,maxX);

            TH2F * ratio2DEta2050 = new TH2F("ratio2DEta2050","ratio2DEta2050",1000,0.25,70.,1000,0.,5);
            SetStyleHistoTH2ForGraphs(ratio2DEta2050, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Fit}", 0.85*textsizeLabelscolYEta2050, textsizeLabelscolYEta2050, 0.85*textsizeLabelscolYEta2050,textsizeLabelscolYEta2050, 1,0.2/(textsizeFaccolYEta2050*margincolYR), 510, 505);
            ratio2DEta2050->GetXaxis()->SetMoreLogLabels(kTRUE);
            ratio2DEta2050->GetXaxis()->SetNoExponent(kTRUE);
            ratio2DEta2050->GetXaxis()->SetRangeUser(minX,maxX);
            ratio2DEta2050->GetYaxis()->SetRangeUser(0.08,2.8);

            DrawGammaSetMarkerTGraphAsym(graphRatioBegunEQPi0_0010, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
            graphRatioBegunEQPi0_0010->SetLineStyle(9);
            DrawGammaSetMarkerTGraphAsym(graphRatioBegunEQPi0_2050, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
            graphRatioBegunEQPi0_2050->SetLineStyle(9);
            DrawGammaSetMarkerTGraphAsym(graphRatioBegunEQEta_0010, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
            graphRatioBegunEQEta_0010->SetLineStyle(9);
            DrawGammaSetMarkerTGraphAsym(graphRatioBegunEQEta_2050, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
            graphRatioBegunEQEta_2050->SetLineStyle(9);

            TGraphAsymmErrors* copydummyEtaEQ= (TGraphAsymmErrors*)graphRatioBegunEQEta_0010->Clone();
            DrawGammaSetMarkerTGraphAsym(copydummyEtaEQ, 0, 0, colorBegunRatio, colorBegunRatio, 5, kTRUE, colorBegunRatio);
            copydummyEtaEQ->SetLineStyle(9);

            padPi0YieldModelRatio0010->cd();
            padPi0YieldModelRatio0010->SetLogx(1);
            ratio2DPi0010->DrawCopy();
                DrawGammaLines(minX,maxX,1., 1.,0.1,kGray);

                while(graphRatioLowPtPi0_0010->GetX()[0] < .95)
                graphRatioLowPtPi0_0010->RemovePoint(0);
                graphRatioLowPtPi0_0010->Draw("l,same");

                while(graphRatioBegunEQPi0_0010->GetX()[0] < .95)
                graphRatioBegunEQPi0_0010->RemovePoint(0);
                graphRatioBegunEQPi0_0010->Draw("l,same");

                while(graphRatioEPOSToFitPi0_0010->GetX()[0] < .95)
                graphRatioEPOSToFitPi0_0010->RemovePoint(0);
                graphRatioEPOSToFitPi0_0010->Draw("same,l,x");

                graphRatioPi0InvYieldFitStatPbPb2760GeV_0010->Draw("p,same");
                graphRatioPi0InvYieldFitSystPbPb2760GeV_0010->Draw("E2same");

                TLatex *labelSystemRatioPlots = new TLatex(0.18,0.78,collisionSystem2760GeV.Data());
                SetStyleTLatex( labelSystemRatioPlots, 0.1,4);
                labelSystemRatioPlots->Draw();
                if(thesisPlotting){
                    TLatex *thesisLabel = new TLatex(0.2,0.64,thisthesis.Data());
                    SetStyleTLatex( thesisLabel,0.12,4);
//                     thesisLabel->Draw();
                }
                TLegend* legen4RowOnlyRatios_pad1 = GetAndSetLegend2(0.7, 0.08, 0.9, 0.19, 0.85*textSizeLabelsPixel);
                legen4RowOnlyRatios_pad1->SetMargin(0.2);
        //           legen4RowOnlyRatios_pad1->SetHeader(collisionSystem2760GeV.Data());
                legen4RowOnlyRatios_pad1->AddEntry(graphRatioPi0InvYieldFitSystPbPb2760GeV_0010,"#pi^{0},  0#font[122]{-}10%","pf");
                legen4RowOnlyRatios_pad1->Draw();

                TLegend* legen4RowsOnlyRatiosTheory = GetAndSetLegend2(0.68, 0.85-3*0.08, 0.97, 0.85, 0.83* textSizeLabelsPixel);
                legen4RowsOnlyRatiosTheory->AddEntry(graphRatioCombLowPtPbPb2760GeV_0010,"NEQ SHM","l");
                legen4RowsOnlyRatiosTheory->AddEntry(copydummyEtaEQ,"EQ SHM","l");
                legen4RowsOnlyRatiosTheory->AddEntry(graphRatioEPOSToFitPbPb2760GeV_0010,"EPOS","l");
                legen4RowsOnlyRatiosTheory->Draw();

                ratio2DPi0010->Draw("axis,same");
            padPi0YieldModelRatio2050->cd();
            padPi0YieldModelRatio2050->SetLogx(1);
            ratio2DPi2050->DrawCopy();
                DrawGammaLines(minX,maxX,1., 1.,0.1,kGray);

                while(graphRatioLowPtPi0_2050->GetX()[0] < .95)
                graphRatioLowPtPi0_2050->RemovePoint(0);
                graphRatioLowPtPi0_2050->Draw("l,same");

                while(graphRatioBegunEQPi0_2050->GetX()[0] < .95)
                graphRatioBegunEQPi0_2050->RemovePoint(0);
                graphRatioBegunEQPi0_2050->Draw("l,same");

                while(graphRatioEPOSToFitPi0_2050->GetX()[0] < .95)
                graphRatioEPOSToFitPi0_2050->RemovePoint(0);
                graphRatioEPOSToFitPi0_2050->Draw("same,l,x");

                graphRatioPi0InvYieldFitStatPbPb2760GeV_2050->Draw("p,same");
                graphRatioPi0InvYieldFitSystPbPb2760GeV_2050->Draw("E2same");

                TLegend* legen4RowOnlyRatios_pad2 = GetAndSetLegend2(0.7, 0.08, 0.9, 0.19, 0.85*textSizeLabelsPixel);
                legen4RowOnlyRatios_pad2->SetMargin(0.2);
        //           legen4RowOnlyRatios_pad2->SetHeader(collisionSystem2760GeV.Data());
                legen4RowOnlyRatios_pad2->AddEntry(graphRatioPi0InvYieldFitSystPbPb2760GeV_0010,"#pi^{0},  20#font[122]{-}50%","pf");
                legen4RowOnlyRatios_pad2->Draw();

                ratio2DPi2050->Draw("axis,same");
            padEtaYieldModelRatio0010->cd();
            padEtaYieldModelRatio0010->SetLogx(1);
                ratio2DEta0010->DrawCopy();
                DrawGammaLines(minX,maxX,1., 1.,0.1,kGray);

                while(graphRatioLowPtEta_0010->GetX()[0] < .95)
                graphRatioLowPtEta_0010->RemovePoint(0);
                graphRatioLowPtEta_0010->Draw("l,same");

                while(graphRatioBegunEQEta_0010->GetX()[0] < .95)
                graphRatioBegunEQEta_0010->RemovePoint(0);
                graphRatioBegunEQEta_0010->Draw("l,same");

                while(graphRatioEPOSToFitEta_0010->GetX()[0] < .95)
                graphRatioEPOSToFitEta_0010->RemovePoint(0);
                graphRatioEPOSToFitEta_0010->Draw("same,l,x");

                graphRatioEtaInvYieldFitStatPbPb2760GeV_0010->Draw("p,same");
                graphRatioEtaInvYieldFitSystPbPb2760GeV_0010->Draw("E2same");

                TLegend* legen4RowOnlyRatios_pad3 = GetAndSetLegend2(0.7, 0.08, 0.9, 0.19, 0.85*textSizeLabelsPixel);
                legen4RowOnlyRatios_pad3->SetMargin(0.2);
        //           legen4RowOnlyRatios_pad3->SetHeader(collisionSystem2760GeV.Data());
                legen4RowOnlyRatios_pad3->AddEntry(graphRatioEtaInvYieldFitSystPbPb2760GeV_0010,"#eta,  0#font[122]{-}10%","pf");
                legen4RowOnlyRatios_pad3->Draw();

                ratio2DEta0010->Draw("axis,same");
            padEtaYieldModelRatio2050->cd();
            padEtaYieldModelRatio2050->SetLogx(1);
            ratio2DEta2050->DrawCopy();
                DrawGammaLines(minX,maxX,1., 1.,0.1,kGray);

                while(graphRatioLowPtEta_2050->GetX()[0] < .95)
                graphRatioLowPtEta_2050->RemovePoint(0);
                graphRatioLowPtEta_2050->Draw("l,same");

                while(graphRatioBegunEQEta_2050->GetX()[0] < .95)
                graphRatioBegunEQEta_2050->RemovePoint(0);
                graphRatioBegunEQEta_2050->Draw("l,same");

                while(graphRatioEPOSToFitEta_2050->GetX()[0] < .95)
                graphRatioEPOSToFitEta_2050->RemovePoint(0);
                graphRatioEPOSToFitEta_2050->Draw("same,l,x");

                graphRatioEtaInvYieldFitStatPbPb2760GeV_2050->Draw("p,same");
                graphRatioEtaInvYieldFitSystPbPb2760GeV_2050->Draw("E2same");

                TLegend* legen4RowOnlyRatios_pad4 = GetAndSetLegend2(0.7, 0.29, 0.9, 0.38, 0.85*textSizeLabelsPixel);
                legen4RowOnlyRatios_pad4->SetMargin(0.2);
        //           legen4RowOnlyRatios_pad4->SetHeader(collisionSystem2760GeV.Data());
                legen4RowOnlyRatios_pad4->AddEntry(graphRatioEtaInvYieldFitSystPbPb2760GeV_0010,"#eta,  20#font[122]{-}50%","pf");
                legen4RowOnlyRatios_pad4->Draw();

                ratio2DEta2050->Draw("axis,same");
            canvasInvYieldOnlyRatios->Update();
            canvasInvYieldOnlyRatios->SaveAs(Form("%s/YieldCombinedLHC11h_DataWithModelsOnlyRatiosColumn.%s",paperPlots.Data(),suffix.Data()));
        }
    }


    Double_t arrayBoundariesX1_ratioCh[3];
    Double_t arrayBoundariesY1_ratioCh[2];
    Double_t relativeMarginsX_ratioCh[3];
    Double_t relativeMarginsY_ratioCh[3];
    ReturnCorrectValuesForCanvasScaling(1000,500, 2, 1,0.06, 0.005, 0.005,0.09,arrayBoundariesX1_ratioCh,arrayBoundariesY1_ratioCh,relativeMarginsX_ratioCh,relativeMarginsY_ratioCh);
    TCanvas* canvas2PartNeutralToChargedRatio = new TCanvas("canvas2PartNeutralToChargedRatio","",0,0,1000,500);  // gives the page size

    TPad* pad2PartNeutralToChargedRatio1 = new TPad("pad2PartNeutralToChargedRatio1", "", arrayBoundariesX1_ratioCh[0], arrayBoundariesY1_ratioCh[1], arrayBoundariesX1_ratioCh[1], arrayBoundariesY1_ratioCh[0],-1, -1, -2);
    DrawGammaPadSettings( pad2PartNeutralToChargedRatio1, relativeMarginsX_ratioCh[0], relativeMarginsX_ratioCh[1], relativeMarginsY_ratioCh[0], relativeMarginsY_ratioCh[2]);
    pad2PartNeutralToChargedRatio1->Draw();

    TPad* pad2PartNeutralToChargedRatio2 = new TPad("pad2PartNeutralToChargedRatio2", "", arrayBoundariesX1_ratioCh[1], arrayBoundariesY1_ratioCh[1], arrayBoundariesX1_ratioCh[2], arrayBoundariesY1_ratioCh[0],-1, -1, -2);
    DrawGammaPadSettings( pad2PartNeutralToChargedRatio2, relativeMarginsX_ratioCh[1], relativeMarginsX_ratioCh[2], relativeMarginsY_ratioCh[0], relativeMarginsY_ratioCh[2]);
    pad2PartNeutralToChargedRatio2->Draw();

    Double_t margin = relativeMarginsX_ratioCh[0]*0.8*1000;
    Double_t textsizeLabels1 = 0;
    Double_t textsizeFac1 = 0;
    Double_t textsizeLabels2 = 0;
    Double_t textsizeFac2 = 0;

    ReturnCorrectValuesTextSize(pad2PartNeutralToChargedRatio1,textsizeLabels1, textsizeFac1, 22, margin);
    ReturnCorrectValuesTextSize(pad2PartNeutralToChargedRatio2,textsizeLabels2, textsizeFac2, 22, margin);

    TH2F * histo2DRatioChargeToNeutral1 = new TH2F("histo2DRatioChargeToNeutral1","histo2DRatioChargeToNeutral1",1000,0.25,70.,1000,0.2,4. );
    if(meson.CompareTo("Pi0")==0){
        SetStyleHistoTH2ForGraphs(histo2DRatioChargeToNeutral1, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}",0.85*textsizeLabels1, textsizeLabels1,
                                      0.85*textsizeLabels1, textsizeLabels1, 0.8,0.25/(textsizeFac1*margin), 512, 505);
    } else if(meson.CompareTo("Eta")==0){
        SetStyleHistoTH2ForGraphs(histo2DRatioChargeToNeutral1, "#it{p}_{T} (GeV/#it{c})","#eta/K^{#pm}",0.85*textsizeLabels1, textsizeLabels1,
                                      0.85*textsizeLabels1, textsizeLabels1, 0.8,0.25/(textsizeFac1*margin), 512, 505);
    }
    histo2DRatioChargeToNeutral1->GetYaxis()->SetRangeUser(0.,2.1);
    histo2DRatioChargeToNeutral1->GetXaxis()->SetLabelOffset(-0.0105);
    histo2DRatioChargeToNeutral1->GetXaxis()->SetRangeUser(minX,maxX);

    TH2F* histo2DRatioChargeToNeutral2 = new TH2F("histo2DRatioChargeToNeutral2","histo2DRatioChargeToNeutral2",1000,0.25,70.,1000,0.2,4.    );
    if(meson.CompareTo("Pi0")==0){
        SetStyleHistoTH2ForGraphs(histo2DRatioChargeToNeutral2, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textsizeLabels2, textsizeLabels2,
                                  0.85*textsizeLabels2, textsizeLabels2, 0.8,0.25/(textsizeFac2*margin), 512, 505);
    } else if(meson.CompareTo("Eta")==0){
        SetStyleHistoTH2ForGraphs(histo2DRatioChargeToNeutral2, "#it{p}_{T} (GeV/#it{c})","#eta/K^{#pm}", 0.85*textsizeLabels2, textsizeLabels2,
                                  0.85*textsizeLabels2, textsizeLabels2, 0.8,0.25/(textsizeFac2*margin), 512, 505);
    }
    histo2DRatioChargeToNeutral2->GetXaxis()->SetLabelOffset(-0.0105);
    histo2DRatioChargeToNeutral2->GetYaxis()->SetRangeUser(0.,2.1);
    histo2DRatioChargeToNeutral2->GetXaxis()->SetRangeUser(minX,maxX);

    pad2PartNeutralToChargedRatio1->cd();
    pad2PartNeutralToChargedRatio1->SetLogx();
    histo2DRatioChargeToNeutral1->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioNeutralToChargedComb_0010, markerStyleComb , markerSizeComb , colorComb, colorComb);
        graphRatioNeutralToChargedComb_0010->Draw("E1psame");

        TLatex *labelPi0CompChargedPionsPbPb0010 = new TLatex(0.16,0.9,collisionSystemPbPb0010.Data());
        SetStyleTLatex( labelPi0CompChargedPionsPbPb0010, 0.85*textsizeLabels1,4);
        labelPi0CompChargedPionsPbPb0010->Draw();
        DrawGammaLines(minX,maxX, 1, 1 ,1, kGray, 2);

        TLegend* legendPi0CombChargedPionsPbPb0010 = new TLegend(0.12,0.78,0.5,0.88);//0.12,0.68,0.95,0.88);
        legendPi0CombChargedPionsPbPb0010->SetFillColor(0);
        legendPi0CombChargedPionsPbPb0010->SetLineColor(0);
//         legendPi0CombChargedPionsPbPb0010->SetNColumns(2);
        legendPi0CombChargedPionsPbPb0010->SetTextSize(FontSize);
        if(meson.CompareTo("Pi0")==0){
          legendPi0CombChargedPionsPbPb0010->AddEntry(graphRatioNeutralToChargedComb_0010,"#pi^{0}/#pi^{#pm} Comb","p");
        } else if(meson.CompareTo("Eta")==0){
          legendPi0CombChargedPionsPbPb0010->AddEntry(graphRatioNeutralToChargedComb_0010,"#eta/K^{#pm} Comb","p");
          TLine *linemTOnset = new TLine(4.,4.,0.1,1.5);
          linemTOnset->SetLineColor(kGray);
          DrawGammaLines(4., 4 , 0.2, 1.5 ,1, kGray, 1);
          legendPi0CombChargedPionsPbPb0010->AddEntry(linemTOnset,"#it{m}_{T} scaling onset","l");
        }

        legendPi0CombChargedPionsPbPb0010->Draw();

    histo2DRatioChargeToNeutral1->Draw("axis,same");
    pad2PartNeutralToChargedRatio1->Update();
    pad2PartNeutralToChargedRatio2->cd();
    pad2PartNeutralToChargedRatio2->SetLogx();
    histo2DRatioChargeToNeutral2->DrawCopy();

        TGraphErrors *graphRatioNeutralToChargedComb_0010_clone = (TGraphErrors*)graphRatioNeutralToChargedComb_0010->Clone();
        DrawGammaSetMarkerTGraphErr(graphRatioNeutralToChargedComb_0010_clone, markerStyleComb , markerSizeComb , kGray+1, kGray+1);
        graphRatioNeutralToChargedComb_0010_clone->Draw("E1psame");

        if(meson.CompareTo("Pi0")==0){
            DrawGammaSetMarkerTGraphErr(graphRatioPi0PCMToChargedPion_0010, markerStyleDet[0], markerSizeDet[0]*0.4 , colorDet[0], colorDet[0]);
            graphRatioPi0PCMToChargedPion_0010->Draw("E1psame");

            DrawGammaSetMarkerTGraphErr(graphRatioPi0PHOSToChargedPion_0010, markerStyleDet[1] , markerSizeDet[1]*0.4 , colorDet[1], colorDet[1]);
            graphRatioPi0PHOSToChargedPion_0010->Draw("E1psame");

            DrawGammaSetMarkerTGraphErr(graphRatioPi0EMCalToChargedPion_0010, markerStyleDet[2] , markerSizeDet[2]*0.4, colorDet[2], colorDet[2]);
            graphRatioPi0EMCalToChargedPion_0010->Draw("E1psame");

        } else if(meson.CompareTo("Eta")==0){
            DrawGammaSetMarkerTGraphErr(graphRatioEtaPCMToChargedKaon_0010, markerStyleDet[0], markerSizeDet[0]*0.4, colorDet[0], colorDet[0]);
            graphRatioEtaPCMToChargedKaon_0010->Draw("E1psame");

            DrawGammaSetMarkerTGraphErr(graphRatioEtaEMCalToChargedKaon_0010, markerStyleDet[2] , markerSizeDet[2]*0.4, colorDet[2], colorDet[2]);
            graphRatioEtaEMCalToChargedKaon_0010->Draw("E1psame");

        }
        DrawGammaLines(minX,maxX, 1, 1 ,1,kGray, 2);

        TLegend* legendPi0SingleMeasChargedPionsPbPb0010 = new TLegend(0.05,0.78,0.95,0.92);
        legendPi0SingleMeasChargedPionsPbPb0010->SetFillColor(0);
        legendPi0SingleMeasChargedPionsPbPb0010->SetLineColor(0);
        legendPi0SingleMeasChargedPionsPbPb0010->SetNColumns(2);
        legendPi0SingleMeasChargedPionsPbPb0010->SetTextSize(FontSize);
        if(meson.CompareTo("Pi0")==0){
//           legendPi0SingleMeasChargedPionsPbPb0010->AddEntry(graphRatioNeutralToChargedComb_0010,"#pi^{0}/#pi^{#pm} Comb","p");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry(graphRatioNeutralToChargedComb_0010_clone,"#pi^{0}/#pi^{#pm} Comb","p");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry(graphRatioPi0PCMToChargedPion_0010 ,"#pi^{0}/#pi^{#pm} PCM","p");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry((TObject*)0,"","");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry(graphRatioPi0PHOSToChargedPion_0010,"#pi^{0}/#pi^{#pm} PHOS","p");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry((TObject*)0,"","");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry(graphRatioPi0EMCalToChargedPion_0010,"#pi^{0}/#pi^{#pm} EMCal","p");
        } else if(meson.CompareTo("Eta")==0){
//           legendPi0SingleMeasChargedPionsPbPb0010->AddEntry(graphRatioNeutralToChargedComb_0010,"#eta/K^{#pm} Comb","p");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry(graphRatioNeutralToChargedComb_0010_clone,"#eta/K^{#pm} Comb","p");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry(graphRatioEtaPCMToChargedKaon_0010,"#eta/K^{#pm} PCM","p");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry((TObject*)0,"","");
          legendPi0SingleMeasChargedPionsPbPb0010->AddEntry(graphRatioEtaEMCalToChargedKaon_0010,"#eta/K^{#pm} EMCal","p");
        }
        legendPi0SingleMeasChargedPionsPbPb0010->Draw();

    histo2DRatioChargeToNeutral2->Draw("axis,same");
    pad2PartNeutralToChargedRatio2->Update();
    canvas2PartNeutralToChargedRatio->Update();
    canvas2PartNeutralToChargedRatio->SaveAs(Form("%s/%s_RatioNeutralToCharged_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));

    pad2PartNeutralToChargedRatio1->cd();
    pad2PartNeutralToChargedRatio1->SetLogx();
    histo2DRatioChargeToNeutral1->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioNeutralToChargedComb_2050, markerStyleComb , markerSizeComb , colorComb, colorComb);
        graphRatioNeutralToChargedComb_2050->Draw("E1psame");

        TLatex *labelPi0CompChargedPionsPbPb2050 = new TLatex(0.16,0.9,collisionSystemPbPb2050.Data());
        SetStyleTLatex( labelPi0CompChargedPionsPbPb2050, 0.85*textsizeLabels1,4);
        labelPi0CompChargedPionsPbPb2050->Draw();
        DrawGammaLines(minX,maxX, 1, 1 ,1, kGray, 2);

        TLegend* legendPi0CombChargedPionsPbPb2050 = new TLegend(0.12,0.78,0.5,0.88);//0.12,0.68,0.95,0.88);
        legendPi0CombChargedPionsPbPb2050->SetFillColor(0);
        legendPi0CombChargedPionsPbPb2050->SetLineColor(0);
//         legendPi0CombChargedPionsPbPb2050->SetNColumns(2);
        legendPi0CombChargedPionsPbPb2050->SetTextSize(FontSize);
        if(meson.CompareTo("Pi0")==0){
          legendPi0CombChargedPionsPbPb2050->AddEntry(graphRatioNeutralToChargedComb_2050,"#pi^{0}/#pi^{#pm} Comb","p");
        } else if(meson.CompareTo("Eta")==0){
          legendPi0CombChargedPionsPbPb2050->AddEntry(graphRatioNeutralToChargedComb_2050,"#eta/K^{#pm} Comb","p");
          TLine *linemTOnset = new TLine(4.,4.,0.1,1.5);
          linemTOnset->SetLineColor(kGray);
          DrawGammaLines(4., 4 , 0.2, 1.5 ,1, kGray, 1);
          legendPi0CombChargedPionsPbPb2050->AddEntry(linemTOnset,"#it{m}_{T} scaling onset","l");
        }

        legendPi0CombChargedPionsPbPb2050->Draw();

    histo2DRatioChargeToNeutral1->Draw("axis,same");
    pad2PartNeutralToChargedRatio1->Update();
    pad2PartNeutralToChargedRatio2->cd();
    pad2PartNeutralToChargedRatio2->SetLogx();
    histo2DRatioChargeToNeutral2->DrawCopy();

        TGraphErrors *graphRatioNeutralToChargedComb_2050_clone = (TGraphErrors*)graphRatioNeutralToChargedComb_2050->Clone();
        DrawGammaSetMarkerTGraphErr(graphRatioNeutralToChargedComb_2050_clone, markerStyleComb , markerSizeComb , kGray+1, kGray+1);
        graphRatioNeutralToChargedComb_2050_clone->Draw("E1psame");

        if(meson.CompareTo("Pi0")==0){
            DrawGammaSetMarkerTGraphErr(graphRatioPi0PCMToChargedPion_2050, markerStyleDet[0], markerSizeDet[0]*0.4 , colorDet[0], colorDet[0]);
            graphRatioPi0PCMToChargedPion_2050->Draw("E1psame");

            DrawGammaSetMarkerTGraphErr(graphRatioPi0EMCalToChargedPion_2050, markerStyleDet[2] , markerSizeDet[2]*0.4, colorDet[2], colorDet[2]);
            graphRatioPi0EMCalToChargedPion_2050->Draw("E1psame");

        } else if(meson.CompareTo("Eta")==0){
            DrawGammaSetMarkerTGraphErr(graphRatioEtaPCMToChargedKaon_2050, markerStyleDet[0], markerSizeDet[0]*0.4, colorDet[0], colorDet[0]);
            graphRatioEtaPCMToChargedKaon_2050->Draw("E1psame");

            DrawGammaSetMarkerTGraphErr(graphRatioEtaEMCalToChargedKaon_2050, markerStyleDet[2] , markerSizeDet[2]*0.4, colorDet[2], colorDet[2]);
            graphRatioEtaEMCalToChargedKaon_2050->Draw("E1psame");

        }
        DrawGammaLines(minX,maxX, 1, 1 ,1,kGray, 2);

        TLegend* legendPi0SingleMeasChargedPionsPbPb2050 = new TLegend(0.05,0.78,0.95,0.92);
        legendPi0SingleMeasChargedPionsPbPb2050->SetFillColor(0);
        legendPi0SingleMeasChargedPionsPbPb2050->SetLineColor(0);
        legendPi0SingleMeasChargedPionsPbPb2050->SetNColumns(2);
        legendPi0SingleMeasChargedPionsPbPb2050->SetTextSize(FontSize);
        if(meson.CompareTo("Pi0")==0){
//           legendPi0SingleMeasChargedPionsPbPb2050->AddEntry(graphRatioNeutralToChargedComb_2050,"#pi^{0}/#pi^{#pm} Comb","p");
          legendPi0SingleMeasChargedPionsPbPb2050->AddEntry(graphRatioNeutralToChargedComb_2050_clone,"#pi^{0}/#pi^{#pm} Comb","p");
          legendPi0SingleMeasChargedPionsPbPb2050->AddEntry(graphRatioPi0PCMToChargedPion_2050 ,"#pi^{0}/#pi^{#pm} PCM","p");
          legendPi0SingleMeasChargedPionsPbPb2050->AddEntry((TObject*)0,"","");
          legendPi0SingleMeasChargedPionsPbPb2050->AddEntry(graphRatioPi0EMCalToChargedPion_2050,"#pi^{0}/#pi^{#pm} EMCal","p");
        } else if(meson.CompareTo("Eta")==0){
//           legendPi0SingleMeasChargedPionsPbPb2050->AddEntry(graphRatioNeutralToChargedComb_2050,"#eta/K^{#pm} Comb","p");
          legendPi0SingleMeasChargedPionsPbPb2050->AddEntry(graphRatioNeutralToChargedComb_2050_clone,"#eta/K^{#pm} Comb","p");
          legendPi0SingleMeasChargedPionsPbPb2050->AddEntry(graphRatioEtaPCMToChargedKaon_2050,"#eta/K^{#pm} PCM","p");
          legendPi0SingleMeasChargedPionsPbPb2050->AddEntry((TObject*)0,"","");
          legendPi0SingleMeasChargedPionsPbPb2050->AddEntry(graphRatioEtaEMCalToChargedKaon_2050,"#eta/K^{#pm} EMCal","p");
        }
        legendPi0SingleMeasChargedPionsPbPb2050->Draw();

    histo2DRatioChargeToNeutral2->Draw("axis,same");
    pad2PartNeutralToChargedRatio2->Update();
    canvas2PartNeutralToChargedRatio->Update();
//     canvas2PartNeutralToChargedRatio->SaveAs(Form("%s/%s_RatioNeutralToCharged_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    delete pad2PartNeutralToChargedRatio1;
    delete pad2PartNeutralToChargedRatio2;
    delete canvas2PartNeutralToChargedRatio;

    TCanvas* canvasTemp = new TCanvas("canvasTemp","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasTemp, 0.09, 0.02, 0.035, 0.09);
    canvasTemp->SetLogx();

      TH2F * histo2DTemp = new TH2F("histo2DTemp","histo2DTemp",1000,0.25,70.,1000,0,4.);
      if(meson.CompareTo("Pi0")==0){
          SetStyleHistoTH2ForGraphs(histo2DTemp, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.035,0.04, 0.035,0.04, 1,1.);
      } else if(meson.CompareTo("Eta")==0){
          SetStyleHistoTH2ForGraphs(histo2DTemp, "#it{p}_{T} (GeV/#it{c})","#eta/K^{#pm}", 0.035,0.04, 0.035,0.04, 1,1.);
      }

    histo2DTemp->GetYaxis()->SetRangeUser(0.2,2.1);
    histo2DTemp->GetXaxis()->SetRangeUser(minX,maxX);
      histo2DTemp->DrawCopy();

        if(meson.CompareTo("Pi0")==0){
            DrawGammaSetMarkerTGraphErr(graphRatioPi0PCMToChargedPion_2050, markerStyleDet[0], markerSizeDet[0]*0.8, colorDet[0], colorDet[0]);
            graphRatioPi0PCMToChargedPion_2050->Draw("E1psame");

            DrawGammaSetMarkerTGraphErr(graphRatioPi0EMCalToChargedPion_2050, markerStyleDet[2] , markerSizeDet[2]*0.8, colorDet[2], colorDet[2]);
            graphRatioPi0EMCalToChargedPion_2050->Draw("E1psame");

        } else if(meson.CompareTo("Eta")==0){
            DrawGammaSetMarkerTGraphErr(graphRatioEtaPCMToChargedKaon_2050, markerStyleDet[0], markerSizeDet[0]*0.8, colorDet[0], colorDet[0]);
            graphRatioEtaPCMToChargedKaon_2050->Draw("E1psame");

            DrawGammaSetMarkerTGraphErr(graphRatioEtaEMCalToChargedKaon_2050, markerStyleDet[2] , markerSizeDet[2]*0.8, colorDet[2], colorDet[2]);
            graphRatioEtaEMCalToChargedKaon_2050->Draw("E1psame");

        }
        DrawGammaLines(minX,maxX , 1, 1 ,1,kGray, 2);

        TLegend* legendSingleMeasChargedPionsPbPb2050 = new TLegend(0.12,0.75,0.52,0.92);
        legendSingleMeasChargedPionsPbPb2050->SetFillColor(0);
        legendSingleMeasChargedPionsPbPb2050->SetLineColor(0);
//         legendSingleMeasChargedPionsPbPb2050->SetNColumns(2);
        legendSingleMeasChargedPionsPbPb2050->SetTextSize(FontSize);
        legendSingleMeasChargedPionsPbPb2050->SetHeader(collisionSystemPbPb2040.Data());
        if(meson.CompareTo("Pi0")==0){
          legendSingleMeasChargedPionsPbPb2050->AddEntry(graphRatioPi0PCMToChargedPion_2050 ,"#pi^{0}/#pi^{#pm} PCM","p");
          legendSingleMeasChargedPionsPbPb2050->AddEntry(graphRatioPi0EMCalToChargedPion_2050,"#pi^{0}/#pi^{#pm} EMCal","p");
          legendSingleMeasChargedPionsPbPb2050->AddEntry((TObject*)0,"(20#font[122]{-}50% for EMCal)","");
        } else if(meson.CompareTo("Eta")==0){
          legendSingleMeasChargedPionsPbPb2050->AddEntry(graphRatioEtaPCMToChargedKaon_2050,"#eta/K^{#pm} PCM","p");
          legendSingleMeasChargedPionsPbPb2050->AddEntry(graphRatioEtaEMCalToChargedKaon_2050,"#eta/K^{#pm} EMCal","p");
          legendSingleMeasChargedPionsPbPb2050->AddEntry((TObject*)0,"(20#font[122]{-}50% for EMCal)","");
        }
        legendSingleMeasChargedPionsPbPb2050->Draw();

    canvasTemp->SaveAs(Form("%s/%s_RatioNeutralToCharged_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    //**********************************************************************************************************************//
    //***********************************************     Raa combination    ***********************************************//
    //**********************************************************************************************************************//
    // RAA calc and plotting
    Bool_t quiet = kTRUE;
    if(meson.CompareTo("Pi0")==0){

          CalcRaa(    graphInvSectionPCMStatPi02760GeVforRAA, graphInvSectionPCMSysPi02760GeVforRAA,graphInvSectionCombStatPi02760GeVforRAA, fitInvCrossSectionTsallisPi0Comb2760GeV,
                          graphPCMInvYieldStatPbPb2760GeVYShifted_0010, graphPCMInvYieldSysPbPb2760GeVYShifted_0010,
                          &graphRAAPCM0010, &graphRAASysPCM0010,
                          nColl0010, nCollErr0010,"Pi0",8.,0,"h",quiet);

          CalcRaa(    graphInvSectionPHOSStatPi02760GeVforRAA, graphInvSectionPHOSSysPi02760GeVforRAA,graphInvSectionCombStatPi02760GeVforRAA, fitInvCrossSectionTsallisPi0Comb2760GeV,
                          graphPHOSInvYieldStatPbPb2760GeVYShifted_0010, graphPHOSInvYieldSysPbPb2760GeVYShifted_0010,
                          &graphRAAPHOS0010, &graphRAASysPHOS0010,
                          nColl0010, nCollErr0010,"Pi0", 12.,0,"h",quiet);

          CalcRaa(    graphInvSectionPCMStatPi02760GeVforRAA, graphInvSectionPCMSysPi02760GeVforRAA,graphInvSectionCombStatPi02760GeVforRAA, fitInvCrossSectionTsallisPi0Comb2760GeV,
                          graphPCMInvYieldStatPbPb2760GeVYShifted_2050, graphPCMInvYieldSysPbPb2760GeVYShifted_2050,
                          &graphRAAPCM2050, &graphRAASysPCM2050,
                          nColl2050, nCollErr2050,"Pi0",8.,0.,"h",quiet);

          CalcRaa(    graphInvSectionEMCalStatPi02760GeVforRAA, graphInvSectionEMCalSysPi02760GeVforRAA,graphInvSectionCombStatPi02760GeVforRAA, fitInvCrossSectionTsallisPi0Comb2760GeV,
                          graphEMCalInvYieldStatPbPb2760GeVYShifted_0010, graphEMCalInvYieldSysPbPb2760GeVYShifted_0010,
                          &graphRAAEMCal0010, &graphRAASysEMCal0010,
                          nColl0010, nCollErr0010,"Pi0",20.,0,"h",quiet);

          CalcRaa(    graphInvSectionEMCalStatPi02760GeVforRAA, graphInvSectionEMCalSysPi02760GeVforRAA,graphInvSectionCombStatPi02760GeVforRAA, fitInvCrossSectionTsallisPi0Comb2760GeV,
                          graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, graphEMCalInvYieldSysPbPb2760GeVYShifted_2050,
                          &graphRAAEMCal2050, &graphRAASysEMCal2050,
                          nColl2050, nCollErr2050,"Pi0",20.,0.,"h",quiet);

            if((thesisPlotting || PaperPi0)){

                while(graphRAAPCM0010->GetX()[0] < /* 0.8 */ 1.)
                graphRAAPCM0010->RemovePoint(0);
                while(graphRAASysPCM0010->GetX()[0] < /* 0.8 */ 1.)
                graphRAASysPCM0010->RemovePoint(0);
                while(graphRAAPCM2050->GetX()[0] < /* 0.8 */ 1.)
                graphRAAPCM2050->RemovePoint(0);
                while(graphRAASysPCM2050->GetX()[0] < /* 0.8 */ 1.)
                graphRAASysPCM2050->RemovePoint(0);
            }

            histoRAAStatPCM0010       = GraphAsymErrorsToHist_withErrors(graphRAAPCM0010,"histoRAAStatPCM0010");
            histoRAAStatPHOS0010       = GraphAsymErrorsToHist_withErrors(graphRAAPHOS0010,"histoRAAStatPHOS0010");
            histoRAAStatEMCal0010       = GraphAsymErrorsToHist_withErrors(graphRAAEMCal0010,"histoRAAStatEMCal0010");
            histoRAAStatPCM2050       = GraphAsymErrorsToHist_withErrors(graphRAAPCM2050,"histoRAAStatPCM2050");
            histoRAAStatEMCal2050       = GraphAsymErrorsToHist_withErrors(graphRAAEMCal2050,"histoRAAStatEMCal2050");


            Double_t nCollCSC = nColl0010/nColl2050;
            CalcRaa(    graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, graphEMCalInvYieldSysPbPb2760GeVYShifted_2050,graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, fitInvCrossSectionTsallisPi0Comb2760GeV,
                            graphEMCalInvYieldStatPbPb2760GeVYShifted_0010, graphEMCalInvYieldSysPbPb2760GeVYShifted_0010,
                            &graphRCSCEMCal, &graphRCSCSysEMCal,
                            nCollCSC, nCollErr2050,"Pi0",20.,0.,"h",quiet);

//           CalcRcp(    graphEMCalInvYieldStatPbPb2760GeVYShifted_0010, graphEMCalInvYieldSysPbPb2760GeVYShifted_0010,
//                       graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, graphEMCalInvYieldSysPbPb2760GeVYShifted_2050,
//                           &graphRCSCEMCal2050, &graphRCSCSysEMCal2050,
//                           nColl0010, nColl2050,"Pi0");


//           return;
    } else if(meson.CompareTo("Eta")==0){

      //for Preliminary:
//       graphPCMInvYieldSysPbPb2760GeVYShifted_0010->RemovePoint(0);
//       graphPCMInvYieldSysPbPb2760GeVYShifted_2050->RemovePoint(0);
      ////////////////
          CalcRaa(    graphInvSectionPCMStatEta2760GeVforRAA, graphInvSectionPCMSysEta2760GeVforRAA,graphInvSectionCombStatEta2760GeVforRAA, fitInvCrossSectionTsallisEtaComb2760GeV,
                          graphPCMInvYieldStatPbPb2760GeVYShifted_2050, graphPCMInvYieldSysPbPb2760GeVYShifted_2050,
                          &graphRAAPCM2050, &graphRAASysPCM2050,
                          nColl2050, nCollErr2050,"Eta",6.,0.,"h",quiet);

          CalcRaa(    graphInvSectionEMCalStatEta2760GeVforRAA, graphInvSectionEMCalSysEta2760GeVforRAA,graphInvSectionCombStatEta2760GeVforRAA, fitInvCrossSectionTsallisEtaComb2760GeV,
                          graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, graphEMCalInvYieldSysPbPb2760GeVYShifted_2050,
                          &graphRAAEMCal2050, &graphRAASysEMCal2050,
                          nColl2050, nCollErr2050,"Eta",20.,0.,"h",quiet);

          CalcRaa(    graphInvSectionPCMStatEta2760GeVforRAA, graphInvSectionPCMSysEta2760GeVforRAA,graphInvSectionCombStatEta2760GeVforRAA, fitInvCrossSectionTsallisEtaComb2760GeV,
                          graphPCMInvYieldStatPbPb2760GeVYShifted_0010, graphPCMInvYieldSysPbPb2760GeVYShifted_0010,
                          &graphRAAPCM0010, &graphRAASysPCM0010,
                          nColl0010, nCollErr0010,"Eta",6.,0,"h",quiet);

          CalcRaa(    graphInvSectionEMCalStatEta2760GeVforRAA, graphInvSectionEMCalSysEta2760GeVforRAA,graphInvSectionCombStatEta2760GeVforRAA, fitInvCrossSectionTsallisEtaComb2760GeV,
                          graphEMCalInvYieldStatPbPb2760GeVYShifted_0010, graphEMCalInvYieldSysPbPb2760GeVYShifted_0010,
                          &graphRAAEMCal0010, &graphRAASysEMCal0010,
                          nColl0010, nCollErr0010,"Eta",20.,0,"h",quiet);

          histoRAAStatPCM0010       = GraphAsymErrorsToHist_withErrors(graphRAAPCM0010,"histoRAAStatPCM0010");
          histoRAAStatEMCal0010       = GraphAsymErrorsToHist_withErrors(graphRAAEMCal0010,"histoRAAStatEMCal0010");
          histoRAAStatPCM2050       = GraphAsymErrorsToHist_withErrors(graphRAAPCM2050,"histoRAAStatPCM2050");
          histoRAAStatEMCal2050       = GraphAsymErrorsToHist_withErrors(graphRAAEMCal2050,"histoRAAStatEMCal2050");

          Double_t nCollCSC = nColl0010/nColl2050;
          CalcRaa(    graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, graphEMCalInvYieldSysPbPb2760GeVYShifted_2050,graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, fitInvCrossSectionTsallisPi0Comb2760GeV,
                          graphEMCalInvYieldStatPbPb2760GeVYShifted_0010, graphEMCalInvYieldSysPbPb2760GeVYShifted_0010,
                          &graphRCSCEMCal, &graphRCSCSysEMCal,
                          nCollCSC, nCollErr2050,"Eta",20.,0.,"h",quiet);

    }


//     graphRAAPCM0010->Print();
//     graphRAAPHOS0010->Print();
//     graphRAAEMCal0010->Print();
//     graphRAASysPCM0010->Print();
//     graphRAASysPHOS0010->Print();
//     graphRAASysEMCal0010->Print();
//     return;

    TCanvas* canvasRAAMeasurements = new TCanvas("canvasRAAMeasurements","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRAAMeasurements, 0.09, 0.02, 0.035, 0.09);
//  canvasRAAMeasurements->SetLogx();

      TH2F * histo2DRAADummy = new TH2F("histo2DRAADummy","histo2DRAADummy",1000,0.,21.,1000,0,1.2);
      SetStyleHistoTH2ForGraphs(histo2DRAADummy, "#it{p}_{T} (GeV/#it{c})","R_{AA}", 0.032,0.04, 0.04,0.04, 1,1.);
      histo2DRAADummy->DrawCopy();

      DrawGammaLines(0., 21. , 1., 1.,0.5,   kGray);

      DrawGammaSetMarkerTGraphAsym(graphRAASysPCM0010, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
      graphRAASysPCM0010->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphRAAPCM0010, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAPCM0010);
      graphRAAPCM0010->Draw("p,same");

      if(meson.CompareTo("Pi0")==0){
          DrawGammaSetMarkerTGraphAsym(graphRAASysPHOS0010, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
          graphRAASysPHOS0010->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(graphRAAPHOS0010, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAPHOS0010);
          graphRAAPHOS0010->Draw("p,same");
      }

      DrawGammaSetMarkerTGraphAsym(graphRAASysEMCal0010,markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
      graphRAASysEMCal0010->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphRAAEMCal0010, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAEMCal0010);
      graphRAAEMCal0010->Draw("p,same");


        // ****************************** Definition of the Legend ******************************************
        // **************** Row def ************************
        Double_t rowsLegendOnlyPi0Ratio2[5]      = {0.75,0.7,0.66,0.62,0.59};
        Double_t rowsLegendOnlyPi0RatioAbs2[5]   = {0.5,.85,0.8,0.75,0.9};
        Double_t columnsLegendOnlyPi0Ratio2[3]   = {0.14,0.37, 0.42};
        Double_t columnsLegendOnlyPi0RatioAbs2[3]= {7.,7.2,8.3};
        Double_t lengthBox2                      = 0.5/2;
        Double_t heightBox2                      = 0.08/4;
        // ****************** first Column **************************************************
        TLatex *textPCMOnlyLHC11h;
        if(thesisPlotting)
            textPCMOnlyLHC11h = new TLatex(columnsLegendOnlyPi0Ratio2[0],rowsLegendOnlyPi0Ratio2[1],"PCM (this thesis)");
        else
            textPCMOnlyLHC11h = new TLatex(columnsLegendOnlyPi0Ratio2[0],rowsLegendOnlyPi0Ratio2[1],"PCM");

        SetStyleTLatex( textPCMOnlyLHC11h, 0.85*40,4);
        textPCMOnlyLHC11h->SetTextFont(43);
        textPCMOnlyLHC11h->Draw();
        TLatex *textPHOSOnlyLHC11h = new TLatex(columnsLegendOnlyPi0Ratio2[0],rowsLegendOnlyPi0Ratio2[3],"PHOS");
        SetStyleTLatex( textPHOSOnlyLHC11h, 0.85*40,4);
        textPHOSOnlyLHC11h->SetTextFont(43);
        if(meson.CompareTo("Pi0")==0){
          textPHOSOnlyLHC11h->Draw();
        }
        TLatex *textEMCalOnlyLHC11h = new TLatex(columnsLegendOnlyPi0Ratio2[0],rowsLegendOnlyPi0Ratio2[2],"EMCal");//" (*)");
        SetStyleTLatex( textEMCalOnlyLHC11h,  0.85*40,4);
        textEMCalOnlyLHC11h->SetTextFont(43);
        textEMCalOnlyLHC11h->Draw();
        // ****************** second Column *************************************************
        TLatex *textStatOnlyLHC11h = new TLatex(columnsLegendOnlyPi0Ratio2[1],rowsLegendOnlyPi0Ratio2[0] ,"stat");
        SetStyleTLatex( textStatOnlyLHC11h, 0.85*40,4);
        textStatOnlyLHC11h->SetTextFont(43);
        textStatOnlyLHC11h->Draw();
        TLatex *textSysOnlyLHC11h = new TLatex(columnsLegendOnlyPi0Ratio2[2] ,rowsLegendOnlyPi0Ratio2[0],"syst");
        SetStyleTLatex( textSysOnlyLHC11h, 0.85*40,4);
        textSysOnlyLHC11h->SetTextFont(43);
        textSysOnlyLHC11h->Draw();

        TMarker* markerPCMOnlyLHC11h = CreateMarkerFromGraph(graphRAASysPCM0010,columnsLegendOnlyPi0Ratio2[1] ,rowsLegendOnlyPi0Ratio2[1],1);
        markerPCMOnlyLHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs2[1] ,rowsLegendOnlyPi0RatioAbs2[1]);
        TBox* boxPCMOnlyPi0 = CreateBoxFromGraph(graphRAASysPCM0010, columnsLegendOnlyPi0RatioAbs2[2]-1.5*lengthBox2 , rowsLegendOnlyPi0RatioAbs2[1]- heightBox2,
                                                        columnsLegendOnlyPi0RatioAbs2[2]+ 2*lengthBox2, rowsLegendOnlyPi0RatioAbs2[1]+ heightBox2);
        boxPCMOnlyPi0->Draw("l");

        if(meson.CompareTo("Pi0")==0){
            TMarker* markerPHOSOnlyLHC11h = CreateMarkerFromGraph(graphRAASysPHOS0010, columnsLegendOnlyPi0Ratio2[1] ,rowsLegendOnlyPi0Ratio2[3],1);
            markerPHOSOnlyLHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs2[1] ,rowsLegendOnlyPi0RatioAbs2[3]);
            TBox* boxPHOSOnly = CreateBoxFromGraph(graphRAASysPHOS0010, columnsLegendOnlyPi0RatioAbs2[2]-1.5*lengthBox2 , rowsLegendOnlyPi0RatioAbs2[3]- heightBox2,
                                                            columnsLegendOnlyPi0RatioAbs2[2]+ 2*lengthBox2, rowsLegendOnlyPi0RatioAbs2[3]+ heightBox2);
            boxPHOSOnly->Draw("l");
        }

        TMarker* markerEMCalOnlyLHC11h = CreateMarkerFromGraph(graphRAASysEMCal0010, columnsLegendOnlyPi0Ratio2[1] ,rowsLegendOnlyPi0Ratio2[2],1);
        markerEMCalOnlyLHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs2[1] ,rowsLegendOnlyPi0RatioAbs2[2]);
        TBox* boxEMCalOnly = CreateBoxFromGraph(graphRAASysEMCal0010, columnsLegendOnlyPi0RatioAbs2[2]-1.5*lengthBox2 , rowsLegendOnlyPi0RatioAbs2[2]- heightBox2,
                                                        columnsLegendOnlyPi0RatioAbs2[2]+ 2*lengthBox2, rowsLegendOnlyPi0RatioAbs2[2]+ heightBox2);
        boxEMCalOnly->Draw("l");

        labelEnergyIndMeasRAA0010->Draw();
        labelDetSysIndMeasRAA->Draw();

        TLatex *textEMCalNoteRaa = new TLatex(0.65,0.75,"(*) ALICE preliminary input");
        SetStyleTLatex( textEMCalNoteRaa,  0.8*40,4);
        textEMCalNoteRaa->SetTextFont(43);
        TLatex *textEMCalNoteRaa1 = new TLatex(0.65,0.71,"arXiv:1609.06106");
        SetStyleTLatex( textEMCalNoteRaa1,  0.8*40,4);
        textEMCalNoteRaa1->SetTextFont(43);
        TLatex *textEMCalNoteRaa2 = new TLatex(0.65,0.67,"analysis by A. Morreale");
        SetStyleTLatex( textEMCalNoteRaa2,  0.8*40,4);
        textEMCalNoteRaa2->SetTextFont(43);
//         textEMCalNoteRaa->Draw();
//         textEMCalNoteRaa1->Draw();
//         textEMCalNoteRaa2->Draw();

        histo2DRAADummy->Draw("axis,same");
    canvasRAAMeasurements->Update();
    canvasRAAMeasurements->SaveAs(Form("%s/%s_RAA_IndividualMeasLHC11h_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRAAMeasurements->SaveAs(Form("%s/%s_RAA_IndividualMeasLHC11h_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasRAAMeasurements->cd();
      histo2DRAADummy->DrawCopy();

      DrawGammaLines(0., 21 , 1., 1.,0.5,   kGray);

      DrawGammaSetMarkerTGraphAsym(graphRAASysPCM2050, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
      graphRAASysPCM2050->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphRAAPCM2050, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAPCM2050);
      graphRAAPCM2050->Draw("p,same");

      DrawGammaSetMarkerTGraphAsym(graphRAASysEMCal2050,markerStyleDet[2] ,markerStyleDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
      graphRAASysEMCal2050->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphRAAEMCal2050, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAEMCal2050);
      graphRAAEMCal2050->Draw("p,same");

        labelEnergyIndMeasRAA2050->Draw();
        labelDetSysIndMeasRAA->Draw();

      textPCMOnlyLHC11h->Draw();
      textEMCalOnlyLHC11h->Draw();

      textStatOnlyLHC11h->Draw();
      textSysOnlyLHC11h->Draw();

        markerPCMOnlyLHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs2[1] ,rowsLegendOnlyPi0RatioAbs2[1]);
        markerEMCalOnlyLHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs2[1] ,rowsLegendOnlyPi0RatioAbs2[2]);
        boxPCMOnlyPi0->Draw("l");
        boxEMCalOnly->Draw("l");

//         textEMCalNoteRaa->Draw();
//         textEMCalNoteRaa1->Draw();
//         textEMCalNoteRaa2->Draw();

        histo2DRAADummy->Draw("axis,same");

    canvasRAAMeasurements->Update();
    canvasRAAMeasurements->SaveAs(Form("%s/%s_RAA_IndividualMeasLHC11h_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRAAMeasurements->SaveAs(Form("%s/%s_RAA_IndividualMeasLHC11h_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));



    // direct combination procedure
    TGraphAsymmErrors *graphDirectCombRAATotPbPb2760GeV_0010 = NULL;
    TGraphAsymmErrors *graphDirectCombRAAStatPbPb2760GeV_0010 = NULL;
    TGraphAsymmErrors *graphDirectCombRAASysPbPb2760GeV_0010 = NULL;
    TGraphAsymmErrors *graphDirectCombRAATotPbPb2760GeV_2050 = NULL;
    TGraphAsymmErrors *graphDirectCombRAAStatPbPb2760GeV_2050 = NULL;
    TGraphAsymmErrors *graphDirectCombRAASysPbPb2760GeV_2050 = NULL;
//     if(meson.CompareTo("Pi0")==0){
//       graphDirectCombRAATotPbPb2760GeV_0010 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeV_0010,     graphPCMPi0RAASysPbPb2760GeV_0010,
//                                                                  graphEMCalPi0RAAStatPbPb2760GeV_0010,   graphEMCalPi0RAASysPbPb2760GeV_0010,
//                                                                  graphDirectCombRAAStatPbPb2760GeV_0010, graphDirectCombRAASysPbPb2760GeV_0010,
//                                                                  xPtLimitsPi0, 24, 0, 1, 16, kFALSE);
//
//       graphDirectCombRAATotPbPb2760GeV_2050 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeV_2050,     graphPCMPi0RAASysPbPb2760GeV_2050,
//                                                                  graphEMCalPi0RAAStatPbPb2760GeV_2050,   graphEMCalPi0RAASysPbPb2760GeV_2050,
//                                                                  graphDirectCombRAAStatPbPb2760GeV_2050, graphDirectCombRAASysPbPb2760GeV_2050,
//                                                                  xPtLimitsPi0, 24, 0, 1, 16, kFALSE);
//
//     } else if(meson.CompareTo("Eta")==0){
//
//       graphDirectCombRAATotPbPb2760GeV_0010 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeV_0010,     graphPCMEtaRAASysPbPb2760GeV_0010,
//                                                                  graphEMCalEtaRAAStatPbPb2760GeV_0010,   graphEMCalEtaRAASysPbPb2760GeV_0010,
//                                                                  graphDirectCombRAAStatPbPb2760GeV_0010, graphDirectCombRAASysPbPb2760GeV_0010,
//                                                                  xPtLimitsEta, 14, 0, 2, 6, kFALSE);
//
//       graphDirectCombRAATotPbPb2760GeV_2050 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeV_2050,     graphPCMEtaRAASysPbPb2760GeV_2050,
//                                                                  graphEMCalEtaRAAStatPbPb2760GeV_2050,   graphEMCalEtaRAASysPbPb2760GeV_2050,
//                                                                  graphDirectCombRAAStatPbPb2760GeV_2050, graphDirectCombRAASysPbPb2760GeV_2050,
//                                                                  xPtLimitsEta, 14, 0, 2, 6, kFALSE);
//
//     }


    for (Int_t i = 0; i< 11; i++){

        statErrorCollectionforRAALHC11h_0010[i] = NULL;
        statErrorCollectionforRAALHC11h_2050[i] = NULL;

        sysErrorCollectionforRAALHC11h_0010[i] = NULL;
        sysErrorCollectionforRAALHC11h_2050[i] = NULL;

    }
    if(meson.CompareTo("Pi0")==0){

        statErrorCollectionforRAALHC11h_0010[0] = (TH1D*)histoRAAStatPCM0010->Clone("statErrPCMPi0forRAA_0010");
        statErrorCollectionforRAALHC11h_0010[1] = (TH1D*)histoRAAStatPHOS0010->Clone("statErrPCMPi0forRAA_0010");
        statErrorCollectionforRAALHC11h_0010[2] = (TH1D*)histoRAAStatEMCal0010->Clone("statErrEMCalPi0forRAA_0010");

        statErrorCollectionforRAALHC11h_2050[0] = (TH1D*)histoRAAStatPCM2050->Clone("statErrPCMPi0forRAA_2050");
        statErrorCollectionforRAALHC11h_2050[2] = (TH1D*)histoRAAStatEMCal2050->Clone("statErrEMCalPi0forRAA_2050");


        sysErrorCollectionforRAALHC11h_0010[0] = (TGraphAsymmErrors*)graphRAASysPCM0010->Clone("sysErrPCMPi0forRAA_0010");
        sysErrorCollectionforRAALHC11h_0010[1] = (TGraphAsymmErrors*)graphRAASysPHOS0010->Clone("sysErrPCMPi0forRAA_0010");
        sysErrorCollectionforRAALHC11h_0010[2] = (TGraphAsymmErrors*)graphRAASysEMCal0010->Clone("sysErrEMCalPi0forRAA_0010");

        sysErrorCollectionforRAALHC11h_2050[0] = (TGraphAsymmErrors*)graphRAASysPCM2050->Clone("sysErrPCMPi0forRAA_2050");
        sysErrorCollectionforRAALHC11h_2050[2] = (TGraphAsymmErrors*)graphRAASysEMCal2050->Clone("sysErrEMCalPi0forRAA_2050");

    } else if(meson.CompareTo("Eta")==0) {

        statErrorCollectionforRAALHC11h_0010[0] = (TH1D*)histoRAAStatPCM0010->Clone("statErrPCMEtaforRAA_0010");
        statErrorCollectionforRAALHC11h_0010[2] = (TH1D*)histoRAAStatEMCal0010->Clone("statErrEMCalEtaforRAA_0010");

        statErrorCollectionforRAALHC11h_2050[0] = (TH1D*)histoRAAStatPCM2050->Clone("statErrPCMEtaforRAA_2050");
        statErrorCollectionforRAALHC11h_2050[2] = (TH1D*)histoRAAStatEMCal2050->Clone("statErrEMCalEtaforRAA_2050");


        sysErrorCollectionforRAALHC11h_0010[0] = (TGraphAsymmErrors*)graphRAASysPCM0010->Clone("sysErrPCMEtaforRAA_0010");
        sysErrorCollectionforRAALHC11h_0010[2] = (TGraphAsymmErrors*)graphRAASysEMCal0010->Clone("sysErrEMCalEtaforRAA_0010");

        sysErrorCollectionforRAALHC11h_2050[0] = (TGraphAsymmErrors*)graphRAASysPCM2050->Clone("sysErrPCMEtaforRAA_2050");
        sysErrorCollectionforRAALHC11h_2050[2] = (TGraphAsymmErrors*)graphRAASysEMCal2050->Clone("sysErrEMCalEtaforRAA_2050");

    }

    cout << " \n\nCombining RAA for " << meson.Data() << endl;

    cout << __LINE__ << endl;
    Int_t npoitRAA = 23;
    if((thesisPlotting || PaperPi0) && meson.CompareTo("Pi0")==0){
        offSetsPi0RAA[0] = 4;
        offSetsPi0RAA[1] = 4;
        offSetsPi0RAA[2] = 15;
        offSetsPi0RAASys[0] = 4;
        offSetsPi0RAASys[1] = 4;
        offSetsPi0RAASys[2] = 15;

    }

    if(meson.CompareTo("Pi0")==0){
        graphCombRAATotPbPb2760GeV_0010 = CombinePtPointsSpectraFullCorrMat( statErrorCollectionforRAALHC11h_0010,  sysErrorCollectionforRAALHC11h_0010,
                                                                                                  xPtLimitsPi0, npoitRAA, offSetsPi0RAA, offSetsPi0RAASys,
                                                                                                  graphCombRAAStatPbPb2760GeV_0010, graphCombRAASysPbPb2760GeV_0010,
                                                                                                    Form("%s/0010LHC11h_WeightingRAA%s.dat",outputDir.Data(),meson.Data()));
        graphCombRAATotPbPb2760GeV_2050 = CombinePtPointsSpectraFullCorrMat( statErrorCollectionforRAALHC11h_2050,  sysErrorCollectionforRAALHC11h_2050,
                                                                                                  xPtLimitsPi0, npoitRAA, offSetsPi0RAA, offSetsPi0RAASys,
                                                                                                  graphCombRAAStatPbPb2760GeV_2050, graphCombRAASysPbPb2760GeV_2050,
                                                                                                    Form("%s/2050LHC11h_WeightingRAA%s.dat",outputDir.Data(),meson.Data()));
        while(graphCombRAATotPbPb2760GeV_0010->GetY()[0] < 1e-14){
          graphCombRAATotPbPb2760GeV_0010->RemovePoint(0);
          graphCombRAAStatPbPb2760GeV_0010->RemovePoint(0);
          graphCombRAASysPbPb2760GeV_0010->RemovePoint(0);
          graphCombRAATotPbPb2760GeV_2050->RemovePoint(0);
          graphCombRAAStatPbPb2760GeV_2050->RemovePoint(0);
          graphCombRAASysPbPb2760GeV_2050->RemovePoint(0);
        }
    } else if(meson.CompareTo("Eta")==0){
        graphCombRAATotPbPb2760GeV_0010 = CombinePtPointsSpectraFullCorrMat( statErrorCollectionforRAALHC11h_0010,  sysErrorCollectionforRAALHC11h_0010,
                                                                                                  xPtLimitsEta, 13, offSetsEtaRAA, offSetsEtaRAASys,
                                                                                                  graphCombRAAStatPbPb2760GeV_0010, graphCombRAASysPbPb2760GeV_0010,
                                                                                                    Form("%s/0010LHC11h_WeightingRAA%s.dat",outputDir.Data(),meson.Data()));

        graphCombRAATotPbPb2760GeV_2050 = CombinePtPointsSpectraFullCorrMat( statErrorCollectionforRAALHC11h_2050,  sysErrorCollectionforRAALHC11h_2050,
                                                                                                  xPtLimitsEta, 13, offSetsEtaRAA, offSetsEtaRAASys,
                                                                                                  graphCombRAAStatPbPb2760GeV_2050, graphCombRAASysPbPb2760GeV_2050,
                                                                                                    Form("%s/2050LHC11h_WeightingRAA%s.dat",outputDir.Data(),meson.Data()));

        while(graphCombRAATotPbPb2760GeV_0010->GetY()[0] < 1e-14){
          graphCombRAATotPbPb2760GeV_0010->RemovePoint(0);
          graphCombRAAStatPbPb2760GeV_0010->RemovePoint(0);
          graphCombRAASysPbPb2760GeV_0010->RemovePoint(0);
          graphCombRAATotPbPb2760GeV_2050->RemovePoint(0);
          graphCombRAAStatPbPb2760GeV_2050->RemovePoint(0);
          graphCombRAASysPbPb2760GeV_2050->RemovePoint(0);
        }
    }

//     graphCombRAATotPbPb2760GeV_0010->Print();
//     graphCombRAAStatPbPb2760GeV_0010->Print();
//     graphCombRAASysPbPb2760GeV_0010->Print();
//     graphCombRAATotPbPb2760GeV_2050->Print();
//     graphCombRAAStatPbPb2760GeV_2050->Print();
//     graphCombRAASysPbPb2760GeV_2050->Print();
//     for(Int_t i=0; i<graphCombRAASysPbPb2760GeV_0010->GetN()-1;i++){
//         cout << "x: " << graphCombRAASysPbPb2760GeV_0010->GetX()[i] << " y: " << graphCombRAASysPbPb2760GeV_0010->GetY()[i] << " err: " << graphCombRAASysPbPb2760GeV_0010->GetEYlow()[i] << " % err: " << graphCombRAASysPbPb2760GeV_0010->GetEYlow()[i]/graphCombRAASysPbPb2760GeV_0010->GetY()[i]*100 << endl;
//     }
//     for(Int_t i=0; i<graphCombRAASysPbPb2760GeV_2050->GetN()-1;i++){
//         cout << "x: " << graphCombRAASysPbPb2760GeV_2050->GetX()[i] << " y: " << graphCombRAASysPbPb2760GeV_2050->GetY()[i] << " err: " << graphCombRAASysPbPb2760GeV_2050->GetEYlow()[i] << " % err: " << graphCombRAASysPbPb2760GeV_2050->GetEYlow()[i]/graphCombRAASysPbPb2760GeV_2050->GetY()[i]*100 << endl;
//     }
//     return;


    TCanvas* canvasRAAcombo = new TCanvas("canvasRAAcombo","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRAAcombo, 0.08, 0.02, 0.035, 0.09);

        TH2F * histo2DRAAcombo;
        histo2DRAAcombo = new TH2F("histo2DRAAcombo","histo2DRAAcombo",11000,0.,70.,1000,0.,2.);
        SetStyleHistoTH2ForGraphs(histo2DRAAcombo, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,.9);
        histo2DRAAcombo->GetYaxis()->SetRangeUser(0.01,1.4);
        histo2DRAAcombo->GetXaxis()->SetRangeUser(0.,21);
        histo2DRAAcombo->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeV_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
        graphCombRAASysPbPb2760GeV_0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeV_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeV_0010);
        graphCombRAAStatPbPb2760GeV_0010->Draw("p,same");

        DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeV_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050, 2, kTRUE);
        graphCombRAASysPbPb2760GeV_2050->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeV_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeV_2050);
        graphCombRAAStatPbPb2760GeV_2050->Draw("p,same");

        labelSystRaa->Draw();

        TLegend* legendRAAcombo2 = new TLegend(0.595,0.73,0.91,0.88);
        legendRAAcombo2->SetFillColor(0);
        legendRAAcombo2->SetLineColor(0);
        legendRAAcombo2->SetTextFont(42);
        legendRAAcombo2->SetMargin(0.17);
        legendRAAcombo2->SetTextSize(FontSize);
        legendRAAcombo2->SetHeader(collisionSystem2760GeV.Data());
        legendRAAcombo2->AddEntry(graphCombRAASysPbPb2760GeV_0010,Form("%s ",cent0010.Data()),"fp");
        legendRAAcombo2->AddEntry(graphCombRAASysPbPb2760GeV_2050,Form("%s ",cent2050.Data()),"fp");
        legendRAAcombo2->Draw();

        DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);
        boxErrorNorm0010->Draw();
        boxErrorNorm2050->Draw();

    canvasRAAcombo->SaveAs(Form("%s/%s_RAAcombined_DataOnly.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRAAcombo->SaveAs(Form("%s/%s_RAAcombined_DataOnly.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    canvasRAAcombo->cd();
    histo2DRAAcombo->DrawCopy();

    // Direct combination
//       DrawGammaSetMarkerTGraphAsym(graphDirectCombRAASysPbPb2760GeV_0010, markerStyle0010, markerSizeComb, kBlack , kBlack, 2, kTRUE);
//       graphDirectCombRAASysPbPb2760GeV_0010->Draw("E2same");
//       DrawGammaSetMarkerTGraphAsym(graphDirectCombRAAStatPbPb2760GeV_0010, markerStyle0010, markerSizeComb, kBlack , kBlack);
//       if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphDirectCombRAAStatPbPb2760GeV_0010);
//       graphDirectCombRAAStatPbPb2760GeV_0010->Draw("p,same");
//
//       DrawGammaSetMarkerTGraphAsym(graphDirectCombRAASysPbPb2760GeV_2050, markerStyle2050, markerSizeComb, kGray+2, kGray+2, 2, kTRUE);
//       graphDirectCombRAASysPbPb2760GeV_2050->Draw("E2same");
//       DrawGammaSetMarkerTGraphAsym(graphDirectCombRAAStatPbPb2760GeV_2050, markerStyle2050, markerSizeComb, kGray+2 , kGray+2);
//       if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphDirectCombRAAStatPbPb2760GeV_2050);
//       graphDirectCombRAAStatPbPb2760GeV_2050->Draw("p,same");

    // Calculation raa -> combination
      DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeV_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
      graphCombRAASysPbPb2760GeV_0010->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeV_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeV_0010);
      graphCombRAAStatPbPb2760GeV_0010->Draw("p,same");

      DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeV_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050, 2, kTRUE);
      graphCombRAASysPbPb2760GeV_2050->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeV_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeV_2050);
      graphCombRAAStatPbPb2760GeV_2050->Draw("p,same");

        TLegend* legendRAAcombo3 = new TLegend(0.5,0.7,0.91,0.88);
        legendRAAcombo3->SetFillColor(0);
        legendRAAcombo3->SetLineColor(0);
        legendRAAcombo3->SetTextFont(42);
        legendRAAcombo3->SetMargin(0.17);
        legendRAAcombo3->SetTextSize(FontSize);
        legendRAAcombo3->SetHeader(collisionSystem2760GeV.Data());
        legendRAAcombo3->AddEntry(graphCombRAASysPbPb2760GeV_0010,Form("%s ",cent0010.Data()),"fp");
//         legendRAAcombo3->AddEntry(graphDirectCombRAASysPbPb2760GeV_0010,Form("%s direct combo",cent0010.Data()),"fp");
        legendRAAcombo3->AddEntry(graphCombRAASysPbPb2760GeV_2050,Form("%s ",cent2050.Data()),"fp");
//         legendRAAcombo3->AddEntry(graphDirectCombRAASysPbPb2760GeV_2050,Form("%s direct combo",cent2050.Data()),"fp");
        legendRAAcombo3->Draw();

    canvasRAAcombo->SaveAs(Form("%s/%s_RAA_CombAndDirectComb.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//     canvasRAAcombo->SaveAs(Form("%s/%s_RAA_CombAndDirectComb.%s",paperPlots.Data(),meson.Data(),suffix.Data()));

    if(MesonInput){

      graphRAAPi0StatBothMeson_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphRAAPi0CombPbPb2760GeVStatErr_0010");
      graphRAAPi0SysBothMeson_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphRAAPi0CombPbPb2760GeVSysErr_0010");
      graphRAAPi0StatBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphRAAPi0CombPbPb2760GeVStatErr_2050");
      graphRAAPi0SysBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphRAAPi0CombPbPb2760GeVSysErr_2050");

      graphRAAEtaStatBothMeson_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphRAAEtaCombPbPb2760GeVStatErr_0010");
      graphRAAEtaSysBothMeson_0010= (TGraphAsymmErrors*)MesonInput->Get("graphRAAEtaCombPbPb2760GeVSysErr_0010");
      graphRAAEtaStatBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphRAAEtaCombPbPb2760GeVStatErr_2050");
      graphRAAEtaSysBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphRAAEtaCombPbPb2760GeVSysErr_2050");

      graphPi0RCSCSysEMCal = (TGraphAsymmErrors*)MesonInput->Get("graphPi0RCSCSysEMCal");
      graphPi0RCSCEMCal = (TGraphAsymmErrors*)MesonInput->Get("graphPi0RCSCEMCal");
      graphEtaRCSCSysEMCal = (TGraphAsymmErrors*)MesonInput->Get("graphEtaRCSCSysEMCal");
      graphEtaRCSCEMCal = (TGraphAsymmErrors*)MesonInput->Get("graphEtaRCSCEMCal");

      if(graphRAAPi0StatBothMeson_0010 && graphRAAEtaStatBothMeson_0010 && graphRAAPi0StatBothMeson_2050 && graphRAAEtaStatBothMeson_2050){

        TCanvas* canvasRAAcomboPi0andEta = new TCanvas("canvasRAAcomboPi0andEta","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasRAAcomboPi0andEta, 0.08, 0.02, 0.035, 0.09);
    //      canvasRAAcomboPi0andEta->SetLogx();

          TH2F * histo2DRAAcomboPi0andEta = new TH2F("histo2DRAAcomboPi0andEta","histo2DRAAcomboPi0andEta",11000,0.,21.,1000,-0.5,2.);
          SetStyleHistoTH2ForGraphs(histo2DRAAcomboPi0andEta, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, .95,0.9);
          histo2DRAAcomboPi0andEta->GetYaxis()->SetRangeUser(0.,1.38);
//           histo2DRAAcomboPi0andEta->GetXaxis()->SetRangeUser(0.0001,21);
          histo2DRAAcomboPi0andEta->Draw("copy");

          boxErrorNorm0010Only->Draw();

//         TGraphAsymmErrors* graphCombRAAStatPbPbPbPb2760GeVPi0_0010 = NULL;
//         TGraphAsymmErrors* graphDirectCombRAASysPbPb2760GeVPi0_0010 = NULL;
//         TGraphAsymmErrors* graphCombRAATotPbPb2760GeVPi0_0010 = NULL;
//
//         TGraphAsymmErrors* graphCombRAAStatPbPbPbPb2760GeVPi0_2050 = NULL;
//         TGraphAsymmErrors* graphDirectCombRAASysPbPb2760GeVPi0_2050 = NULL;
//         TGraphAsymmErrors* graphCombRAATotPbPb2760GeVPi0_2050 = NULL;
//
//         TGraphAsymmErrors* graphCombRAAStatPbPbPbPb2760GeVEta_0010 = NULL;
//         TGraphAsymmErrors* graphDirectCombRAASysPbPb2760GeVEta_0010 = NULL;
//         TGraphAsymmErrors* graphCombRAATotPbPb2760GeVEta_0010 = NULL;
//
//         TGraphAsymmErrors* graphCombRAAStatPbPbPbPb2760GeVEta_2050 = NULL;
//         TGraphAsymmErrors* graphDirectCombRAASysPbPb2760GeVEta_2050 = NULL;
//         TGraphAsymmErrors* graphCombRAATotPbPb2760GeVEta_2050 = NULL;

//      TGraphAsymmErrors *graphChargedKaonRAAStatandSyst0010 = CalculateCombinedSysAndStatError(graphChargedKaonRAA0010,graphChargedKaonRAASys0010);
//      DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAAStatandSyst0010, 25,1.6, colorCharged,colorCharged, 0.1, kFALSE);
        DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0010, 27,1.5, colorCharged,colorCharged, 0.1, kFALSE);
        DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0010,27,1.5, colorCharged,colorCharged, 1, kTRUE);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAA0010);
        graphChargedKaonRAASys0010->Draw("2same");
        graphChargedKaonRAA0010->Draw("p,same");
//      graphChargedKaonRAAStatandSyst0010->Draw("p,same");
//
//         //=============Combined Pi0
//         DrawGammaSetMarkerTGraphAsym(graphDirectCombRAASysPbPb2760GeVPi0_0010, 24, markerSizeComb, colorPCM0010, colorPCM0010, widthLinesBoxes, kTRUE);
//         graphDirectCombRAASysPbPb2760GeVPi0_0010->Draw("E2same");
//         DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPbPbPb2760GeVPi0_0010, 24, markerSizeComb, colorPCM0010, colorPCM0010);
//         if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPbPbPb2760GeVPi0_0010);
//         graphCombRAAStatPbPbPbPb2760GeVPi0_0010->Draw("p,same");
//
//         //=============Combined Eta
//         DrawGammaSetMarkerTGraphAsym(graphDirectCombRAASysPbPb2760GeVEta_0010, markerStyle0010, markerSizeComb, colorPCM0010, colorPCM0010, widthLinesBoxes, kTRUE);
//         graphDirectCombRAASysPbPb2760GeVEta_0010->Draw("E2same");
//         DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPbPbPb2760GeVEta_0010, markerStyle0010, markerSizeComb, colorPCM0010 , colorPCM0010);
//         if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPbPbPb2760GeVEta_0010);
//         graphCombRAAStatPbPbPbPb2760GeVEta_0010->Draw("p,same");

        DrawGammaSetMarkerTGraphAsym(graphRAAPi0SysBothMeson_0010, 24, markerSizeComb, colorPCM0010, colorPCM0010, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRAAPi0StatBothMeson_0010, 24, markerSizeComb, colorPCM0010, colorPCM0010);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAPi0StatBothMeson_0010);
        graphRAAPi0SysBothMeson_0010->Draw("E2same");
        graphRAAPi0StatBothMeson_0010->Draw("p,same");

        DrawGammaSetMarkerTGraphAsym(graphRAAEtaSysBothMeson_0010, markerStyle0010, markerSizeComb, colorPCM0010, colorPCM0010, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRAAEtaStatBothMeson_0010, markerStyle0010, markerSizeComb, colorPCM0010 , colorPCM0010);
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAEtaStatBothMeson_0010);
        graphRAAEtaSysBothMeson_0010->Draw("E2same");
        graphRAAEtaStatBothMeson_0010->Draw("p,same");

//         TLegend* legendRAAcomboPi0andEta0010 = new TLegend(0.5,0.73,0.9,0.93);
//         TLegend* legendRAAcomboPi0andEta0010 = new TLegend(0.54,0.78,0.9,0.88);
//         legendRAAcomboPi0andEta0010->SetFillColor(0);
//         legendRAAcomboPi0andEta0010->SetLineColor(0);
//         legendRAAcomboPi0andEta0010->SetTextFont(42);
// //      legendRAAcomboPi0andEta0010->SetMargin(0.17);
//         legendRAAcomboPi0andEta0010->SetNColumns(3);
//         legendRAAcomboPi0andEta0010->SetTextSize(FontSize);
// //         legendRAAcomboPi0andEta0010->SetHeader(Form("%s",collisionSystemPbPb0010.Data()));
//         legendRAAcomboPi0andEta0010->AddEntry(graphRAAPi0SysBothMeson_0010,"#pi^{0}","pf");
//         legendRAAcomboPi0andEta0010->AddEntry((TObject*)0,"","");
//         legendRAAcomboPi0andEta0010->AddEntry(graphChargedKaonRAA0010,"K^{#pm}","fp");
//         legendRAAcomboPi0andEta0010->AddEntry((TObject*)0,"","");
//         legendRAAcomboPi0andEta0010->AddEntry(graphRAAEtaSysBothMeson_0010,"#eta","fp");
// //         legendRAAcomboPi0andEta0010->AddEntry((TObject*)0,"","");
//         legendRAAcomboPi0andEta0010->AddEntry((TObject*)0,"PLB 736 (2014)","");
//         legendRAAcomboPi0andEta0010->Draw();

          TLegend* legendRAAcomboPi0andEta0010 = new TLegend(0.54,0.79,0.77,0.88);
          legendRAAcomboPi0andEta0010->SetFillColor(0);
          legendRAAcomboPi0andEta0010->SetLineColor(0);
          legendRAAcomboPi0andEta0010->SetTextFont(42);
          legendRAAcomboPi0andEta0010->SetTextSize(0.037);
          legendRAAcomboPi0andEta0010->AddEntry(graphRAAPi0SysBothMeson_0010,"#pi^{0}","pf");
          legendRAAcomboPi0andEta0010->AddEntry(graphRAAEtaSysBothMeson_0010,"#eta","fp");
          legendRAAcomboPi0andEta0010->Draw();
          TLegend* legendRAAKaon0010 = new TLegend(0.68,0.79,0.88,0.88);
          legendRAAKaon0010->SetFillColor(0);
          legendRAAKaon0010->SetLineColor(0);
          legendRAAKaon0010->SetTextFont(42);
          legendRAAKaon0010->SetTextSize(0.037);
          legendRAAKaon0010->AddEntry(graphChargedKaonRAA0010,"K^{#pm}","fp");
          legendRAAKaon0010->AddEntry((TObject*)0,"PLB 736 (2014)","");
          legendRAAKaon0010->Draw();

        TLatex *labelEnergyRAAcomboPi0andEta0010 = new TLatex(0.54,0.9,collisionSystemPbPb0010.Data());
        SetStyleTLatex( labelEnergyRAAcomboPi0andEta0010, 0.035,4);
        labelEnergyRAAcomboPi0andEta0010->Draw();
        DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);

      canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEta_0010.%s",outputDir.Data(),suffix.Data()));

      canvasRAAcomboPi0andEta->cd();
      histo2DRAAcomboPi0andEta->Draw("copy");
      boxErrorNorm2050EPOnly->Draw();

//           DrawGammaSetMarkerTGraphAsym(graphDirectCombRAASysPbPb2760GeVPi0_2050, 25, markerSizeComb, colorEMCal2050, colorEMCal2050, widthLinesBoxes, kTRUE);
//           graphDirectCombRAASysPbPb2760GeVPi0_2050->Draw("E2same");
//           DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPbPbPb2760GeVPi0_2050, 25, markerSizeComb, colorEMCal2050, colorEMCal2050);
//           if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPbPbPb2760GeVPi0_2050);
//           graphCombRAAStatPbPbPbPb2760GeVPi0_2050->Draw("p,same");
//
//           DrawGammaSetMarkerTGraphAsym(graphDirectCombRAASysPbPb2760GeVEta_2050, markerStyle2050, markerSizeComb, colorEMCal2050, colorEMCal2050, widthLinesBoxes, kTRUE);
//           graphDirectCombRAASysPbPb2760GeVEta_2050->Draw("E2same");
//           DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPbPbPb2760GeVEta_2050, markerStyle2050, markerSizeComb, colorEMCal2050 , colorEMCal2050);
//           if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPbPbPb2760GeVEta_2050);
//           graphCombRAAStatPbPbPbPb2760GeVEta_2050->Draw("p,same");

          DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA2040, 27,1.5, colorCharged,colorCharged, 0.1, kFALSE);
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAA2040);
          DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys2040,27,1.5, colorCharged,colorCharged, 1, kTRUE);
          graphChargedKaonRAASys2040->Draw("2same");
          graphChargedKaonRAA2040->Draw("p,same");

          DrawGammaSetMarkerTGraphAsym(graphRAAPi0SysBothMeson_2050, 25, markerSizeComb, colorEMCal2050, colorEMCal2050, widthLinesBoxes, kTRUE);
          graphRAAPi0SysBothMeson_2050->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(graphRAAPi0StatBothMeson_2050, 25, markerSizeComb, colorEMCal2050, colorEMCal2050);
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAPi0StatBothMeson_2050);
          graphRAAPi0StatBothMeson_2050->Draw("p,same");

          DrawGammaSetMarkerTGraphAsym(graphRAAEtaSysBothMeson_2050, markerStyle2050, markerSizeComb, colorEMCal2050, colorEMCal2050, widthLinesBoxes, kTRUE);
          graphRAAEtaSysBothMeson_2050->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(graphRAAEtaStatBothMeson_2050, markerStyle2050, markerSizeComb, colorEMCal2050 , colorEMCal2050);
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAEtaStatBothMeson_2050);
          graphRAAEtaStatBothMeson_2050->Draw("p,same");

          TLegend* legendRAAcomboPi0andEta2050 = new TLegend(0.57,0.75,0.8,0.88);
          legendRAAcomboPi0andEta2050->SetFillColor(0);
          legendRAAcomboPi0andEta2050->SetLineColor(0);
//           legendRAAcomboPi0andEta2050->SetNColumns(2);
//           legendRAAcomboPi0andEta2050->SetMargin(0.17);
          legendRAAcomboPi0andEta2050->SetTextFont(42);
          legendRAAcomboPi0andEta2050->SetTextSize(0.037);
//           legendRAAcomboPi0andEta2050->SetHeader(Form("%s",collisionSystemPbPb2050.Data()));
//           legendRAAcomboPi0andEta2050->SetHeader(Form("%s",collisionSystem2760GeV.Data()));
          legendRAAcomboPi0andEta2050->SetHeader("20#font[122]{-}50%");
//           legendRAAcomboPi0andEta2050->AddEntry((TObject*)0,"20#font[122]{-}50%","");
//           legendRAAcomboPi0andEta2050->AddEntry((TObject*)0,"20#font[122]{-}40%","");
          legendRAAcomboPi0andEta2050->AddEntry(graphRAAPi0SysBothMeson_2050,"#pi^{0}","pf");
//           legendRAAcomboPi0andEta2050->AddEntry(graphChargedKaonRAA2040,"K^{#pm}","fp");
          legendRAAcomboPi0andEta2050->AddEntry(graphRAAEtaSysBothMeson_2050,"#eta","fp");
//           legendRAAcomboPi0andEta2050->AddEntry((TObject*)0,"PLB 736 (2014)","");
          legendRAAcomboPi0andEta2050->Draw();
          TLegend* legendRAAKaon2040 = new TLegend(0.7,0.75,0.9,0.88);
          legendRAAKaon2040->SetFillColor(0);
          legendRAAKaon2040->SetLineColor(0);
//           legendRAAKaon2040->SetMargin(0.17);
          legendRAAKaon2040->SetTextFont(42);
          legendRAAKaon2040->SetTextSize(0.037);
          legendRAAKaon2040->SetHeader("20#font[122]{-}40%");
          legendRAAKaon2040->AddEntry(graphChargedKaonRAA2040,"K^{#pm}","fp");
//           legendRAAKaon2040->AddEntry((TObject*)0,"","");
          legendRAAKaon2040->AddEntry((TObject*)0,"PLB 736 (2014)","");
          legendRAAKaon2040->Draw();

          TLatex *labelEnergyRAAcomboPi0andEta = new TLatex(0.57,0.9,collisionSystem2760GeV.Data());
          SetStyleTLatex( labelEnergyRAAcomboPi0andEta, 0.035,4);
          labelEnergyRAAcomboPi0andEta->Draw();
          DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);


        canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEta_2050.%s",outputDir.Data(),suffix.Data()));
//         canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEta_2050.%s",paperPlots.Data(),suffix.Data()));

        if(graphPi0RCSCSysEMCal && graphPi0RCSCEMCal && graphEtaRCSCSysEMCal && graphEtaRCSCEMCal){

            canvasRAAcomboPi0andEta->cd();
              TH2F * histo2DRCPcomboPi0andEta = new TH2F("histo2DRCPcomboPi0andEta","histo2DRCPcomboPi0andEta",11000,0.,70.,1000,-0.5,2.);
              SetStyleHistoTH2ForGraphs(histo2DRCPcomboPi0andEta, "#it{p}_{T} (GeV/#it{c})","#it{R} #frac{Cent}{SemiC}  ",0.035,0.035, 0.035,0.035, .95,0.9);
              histo2DRCPcomboPi0andEta->GetYaxis()->SetRangeUser(0.03,1.2);
              histo2DRCPcomboPi0andEta->GetXaxis()->SetRangeUser(2.,21);
              histo2DRCPcomboPi0andEta->DrawCopy();

                DrawGammaSetMarkerTGraphAsym(graphEtaRCSCSysEMCal, markerStyle0010, markerSizeComb, kOrange+2, kOrange+2, widthLinesBoxes, kTRUE);
                graphEtaRCSCSysEMCal->Draw("E2same");
                DrawGammaSetMarkerTGraphAsym(graphEtaRCSCEMCal, markerStyle0010, markerSizeComb, kOrange+2, kOrange+2);
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphEtaRCSCEMCal);
                graphEtaRCSCEMCal->Draw("p,same");

                DrawGammaSetMarkerTGraphAsym(graphPi0RCSCSysEMCal,markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
                graphPi0RCSCSysEMCal->Draw("E2same");
                DrawGammaSetMarkerTGraphAsym(graphPi0RCSCEMCal, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphPi0RCSCEMCal);
                graphPi0RCSCEMCal->Draw("p,same");

        //         ratioRCSC->Draw("same");
    //             labelSystRaa->Draw();

                TLegend* legendRCSC = new TLegend(0.6,0.8,0.8,0.92);
                legendRCSC->SetFillColor(0);
                legendRCSC->SetLineColor(0);
                legendRCSC->SetTextFont(42);
//                 legendRCSC->SetMargin(0.17);
//                 legendRCSC->SetNColumns(2);
                legendRCSC->SetTextSize(FontSize);
                legendRCSC->SetHeader(collisionSystem2760GeV.Data());
                legendRCSC->AddEntry(graphPi0RCSCSysEMCal,"#pi^{0}","pf");
                legendRCSC->AddEntry(graphEtaRCSCSysEMCal,"#eta","pf");
                legendRCSC->Draw();
//               DrawGammaLines(2., 20.5 , 1., 1.,0.5,   kGray);

            canvasRAAcomboPi0andEta->Update();
            canvasRAAcomboPi0andEta->SaveAs(Form("%s/RCSC_EMCal.%s",outputDir.Data(),suffix.Data()));
        }

      canvasRAAcomboPi0andEta->cd();
//       canvasRAAcomboPi0andEta->SetLogx();
      histo2DRAAcomboPi0andEta->Draw("copy");
      histo2DRAAcomboPi0andEta->GetYaxis()->SetRangeUser(0.,1.2);
      histo2DRAAcomboPi0andEta->GetXaxis()->SetRangeUser(0.,21);

          TBox* boxErrorNorm0010OnlyV2      = CreateBoxConv(kRed-7, 0.25, 1.-normErr0010 , 0.5, 1.+normErr0010);
          boxErrorNorm0010OnlyV2->Draw();
          DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0010, 33,2.3, colorCharged,colorCharged, 0.1, kFALSE);
          DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0010,33,2.3, colorCharged,colorCharged, 2, kTRUE);

          DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA0010, 27,2.3, colorCharged,colorCharged, 0.1, kFALSE);
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedPionRAA0010);
          DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys0010,27,2.3, colorCharged,colorCharged, 1, kTRUE);
//           graphChargedKaonRAASys0010->SetFillStyle(3005);
//           graphChargedKaonRAASys0010->SetFillColor(colorCharged+1);
          graphChargedPionRAASys0010->SetFillStyle(3001);
          graphChargedPionRAASys0010->SetFillColor(colorCharged);

          graphChargedKaonRAASys0010->Draw("E2same");
          graphChargedKaonRAA0010->Draw("p,same");
          graphChargedPionRAASys0010->Draw("E2same");
          graphChargedPionRAA0010->Draw("p,same");
          graphRAAPi0StatBothMeson_0010->Draw("p,same");
          graphRAAPi0SysBothMeson_0010->Draw("E2same");
          graphRAAPi0SysBothMeson_0010->SetFillStyle(3005);
          graphRAAPi0SysBothMeson_0010->SetFillColor(colorPCM0010);
          graphRAAEtaStatBothMeson_0010->Draw("p,same");
          graphRAAEtaSysBothMeson_0010->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(graphRAAEtaSysBothMeson_0010,20,2.3, colorPCM0010,colorPCM0010, 2, kTRUE);

//           graphChargedKaonRAATotErr0010->Draw("p,same");
//             if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAATotErr0010);
//           graphChargedPionRAATotErr0010->Draw("p,same");
//             if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedPionRAATotErr0010);
//             TGraphAsymmErrors* grapPi0RAATotErr0010 = AddErrorsOfGraphsQuadratically(graphRAAPi0StatBothMeson_0010,graphRAAPi0SysBothMeson_0010);
//             DrawGammaSetMarkerTGraphAsym(grapPi0RAATotErr0010, markerStyle0010+4, markerSizeComb, colorPCM0010, colorPCM0010);
//             TGraphAsymmErrors* graphEtaRAATotErr0010 = AddErrorsOfGraphsQuadratically(graphRAAEtaStatBothMeson_0010,graphRAAEtaSysBothMeson_0010);
//             DrawGammaSetMarkerTGraphAsym(graphEtaRAATotErr0010, markerStyle0010, markerSizeComb, colorPCM0010 , colorPCM0010);
//           grapPi0RAATotErr0010->Draw("p,same");
//           graphEtaRAATotErr0010->Draw("p,same");

          TLegend* legendRAAKaonPion00101 = new TLegend(0.55,0.87-(0.037*2.7),0.81,0.87);//0.52,0.87-(0.037*1.3),0.8,0.87);
          legendRAAKaonPion00101->SetFillColor(0);
          legendRAAKaonPion00101->SetLineColor(0);
          legendRAAKaonPion00101->SetTextFont(42);
          legendRAAKaonPion00101->SetMargin(0.35);
          legendRAAKaonPion00101->SetTextSize(0.037);
          legendRAAKaonPion00101->SetNColumns(2);
          legendRAAKaonPion00101->SetHeader(cent0010.Data());
//           legendRAAKaonPion00101->AddEntry(grapPi0RAATotErr0010,"#pi^{0}","pl");
//           legendRAAKaonPion00101->AddEntry(graphEtaRAATotErr0010,"#eta","lp");
          legendRAAKaonPion00101->AddEntry(graphRAAPi0SysBothMeson_0010,"#pi^{0}","pf");
          legendRAAKaonPion00101->AddEntry(graphRAAEtaSysBothMeson_0010,"#eta","fp");
          legendRAAKaonPion00101->Draw();
          TLegend* legendRAAKaonPion0010 = new TLegend(0.55,0.77-(0.037*2.7),0.86,0.77);//0.515,0.812-(2.5*0.037),0.86,0.812);
          legendRAAKaonPion0010->SetFillColor(0);
          legendRAAKaonPion0010->SetLineColor(0);
          legendRAAKaonPion0010->SetTextFont(42);
          legendRAAKaonPion0010->SetMargin(0.3);
          legendRAAKaonPion0010->SetTextSize(0.037);
          legendRAAKaonPion0010->SetNColumns(2);
          legendRAAKaonPion0010->SetHeader("PLB 736 (2014)");
//           legendRAAKaonPion0010->AddEntry((TObject*)0,"","");
//           legendRAAKaonPion0010->AddEntry(graphChargedPionRAATotErr0010,"#pi^{#pm}","lp");
//           legendRAAKaonPion0010->AddEntry(graphChargedKaonRAATotErr0010,"K^{#pm}","lp");
          legendRAAKaonPion0010->AddEntry(graphChargedPionRAASys0010,"#pi^{#pm}","fp");
          legendRAAKaonPion0010->AddEntry(graphChargedKaonRAASys0010,"K^{#pm}","fp");
          legendRAAKaonPion0010->Draw();
          TLatex *labelEnergyRAAcomboPi0andEta0010andCharged = new TLatex(0.55,0.88,collisionSystem2760GeV.Data());
          SetStyleTLatex( labelEnergyRAAcomboPi0andEta0010andCharged, 0.037,4);
          labelEnergyRAAcomboPi0andEta0010andCharged->Draw();


//           TLegend* legendRAAKaonPion00101 = new TLegend(0.5,0.76,0.83,0.81);
//           legendRAAKaonPion00101->SetFillColor(0);
//           legendRAAKaonPion00101->SetLineColor(0);
//           legendRAAKaonPion00101->SetTextFont(42);
//           legendRAAKaonPion00101->SetTextSize(0.037);
//           legendRAAKaonPion00101->SetNColumns(2);
// //           legendRAAKaonPion00101->SetHeader("");
//           legendRAAKaonPion00101->AddEntry(graphRAAPi0SysBothMeson_0010,"#pi^{0}","pf");
//           legendRAAKaonPion00101->AddEntry(graphRAAEtaSysBothMeson_0010,"#eta","fp");
//           legendRAAKaonPion00101->Draw();
//           TLegend* legendRAAKaonPion0010 = new TLegend(0.5,0.65,0.9,0.76);
//           legendRAAKaonPion0010->SetFillColor(0);
//           legendRAAKaonPion0010->SetLineColor(0);
//           legendRAAKaonPion0010->SetTextFont(42);
//           legendRAAKaonPion0010->SetTextSize(0.037);
//           legendRAAKaonPion0010->SetNColumns(2);
//           legendRAAKaonPion0010->SetHeader("PLB 736 (2014)");
// //           legendRAAKaonPion0010->AddEntry((TObject*)0,"","");
//           legendRAAKaonPion0010->AddEntry(graphChargedPionRAASys0010,"#pi^{#pm}","fp");
//           legendRAAKaonPion0010->AddEntry(graphChargedKaonRAASys0010,"K^{#pm}","fp");
//           legendRAAKaonPion0010->Draw();
//           TLatex *labelEnergyRAAcomboPi0andEta0010andCharged = new TLatex(0.5,0.835,collisionSystemPbPb0010.Data());
//           SetStyleTLatex( labelEnergyRAAcomboPi0andEta0010andCharged, 0.037,4);
//           labelEnergyRAAcomboPi0andEta0010andCharged->Draw();
          DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);

        if(thesisPlotting){
          TLatex *thesisLabel = new TLatex(0.5,0.88,thisthesis.Data());
          SetStyleTLatex( thesisLabel,FontSize,4);
//           thesisLabel->Draw();
        }

          histo2DRAAcomboPi0andEta->Draw("same,axis");
      canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEtaandCharged_0010.%s",outputDir.Data(),suffix.Data()));
      canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEtaandCharged_0010.%s",paperPlots.Data(),suffix.Data()));

      canvasRAAcomboPi0andEta->cd();
      histo2DRAAcomboPi0andEta->Draw("copy");
          TBox* boxErrorNorm2050OnlyV2      = CreateBoxConv(colorEMCal2050, 0.25, 1.-normErr2050 , 0.5, 1.+normErr2050);
          boxErrorNorm2050OnlyV2->Draw();
          DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA2040, 33,2.5, colorCharged,colorCharged, 0.1, kFALSE);
          DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys2040,33,2.5, colorCharged,colorCharged, 2, kTRUE);

          DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA2040, 27,2.5, colorCharged,colorCharged, 0.1, kFALSE);
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedPionRAA2040);
          DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys2040,27,2.5, colorCharged,colorCharged, 1, kTRUE);

//           graphChargedKaonRAASys2040->SetFillStyle(3005);
//           graphChargedKaonRAASys2040->SetFillColor(colorCharged+1);
          graphChargedPionRAASys2040->SetFillStyle(3001);
          graphChargedPionRAASys2040->SetFillColor(colorCharged);

          graphChargedKaonRAASys2040->Draw("E2same");
          graphChargedKaonRAA2040->Draw("p,same");
          graphChargedPionRAASys2040->Draw("E2same");
          graphChargedPionRAA2040->Draw("p,same");
          graphRAAPi0SysBothMeson_2050->Draw("E2same");
          graphRAAPi0SysBothMeson_2050->SetFillStyle(3005);
          graphRAAPi0SysBothMeson_2050->SetFillColor(colorEMCal2050);
          graphRAAPi0StatBothMeson_2050->Draw("p,same");
          graphRAAEtaSysBothMeson_2050->Draw("E2same");
          graphRAAEtaStatBothMeson_2050->Draw("p,same");
          DrawGammaSetMarkerTGraphAsym(graphRAAEtaSysBothMeson_2050,21,2.3, colorEMCal2050,colorEMCal2050, 2, kTRUE);
//           graphChargedKaonRAATotErr2040->Draw("p,same");
//             if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAATotErr2040);
//             graphChargedPionRAATotErr2040->Draw("p,same");
//             if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedPionRAATotErr2040);
//             TGraphAsymmErrors* graphPi0RAATotErr2050 = AddErrorsOfGraphsQuadratically(graphRAAPi0StatBothMeson_2050,graphRAAPi0SysBothMeson_2050);
//             DrawGammaSetMarkerTGraphAsym(graphPi0RAATotErr2050, markerStyle2050+4, markerSizeComb, colorEMCal2050 , colorEMCal2050);
//             TGraphAsymmErrors* graphEtaRAATotErr2050 = AddErrorsOfGraphsQuadratically(graphRAAEtaStatBothMeson_2050,graphRAAEtaSysBothMeson_2050);
//             DrawGammaSetMarkerTGraphAsym(graphEtaRAATotErr2050, markerStyle2050, markerSizeComb, colorEMCal2050 , colorEMCal2050);
//           graphPi0RAATotErr2050->Draw("p,same");
//           graphEtaRAATotErr2050->Draw("p,same");

          TLegend* legendRAAKaonPion20501 = new TLegend(0.55,0.87-(0.037*2.7),0.81,0.87);
          legendRAAKaonPion20501->SetFillColor(0);
          legendRAAKaonPion20501->SetLineColor(0);
          legendRAAKaonPion20501->SetTextFont(42);
          legendRAAKaonPion20501->SetTextSize(0.037);
          legendRAAKaonPion20501->SetNColumns(2);
          legendRAAKaonPion20501->SetMargin(0.35);
          legendRAAKaonPion20501->SetHeader("20#font[122]{-}50%");
//           legendRAAKaonPion20501->AddEntry(graphPi0RAATotErr2050,"#pi^{0}","pl");
//           legendRAAKaonPion20501->AddEntry(graphEtaRAATotErr2050,"#eta","lp");
          legendRAAKaonPion20501->AddEntry(graphRAAPi0SysBothMeson_2050,"#pi^{0}","pf");
          legendRAAKaonPion20501->AddEntry(graphRAAEtaSysBothMeson_2050,"#eta","fp");
          legendRAAKaonPion20501->Draw();
          TLegend* legendRAAKaonPion2050 = new TLegend(0.55,0.77-(0.037*2.7),0.86,0.77);
          legendRAAKaonPion2050->SetFillColor(0);
          legendRAAKaonPion2050->SetLineColor(0);
          legendRAAKaonPion2050->SetTextFont(42);
          legendRAAKaonPion2050->SetTextSize(0.037);
          legendRAAKaonPion2050->SetNColumns(2);
          legendRAAKaonPion2050->SetMargin(0.3);
          legendRAAKaonPion2050->SetHeader("20#font[122]{-}40%, PLB 736 (2014)");
          legendRAAKaonPion2050->AddEntry(graphChargedPionRAASys2040,"#pi^{#pm}","fp");
          legendRAAKaonPion2050->AddEntry(graphChargedKaonRAASys2040,"K^{#pm}","fp");
//           legendRAAKaonPion2050->AddEntry(graphChargedPionRAATotErr2040,"#pi^{#pm}","lp");
//           legendRAAKaonPion2050->AddEntry(graphChargedKaonRAATotErr2040,"K^{#pm}","lp");
          legendRAAKaonPion2050->Draw();
          TLatex *labelEnergyRAAcomboPi0andEtaanfCharged = new TLatex(0.55,0.88,collisionSystem2760GeV.Data());
          SetStyleTLatex( labelEnergyRAAcomboPi0andEtaanfCharged, 0.037,4);
          labelEnergyRAAcomboPi0andEtaanfCharged->Draw();
          DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);

          histo2DRAAcomboPi0andEta->Draw("same,axis");
      canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEtaandCharged_2050.%s",outputDir.Data(),suffix.Data()));
      canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEtaandCharged_2050.%s",paperPlots.Data(),suffix.Data()));


        Double_t arrayBoundariesX_RAA4pads[3];
        Double_t arrayBoundariesY_RAA4pads[3];
        Double_t relativeMarginsX_RAA4pads[3];
        Double_t relativeMarginsY_RAA4pads[3];
        Double_t marginRatio = 0.16*2400;
        Double_t textsizeLabelsRatioUp = 0;
        Double_t textsizeFacRatioUp = 0;
        Double_t textsizeLabelsRatioDown = 0;
        Double_t textsizeFacRatioDown = 0;

        TCanvas* canvasRAAwithModels_4pads = new TCanvas("canvasRAAwithModels_4pads","",0,0,2400,2000);  // gives the page size
        ReturnCorrectValuesForCanvasScaling(2400,2000, 2, 2,0.05, 0.008, 0.008,0.06,arrayBoundariesX_RAA4pads,arrayBoundariesY_RAA4pads,relativeMarginsX_RAA4pads,relativeMarginsY_RAA4pads);

        TPad* padRAAwithModels_pad1 = new TPad("", "", arrayBoundariesX_RAA4pads[0], arrayBoundariesY_RAA4pads[2], arrayBoundariesX_RAA4pads[1], arrayBoundariesY_RAA4pads[1],-1, -1, -2);
        DrawGammaPadSettings( padRAAwithModels_pad1, relativeMarginsX_RAA4pads[0], relativeMarginsX_RAA4pads[1], relativeMarginsY_RAA4pads[1], relativeMarginsY_RAA4pads[2]);
            padRAAwithModels_pad1->Draw();
        TPad* padRAAwithModels_pad2 = new TPad("", "", arrayBoundariesX_RAA4pads[0], arrayBoundariesY_RAA4pads[1], arrayBoundariesX_RAA4pads[1], arrayBoundariesY_RAA4pads[0],-1, -1, -2);
        DrawGammaPadSettings( padRAAwithModels_pad2, relativeMarginsX_RAA4pads[0], relativeMarginsX_RAA4pads[1], relativeMarginsY_RAA4pads[0], relativeMarginsY_RAA4pads[1]);
            padRAAwithModels_pad2->Draw();
        TPad* padRAAwithModels_pad3 = new TPad("", "", arrayBoundariesX_RAA4pads[1], arrayBoundariesY_RAA4pads[2], arrayBoundariesX_RAA4pads[2], arrayBoundariesY_RAA4pads[1],-1, -1, -2);
        DrawGammaPadSettings( padRAAwithModels_pad3, relativeMarginsX_RAA4pads[1], relativeMarginsX_RAA4pads[2], relativeMarginsY_RAA4pads[1], relativeMarginsY_RAA4pads[2]);
            padRAAwithModels_pad3->Draw();
        TPad* padRAAwithModels_pad4 = new TPad("", "", arrayBoundariesX_RAA4pads[1], arrayBoundariesY_RAA4pads[1], arrayBoundariesX_RAA4pads[2], arrayBoundariesY_RAA4pads[0],-1, -1, -2);
        DrawGammaPadSettings( padRAAwithModels_pad4, relativeMarginsX_RAA4pads[1], relativeMarginsX_RAA4pads[2], relativeMarginsY_RAA4pads[0], relativeMarginsY_RAA4pads[1]);
            padRAAwithModels_pad4->Draw();
        if (padRAAwithModels_pad2->XtoPixel(padRAAwithModels_pad2->GetX2()) < padRAAwithModels_pad2->YtoPixel(padRAAwithModels_pad2->GetY1())){
            textsizeLabelsRatioUp = (Double_t)50/padRAAwithModels_pad2->XtoPixel(padRAAwithModels_pad2->GetX2()) ;
            textsizeFacRatioUp = (Double_t)1./padRAAwithModels_pad2->XtoPixel(padRAAwithModels_pad2->GetX2()) ;
        } else {
            textsizeLabelsRatioUp = (Double_t)50/padRAAwithModels_pad2->YtoPixel(padRAAwithModels_pad2->GetY1());
            textsizeFacRatioUp = (Double_t)1./padRAAwithModels_pad2->YtoPixel(padRAAwithModels_pad2->GetY1());
        }
        if (padRAAwithModels_pad1->XtoPixel(padRAAwithModels_pad1->GetX2()) < padRAAwithModels_pad1->YtoPixel(padRAAwithModels_pad1->GetY1())){
            textsizeLabelsRatioDown = (Double_t)50/padRAAwithModels_pad1->XtoPixel(padRAAwithModels_pad1->GetX2()) ;
            textsizeFacRatioDown = (Double_t)1./padRAAwithModels_pad1->XtoPixel(padRAAwithModels_pad1->GetX2()) ;
        } else {
            textsizeLabelsRatioDown = (Double_t)50/padRAAwithModels_pad1->YtoPixel(padRAAwithModels_pad1->GetY1());
            textsizeFacRatioDown = (Double_t)1./padRAAwithModels_pad1->YtoPixel(padRAAwithModels_pad1->GetY1());
        }

        TH2F * histo2DRAATopRowPad = new TH2F("histo2DRAATopRowPad","histo2DRAATopRowPad",1000,0.0001,30.,1000,-0.05,10.);
            SetStyleHistoTH2ForGraphs(histo2DRAATopRowPad, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}  ", 0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp, 0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp, 1,0.3/(textsizeFacRatioUp*marginRatio), 512, 505);
            histo2DRAATopRowPad->GetXaxis()->SetRangeUser(0.1,21);
            histo2DRAATopRowPad->GetYaxis()->SetRangeUser(0.01,1.4);
        TH2F * histo2DRAABottomRowPad = new TH2F("histo2DRAABottomRowPad","histo2DRAABottomRowPad",1000,0.0001,30.,1000,-0.05,10.);
            SetStyleHistoTH2ForGraphs(histo2DRAABottomRowPad, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}  ", 0.85*textsizeLabelsRatioDown, textsizeLabelsRatioDown,  0.85*textsizeLabelsRatioDown, textsizeLabelsRatioDown, 1,0.3/(textsizeFacRatioDown*marginRatio), 512, 505);
            histo2DRAABottomRowPad->GetXaxis()->SetRangeUser(0.1,21);
            histo2DRAABottomRowPad->GetYaxis()->SetRangeUser(0.01,1.4);

        padRAAwithModels_pad1->cd();
//         padRAAwithModels_pad1->SetLogx();
        histo2DRAABottomRowPad->DrawCopy();

            boxErrorNorm0010_Single->Draw();
            DrawGammaLines(0.1, 21, 1, 1 ,1,kGray, 2);

            DrawGammaSetMarkerTGraphAsym(graphEtaRAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
            graphEtaRAAJetQuenching_0010->Draw("3,same");

            DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
            gWHDG_Eta_Raa_0010->SetFillStyle(fillStyleWHDG);
            gWHDG_Eta_Raa_0010->SetFillColor(colorWHDG0005);
            gWHDG_Eta_Raa_0010->Draw("3 same");

            DrawGammaSetMarkerTGraphAsym(graphRAAEtaSysBothMeson_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010, widthLinesBoxes, kTRUE);
            graphRAAEtaSysBothMeson_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphRAAEtaStatBothMeson_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);
            graphRAAEtaStatBothMeson_0010->Draw("p,same");

            TLegend* legendRAAwithModels_eta0010 = new TLegend(0.19,0.86,0.55,0.96);
            legendRAAwithModels_eta0010->SetFillColor(0);
            legendRAAwithModels_eta0010->SetLineColor(0);
            legendRAAwithModels_eta0010->SetTextFont(42);
            legendRAAwithModels_eta0010->SetMargin(0.17);
            legendRAAwithModels_eta0010->SetTextSize(0.85*textsizeLabelsRatioDown);
            legendRAAwithModels_eta0010->SetHeader(collisionSystem2760GeV.Data());
            legendRAAwithModels_eta0010->AddEntry(graphRAAEtaSysBothMeson_0010,"#eta,  0#font[122]{-}10%","fp");
            legendRAAwithModels_eta0010->Draw();

            if(thesisPlotting){
                TLatex *thesisLabel = new TLatex(0.62,0.7,thisthesis.Data());
                SetStyleTLatex( thesisLabel,0.85*textsizeLabelsRatioDown,4);
//                 thesisLabel->Draw();
            }

        histo2DRAABottomRowPad->Draw("axis,same");
        padRAAwithModels_pad1->Update();
        padRAAwithModels_pad2->cd();
//         padRAAwithModels_pad2->SetLogx();
        histo2DRAATopRowPad->DrawCopy();

            boxErrorNorm0010_Single->Draw();
            DrawGammaLines(0.1, 21, 1, 1 ,1,kGray, 2);

            DrawGammaSetMarkerTGraphAsym(graphPi0RAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
            graphPi0RAAJetQuenching_0010->Draw("3,same");

            DrawGammaSetMarkerTGraphAsym(gWHDG_Pi0_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
            gWHDG_Pi0_Raa_0010->SetFillStyle(fillStyleWHDG);
            gWHDG_Pi0_Raa_0010->SetFillColor(colorWHDG0005);
            gWHDG_Pi0_Raa_0010->Draw("3 same");

            DrawGammaSetMarkerTGraphAsym(graphPi0Djordjevic_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorVitevBas0005,colorVitevBas0005,widthLinesBoxes, kTRUE);
            graphPi0Djordjevic_0010->SetFillStyle(fillStyleVitev);
            graphPi0Djordjevic_0010->SetFillColor(colorVitevBas0005);
            graphPi0Djordjevic_0010->Draw("3 same");

            TGraphAsymmErrors *copyWHDG = (TGraphAsymmErrors*) gWHDG_Pi0_Raa_0010->Clone();
            DrawGammaSetMarkerTGraphAsym(copyWHDG, markerStyleCommmonSpectrum0010,markerSizeSpectrum, kBlack, kBlack,widthLinesBoxes, kTRUE);
            copyWHDG->SetFillStyle(fillStyleWHDG);
            copyWHDG->SetFillColor(kBlack);

            TGraphAsymmErrors *copyDj = (TGraphAsymmErrors*) graphPi0Djordjevic_0010->Clone();
            DrawGammaSetMarkerTGraphAsym(copyDj, markerStyleCommmonSpectrum0010,markerSizeSpectrum, kBlack, kBlack,widthLinesBoxes, kTRUE);
            copyDj->SetFillStyle(fillStyleVitev);
            copyDj->SetFillColor(kBlack);

            DrawGammaSetMarkerTGraphAsym(graphRAAPi0SysBothMeson_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010, widthLinesBoxes, kTRUE);
            graphRAAPi0SysBothMeson_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphRAAPi0StatBothMeson_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010);
            graphRAAPi0StatBothMeson_0010->Draw("p,same");

            TLegend* legendRAAwithModels_pi00010 = new TLegend(0.19,0.72,0.89,0.95);
            legendRAAwithModels_pi00010->SetFillColor(0);
            legendRAAwithModels_pi00010->SetLineColor(0);
            legendRAAwithModels_pi00010->SetTextFont(42);
            legendRAAwithModels_pi00010->SetMargin(0.17);
            legendRAAwithModels_pi00010->SetNColumns(2);
            legendRAAwithModels_pi00010->SetTextSize(0.85*textsizeLabelsRatioUp);
            legendRAAwithModels_pi00010->SetHeader(collisionSystem2760GeV.Data());
            legendRAAwithModels_pi00010->AddEntry(graphRAAPi0SysBothMeson_0010,"#pi^{0},  0#font[122]{-}10%","fp");
            legendRAAwithModels_pi00010->AddEntry(graphPi0RAAJetQuenching_0010,"0#font[122]{-}10% NLO DCZW","f");
            legendRAAwithModels_pi00010->AddEntry((TObject*)0,"","");
            legendRAAwithModels_pi00010->AddEntry(copyWHDG,"WHDG","f");
            legendRAAwithModels_pi00010->AddEntry((TObject*)0,"","");
            legendRAAwithModels_pi00010->AddEntry(copyDj,"Djordjevic","f");
            legendRAAwithModels_pi00010->Draw();

            if(thesisPlotting){
                TLatex *thesisLabel = new TLatex(0.62,0.65,thisthesis.Data());
                SetStyleTLatex( thesisLabel,0.85*textsizeLabelsRatioUp,4);
//                 thesisLabel->Draw();
            }

        histo2DRAATopRowPad->Draw("axis,same");
        padRAAwithModels_pad2->Update();
        padRAAwithModels_pad3->cd();
//         padRAAwithModels_pad3->SetLogx();
        histo2DRAABottomRowPad->DrawCopy();

            boxErrorNorm2050_Single->Draw();
            DrawGammaLines(0.1, 21, 1, 1 ,1,kGray, 2);

            DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
            gWHDG_Eta_Raa_2050->SetFillStyle(fillStyleWHDG);
            gWHDG_Eta_Raa_2050->SetFillColor(colorWHDG4060);
            gWHDG_Eta_Raa_2050->Draw("3 same");

            DrawGammaSetMarkerTGraphAsym(graphRAAEtaSysBothMeson_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050, widthLinesBoxes, kTRUE);
            graphRAAEtaSysBothMeson_2050->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphRAAEtaStatBothMeson_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050);
            graphRAAEtaStatBothMeson_2050->Draw("p,same");

            TLegend* legendRAAwithModels_eta2050 = new TLegend(0.1,0.86,0.47,0.96);
            legendRAAwithModels_eta2050->SetFillColor(0);
            legendRAAwithModels_eta2050->SetLineColor(0);
            legendRAAwithModels_eta2050->SetTextFont(42);
            legendRAAwithModels_eta2050->SetMargin(0.17);
            legendRAAwithModels_eta2050->SetTextSize(0.85*textsizeLabelsRatioDown);
            legendRAAwithModels_eta2050->SetHeader(collisionSystem2760GeV.Data());
            legendRAAwithModels_eta2050->AddEntry(graphRAAEtaSysBothMeson_2050,"#eta, 20#font[122]{-}50%","fp");
            legendRAAwithModels_eta2050->Draw();

            if(thesisPlotting){
                TLatex *thesisLabel = new TLatex(0.6,0.7,thisthesis.Data());
                SetStyleTLatex( thesisLabel,0.85*textsizeLabelsRatioDown,4);
//                 thesisLabel->Draw();
            }

        histo2DRAABottomRowPad->Draw("axis,same");
        padRAAwithModels_pad3->Update();
        padRAAwithModels_pad4->cd();
//         padRAAwithModels_pad4->SetLogx();
        histo2DRAATopRowPad->DrawCopy();

            boxErrorNorm2050_Single->Draw();
            DrawGammaLines(0.1, 21, 1, 1 ,1, kGray, 2);

            DrawGammaSetMarkerTGraphAsym(gWHDG_Pi0_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
            gWHDG_Pi0_Raa_2050->SetFillStyle(fillStyleWHDG);
            gWHDG_Pi0_Raa_2050->SetFillColor(colorWHDG4060);
            gWHDG_Pi0_Raa_2050->Draw("3 same");

            DrawGammaSetMarkerTGraphAsym(graphPi0Djordjevic_2050, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorVitevBas4060,colorVitevBas4060,widthLinesBoxes, kTRUE);
            graphPi0Djordjevic_2050->SetFillStyle(fillStyleVitev);
            graphPi0Djordjevic_2050->SetFillColor(colorVitevBas4060);
            graphPi0Djordjevic_2050->Draw("3 same");

            DrawGammaSetMarkerTGraphAsym(graphRAAPi0SysBothMeson_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050, widthLinesBoxes, kTRUE);
            graphRAAPi0SysBothMeson_2050->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphRAAPi0StatBothMeson_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050);
            graphRAAPi0StatBothMeson_2050->Draw("p,same");

            TLegend* legendRAAwithModels_pi02050 = new TLegend(0.1,0.83,0.47,0.94);
            legendRAAwithModels_pi02050->SetFillColor(0);
            legendRAAwithModels_pi02050->SetLineColor(0);
            legendRAAwithModels_pi02050->SetTextFont(42);
            legendRAAwithModels_pi02050->SetMargin(0.17);
            legendRAAwithModels_pi02050->SetTextSize(0.85*textsizeLabelsRatioUp);
            legendRAAwithModels_pi02050->SetHeader(collisionSystem2760GeV.Data());
            legendRAAwithModels_pi02050->AddEntry(graphRAAPi0SysBothMeson_2050,"#pi^{0}, 20#font[122]{-}50%","fp");
            legendRAAwithModels_pi02050->Draw();

            if(thesisPlotting){
                TLatex *thesisLabel = new TLatex(0.6,0.65,thisthesis.Data());
                SetStyleTLatex( thesisLabel,0.85*textsizeLabelsRatioUp,4);
//                 thesisLabel->Draw();
            }

        histo2DRAATopRowPad->Draw("axis,same");
        padRAAwithModels_pad4->Update();
        canvasRAAwithModels_4pads->Update();
        canvasRAAwithModels_4pads->SaveAs(Form("%s/RAA_TheoryModels.%s",outputDir.Data(),suffix.Data()));
        canvasRAAwithModels_4pads->SaveAs(Form("%s/RAA_TheoryModels.%s",paperPlots.Data(),suffix.Data()));
        delete padRAAwithModels_pad1;
        delete padRAAwithModels_pad3;
        delete padRAAwithModels_pad2;
        delete padRAAwithModels_pad4;
        delete canvasRAAwithModels_4pads;
      }
    }


    canvasRAAcombo->cd();
    histo2DRAAcombo->DrawCopy();

        TLegend* legendRAATheoryPbPb = new TLegend(0.53,0.5,0.95,0.68);
        legendRAATheoryPbPb->SetFillColor(0);
        legendRAATheoryPbPb->SetLineColor(0);
        legendRAATheoryPbPb->SetTextFont(42);
        legendRAATheoryPbPb->SetTextSize(FontSize);
//         legendRAATheoryPbPb->SetMargin(0.2);
        legendRAATheoryPbPb->SetNColumns(2);
        legendRAATheoryPbPb->AddEntry((TObject*)0,cent0010.Data(),"");
        legendRAATheoryPbPb->AddEntry((TObject*)0,cent2050.Data(),"");

        if(meson.CompareTo("Pi0")==0){

            DrawGammaSetMarkerTGraphAsym(graphPi0RAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
            graphPi0RAAJetQuenching_0010->Draw("3,same");

            DrawGammaSetMarkerTGraphAsym(gWHDG_Pi0_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
            gWHDG_Pi0_Raa_0010->SetFillStyle(fillStyleWHDG);
            gWHDG_Pi0_Raa_0010->SetFillColor(colorWHDG0005);
            gWHDG_Pi0_Raa_0010->Draw("3 same");

            DrawGammaSetMarkerTGraphAsym(gWHDG_Pi0_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
            gWHDG_Pi0_Raa_2050->SetFillStyle(fillStyleWHDG);
            gWHDG_Pi0_Raa_2050->SetFillColor(colorWHDG4060);
            gWHDG_Pi0_Raa_2050->Draw("3 same");

            DrawGammaSetMarkerTGraphAsym(graphPi0Djordjevic_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorVitevBas0005, colorVitevBas0005,widthLinesBoxes, kTRUE);
            graphPi0Djordjevic_0010->SetFillStyle(fillStyleVitev);
            graphPi0Djordjevic_0010->SetFillColor(colorVitevBas0005);
            graphPi0Djordjevic_0010->Draw("3 same");

            DrawGammaSetMarkerTGraphAsym(graphPi0Djordjevic_2050, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorVitevBas4060,colorVitevBas4060,widthLinesBoxes, kTRUE);
            graphPi0Djordjevic_2050->SetFillStyle(fillStyleVitev);
            graphPi0Djordjevic_2050->SetFillColor(colorVitevBas4060);
            graphPi0Djordjevic_2050->Draw("3 same");

            legendRAATheoryPbPb->AddEntry(graphPi0RAAJetQuenching_0010,"NLO DCZW","f");
            legendRAATheoryPbPb->AddEntry((TObject*)0,"","");
            legendRAATheoryPbPb->AddEntry(gWHDG_Pi0_Raa_0010,"WHDG","f");
            legendRAATheoryPbPb->AddEntry(gWHDG_Pi0_Raa_2050,"WHDG","f");
            legendRAATheoryPbPb->AddEntry(graphPi0Djordjevic_0010,"Djordjevic","f");
            legendRAATheoryPbPb->AddEntry(graphPi0Djordjevic_2050,"Djordjevic","f");


        } else if(meson.CompareTo("Eta")==0){

            DrawGammaSetMarkerTGraphAsym(graphEtaRAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
            graphEtaRAAJetQuenching_0010->Draw("3,same");

            DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
            gWHDG_Eta_Raa_0010->SetFillStyle(fillStyleWHDG);
            gWHDG_Eta_Raa_0010->SetFillColor(colorWHDG0005);
            gWHDG_Eta_Raa_0010->Draw("3 same");

            DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
            gWHDG_Eta_Raa_2050->SetFillStyle(fillStyleWHDG);
            gWHDG_Eta_Raa_2050->SetFillColor(colorWHDG4060);
            gWHDG_Eta_Raa_2050->Draw("3 same");

            legendRAATheoryPbPb->AddEntry(graphEtaRAAJetQuenching_0010,"NLO DCZW","f");
            legendRAATheoryPbPb->AddEntry((TObject*)0,"","");
            legendRAATheoryPbPb->AddEntry(gWHDG_Eta_Raa_0010,"WHDG","f");
            legendRAATheoryPbPb->AddEntry(gWHDG_Eta_Raa_2050,"WHDG","f");
            legendRAATheoryPbPb->AddEntry((TObject*)0,"","");
            legendRAATheoryPbPb->AddEntry((TObject*)0,"","");


        }
        legendRAATheoryPbPb->Draw();

        labelSystRaa->Draw();
        legendRAAcombo2->Draw();
        boxErrorNorm0010->Draw();
        boxErrorNorm2050->Draw();
        DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);


//         graphDirectCombRAASysPbPb2760GeV_0010->Draw("E2same");
//         if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPbPbPb2760GeV_0010);
//         graphCombRAAStatPbPbPbPb2760GeV_0010->Draw("p,same");
//         graphDirectCombRAASysPbPb2760GeV_2050->Draw("E2same");
//         if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPbPbPb2760GeV_2050);
//         graphCombRAAStatPbPbPbPb2760GeV_2050->Draw("p,same");

        graphCombRAASysPbPb2760GeV_0010->Draw("E2same");
        graphCombRAAStatPbPb2760GeV_0010->Draw("p,same");
        graphCombRAASysPbPb2760GeV_2050->Draw("E2same");
        graphCombRAAStatPbPb2760GeV_2050->Draw("p,same");

            TLatex *thesisLabel = new TLatex(0.12,0.89,thisthesis.Data());
            SetStyleTLatex( thesisLabel,FontSize,4);
        if(thesisPlotting){
            thesisLabel->Draw();
        }


    canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithTheoryModels.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithTheoryModels.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));



    canvasRAAcombo->cd();
    histo2DRAAcombo->DrawCopy();

        labelSystRaa->Draw();

        TLegend* legendRAAcomboCharged = new TLegend(0.6,0.73,0.91,0.93);
        legendRAAcomboCharged->SetFillColor(0);
        legendRAAcomboCharged->SetLineColor(0);
        legendRAAcomboCharged->SetTextFont(42);
        legendRAAcomboCharged->SetMargin(0.17);
        legendRAAcomboCharged->SetTextSize(FontSize);
        legendRAAcomboCharged->SetHeader(collisionSystem2760GeV.Data());

        if(meson.CompareTo("Pi0")==0){

            legendRAAcomboCharged->AddEntry(graphCombRAASysPbPb2760GeV_0010,Form("#pi^{0} - %s ",cent0010.Data()),"fp");
            DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA0010, 24,2, colorCharged,colorCharged, 0.1, kFALSE);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedPionRAA0010);
            DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys0010,24,2, colorCharged,colorCharged, 1, kTRUE);
            graphChargedPionRAASys0010->Draw("2same");
            graphChargedPionRAA0010->Draw("p,same");
            legendRAAcomboCharged->AddEntry(graphChargedPionRAA0010,Form("#pi^{#pm} - %s ",cent0010.Data()),"fp");
            legendRAAcomboCharged->AddEntry((TObject*)0,"PLB 736 (2014)","");

        } else if(meson.CompareTo("Eta")==0){

            legendRAAcomboCharged->AddEntry(graphCombRAASysPbPb2760GeV_0010,Form("#eta - %s ",cent0010.Data()),"fp");
            DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0010, 24,2, colorCharged,colorCharged, 0.1, kFALSE);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAA0010);
            DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0010,24,2, colorCharged,colorCharged, 1, kTRUE);
            graphChargedKaonRAASys0010->Draw("2same");
            graphChargedKaonRAA0010->Draw("p,same");
            legendRAAcomboCharged->AddEntry(graphChargedKaonRAA0010,Form("K^{#pm} - %s ",cent0010.Data()),"fp");
            legendRAAcomboCharged->AddEntry((TObject*)0,"PLB 736 (2014)","");

        }

        graphCombRAASysPbPb2760GeV_0010->Draw("E2same");
        graphCombRAAStatPbPb2760GeV_0010->Draw("p,same");

        boxErrorNorm0010Only->Draw();
                    thesisLabel->Draw();

        DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);
        legendRAAcomboCharged->Draw();

    canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithCharged_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithCharged_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    canvasRAAcombo->cd();
    histo2DRAAcombo->DrawCopy();

        labelSystRaa->Draw();

        TLegend* legendRAAcomboCharged2050 = new TLegend(0.6,0.73,0.91,0.93);
        legendRAAcomboCharged2050->SetFillColor(0);
        legendRAAcomboCharged2050->SetLineColor(0);
        legendRAAcomboCharged2050->SetTextFont(42);
        legendRAAcomboCharged2050->SetMargin(0.17);
        legendRAAcomboCharged2050->SetTextSize(FontSize);
        legendRAAcomboCharged2050->SetHeader(collisionSystem2760GeV.Data());

        if(meson.CompareTo("Pi0")==0){

            legendRAAcomboCharged2050->AddEntry(graphCombRAASysPbPb2760GeV_2050,Form("#pi^{0} - %s ",cent2050.Data()),"fp");
            DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA2040, 25,2, colorCharged,colorCharged, 0.1, kFALSE);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedPionRAA2040);
            DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys2040,25,2, colorCharged,colorCharged, 1, kTRUE);
            graphChargedPionRAASys2040->Draw("2same");
            graphChargedPionRAA2040->Draw("p,same");
            legendRAAcomboCharged2050->AddEntry(graphChargedPionRAA2040,"#pi^{#pm} - 20-40%","fp");
            legendRAAcomboCharged2050->AddEntry((TObject*)0,"PLB 736 (2014)","");

        } else if(meson.CompareTo("Eta")==0){

            legendRAAcomboCharged2050->AddEntry(graphCombRAASysPbPb2760GeV_2050,Form("#eta - %s ",cent2050.Data()),"fp");
            DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA2040, 25,2,  colorCharged,colorCharged, 0.1, kFALSE);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAA2040);
            DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys2040,25,2, colorCharged,colorCharged, 1, kTRUE);
            graphChargedKaonRAASys2040->Draw("2same");
            graphChargedKaonRAA2040->Draw("p,same");
            legendRAAcomboCharged2050->AddEntry(graphChargedKaonRAA2040,"K^{#pm} - 20-40%","fp");
            legendRAAcomboCharged2050->AddEntry((TObject*)0,"PLB 736 (2014)","");

        }

        graphCombRAASysPbPb2760GeV_2050->Draw("E2same");
        graphCombRAAStatPbPb2760GeV_2050->Draw("p,same");

                    thesisLabel->Draw();

        boxErrorNorm2050Only->Draw();
        DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);
        legendRAAcomboCharged2050->Draw();

    canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithCharged_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithCharged_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


    Int_t textSizeLabelsPixelRAA = 50;
    Double_t marginRAA = 0.14*1200;
    Double_t textsizeLabelsRAA = 0;
    Double_t textsizeFacRAA = 0;
    if (canvasRAAcombo->XtoPixel(canvasRAAcombo->GetX2()) < canvasRAAcombo->YtoPixel(canvasRAAcombo->GetY1())){
        textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAAcombo->XtoPixel(canvasRAAcombo->GetX2()) ;
        textsizeFacRAA = (Double_t)1./canvasRAAcombo->XtoPixel(canvasRAAcombo->GetX2()) ;
    } else {
        textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAAcombo->YtoPixel(canvasRAAcombo->GetY1());
        textsizeFacRAA = (Double_t)1./canvasRAAcombo->YtoPixel(canvasRAAcombo->GetY1());
    }

    if(meson.CompareTo("Pi0")==0){
        canvasRAAcombo->cd();
            TH2F * histo2DRAAcomboPHENIX = new TH2F("histo2DRAAcomboPHENIX","histo2DRAAcomboPHENIX",11000,0.,70.,1000,-0.5,2.);
            SetStyleHistoTH2ForGraphs(histo2DRAAcomboPHENIX, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,.92);
            histo2DRAAcomboPHENIX->GetYaxis()->SetRangeUser(0.,2.);
            histo2DRAAcomboPHENIX->GetXaxis()->SetRangeUser(0.,21.);
            histo2DRAAcomboPHENIX->Draw("copy");

            DrawGammaLines(0., 21. , 1, 1 ,1,kGray,2);

            while(graphCombRAASysPbPb2760GeV_0010->GetX()[0]<1)graphCombRAASysPbPb2760GeV_0010->RemovePoint(0);
            while(graphCombRAAStatPbPb2760GeV_0010->GetX()[0]<1)graphCombRAAStatPbPb2760GeV_0010->RemovePoint(0);
            graphCombRAASysPbPb2760GeV_0010->Draw("E2same");
            graphCombRAAStatPbPb2760GeV_0010->SetMarkerSize(2.7);
            graphCombRAAStatPbPb2760GeV_0010->Draw("p,same");



            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVPi0RAA_0010, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
            DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVPi0RAA_0010, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
            DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVPi0RAA_0010, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);
            DrawGammaSetMarkerTGraphErr(graphWA98_17_3GeVPi0RAA_0013, markerStyleWA98,markerSizeWA98, kGray+2 , kGray+2);

            graphPHENIX200GeVPi0RAA_0010->Draw("p,same");
            graphPHENIX39GeVPi0RAA_0010->Draw("p,same");
            graphPHENIX62GeVPi0RAA_0010->Draw("p,same");
            graphWA98_17_3GeVPi0RAA_0013->Draw("p,same");

            TLatex *labelRAAALICEPbPb0010 = new TLatex(0.35,0.9,"#pi^{0} ALICE 0#font[122]{-}10% Pb#font[122]{-}Pb (2011)");
            SetStyleTLatex( labelRAAALICEPbPb0010, 0.85*textsizeLabelsRAA,4);
            labelRAAALICEPbPb0010->Draw();
            TLegend* legendRAASinglePbPb0010 = new TLegend(0.35,0.84,0.65,0.885);
            legendRAASinglePbPb0010->SetFillColor(0);
            legendRAASinglePbPb0010->SetLineColor(0);
            legendRAASinglePbPb0010->SetNColumns(1);
            legendRAASinglePbPb0010->SetTextFont(42);
            legendRAASinglePbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
            legendRAASinglePbPb0010->SetMargin(0.17);
            legendRAASinglePbPb0010->AddEntry(graphCombRAASysPbPb2760GeV_0010,"0#font[122]{-}10% #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
            legendRAASinglePbPb0010->Draw();

            TLatex *labelRAAPHENIXPbPb0010 = new TLatex(0.35,0.79,"#pi^{0} PHENIX 0#font[122]{-}10% Au#font[122]{-}Au");
            SetStyleTLatex( labelRAAPHENIXPbPb0010, 0.85*textsizeLabelsRAA,4);
            labelRAAPHENIXPbPb0010->Draw();

            TLegend* legendRAARHICPbPb0010 = new TLegend(0.35,0.66,0.95,0.78);
            legendRAARHICPbPb0010->SetFillColor(0);
            legendRAARHICPbPb0010->SetLineColor(0);
            legendRAARHICPbPb0010->SetNColumns(2);
            legendRAARHICPbPb0010->SetTextFont(42);
            legendRAARHICPbPb0010->SetMargin(0.17);
            legendRAARHICPbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
            legendRAARHICPbPb0010->AddEntry(graphPHENIX200GeVPi0RAA_0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendRAARHICPbPb0010->AddEntry(graphPHENIX62GeVPi0RAA_0010,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
            legendRAARHICPbPb0010->AddEntry(graphPHENIX39GeVPi0RAA_0010,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
            legendRAARHICPbPb0010->Draw();

            TLatex *labelRAAWA98PbPb0010 = new TLatex(0.35,0.61,"#pi^{0} WA98     0#font[122]{-}13% Pb#font[122]{-}Pb");
            SetStyleTLatex( labelRAAWA98PbPb0010, 0.85*textsizeLabelsRAA,4);
            labelRAAWA98PbPb0010->Draw();

            TLegend* legendRAASPSPbPb0010 = new TLegend(0.35,0.55,0.95,0.59);
            legendRAASPSPbPb0010->SetFillColor(0);
            legendRAASPSPbPb0010->SetLineColor(0);
            legendRAASPSPbPb0010->SetNColumns(2);
            legendRAASPSPbPb0010->SetTextFont(42);
            legendRAASPSPbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
            legendRAASPSPbPb0010->SetMargin(0.17);
            legendRAASPSPbPb0010->AddEntry(graphWA98_17_3GeVPi0RAA_0013,"#sqrt{#it{s}_{_{NN}}} = 17.3 GeV","p");
            legendRAASPSPbPb0010->Draw();

            boxErrorNorm0010_Single->Draw();
            histo2DRAAcomboPHENIX->Draw("axis,same");

        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_asPaper_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_asPaper_0010.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_asPaper_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));


        canvasRAAcombo->cd();
        histo2DRAAcomboPHENIX->DrawCopy();

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVPi0RAA_2040, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
            DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVPi0RAA_2040, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
            DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVPi0RAA_2040, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);

            graphPHENIX200GeVPi0RAA_2040->Draw("p,same");
            graphPHENIX39GeVPi0RAA_2040->Draw("p,same");
            graphPHENIX62GeVPi0RAA_2040->Draw("p,same");

            while(graphCombRAASysPbPb2760GeV_2050->GetX()[0]<1)graphCombRAASysPbPb2760GeV_2050->RemovePoint(0);
            while(graphCombRAAStatPbPb2760GeV_2050->GetX()[0]<1)graphCombRAAStatPbPb2760GeV_2050->RemovePoint(0);
            graphCombRAASysPbPb2760GeV_2050->Draw("E2same");
            graphCombRAAStatPbPb2760GeV_2050->Draw("p,same");

            TLatex *labelRAAALICEPbPb2040 = new TLatex(0.5,0.87,"#pi^{0} ALICE Pb#font[122]{-}Pb (2011)");
            SetStyleTLatex( labelRAAALICEPbPb2040, 0.85*textsizeLabelsRAA,4);
            labelRAAALICEPbPb2040->Draw();
            TLegend* legendRAASinglePbPb2040 = new TLegend(0.5,0.81,0.83,0.85);
            legendRAASinglePbPb2040->SetFillColor(0);
            legendRAASinglePbPb2040->SetLineColor(0);
            legendRAASinglePbPb2040->SetNColumns(1);
            legendRAASinglePbPb2040->SetTextFont(42);
            legendRAASinglePbPb2040->SetTextSize(0.85*textsizeLabelsRAA);
            legendRAASinglePbPb2040->SetMargin(0.17);
            legendRAASinglePbPb2040->AddEntry(graphCombRAASysPbPb2760GeV_2050,"20#font[122]{-}50% #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
            legendRAASinglePbPb2040->Draw();

            TLatex *labelRAAPHENIXPbPb2040 = new TLatex(0.5,0.75,"#pi^{0} PHENIX 20-40% Au#font[122]{-}Au");
            SetStyleTLatex( labelRAAPHENIXPbPb2040, 0.85*textsizeLabelsRAA,4);
            labelRAAPHENIXPbPb2040->Draw();

            TLegend* legendRAARHICPbPb2040 = new TLegend(0.5,0.55,0.84,0.73);
            legendRAARHICPbPb2040->SetFillColor(0);
            legendRAARHICPbPb2040->SetLineColor(0);
        //  legendRAARHICPbPb2040->SetNColumns(2);
            legendRAARHICPbPb2040->SetTextFont(42);
            legendRAARHICPbPb2040->SetTextSize(0.85*textsizeLabelsRAA);
            legendRAARHICPbPb2040->SetMargin(0.17);
            legendRAARHICPbPb2040->AddEntry(graphPHENIX200GeVPi0RAA_2040,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendRAARHICPbPb2040->AddEntry(graphPHENIX62GeVPi0RAA_2040,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
            legendRAARHICPbPb2040->AddEntry(graphPHENIX39GeVPi0RAA_2040,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
            legendRAARHICPbPb2040->Draw();

    //      boxErrorNorm2040_Single->Draw();
            boxErrorNorm2050_Single->Draw();
            DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray,2);

        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_asPaper_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_asPaper_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    } else if(meson.CompareTo("Eta")==0){

        canvasRAAcombo->cd();
        histo2DRAAcombo->GetYaxis()->SetRangeUser(0.,1.2);
        histo2DRAAcombo->DrawCopy();

            graphCombRAAStatPbPb2760GeV_0010->SetMarkerSize(2.7);
            graphCombRAASysPbPb2760GeV_0010->Draw("E2same");
            graphCombRAAStatPbPb2760GeV_0010->Draw("p,same");

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_0010, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
            graphPHENIX200GeVEtaRAA_0010->Draw("p,same");

            TLegend* legendEtaRAAcomp =new TLegend(0.5,0.83,0.8,0.93);
            legendEtaRAAcomp->SetFillColor(0);
            legendEtaRAAcomp->SetLineColor(0);
            legendEtaRAAcomp->SetTextFont(42);
            legendEtaRAAcomp->SetMargin(0.17);
            legendEtaRAAcomp->SetTextSize(0.85*textsizeLabelsRAA);
            legendEtaRAAcomp->SetHeader("#eta ALICE 0#font[122]{-}10% Pb#font[122]{-}Pb (2011)");
            legendEtaRAAcomp->AddEntry(graphCombRAASysPbPb2760GeV_0010,"#sqrt{#it{s}_{_{NN}}} = 2.76 TeV", "pf");
            legendEtaRAAcomp->Draw();
            TLegend* legendEtaRAAcompPHENIX1 =new TLegend(0.5,0.7,0.8,0.8);
            legendEtaRAAcompPHENIX1->SetFillColor(0);
            legendEtaRAAcompPHENIX1->SetLineColor(0);
            legendEtaRAAcompPHENIX1->SetTextFont(42);
            legendEtaRAAcompPHENIX1->SetMargin(0.17);
            legendEtaRAAcompPHENIX1->SetTextSize(0.85*textsizeLabelsRAA);
            legendEtaRAAcompPHENIX1->SetHeader("#eta PHENIX 0#font[122]{-}10% Au#font[122]{-}Au");
            legendEtaRAAcompPHENIX1->AddEntry(graphPHENIX200GeVEtaRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendEtaRAAcompPHENIX1->Draw();

            boxErrorNorm0010_Single->Draw();
            DrawGammaLines(0., 21 , 1, 1 ,1,kGray,2);
            histo2DRAAcombo->Draw("axis,same");

        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_0010.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_0010.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

        canvasRAAcombo->cd();
        histo2DRAAcombo->DrawCopy();
            DrawGammaLines(0., 21 , 1, 1 ,1,kGray,2);

            graphCombRAASysPbPb2760GeV_0010->Draw("E2same");
            graphCombRAAStatPbPb2760GeV_0010->Draw("p,same");

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_0020, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
            graphPHENIX200GeVEtaRAA_0020->Draw("p,same");

            legendEtaRAAcomp->Draw();
            TLegend* legendEtaRAAcompPHENIX2 =new TLegend(0.5,0.72,0.8,0.82);
            legendEtaRAAcompPHENIX2->SetFillColor(0);
            legendEtaRAAcompPHENIX2->SetLineColor(0);
            legendEtaRAAcompPHENIX2->SetTextFont(42);
            legendEtaRAAcompPHENIX2->SetMargin(0.17);
            legendEtaRAAcompPHENIX2->SetTextSize(0.85*textsizeLabelsRAA);
            legendEtaRAAcompPHENIX2->SetHeader("#eta PHENIX 0#font[122]{-}20% Au#font[122]{-}Au");
            legendEtaRAAcompPHENIX2->AddEntry(graphPHENIX200GeVEtaRAA_0020,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendEtaRAAcompPHENIX2->Draw();

            boxErrorNorm0010_Single->Draw();

        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_0010V2.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//         canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_0010V2.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_0010V2.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

        canvasRAAcombo->cd();
        histo2DRAAcombo->Draw("copy");

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_2060, 25,2, colorPhenix,colorPhenix, 0.1, kFALSE);
            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_2060,25,2, colorPhenix,colorPhenix, 1, kTRUE, colorPhenix);
            graphPHENIX200GeVEtaRAA_2060->SetFillStyle(3002);
            graphPHENIX200GeVEtaRAA_2060->Draw("2same");
            graphPHENIX200GeVEtaRAA_2060->Draw("p,same");

            graphCombRAASysPbPb2760GeV_2050->Draw("E2same");
            graphCombRAAStatPbPb2760GeV_2050->Draw("p,same");

            TLegend* legendEtaRAAcomp1 = new TLegend(0.48,0.73,0.8,0.93);
            legendEtaRAAcomp1->SetFillColor(0);
            legendEtaRAAcomp1->SetLineColor(0);
            legendEtaRAAcomp1->SetTextFont(42);
            legendEtaRAAcomp1->SetMargin(0.17);
            legendEtaRAAcomp1->SetTextSize(FontSize);
            legendEtaRAAcomp1->SetHeader("#eta - semicentral");
            legendEtaRAAcomp1->AddEntry(graphCombRAASysPbPb2760GeV_2050,"20#font[122]{-}50%, Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV", "pf");
            legendEtaRAAcomp1->AddEntry(graphPHENIX200GeVEtaRAA_2060,"20-60%, Au-Au #sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendEtaRAAcomp1->Draw();

            boxErrorNorm2050Only->Draw();
            DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);

        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_2050V2.%s",outputDir.Data(),meson.Data(),suffix.Data()));
        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_2050V2.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

        canvasRAAcombo->cd();
        histo2DRAAcombo->Draw("copy");

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_2040, 25,2, colorPhenix,colorPhenix, 0.1, kFALSE);
            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_2040,24,2, colorPhenix,colorPhenix, 1, kTRUE,colorPhenix);
            graphPHENIX200GeVEtaRAA_2040->SetFillStyle(3003);
            graphPHENIX200GeVEtaRAA_2040->Draw("2same");
            graphPHENIX200GeVEtaRAA_2040->Draw("p,same");

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaHighPtRAA_2060, 25,2, colorPhenix,colorPhenix, 0.1, kFALSE);
            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaHighPtRAA_2060,25,2, colorPhenix,colorPhenix, 1, kTRUE, colorPhenix);
            graphPHENIX200GeVEtaHighPtRAA_2060->SetFillStyle(3002);
            graphPHENIX200GeVEtaHighPtRAA_2060->Draw("2same");
            graphPHENIX200GeVEtaHighPtRAA_2060->Draw("p,same");

            graphCombRAASysPbPb2760GeV_2050->Draw("E2same");
            graphCombRAAStatPbPb2760GeV_2050->Draw("p,same");

            TLegend* legendEtaRAAcomp2 = new TLegend(0.48,0.73,0.8,0.93);
            legendEtaRAAcomp2->SetFillColor(0);
            legendEtaRAAcomp2->SetLineColor(0);
            legendEtaRAAcomp2->SetTextFont(42);
            legendEtaRAAcomp2->SetMargin(0.17);
            legendEtaRAAcomp2->SetTextSize(FontSize);
            legendEtaRAAcomp2->SetHeader("#eta - semicentral");
            legendEtaRAAcomp2->AddEntry(graphCombRAASysPbPb2760GeV_2050,"20#font[122]{-}50%, Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV", "pf");
            legendEtaRAAcomp2->AddEntry(graphPHENIX200GeVEtaRAA_2040,"20-40%, Au-Au #sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendEtaRAAcomp2->AddEntry(graphPHENIX200GeVEtaHighPtRAA_2060,"20-60%, Au-Au #sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendEtaRAAcomp2->Draw();

            boxErrorNorm2050Only->Draw();
            DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);

        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
        canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_2050.%s",PubNotePlots.Data(),meson.Data(),suffix.Data()));

    }


    //**********************************************************************************************************************//
    //********************************************   Eta to Pi0 combination    *********************************************//
    //**********************************************************************************************************************//
    if(meson.CompareTo("Eta")==0){

        //eta to pi0 ratio
        graphCombEtatoPi0TotPbPb2760GeV_0010 = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEtatoPi0LHC11h_0010,  sysErrorCollectionEtatoPi0LHC11h_0010,
                                                                                                    xPtLimitsEta, /*17*/13, offSetsEtaPi0Ratio, offSetsEtaPi0RatioSys,
                                                                                                    graphCombEtatoPi0StatPbPb2760GeV_0010, graphCombEtatoPi0SysPbPb2760GeV_0010, "weightEtatoPi0_0010.dat");
        graphCombEtatoPi0TotPbPb2760GeV_2050 = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEtatoPi0LHC11h_2050,  sysErrorCollectionEtatoPi0LHC11h_2050,
                                                                                                    xPtLimitsEta, /*17*/13, offSetsEtaPi0Ratio, offSetsEtaPi0RatioSys,
                                                                                                    graphCombEtatoPi0StatPbPb2760GeV_2050, graphCombEtatoPi0SysPbPb2760GeV_2050, "weightEtatoPi0_2050.dat");


        while(graphCombEtatoPi0TotPbPb2760GeV_0010->GetY()[0] == 0){
            graphCombEtatoPi0TotPbPb2760GeV_0010->RemovePoint(0);
            graphCombEtatoPi0StatPbPb2760GeV_0010->RemovePoint(0);
            graphCombEtatoPi0SysPbPb2760GeV_0010->RemovePoint(0);

            graphCombEtatoPi0TotPbPb2760GeV_2050->RemovePoint(0);
            graphCombEtatoPi0StatPbPb2760GeV_2050->RemovePoint(0);
            graphCombEtatoPi0SysPbPb2760GeV_2050->RemovePoint(0);

        }
//       graphCombEtatoPi0TotPbPb2760GeV_0010->Print();
//       graphCombEtatoPi0StatPbPb2760GeV_0010->Print();
//       graphCombEtatoPi0SysPbPb2760GeV_0010->Print();
//       graphCombEtatoPi0TotPbPb2760GeV_2050->Print();
//       graphCombEtatoPi0StatPbPb2760GeV_2050->Print();
//       graphCombEtatoPi0SysPbPb2760GeV_2050->Print();
//       return;


        TCanvas* canvasEtatoPi0combo = new TCanvas("canvasEtatoPi0combo","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.09, 0.03, 0.035, 0.1);
        canvasEtatoPi0combo->SetLogx();

            TH2F * histo2DEtatoPi0combo = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",11000,minPtRange,maxPtRange,1000,0.01,1.2);
            SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}",0.035,0.04, 0.035,0.04, 1.,1.);
            histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
            histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
            histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.01,1.2);
            histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphPCMEtatoPi0Sys2760GeV_0010, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
            graphPCMEtatoPi0Sys2760GeV_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphPCMEtatoPi0Stat2760GeV_0010, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphPCMEtatoPi0Stat2760GeV_0010);
            graphPCMEtatoPi0Stat2760GeV_0010->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Sys2760GeV_0010,markerStyleDet[2] ,markerStyleDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
            graphEMCalEtatoPi0Sys2760GeV_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Stat2760GeV_0010, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphEMCalEtatoPi0Stat2760GeV_0010);
            graphEMCalEtatoPi0Stat2760GeV_0010->Draw("p,same");

            // ****************************** Definition of the Legend ******************************************
            // **************** Row def ************************
            Double_t rowsLegendOnlyEtaToPi0Ratio[3]      = {0.85,0.8,0.77};
            Double_t rowsLegendOnlyEtaToPi0RatioAbs[3]   = {0.95,1.,.95};
            Double_t columnsLegendOnlyEtaToPi0Ratio[3]   = {0.15,0.4, 0.45};
            Double_t columnsLegendOnlyEtaToPi0RatioAbs[3]= {0.15,2.,2.6};
            Double_t lengthBox3                      = 0.3/2.;
            Double_t heightBox3                      = 0.08/4;
            textSizeLabelsPixel = 40;
            // ****************** first Column **************************************************
            TLatex *textPCMOnlyEtaToPi0RatioLHC11h = new TLatex(columnsLegendOnlyEtaToPi0Ratio[0],rowsLegendOnlyEtaToPi0Ratio[1],"PCM");//" (this thesis)");
            SetStyleTLatex( textPCMOnlyEtaToPi0RatioLHC11h, 0.85*textSizeLabelsPixel,4);
            textPCMOnlyEtaToPi0RatioLHC11h->SetTextFont(43);
            textPCMOnlyEtaToPi0RatioLHC11h->Draw();
            TLatex *textEMCalOnlyEtaToPi0RatioLHC11h = new TLatex(columnsLegendOnlyEtaToPi0Ratio[0],rowsLegendOnlyEtaToPi0Ratio[2],"EMCal");
            SetStyleTLatex( textEMCalOnlyEtaToPi0RatioLHC11h,0.85*textSizeLabelsPixel,4);
            textEMCalOnlyEtaToPi0RatioLHC11h->SetTextFont(43);
            textEMCalOnlyEtaToPi0RatioLHC11h->Draw();
            // ****************** second Column *************************************************
            TLatex *textStatOnlyEtaToPi0RatioLHC11h = new TLatex(columnsLegendOnlyEtaToPi0Ratio[1],rowsLegendOnlyEtaToPi0Ratio[0] ,"stat");
            SetStyleTLatex( textStatOnlyEtaToPi0RatioLHC11h, 0.85*textSizeLabelsPixel,4);
            textStatOnlyEtaToPi0RatioLHC11h->SetTextFont(43);
            textStatOnlyEtaToPi0RatioLHC11h->Draw();
            TLatex *textSysOnlyEtaToPi0RatioLHC11h = new TLatex(columnsLegendOnlyEtaToPi0Ratio[2] ,rowsLegendOnlyEtaToPi0Ratio[0],"syst");
            SetStyleTLatex( textSysOnlyEtaToPi0RatioLHC11h, 0.85*textSizeLabelsPixel,4);
            textSysOnlyEtaToPi0RatioLHC11h->SetTextFont(43);
            textSysOnlyEtaToPi0RatioLHC11h->Draw();

            TMarker* markerPCMPi0OnlyRatioPi0LHC11h = CreateMarkerFromGraph(graphPCMEtatoPi0Sys2760GeV_0010,columnsLegendOnlyEtaToPi0Ratio[1] ,rowsLegendOnlyEtaToPi0Ratio[1],0.7);
            markerPCMPi0OnlyRatioPi0LHC11h->DrawMarker(columnsLegendOnlyEtaToPi0RatioAbs[1] ,rowsLegendOnlyEtaToPi0RatioAbs[1]);
            TBox* boxPCMPi0OnlyRatioPi0 = CreateBoxFromGraph(graphPCMEtatoPi0Sys2760GeV_0010, columnsLegendOnlyEtaToPi0RatioAbs[2]-1.5*lengthBox3 , rowsLegendOnlyEtaToPi0RatioAbs[1]- heightBox3, columnsLegendOnlyEtaToPi0RatioAbs[2]+ 2*lengthBox3, rowsLegendOnlyEtaToPi0RatioAbs[1]+ heightBox3);
            boxPCMPi0OnlyRatioPi0->Draw("l");

            TMarker* markerEMCalPi0OnlyRatioPi0LHC11h = CreateMarkerFromGraph(graphEMCalEtatoPi0Sys2760GeV_0010, columnsLegendOnlyEtaToPi0Ratio[1] ,rowsLegendOnlyEtaToPi0Ratio[2],0.2);
            markerEMCalPi0OnlyRatioPi0LHC11h->DrawMarker(columnsLegendOnlyEtaToPi0RatioAbs[1] ,rowsLegendOnlyEtaToPi0RatioAbs[2]);
            TBox* boxEMCalPi0OnlyRatioPi0 = CreateBoxFromGraph(graphEMCalEtatoPi0Sys2760GeV_0010, columnsLegendOnlyEtaToPi0RatioAbs[2]-1.5*lengthBox3 , rowsLegendOnlyEtaToPi0RatioAbs[2]- heightBox3, columnsLegendOnlyEtaToPi0RatioAbs[2]+ 2*lengthBox3, rowsLegendOnlyEtaToPi0RatioAbs[2]+ heightBox3);
            boxEMCalPi0OnlyRatioPi0->Draw("l");

            TLatex *labelEtaToPi0RatioEnergy = new TLatex(0.15,0.9,collisionSystemPbPb0010.Data());
            SetStyleTLatex( labelEtaToPi0RatioEnergy, 0.85*textSizeLabelsPixel,4);
            labelEtaToPi0RatioEnergy->SetTextFont(43);
            labelEtaToPi0RatioEnergy->Draw();

            TLatex *textEMCalNoteEtaToPi0 = new TLatex(0.65,0.9,"(*) ALICE preliminary input");
            SetStyleTLatex( textEMCalNoteEtaToPi0,  0.8*textSizeLabelsPixel,4);
            textEMCalNoteEtaToPi0->SetTextFont(43);
            TLatex *textEMCalNoteEtaToPi01 = new TLatex(0.65,0.86,"arXiv:1609.06106");
            SetStyleTLatex( textEMCalNoteEtaToPi01,  0.8*textSizeLabelsPixel,4);
            textEMCalNoteEtaToPi01->SetTextFont(43);
            TLatex *textEMCalNoteEtaToPi02 = new TLatex(0.65,0.82,"analysis by A. Morreale");
            SetStyleTLatex( textEMCalNoteEtaToPi02,  0.8*textSizeLabelsPixel,4);
            textEMCalNoteEtaToPi02->SetTextFont(43);
//             textEMCalNoteEtaToPi0->Draw();
//             textEMCalNoteEtaToPi01->Draw();
//             textEMCalNoteEtaToPi02->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_IndividualMeas_0010.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_IndividualMeas_0010.%s",PubNotePlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphPCMEtatoPi0Sys2760GeV_2050, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
            graphPCMEtatoPi0Sys2760GeV_2050->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphPCMEtatoPi0Stat2760GeV_2050, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphPCMEtatoPi0Stat2760GeV_2050);
            graphPCMEtatoPi0Stat2760GeV_2050->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Sys2760GeV_2050,markerStyleDet[2] ,markerStyleDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
            graphEMCalEtatoPi0Sys2760GeV_2050->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Stat2760GeV_2050, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphEMCalEtatoPi0Stat2760GeV_2050);
            graphEMCalEtatoPi0Stat2760GeV_2050->Draw("p,same");

            TLatex *labelEtaToPi0RatioEnergy2 = new TLatex(0.15,0.9,collisionSystemPbPb2050.Data());
            SetStyleTLatex( labelEtaToPi0RatioEnergy2, 0.85*textSizeLabelsPixel,4);
            labelEtaToPi0RatioEnergy2->SetTextFont(43);
            labelEtaToPi0RatioEnergy2->Draw();

            textPCMOnlyEtaToPi0RatioLHC11h->Draw();
            textEMCalOnlyEtaToPi0RatioLHC11h->Draw();
            textStatOnlyEtaToPi0RatioLHC11h->Draw();
            textSysOnlyEtaToPi0RatioLHC11h->Draw();

            TMarker* markerPCMEtaToPi0RatioLHC11h_2050 = CreateMarkerFromGraph(graphPCMEtatoPi0Sys2760GeV_2050,columnsLegendOnlyEtaToPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],0.7);
            markerPCMEtaToPi0RatioLHC11h_2050->DrawMarker(columnsLegendOnlyEtaToPi0RatioAbs[1] ,rowsLegendOnlyEtaToPi0RatioAbs[1]);
            TBox* boxPCMEtaToPi0Ratio_2050 = CreateBoxFromGraph(graphPCMEtatoPi0Sys2760GeV_2050, columnsLegendOnlyEtaToPi0RatioAbs[2]-1.5*lengthBox3 , rowsLegendOnlyEtaToPi0RatioAbs[1]- heightBox3, columnsLegendOnlyEtaToPi0RatioAbs[2]+ 2*lengthBox3, rowsLegendOnlyEtaToPi0RatioAbs[1]+ heightBox3);
            boxPCMEtaToPi0Ratio_2050->Draw("l");

            TMarker* markerEMCalEtaToPi0RatioLHC11h_2050 = CreateMarkerFromGraph(graphEMCalEtatoPi0Sys2760GeV_2050, columnsLegendOnlyEtaToPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[2],0.2);
            markerEMCalEtaToPi0RatioLHC11h_2050->DrawMarker(columnsLegendOnlyEtaToPi0RatioAbs[1] ,rowsLegendOnlyEtaToPi0RatioAbs[2]);
            TBox* boxEMCalEtaToPi0Ratio_2050 = CreateBoxFromGraph(graphEMCalEtatoPi0Sys2760GeV_2050, columnsLegendOnlyEtaToPi0RatioAbs[2]-1.5*lengthBox3 , rowsLegendOnlyEtaToPi0RatioAbs[2]- heightBox3, columnsLegendOnlyEtaToPi0RatioAbs[2]+ 2*lengthBox3, rowsLegendOnlyEtaToPi0RatioAbs[2]+ heightBox3);
            boxEMCalEtaToPi0Ratio_2050->Draw("l");

//             textEMCalNoteEtaToPi0->Draw();
//             textEMCalNoteEtaToPi01->Draw();
//             textEMCalNoteEtaToPi02->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_IndividualMeas_2050.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_IndividualMeas_2050.%s",PubNotePlots.Data(),suffix.Data()));


//       canvasEtatoPi0combo->cd();
//           histo2DEtatoPi0combo->Draw("copy");
//           histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(3.,30.);
//
//             DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Sys2760GeV_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
//             DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Stat2760GeV_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);
//             graphEMCalEtatoPi0Sys2760GeV_0010->Draw("E2same");
//             graphEMCalEtatoPi0Stat2760GeV_0010->Draw("p,same");
//
//             DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Sys2760GeV_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050, 2, kTRUE);
//             DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Stat2760GeV_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050);
//             graphEMCalEtatoPi0Sys2760GeV_2050->Draw("E2same");
//             graphEMCalEtatoPi0Stat2760GeV_2050->Draw("p,same");
//
//             TLegend* legendEtatoPi0EMCalOnly = new TLegend(0.12,0.76,0.53,0.92);
//             legendEtatoPi0EMCalOnly->SetFillColor(0);
//             legendEtatoPi0EMCalOnly->SetLineColor(0);
//             legendEtatoPi0EMCalOnly->SetTextFont(42);
//             legendEtatoPi0EMCalOnly->SetTextSize(0.037);
//             legendEtatoPi0EMCalOnly->SetMargin(0.17);
//             legendEtatoPi0EMCalOnly->SetHeader(Form("EMCal only, %s",collisionSystem2760GeV.Data()));
//             legendEtatoPi0EMCalOnly->AddEntry(graphEMCalEtatoPi0Sys2760GeV_0010,Form("  %s",cent0010.Data()),"fp");
//             legendEtatoPi0EMCalOnly->AddEntry(graphEMCalEtatoPi0Sys2760GeV_2050,Form("%s",cent2050.Data()),"fp");
//             legendEtatoPi0EMCalOnly->Draw();
//
//       canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_EMCalOnly.%s",outputDir.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0SysPbPb2760GeV_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0StatPbPb2760GeV_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_0010);
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0SysPbPb2760GeV_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050, 2, kTRUE);
            graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0StatPbPb2760GeV_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050);
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_2050);
            graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

            TF1* linEtaToPi0 = new TF1("linEtaToPi0","[0]",minPtRange,maxPtRange);
            linEtaToPi0->SetLineColor(kBlack);
            graphCombEtatoPi0StatPbPb2760GeV_0010->Fit(linEtaToPi0,"NRMEX0+","",3,20.);
            graphCombEtatoPi0SysPbPb2760GeV_0010->Fit(linEtaToPi0,"NRMEX0+","",3,20.);
//             linEtaToPi0->Draw("same");

            TLegend* legendEtatoPi0combo_onlyPbPb = new TLegend(0.12,0.76,0.53,0.92);
            legendEtatoPi0combo_onlyPbPb->SetFillColor(0);
            legendEtatoPi0combo_onlyPbPb->SetLineColor(0);
            legendEtatoPi0combo_onlyPbPb->SetTextFont(42);
            legendEtatoPi0combo_onlyPbPb->SetTextSize(0.037);
            legendEtatoPi0combo_onlyPbPb->SetMargin(0.17);
            legendEtatoPi0combo_onlyPbPb->SetHeader(collisionSystem2760GeV.Data());
            legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,Form("  %s",cent0010.Data()),"fp");
            legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("%s",cent2050.Data()),"fp");
            legendEtatoPi0combo_onlyPbPb->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_DataOnly.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_DataOnly.%s",paperPlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaToPi0Ratio_0020, markerStylePHENIX200GeV, markerSizePHENIX200GeV, kBlack , kBlack);
            graphPHENIX200GeVEtaToPi0Ratio_0020->Draw("p,same");

            TLegend* legendEtatoPi0combo_onlyPbPbWithPHENIX = new TLegend(0.12,0.76,0.53,0.92);
            legendEtatoPi0combo_onlyPbPbWithPHENIX->SetFillColor(0);
            legendEtatoPi0combo_onlyPbPbWithPHENIX->SetLineColor(0);
            legendEtatoPi0combo_onlyPbPbWithPHENIX->SetTextFont(42);
            legendEtatoPi0combo_onlyPbPbWithPHENIX->SetTextSize(0.037);
            legendEtatoPi0combo_onlyPbPbWithPHENIX->SetMargin(0.17);
            legendEtatoPi0combo_onlyPbPbWithPHENIX->SetHeader("This thesis"); //collisionSystem2760GeV.Data());
            legendEtatoPi0combo_onlyPbPbWithPHENIX->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,Form("ALICE %s",collisionSystemPbPb0010.Data()),"fp");
            legendEtatoPi0combo_onlyPbPbWithPHENIX->AddEntry(graphPHENIX200GeVEtaToPi0Ratio_0020,"PHENIX 0#font[122]{-}20% Au#font[122]{-}Au, #sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
            legendEtatoPi0combo_onlyPbPbWithPHENIX->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_DataWithPHENIX_0010.%s",outputDir.Data(),suffix.Data()));
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_DataWithPHENIX_0010.%s",paperPlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

            DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaToPi0Ratio_2060, markerStylePHENIX200GeV, markerSizePHENIX200GeV, kBlack , kBlack);
            graphPHENIX200GeVEtaToPi0Ratio_2060->Draw("p,same");

            TLegend* legendEtatoPi0combo_onlyPbPbWithPHENIX2 = new TLegend(0.12,0.76,0.53,0.92);
            legendEtatoPi0combo_onlyPbPbWithPHENIX2->SetFillColor(0);
            legendEtatoPi0combo_onlyPbPbWithPHENIX2->SetLineColor(0);
            legendEtatoPi0combo_onlyPbPbWithPHENIX2->SetTextFont(42);
            legendEtatoPi0combo_onlyPbPbWithPHENIX2->SetTextSize(0.037);
            legendEtatoPi0combo_onlyPbPbWithPHENIX2->SetMargin(0.17);
            legendEtatoPi0combo_onlyPbPbWithPHENIX2->SetHeader("This thesis"); //collisionSystem2760GeV.Data());
            legendEtatoPi0combo_onlyPbPbWithPHENIX2->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("ALICE %s",collisionSystemPbPb2050.Data()),"fp");
            legendEtatoPi0combo_onlyPbPbWithPHENIX2->AddEntry(graphPHENIX200GeVEtaToPi0Ratio_2060,"PHENIX 20#font[122]{-}60% Au#font[122]{-}Au, #sqrt{#it{s}_{_{NN}}}= 200 GeV","p");
            legendEtatoPi0combo_onlyPbPbWithPHENIX2->Draw();


        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_DataWithPHENIX_2050.%s",outputDir.Data(),suffix.Data()));
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_DataWithPHENIX_2050.%s",paperPlots.Data(),suffix.Data()));


        TH1D *dummy = new TH1D("dummy","dummy",13, xPt);
        TFile *ppInput7TeV = new TFile("/home/admin1/leardini/newSoftware/cocktail_input/pp/pp7TeV_data_PCMResultsFullCorrection_PP_NoBinShifting_current20161117.root");
        TDirectoryFile *Pi07TeV = (TDirectoryFile*)ppInput7TeV->Get("Pi07TeV");
        TDirectoryFile *Pi0Fits = (TDirectoryFile*)Pi07TeV->Get("Fits");
        TF1 *pi0LevyFit7TeV = (TF1*)Pi0Fits->Get("fitPtLevyPi0");
        TF1 *mTScaledEtaFromPi07TeV = (TF1*)MtScaledParam(pi0LevyFit7TeV, 221, 0.476);
        TF1* etapi0Ratio7TeV  = DivideTF1(mTScaledEtaFromPi07TeV, pi0LevyFit7TeV, "etapi0Ratio7TeV");

        TFile *ppOldInput7TeV = new TFile("/home/admin1/leardini/newSoftware/cocktail_input/pp/pp7TeV_data_PCMResultsFullCorrection_PP_NoBinShifting_current20161117.root");
        TDirectoryFile *Pi0Old7TeV = (TDirectoryFile*)ppOldInput7TeV->Get("Pi07TeV");
        TDirectoryFile *Pi0OldFits = (TDirectoryFile*)Pi0Old7TeV->Get("Fits");
        TF1 *pi0OldLevyFit7TeV = (TF1*)Pi0OldFits->Get("fitPtLevyPi0");
        TF1 *mTScaledEtaFromOldPi07TeV = (TF1*)MtScaledParam(pi0OldLevyFit7TeV, 221, 0.476);
        TF1* etapi0RatioOld7TeV  = DivideTF1(mTScaledEtaFromOldPi07TeV, pi0OldLevyFit7TeV, "etapi0RatioOld7TeV");

        TF1 *mTScaledEtaFromPi02760GeV = (TF1*)MtScaledParam(fitInvCrossSectionTsallisPi0Comb2760GeV, 221, 0.476);
        TF1* etapi0Ratio2760GeV  = DivideTF1(mTScaledEtaFromPi02760GeV, fitInvCrossSectionTsallisPi0Comb2760GeV, "etapi0Ratio2760GeV");

//       cout << "PCM Eta for pPb" << endl;
//       TString nameFilepPb = "LHC11hExternalInputs/data_PCMResults_pPb_20150624_standard_dc4.root";
//       TFile* filePCMpPb                   = new TFile(nameFilepPb.Data());
//       TDirectory* fEtaToPi0pPbContainer   = (TDirectory*) filePCMpPb->GetDirectory("Eta_pPb_5.023TeV_0#font[122]{-}100%");
//       TH1D* histoPCMEtaToPi0RatiopPb      = (TH1D*)fEtaToPi0pPbContainer->Get("EtatoPi0Ratio");
//       TGraphAsymmErrors* graphPCMEtaToPi0RatioSysErrpPb=    (TGraphAsymmErrors*)fEtaToPi0pPbContainer->Get("EtatoPi0RatioSys");
//       DrawGammaSetMarker(histoPCMEtaToPi0RatiopPb, 20, 1.5, kBlue+2, kBlue+2);
//       DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioSysErrpPb, 21, 1.5,  kBlue+2, kBlue+2, widthLinesBoxes, kTRUE);

        TFile *PbPbCoktailInput = new TFile("50100013_00200009247602008850404000_0652501500000000/PbPb_2.76TeV/GammaCocktail_0.85_50100013_00200009247602008850404000_0652501500000000.root");
        TF1* paramPi0PbPb2760GeV     = (TF1*)PbPbCoktailInput->Get("111_pt");
        TF1* paramEtaPbPb2760GeV  = (TF1*)PbPbCoktailInput->Get("221_pt");
        TF1* etapi0RatioFromParamPbPb2760GeV  = DivideTF1(paramEtaPbPb2760GeV, paramPi0PbPb2760GeV, "etapi0RatioPbPb2760GTeV");

        TLatex *thesisLabelEtaPi0Ratio = new TLatex(0.12,0.9,thisthesis.Data());
        SetStyleTLatex( thesisLabelEtaPi0Ratio,FontSize,4);
        if(MesonInput){
            canvasEtatoPi0combo->cd();
                histo2DEtatoPi0combo->Draw("copy");

                cout << __LINE__ << endl;
                TF1 *pi0fitPbPb_0010 = (TF1*)MesonInput->Get("FitToYieldPi0_0010");
                TF1 *mTScaledEtaFromPi0_0010 = (TF1*)MesonInput->Get("mTScaledEtaFromPi0_0010");
                TF1* etapi0RatioPbPb2760GeV_0010  = DivideTF1(mTScaledEtaFromPi0_0010, pi0fitPbPb_0010, "etapi0RatioPbPb2760GeV_0010");

                graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_0010);
                graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

                etapi0RatioFromMtScalingCombPbPb2760GeV_0010->SetLineColor(colorCombo0010+2);
                etapi0RatioFromMtScalingCombPbPb2760GeV_0010->SetLineStyle(2);
                etapi0RatioFromMtScalingCombPbPb2760GeV_0010->SetLineWidth(4);
                etapi0RatioFromMtScalingCombPbPb2760GeV_0010->SetRange(.95,30.);
                etapi0RatioFromMtScalingCombPbPb2760GeV_0010->Draw("l,same");
//                 etapi0RatioPbPb2760GeV_0010->SetLineColor(colorCombo0010+2);
//                 etapi0RatioPbPb2760GeV_0010->SetLineStyle(2);
//                 etapi0RatioPbPb2760GeV_0010->SetLineWidth(4);
//                 etapi0RatioPbPb2760GeV_0010->Draw("c,histo,same");

                TLegend* legendEtatoPi0combo_onlyPbPb = new TLegend(0.12,0.7,0.53,0.88);
                legendEtatoPi0combo_onlyPbPb->SetFillColor(0);
                legendEtatoPi0combo_onlyPbPb->SetLineColor(0);
                legendEtatoPi0combo_onlyPbPb->SetTextFont(42);
                legendEtatoPi0combo_onlyPbPb->SetTextSize(0.04);
                legendEtatoPi0combo_onlyPbPb->SetMargin(0.17);
                legendEtatoPi0combo_onlyPbPb->SetHeader(collisionSystem2760GeV.Data());
                legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,Form("  %s",cent0010.Data()),"fp");
    //             legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("%s",cent2050.Data()),"fp");
                legendEtatoPi0combo_onlyPbPb->AddEntry(etapi0RatioFromMtScalingCombPbPb2760GeV_0010,"#eta from #it{m}_{T} scaled #pi^{0}","l");
                legendEtatoPi0combo_onlyPbPb->Draw();

                if(thesisPlotting) thesisLabelEtaPi0Ratio->Draw();

            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_mtscaled_DataOnly_0010.%s",outputDir.Data(),suffix.Data()));
            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_mtscaled_DataOnly_0010.%s",PubNotePlots.Data(),suffix.Data()));

            canvasEtatoPi0combo->cd();
                histo2DEtatoPi0combo->Draw("copy");

                cout << __LINE__ << endl;
                TF1 *pi0fitPbPb_2050 = (TF1*)MesonInput->Get("FitToYieldPi0_2050");
                TF1 *mTScaledEtaFromPi0_2050 = (TF1*)MesonInput->Get("mTScaledEtaFromPi0_0010");
                TF1* etapi0RatioPbPb2760GeV_2050  = DivideTF1(mTScaledEtaFromPi0_2050, pi0fitPbPb_2050, "etapi0RatioPbPb2760GeV_2050");

                graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_2050);
                graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

//                 etapi0RatioPbPb2760GeV_2050->SetLineColor(colorCombo0010+2);
//                 etapi0RatioPbPb2760GeV_2050->SetLineStyle(2);
//                 etapi0RatioPbPb2760GeV_2050->SetLineWidth(4);
//                 etapi0RatioPbPb2760GeV_2050->Draw("c,histo,same");
                etapi0RatioFromMtScalingCombPbPb2760GeV_2050->SetLineColor(colorCombo2050+2);
                etapi0RatioFromMtScalingCombPbPb2760GeV_2050->SetLineStyle(2);
                etapi0RatioFromMtScalingCombPbPb2760GeV_2050->SetLineWidth(4);
                etapi0RatioFromMtScalingCombPbPb2760GeV_2050->SetRange(.95,30.);
                etapi0RatioFromMtScalingCombPbPb2760GeV_2050->Draw("l,same");

                TLegend* legendEtatoPi0combo_onlyPbPbSC = new TLegend(0.12,0.7,0.53,0.88);
                legendEtatoPi0combo_onlyPbPbSC->SetFillColor(0);
                legendEtatoPi0combo_onlyPbPbSC->SetLineColor(0);
                legendEtatoPi0combo_onlyPbPbSC->SetTextFont(42);
                legendEtatoPi0combo_onlyPbPbSC->SetTextSize(0.04);
                legendEtatoPi0combo_onlyPbPbSC->SetMargin(0.17);
                legendEtatoPi0combo_onlyPbPbSC->SetHeader(collisionSystem2760GeV.Data());
                legendEtatoPi0combo_onlyPbPbSC->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("  %s",cent2050.Data()),"fp");
    //             legendEtatoPi0combo_onlyPbPbSC->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("%s",cent2050.Data()),"fp");
                legendEtatoPi0combo_onlyPbPbSC->AddEntry(etapi0RatioFromMtScalingCombPbPb2760GeV_2050,"#eta from #it{m}_{T} scaled #pi^{0}","l");
                legendEtatoPi0combo_onlyPbPbSC->Draw();

            if(thesisPlotting)          thesisLabelEtaPi0Ratio->Draw();

            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_mtscaled_DataOnly_2050.%s",outputDir.Data(),suffix.Data()));
            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_mtscaled_DataOnly_2050.%s",PubNotePlots.Data(),suffix.Data()));

            canvasEtatoPi0combo->cd();
                histo2DEtatoPi0combo->Draw("copy");

                if(thesisPlotting) thesisLabelEtaPi0Ratio->Draw();

                graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_0010);
                graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

//                 etapi0RatioPbPb2760GeV_0010->Draw("c,histo,same");
                legendEtatoPi0combo_onlyPbPb->Draw();
                etapi0RatioFromMtScalingCombPbPb2760GeV_0010->Draw("l,same");

                graphRatioEtaToPi0Comb2760GeVStatErr->Draw("p,same");
                graphRatioEtaToPi0Comb2760GeVSysErr->Draw("E2same");

                etapi0Ratio2760GeV->SetLineColor(kBlack);
                etapi0Ratio2760GeV->SetLineStyle(3);
                etapi0Ratio2760GeV->SetLineWidth(2);
                etapi0Ratio2760GeV->Draw("c,histo,same");

                TLegend* legendEtatoPi0combo_withPP2760GeV = new TLegend(0.55,0.15,0.95,0.3/*26*/);
                legendEtatoPi0combo_withPP2760GeV->SetFillColor(0);
                legendEtatoPi0combo_withPP2760GeV->SetLineColor(0);
                legendEtatoPi0combo_withPP2760GeV->SetTextFont(42);
                legendEtatoPi0combo_withPP2760GeV->SetTextSize(0.037);
                legendEtatoPi0combo_withPP2760GeV->SetMargin(0.17);
                legendEtatoPi0combo_withPP2760GeV->SetHeader(collisionSystemPP2760GeV.Data());
                legendEtatoPi0combo_withPP2760GeV->AddEntry(graphRatioEtaToPi0Comb2760GeVSysErr,"EPJC 77 (2017) 339"/*"arXiv: 1702.00917"*/,"fp");
                legendEtatoPi0combo_withPP2760GeV->AddEntry(etapi0Ratio2760GeV,"#eta from #it{m}_{T} scaled #pi^{0}","l");
                legendEtatoPi0combo_withPP2760GeV->Draw();

            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_mtscaled_DataOnlyWithPP2760GeV.%s",outputDir.Data(),suffix.Data()));
            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_mtscaled_DataOnlyWithPP2760GeV.%s",PubNotePlots.Data(),suffix.Data()));

            canvasEtatoPi0combo->cd();
                histo2DEtatoPi0combo->Draw("copy");

                graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
                if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_0010);
                graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");
                etapi0RatioPbPb2760GeV_0010->Draw("c,histo,same");

                graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
                graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");

                etapi0RatioOld7TeV->SetLineColor(kBlack);
                etapi0RatioOld7TeV->SetLineStyle(3);
                etapi0RatioOld7TeV->SetLineWidth(2);
                etapi0RatioOld7TeV->Draw("c,histo,same");

    //             etapi0Ratio7TeV->SetLineColor(kBlack);
    //             etapi0Ratio7TeV->SetLineStyle(3);
    //             etapi0Ratio7TeV->SetLineWidth(2);
    //             etapi0Ratio7TeV->Draw("c,histo,same");

    //             DrawGammaSetMarker(cocktailEtaToPi0Ratio_MtScaledRebinned, 2, 0, kBlack, kBlack);
    //             cocktailEtaToPi0Ratio7TeVRebined->GetXaxis()->SetRangeUser(0.4,20);
    //             cocktailEtaToPi0Ratio_K0ScaledRebinned->GetXaxis()->SetRangeUser(0.4,20);
    //             cocktailEtaToPi0Ratio_MtScaledRebinned->GetXaxis()->SetRangeUser(0.4,20);
    //             cocktailEtaToPi0Ratio7TeVRebined->Draw("same,hist,c");
    //             cocktailEtaToPi0Ratio_K0ScaledRebinned->Draw("same,hist,c");
    //             cocktailEtaToPi0Ratio_MtScaledRebinned->Draw("same,hist,c");

                legendEtatoPi0combo_onlyPbPb->Draw();

                if(thesisPlotting) thesisLabelHighRight->Draw();

                TLegend* legendEtatoPi0combo_withPP = new TLegend(0.55,0.15,0.95,0.3/*26*/);
                legendEtatoPi0combo_withPP->SetFillColor(0);
                legendEtatoPi0combo_withPP->SetLineColor(0);
                legendEtatoPi0combo_withPP->SetTextFont(42);
                legendEtatoPi0combo_withPP->SetTextSize(0.037);
                legendEtatoPi0combo_withPP->SetMargin(0.17);
                legendEtatoPi0combo_withPP->SetHeader(collisionSystemPP7TeV.Data());
                legendEtatoPi0combo_withPP->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,"PLB 717 (2012) 162","fp");//"Phys. Lett. B 717 (2012) 162-172","fp");
    //             legendEtatoPi0combo_withPP->AddEntry(etapi0RatioOld7TeV,"#eta from #it{m}_{T} scaled #pi^{0}","l");
                legendEtatoPi0combo_withPP->Draw();

            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_mtscaled_DataOnlyWithPP7TeV.%s",outputDir.Data(),suffix.Data()));
            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_mtscaled_DataOnlyWithPP7TeV.%s",PubNotePlots.Data(),suffix.Data()));


            canvasEtatoPi0combo->cd();
                histo2DEtatoPi0combo->Draw("copy");
                histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.,1.15);

                graphChargedRatioKaonToPionSys0010->Draw("2same");
                graphChargedRatioKaonToPion0010->Draw("p,same");

                graphRatioEtaToPi0Comb2760GeVStatErr->Draw("p,same");
                graphRatioEtaToPi0Comb2760GeVSysErr->Draw("E2same");

//                 etapi0RatioPbPb2760GeV_0010->Draw("c,histo,same");
                etapi0RatioFromMtScalingCombPbPb2760GeV_0010->Draw("l,same");

                TLegend* legendChargedRatio3 = new TLegend(0.12,0.63,0.48,0.88);
                legendChargedRatio3->SetFillColor(0);
                legendChargedRatio3->SetLineColor(0);
                legendChargedRatio3->SetTextFont(42);
                legendChargedRatio3->SetTextSize(0.037);
                legendChargedRatio3->SetMargin(0.17);
                legendChargedRatio3->SetHeader(collisionSystemPbPb0010.Data());
                legendChargedRatio3->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,"#eta/#pi^{0}","fp");
                legendChargedRatio3->AddEntry(etapi0RatioFromMtScalingCombPbPb2760GeV_0010,"#eta from #it{m}_{T} scaled #pi^{0}","l");
                legendChargedRatio3->AddEntry(graphChargedRatioKaonToPionSys0010,"K^{#pm}/#pi^{#pm}","fp");
                legendChargedRatio3->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
                legendChargedRatio3->Draw();

                TLegend* legendEtatoPi0combo_withPP2760GeVmt = GetAndSetLegend(0.6,0.15,2);// new TLegend(0.55,0.15,0.95,0.3/*26*/);
                legendEtatoPi0combo_withPP2760GeVmt->SetFillColor(0);
                legendEtatoPi0combo_withPP2760GeVmt->SetLineColor(0);
                legendEtatoPi0combo_withPP2760GeVmt->SetTextFont(42);
                legendEtatoPi0combo_withPP2760GeVmt->SetTextSize(0.037);
                legendEtatoPi0combo_withPP2760GeVmt->SetMargin(0.17);
                legendEtatoPi0combo_withPP2760GeVmt->SetHeader(collisionSystemPP2760GeV.Data());
                legendEtatoPi0combo_withPP2760GeVmt->AddEntry(graphRatioEtaToPi0Comb2760GeVSysErr,"EPJC 77 (2017) 339"/*"arXiv: 1702.00917"*/,"fp");
    //             legendEtatoPi0combo_withPP2760GeVmt->AddEntry(etapi0Ratio2760GeV,"#eta from #it{m}_{T} scaled #pi^{0}","l");
                legendEtatoPi0combo_withPP2760GeVmt->Draw();

//                 if(thesisPlotting)          thesisLabelEtaPi0Ratio->Draw();

                graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
                graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_withPP2760GeVandKaonsToPionsMt_0010.%s",outputDir.Data(),suffix.Data()));
//             canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_withPP2760GeVandKaonsToPionsMt_0010.%s",paperPlots.Data(),suffix.Data()));

            canvasEtatoPi0combo->cd();
                histo2DEtatoPi0combo->Draw("copy");
                histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.,1.15);

                DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPion2040, 25,markerSizeComb, colorCharged, colorCharged);
                DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPionSys2040,25,markerSizeComb, colorCharged, colorCharged, 1, kTRUE);
                graphChargedRatioKaonToPionSys2040->Draw("2same");
                graphChargedRatioKaonToPion2040->Draw("p,same");

                graphRatioEtaToPi0Comb2760GeVStatErr->Draw("p,same");
                graphRatioEtaToPi0Comb2760GeVSysErr->Draw("E2same");

//                 etapi0RatioPbPb2760GeV_2050->Draw("c,histo,same");
                etapi0RatioFromMtScalingCombPbPb2760GeV_2050->Draw("l,same");

                TLegend* legendChargedRatio3SC = new TLegend(0.12,0.63,0.48,0.88);
                legendChargedRatio3SC->SetFillColor(0);
                legendChargedRatio3SC->SetLineColor(0);
                legendChargedRatio3SC->SetTextFont(42);
                legendChargedRatio3SC->SetTextSize(0.037);
                legendChargedRatio3SC->SetMargin(0.17);
                legendChargedRatio3SC->SetHeader(collisionSystem2760GeV.Data());
                legendChargedRatio3SC->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("%s #eta/#pi^{0}",cent2050.Data()),"fp");
                legendChargedRatio3SC->AddEntry(etapi0RatioFromMtScalingCombPbPb2760GeV_2050,"#eta from #it{m}_{T} scaled #pi^{0}","l");
                legendChargedRatio3SC->AddEntry(graphChargedRatioKaonToPionSys2040,Form("%s #eta/#pi^{0}",cent2040.Data()),"fp");
                legendChargedRatio3SC->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
                legendChargedRatio3SC->Draw();

                legendEtatoPi0combo_withPP2760GeVmt->Draw();

                if(thesisPlotting)          thesisLabelEtaPi0Ratio->Draw();

                graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
                graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

            canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_withPP2760GeVandKaonsToPionsMt_2050.%s",outputDir.Data(),suffix.Data()));
//             canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_withPP2760GeVandKaonsToPionsMt_2050.%s",paperPlots.Data(),suffix.Data()));

            canvasEtatoPi0combo->cd();
                    histo2DEtatoPi0combo->Draw("copy");


                    DrawGammaSetMarkerTGraphAsym(graphEtaToPi0JetQuenching_0010, 0, 0, colorNLO0010,colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
                    graphEtaToPi0JetQuenching_0010->Draw("3,same");

                    graphEPOSEtaToPi0_0010->GetXaxis()->SetRangeUser(0,13);
                    graphEPOSEtaToPi0_0010->Draw("same,l,x");

                    DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaToPi0LowPt_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE,colorCracow0010);
                    TheoryCracowEtaToPi0LowPt_0010->Draw("l,same");
//                     DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaToPi0LowPt_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE,colorCracow0010);
                    TheoryBegunEQEtaToPi0_0010->Draw("l,same");

                    graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
                    graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

        //          graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
        //          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_2050);
        //          graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

        //           labelEtaToPi0Energy->Draw();

                    TLegend* legendEtatoPi0combo_onlyPbPb3;
                    if(thesisPlotting)legendEtatoPi0combo_onlyPbPb3 = new TLegend(0.13,0.78,0.53,0.9);
                    else legendEtatoPi0combo_onlyPbPb3 = new TLegend(0.13,0.83,0.53,0.93);
                    legendEtatoPi0combo_onlyPbPb3->SetFillColor(0);
                    legendEtatoPi0combo_onlyPbPb3->SetLineColor(0);
                    legendEtatoPi0combo_onlyPbPb3->SetTextFont(42);
                    legendEtatoPi0combo_onlyPbPb3->SetTextSize(0.037);
                    legendEtatoPi0combo_onlyPbPb3->SetMargin(0.17);
                    legendEtatoPi0combo_onlyPbPb3->SetHeader(collisionSystem2760GeV.Data());
                    legendEtatoPi0combo_onlyPbPb3->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,Form("  %s",cent0010.Data()),"fp");
                    legendEtatoPi0combo_onlyPbPb3->AddEntry(etapi0RatioPbPb2760GeV_0010,"#eta from #it{m}_{T} scaled #pi^{0}","l");
                    legendEtatoPi0combo_onlyPbPb3->Draw();

                    TLegend* legendRatioALICE3 = new TLegend(0.58,0.83,0.94,0.93); //0.12,0.62,0.5,0.71);
                    legendRatioALICE3->SetFillColor(0);
                    legendRatioALICE3->SetLineColor(0);
                    legendRatioALICE3->SetTextFont(42);
                    legendRatioALICE3->SetTextSize(0.037);
                    legendRatioALICE3->SetMargin(0.17);
            //           legendRatioALICE2->SetHeader("NLO DCZW (#tau_{0} = 0.6 fm)");
                    legendRatioALICE3->AddEntry(graphEtaToPi0JetQuenching_0010,"0#font[122]{-}10% NLO DCZW","f");
                    legendRatioALICE3->AddEntry((TObject*)0,"PLB 750 (2015)" /*"arXiv:1506.00838"*/,"");
                    legendRatioALICE3->Draw();

                    TLegend* legendRatioALICE4A = new TLegend(0.54,0.23,0.76,0.33);
                    legendRatioALICE4A->SetFillColor(0);
                    legendRatioALICE4A->SetLineColor(0);
                    legendRatioALICE4A->SetTextFont(42);
                    legendRatioALICE4A->SetTextSize(0.037);
                    legendRatioALICE4A->SetMargin(0.28);
            //           legendRatioALICE3A->SetHeader("EPOS");
                    legendRatioALICE4A->AddEntry(graphEPOSEtaToPi0_0010,"0#font[122]{-}10% EPOS", "l");
                    legendRatioALICE4A->AddEntry((TObject*)0,"PRC 85, 064907 (2012)","");
                    legendRatioALICE4A->Draw();
                    TLegend* legendRatioALICE4B = new TLegend(0.54,0.12,0.92,0.22);
                    legendRatioALICE4B->SetFillColor(0);
                    legendRatioALICE4B->SetLineColor(0);
                    legendRatioALICE4B->SetTextFont(42);
                    legendRatioALICE4B->SetTextSize(0.037);
                    legendRatioALICE4B->SetMargin(0.17);
            //           legendRatioALICE3B->SetHeader("NEQ SHM");
                    legendRatioALICE4B->AddEntry(TheoryCracowEtaToPi0LowPt_0010,"0#font[122]{-}10% NEQ SHM","l");
                    legendRatioALICE4B->AddEntry((TObject*)0,"PRC 90, 014906 (2014)","");
                    legendRatioALICE4B->Draw();

                    etapi0RatioPbPb2760GeV_0010->Draw("c,histo,same");

        //             if(thesisPlotting)          thesisLabelHighLeft2->Draw();

                canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsAndMt_0010.%s",outputDir.Data(),suffix.Data()));
//                 canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsAndMt_0010.%s",paperPlots.Data(),suffix.Data()));
                canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsAndMt_0010.%s",PubNotePlots.Data(),suffix.Data()));

        }

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_0010);
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_2050);
            graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

            graphRatioEtaToPi0Comb2760GeVStatErr->Draw("p,same");
            graphRatioEtaToPi0Comb2760GeVSysErr->Draw("E2same");

//          TLegend* legendEtatoPi0combo_onlyPbPbwithPP = new TLegend(0.12,0.73,0.53,0.92);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetFillColor(0);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetLineColor(0);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetTextFont(42);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetTextSize(0.04);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetMargin(0.17);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetHeader(collisionSystem2760GeV.Data());
//          legendEtatoPi0combo_onlyPbPbwithPP->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,Form("  %s",cent0010.Data()),"fp");
//          legendEtatoPi0combo_onlyPbPbwithPP->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("%s",cent2050.Data()),"fp");
//          legendEtatoPi0combo_onlyPbPbwithPP->Draw();
            legendEtatoPi0combo_onlyPbPb->Draw();

            if(thesisPlotting)          thesisLabelEtaPi0Ratio->Draw();

            TLegend* legendEtatoPi0combo_withPP276GeV = new TLegend(0.55,0.15,0.95,0.26);
            legendEtatoPi0combo_withPP276GeV->SetFillColor(0);
            legendEtatoPi0combo_withPP276GeV->SetLineColor(0);
            legendEtatoPi0combo_withPP276GeV->SetTextFont(42);
            legendEtatoPi0combo_withPP276GeV->SetTextSize(0.037);
            legendEtatoPi0combo_withPP276GeV->SetMargin(0.17);
            legendEtatoPi0combo_withPP276GeV->SetHeader(collisionSystemPP2760GeV.Data());
            legendEtatoPi0combo_withPP276GeV->AddEntry(graphRatioEtaToPi0Comb2760GeVSysErr,"EPJC 77 (2017) 339"/*"arXiv: 1702.00917"*/,"fp");
            legendEtatoPi0combo_withPP276GeV->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_combined_DataOnlyWithPP276GeV.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_combined_DataOnlyWithPP276GeV.%s",PubNotePlots.Data(),suffix.Data()));


        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_0010);
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0StatPbPb2760GeV_2050);
            graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

            graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
            graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");

//          TLegend* legendEtatoPi0combo_onlyPbPbwithPP = new TLegend(0.12,0.73,0.53,0.92);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetFillColor(0);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetLineColor(0);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetTextFont(42);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetTextSize(0.04);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetMargin(0.17);
//          legendEtatoPi0combo_onlyPbPbwithPP->SetHeader(collisionSystem2760GeV.Data());
//          legendEtatoPi0combo_onlyPbPbwithPP->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,Form("  %s",cent0010.Data()),"fp");
//          legendEtatoPi0combo_onlyPbPbwithPP->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("%s",cent2050.Data()),"fp");
//          legendEtatoPi0combo_onlyPbPbwithPP->Draw();
            legendEtatoPi0combo_onlyPbPb->Draw();

            TLegend* legendEtatoPi0combo_withPP = new TLegend(0.55,0.15,0.95,0.26);
            legendEtatoPi0combo_withPP->SetFillColor(0);
            legendEtatoPi0combo_withPP->SetLineColor(0);
            legendEtatoPi0combo_withPP->SetTextFont(42);
            legendEtatoPi0combo_withPP->SetTextSize(0.037);
            legendEtatoPi0combo_withPP->SetMargin(0.17);
            legendEtatoPi0combo_withPP->SetHeader(collisionSystemPP7TeV.Data());
            legendEtatoPi0combo_withPP->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,"PLB 717 (2012) 162","fp");
            legendEtatoPi0combo_withPP->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_combined_DataOnlyWithPP7TeV.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_combined_DataOnlyWithPP7TeV.%s",PubNotePlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphEtaToPi0JetQuenching_0010, 0, 0, colorNLO0010,colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
            graphEtaToPi0JetQuenching_0010->Draw("3,same");

            graphEPOSEtaToPi0_0010->GetXaxis()->SetRangeUser(0,13);
            graphEPOSEtaToPi0_0010->Draw("same,l,x");

            DrawGammaSetMarkerTGraphAsym(TheoryBegunEQEtaToPi0_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE,colorCracow0010);
            TheoryBegunEQEtaToPi0_0010->SetLineStyle(9);
            TheoryBegunEQEtaToPi0_0010->Draw("l,same");

            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            TLegend* legendEtatoPi0combo_onlyPbPb2;
//             if(thesisPlotting)legendEtatoPi0combo_onlyPbPb2 = new TLegend(0.13,0.8,0.53,0.9);
//             else
                legendEtatoPi0combo_onlyPbPb2 = new TLegend(0.13,0.83,0.53,0.93);
            legendEtatoPi0combo_onlyPbPb2->SetFillColor(0);
            legendEtatoPi0combo_onlyPbPb2->SetLineColor(0);
            legendEtatoPi0combo_onlyPbPb2->SetTextFont(42);
            legendEtatoPi0combo_onlyPbPb2->SetTextSize(0.037);
            legendEtatoPi0combo_onlyPbPb2->SetMargin(0.17);
            legendEtatoPi0combo_onlyPbPb2->SetHeader(collisionSystem2760GeV.Data());
            legendEtatoPi0combo_onlyPbPb2->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,Form("  %s",cent0010.Data()),"fp");
            legendEtatoPi0combo_onlyPbPb2->Draw();

            TLegend* legendRatioALICE2 = new TLegend(0.58,0.83,0.94,0.93); //0.12,0.62,0.5,0.71);
            legendRatioALICE2->SetFillColor(0);
            legendRatioALICE2->SetLineColor(0);
            legendRatioALICE2->SetTextFont(42);
            legendRatioALICE2->SetTextSize(0.037);
            legendRatioALICE2->SetMargin(0.17);
    //           legendRatioALICE2->SetHeader("NLO DCZW (#tau_{0} = 0.6 fm)");
            legendRatioALICE2->AddEntry(graphEtaToPi0JetQuenching_0010,"0#font[122]{-}10% NLO DCZW","f");
            legendRatioALICE2->AddEntry((TObject*)0,"PLB 750 (2015)" /*"arXiv:1506.00838"*/,"");
            legendRatioALICE2->Draw();

            TLegend* legendRatioALICE3A = new TLegend(0.54,0.23,0.76,0.33);
            legendRatioALICE3A->SetFillColor(0);
            legendRatioALICE3A->SetLineColor(0);
            legendRatioALICE3A->SetTextFont(42);
            legendRatioALICE3A->SetTextSize(0.037);
            legendRatioALICE3A->SetMargin(0.28);
    //           legendRatioALICE3A->SetHeader("EPOS");
            legendRatioALICE3A->AddEntry(graphEPOSEtaToPi0_0010,"0#font[122]{-}10% EPOS", "l");
            legendRatioALICE3A->AddEntry((TObject*)0,"PRC 85, 064907 (2012)","");
            legendRatioALICE3A->Draw();

            TLegend* legendRatioALICE3BEQ = new TLegend(0.54,0.12,0.92,0.22);
            legendRatioALICE3BEQ->SetFillColor(0);
            legendRatioALICE3BEQ->SetLineColor(0);
            legendRatioALICE3BEQ->SetTextFont(42);
            legendRatioALICE3BEQ->SetTextSize(0.037);
            legendRatioALICE3BEQ->SetMargin(0.17);
    //           legendRatioALICE3BEQ->SetHeader("NEQ SHM");
//             legendRatioALICE3BEQ->AddEntry(TheoryCracowEtaToPi0LowPt_0010,"0#font[122]{-}10% NEQ SHM","l");
            legendRatioALICE3BEQ->AddEntry(TheoryBegunEQEtaToPi0_0010,"0#font[122]{-}10% EQ SHM","l");
            legendRatioALICE3BEQ->AddEntry((TObject*)0,"PRC 90, 014906 (2014)","");
            legendRatioALICE3BEQ->Draw();

            if(thesisPlotting)          thesisLabelHighLeft2->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsEQ_0010.%s",outputDir.Data(),suffix.Data()));
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsEQ_0010.%s",paperPlots.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsEQ_0010.%s",PubNotePlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphEtaToPi0JetQuenching_0010, 0, 0, colorNLO0010,colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
            graphEtaToPi0JetQuenching_0010->Draw("3,same");

            graphEPOSEtaToPi0_0010->GetXaxis()->SetRangeUser(0,13);
            graphEPOSEtaToPi0_0010->Draw("same,l,x");

            DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaToPi0LowPt_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE,colorCracow0010);
            TheoryCracowEtaToPi0LowPt_0010->Draw("l,same");

            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            legendEtatoPi0combo_onlyPbPb2->Draw();
            legendRatioALICE2->Draw();

            legendRatioALICE3A->Draw();

            TLegend* legendRatioALICE3BNEQ = new TLegend(0.54,0.12,0.92,0.22);
            legendRatioALICE3BNEQ->SetFillColor(0);
            legendRatioALICE3BNEQ->SetLineColor(0);
            legendRatioALICE3BNEQ->SetTextFont(42);
            legendRatioALICE3BNEQ->SetTextSize(0.037);
            legendRatioALICE3BNEQ->SetMargin(0.17);
            legendRatioALICE3BNEQ->AddEntry(TheoryCracowEtaToPi0LowPt_0010,"0#font[122]{-}10% NEQ SHM","l");
            legendRatioALICE3BNEQ->AddEntry((TObject*)0,"PRC 90, 014906 (2014)","");
            legendRatioALICE3BNEQ->Draw();

            if(thesisPlotting)          thesisLabelHighLeft2->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsNEQ_0010.%s",outputDir.Data(),suffix.Data()));
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsNEQ_0010.%s",paperPlots.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsNEQ_0010.%s",PubNotePlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphEtaToPi0JetQuenching_0010, 0, 0, colorNLO0010,colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
            graphEtaToPi0JetQuenching_0010->Draw("3,same");

            while(graphEPOSEtaToPi0_0010->GetX()[graphEPOSEtaToPi0_0010->GetN()-1]>17.)graphEPOSEtaToPi0_0010->RemovePoint(graphEPOSEtaToPi0_0010->GetN()-1);
            graphEPOSEtaToPi0_0010->Draw("same,l,x");

            DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaToPi0LowPt_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE,colorCracow0010);
            TheoryCracowEtaToPi0LowPt_0010->Draw("l,same");
            DrawGammaSetMarkerTGraphAsym(TheoryBegunEQEtaToPi0_0010, 0, 0, colorBegun0010,colorBegun0010, 5, kTRUE,colorBegun0010);
            TheoryBegunEQEtaToPi0_0010->SetLineStyle(9);
            TheoryBegunEQEtaToPi0_0010->Draw("l,same");


            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            legendEtatoPi0combo_onlyPbPb2->Draw();
            legendRatioALICE2->Draw();
            legendRatioALICE3A->Draw();

            TLegend* legendRatioSHM = new TLegend(0.5,0.12,0.95,0.22);
            legendRatioSHM->SetFillColor(0);
            legendRatioSHM->SetLineColor(0);
            legendRatioSHM->SetTextFont(42);
            legendRatioSHM->SetTextSize(0.037);
            legendRatioSHM->SetMargin(0.17);
            legendRatioSHM->SetNColumns(2);
            legendRatioSHM->SetHeader("SHM - PRC 90, 014906 (2014)");
            legendRatioSHM->AddEntry(TheoryCracowEtaToPi0LowPt_0010,"0#font[122]{-}10% NEQ","l");
            legendRatioSHM->AddEntry(TheoryBegunEQEtaToPi0_0010,"0#font[122]{-}10% EQ","l");
            legendRatioSHM->Draw();

//             if(thesisPlotting)          thesisLabelHighLeft2->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModels_0010.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModels_0010.%s",paperPlots.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModels_0010.%s",PubNotePlots.Data(),suffix.Data()));


        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphErr(graphEPOSEtaToPi0_2050,0,0,colorEPOS2050,colorEPOS2050,2,kTRUE,colorEPOS2050);
            graphEPOSEtaToPi0_2050->SetLineStyle(5);
            graphEPOSEtaToPi0_2050->GetXaxis()->SetRangeUser(0,13);
            graphEPOSEtaToPi0_2050->Draw("same,l,x");

            DrawGammaSetMarkerTGraphAsym(TheoryBegunEQEtaToPi0_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
            TheoryBegunEQEtaToPi0_2050->SetLineStyle(9);
            TheoryBegunEQEtaToPi0_2050->Draw("l,same");

            graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

            TLegend* legendEtatoPi0combo_onlyPbPb3;
            if(thesisPlotting)legendEtatoPi0combo_onlyPbPb3 = new TLegend(0.13,0.8,0.53,0.9);
            else legendEtatoPi0combo_onlyPbPb3 = new TLegend(0.13,0.83,0.53,0.93);
            legendEtatoPi0combo_onlyPbPb3->SetFillColor(0);
            legendEtatoPi0combo_onlyPbPb3->SetLineColor(0);
            legendEtatoPi0combo_onlyPbPb3->SetTextFont(42);
            legendEtatoPi0combo_onlyPbPb3->SetTextSize(0.037);
            legendEtatoPi0combo_onlyPbPb3->SetMargin(0.17);
            legendEtatoPi0combo_onlyPbPb3->SetHeader(collisionSystem2760GeV.Data());
            legendEtatoPi0combo_onlyPbPb3->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("%s ",cent2050.Data()),"fp");
            legendEtatoPi0combo_onlyPbPb3->Draw();

            TLegend* legendRatioEtaPi0RatioEpos = new TLegend(0.54,0.23,0.76,0.33);
            legendRatioEtaPi0RatioEpos->SetFillColor(0);
            legendRatioEtaPi0RatioEpos->SetLineColor(0);
            legendRatioEtaPi0RatioEpos->SetTextFont(42);
            legendRatioEtaPi0RatioEpos->SetTextSize(0.037);
            legendRatioEtaPi0RatioEpos->SetMargin(0.28);
    //           legendRatioEtaPi0RatioEpos->SetHeader("EPOS");
            legendRatioEtaPi0RatioEpos->AddEntry(graphEPOSEtaToPi0_2050,"20#font[122]{-}50% EPOS", "l");
            legendRatioEtaPi0RatioEpos->AddEntry((TObject*)0,"PRC 85, 064907 (2012)","");
            legendRatioEtaPi0RatioEpos->Draw();
            TLegend* legendEtaPi0RatioEQ = new TLegend(0.54,0.12,0.92,0.22);
            legendEtaPi0RatioEQ->SetFillColor(0);
            legendEtaPi0RatioEQ->SetLineColor(0);
            legendEtaPi0RatioEQ->SetTextFont(42);
            legendEtaPi0RatioEQ->SetTextSize(0.037);
            legendEtaPi0RatioEQ->SetMargin(0.17);
            legendEtaPi0RatioEQ->AddEntry(TheoryBegunEQEtaToPi0_2050,"20#font[122]{-}50% EQ SHM","l");
            legendEtaPi0RatioEQ->AddEntry((TObject*)0,"PRC 90, 014906 (2014)","");
            legendEtaPi0RatioEQ->Draw();

            if(thesisPlotting)          thesisLabelHighLeft2->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsEQ_2050.%s",outputDir.Data(),suffix.Data()));
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsEQ_2050.%s",PubNotePlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphErr(graphEPOSEtaToPi0_2050,0,0,colorEPOS2050,colorEPOS2050,2,kTRUE,colorEPOS2050);
            graphEPOSEtaToPi0_2050->SetLineStyle(5);
            graphEPOSEtaToPi0_2050->GetXaxis()->SetRangeUser(0,13);
            graphEPOSEtaToPi0_2050->Draw("same,l,x");

            DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaToPi0LowPt_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
            TheoryCracowEtaToPi0LowPt_2050->Draw("l,same");

            graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

            legendEtatoPi0combo_onlyPbPb3->Draw();
            legendRatioEtaPi0RatioEpos->Draw();
            TLegend* legendEtaPi0RatioNEQ = new TLegend(0.54,0.12,0.92,0.22);
            legendEtaPi0RatioNEQ->SetFillColor(0);
            legendEtaPi0RatioNEQ->SetLineColor(0);
            legendEtaPi0RatioNEQ->SetTextFont(42);
            legendEtaPi0RatioNEQ->SetTextSize(0.037);
            legendEtaPi0RatioNEQ->SetMargin(0.17);
            legendEtaPi0RatioNEQ->AddEntry(TheoryCracowEtaToPi0LowPt_2050,"20#font[122]{-}50% NEQ SHM","l");
            legendEtaPi0RatioNEQ->AddEntry((TObject*)0,"PRC 90, 014906 (2014)","");
            legendEtaPi0RatioNEQ->Draw();

            if(thesisPlotting)          thesisLabelHighLeft2->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsNEQ_2050.%s",outputDir.Data(),suffix.Data()));
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModelsNEQ_2050.%s",PubNotePlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphErr(graphEPOSEtaToPi0_2050,0,0,colorEPOS2050,colorEPOS2050,2,kTRUE,colorEPOS2050);
            graphEPOSEtaToPi0_2050->SetLineStyle(5);
            graphEPOSEtaToPi0_2050->GetXaxis()->SetRangeUser(0,13);
            graphEPOSEtaToPi0_2050->Draw("same,l,x");

            DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaToPi0LowPt_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
            TheoryCracowEtaToPi0LowPt_2050->Draw("l,same");
            DrawGammaSetMarkerTGraphAsym(TheoryBegunEQEtaToPi0_2050, 0, 0, colorBegun2050,colorBegun2050, 5, kTRUE,colorBegun2050);
            TheoryBegunEQEtaToPi0_2050->SetLineStyle(9);
            TheoryBegunEQEtaToPi0_2050->Draw("l,same");

            graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

            legendEtatoPi0combo_onlyPbPb3->Draw();
            legendRatioEtaPi0RatioEpos->Draw();
            TLegend* legendEtaPi0RatioSHM = new TLegend(0.5,0.12,0.95,0.22);
            legendEtaPi0RatioSHM->SetFillColor(0);
            legendEtaPi0RatioSHM->SetLineColor(0);
            legendEtaPi0RatioSHM->SetTextFont(42);
            legendEtaPi0RatioSHM->SetTextSize(0.037);
            legendEtaPi0RatioSHM->SetMargin(0.17);
            legendEtaPi0RatioSHM->SetNColumns(2);
            legendEtaPi0RatioSHM->SetHeader("SHM - PRC 90, 014906 (2014)");
            legendEtaPi0RatioSHM->AddEntry(TheoryCracowEtaToPi0LowPt_2050,"20#font[122]{-}50% NEQ","l");
            legendEtaPi0RatioSHM->AddEntry(TheoryBegunEQEtaToPi0_2050,"20#font[122]{-}50% EQ","l");
            legendEtaPi0RatioSHM->Draw();

            if(thesisPlotting)          thesisLabelHighLeft2->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModels_2050.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_WithModels_2050.%s",PubNotePlots.Data(),suffix.Data()));


        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");
            histo2DEtatoPi0combo->GetYaxis()->SetTitleOffset(1.1);
            histo2DEtatoPi0combo->GetYaxis()->SetTitle("Particle ratio");

            graphChargedRatioKaonToPionSys0010->Draw("2same");
            graphChargedRatioKaonToPion0010->Draw("p,same");

            TLegend* legendChargedRatio = new TLegend(0.12,0.73,0.95,0.89);
            legendChargedRatio->SetFillColor(0);
            legendChargedRatio->SetLineColor(0);
            legendChargedRatio->SetTextFont(42);
            legendChargedRatio->SetTextSize(0.037);
            legendChargedRatio->SetMargin(0.17);
            legendChargedRatio->SetHeader(collisionSystem2760GeV.Data());
            legendChargedRatio->SetNColumns(2);
            legendChargedRatio->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,Form("#eta/#pi^{0},   %s",cent0010.Data()),"fp");
            legendChargedRatio->AddEntry(graphChargedRatioKaonToPionSys0010,Form("K^{#pm}/#pi^{#pm}, %s",cent0010.Data()),"fp");
            legendChargedRatio->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("#eta/#pi^{0}, %s",cent2050.Data()),"fp");
            legendChargedRatio->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
            legendChargedRatio->Draw();

            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");
                        if(thesisPlotting)          thesisLabelHighLeft2->Draw();


        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_withKaonsToPions.%s",outputDir.Data(),suffix.Data()));
    //       canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_withKaonsToPions.%s",paperPlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");
            histo2DEtatoPi0combo->GetYaxis()->SetTitleOffset(1.1);
            histo2DEtatoPi0combo->GetYaxis()->SetTitle("Particle ratio");

            DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPion0010, 24,markerSizeComb, colorCharged, colorCharged);
            DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPionSys0010,24,markerSizeComb, colorCharged, colorCharged, 1, kTRUE);
            graphChargedRatioKaonToPionSys0010->Draw("2same");
            graphChargedRatioKaonToPion0010->Draw("p,same");

            TLegend* legendChargedRatio0010 = new TLegend(0.12,0.72,0.5,0.9);
            legendChargedRatio0010->SetFillColor(0);
            legendChargedRatio0010->SetLineColor(0);
            legendChargedRatio0010->SetTextFont(42);
            legendChargedRatio0010->SetTextSize(0.037);
            legendChargedRatio0010->SetMargin(0.17);
            legendChargedRatio0010->SetHeader(collisionSystem2760GeV.Data());
    //           legendChargedRatio0010->SetNColumns(2);
            legendChargedRatio0010->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,Form("#eta/#pi^{0},   %s",cent0010.Data()),"fp");
            legendChargedRatio0010->AddEntry(graphChargedRatioKaonToPionSys0010,Form("K^{#pm}/#pi^{#pm}, %s",cent0010.Data()),"fp");
            legendChargedRatio0010->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
            legendChargedRatio0010->Draw();

            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

            if(thesisPlotting)          thesisLabelHighLeft2->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_withKaonsToPions_0010.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_withKaonsToPions_0010.%s",PubNotePlots.Data(),suffix.Data()));

        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");
            histo2DEtatoPi0combo->GetYaxis()->SetTitleOffset(1.1);
            histo2DEtatoPi0combo->GetYaxis()->SetTitle("Particle ratio");

            DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPion2040, 25,markerSizeComb, colorCharged, colorCharged);
            DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPionSys2040,25,markerSizeComb, colorCharged, colorCharged, 1, kTRUE);
            graphChargedRatioKaonToPionSys2040->Draw("2same");
            graphChargedRatioKaonToPion2040->Draw("p,same");

            TLegend* legendChargedRatio2050 = new TLegend(0.12,0.72,0.5,0.9);
            legendChargedRatio2050->SetFillColor(0);
            legendChargedRatio2050->SetLineColor(0);
            legendChargedRatio2050->SetTextFont(42);
            legendChargedRatio2050->SetTextSize(0.037);
            legendChargedRatio2050->SetMargin(0.17);
            legendChargedRatio2050->SetHeader(collisionSystem2760GeV.Data());
    //           legendChargedRatio2050->SetNColumns(2);
            legendChargedRatio2050->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_2050,Form("#eta/#pi^{0},   %s",cent2050.Data()),"fp");
            legendChargedRatio2050->AddEntry(graphChargedRatioKaonToPionSys2040,Form("K^{#pm}/#pi^{#pm}, %s",cent2040.Data()),"fp");
            legendChargedRatio2050->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
            legendChargedRatio2050->Draw();

            if(thesisPlotting)          thesisLabelHighLeft2->Draw();
            graphCombEtatoPi0SysPbPb2760GeV_2050->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_2050->Draw("p,same");

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_withKaonsToPions_2050.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_withKaonsToPions_2050.%s",PubNotePlots.Data(),suffix.Data()));


        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");
            histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.,1.15);

            graphChargedRatioKaonToPionSys0010->Draw("2same");
            graphChargedRatioKaonToPion0010->Draw("p,same");

            graphRatioEtaToPi0Comb2760GeVStatErr->Draw("p,same");
            graphRatioEtaToPi0Comb2760GeVSysErr->Draw("E2same");

            TLegend* legendChargedRatio2 = new TLegend(0.12,0.72,0.48,0.92);//0.12,0.7,0.48,0.9);
            legendChargedRatio2->SetFillColor(0);
            legendChargedRatio2->SetLineColor(0);
            legendChargedRatio2->SetTextFont(42);
            legendChargedRatio2->SetTextSize(0.037);
            legendChargedRatio2->SetMargin(0.17);
            legendChargedRatio2->SetHeader(collisionSystemPbPb0010.Data());
            legendChargedRatio2->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,"#eta/#pi^{0}","fp");
    //             legendChargedRatio2->AddEntry((TObject*)0,"","");
            legendChargedRatio2->AddEntry(graphChargedRatioKaonToPionSys0010,"K^{#pm}/#pi^{#pm}","fp");
            legendChargedRatio2->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
            legendChargedRatio2->Draw();

            TLegend* legendEtatoPi0combo_withPP2760GeV = GetAndSetLegend(0.6,0.14,3.);// new TLegend(0.55,0.15,0.95,0.3/*26*/);
            legendEtatoPi0combo_withPP2760GeV->SetFillColor(0);
            legendEtatoPi0combo_withPP2760GeV->SetLineColor(0);
            legendEtatoPi0combo_withPP2760GeV->SetTextFont(42);
            legendEtatoPi0combo_withPP2760GeV->SetTextSize(0.037);
            legendEtatoPi0combo_withPP2760GeV->SetMargin(0.17);
            legendEtatoPi0combo_withPP2760GeV->SetHeader(collisionSystemPP2760GeV.Data());
            legendEtatoPi0combo_withPP2760GeV->AddEntry(graphRatioEtaToPi0Comb2760GeVSysErr,"#eta/#pi^{0}"/*"arXiv: 1702.00917"*/,"fp");
            legendEtatoPi0combo_withPP2760GeV->AddEntry((TObject*)0,"EPJC 77 (2017) 339","");
    //           legendEtatoPi0combo_withPP2760GeV->AddEntry(etapi0Ratio2760GeV,"#eta from #it{m}_{T} scaled #pi^{0}","l");
            legendEtatoPi0combo_withPP2760GeV->Draw();

//             if(thesisPlotting)          thesisLabelHighLeft2->Draw();

            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_withPP2760GeVandKaonsToPions.%s",outputDir.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_withPP2760GeVandKaonsToPions.%s",paperPlots.Data(),suffix.Data()));


        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");
            histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.,1.15);

            graphChargedRatioKaonToPionSys0010->Draw("2same");
            graphChargedRatioKaonToPion0010->Draw("p,same");

            graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
            graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");

            legendChargedRatio2->Draw();

            legendEtatoPi0combo_withPP->Draw();

            graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
            graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_withPP7TeVandKaonsToPions.%s",outputDir.Data(),suffix.Data()));
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_withPP7TeVandKaonsToPions.%s",paperPlots.Data(),suffix.Data()));

    }


    //**********************************************************************************************************************//
    //*************************************************   Theory plots    **************************************************//
    //**********************************************************************************************************************//
    // theory only plotting

//     // Cracow Low pt model ---- NEQ SHM
//     TCanvas* canvasTheoryCracow = new TCanvas("canvasTheoryCracow","",200,10,1350,1350*1.15);  // gives the page size
//     DrawGammaCanvasSettings( canvasTheoryCracow, 0.16, 0.02, 0.02, 0.09);
//     canvasTheoryCracow->SetLogx();
//     canvasTheoryCracow->SetLogy();
//
//       TH2F * histo2DCracowYields = new TH2F("histo2DCracowYields","histo2DCracowYields",11000,0.09,5.,1000,1e-2,1e4);
//       SetStyleHistoTH2ForGraphs(histo2DCracowYields, "#it{p}_{T} (GeV/#it{c})","",0.035,0.04, 0.035,0.04, 1.,1.65);
//       histo2DCracowYields->DrawCopy();
//
//       TheoryCracowPi0LowPt_0010->Draw("l,same");
//       TheoryCracowEtaLowPt_0010->Draw("l,same");
//
//       TheoryCracowChargedPionLowPt_0010->Draw("l,same");
//       TheoryCracowChargedKaonLowPt_0010->Draw("l,same");
//
//       TLegend* legYield = new TLegend(0.45,0.85,0.95,0.95);
//       legYield->SetFillColor(0);
//       legYield->SetLineColor(0);
//       legYield->SetTextFont(42);
//       legYield->SetTextSize(FontSize);
//       legYield->SetNColumns(2);
//       legYield->SetHeader("NEQ SHM - PRC 90, 014906 (2014)");
//       legYield->AddEntry(TheoryCracowChargedPionLowPt_0010,"#pi^{#pm}","l");
//       legYield->AddEntry(TheoryCracowChargedKaonLowPt_0010,"K^{#pm}","l");
//       legYield->AddEntry(TheoryCracowPi0LowPt_0010,"#pi^{0}","l");
//       legYield->AddEntry(TheoryCracowEtaLowPt_0010,"#eta","l");
//       legYield->Draw();
//
//     canvasTheoryCracow->SaveAs(Form("%s/CracowTheoryWithCharged.%s",outputDir.Data(),suffix.Data()));
//
//     TCanvas* canvasTheoryCracowRatio = new TCanvas("canvasTheoryCracowRatio","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
//     DrawGammaCanvasSettings( canvasTheoryCracowRatio, 0.09, 0.03, 0.035, 0.1);
//     canvasTheoryCracowRatio->SetLogx();
//
//       TH2F * histo2DTheoryCracowRatio = new TH2F("histo2DTheoryCracowRatio","histo2DTheoryCracowRatio",11000,0.09,10.,1000,0.01,2.);
//       SetStyleHistoTH2ForGraphs(histo2DTheoryCracowRatio, "#it{p}_{T} (GeV/#it{c})","Particle ratio",0.035,0.04, 0.035,0.04, 1.2,1.);
//       histo2DTheoryCracowRatio->Draw("copy");
//       histo2DTheoryCracowRatio->GetYaxis()->SetRangeUser(0.,1.);
//       histo2DTheoryCracowRatio->GetXaxis()->SetRangeUser(0.,7.);
//
//
//         if(meson.CompareTo("Eta")==0){
//           graphChargedRatioKaonToPionSys0010->Draw("2same");
//           graphChargedRatioKaonToPion0010->Draw("p,same");
//           graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
//           graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");
//         }
//
//         DrawGammaSetMarkerTGraphAsym(TheoryCracowChargedKaonToPionLowPt_0010, 0, 0, kOrange+1,kOrange+1, 5, kTRUE, kOrange+1);
//         TheoryCracowChargedKaonToPionLowPt_0010->Draw("l,same");
//
//         DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaToPi0_0010,0, 0, kRed+1,kRed+1, 5, kTRUE, kRed+1);
//         TheoryCracowEtaToPi0_0010->Draw("l,same");
//
//         DrawGammaSetMarkerTGraphAsym(TheoryCracowChargedEtaToKaonLowPt_0010, 0, 0, kGreen+1,kGreen+1, 5, kTRUE,kGreen+1);
//         TheoryCracowChargedEtaToKaonLowPt_0010->Draw("l,same");
//
//         TLegend* legRatio = new TLegend(0.15,0.6,0.5,0.93);
//         legRatio->SetFillColor(0);
//         legRatio->SetLineColor(0);
//         legRatio->SetTextFont(42);
//         legRatio->SetTextSize(FontSize);
//         legRatio->SetMargin(0.17);
//         legRatio->SetHeader("NEQ SHM - PRC 90, 014906 (2014)");
//         legRatio->AddEntry(TheoryCracowChargedKaonToPionLowPt_0010,"K^{#pm}/#pi^{#pm}","l");
//         legRatio->AddEntry(TheoryCracowEtaToPi0_0010,"#eta/#pi^{0}","l");
//         legRatio->AddEntry(TheoryCracowChargedEtaToKaonLowPt_0010,"#eta/K^{#pm}","l");
//         if(meson.CompareTo("Eta")==0){
//           legRatio->AddEntry((TObject*)0,"0#font[122]{-}10% ALICE ","");
//           legRatio->AddEntry(graphChargedRatioKaonToPionSys0010,"K^{#pm}/#pi^{#pm}","fp");
//           legRatio->AddEntry(graphCombEtatoPi0SysPbPb2760GeV_0010,"ALICE #eta/#pi^{0}","fp");
//         }
//         legRatio->Draw();
//
//     canvasTheoryCracowRatio->SaveAs(Form("%s/CracowTheoryWithChargedRatio.%s",outputDir.Data(),suffix.Data()));
//     canvasTheoryCracowRatio->SaveAs(Form("%s/CracowTheoryWithChargedRatio.%s",PubNotePlots.Data(),suffix.Data()));
//
//     canvasTheoryCracowRatio->cd();
//     histo2DTheoryCracowRatio->Draw("copy");
//     histo2DTheoryCracowRatio->GetYaxis()->SetRangeUser(0.8,1.4);
//     histo2DTheoryCracowRatio->GetXaxis()->SetRangeUser(0.,5.);
//
//       DrawGammaSetMarkerTGraphAsym(TheoryCracowPi0ToChargedPionLowPt_0010,0, 0, kCyan+2,kCyan+2, 5, kTRUE, kCyan+2);
//       TheoryCracowPi0ToChargedPionLowPt_0010->Draw("l,same");
//       DrawGammaLines(0.09, 5. , 1., 1.,1.,kGray);
//
//       TLegend* legRatioPi = new TLegend(0.15,0.77,0.5,0.93);
//       legRatioPi->SetFillColor(0);
//       legRatioPi->SetLineColor(0);
//       legRatioPi->SetTextFont(42);
//       legRatioPi->SetTextSize(FontSize);
//       legRatioPi->SetHeader("NEQ SHM - PRC 90, 014906 (2014)");
//       legRatioPi->AddEntry(TheoryCracowPi0ToChargedPionLowPt_0010,"#pi^{0}/#pi^{#pm}","l");
//       legRatioPi->Draw();
//
//     canvasTheoryCracowRatio->SaveAs(Form("%s/CracowTheoryPi0ToChargedPions.%s",outputDir.Data(),suffix.Data()));
//
//
//
//     // EPOS
//     DrawGammaSetMarkerTGraph(gEPOS_Spec_0005, 20,1.2, colorEPOS0005, colorEPOS0005,2);
//     DrawGammaSetMarkerTGraph(gEPOS_Spec_0510, 20,1.2, colorEPOS0510, colorEPOS0510,2);
//     DrawGammaSetMarkerTGraphErr(gEPOS_Spec_0010, 20,1.2, colorEPOS0005, colorEPOS0005,2);
//     TGraph *gEPOS_SpecScaled_2040 = ScaleGraph(gEPOS_Spec_2040,0.1);
//     DrawGammaSetMarkerTGraph(gEPOS_SpecScaled_2040, 21,1.2, colorEPOS2040, colorEPOS2040,2);
// //       DrawGammaSetMarkerTGraph(graphEPOS3Pi0_0010, 24,1.2, kRed, kRed,2);
//
//     TCanvas* canvasTheoryEPOS = new TCanvas("canvasTheoryEPOS","",200,10,1350,1350*1.15);  // gives the page size
//     DrawGammaCanvasSettings( canvasTheoryEPOS, 0.16, 0.02, 0.02, 0.09);
//     canvasTheoryEPOS->SetLogx();
//     canvasTheoryEPOS->SetLogy();
//
//         TH2F * histo2DTheoryEPOS = new TH2F("histo2DTheoryEPOS","histo2DTheoryEPOS",11000,0.01,20.,1000,1e-7,1e4);
//         SetStyleHistoTH2ForGraphs(histo2DTheoryEPOS, "#it{p}_{T} (GeV/#it{c})","",0.035,0.04, 0.035,0.04, 1.,1.);
//         histo2DTheoryEPOS->DrawCopy();
//
//   //         DrawGammaSetMarkerTGraphAsym(graphEPOS3Pi0_0010,1,2,kRed+2,kRed+2);
//   //         graphEPOS3Pi0_0010->SetLineWidth(2);
//   //         graphEPOS3Pi0_0010->Draw("same,histo");
//
//         graphEPOSpred_0010->Draw("l,same");
//
//         TH1D *graphEPOSpredScaledforPlot_2050 = (TH1D*)graphEPOSpred_2050->Clone("graphEPOSpredScaledforPlot_2050");
//         graphEPOSpredScaledforPlot_2050->Scale(1./10);
//         graphEPOSpredScaledforPlot_2050->SetLineWidth(2);
//         graphEPOSpredScaledforPlot_2050->Draw("l,same");
//
//   //         gEPOS_Spec_0005->Draw("same,p");
//   //         gEPOS_Spec_0510->Draw("same,p");
//         gEPOS_Spec_0010->Draw("same,p");
//         gEPOS_SpecScaled_2040->Draw("same,p");
//
//   //         graphEPOS3Pi0_0010->Draw("same,p");
//
//   //         TF1* fitToEPOS_0010 = (TF1*)fitBylinkinPbPb2760GeVPtLHC11h_0010->Clone();
//   //         graphEPOSpred_0010->Fit(fitToEPOS_0010,"NRMEX0+","",0.,20.);
//   //         fitToEPOS_0010->SetLineColor(kBlue);
//   //         fitToEPOS_0010->Draw("same");
//
//         TLegend* legEPOS = new TLegend(0.2,0.15,0.35,0.35);
//         legEPOS->SetFillColor(0);
//         legEPOS->SetLineColor(0);
//         legEPOS->SetTextFont(42);
//         legEPOS->SetTextSize(FontSize);
//     //       legEPOS->AddEntry(graphEPOS3Pi0_0010,"EPOS3 new - 0#font[122]{-}10%","pl");
//         legEPOS->AddEntry(graphEPOSpred_0010,"EPOS2 new - 0#font[122]{-}10%","l");
//     //       legEPOS->AddEntry(gEPOS_Spec_0005,"EPOS old - 0-5%","lp");
//         legEPOS->AddEntry(gEPOS_Spec_0010,"EPOS (from published pi0 paper) - 0#font[122]{-}10%","pl");
//         legEPOS->AddEntry((TObject*)0,"semicentral models scaled by 10:","");
//         legEPOS->AddEntry(graphEPOSpredScaledforPlot_2050,"EPOS2 new - 20#font[122]{-}50%","l");
//         legEPOS->AddEntry(gEPOS_SpecScaled_2040,"EPOS (from published pi0 paper) - 20-40%","pl");
//         legEPOS->Draw();
//
//     canvasTheoryEPOS->SaveAs(Form("%s/Pi0_EPOSTheoryComparison.%s",outputDir.Data(),suffix.Data()));
//
//     canvasTheoryEPOS->cd();
//         histo2DTheoryEPOS->DrawCopy();
//
//         graphEPOSpred_0010->Draw("l,same");
//
//         DrawGammaSetMarker(histoUrQMDPionEPOS_0010,1,2,kGray+1,kGray+1);
//         histoUrQMDPionEPOS_0010->SetLineWidth(2);
//         DrawGammaSetMarker(histoNOUrQMDPionEPOS_0010,1,2,kBlack,kBlack);
//         histoNOUrQMDPionEPOS_0010->SetLineWidth(2);
//
//         DrawGammaSetMarker(histoUrQMDKaonEPOS_0010,1,2,kGray+1,kGray+1);
//         histoUrQMDKaonEPOS_0010->SetLineWidth(2);
//         DrawGammaSetMarker(histoNOUrQMDKaonEPOS_0010,1,2,kBlack,kBlack);
//         histoNOUrQMDKaonEPOS_0010->SetLineWidth(2);
//
// //         DrawGammaSetMarkerTGraph(graphEPOS3Eta_0010, 1,2, kRed, kRed);
//
//         if(meson.CompareTo("Pi0")==0){
// //           graphEPOS3Pi0_0010->Draw("same,l");
//           histoUrQMDPionEPOS_0010->Draw("l,same");
//           histoNOUrQMDPionEPOS_0010->Draw("l,same");
//
//         } else if(meson.CompareTo("Eta")==0){
// //           graphEPOS3Eta_0010->Draw("same,l");
//           histoUrQMDKaonEPOS_0010->Draw("l,same");
//           histoNOUrQMDKaonEPOS_0010->Draw("l,same");
//         }
//
//         TLegend* legEPOS2 = new TLegend(0.2,0.15,0.35,0.35);
//         legEPOS2->SetFillColor(0);
//         legEPOS2->SetLineColor(0);
//         legEPOS2->SetTextFont(42);
//         legEPOS2->SetTextSize(FontSize);
//         if(meson.CompareTo("Pi0")==0){
//           legEPOS2->AddEntry(histoNOUrQMDPionEPOS_0010,"EPOS #pi^{#pm} no UrQMD  - 0#font[122]{-}10%","l");
//           legEPOS2->AddEntry(histoUrQMDPionEPOS_0010,"EPOS #pi^{#pm} with UrQMD  - 0#font[122]{-}10%","l");
//   //         legEPOS2->AddEntry(graphEPOS3Pi0_0010,"EPOS3 #pi^{0} - 0#font[122]{-}10%","pl");
//           legEPOS2->AddEntry(graphEPOSpred_0010,"EPOS2 #pi^{0} - 0#font[122]{-}10%","l");
//         } else if(meson.CompareTo("Eta")==0){
//           legEPOS2->AddEntry(histoNOUrQMDKaonEPOS_0010,"EPOS K^{#pm} no UrQMD  - 0#font[122]{-}10%","l");
//           legEPOS2->AddEntry(histoUrQMDKaonEPOS_0010,"EPOS K^{#pm} with UrQMD  - 0#font[122]{-}10%","l");
//   //         legEPOS2->AddEntry(graphEPOS3Eta_0010,"EPOS3 #eta - 0#font[122]{-}10%","pl");
//           legEPOS2->AddEntry(graphEPOSpred_0010,"EPOS2 #eta - 0#font[122]{-}10%","l");
//         }
//   //       legEPOS2->AddEntry(gEPOS_Spec_0010,"EPOS (from published pi0 paper) - 0#font[122]{-}10%","pl");
//         legEPOS2->Draw();
//
//     canvasTheoryEPOS->SaveAs(Form("%s/%s_EPOSTheoryChargeComparison.%s",outputDir.Data(),meson.Data(),suffix.Data()));
//
//
//     TCanvas* canvasTheoryEPOSRatio = new TCanvas("canvasTheoryEPOSRatio","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
//     DrawGammaCanvasSettings( canvasTheoryEPOSRatio, 0.09, 0.03, 0.035, 0.1);
//     canvasTheoryEPOSRatio->SetLogx();
//
//       TH2F * histo2DTheoryEPOSRatio = new TH2F("histo2DTheoryEPOSRatio","histo2DTheoryEPOSRatio",11000,0.09,11.,1000,0.01,10.);
//       SetStyleHistoTH2ForGraphs(histo2DTheoryEPOSRatio, "#it{p}_{T} (GeV/#it{c})","EPOS particle ratios",0.035,0.04, 0.035,0.04, 1.2,1.);
//       histo2DTheoryEPOSRatio->Draw("copy");
//       histo2DTheoryEPOSRatio->GetYaxis()->SetRangeUser(0.,2.);
//
//         TH1D *histoEtaPi0ratio = (TH1D*)graphEPOSEta_0010->Clone("histoEtaPi0ratio");
//         histoEtaPi0ratio->Divide(graphEPOSEta_0010,graphEPOSPi0_0010,1.,1.,"");
//         DrawGammaSetMarker(histoEtaPi0ratio,1,2,kRed+1,kRed+1);
//         histoEtaPi0ratio->SetLineWidth(2);
//         histoEtaPi0ratio->SetLineStyle(1);
//         histoEtaPi0ratio->Draw("same,l");
//
//         TH1D *histoKaonsPionsratio = (TH1D*)histoNOUrQMDKaonEPOS_0010->Clone("histoKaonsPionsratio");
//         histoKaonsPionsratio->Divide(histoNOUrQMDKaonEPOS_0010,histoNOUrQMDPionEPOS_0010,1.,1.,"");
//         DrawGammaSetMarker(histoKaonsPionsratio,1,2,kOrange+1,kOrange+1);
//         histoKaonsPionsratio->SetLineWidth(2);
//         histoKaonsPionsratio->SetLineStyle(1);
//         histoKaonsPionsratio->Draw("same,l");
//
//         TGraphErrors* a  = NULL;
//         TGraphErrors* b  = NULL;
//         TGraphErrors* c  = NULL;
//         TGraphErrors* d  = NULL;
//         TGraphErrors* histoPi0Pionsratio = CalculateRatioBetweenSpectraWithDifferentBinning(
//                                                     graphEPOSPi0_0010, graphEPOSPi0_0010,
//                                                     histoNOUrQMDPionEPOS_0010, histoNOUrQMDPionEPOS_0010,
//                                                     kTRUE,  kTRUE, &a, &b, &c, &d);
//
//         DrawGammaSetMarkerTGraph(histoPi0Pionsratio, 24,2, kBlue, kBlue,2);
//         histoPi0Pionsratio->Draw("same,l");
//
//         a  = NULL;
//         b  = NULL;
//         c  = NULL;
//         d  = NULL;
//         TGraphErrors* histoEtaKaonratio = CalculateRatioBetweenSpectraWithDifferentBinning(
//                                                     graphEPOSEta_0010, graphEPOSEta_0010,
//                                                     histoNOUrQMDKaonEPOS_0010, histoNOUrQMDKaonEPOS_0010,
//                                                     kTRUE,  kTRUE, &a, &b, &c, &d);
//         DrawGammaSetMarkerTGraph(histoEtaKaonratio, 24,2, kGreen+1, kGreen+1,2);
//         histoEtaKaonratio->Draw("same,l");
//
//         TLegend* legEPOSRatio = new TLegend(0.2,0.77,0.35,0.93);
//         legEPOSRatio->SetFillColor(0);
//         legEPOSRatio->SetLineColor(0);
//         legEPOSRatio->SetTextFont(42);
//         legEPOSRatio->SetTextSize(FontSize);
//         legEPOSRatio->AddEntry(histoEtaPi0ratio,"#eta/#pi^{0}","l");
//         legEPOSRatio->AddEntry(histoKaonsPionsratio,"K^{#pm}/#pi^{#pm}","l");
//         legEPOSRatio->AddEntry(histoPi0Pionsratio,"#pi^{0}/#pi^{#pm}","l");
//         legEPOSRatio->AddEntry(histoEtaKaonratio,"#eta/K^{#pm}","l");
//         legEPOSRatio->Draw();
//
//     canvasTheoryEPOSRatio->SaveAs(Form("%s/EPOSTheoryChargedRatios.%s",outputDir.Data(),suffix.Data()));
//     canvasTheoryEPOSRatio->SaveAs(Form("%s/EPOSTheoryChargedRatios.%s",PubNotePlots.Data(),suffix.Data()));
//
//     canvasTheoryEPOSRatio->cd();
//       histo2DTheoryEPOSRatio->Draw("copy");
//       histo2DTheoryEPOSRatio->GetYaxis()->SetRangeUser(0.,1.);
//       histo2DTheoryEPOSRatio->GetXaxis()->SetRangeUser(0.4,11.);
//
//         graphChargedRatioKaonToPionSys0010->Draw("2same");
//         graphChargedRatioKaonToPion0010->Draw("p,same");
//
//         if(meson.CompareTo("Eta")==0){
//           graphCombEtatoPi0SysPbPb2760GeV_0010->Draw("E2same");
//           graphCombEtatoPi0StatPbPb2760GeV_0010->Draw("p,same");
//         }
//
//         histoEtaPi0ratio->Draw("same,l");
//         histoKaonsPionsratio->Draw("same,l");
//
//         TLegend* legEPOSRatioWData = new TLegend(0.2,0.77,0.35,0.93);
//         legEPOSRatioWData->SetFillColor(0);
//         legEPOSRatioWData->SetLineColor(0);
//         legEPOSRatioWData->SetTextFont(42);
//         legEPOSRatioWData->SetTextSize(FontSize);
//         legEPOSRatioWData->AddEntry(histoEtaPi0ratio,"#eta/#pi^{0}","l");
//         legEPOSRatioWData->AddEntry(histoKaonsPionsratio,"K^{#pm}/#pi^{#pm}","l");
//         legEPOSRatioWData->AddEntry(histoEtaKaonratio,"#eta/K^{#pm}","l");
//         legEPOSRatioWData->Draw();
//
//     canvasTheoryEPOSRatio->SaveAs(Form("%s/EPOSTheoryChargedRatiosWithData.%s",outputDir.Data(),suffix.Data()));
//     canvasTheoryEPOSRatio->SaveAs(Form("%s/EPOSTheoryChargedRatiosWithData.%s",PubNotePlots.Data(),suffix.Data()));

/*
    canvasTheoryEPOSRatio->cd();
      histo2DTheoryEPOSRatio->Draw("copy");
      histo2DTheoryEPOSRatio->GetYaxis()->SetRangeUser(0., 3.);
      a  = NULL;
      b  = NULL;
      c  = NULL;
      d  = NULL;

      TGraphErrors *graphEPOS2old_0010 = (TGraphErrors*)gEPOS_Spec_0010->Clone();
      TGraphErrors *graphEPOS2new_0010 = new TGraphErrors(graphEPOSpred_0010);
      graphEPOS2old_0010->Print();
      graphEPOS2new_0010->Print();
//       Double_t *yPi0EPOS2old = graphEPOS2old_0010->GetY();
//       Double_t *yPi0EPOS2new = graphEPOS2new_0010->GetY();
//       Double_t *xPi0EPOS2new = graphEPOS2new_0010->GetX();
//       TGraphErrors* graphEPOS2old2new = (TGraphErrors*)graphEPOS2new_0010->Clone("graphEPOS2old2new");
//       Int_t binnew = graphEPOS2old2new->GetN();
//       Double_t yPi0EPOS2Ratioold2new[binnew];
//
//       for(Int_t a = 0; a<graphEPOS2old2new->GetN();a++){
//
//         cout << "from old " << yPi0EPOS2old[a] << endl;
//         cout << "from new " << yPi0EPOS2new[a] << endl;
// //         yPi0EPOS2Ratioold2new[a] = yPi0EPOS2old[a]/yPi0EPOS2new[a];
// //         graphEPOS2old2new->SetPoint(a,xPi0EPOS2new[a],yPi0EPOS2Ratioold2new[a]);
//
//       }
//       return;
      TGraphErrors* histoPi0EPOS2 = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                  graphEPOS2new_0010, graphEPOS2new_0010,
                                                  graphEPOS2old_0010, graphEPOS2old_0010,
                                                  kFALSE,  kFALSE, &a, &b, &c, &d);
      DrawGammaSetMarkerTGraphErr(histoPi0EPOS2, 24,1.2, kBlack, kBlack,2);
      histoPi0EPOS2->Draw("same,p");
      histoPi0EPOS2->Print();
//       return;

      a  = NULL;
      b  = NULL;
      c  = NULL;
      d  = NULL;
      TGraphErrors* histoPi0EPOS3to2new = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                  graphEPOSpred_0010, graphEPOSpred_0010,
                                                  graphEPOS3Pi0_0010, graphEPOS3Pi0_0010,
                                                  kFALSE,  kFALSE, &a, &b, &c, &d);
      DrawGammaSetMarkerTGraphErr(histoPi0EPOS3to2new, 24,1.2, kBlue, kBlue,2);
      histoPi0EPOS3to2new->Draw("same,l");
      histoPi0EPOS3to2new->Print();
//       return;


      a  = NULL;
      b  = NULL;
      c  = NULL;
      d  = NULL;
      TGraphErrors* histoPi0EPOS3to2old = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                  graphEPOS2old_0010, graphEPOS2old_0010,
                                                  graphEPOS3Pi0_0010, graphEPOS3Pi0_0010,
                                                  kFALSE,  kFALSE, &a, &b, &c, &d);
      DrawGammaSetMarkerTGraphErr(histoPi0EPOS3to2old, 24,1.2, kRed, kRed,2);
      histoPi0EPOS3to2old->Draw("same,l");
      histoPi0EPOS3to2old->Print();
//       return;
*/

// //       a  = NULL;
// //       b  = NULL;
// //       c  = NULL;
// //       d  = NULL;
// //       TGraphErrors* histoPi0EPOS2_2040to2050 = CalculateRatioBetweenSpectraWithDifferentBinning(
// //                                                   gEPOS_SpecScaled_2040, gEPOS_SpecScaled_2040,
// //                                                   graphEPOSpred_2050, graphEPOSpred_2050,
// //                                                   kTRUE,  kTRUE, &a, &b, &c, &d);
// //       DrawGammaSetMarkerTGraphErr(histoPi0EPOS2_2040to2050, 24,1.2, kGreen+2, kGreen+2,2);
// //       histoPi0EPOS2_2040to2050->Draw("same,l");
//
//
//       TLegend* legEPOSratio = new TLegend(0.2,0.15,0.35,0.35);
//       legEPOSratio->SetFillColor(0);
//       legEPOSratio->SetLineColor(0);
//       legEPOSratio->SetTextFont(42);
//       legEPOSratio->SetTextSize(FontSize);
//       legEPOSratio->AddEntry(histoPi0EPOS2,"EPOS2 old/new - 0#font[122]{-}10%","lp");
//       legEPOSratio->AddEntry(histoPi0EPOS3to2new,"EPOS3/EPOS2(new) - 0#font[122]{-}10%","pl");
//       legEPOSratio->AddEntry(histoPi0EPOS3to2old,"EPOS3/EPOS2(old) - 0#font[122]{-}10%","pl");
// //       legEPOSratio->AddEntry(histoPi0EPOS2_2040to2050,"EPOS2 old 20-40% /new 20#font[122]{-}50% ","l");
//       legEPOSratio->Draw();
//
//
//     canvasTheoryEPOSRatio->SaveAs(Form("%s/EPOSTheoryRatios.%s",outputDir.Data(),suffix.Data()));


//     Double_t arrayBoundariesX1_MassWidth[3];
//     Double_t arrayBoundariesY1_MassWidth[3];
//     Double_t relativeMarginsX_RAA4pads_MassWidth[3];
//     Double_t relativeMarginsY_RAA4pads_MassWidth[3];
//     textSizeLabelsPixel = 90;
//     ReturnCorrectValuesForCanvasScaling(2400,1600, 2, 2,0.25, 0.2, 0.003,0.2,arrayBoundariesX1_MassWidth,arrayBoundariesY1_MassWidth,relativeMarginsX_RAA4pads_MassWidth,relativeMarginsY_RAA4pads_MassWidth);
//
//     TCanvas * canvas4PartMassWidth = new TCanvas("canvas4PartMassWidth","",0,0,2400,1600);  // gives the page size
// //     DrawGammaCanvasSettings( canvas4PartMassWidth, 0.13, 0.0, 0.02, 0.1);
//     canvas4PartMassWidth->cd();
//
//     TPad* pad4PartMassWidth1 = new TPad("pad4PartMassWidth1", "", arrayBoundariesX_RAA4pads[0], arrayBoundariesY_RAA4pads[1], arrayBoundariesX_RAA4pads[1], arrayBoundariesY_RAA4pads[0],-1, -1, -2);
//     DrawGammaPadSettings( pad4PartMassWidth1, relativeMarginsX_RAA4pads[0], relativeMarginsX_RAA4pads[1], relativeMarginsY_RAA4pads[0], relativeMarginsY_RAA4pads[1]);
//     pad4PartMassWidth1->Draw();
//
//     TPad* pad4PartMassWidth2 = new TPad("pad4PartMassWidth2", "", arrayBoundariesX_RAA4pads[0], arrayBoundariesY_RAA4pads[2], arrayBoundariesX_RAA4pads[1], arrayBoundariesY_RAA4pads[1],-1, -1, -2);
//     DrawGammaPadSettings( pad4PartMassWidth2, relativeMarginsX_RAA4pads[0], relativeMarginsX_RAA4pads[1], relativeMarginsY_RAA4pads[1], relativeMarginsY_RAA4pads[2]);
//     pad4PartMassWidth2->Draw();
//
//     TPad* pad4PartMassWidth3 = new TPad("pad4PartMassWidth3", "", arrayBoundariesX_RAA4pads[1], arrayBoundariesY_RAA4pads[1], arrayBoundariesX_RAA4pads[2], arrayBoundariesY_RAA4pads[0],-1, -1, -2);
//     DrawGammaPadSettings( pad4PartMassWidth3, relativeMarginsX_RAA4pads[1], relativeMarginsX_RAA4pads[1], relativeMarginsY_RAA4pads[0], relativeMarginsY_RAA4pads[1]);
//     pad4PartMassWidth3->Draw();
//
//     TPad* pad4PartMassWidth4 = new TPad("pad4PartMassWidth4", "", arrayBoundariesX_RAA4pads[1], arrayBoundariesY_RAA4pads[2], arrayBoundariesX_RAA4pads[2], arrayBoundariesY_RAA4pads[1],-1, -1, -2);
//     DrawGammaPadSettings( pad4PartMassWidth4, relativeMarginsX_RAA4pads[1], relativeMarginsX_RAA4pads[1], relativeMarginsY_RAA4pads[1], relativeMarginsY_RAA4pads[2]);
//     pad4PartMassWidth4->Draw();
//
//
//
//     Int_t textSizeLabelsPixelMass = 50;
//     Double_t marginMass = 0.13*2400;
//     Double_t textsizeLabelsMass = 0;
//     Double_t textsizeFacMass = 0;
//     Double_t textsizeLabelsWidth = 0;
//     Double_t textsizeFacWidth = 0;
//
//     if (pad4PartMassWidth1->XtoPixel(pad4PartMassWidth1->GetX2()) < pad4PartMassWidth1->YtoPixel(pad4PartMassWidth1->GetY1())){
//         textsizeLabelsWidth = (Double_t)textSizeLabelsPixelMass/pad4PartMassWidth1->XtoPixel(pad4PartMassWidth1->GetX2()) ;
//         textsizeFacWidth = (Double_t)1./pad4PartMassWidth1->XtoPixel(pad4PartMassWidth1->GetX2()) ;
//     } else {
//         textsizeLabelsWidth = (Double_t)textSizeLabelsPixelMass/pad4PartMassWidth1->YtoPixel(pad4PartMassWidth1->GetY1());
//         textsizeFacWidth = (Double_t)1./pad4PartMassWidth1->YtoPixel(pad4PartMassWidth1->GetY1());
//     }
//     if (pad4PartMassWidth2->XtoPixel(pad4PartMassWidth2->GetX2()) < pad4PartMassWidth2->YtoPixel(pad4PartMassWidth2->GetY1())){
//         textsizeLabelsMass = (Double_t)textSizeLabelsPixelMass/pad4PartMassWidth2->XtoPixel(pad4PartMassWidth2->GetX2()) ;
//         textsizeFacMass = (Double_t)1./pad4PartMassWidth2->XtoPixel(pad4PartMassWidth2->GetX2()) ;
//     } else {
//         textsizeLabelsMass = (Double_t)textSizeLabelsPixelMass/pad4PartMassWidth2->YtoPixel(pad4PartMassWidth2->GetY1());
//         textsizeFacMass = (Double_t)1./pad4PartMassWidth2->YtoPixel(pad4PartMassWidth2->GetY1());
//     }
//
//     cout << textsizeLabelsMass << endl;
//
//     TPad* padFWHMLegend1 = new TPad("padFWHMLegend1", "", 0.07, 0.815, 0.35, 0.93,-1, -1, -2);
//     DrawGammaPadSettings( padFWHMLegend1, 0., 0., 0., 0.);
//     padFWHMLegend1->Draw();
//
//
//     TH2D *histo2DPi0FWHM;
//     histo2DPi0FWHM = new TH2D("histo2DPi0FWHM", "histo2DPi0FWHM", 20,0.35,20. ,1000.,-30,40);
//     SetStyleHistoTH2ForGraphs(histo2DPi0FWHM, "#it{p}_{T} (GeV/#it{c})","peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth,textsizeLabelsWidth, 0.85*textsizeLabelsWidth,textsizeLabelsWidth, 1,0.3/(textsizeFacWidth*marginMass), 515, 504);
//     histo2DPi0FWHM->GetYaxis()->SetRangeUser(-1.,18);
//     histo2DPi0FWHM->GetYaxis()->SetLabelOffset(0.01);
//     histo2DPi0FWHM->GetYaxis()->SetLabelFont(42);
//     histo2DPi0FWHM->GetXaxis()->SetLabelFont(42);
//
//     TH2D *histo2DPi0Mass;
//     histo2DPi0Mass = new TH2D("histo2DPi0Mass", "histo2DPi0Mass", 20,0.35,20. ,1000.,125.,150);
//     SetStyleHistoTH2ForGraphs(histo2DPi0Mass, "#it{p}_{T} (GeV/#it{c})","peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass,textsizeLabelsMass, 0.85*textsizeLabelsMass,textsizeLabelsMass, 0.9,0.3/(textsizeFacMass*marginMass), 515, 510);
//     histo2DPi0Mass->GetYaxis()->SetRangeUser(128.,143.5);
//     histo2DPi0Mass->GetXaxis()->SetLabelOffset(-0.02);
//     histo2DPi0Mass->GetYaxis()->SetLabelOffset(0.01);
//     histo2DPi0Mass->GetYaxis()->SetLabelFont(42);
//     histo2DPi0Mass->GetXaxis()->SetLabelFont(42);
//
//     pad4PartMassWidth1->cd();
//     pad4PartMassWidth1->SetLogx();
//     histo2DPi0FWHM->DrawCopy();
//
// //     DrawGammaSetMarker(histoPHOSWidthDataPP, markerStylePHOS, markerSizeMass, colorPHOSMass, colorPHOSMass);
// //     histoPHOSWidthDataPP->DrawCopy("same,p");
// //     DrawGammaSetMarker(histoPHOSWidthMCPP, markerStylePHOSMC, markerSizeMass, colorPHOSMCMass , colorPHOSMCMass);
// //     histoPHOSWidthMCPP->DrawCopy("same,p");
// //
// //     DrawGammaSetMarker(histoPCMWidthDataPP, markerStyleConv, markerSizeMass, colorConv, colorConv);
// //     histoPCMWidthDataPP->DrawCopy("same,p");
// //     DrawGammaSetMarker(histoPCMWidthMCPP, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);
// //     histoPCMWidthMCPP->DrawCopy("same,p");
// //
// //     DrawGammaSetMarker(histoPHOSWidthData0010, markerStylePHOS, markerSizeMass, colorPHOSMass, colorPHOSMass);
// //     DrawGammaSetMarker(histoPHOSWidthMC0010, markerStylePHOSMC, markerSizeMass, colorPHOSMCMass , colorPHOSMCMass);
// //
// //     TLatex *labelMassPi0PP = new TLatex(0.2,0.88,collisionSystemPP.Data());
// //     SetStyleTLatex( labelMassPi0PP, 0.85*textsizeLabelsWidth,4);
// //     labelMassPi0PP->Draw();
// //     TLatex *labelLegendAMass = new TLatex(0.92,0.88,"a)");
// //     SetStyleTLatex( labelLegendAMass,0.85*textsizeLabelsWidth,4);
// //     labelLegendAMass->Draw();
// //
// //     //********************************** Defintion of the Legend **************************************************
// //     Double_t columnsLegendFWHM[4]   = {0.,0.2,0.37,0.55};
// //     Double_t rowsLegendFWHM[3]      = {0.66,0.33,0.0};
// //     //******************* Text sizes *******************
// //     Size_t textSizeLeftColumnFWHM   = 0.301;
// //     Size_t textSizeTopRowFWHM   = 0.301;
// //     Size_t textSizeSecondRowFWHM    = 0.301;
// //     //******************* Offsets ***********************
// //     Double_t offsetMarkerXFWHM  = 0.07;
// //     Double_t offsetMarkerYFWHM  = 0.07;
// //     //****************** Scale factors ******************
// //     Double_t scaleMarkerFWHM        = 1.;
// //
// //     padFWHMLegend1->cd();
// //     //****************** first Column **************************************************
// //     TLatex *textFWHMCTS = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM[1],"PCM");
// //     SetStyleTLatex( textFWHMCTS, textSizeLeftColumnFWHM,4);
// //     textFWHMCTS->Draw();
// //     TLatex *textFWHMPHOS = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM[2],"PHOS");
// //     SetStyleTLatex( textFWHMPHOS, textSizeLeftColumnFWHM,4);
// //     textFWHMPHOS->Draw();
// //
// //     //****************** second Column *************************************************
// //     TLatex *textFWHMData2 = new TLatex(columnsLegendFWHM[1],rowsLegendFWHM[0] ,"Data");
// //     SetStyleTLatex( textFWHMData2, textSizeTopRowFWHM ,4);
// //     textFWHMData2->Draw();
// //     TLatex *textFWHMMC2 = new TLatex(columnsLegendFWHM[2] ,rowsLegendFWHM[0],"MC");
// //     SetStyleTLatex( textFWHMMC2, textSizeTopRowFWHM,4);
// //     textFWHMMC2->Draw();
// //
// //     TMarker* markerCTSPi0FWHM = CreateMarkerFromHisto(histoPCMWidthDataPP,columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
// //     markerCTSPi0FWHM->DrawMarker(columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM);
// //     TMarker* markerPHOSPi0FWHM = CreateMarkerFromHisto(histoPHOSWidthData0010,columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
// //     markerPHOSPi0FWHM->DrawMarker(columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM);
// //
// //     TMarker* markerCTSPi0FWHMMC = CreateMarkerFromHisto(histoPCMWidthMCPP,columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
// //     markerCTSPi0FWHMMC->DrawMarker(columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM);
// //     TMarker* markerPHOSPi0FWHMMC = CreateMarkerFromHisto(histoPHOSWidthMC0010,columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
// //     markerPHOSPi0FWHMMC->DrawMarker(columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM);
// //
// //     TLatex *textWidthConv2 = new TLatex(columnsLegendFWHM[3],rowsLegendFWHM[1] ,"FWHM/2.35");
// //     SetStyleTLatex( textWidthConv2, textSizeSecondRowFWHM,4);
// //     textWidthConv2->Draw();
// //     TLatex *textWidthPHOS2 = new TLatex(columnsLegendFWHM[3] ,rowsLegendFWHM[2],"#sigma");
// //     SetStyleTLatex( textWidthPHOS2, textSizeSecondRowFWHM,4);
// //     textWidthPHOS2->Draw();
//
//
//     pad4PartMassWidth1->cd();
//     pad4PartMassWidth1->SetLogx();
//     histo2DPi0FWHM->DrawCopy();
//
//     DrawGammaSetMarker(histoPCMPi0FWHMMeVPbPb2760GeV_0010, markerStyleConv, markerSizeMass, colorConv, colorConv);
//     histoPCMPi0FWHMMeVPbPb2760GeV_0010->DrawCopy("same,p");
//     DrawGammaSetMarker(histoPCMPi0TrueFWHMMeVPbPb2760GeV_0010, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);
//     histoPCMPi0TrueFWHMMeVPbPb2760GeV_0010->DrawCopy("same,p");
//
//     TLatex *labelMassPi0PbPb0005 = new TLatex(0.05,0.88,collisionSystemPbPb0010.Data());
//     SetStyleTLatex( labelMassPi0PbPb0005, 0.85*textsizeLabelsWidth,4);
//     labelMassPi0PbPb0005->Draw();
//     TLatex *labelLegendBMass = new TLatex(0.89,0.88,"c)");
//     SetStyleTLatex( labelLegendBMass, 0.85*textsizeLabelsWidth,4);
//     labelLegendBMass->Draw();
//
//     pad4PartMassWidth1->Update();
//     pad4PartMassWidth2->cd();
//     pad4PartMassWidth2->SetLogx();
//     histo2DPi0Mass->DrawCopy();
//
//     DrawGammaSetMarker(histoPCMPi0MassPbPb2760GeV_0010 , markerStyleConv, markerSizeMass, colorConv, colorConv);
//     histoPCMPi0MassPbPb2760GeV_0010->DrawCopy("same,p");
//     DrawGammaSetMarker(histoPCMPi0TrueMassPbPb2760GeV_0010 , markerStyleConvMC , markerSizeMass, colorConvMC, colorConvMC);
//     histoPCMPi0TrueMassPbPb2760GeV_0010->DrawCopy("same,p");
//
//     DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,colorConv);
//     TLatex *labelLegendFMass = new TLatex(0.89,0.9,"f)");
//     SetStyleTLatex( labelLegendFMass, 0.85*textsizeLabelsMass,4);
//     labelLegendFMass->Draw();
//
//     pad4PartMassWidth2->Update();
//
//     pad4PartMassWidth3->cd();
//     pad4PartMassWidth3->SetLogx();
//     histo2DPi0FWHM->DrawCopy();
//
//     DrawGammaSetMarker(histoPCMEtaFWHMMeVPbPb2760GeV_0010, markerStyleConv, markerSizeMass, colorConv, colorConv);
//     histoPCMEtaFWHMMeVPbPb2760GeV_0010->DrawCopy("same,p");
//     DrawGammaSetMarker(histoPCMEtaTrueFWHMMeVPbPb2760GeV_0010, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);
//     histoPCMEtaTrueFWHMMeVPbPb2760GeV_0010->DrawCopy("same,p");
//
//     TLatex *labelMassPi0PbPb6080 = new TLatex(0.05,0.88,collisionSystemCent6080.Data());
//     SetStyleTLatex( labelMassPi0PbPb6080, 0.85*textsizeLabelsWidth,4);
//     labelMassPi0PbPb6080->Draw();
//     TLatex *labelLegendCMass = new TLatex(0.91,0.88,"b)");
//     SetStyleTLatex( labelLegendCMass, 0.85*textsizeLabelsWidth,4);
//     labelLegendCMass->Draw();
//
//     pad4PartMassWidth3->Update();
//     pad4PartMassWidth4->cd();
//     pad4PartMassWidth4->SetLogx();
//     histo2DPi0Mass->DrawCopy();
//
//     DrawGammaSetMarker(histoPCMPi0MassPbPb2760GeV_0010 , markerStyleConv, markerSizeMass, colorConv, colorConv);
//     histoPCMPi0MassPbPb2760GeV_0010->DrawCopy("same,p");
//     DrawGammaSetMarker(histoPCMPi0TrueMassPbPb2760GeV_0010 , markerStyleConvMC , markerSizeMass, colorConvMC, colorConvMC);
//     histoPCMPi0TrueMassPbPb2760GeV_0010->DrawCopy("same,p");
//
//     DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,colorConv);
//     TLatex *labelLegendEMass = new TLatex(0.91,0.9,"e)");
//     SetStyleTLatex( labelLegendEMass, 0.85*textsizeLabelsMass,4);
//     labelLegendEMass->Draw();
//
//     pad4PartMassWidth4->Update();
//
//     canvas4PartMassWidth->Update();
//     canvas4PartMassWidth->SaveAs(Form("%s/MassWidth_4Parted_Paper.%s",paperPlots.Data(),suffix.Data()));
//     canvas4PartMassWidth->SaveAs(Form("%s/MassWidth_4Parted_Paper.%s",outputDir.Data(),suffix.Data()));
//     delete pad4PartMassWidth1;
//     delete pad4PartMassWidth2;
//     delete pad4PartMassWidth3;
//     delete pad4PartMassWidth4;
//     delete canvas4PartMassWidth;



    cout << "Parameter Bylinkin fit 0010: " << endl;
    cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_0010)<< endl << endl;
    cout << "Parameter Bylinkin fit 2050: " << endl;
    cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_2050)<< endl << endl;
//     cout << "Parameter low component Bylinkin fit 0010: " << endl;
//     cout << WriteParameterToFile(fitLowPtBylinkin_0010)<< endl << endl;
//     cout << "Parameter high component Bylinkin fit 0010: " << endl;
//     cout << WriteParameterToFile(fitHighPtBylinkin_0010)<< endl << endl;
//     cout << "Parameter low component Bylinkin fit 2050: " << endl;
//     cout << WriteParameterToFile(fitLowPtBylinkin_2050)<< endl << endl;
//     cout << "Parameter high component Bylinkin fit 2050: " << endl;
//     cout << WriteParameterToFile(fitHighPtBylinkin_2050)<< endl << endl;
//     cout << "low 0010: " << fitLowPtBylinkin_0010->Eval(0) << "  high 0010: " << fitHighPtBylinkin_0010->Eval(0) << endl;
//     cout << "low 2050: " << fitLowPtBylinkin_2050->Eval(0) << "  high 2050: " << fitHighPtBylinkin_2050->Eval(0) << endl;

//        cout << WriteParameterToFile(fitQCDInvYieldPbPb2760GeV_0010)<< endl << endl;
//        cout << WriteParameterToFile(fitQCDInvYieldPbPb2760GeV_2050)<< endl << endl;





    if(meson.CompareTo("Pi0")==0 && (thesisPlotting || PaperPi0)){

        while(graphPCMInvYieldStatPbPb2760GeVUnShifted_0010->GetX()[0] < /* 0.8 */ 1.)
          graphPCMInvYieldStatPbPb2760GeVUnShifted_0010->RemovePoint(0);
        while(graphPCMInvYieldSysPbPb2760GeVUnShifted_0010->GetX()[0] < /* 0.8 */ 1.)
          graphPCMInvYieldSysPbPb2760GeVUnShifted_0010->RemovePoint(0);
        while(graphPCMInvYieldStatPbPb2760GeVUnShifted_2050->GetX()[0] < /* 0.8 */ 1.)
          graphPCMInvYieldStatPbPb2760GeVUnShifted_2050->RemovePoint(0);
        while(graphPCMInvYieldSysPbPb2760GeVUnShifted_2050->GetX()[0] < /* 0.8 */ 1.)
          graphPCMInvYieldSysPbPb2760GeVUnShifted_2050->RemovePoint(0);
        while(graphPCMInvYieldStatPbPb2760GeVUnShifted_2040->GetX()[0] < /* 0.8 */ 1.)
          graphPCMInvYieldStatPbPb2760GeVUnShifted_2040->RemovePoint(0);
        while(graphPCMInvYieldSysPbPb2760GeVUnShifted_2040->GetX()[0] < /* 0.8 */ 1.)
          graphPCMInvYieldSysPbPb2760GeVUnShifted_2040->RemovePoint(0);
    }

    //******************************************************************************//
    //************************* Saving of final results ****************************//
    //******************************************************************************//
    TFile *fCombResults = new TFile(Form("%s/CombinedResultsPaperPbPb2760GeV_%s.root", outputDir.Data(),dateForOutput.Data()), "UPDATE");

        //x shifted spectra
        graphCombInvYieldTotPbPb2760GeV_0010->Write(Form("graphInvYield%sCombPbPb2760GeV_0010",meson.Data()));
        graphCombInvYieldStatPbPb2760GeV_0010->Write(Form("graphInvYield%sCombPbPb2760GeVStatErr_0010",meson.Data()));
        graphCombInvYieldSysPbPb2760GeV_0010->Write(Form("graphInvYield%sCombPbPb2760GeVSysErr_0010",meson.Data()));
        graphCombInvYieldTotPbPb2760GeV_2050->Write(Form("graphInvYield%sCombPbPb2760GeV_2050",meson.Data()));
        graphCombInvYieldStatPbPb2760GeV_2050->Write(Form("graphInvYield%sCombPbPb2760GeVStatErr_2050",meson.Data()));
        graphCombInvYieldSysPbPb2760GeV_2050->Write(Form("graphInvYield%sCombPbPb2760GeVSysErr_2050",meson.Data()));

        graphPCMPi0InvYieldStatPbPb2760GeV_0010->Write("graphInvYieldPi0PCMPbPb2760GeVStatErr_0010");
        graphPCMPi0InvYieldSysPbPb2760GeV_0010->Write("graphInvYieldPi0PCMPbPb2760GeVSysErr_0010");
        graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Write("graphInvYieldPi0EMCalPbPb2760GeVStatErr_0010");
        graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Write("graphInvYieldPi0EMCalPbPb2760GeVSysErr_0010");
        graphPHOSPi0InvYieldStatPbPb2760GeV_0010->Write("graphInvYieldPi0PHOSPbPb2760GeVStatErr_0010");
        graphPHOSPi0InvYieldSysPbPb2760GeV_0010->Write("graphInvYieldPi0PHOSPbPb2760GeVSysErr_0010");
        graphPCMPi0InvYieldStatPbPb2760GeV_2050->Write("graphInvYieldPi0PCMPbPb2760GeVStatErr_2050");
        graphPCMPi0InvYieldSysPbPb2760GeV_2050->Write("graphInvYieldPi0PCMPbPb2760GeVSysErr_2050");
        graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Write("graphInvYieldPi0EMCalPbPb2760GeVStatErr_2050");
        graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Write("graphInvYieldPi0EMCalPbPb2760GeVSysErr_2050");

        graphPCMEtaInvYieldStatPbPb2760GeV_0010->Write("graphInvYieldEtaPCMPbPb2760GeVStatErr_0010");
        graphPCMEtaInvYieldSysPbPb2760GeV_0010->Write("graphInvYieldEtaPCMPbPb2760GeVSysErr_0010");
        graphEMCalEtaInvYieldStatPbPb2760GeV_0010->Write("graphInvYieldEtaEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Write("graphInvYieldEtaEMCalPbPb2760GeVSysErr_0010");
        graphPCMEtaInvYieldStatPbPb2760GeV_2050->Write("graphInvYieldEtaPCMPbPb2760GeVStatErr_2050");
        graphPCMEtaInvYieldSysPbPb2760GeV_2050->Write("graphInvYieldEtaPCMPbPb2760GeVSysErr_2050");
        graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Write("graphInvYieldEtaEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Write("graphInvYieldEtaEMCalPbPb2760GeVSysErr_2050");

        //y shifted spectra
        graphCombInvYieldSysPbPb2760GeVYShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVSysErrYshifted_0010",meson.Data()));
        graphCombInvYieldStatPbPb2760GeVYShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVStatErrYshifted_0010",meson.Data()));
        graphCombInvYieldSysPbPb2760GeVYShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVSysErrYshifted_2050",meson.Data()));
        graphCombInvYieldStatPbPb2760GeVYShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVStatErrYshifted_2050",meson.Data()));

        fitBylinkinPbPb2760GeVPtLHC11h_0010->Write(Form("FitToYield%s_0010",meson.Data()));
        fitBylinkinPbPb2760GeVPtLHC11h_2050->Write(Form("FitToYield%s_2050",meson.Data()));

        if(mTScaledEtaFromPi0)mTScaledEtaFromPi0->Write("mTScaledEtaFromPi0_0010");

        graphCombRAASysPbPb2760GeV_0010->Write(Form("graphRAA%sCombPbPb2760GeVSysErr_0010",meson.Data()));
        graphCombRAAStatPbPb2760GeV_0010->Write(Form("graphRAA%sCombPbPb2760GeVStatErr_0010",meson.Data()));
        graphCombRAASysPbPb2760GeV_2050->Write(Form("graphRAA%sCombPbPb2760GeVSysErr_2050",meson.Data()));
        graphCombRAAStatPbPb2760GeV_2050->Write(Form("graphRAA%sCombPbPb2760GeVStatErr_2050",meson.Data()));

        graphPCMEtaRAASysPbPb2760GeV_2040->Write("graphPCMEtaRAASysPbPb2760GeV_2040");
        graphPCMEtaRAAStatPbPb2760GeV_2040->Write("graphPCMEtaRAAStatPbPb2760GeV_2040");
        graphPCMPi0RAASysPbPb2760GeV_2040->Write("graphPCMPi0RAASysPbPb2760GeV_2040");
        graphPCMPi0RAAStatPbPb2760GeV_2040->Write("graphPCMPi0RAAStatPbPb2760GeV_2040");

        graphRCSCEMCal->Write(Form("graph%sRCSCEMCal",meson.Data()));
        graphRCSCSysEMCal->Write(Form("graph%sRCSCSysEMCal",meson.Data()));

        if(graphCombEtatoPi0SysPbPb2760GeV_0010){
            while(graphCombEtatoPi0SysPbPb2760GeV_0010->GetY()[0] < 1e-14)graphCombEtatoPi0SysPbPb2760GeV_0010->RemovePoint(0);
            graphCombEtatoPi0SysPbPb2760GeV_0010->Write("graphEtaToPi0RatioCombPbPb2760GeVSysErr_0010");
//             graphCombEtatoPi0SysPbPb2760GeV_0010->Print();
        }
        if(graphCombEtatoPi0StatPbPb2760GeV_0010){
            while(graphCombEtatoPi0StatPbPb2760GeV_0010->GetY()[0] < 1e-14)graphCombEtatoPi0StatPbPb2760GeV_0010->RemovePoint(0);
            graphCombEtatoPi0StatPbPb2760GeV_0010->Write("graphEtaToPi0RatioCombPbPb2760GeVStatErr_0010");
//             graphCombEtatoPi0StatPbPb2760GeV_0010->Print();
        }
        if(graphCombEtatoPi0SysPbPb2760GeV_2050){
            while(graphCombEtatoPi0SysPbPb2760GeV_2050->GetY()[0] < 1e-14)graphCombEtatoPi0SysPbPb2760GeV_2050->RemovePoint(0);
            graphCombEtatoPi0SysPbPb2760GeV_2050->Write("graphEtaToPi0RatioCombPbPb2760GeVSysErr_2050");
            graphCombEtatoPi0SysPbPb2760GeV_2050->Print();
        }
        if(graphCombEtatoPi0StatPbPb2760GeV_2050){
            while(graphCombEtatoPi0StatPbPb2760GeV_2050->GetY()[0] < 1e-14)graphCombEtatoPi0StatPbPb2760GeV_2050->RemovePoint(0);
            graphCombEtatoPi0StatPbPb2760GeV_2050->Write("graphEtaToPi0RatioCombPbPb2760GeVStatErr_2050");
//             graphCombEtatoPi0StatPbPb2760GeV_2050->Print();
        }
        graphPCMEtatoPi0Stat2760GeV_0010->Write("graphPCMEtatoPi0Stat2760GeV_0010");
        graphPCMEtatoPi0Sys2760GeV_0010->Write("graphPCMEtatoPi0Sys2760GeV_0010");
        graphPCMEtatoPi0Stat2760GeV_2050->Write("graphPCMEtatoPi0Stat2760GeV_2050");
        graphPCMEtatoPi0Sys2760GeV_2050->Write("graphPCMEtatoPi0Sys2760GeV_2050");

        graphEPOSpred_0010->Write(Form("EPOSprediction%s_0010",meson.Data()));
        graphEPOSpred_2050->Write(Form("EPOSprediction%s_2050",meson.Data()));
        TheoryCracowLowPt_0010->Write(Form("Cracowprediction%s_0010",meson.Data()));
        TheoryCracowLowPt_2050->Write(Form("Cracowprediction%s_2050",meson.Data()));
        TheoryBegunEQ_0010->Write(Form("BegunEQprediction%s_0010",meson.Data()));
        TheoryBegunEQ_2050->Write(Form("BegunEQprediction%s_2050",meson.Data()));

        fCombResults->mkdir("UnshiftedSpectra");
        TDirectoryFile* directoryNoShift = (TDirectoryFile*)fCombResults->Get("UnshiftedSpectra");
        fCombResults->cd("UnshiftedSpectra");
            //unshifted spectra
            graphCombInvYieldTotPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVNoShift_0010",meson.Data()));
            graphCombInvYieldStatPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVStatErrNoShift_0010",meson.Data()));
            graphCombInvYieldSysPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVSysErrNoShift_0010",meson.Data()));
            graphCombInvYieldTotPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVNoShift_2050",meson.Data()));
            graphCombInvYieldStatPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVStatErrNoShift_2050",meson.Data()));
            graphCombInvYieldSysPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVSysErrNoShift_2050",meson.Data()));

            if(graphPCMInvYieldStatPbPb2760GeVUnShifted_0010)graphPCMInvYieldStatPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sPCMPbPb2760GeVStatErrNoShift_0010",meson.Data()));
            if(graphPCMInvYieldSysPbPb2760GeVUnShifted_0010)graphPCMInvYieldSysPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sPCMPbPb2760GeVSysErrNoShift_0010",meson.Data()));
            if(graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010)graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010->Write(Form("graphInvYield%sEMCalPbPb2760GeVStatErrNoShift_0010",meson.Data()));
            if(graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010)graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010->Write(Form("graphInvYield%sEMCalPbPb2760GeVSysErrNoShift_0010",meson.Data()));
            if(graphPHOSInvYieldStatPbPb2760GeVUnshifted_0010)graphPHOSInvYieldStatPbPb2760GeVUnshifted_0010->Write(Form("graphInvYield%sPHOSPbPb2760GeVStatErrNoShift_0010",meson.Data()));
            if(graphPHOSInvYieldSysPbPb2760GeVUnshifted_0010)graphPHOSInvYieldSysPbPb2760GeVUnshifted_0010->Write(Form("graphInvYield%sPHOSPbPb2760GeVSysErrNoShift_0010",meson.Data()));

            if(graphPCMInvYieldStatPbPb2760GeVUnShifted_2040)graphPCMInvYieldStatPbPb2760GeVUnShifted_2040->Write(Form("graphInvYield%sPCMPbPb2760GeVStatErrNoShift_2040",meson.Data()));
            if(graphPCMInvYieldSysPbPb2760GeVUnShifted_2040)graphPCMInvYieldSysPbPb2760GeVUnShifted_2040->Write(Form("graphInvYield%sPCMPbPb2760GeVSysErrNoShift_2040",meson.Data()));

            if(graphPCMInvYieldStatPbPb2760GeVUnShifted_2050)graphPCMInvYieldStatPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sPCMPbPb2760GeVStatErrNoShift_2050",meson.Data()));
            if(graphPCMInvYieldSysPbPb2760GeVUnShifted_2050)graphPCMInvYieldSysPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sPCMPbPb2760GeVSysErrNoShift_2050",meson.Data()));
            if(graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050)graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050->Write(Form("graphInvYield%sEMCalPbPb2760GeVStatErrNoShift_2050",meson.Data()));
            if(graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050)graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050->Write(Form("graphInvYield%sEMCalPbPb2760GeVSysErrNoShift_2050",meson.Data()));

            if(graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_0010)graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_0010->Write(Form("UnshiftedSysforRAAPCM%s_0010",meson.Data()));
            if(graphPHOSInvYieldSysPbPb2760GeVforRAAUnshifted_0010)graphPHOSInvYieldSysPbPb2760GeVforRAAUnshifted_0010->Write(Form("UnshiftedSysforRAAPHOS%s_0010",meson.Data()));
            if(graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_0010)graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_0010->Write(Form("UnshiftedSysforRAAEMCal%s_0010",meson.Data()));
            if(graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_2050)graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_2050->Write(Form("UnshiftedSysWOMatPCM%s_2050",meson.Data()));
            if(graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_2050)graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_2050->Write(Form("UnshiftedSysforRAAEMCal%s_2050",meson.Data()));


    fCombResults->Close();


}
