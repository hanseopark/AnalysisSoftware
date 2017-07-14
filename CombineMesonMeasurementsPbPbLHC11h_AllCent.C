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
// #include "CombineMesonMeasurementsPbPbLHC11hV1.h"
#include "CombineMesonMeasurementsPbPbLHC11h_AllCent.h"

extern TRandom* gRandom;
extern TBenchmark*  gBenchmark;
extern TSystem* gSystem;
extern TMinuit*     gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //  TString name;
};

void CombineMesonMeasurementsPbPbLHC11h_AllCent(TString meson = "Eta",
                                                TString fileNamePCM = "",
                                                TString suffix = "pdf",
                                                TString thisthesis="This thesis",
                                                Bool_t noXerrorBars = kTRUE,
                                                TString bWCorrection="X"){

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput           = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;

    Bool_t thesisPlotting;
    if(thisthesis.CompareTo("")==0) thesisPlotting = kFALSE;
    else thesisPlotting = kTRUE;
    
    Double_t maxX, minX;
    if(thesisPlotting){
        maxX = 35;
        minX = 0.5;
    } else {
        maxX = 35;
        minX = 0.25;
    }

    //___________________________________ Labels definition _____________________________________________

    Double_t FontSize = 0.035;
    TLatex *labelSystOnlyPbPb;
    if(meson.CompareTo("Pi0")==0){     labelSystOnlyPbPb= new TLatex(0.75,0.93,"#pi^{0} #rightarrow #gamma#gamma");}
    else if(meson.CompareTo("Eta")==0){labelSystOnlyPbPb= new TLatex(0.75,0.93,"#eta #rightarrow #gamma#gamma");}
    SetStyleTLatex(labelSystOnlyPbPb,FontSize,4);

    TLatex *labelSyst;
    if(meson.CompareTo("Pi0")==0){     labelSyst= new TLatex(0.6,0.93,"#pi^{0} #rightarrow #gamma#gamma");}
    else if(meson.CompareTo("Eta")==0){labelSyst= new TLatex(0.6,0.93,"#eta #rightarrow #gamma#gamma");}
    SetStyleTLatex(labelSyst,FontSize,4);

    TLatex *labelEnergyInvYieldSectionPi0LHC11h = new TLatex(0.6,0.88,collisionSystem2760GeV.Data());
    SetStyleTLatex( labelEnergyInvYieldSectionPi0LHC11h,FontSize,4);
    
    TLatex *labelDetSysInvYieldSectionPi0LHC11h;
    if(meson.CompareTo("Pi0")==0){
        labelDetSysInvYieldSectionPi0LHC11h= new TLatex(0.6,0.84,"#pi^{0} #rightarrow #gamma#gamma");
    } else if(meson.CompareTo("Eta")==0){
        labelDetSysInvYieldSectionPi0LHC11h= new TLatex(0.6,0.84,"#eta #rightarrow #gamma#gamma");
    }
    SetStyleTLatex( labelDetSysInvYieldSectionPi0LHC11h, FontSize,4);

    TLatex *labelEnergyInvYieldSectionPi0LHC11hnoPrelim = new TLatex(0.6,0.92,collisionSystem2760GeV.Data());
    SetStyleTLatex( labelEnergyInvYieldSectionPi0LHC11hnoPrelim, 0.035,4);
    
    TLatex *labelDetSysInvYieldSectionPi0LHC11hnoPrelim;
    if(meson.CompareTo("Pi0")==0){
        labelDetSysInvYieldSectionPi0LHC11hnoPrelim= new TLatex(0.6,0.88,"#pi^{0} #rightarrow #gamma#gamma");
    } else if(meson.CompareTo("Eta")==0){
        labelDetSysInvYieldSectionPi0LHC11hnoPrelim= new TLatex(0.6,0.88,"#eta #rightarrow #gamma#gamma");
    }
    SetStyleTLatex( labelDetSysInvYieldSectionPi0LHC11hnoPrelim, 0.035,4);

    TLatex *labelDetSysInvYieldSectionPi0LHC11hwithPP;
    if(meson.CompareTo("Pi0")==0){
        labelDetSysInvYieldSectionPi0LHC11hwithPP= new TLatex(0.62,0.88,"#pi^{0} #rightarrow #gamma#gamma");
    } else if(meson.CompareTo("Eta")==0){
        labelDetSysInvYieldSectionPi0LHC11hwithPP= new TLatex(0.62,0.88,"#eta #rightarrow #gamma#gamma");
    }
    SetStyleTLatex( labelDetSysInvYieldSectionPi0LHC11hwithPP, 0.035,4);

    TLatex *labelEtaToPi0Energy = new TLatex(0.12,0.89,collisionSystem2760GeV.Data());
    SetStyleTLatex( labelEtaToPi0Energy,0.04,4);

//     GetColorDefaultColor
//     GetDefaultMarkerStyle
//     GetDefaultMarkerSize

    TString cent[5] = {"0-5%","5-10%","0-10%","20-40%","20-50%"};
    TString labelcent[5] = {"0#font[122]{-}5%","5#font[122]{-}10%","0#font[122]{-}10%","20#font[122]{-}40%","20#font[122]{-}50%"};    
    TString collisionSystPbPb[5] = { "0#font[122]{-}5% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","5#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","0#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV", "20#font[122]{-}40% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","20#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV"};

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
    Double_t mass;
    if (meson.CompareTo("Pi0")==0){
        mass = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    } else if (meson.CompareTo("Eta")==0){
        mass = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    }
    TString energy= "PbPb_2.76TeV";

    Double_t normErr0005 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0005/tAA0005),2)+pow((commonCentralityErr0005/100),2),0.5);
    Double_t normErr0510 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0510/tAA0510),2)+pow((commonCentralityErr0510/100),2),0.5);
    Double_t normErr0010 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0010/tAA0010),2)+pow((commonCentralityErr0010/100),2),0.5);
    Double_t normErr2040 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2040/tAA2040),2)+pow((commonCentralityErr2040/100),2),0.5);
    Double_t normErr2050 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2050/tAA2050),2)+pow((commonCentralityErr2040/100),2),0.5);

    TBox* boxErrorNorm0005 = CreateBoxConv(GetColorDefaultColor(energy.Data(),"",cent[0]), 0.01, 1.-normErr0005 , 0.17, 1.+normErr0005);
    TBox* boxErrorNorm0510 = CreateBoxConv(GetColorDefaultColor(energy.Data(),"",cent[1]), 0.17, 1.-normErr0510 , 0.29, 1.+normErr0510);
    TBox* boxErrorNorm0010 = CreateBoxConv(GetColorDefaultColor(energy.Data(),"",cent[2]), 0.29, 1.-normErr0010 , 0.31, 1.+normErr0010);
    TBox* boxErrorNorm2040 = CreateBoxConv(GetColorDefaultColor(energy.Data(),"",cent[3]), 0.31, 1.-normErr2040 , 0.43, 1.+normErr2040);
    TBox* boxErrorNorm2050 = CreateBoxConv(GetColorDefaultColor(energy.Data(),"",cent[4]), 0.43, 1.-normErr2050 , 0.55, 1.+normErr2050);

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
    //*********************************************     Theory    *********************************************//
    //*********************************************************************************************************//
    TFile* fileTheoryGraphs = new TFile(fileNameTheory);
    TGraphErrors *Vitev_Bas_Raa_0005 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA0005");
    TGraphErrors *Vitev_Bas_Raa_0510 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA0510");
    TGraphErrors *Vitev_Bas_Raa_2040 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA2040");
    TGraphAsymmErrors* gWHDG_Raa_0510 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0510");
    TGraphAsymmErrors* gWHDG_Raa_0005 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0005");
    TGraphAsymmErrors* gWHDG_Raa_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0010");
    TGraphAsymmErrors* gWHDG_Raa_2040 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA2040");
    TGraphAsymmErrors* gWHDG_Raa_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA2050");

    TGraphAsymmErrors* graphPi0RAAJetQuenching_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphPi0RAAJetQuenching_0010");
    TGraphAsymmErrors* graphEtaRAAJetQuenching_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphEtaRAAJetQuenching_0010");

    gWHDG_Eta_Raa_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGetaRAA0010");
    gWHDG_Eta_Raa_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGetaRAA2050");

    TGraphAsymmErrors* graphPi0Djordjevic_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphPi0Djordjevic_0010");
    TGraphAsymmErrors* graphPi0Djordjevic_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphPi0Djordjevic_2050");


    //*********************************************************************************************************//
    //********************************        Other experiment inputs     *************************************//
    //*********************************************************************************************************//
    TFile* fileDataOtherEnergies = new TFile(fileNameDataOtherEnergyInput);
	TGraphErrors* graphPHENIX200GeVInvYield_0010       = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVInvYield_0010");
	TGraphErrors* graphPHENIX200GeVInvYield_2040       = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVInvYield_2040");
	TGraphErrors* graphPHENIX200GeVInvYield_2050       = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVInvYield_2050");
	TGraphErrors* graphPHENIX200GeVEtaInvYield_0020    = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVInvYield_0010");
	TGraphErrors* graphPHENIX200GeVEtaInvYield_2060    = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVInvYield_2050");

    TGraphErrors* graphPHENIX39GeVInvYield_0010        = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVInvYield_0010");
	TGraphErrors* graphPHENIX39GeVInvYield_2040        = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVInvYield_2040");
	TGraphErrors* graphPHENIX62GeVInvYield_0010        = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVInvYield_0010");
	TGraphErrors* graphPHENIX62GeVInvYield_2040        = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVInvYield_2040");
    
    graphWA98_17_3GeVPi0RAA_0013                       = (TGraphErrors*)fileDataOtherEnergies->Get("graphWA98RAA_0013");
    graphPHENIX200GeVPi0RAA_0010                       = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_0010");
    graphPHENIX200GeVPi0RAA_2040                       = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_2040");
    graphPHENIX39GeVPi0RAA_0010                        = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_0010");
    graphPHENIX39GeVPi0RAA_2040                        = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_2040");
    graphPHENIX62GeVPi0RAA_0010                        = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_0010");
    graphPHENIX62GeVPi0RAA_2040                        = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_2040");
    
    graphPHENIX200GeVEtaToPi0Ratio_0020                = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaToPi0Ratio_0020");
    graphPHENIX200GeVEtaToPi0Ratio_2060                = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaToPi0Ratio_2060");

    graphPHENIX200GeVEtaRAA_0020                       = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_0020");
    graphPHENIX200GeVEtaRAA_2060                       = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_2060");
    graphPHENIX200GeVEtaRAA_0010                       = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_0010");
    graphPHENIX200GeVEtaRAA_2040                       = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_2040");
    graphPHENIX200GeVEtaHighPtRAA_2060                 = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAAHighPt_2060");


    TFile* fileDataALICE = new TFile(fileNameALICEData);
    //*********************************************************************************************************//
    //**********************************     Charged pions and kaons     **************************************//
    //*********************************************************************************************************//
    TDirectoryFile* directoryChargedPbPb = (TDirectoryFile*)fileDataALICE->Get("ChargedParticles_PbPb_2.76TeV");

    histoChargedPionSpectraStat0005   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraStat0005");
    histoChargedPionSpectraSyst0005   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraSyst0005");
    histoChargedPionSpectraStat0510   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraStat0510");
    histoChargedPionSpectraSyst0510   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraSyst0510");
    histoChargedPionSpectraStat0010   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraStat0010");
    histoChargedPionSpectraSyst0010   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraSyst0010");
    histoChargedPionSpectraStat2040   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraStat2040");
    histoChargedPionSpectraSyst2040   = (TH1D*)directoryChargedPbPb->Get("histoChargedPionSpectraSyst2040");

    histoChargedKaonSpectraStat0005   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraStat0005");
    histoChargedKaonSpectraSyst0005   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraSyst0005");
    histoChargedKaonSpectraStat0510   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraStat0510");
    histoChargedKaonSpectraSyst0510   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraSyst0510");
    histoChargedKaonSpectraStat0010   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraStat0010");
    histoChargedKaonSpectraSyst0010   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraSyst0010");
    histoChargedKaonSpectraStat2040   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraStat2040");
    histoChargedKaonSpectraSyst2040   = (TH1D*)directoryChargedPbPb->Get("histoChargedKaonSpectraSyst2040");

    graphChargedRatioKaonToPion0005      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPion0005");
    graphChargedRatioKaonToPionSys0005   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPionSys0005");
    graphChargedRatioKaonToPion0510      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPion0510");
    graphChargedRatioKaonToPionSys0510   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPionSys0510");
    graphChargedRatioKaonToPion0010      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPion0010");
    graphChargedRatioKaonToPionSys0010   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPionSys0010");
    graphChargedRatioKaonToPion2040      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPion2040");
    graphChargedRatioKaonToPionSys2040   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPionSys2040");

    graphChargedPionRAA0005      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAA0005");
    graphChargedPionRAASys0005   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAASys0005");
    graphChargedPionRAA0510      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAA0510");
    graphChargedPionRAASys0510   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAASys0510");
    graphChargedPionRAA0010      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAA0010");
    graphChargedPionRAASys0010   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAASys0010");
    graphChargedPionRAA2040      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAA2040");
    graphChargedPionRAASys2040   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedPionRAASys2040");

    graphChargedKaonRAA0005      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAA0005");
    graphChargedKaonRAASys0005   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAASys0005");
    graphChargedKaonRAA0510      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAA0510");
    graphChargedKaonRAASys0510   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAASys0510");
    graphChargedKaonRAA0010      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAA0010");
    graphChargedKaonRAASys0010   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAASys0010");
    graphChargedKaonRAA2040      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAA2040");
    graphChargedKaonRAASys2040   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedKaonRAASys2040");


    TGraphAsymmErrors* graphChargedPionRAA[5];
    TGraphAsymmErrors* graphChargedPionRAASys[5];
    TGraphAsymmErrors* graphChargedKaonRAA[5];
    TGraphAsymmErrors* graphChargedKaonRAASys[5];
    TGraphAsymmErrors* graphChargedKaonToPion[5];
    TGraphAsymmErrors* graphChargedKaonToPionSys[5];
    for(Int_t c=0; c<5; c++){
      graphChargedPionRAA[c] = NULL;
      graphChargedPionRAASys[c] = NULL;
      graphChargedKaonRAA[c] = NULL;
      graphChargedKaonRAASys[c] = NULL;
      graphChargedKaonToPion[c] = NULL;
      graphChargedKaonToPionSys[c] = NULL;
    }
    graphChargedPionRAA[0] = (TGraphAsymmErrors*)graphChargedPionRAA0005->Clone("graphChargedPionRAA");
    graphChargedPionRAASys[0] = (TGraphAsymmErrors*)graphChargedPionRAASys0005->Clone("graphChargedPionRAASys");
    graphChargedKaonRAA[0] = (TGraphAsymmErrors*)graphChargedKaonRAA0005->Clone("graphChargedKaonRAA");
    graphChargedKaonRAASys[0] =  (TGraphAsymmErrors*)graphChargedKaonRAASys0005->Clone("graphChargedKaonRAASys");
    graphChargedKaonToPion[0] = (TGraphAsymmErrors*)graphChargedRatioKaonToPion0005->Clone("graphChargedKaonToPion");
    graphChargedKaonToPionSys[0] =  (TGraphAsymmErrors*)graphChargedRatioKaonToPionSys0005->Clone("graphChargedKaonToPionSys");

    graphChargedPionRAA[1] = (TGraphAsymmErrors*)graphChargedPionRAA0510->Clone("graphChargedPionRAA");
    graphChargedPionRAASys[1] = (TGraphAsymmErrors*)graphChargedPionRAASys0510->Clone("graphChargedPionRAASys");
    graphChargedKaonRAA[1] = (TGraphAsymmErrors*)graphChargedKaonRAA0510->Clone("graphChargedKaonRAA");
    graphChargedKaonRAASys[1] =  (TGraphAsymmErrors*)graphChargedKaonRAASys0510->Clone("graphChargedKaonRAASys");
    graphChargedKaonToPion[1] = (TGraphAsymmErrors*)graphChargedRatioKaonToPion0510->Clone("graphChargedKaonToPion");
    graphChargedKaonToPionSys[1] =  (TGraphAsymmErrors*)graphChargedRatioKaonToPionSys0510->Clone("graphChargedKaonToPionSys");

    graphChargedPionRAA[2] = (TGraphAsymmErrors*)graphChargedPionRAA0010->Clone("graphChargedPionRAA");
    graphChargedPionRAASys[2] = (TGraphAsymmErrors*)graphChargedPionRAASys0010->Clone("graphChargedPionRAASys");
    graphChargedKaonRAA[2] = (TGraphAsymmErrors*)graphChargedKaonRAA0010->Clone("graphChargedKaonRAA");
    graphChargedKaonRAASys[2] =  (TGraphAsymmErrors*)graphChargedKaonRAASys0010->Clone("graphChargedKaonRAASys");
    graphChargedKaonToPion[2] = (TGraphAsymmErrors*)graphChargedRatioKaonToPion0010->Clone("graphChargedKaonToPion");
    graphChargedKaonToPionSys[2] =  (TGraphAsymmErrors*)graphChargedRatioKaonToPionSys0010->Clone("graphChargedKaonToPionSys");

    graphChargedPionRAA[3] = (TGraphAsymmErrors*)graphChargedPionRAA2040->Clone("graphChargedPionRAA");
    graphChargedPionRAASys[3] = (TGraphAsymmErrors*)graphChargedPionRAASys2040->Clone("graphChargedPionRAASys");
    graphChargedKaonRAA[3] = (TGraphAsymmErrors*)graphChargedKaonRAA2040->Clone("graphChargedKaonRAA");
    graphChargedKaonRAASys[3] =  (TGraphAsymmErrors*)graphChargedKaonRAASys2040->Clone("graphChargedKaonRAASys");
    graphChargedKaonToPion[3] = (TGraphAsymmErrors*)graphChargedRatioKaonToPion2040->Clone("graphChargedKaonToPion");
    graphChargedKaonToPionSys[3] =  (TGraphAsymmErrors*)graphChargedRatioKaonToPionSys2040->Clone("graphChargedKaonToPionSys");

    graphChargedPionRAA[4] = (TGraphAsymmErrors*)graphChargedPionRAA2040->Clone("graphChargedPionRAA");
    graphChargedPionRAASys[4] = (TGraphAsymmErrors*)graphChargedPionRAASys2040->Clone("graphChargedPionRAASys");
    graphChargedKaonRAA[4] = (TGraphAsymmErrors*)graphChargedKaonRAA2040->Clone("graphChargedKaonRAA");
    graphChargedKaonRAASys[4] =  (TGraphAsymmErrors*)graphChargedKaonRAASys2040->Clone("graphChargedKaonRAASys");
    graphChargedKaonToPion[4] = (TGraphAsymmErrors*)graphChargedRatioKaonToPion2040->Clone("graphChargedKaonToPion");
    graphChargedKaonToPionSys[4] =  (TGraphAsymmErrors*)graphChargedRatioKaonToPionSys2040->Clone("graphChargedKaonToPionSys");


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

        histoPCMPi0InvYieldPbPb2760GeV_0005         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0PCMPbPb2760GeVStatErr_0005");
        graphPCMPi0InvYieldStatPbPb2760GeV_0005     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_0005");
        graphPCMPi0InvYieldSysPbPb2760GeV_0005      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_0005");
        graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0005 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_yShifted_0005");
        graphPCMPi0InvYieldSysWOMat2760GeV_0005     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_forRAA_0005");
        histoPCMPi0RAAStatPbPb2760GeV_0005          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAPi0PCMPbPb2760GeVStatErr_0005");
        graphPCMPi0RAAStatPbPb2760GeV_0005          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVStatErr_0005");
        graphPCMPi0RAASysPbPb2760GeV_0005           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVSysErr_0005");

        histoPCMPi0InvYieldPbPb2760GeV_0510         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0PCMPbPb2760GeVStatErr_0510");
        graphPCMPi0InvYieldStatPbPb2760GeV_0510     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_0510");
        graphPCMPi0InvYieldSysPbPb2760GeV_0510      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_0510");
        graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0510 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_yShifted_0510");
        graphPCMPi0InvYieldSysWOMat2760GeV_0510     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_forRAA_0510");
        histoPCMPi0RAAStatPbPb2760GeV_0510          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAPi0PCMPbPb2760GeVStatErr_0510");
        graphPCMPi0RAAStatPbPb2760GeV_0510          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVStatErr_0510");
        graphPCMPi0RAASysPbPb2760GeV_0510           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVSysErr_0510");

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

        histoPCMPi0InvYieldPbPb2760GeV_2040         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0PCMPbPb2760GeVStatErr_2040");
        graphPCMPi0InvYieldStatPbPb2760GeV_2040     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_2040");
        graphPCMPi0InvYieldSysPbPb2760GeV_2040      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_2040");
        graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2040 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVStatErr_yShifted_2040");
        graphPCMPi0InvYieldSysWOMat2760GeV_2040     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0PCMPbPb2760GeVSysErr_forRAA_2040");
        histoPCMPi0RAAStatPbPb2760GeV_2040          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAPi0PCMPbPb2760GeVStatErr_2040");
        graphPCMPi0RAAStatPbPb2760GeV_2040          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVStatErr_2040");
        graphPCMPi0RAASysPbPb2760GeV_2040           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0PCMPbPb2760GeVSysErr_2040");

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

        histoPCMEtaInvYieldPbPb2760GeV_0005         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaPCMPbPb2760GeVStatErr_0005");
        graphPCMEtaInvYieldStatPbPb2760GeV_0005     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_0005");
        graphPCMEtaInvYieldSysPbPb2760GeV_0005      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_0005");
        graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0005 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_yShifted_0005");
        graphPCMEtaInvYieldSysWOMat2760GeV_0005     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_forRAA_0005");
        histoPCMEtaRAAStatPbPb2760GeV_0005          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAEtaPCMPbPb2760GeVStatErr_0005");
        graphPCMEtaRAAStatPbPb2760GeV_0005          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVStatErr_0005");
        graphPCMEtaRAASysPbPb2760GeV_0005           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVSysErr_0005");
        histoPCMEtatoPi0Stat2760GeV_0005            = (TH1D*)directoryNeutralMesonPbPb->Get("histoEtaToPi0RatioPCMPbPb2760GeVStatErr_0005");
        graphPCMEtatoPi0Stat2760GeV_0005            = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVStatErr_0005");
        graphPCMEtatoPi0Sys2760GeV_0005             = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVSysErr_0005");

        histoPCMEtaInvYieldPbPb2760GeV_0510         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaPCMPbPb2760GeVStatErr_0510");
        graphPCMEtaInvYieldStatPbPb2760GeV_0510     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_0510");
        graphPCMEtaInvYieldSysPbPb2760GeV_0510      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_0510");
        graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0510 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_yShifted_0510");
        graphPCMEtaInvYieldSysWOMat2760GeV_0510     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_forRAA_0510");
        histoPCMEtaRAAStatPbPb2760GeV_0510          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAEtaPCMPbPb2760GeVStatErr_0510");
        graphPCMEtaRAAStatPbPb2760GeV_0510          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVStatErr_0510");
        graphPCMEtaRAASysPbPb2760GeV_0510           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVSysErr_0510");
        histoPCMEtatoPi0Stat2760GeV_0510            = (TH1D*)directoryNeutralMesonPbPb->Get("histoEtaToPi0RatioPCMPbPb2760GeVStatErr_0510");
        graphPCMEtatoPi0Stat2760GeV_0510            = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVStatErr_0510");
        graphPCMEtatoPi0Sys2760GeV_0510             = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVSysErr_0510");

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

        histoPCMEtaInvYieldPbPb2760GeV_2040         = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaPCMPbPb2760GeVStatErr_2040");
        graphPCMEtaInvYieldStatPbPb2760GeV_2040     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_2040");
        graphPCMEtaInvYieldSysPbPb2760GeV_2040      = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_2040");
        graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2040 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVStatErr_yShifted_2040");
        graphPCMEtaInvYieldSysWOMat2760GeV_2040     = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaPCMPbPb2760GeVSysErr_forRAA_2040");
        histoPCMEtaRAAStatPbPb2760GeV_2040          = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAEtaPCMPbPb2760GeVStatErr_2040");
        graphPCMEtaRAAStatPbPb2760GeV_2040          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVStatErr_2040");
        graphPCMEtaRAASysPbPb2760GeV_2040           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaPCMPbPb2760GeVSysErr_2040");
        histoPCMEtatoPi0Stat2760GeV_2040            = (TH1D*)directoryNeutralMesonPbPb->Get("histoEtaToPi0RatioPCMPbPb2760GeVStatErr_2040");
        graphPCMEtatoPi0Stat2760GeV_2040            = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVStatErr_2040");
        graphPCMEtatoPi0Sys2760GeV_2040             = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioPCMPbPb2760GeVSysErr_2040");
        
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

//         histoEMCalPi0InvYieldStatPbPb2760GeV_0010   = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0EMCalPbPb2760GeVStatErr_0010");
//         graphEMCalPi0InvYieldStatPbPb2760GeV_0010   = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVStatErr_0010");
//         graphEMCalPi0InvYieldSysPbPb2760GeV_0010    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErr_0010");
//         graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_0010 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErr_forRAA_0010");
//         histoEMCalPi0RAAStatPbPb2760GeV_0010        = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAPi0EMCalPbPb2760GeVStatErr_0010");
//         graphEMCalPi0RAAStatPbPb2760GeV_0010        = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0EMCalPbPb2760GeVStatErr_0010");
//         graphEMCalPi0RAASysPbPb2760GeV_0010         = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0EMCalPbPb2760GeVSysErr_0010");
// 
//         histoEMCalPi0InvYieldStatPbPb2760GeV_2050   = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldPi0EMCalPbPb2760GeVStatErr_2050");
//         graphEMCalPi0InvYieldStatPbPb2760GeV_2050   = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVStatErr_2050");
//         graphEMCalPi0InvYieldSysPbPb2760GeV_2050    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErr_2050");
//         graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_2050 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErr_forRAA_2050");
//         histoEMCalPi0RAAStatPbPb2760GeV_2050        = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAPi0EMCalPbPb2760GeVStatErr_2050");
//         graphEMCalPi0RAAStatPbPb2760GeV_2050        = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0EMCalPbPb2760GeVStatErr_2050");
//         graphEMCalPi0RAASysPbPb2760GeV_2050         = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAPi0EMCalPbPb2760GeVSysErr_2050");
// 
//         histoEMCalEtaInvYieldStatPbPb2760GeV_0010   = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaEMCalPbPb2760GeVStatErr_0010");
//         graphEMCalEtaInvYieldStatPbPb2760GeV_0010   = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVStatErr_0010");
//         graphEMCalEtaInvYieldSysPbPb2760GeV_0010    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErr_0010");
//         graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErr_forRAA_0010");
//         histoEMCalEtaRAAStatPbPb2760GeV_0010        = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAEtaEMCalPbPb2760GeVStatErr_0010");
//         graphEMCalEtaRAAStatPbPb2760GeV_0010        = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaEMCalPbPb2760GeVStatErr_0010");
//         graphEMCalEtaRAASysPbPb2760GeV_0010         = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaEMCalPbPb2760GeVSysErr_0010");
//         histoEMCalEtatoPi0Stat2760GeV_0010          = (TH1D*)directoryNeutralMesonPbPb->Get("histoEtaToPi0RatioEMCalPbPb2760GeVStatErr_0010");
//         graphEMCalEtatoPi0Stat2760GeV_0010          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioEMCalPbPb2760GeVStatErr_0010");
//         graphEMCalEtatoPi0Sys2760GeV_0010           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioEMCalPbPb2760GeVSysErr_0010");
// 
//         histoEMCalEtaInvYieldStatPbPb2760GeV_2050   = (TH1D*)directoryNeutralMesonPbPb->Get("histoInvYieldEtaEMCalPbPb2760GeVStatErr_2050");
//         graphEMCalEtaInvYieldStatPbPb2760GeV_2050   = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVStatErr_2050");
//         graphEMCalEtaInvYieldSysPbPb2760GeV_2050    = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErr_2050");
//         graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_2050 = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErr_forRAA_2050");
//         histoEMCalEtaRAAStatPbPb2760GeV_2050        = (TH1D*)directoryNeutralMesonPbPb->Get("histoRAAEtaEMCalPbPb2760GeVStatErr_2050");
//         graphEMCalEtaRAAStatPbPb2760GeV_2050        = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaEMCalPbPb2760GeVStatErr_2050");
//         graphEMCalEtaRAASysPbPb2760GeV_2050         = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphRAAEtaEMCalPbPb2760GeVSysErr_2050");
//         histoEMCalEtatoPi0Stat2760GeV_2050          = (TH1D*)directoryNeutralMesonPbPb->Get("histoEtaToPi0RatioEMCalPbPb2760GeVStatErr_2050");
//         graphEMCalEtatoPi0Stat2760GeV_2050          = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioEMCalPbPb2760GeVStatErr_2050");
//         graphEMCalEtatoPi0Sys2760GeV_2050           = (TGraphAsymmErrors*)directoryNeutralMesonPbPb->Get("graphEtaToPi0RatioEMCalPbPb2760GeVSysErr_2050");


    //*********************************************************************************************************//
    //***************************************     generalise arrays     ***************************************//
    //*********************************************************************************************************//

    for(Int_t c=0; c<5; c++){
      histoarrayPCMPi0InvYieldPbPb2760GeV[c] = NULL;
      grapharrayPCMPi0InvYieldStatPbPb2760GeV[c] = NULL;
      grapharrayPCMPi0InvYieldSysPbPb2760GeV[c] = NULL;
      grapharrayPCMPi0InvYieldSysWOMat2760GeV[c] = NULL;

      histoarrayPCMEtaInvYieldPbPb2760GeV[c] = NULL;
      grapharrayPCMEtaInvYieldStatPbPb2760GeV[c] = NULL;
      grapharrayPCMEtaInvYieldSysPbPb2760GeV[c] = NULL;
      grapharrayPCMEtaInvYieldSysWOMat2760GeV[c] = NULL;
      histoarrayPCMEtatoPi0Stat2760GeV[c] = NULL;
      grapharrayPCMEtatoPi0Stat2760GeV[c] = NULL;
      grapharrayPCMEtatoPi0Sys2760GeV[c] = NULL;

      grapharrayPCMInvYieldTotPbPb2760GeV[c] = NULL;
      grapharrayPCMInvYieldStatPbPb2760GeV[c] = NULL;
      grapharrayPCMInvYieldSysPbPb2760GeV[c] = NULL;
      grapharrayPCMInvYieldSysWOMat2760GeV[c] = NULL;

      grapharrayPCMPi0RAAStatPbPb2760GeV[c] = NULL;
      grapharrayPCMPi0RAASysPbPb2760GeV[c] = NULL;
      grapharrayPCMEtaRAAStatPbPb2760GeV[c] = NULL;
      grapharrayPCMEtaRAASysPbPb2760GeV[c] = NULL;
      grapharrayPCMRAAStatPbPb2760GeV[c] = NULL;
      grapharrayPCMRAASysPbPb2760GeV[c] = NULL;

    }

    histoarrayPCMPi0InvYieldPbPb2760GeV[0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_0005->Clone("histoPCMPi0InvYieldPbPb2760GeV_0005");
    grapharrayPCMPi0InvYieldStatPbPb2760GeV[0] = new TGraphAsymmErrors(histoarrayPCMPi0InvYieldPbPb2760GeV[0]);
    grapharrayPCMPi0InvYieldSysPbPb2760GeV[0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0005->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_0005");
    grapharrayPCMPi0InvYieldSysWOMat2760GeV[0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysWOMat2760GeV_0005->Clone("graphPCMPi0InvYieldSysWOMat2760GeV_0005");

    histoarrayPCMEtaInvYieldPbPb2760GeV[0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_0005->Clone("histoPCMPi0InvYieldPbPb2760GeV_0005");
    grapharrayPCMEtaInvYieldStatPbPb2760GeV[0] = new TGraphAsymmErrors(histoarrayPCMEtaInvYieldPbPb2760GeV[0]);
    grapharrayPCMEtaInvYieldSysPbPb2760GeV[0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0005->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_0005");
    grapharrayPCMEtaInvYieldSysWOMat2760GeV[0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysWOMat2760GeV_0005->Clone("graphPCMEtaInvYieldSysWOMat2760GeV_0005");
    histoarrayPCMEtatoPi0Stat2760GeV[0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_0005->Clone("histoPCMEtatoPi0Stat2760GeV_0005");
    grapharrayPCMEtatoPi0Stat2760GeV[0] = new TGraphAsymmErrors(histoarrayPCMEtatoPi0Stat2760GeV[0]);
    grapharrayPCMEtatoPi0Sys2760GeV[0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_0005->Clone("graphPCMEtatoPi0Sys2760GeV_0005");

    grapharrayPCMPi0RAAStatPbPb2760GeV[0] = (TGraphAsymmErrors*)graphPCMPi0RAAStatPbPb2760GeV_0005->Clone("grapharrayPCMPi0RAAStatPbPb2760GeV_0005");
    grapharrayPCMPi0RAASysPbPb2760GeV[0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_0005->Clone("grapharrayPCMPi0RAASysPbPb2760GeV_0005");
    grapharrayPCMEtaRAAStatPbPb2760GeV[0] = (TGraphAsymmErrors*)graphPCMEtaRAAStatPbPb2760GeV_0005->Clone("grapharrayPCMEtaRAAStatPbPb2760GeV_0005");
    grapharrayPCMEtaRAASysPbPb2760GeV[0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_0005->Clone("grapharrayPCMEtaRAASysPbPb2760GeV_0005");

    histoarrayPCMPi0InvYieldPbPb2760GeV[1] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_0510->Clone("histoPCMPi0InvYieldPbPb2760GeV_0510");
    grapharrayPCMPi0InvYieldStatPbPb2760GeV[1] = new TGraphAsymmErrors(histoarrayPCMPi0InvYieldPbPb2760GeV[1]);
    grapharrayPCMPi0InvYieldSysPbPb2760GeV[1] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0510->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_0510");
    grapharrayPCMPi0InvYieldSysWOMat2760GeV[1] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysWOMat2760GeV_0510->Clone("graphPCMPi0InvYieldSysWOMat2760GeV_0510");

    histoarrayPCMEtaInvYieldPbPb2760GeV[1] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_0510->Clone("histoPCMPi0InvYieldPbPb2760GeV_0510");
    grapharrayPCMEtaInvYieldStatPbPb2760GeV[1] = new TGraphAsymmErrors(histoarrayPCMEtaInvYieldPbPb2760GeV[1]);
    grapharrayPCMEtaInvYieldSysPbPb2760GeV[1] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0510->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_0510");
    grapharrayPCMEtaInvYieldSysWOMat2760GeV[1] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysWOMat2760GeV_0510->Clone("graphPCMEtaInvYieldSysWOMat2760GeV_0510");
    histoarrayPCMEtatoPi0Stat2760GeV[1] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_0510->Clone("histoPCMEtatoPi0Stat2760GeV_0510");
    grapharrayPCMEtatoPi0Stat2760GeV[1] = new TGraphAsymmErrors(histoarrayPCMEtatoPi0Stat2760GeV[1]);
    grapharrayPCMEtatoPi0Sys2760GeV[1] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_0510->Clone("graphPCMEtatoPi0Sys2760GeV_0510");

    grapharrayPCMPi0RAAStatPbPb2760GeV[1] = (TGraphAsymmErrors*)graphPCMPi0RAAStatPbPb2760GeV_0510->Clone("grapharrayPCMPi0RAAStatPbPb2760GeV_0510");
    grapharrayPCMPi0RAASysPbPb2760GeV[1] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_0510->Clone("grapharrayPCMPi0RAASysPbPb2760GeV_0510");
    grapharrayPCMEtaRAAStatPbPb2760GeV[1] = (TGraphAsymmErrors*)graphPCMEtaRAAStatPbPb2760GeV_0510->Clone("grapharrayPCMEtaRAAStatPbPb2760GeV_0510");
    grapharrayPCMEtaRAASysPbPb2760GeV[1] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_0510->Clone("grapharrayPCMEtaRAASysPbPb2760GeV_0510");

    histoarrayPCMPi0InvYieldPbPb2760GeV[2] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_0010->Clone("histoPCMPi0InvYieldPbPb2760GeV_0010");
    grapharrayPCMPi0InvYieldStatPbPb2760GeV[2] = new TGraphAsymmErrors(histoarrayPCMPi0InvYieldPbPb2760GeV[2]);
    grapharrayPCMPi0InvYieldSysPbPb2760GeV[2] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0010->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_0010");
    grapharrayPCMPi0InvYieldSysWOMat2760GeV[2] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysWOMat2760GeV_0010->Clone("graphPCMPi0InvYieldSysWOMat2760GeV_0010");

    histoarrayPCMEtaInvYieldPbPb2760GeV[2] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_0010->Clone("histoPCMPi0InvYieldPbPb2760GeV_0010");
    grapharrayPCMEtaInvYieldStatPbPb2760GeV[2] = new TGraphAsymmErrors(histoarrayPCMEtaInvYieldPbPb2760GeV[2]);
    grapharrayPCMEtaInvYieldSysPbPb2760GeV[2] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0010->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_0010");
    grapharrayPCMEtaInvYieldSysWOMat2760GeV[2] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysWOMat2760GeV_0010->Clone("graphPCMEtaInvYieldSysWOMat2760GeV_0010");
    histoarrayPCMEtatoPi0Stat2760GeV[2] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_0010->Clone("histoPCMEtatoPi0Stat2760GeV_0010");
    grapharrayPCMEtatoPi0Stat2760GeV[2] = new TGraphAsymmErrors(histoarrayPCMEtatoPi0Stat2760GeV[2]);
    grapharrayPCMEtatoPi0Sys2760GeV[2] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_0010->Clone("graphPCMEtatoPi0Sys2760GeV_0010");

    grapharrayPCMPi0RAAStatPbPb2760GeV[2] = (TGraphAsymmErrors*)graphPCMPi0RAAStatPbPb2760GeV_0010->Clone("grapharrayPCMPi0RAAStatPbPb2760GeV_0010");
    grapharrayPCMPi0RAASysPbPb2760GeV[2] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_0010->Clone("grapharrayPCMPi0RAASysPbPb2760GeV_0010");
    grapharrayPCMEtaRAAStatPbPb2760GeV[2] = (TGraphAsymmErrors*)graphPCMEtaRAAStatPbPb2760GeV_0010->Clone("grapharrayPCMEtaRAAStatPbPb2760GeV_0010");
    grapharrayPCMEtaRAASysPbPb2760GeV[2] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_0010->Clone("grapharrayPCMEtaRAASysPbPb2760GeV_0010");

    histoarrayPCMPi0InvYieldPbPb2760GeV[3] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_2040->Clone("histoPCMPi0InvYieldPbPb2760GeV_2040");
    grapharrayPCMPi0InvYieldStatPbPb2760GeV[3] = new TGraphAsymmErrors(histoarrayPCMPi0InvYieldPbPb2760GeV[3]);
    grapharrayPCMPi0InvYieldSysPbPb2760GeV[3] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2040->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_2040");
    grapharrayPCMPi0InvYieldSysWOMat2760GeV[3] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysWOMat2760GeV_2040->Clone("graphPCMPi0InvYieldSysWOMat2760GeV_2040");

    histoarrayPCMEtaInvYieldPbPb2760GeV[3] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_2040->Clone("histoPCMPi0InvYieldPbPb2760GeV_2040");
    grapharrayPCMEtaInvYieldStatPbPb2760GeV[3] = new TGraphAsymmErrors(histoarrayPCMEtaInvYieldPbPb2760GeV[3]);
    grapharrayPCMEtaInvYieldSysPbPb2760GeV[3] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2040->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_2040");
    grapharrayPCMEtaInvYieldSysWOMat2760GeV[3] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysWOMat2760GeV_2040->Clone("graphPCMEtaInvYieldSysWOMat2760GeV_2040");
    histoarrayPCMEtatoPi0Stat2760GeV[3] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_2040->Clone("histoPCMEtatoPi0Stat2760GeV_2040");
    grapharrayPCMEtatoPi0Stat2760GeV[3] = new TGraphAsymmErrors(histoarrayPCMEtatoPi0Stat2760GeV[3]);
    grapharrayPCMEtatoPi0Sys2760GeV[3] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_2040->Clone("graphPCMEtatoPi0Sys2760GeV_2040");

    grapharrayPCMPi0RAAStatPbPb2760GeV[3] = (TGraphAsymmErrors*)graphPCMPi0RAAStatPbPb2760GeV_2040->Clone("grapharrayPCMPi0RAAStatPbPb2760GeV_2040");
    grapharrayPCMPi0RAASysPbPb2760GeV[3] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_2040->Clone("grapharrayPCMPi0RAASysPbPb2760GeV_2040");
    grapharrayPCMEtaRAAStatPbPb2760GeV[3] = (TGraphAsymmErrors*)graphPCMEtaRAAStatPbPb2760GeV_2040->Clone("grapharrayPCMEtaRAAStatPbPb2760GeV_2040");
    grapharrayPCMEtaRAASysPbPb2760GeV[3] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_2040->Clone("grapharrayPCMEtaRAASysPbPb2760GeV_2040");

    histoarrayPCMPi0InvYieldPbPb2760GeV[4] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_2050->Clone("histoPCMPi0InvYieldPbPb2760GeV_2050");
    grapharrayPCMPi0InvYieldStatPbPb2760GeV[4] = new TGraphAsymmErrors(histoarrayPCMPi0InvYieldPbPb2760GeV[4]);
    grapharrayPCMPi0InvYieldSysPbPb2760GeV[4] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2050->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_2050");
    grapharrayPCMPi0InvYieldSysWOMat2760GeV[4] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysWOMat2760GeV_2050->Clone("graphPCMPi0InvYieldSysWOMat2760GeV_2050");

    histoarrayPCMEtaInvYieldPbPb2760GeV[4] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_2050->Clone("histoPCMPi0InvYieldPbPb2760GeV_2050");
    grapharrayPCMEtaInvYieldStatPbPb2760GeV[4] = new TGraphAsymmErrors(histoarrayPCMEtaInvYieldPbPb2760GeV[4]);
    grapharrayPCMEtaInvYieldSysPbPb2760GeV[4] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2050->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_2050");
    grapharrayPCMEtaInvYieldSysWOMat2760GeV[4] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysWOMat2760GeV_2050->Clone("graphPCMEtaInvYieldSysWOMat2760GeV_2050");
    histoarrayPCMEtatoPi0Stat2760GeV[4] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_2050->Clone("histoPCMEtatoPi0Stat2760GeV_2050");
    grapharrayPCMEtatoPi0Stat2760GeV[4] = new TGraphAsymmErrors(histoarrayPCMEtatoPi0Stat2760GeV[4]);
    grapharrayPCMEtatoPi0Sys2760GeV[4] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_2050->Clone("graphPCMEtatoPi0Sys2760GeV_2050");

    grapharrayPCMPi0RAAStatPbPb2760GeV[4] = (TGraphAsymmErrors*)graphPCMPi0RAAStatPbPb2760GeV_2050->Clone("grapharrayPCMPi0RAAStatPbPb2760GeV_2050");
    grapharrayPCMPi0RAASysPbPb2760GeV[4] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_2050->Clone("grapharrayPCMPi0RAASysPbPb2760GeV_2050");
    grapharrayPCMEtaRAAStatPbPb2760GeV[4] = (TGraphAsymmErrors*)graphPCMEtaRAAStatPbPb2760GeV_2050->Clone("grapharrayPCMEtaRAAStatPbPb2760GeV_2050");
    grapharrayPCMEtaRAASysPbPb2760GeV[4] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_2050->Clone("grapharrayPCMEtaRAASysPbPb2760GeV_2050");

    if(meson.CompareTo("Pi0")==0){

      grapharrayPCMInvYieldStatPbPb2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldStatPbPb2760GeV[0]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[0]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysWOMat2760GeV[0]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[0] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[0],grapharrayPCMInvYieldSysPbPb2760GeV[0]);

      grapharrayPCMRAAStatPbPb2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMPi0RAAStatPbPb2760GeV[0]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMPi0RAASysPbPb2760GeV[0]->Clone("grapharrayPCMRAASysPbPb2760GeV");

      grapharrayPCMInvYieldStatPbPb2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldStatPbPb2760GeV[1]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[1]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysWOMat2760GeV[1]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[1] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[1],grapharrayPCMInvYieldSysPbPb2760GeV[1]);

      grapharrayPCMRAAStatPbPb2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMPi0RAAStatPbPb2760GeV[1]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMPi0RAASysPbPb2760GeV[1]->Clone("grapharrayPCMRAASysPbPb2760GeV");

      grapharrayPCMInvYieldStatPbPb2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldStatPbPb2760GeV[2]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[2]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysWOMat2760GeV[2]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[2] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[2],grapharrayPCMInvYieldSysPbPb2760GeV[2]);

      grapharrayPCMRAAStatPbPb2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMPi0RAAStatPbPb2760GeV[2]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMPi0RAASysPbPb2760GeV[2]->Clone("grapharrayPCMRAASysPbPb2760GeV");

      grapharrayPCMInvYieldStatPbPb2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldStatPbPb2760GeV[3]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[3]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysWOMat2760GeV[3]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[3] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[3],grapharrayPCMInvYieldSysPbPb2760GeV[3]);

      grapharrayPCMRAAStatPbPb2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMPi0RAAStatPbPb2760GeV[3]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMPi0RAASysPbPb2760GeV[3]->Clone("grapharrayPCMRAASysPbPb2760GeV");

      grapharrayPCMInvYieldStatPbPb2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldStatPbPb2760GeV[4]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[4]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysWOMat2760GeV[4]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[4] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[4],grapharrayPCMInvYieldSysPbPb2760GeV[4]);

      grapharrayPCMRAAStatPbPb2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMPi0RAAStatPbPb2760GeV[4]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMPi0RAASysPbPb2760GeV[4]->Clone("grapharrayPCMRAASysPbPb2760GeV");

    } else if(meson.CompareTo("Eta")==0){

      grapharrayPCMInvYieldStatPbPb2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldStatPbPb2760GeV[0]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[0]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysWOMat2760GeV[0]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[0] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[0],grapharrayPCMInvYieldSysPbPb2760GeV[0]);

      grapharrayPCMRAAStatPbPb2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMEtaRAAStatPbPb2760GeV[0]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[0] = (TGraphAsymmErrors*)grapharrayPCMEtaRAASysPbPb2760GeV[0]->Clone("grapharrayPCMRAASysPbPb2760GeV");

      grapharrayPCMInvYieldStatPbPb2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldStatPbPb2760GeV[1]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[1]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysWOMat2760GeV[1]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[1] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[1],grapharrayPCMInvYieldSysPbPb2760GeV[1]);

      grapharrayPCMRAAStatPbPb2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMEtaRAAStatPbPb2760GeV[1]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[1] = (TGraphAsymmErrors*)grapharrayPCMEtaRAASysPbPb2760GeV[1]->Clone("grapharrayPCMRAASysPbPb2760GeV");

      grapharrayPCMInvYieldStatPbPb2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldStatPbPb2760GeV[2]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[2]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysWOMat2760GeV[2]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[2] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[2],grapharrayPCMInvYieldSysPbPb2760GeV[2]);

      grapharrayPCMRAAStatPbPb2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMEtaRAAStatPbPb2760GeV[2]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[2] = (TGraphAsymmErrors*)grapharrayPCMEtaRAASysPbPb2760GeV[2]->Clone("grapharrayPCMRAASysPbPb2760GeV");

      grapharrayPCMInvYieldStatPbPb2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldStatPbPb2760GeV[3]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[3]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysWOMat2760GeV[3]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[3] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[3],grapharrayPCMInvYieldSysPbPb2760GeV[3]);

      grapharrayPCMRAAStatPbPb2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMEtaRAAStatPbPb2760GeV[3]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[3] = (TGraphAsymmErrors*)grapharrayPCMEtaRAASysPbPb2760GeV[3]->Clone("grapharrayPCMRAASysPbPb2760GeV");

      grapharrayPCMInvYieldStatPbPb2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldStatPbPb2760GeV[4]->Clone("grapharrayPCMInvYieldStatPbPb2760GeV");
      grapharrayPCMInvYieldSysPbPb2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[4]->Clone("grapharrayPCMInvYieldSysPbPb2760GeV");
      grapharrayPCMInvYieldSysWOMat2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysWOMat2760GeV[4]->Clone("grapharrayPCMInvYieldSysWOMat2760GeV");
      grapharrayPCMInvYieldTotPbPb2760GeV[4] = CalculateCombinedSysAndStatError(grapharrayPCMInvYieldStatPbPb2760GeV[4],grapharrayPCMInvYieldSysPbPb2760GeV[4]);

      grapharrayPCMRAAStatPbPb2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMEtaRAAStatPbPb2760GeV[4]->Clone("grapharrayPCMRAAStatPbPb2760GeV");
      grapharrayPCMRAASysPbPb2760GeV[4] = (TGraphAsymmErrors*)grapharrayPCMEtaRAASysPbPb2760GeV[4]->Clone("grapharrayPCMRAASysPbPb2760GeV");

    }

    //syst and stat relative errors
    TH1D *statErrorCollectionLHC11h[5][11];
    TH1D *statErrorCollectionEtatoPi0LHC11h[5][11];
    TH1D *statErrorCollectionRaaLHC11h[5][11];
    TGraphAsymmErrors *sysErrorCollectionLHC11h[5][11];
    TGraphAsymmErrors *sysErrorCollectionEtatoPi0LHC11h[5][11];
    TGraphAsymmErrors *sysErrorCollectionRaaLHC11h[5][11];

    for (Int_t i = 0; i< 11; i++){
      for (Int_t c = 0; c< 5; c++){
        statErrorCollectionLHC11h[c][i] = NULL;
        sysErrorCollectionLHC11h[c][i] = NULL;
        statErrorCollectionEtatoPi0LHC11h[c][i] = NULL;
        sysErrorCollectionEtatoPi0LHC11h[c][i] = NULL;
        statErrorCollectionRaaLHC11h[c][i] = NULL;
        sysErrorCollectionRaaLHC11h[c][i] = NULL;
      }
    }
    if(meson.CompareTo("Pi0")==0){
      statErrorCollectionLHC11h[0][0] = (TH1D*)histoarrayPCMPi0InvYieldPbPb2760GeV[0]->Clone("statErrPCMPi0_0005");
      statErrorCollectionLHC11h[1][0] = (TH1D*)histoarrayPCMPi0InvYieldPbPb2760GeV[1]->Clone("statErrPCMPi0_0510");
      statErrorCollectionLHC11h[2][0] = (TH1D*)histoarrayPCMPi0InvYieldPbPb2760GeV[2]->Clone("statErrPCMPi0_0010");
      statErrorCollectionLHC11h[3][0] = (TH1D*)histoarrayPCMPi0InvYieldPbPb2760GeV[3]->Clone("statErrPCMPi0_2040");
      statErrorCollectionLHC11h[4][0] = (TH1D*)histoarrayPCMPi0InvYieldPbPb2760GeV[4]->Clone("statErrPCMPi0_2050");

      sysErrorCollectionLHC11h[0][0] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[0]->Clone("sysErrPCMPi0_0005");
      sysErrorCollectionLHC11h[1][0] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[1]->Clone("sysErrPCMPi0_0510");
      sysErrorCollectionLHC11h[2][0] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[2]->Clone("sysErrPCMPi0_0010");
      sysErrorCollectionLHC11h[3][0] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[3]->Clone("sysErrPCMPi0_2040");
      sysErrorCollectionLHC11h[4][0] = (TGraphAsymmErrors*)grapharrayPCMPi0InvYieldSysPbPb2760GeV[4]->Clone("sysErrPCMPi0_2050");

//       statErrorCollectionRaaLHC11h[0][0] = (TH1D*)histoarrayPCMPi0RAAStatPbPb2760GeV[0]->Clone("statErrPCMPi0RAA_0005");
//       statErrorCollectionRaaLHC11h[1][0] = (TH1D*)histoarrayPCMPi0RAAStatPbPb2760GeV[1]->Clone("statErrPCMPi0RAA_0510");
//       statErrorCollectionRaaLHC11h[2][0] = (TH1D*)histoarrayPCMPi0RAAStatPbPb2760GeV[2]->Clone("statErrPCMPi0RAA_0010");
//       statErrorCollectionRaaLHC11h[3][0] = (TH1D*)histoarrayPCMPi0RAAStatPbPb2760GeV[3]->Clone("statErrPCMPi0RAA_2040");
//       statErrorCollectionRaaLHC11h[4][0] = (TH1D*)histoarrayPCMPi0RAAStatPbPb2760GeV[4]->Clone("statErrPCMPi0RAA_2050");

//       sysErrorCollectionRaaLHC11h[0][0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV[0]->Clone("sysErrPCMPi0Raa_0005");
//       sysErrorCollectionRaaLHC11h[1][0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV[1]->Clone("sysErrPCMPi0Raa_0510");
//       sysErrorCollectionRaaLHC11h[2][0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV[2]->Clone("sysErrPCMPi0Raa_0010");
//       sysErrorCollectionRaaLHC11h[3][0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV[3]->Clone("sysErrPCMPi0Raa_2040");
//       sysErrorCollectionRaaLHC11h[4][0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV[4]->Clone("sysErrPCMPi0Raa_2050");

    } else if(meson.CompareTo("Eta")==0) {

      statErrorCollectionLHC11h[0][0] = (TH1D*)histoarrayPCMEtaInvYieldPbPb2760GeV[0]->Clone("statErrPCMEta_0005");
      statErrorCollectionLHC11h[1][0] = (TH1D*)histoarrayPCMEtaInvYieldPbPb2760GeV[1]->Clone("statErrPCMEta_0510");
      statErrorCollectionLHC11h[2][0] = (TH1D*)histoarrayPCMEtaInvYieldPbPb2760GeV[2]->Clone("statErrPCMEta_0010");
      statErrorCollectionLHC11h[3][0] = (TH1D*)histoarrayPCMEtaInvYieldPbPb2760GeV[3]->Clone("statErrPCMEta_2040");
      statErrorCollectionLHC11h[4][0] = (TH1D*)histoarrayPCMEtaInvYieldPbPb2760GeV[4]->Clone("statErrPCMEta_2050");

      sysErrorCollectionLHC11h[0][0] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[0]->Clone("sysErrPCMEta_0005");
      sysErrorCollectionLHC11h[1][0] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[1]->Clone("sysErrPCMEta_0510");
      sysErrorCollectionLHC11h[2][0] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[2]->Clone("sysErrPCMEta_0010");
      sysErrorCollectionLHC11h[3][0] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[3]->Clone("sysErrPCMEta_2040");
      sysErrorCollectionLHC11h[4][0] = (TGraphAsymmErrors*)grapharrayPCMEtaInvYieldSysPbPb2760GeV[4]->Clone("sysErrPCMEta_2050");

//       statErrorCollectionRaaLHC11h[0][0] = (TH1D*)histoarrayPCMEtaRAAStatPbPb2760GeV[0]->Clone("statErrPCMEtaRAA_0005");
//       statErrorCollectionRaaLHC11h[1][0] = (TH1D*)histoarrayPCMEtaRAAStatPbPb2760GeV[1]->Clone("statErrPCMEtaRAA_0510");
//       statErrorCollectionRaaLHC11h[2][0] = (TH1D*)histoarrayPCMEtaRAAStatPbPb2760GeV[2]->Clone("statErrPCMEtaRAA_0010");
//       statErrorCollectionRaaLHC11h[3][0] = (TH1D*)histoarrayPCMEtaRAAStatPbPb2760GeV[3]->Clone("statErrPCMEtaRAA_2040");
//       statErrorCollectionRaaLHC11h[4][0] = (TH1D*)histoarrayPCMEtaRAAStatPbPb2760GeV[4]->Clone("statErrPCMEtaRAA_2050");

//       sysErrorCollectionRaaLHC11h[0][0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV[0]->Clone("sysErrPCMEtaRaa_0005");
//       sysErrorCollectionRaaLHC11h[1][0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV[1]->Clone("sysErrPCMEtaRaa_0510");
//       sysErrorCollectionRaaLHC11h[2][0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV[2]->Clone("sysErrPCMEtaRaa_0010");
//       sysErrorCollectionRaaLHC11h[3][0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV[3]->Clone("sysErrPCMEtaRaa_2040");
//       sysErrorCollectionRaaLHC11h[4][0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV[4]->Clone("sysErrPCMEtaRaa_2050");

      statErrorCollectionEtatoPi0LHC11h[0][0] = (TH1D*)histoarrayPCMEtatoPi0Stat2760GeV[0]->Clone("statErrPCMEtatoPi0_0005");
      statErrorCollectionEtatoPi0LHC11h[1][0] = (TH1D*)histoarrayPCMEtatoPi0Stat2760GeV[1]->Clone("statErrPCMEtatoPi0_0510");
      statErrorCollectionEtatoPi0LHC11h[2][0] = (TH1D*)histoarrayPCMEtatoPi0Stat2760GeV[2]->Clone("statErrPCMEtatoPi0_0010");
      statErrorCollectionEtatoPi0LHC11h[3][0] = (TH1D*)histoarrayPCMEtatoPi0Stat2760GeV[3]->Clone("statErrPCMEtatoPi0_2040");
      statErrorCollectionEtatoPi0LHC11h[4][0] = (TH1D*)histoarrayPCMEtatoPi0Stat2760GeV[4]->Clone("statErrPCMEtatoPi0_2050");

      sysErrorCollectionEtatoPi0LHC11h[0][0] = (TGraphAsymmErrors*)grapharrayPCMEtatoPi0Sys2760GeV[0]->Clone("sysErrPCMEtatoPi0_0005");
      sysErrorCollectionEtatoPi0LHC11h[1][0] = (TGraphAsymmErrors*)grapharrayPCMEtatoPi0Sys2760GeV[1]->Clone("sysErrPCMEtatoPi0_0510");
      sysErrorCollectionEtatoPi0LHC11h[2][0] = (TGraphAsymmErrors*)grapharrayPCMEtatoPi0Sys2760GeV[2]->Clone("sysErrPCMEtatoPi0_0010");
      sysErrorCollectionEtatoPi0LHC11h[3][0] = (TGraphAsymmErrors*)grapharrayPCMEtatoPi0Sys2760GeV[3]->Clone("sysErrPCMEtatoPi0_2040");
      sysErrorCollectionEtatoPi0LHC11h[4][0] = (TGraphAsymmErrors*)grapharrayPCMEtatoPi0Sys2760GeV[4]->Clone("sysErrPCMEtatoPi0_2050");

    }


    //*********************************************************************************************************************//
    //******************************************** Visualize relative errors **********************************************//
    //*********************************************************************************************************************//
    Int_t nMeasSetLHC11h = 1;
    //relative errors
    TH1D* statErrorRelCollectionLHC11h[5][11];
    TGraphAsymmErrors* sysErrorRelCollectionLHC11h[5][11];
    for (Int_t i = 0; i< 11; i++){
      for (Int_t c = 0; c< 5; c++){
        statErrorRelCollectionLHC11h[c][i] = NULL;
        sysErrorRelCollectionLHC11h[c][i] = NULL;
      }
    }
    for (Int_t i = 0; i < 11; i++){
      for (Int_t c = 0; c< 5; c++){
        if (statErrorCollectionLHC11h[c][i])statErrorRelCollectionLHC11h[c][i] = CalculateRelErrUpTH1D(statErrorCollectionLHC11h[c][i], Form("relativeStatErrorLHC11h_%s_%d", nameMeasGlobal[i].Data(),c));
        if (sysErrorCollectionLHC11h[c][i]) sysErrorRelCollectionLHC11h[c][i] = CalculateRelErrUpAsymmGraph(sysErrorCollectionLHC11h[c][i], Form("relativeSysErrorLHC11h_%s_%d", nameMeasGlobal[i].Data(),c));
      }
    }

    Int_t textSizeLabelsPixel = 900*0.04;
    TCanvas* canvasRelSysErrLHC11h = new TCanvas("canvasRelSysErrLHC11h","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErrLHC11h, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErrLHC11h->SetLogx();

      TH2F * histo2DRelSysErrLHC11h = new TH2F("histo2DRelSysErrLHC11h","histo2DRelSysErrLHC11h",11000,0.23,30.,1000,0,80.5);
      SetStyleHistoTH2ForGraphs(histo2DRelSysErrLHC11h, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
      histo2DRelSysErrLHC11h->GetXaxis()->SetMoreLogLabels();
      histo2DRelSysErrLHC11h->GetXaxis()->SetLabelOffset(-0.01);
      histo2DRelSysErrLHC11h->Draw("copy");

//       if(meson.CompareTo("Pi0")==0){
//         for(Int_t a=0; a<3;a++){
//           for (Int_t c = 0; c< 5; c++){
//             sysErrorRelCollectionLHC11h[c][0]->RemovePoint(0);
//           }
//         }
//       }
      TLegend* legendRelSysErrLHC11h = GetAndSetLegend2(0.12, 0.94-(0.04*5*1.35), 0.5, 0.94, 32);
      for (Int_t i = 0; i < nMeasSetLHC11h; i++){
        for (Int_t c = 0; c< 5; c++){
          DrawGammaSetMarkerTGraph(sysErrorRelCollectionLHC11h[c][nMeasSetLHC11h-1], markerStyleDet[nMeasSetLHC11h-1],
                                      markerSizeDet[nMeasSetLHC11h-1]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
          sysErrorRelCollectionLHC11h[c][nMeasSetLHC11h-1]->Draw("p,same,e0");
          legendRelSysErrLHC11h->AddEntry(sysErrorRelCollectionLHC11h[c][nMeasSetLHC11h-1],cent[c],"p");

        }
      }
      legendRelSysErrLHC11h->Draw();

      TLatex *labelRelSysErrEnergy = new TLatex(0.6,0.89,collisionSystem2760GeV.Data());
      SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel,4);
      labelRelSysErrEnergy->SetTextFont(43);

      TLatex *labelRelSysErrPi0;
      if(meson.CompareTo("Pi0")==0){
          labelRelSysErrPi0= new TLatex(0.6,0.85,"#pi^{0} #rightarrow #gamma#gamma");
      } else if(meson.CompareTo("Eta")==0){
          labelRelSysErrPi0= new TLatex(0.6,0.85,"#eta #rightarrow #gamma#gamma");
      }
      SetStyleTLatex( labelRelSysErrPi0, 0.035,4);
      labelRelSysErrPi0->SetTextFont(43);
      labelRelSysErrEnergy->Draw();
      labelRelSysErrPi0->Draw();

//     canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    TCanvas* canvasRelStatErrLHC11h = new TCanvas("canvasRelStatErrLHC11h","",200,10,1350,900);  // gives the page size
      DrawGammaCanvasSettings( canvasRelStatErrLHC11h, 0.08, 0.02, 0.035, 0.09);
      canvasRelStatErrLHC11h->SetLogx();

      TH2F * histo2DRelStatErrLHC11h = new TH2F("histo2DRelStatErrLHC11h","histo2DRelStatErrLHC11h",11000,0.23,30.,1000,0,80.5);
      SetStyleHistoTH2ForGraphs(histo2DRelStatErrLHC11h, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
      histo2DRelStatErrLHC11h->GetXaxis()->SetMoreLogLabels();
      histo2DRelStatErrLHC11h->GetXaxis()->SetLabelOffset(-0.01);
      histo2DRelStatErrLHC11h->Draw("copy");

      TLegend* legendRelStatErrLHC11h = GetAndSetLegend2(0.32, 0.94-(0.04*5*1.35), 0.5, 0.94, 32);
      for (Int_t i = 0; i < nMeasSetLHC11h; i++){
        for (Int_t c = 0; c< 5; c++){
          DrawGammaSetMarker(statErrorRelCollectionLHC11h[c][nMeasSetLHC11h-1], markerStyleDet[nMeasSetLHC11h-1],
                                      markerSizeDet[nMeasSetLHC11h-1]*0.5,
                                      GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
          statErrorRelCollectionLHC11h[c][nMeasSetLHC11h-1]->Draw("p,same,e0");
          legendRelStatErrLHC11h->AddEntry(statErrorRelCollectionLHC11h[c][nMeasSetLHC11h-1],cent[c],"p");

        }
      }
      legendRelStatErrLHC11h->Draw();

//       TLatex *labelRelStatErrEnergy = new TLatex(0.6,0.89,collisionSystem2760GeV.Data());
//       SetStyleTLatex( labelRelStatErrEnergy,textSizeLabelsPixel, 4);
//       labelRelStatErrEnergy->SetTextFont(43);
//       labelRelStatErrEnergy->Draw();

//       if(meson.CompareTo("Pi0")==0){
//           labelRelStatErrPi0= new TLatex(0.6,0.85,"#pi^{0} #rightarrow #gamma#gamma");
//       } else if(meson.CompareTo("Eta")==0){
//           labelRelStatErrPi0= new TLatex(0.6,0.85,"#eta #rightarrow #gamma#gamma");
//       }
//       SetStyleTLatex( labelRelStatErrPi0,textSizeLabelsPixel, 4);
//       labelRelStatErrPi0->SetTextFont(43);
//       labelRelStatErrPi0->Draw();

      labelRelSysErrEnergy->Draw();
      labelRelSysErrPi0->Draw();

//     canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    //  *********************************************************************************************************************
    //  ************************ Visualize relative total errors of different combination methods Pi0 ***********************
    //  *********************************************************************************************************************
    TGraphAsymmErrors* graphInvYieldsRelStat[5];
    TGraphAsymmErrors* graphInvYieldsRelSys[5];
    TGraphAsymmErrors* graphInvYieldsRelTot[5];
    for (Int_t c = 0; c< 5; c++){
      graphInvYieldsRelStat[c] = CalculateRelErrUpAsymmGraph( grapharrayPCMInvYieldStatPbPb2760GeV[c], Form("relativeStatError_%d",c));
      graphInvYieldsRelSys[c] = CalculateRelErrUpAsymmGraph( grapharrayPCMInvYieldSysPbPb2760GeV[c], Form("relativeSysError_%d",c));
      graphInvYieldsRelTot[c] = CalculateRelErrUpAsymmGraph( grapharrayPCMInvYieldTotPbPb2760GeV[c], Form("relativeTotalError_%d",c));
    }

    TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();

    TH2F * histo2DRelTotErr = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,0.3,30.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelTotErr->GetYaxis()->SetRangeUser(0,50.5);

      labelRelSysErrEnergy->Draw();
      labelRelSysErrPi0->Draw();

    for(Int_t c = 0; c< 5; c++){
      canvasRelTotErr->cd();
      histo2DRelTotErr->Draw("copy");


//       if(meson.CompareTo("Pi0")==0){
//         for(Int_t a=0; a<3;a++){
//           graphInvYieldsRelStat[c]->RemovePoint(0);
//           graphInvYieldsRelSys[c]->RemovePoint(0);
//           graphInvYieldsRelTot[c]->RemovePoint(0);
//         }
//       }

      DrawGammaSetMarkerTGraphAsym(graphInvYieldsRelTot[c], markerStyleComb, markerSizeComb, colorComb , colorComb);
      graphInvYieldsRelTot[c]->Draw("p,same,z");
      DrawGammaSetMarkerTGraphAsym(graphInvYieldsRelStat[c], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
      graphInvYieldsRelStat[c]->Draw("l,x0,same,e1");
      DrawGammaSetMarkerTGraphAsym(graphInvYieldsRelSys[c], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
      graphInvYieldsRelSys[c]->SetLineStyle(7);
      graphInvYieldsRelSys[c]->Draw("l,x0,same,e1");

      TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
      legendRelTotErr3->AddEntry(graphInvYieldsRelTot[c],"tot","p");
      legendRelTotErr3->AddEntry(graphInvYieldsRelStat[c],"stat","l");
      legendRelTotErr3->AddEntry(graphInvYieldsRelSys[c],"sys","l");
      legendRelTotErr3->Draw();

      labelRelSysErrEnergy->Draw();
      labelRelSysErrPi0->Draw();

//       canvasRelTotErr->SaveAs(Form("%s/%s_RelErrdecomp_%d.%s",outputDir.Data(),meson.Data(),c,suffix.Data()));
    }


    //**********************************************************************************************************************//
    //************************************* Calculating bin shifted spectra & fitting **************************************//
    //**********************************************************************************************************************//
    //cloning spectra for shifting
    TGraphAsymmErrors* grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[5];
    TGraphAsymmErrors* grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[5];
    TGraphAsymmErrors* grapharrayPCMInvYieldSysPbPb2760GeVUnshifted[5];
    TGraphAsymmErrors* grapharrayPCMInvYieldSysWOMat2760GeVUnshifted[5];
    TGraphAsymmErrors* grapharrayPCMInvYieldStatPbPb2760GeVYshifted[5];
    TGraphAsymmErrors* grapharrayPCMInvYieldSysPbPb2760GeVYshifted[5];
    TGraphAsymmErrors* grapharrayRatioToFitStatPbPb2760GeV[5];
    TGraphAsymmErrors* grapharrayRatioToFitSysPbPb2760GeV[5];
    for(Int_t c=0; c<5; c++){
      grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[c] = NULL;
      grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[c] = NULL;
      grapharrayPCMInvYieldSysPbPb2760GeVUnshifted[c] = NULL;
      grapharrayPCMInvYieldSysWOMat2760GeVUnshifted[c] = NULL;
      grapharrayPCMInvYieldStatPbPb2760GeVYshifted[c] = NULL;
      grapharrayPCMInvYieldSysPbPb2760GeVYshifted[c] = NULL;
      grapharrayRatioToFitStatPbPb2760GeV[c] = NULL;
      grapharrayRatioToFitSysPbPb2760GeV[c] = NULL;
    }

    for(Int_t c = 0; c< 5; c++){
      grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldTotPbPb2760GeV[c]->Clone(Form("UnshiftedSpectrumStat_%d",c));
      grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldStatPbPb2760GeV[c]->Clone(Form("UnshiftedSpectrumStat_%d",c));
      grapharrayPCMInvYieldSysPbPb2760GeVUnshifted[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldSysPbPb2760GeV[c]->Clone(Form("UnshiftedSpectrumSys_%d",c));
      grapharrayPCMInvYieldSysWOMat2760GeVUnshifted[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldSysWOMat2760GeV[c]->Clone(Form("UnshiftedSpectrumSysWOM_%d",c));
    }

    Double_t fitStart;
    Double_t fitStop;
    if(meson.CompareTo("Pi0")==0){
      fitStart = 0.4;
      fitStop = 14.;
    } else if(meson.CompareTo("Eta")==0){
      fitStart = 1.;
      fitStop = 10.;
    }
    TF1 *fitBylinkinUnshifted;
    TF1 *fitTsallisPbPb2760GeVPtLHC11h;
    TF1 *fitQCDPbPb2760GeVPtLHC11h;
    for(Int_t c = 0; c< 5; c++){

      Double_t paramTCM[5] = {grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[c]->GetY()[0],0.3,grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[c]->GetY()[0],0.3,8};
      //Bylinkin
      fitBylinkinUnshifted = FitObject("tcm","tcm",meson.Data(),grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[c],fitStart,fitStop,paramTCM);
      grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[c]->Fit(fitBylinkinUnshifted,"NRMEX0+","",fitStart,fitStop);
      cout << WriteParameterToFile(fitBylinkinUnshifted) << endl << endl;

      Double_t paramQCD[5] = {24,5.,-20.,2.,20};
      fitQCDPbPb2760GeVPtLHC11h = FitObject("qcd","fitQCD",meson.Data(),grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[c],fitStart,fitStop,paramQCD);
      grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[c]->Fit(fitQCDPbPb2760GeVPtLHC11h,"NRMEX0+","",fitStart,fitStop);
      cout << WriteParameterToFile(fitQCDPbPb2760GeVPtLHC11h) << endl << endl;

      //Calculating binshifts
      TF1* fitFunctionShiftingX;
      if(bWCorrection.CompareTo("X")==0 ){

        //shifting in X
        fitFunctionShiftingX = FitObject("tcmpt","tcmpt",meson.Data(),grapharrayPCMInvYieldTotPbPb2760GeV[c]);
        grapharrayPCMInvYieldTotPbPb2760GeV[c] = ApplyXshift(grapharrayPCMInvYieldTotPbPb2760GeV[c], fitFunctionShiftingX,meson.Data());
        fitFunctionShiftingX = FitObject("tcmpt","tcmpt",meson.Data(),grapharrayPCMInvYieldStatPbPb2760GeV[c]);
        grapharrayPCMInvYieldStatPbPb2760GeV[c] = ApplyXshift(grapharrayPCMInvYieldStatPbPb2760GeV[c], fitFunctionShiftingX,meson.Data());
        fitFunctionShiftingX = FitObject("tcmpt","tcmpt",meson.Data(),grapharrayPCMInvYieldSysPbPb2760GeV[c]);
        grapharrayPCMInvYieldSysPbPb2760GeV[c] = ApplyXshift(grapharrayPCMInvYieldSysPbPb2760GeV[c], fitFunctionShiftingX,meson.Data());

        TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1350,1350*1.15);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2, 0.16, 0.02, 0.02, 0.09);
        canvasDummy2->SetLogx();
        canvasDummy2->SetLogy();

        TH2F * histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",11000,minPtRange-0.2,maxPtRange,1000,minYaxisYields,maxYaxisYields);
            SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
        histo2DDummy2->GetXaxis()->SetMoreLogLabels();
        histo2DDummy2->GetXaxis()->SetLabelOffset(-0.01);
        histo2DDummy2->Draw("copy");

          fitBylinkinUnshifted->SetLineColor(kBlue);
          fitBylinkinUnshifted->Draw("same");
          fitQCDPbPb2760GeVPtLHC11h->SetLineColor(kGreen+2);
          fitQCDPbPb2760GeVPtLHC11h->Draw("same");

          DrawGammaSetMarkerTGraphAsym(grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[c], 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
          grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[c]->Draw("pEsame");
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMInvYieldStatPbPb2760GeV[c], 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
          grapharrayPCMInvYieldStatPbPb2760GeV[c]->Draw("pEsame");

          TLegend* legendXdummy = new TLegend(0.55,0.85,0.85,0.95);
          legendXdummy->SetFillColor(0);
          legendXdummy->SetLineColor(0);
          legendXdummy->SetTextFont(42);
          legendXdummy->SetTextSize(FontSize);
          legendXdummy->AddEntry(grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[c],"combined unshifted","lp");
          legendXdummy->AddEntry(grapharrayPCMInvYieldStatPbPb2760GeV[c],"combined shifted","lp");
          legendXdummy->AddEntry(fitBylinkinUnshifted,"Bylinkin","l");
          legendXdummy->AddEntry(fitQCDPbPb2760GeVPtLHC11h,"QCD","l");
          legendXdummy->Draw();

        canvasDummy2->Update();
//         canvasDummy2->Print(Form("%s/%s_ComparisonXShifted_%d.%s",outputDir.Data(),meson.Data(),c,suffix.Data()));

        //shifting in Y
        TF1* fitFunctionShiftingY = (TF1*)fitFunctionShiftingX->Clone("fitFunctionShiftingY");

        cout << "combined binshift Y" << endl;
        grapharrayPCMInvYieldSysPbPb2760GeVYshifted[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldSysWOMat2760GeVUnshifted[c]->Clone(Form("YShiftedCombSys%d",c));
        grapharrayPCMInvYieldSysPbPb2760GeVYshifted[c] = ApplyYshiftIndividualSpectra(grapharrayPCMInvYieldSysPbPb2760GeVYshifted[c],fitFunctionShiftingY);

        grapharrayPCMInvYieldStatPbPb2760GeVYshifted[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[c]->Clone(Form("YShiftedCombStat%d",c));
        grapharrayPCMInvYieldStatPbPb2760GeVYshifted[c] = ApplyYshiftIndividualSpectra(grapharrayPCMInvYieldStatPbPb2760GeVYshifted[c], fitFunctionShiftingY);

        TCanvas* canvasDummy3 = new TCanvas("canvasDummy3","",200,10,1350,1350*1.15);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy3, 0.16, 0.02, 0.02, 0.09);
        canvasDummy3->SetLogx();
        canvasDummy3->SetLogy();

        TH2F * histo2DDummy3 = new TH2F("histo2DDummy3","histo2DDummy3",11000,minPtRange-0.2,maxPtRange,1000,minYaxisYields,maxYaxisYields);
            SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
        histo2DDummy3->GetXaxis()->SetMoreLogLabels();
        histo2DDummy3->GetXaxis()->SetLabelOffset(-0.01);
        histo2DDummy3->Draw("copy");

          fitBylinkinUnshifted->Draw("same");

          DrawGammaSetMarkerTGraphAsym(grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[c], 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
          grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[c]->Draw("pEsame");
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMInvYieldStatPbPb2760GeVYshifted[c], 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
          grapharrayPCMInvYieldStatPbPb2760GeVYshifted[c]->Draw("pEsame");

          TLegend* legendYdummy = new TLegend(0.55,0.8,0.85,0.95);
          legendYdummy->SetFillColor(0);
          legendYdummy->SetLineColor(0);
          legendYdummy->SetTextFont(42);
          legendYdummy->SetTextSize(FontSize);
          legendYdummy->AddEntry(grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[c],"combined unshifted","lp");
          legendYdummy->AddEntry(grapharrayPCMInvYieldStatPbPb2760GeVYshifted[c],"combined shifted","lp");
          legendYdummy->Draw();

        canvasDummy3->Update();
        canvasDummy3->Print(Form("%s/%s_ComparisonYShifted_%d.%s",outputDir.Data(),meson.Data(),c,suffix.Data()));

      }
    }

    Double_t limitPar = kTRUE;
    TF1 *fitLowPtBylinkin;
    TF1 *fitHighPtBylinkin;
    TF1 *fitBylinkinPbPb2760GeVPtLHC11h[5];
    TF1 *fitQCDInvYieldPbPb2760GeV;
    for(Int_t c = 0; c< 5; c++){
      fitBylinkinPbPb2760GeVPtLHC11h[c] = NULL;

      //Bylinkin fit ---------------------------------------
      Double_t parampart1TCM[2] = {grapharrayPCMInvYieldTotPbPb2760GeV[c]->GetY()[0],0.3};
      fitLowPtBylinkin = FitObject("Lowtcm","Lowtcm",meson.Data(),grapharrayPCMInvYieldTotPbPb2760GeV[c],0.,6.,parampart1TCM,"QNRMEX0+","",limitPar);
      grapharrayPCMInvYieldTotPbPb2760GeV[c]->Fit(fitLowPtBylinkin,"QNRMEX0+","",0.8,2.);
      Double_t parampart2TCM[3] = {grapharrayPCMInvYieldTotPbPb2760GeV[c]->GetY()[3],0.3,8};
//       cout << WriteParameterToFile(fitLowPtBylinkin)<< endl << endl;
      fitHighPtBylinkin = FitObject("Hightcm","Hightcm",meson.Data(),grapharrayPCMInvYieldTotPbPb2760GeV[c],0.,30.,parampart2TCM,"QNRMEX0+","",limitPar);
      grapharrayPCMInvYieldTotPbPb2760GeV[c]->Fit(fitHighPtBylinkin,"QNRMEX0+","",1.,20.);
//       cout << WriteParameterToFile(fitHighPtBylinkin)<< endl << endl;

      Double_t paramTCM[5] = {fitLowPtBylinkin->GetParameter(0),fitLowPtBylinkin->GetParameter(1),fitHighPtBylinkin->GetParameter(0),fitHighPtBylinkin->GetParameter(1),fitHighPtBylinkin->GetParameter(2)};
      fitBylinkinPbPb2760GeVPtLHC11h[c] = FitObject("tcm","tcm",meson.Data(),grapharrayPCMInvYieldTotPbPb2760GeV[c],fitStart,fitStop,paramTCM,"QNRMEX0+","",limitPar);
      grapharrayPCMInvYieldTotPbPb2760GeV[c]->Fit(fitBylinkinPbPb2760GeVPtLHC11h[c],"NRMEX0+","",fitStart,fitStop);
      cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h[c])<< endl << endl;

      //QCD fit ---------------------------------------
      fitQCDInvYieldPbPb2760GeV = FitObject("qcd","fitQCD",meson.Data(),grapharrayPCMInvYieldTotPbPb2760GeV[c],fitStart,fitStop,NULL,"QNRMEX0+");
      grapharrayPCMInvYieldTotPbPb2760GeV[c]->Fit(fitQCDInvYieldPbPb2760GeV,"QNRMEX0+","",fitStart,fitStop);
      cout << WriteParameterToFile(fitQCDInvYieldPbPb2760GeV)<< endl << endl;

      // Ratio fit to data
      grapharrayRatioToFitStatPbPb2760GeV[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldStatPbPb2760GeV[c]->Clone();
      grapharrayRatioToFitStatPbPb2760GeV[c] = CalculateGraphErrRatioToFit(grapharrayRatioToFitStatPbPb2760GeV[c], fitBylinkinPbPb2760GeVPtLHC11h[c]);
      grapharrayRatioToFitSysPbPb2760GeV[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldSysPbPb2760GeV[c]->Clone();
      grapharrayRatioToFitSysPbPb2760GeV[c] = CalculateGraphErrRatioToFit(grapharrayRatioToFitSysPbPb2760GeV[c], fitBylinkinPbPb2760GeVPtLHC11h[c]);

    }
    textSizeLabelsPixel = 48;
    TCanvas* canvasRatioToCombFitLHC11h = new TCanvas("canvasRatioToCombFitLHC11h","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFitLHC11h, 0.12, 0.03, 0.01, 0.11);
    canvasRatioToCombFitLHC11h->SetLogx();

    Double_t textsizeLabelsPP = 0;
    Double_t textsizeFacPP= 0;
    if (canvasRatioToCombFitLHC11h->XtoPixel(canvasRatioToCombFitLHC11h->GetX2()) <canvasRatioToCombFitLHC11h->YtoPixel(canvasRatioToCombFitLHC11h->GetY1()) ){
        textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFitLHC11h->XtoPixel(canvasRatioToCombFitLHC11h->GetX2()) ;
        textsizeFacPP = (Double_t)1./canvasRatioToCombFitLHC11h->XtoPixel(canvasRatioToCombFitLHC11h->GetX2()) ;
    } else {
        textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFitLHC11h->YtoPixel(canvasRatioToCombFitLHC11h->GetY1());
        textsizeFacPP = (Double_t)1./canvasRatioToCombFitLHC11h->YtoPixel(canvasRatioToCombFitLHC11h->GetY1());
    }

    TH2F * histo2DRatioToCombFitLHC11h = new TH2F("histo2DRatioToCombFitLHC11h","histo2DRatioToCombFitLHC11h",1000,0.23,30.,1000,0.2,4.  );
    SetStyleHistoTH2ForGraphs(histo2DRatioToCombFitLHC11h, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit} ", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DRatioToCombFitLHC11h->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToCombFitLHC11h->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRatioToCombFitLHC11h->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DRatioToCombFitLHC11h->Draw("copy");

    DrawGammaLines(0.23, 30. , 1., 1.,0.5,   kGray+2);
    DrawGammaLines(0.23, 30. , 1.1, 1.1,0.5, kGray, 7);
    DrawGammaLines(0.23, 30. , 0.9, 0.9,0.5, kGray, 7);
    TLegend* legendCombToCombFit = new TLegend(0.2,0.75,0.5,0.94);
    legendCombToCombFit->SetFillColor(0);
    legendCombToCombFit->SetLineColor(0);
    legendCombToCombFit->SetTextFont(42);
    legendCombToCombFit->SetTextSize(0.04);

    for(Int_t c = 0; c< 5; c++){
      DrawGammaSetMarkerTGraphAsym(grapharrayRatioToFitSysPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
      grapharrayRatioToFitSysPbPb2760GeV[c]->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(grapharrayRatioToFitStatPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(grapharrayRatioToFitStatPbPb2760GeV[c]);
      grapharrayRatioToFitStatPbPb2760GeV[c]->Draw("p,same,e0");
      legendCombToCombFit->AddEntry(grapharrayRatioToFitSysPbPb2760GeV[c],cent[c],"fp");
    }
    legendCombToCombFit->Draw();

    TLatex *labelRatioToFitEnergy = new TLatex(0.65,0.92,collisionSystem2760GeV.Data());
    SetStyleTLatex( labelRatioToFitEnergy, 0.85*textSizeLabelsPixel,4);
    labelRatioToFitEnergy->SetTextFont(43);
    TLatex *labelRatioToFitPi0;
    if(meson.CompareTo("Pi0")==0){
        labelRatioToFitPi0= new TLatex(0.73,0.87,"#pi^{0} #rightarrow #gamma#gamma");
    } else if(meson.CompareTo("Eta")==0){
        labelRatioToFitPi0= new TLatex(0.73,0.87,"#eta #rightarrow #gamma#gamma");
    }
    SetStyleTLatex( labelRatioToFitPi0, 0.85*textSizeLabelsPixel,4);
    labelRatioToFitPi0->SetTextFont(43);
    labelRatioToFitEnergy->Draw();
    labelRatioToFitPi0->Draw();

    canvasRatioToCombFitLHC11h->SaveAs(Form("%s/%s_RatioOfCombToCombFit_PbPbPbPb2760GeV.%s",outputDir.Data(),meson.Data(),suffix.Data()));

    //**********************************************************************************************************************//
    //***************************************** Plotting Combined Invariant Yields *****************************************//
    //**********************************************************************************************************************//
    // yields plotting

    TCanvas* canvasInvYieldSectionPi0LHC11h = new TCanvas("canvasInvYieldSectionPi0LHC11h","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldSectionPi0LHC11h, 0.16, 0.02, 0.02, 0.09);
    canvasInvYieldSectionPi0LHC11h->SetLogx();
    canvasInvYieldSectionPi0LHC11h->SetLogy();

    TH2F *histo2DInvYieldSectionPi0LHC11h;
    if(meson.CompareTo("Pi0")==0){
        histo2DInvYieldSectionPi0LHC11h = new TH2F("histo2DInvYieldSectionPi0LHC11h","histo2DInvYieldSectionPi0LHC11h",11000,minPtYields,30.,1000,1e-10,maxYaxisYields);
        SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionPi0LHC11h, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.1,1.5);

    } else if(meson.CompareTo("Eta")==0){
        histo2DInvYieldSectionPi0LHC11h = new TH2F("histo2DInvYieldSectionPi0LHC11h","histo2DInvYieldSectionPi0LHC11h",11000,0.6,20.,1000,1e-10,1e2);
        SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionPi0LHC11h, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.1,1.6);
    }
    histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetMoreLogLabels();
    histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetLabelOffset(-0.01);
    histo2DInvYieldSectionPi0LHC11h->Draw("copy");

      TLegend* legendInvYieldSectionPi0LHC11h_onlyPbPb = new TLegend(0.2,0.12,0.53,0.3);
      legendInvYieldSectionPi0LHC11h_onlyPbPb->SetMargin(0.17);
      legendInvYieldSectionPi0LHC11h_onlyPbPb->SetFillColor(0);
      legendInvYieldSectionPi0LHC11h_onlyPbPb->SetLineColor(0);
      legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextFont(42);
      legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextSize(FontSize);
      legendInvYieldSectionPi0LHC11h_onlyPbPb->SetHeader(collisionSystem2760GeV.Data());

      TGraphAsymmErrors* grapharrayPCMInvYieldStatPbPb2760GeVforPlot[5];
      TGraphAsymmErrors* grapharrayPCMInvYieldSysPbPb2760GeVforPlot[5];
      for(Int_t c=0; c<5; c++){
        grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c] = NULL;
        grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c] = NULL;
        
        if(thesisPlotting){
            while(grapharrayPCMInvYieldStatPbPb2760GeV[c]->GetX()[0] <  1.) grapharrayPCMInvYieldStatPbPb2760GeV[c]->RemovePoint(0);
            while(grapharrayPCMInvYieldSysPbPb2760GeV[c]->GetX()[0] <  1.) grapharrayPCMInvYieldSysPbPb2760GeV[c]->RemovePoint(0);
        }
      }
      for(Int_t c=0; c<5; c++){
        grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldStatPbPb2760GeV[c]->Clone("grapharrayPCMInvYieldStatPbPb2760GeVforPlot");
        grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c] = ScaleGraph(grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c],pow(10,-c));
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c]);
        grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldSysPbPb2760GeV[c]->Clone("grapharrayPCMInvYieldSysPbPb2760GeVforPlot");
        grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c] = ScaleGraph(grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c],pow(10,-c));

        DrawGammaSetMarkerTGraphAsym(grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c], markerStyleDet[c], markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
        grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c], markerStyleDet[c], markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
        grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c]->Draw("p,same,e1");
        legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c],Form("%s x 10^{-%d}",cent[c].Data(),c),"fp");
      }
      legendInvYieldSectionPi0LHC11h_onlyPbPb->Draw();
      labelSystOnlyPbPb->Draw();

    canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldLHC11h_DataOnly.%s",outputDir.Data(),meson.Data(),suffix.Data()));

    canvasInvYieldSectionPi0LHC11h->cd();
    histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetRangeUser(minX,maxX);
    histo2DInvYieldSectionPi0LHC11h->Draw("copy");

      TLegend* legendInvYieldForThesis = new TLegend(0.2,0.12,0.53,0.12+3*0.04);
      legendInvYieldForThesis->SetMargin(0.17);
      legendInvYieldForThesis->SetFillColor(0);
      legendInvYieldForThesis->SetLineColor(0);
      legendInvYieldForThesis->SetTextFont(42);
      legendInvYieldForThesis->SetTextSize(FontSize);
      legendInvYieldForThesis->SetHeader(collisionSystem2760GeV.Data());

      for(Int_t c=0; c<5; c++){
          if(c==2 || c==4){
            grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldStatPbPb2760GeV[c]->Clone("grapharrayPCMInvYieldStatPbPb2760GeVforPlot");
            grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c] = ScaleGraph(grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c],pow(10,-c-2));
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c]);
            grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c] = (TGraphAsymmErrors*)grapharrayPCMInvYieldSysPbPb2760GeV[c]->Clone("grapharrayPCMInvYieldSysPbPb2760GeVforPlot");
            grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c] = ScaleGraph(grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c],pow(10,-c-2));

            DrawGammaSetMarkerTGraphAsym(grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c], markerStyleDet[c], markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
            grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c], markerStyleDet[c], markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
            grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c]->Draw("p,same,e1");
            legendInvYieldForThesis->AddEntry(grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c],Form("%s x 10^{-%d}",cent[c].Data(),c),"fp");
          }
      }
      legendInvYieldForThesis->Draw();
      labelSystOnlyPbPb->Draw();

      if(thesisPlotting){
        TLatex *thesisLabel = new TLatex(0.75,0.9,thisthesis.Data());
        SetStyleTLatex( thesisLabel,FontSize,4);
        thesisLabel->Draw();
    }

    canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldForThesis_DataOnly.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    
    canvasInvYieldSectionPi0LHC11h->cd();
    histo2DInvYieldSectionPi0LHC11h->Draw("copy");

      TLegend* legendInvYieldSectionPi0LHC11h_WithFit= new TLegend(0.595,0.75,0.91,0.91); //0.17,0.13,0.5,0.24);
      legendInvYieldSectionPi0LHC11h_WithFit->SetFillColor(0);
      legendInvYieldSectionPi0LHC11h_WithFit->SetMargin(0.17);
      legendInvYieldSectionPi0LHC11h_WithFit->SetLineColor(0);
      legendInvYieldSectionPi0LHC11h_WithFit->SetTextFont(42);
      legendInvYieldSectionPi0LHC11h_WithFit->SetTextSize(FontSize);
      legendInvYieldSectionPi0LHC11h_WithFit->SetHeader(collisionSystem2760GeV.Data());
      TF1 *fitBylinkinPbPb2760GeVPtLHC11hforPlot[5];
      for(Int_t c=0; c<5; c++){
          grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c]->Draw("E2same");
          grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c]->Draw("p,same,e1");

          fitBylinkinPbPb2760GeVPtLHC11hforPlot[c] = NULL;
          fitBylinkinPbPb2760GeVPtLHC11hforPlot[c] = (TF1*)fitBylinkinPbPb2760GeVPtLHC11h[c]->Clone();
          fitBylinkinPbPb2760GeVPtLHC11hforPlot[c]->SetLineStyle(4);
          fitBylinkinPbPb2760GeVPtLHC11hforPlot[c]->SetParameter(0,fitBylinkinPbPb2760GeVPtLHC11hforPlot[c]->GetParameter(0)*pow(10,-c));
          fitBylinkinPbPb2760GeVPtLHC11hforPlot[c]->SetParameter(2,fitBylinkinPbPb2760GeVPtLHC11hforPlot[c]->GetParameter(2)*pow(10,-c));
          fitBylinkinPbPb2760GeVPtLHC11hforPlot[c]->Draw("same");

          legendInvYieldSectionPi0LHC11h_WithFit->AddEntry(grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c],Form("%s x 10^{-%d}",cent[c].Data(),c),"fp");
          if(c==4)legendInvYieldSectionPi0LHC11h_WithFit->AddEntry(fitBylinkinPbPb2760GeVPtLHC11hforPlot[c],"fits to Pb#font[122]{-}Pb","l");
      }
      legendInvYieldSectionPi0LHC11h_WithFit->Draw();
      labelSyst->Draw();

    canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnlyWithFit.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    canvasInvYieldSectionPi0LHC11h->cd();
    TH2F *histo2DInvYieldSectionLHC11hwithPP;
    if(meson.CompareTo("Pi0")==0){
        histo2DInvYieldSectionLHC11hwithPP = new TH2F("histo2DInvYieldSectionLHC11hwithPP","histo2DInvYieldSectionLHC11hwithPP",11000,0.25,30.,1000,1e-12,1e3);
        SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionLHC11hwithPP, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2} ",0.035,0.04, 0.035,0.04, 1.,1.65);

    } else if(meson.CompareTo("Eta")==0){
        histo2DInvYieldSectionLHC11hwithPP = new TH2F("histo2DInvYieldSectionLHC11hwithPP","histo2DInvYieldSectionLHC11hwithPP",11000,0.4,30.,1000,2e-10,1e3);
        SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionLHC11hwithPP, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2} ",0.035,0.04, 0.035,0.04, 1.,1.6);
    }
    histo2DInvYieldSectionLHC11hwithPP->GetXaxis()->SetMoreLogLabels();
    histo2DInvYieldSectionLHC11hwithPP->GetXaxis()->SetLabelOffset(-0.01);
    histo2DInvYieldSectionLHC11hwithPP->Draw("copy");

      TLegend* legendSpectraPP = new TLegend(0.2,0.15,0.54,0.22);
      legendSpectraPP->SetFillColor(0);
      legendSpectraPP->SetLineColor(0);
      legendSpectraPP->SetTextFont(42);
      legendSpectraPP->SetMargin(0.17);
      legendSpectraPP->SetTextSize(FontSize);
      if(meson.CompareTo("Pi0")==0){

        graphInvSectionCombStatPi02760GeVPlot->Draw("p,same,e0");
        graphInvSectionCombSysPi02760GeVPlot->Draw("E2same");

        legendSpectraPP->SetHeader(collisionSystemPP2760GeV.Data());
        legendSpectraPP->AddEntry(graphInvSectionCombSysPi02760GeVPlot,"arXiv:XXXX.XXXX","pf");

      } else if(meson.CompareTo("Eta")==0){

          graphInvSectionCombStatEta2760GeVPlot->Draw("p,same,e0");
          graphInvSectionCombSysEta2760GeVPlot->Draw("E2same");

          legendSpectraPP->SetHeader(collisionSystemPP2760GeV.Data());
          legendSpectraPP->AddEntry(graphInvSectionCombSysEta2760GeVPlot,"arXiv:XXXX.XXXX","pf");

      }
      legendSpectraPP->Draw();

      TLegend* legendInvYieldSectionPi0LHC11h_WitPP = new TLegend(0.595,0.75,0.91,0.91); //0.17,0.13,0.5,0.24);
      legendInvYieldSectionPi0LHC11h_WitPP->SetFillColor(0);
      legendInvYieldSectionPi0LHC11h_WitPP->SetMargin(0.17);
      legendInvYieldSectionPi0LHC11h_WitPP->SetLineColor(0);
      legendInvYieldSectionPi0LHC11h_WitPP->SetTextFont(42);
      legendInvYieldSectionPi0LHC11h_WitPP->SetTextSize(FontSize);
      legendInvYieldSectionPi0LHC11h_WitPP->SetHeader(collisionSystem2760GeV.Data());
      for(Int_t c=0; c<5; c++){
          grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c]->Draw("E2same");
          grapharrayPCMInvYieldStatPbPb2760GeVforPlot[c]->Draw("p,same,e1");
          fitBylinkinPbPb2760GeVPtLHC11hforPlot[c]->Draw("same");
          legendInvYieldSectionPi0LHC11h_WitPP->AddEntry(grapharrayPCMInvYieldSysPbPb2760GeVforPlot[c],Form("%s x 10^{-%d}",cent[c].Data(),c),"fp");
          if(c==4)legendInvYieldSectionPi0LHC11h_WitPP->AddEntry(fitBylinkinPbPb2760GeVPtLHC11hforPlot[c],"fits to Pb#font[122]{-}Pb","l");
      }
      legendInvYieldSectionPi0LHC11h_WitPP->Draw();
      labelSyst->Draw();

    canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly_withPP.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    //**********************************************************************************************************************//
    //***********************************************     Raa combination    ***********************************************//
    //**********************************************************************************************************************//
    // RAA calc and plotting

    Bool_t quiet = kFALSE;
    TGraphAsymmErrors* graphRAAPCM[5];
    TGraphAsymmErrors* graphRAASysPCM[5];

    if(meson.CompareTo("Pi0")==0){

      for(Int_t c=0; c<5; c++){
        CalcRaa(    graphInvSectionPCMStatPi02760GeVforRAA, graphInvSectionPCMSysPi02760GeVforRAA,graphInvSectionCombStatPi02760GeVforRAA, fitInvCrossSectionTsallisPi0Comb2760GeV,
                    grapharrayPCMInvYieldStatPbPb2760GeVYshifted[c], grapharrayPCMInvYieldSysPbPb2760GeVYshifted[c],
                    &graphRAAPCM[c], &graphRAASysPCM[c], nColl[c], nCollErr[c],"Pi0",8.,0,"h",quiet);
      }

    } else if(meson.CompareTo("Eta")==0){

      for(Int_t c=0; c<5; c++){
        CalcRaa(    graphInvSectionPCMStatEta2760GeVforRAA, graphInvSectionPCMSysEta2760GeVforRAA,graphInvSectionCombStatEta2760GeVforRAA, fitInvCrossSectionTsallisEtaComb2760GeV,
                    grapharrayPCMInvYieldStatPbPb2760GeVYshifted[c], grapharrayPCMInvYieldSysPbPb2760GeVYshifted[c],
                    &graphRAAPCM[c], &graphRAASysPCM[c],
                    nColl[c], nCollErr[c],"Eta",6.,0,"h",quiet);
      }
    }

    TCanvas* canvasRAAMeasurements = new TCanvas("canvasRAAMeasurements","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRAAMeasurements, 0.09, 0.02, 0.035, 0.09);
    TH2F * histo2DRAADummy = new TH2F("histo2DRAADummy","histo2DRAADummy",1000,0.,16.,1000,0,1.2);
    SetStyleHistoTH2ForGraphs(histo2DRAADummy, "p_{T} (GeV/c)","R_{AA}", 0.032,0.032, 0.04,0.04, 1,1.);
    histo2DRAADummy->DrawCopy();


      TLegend* legendRAAcombo2 = new TLegend(0.6,0.6,0.91,0.8);
      legendRAAcombo2->SetFillColor(0);
      legendRAAcombo2->SetLineColor(0);
      legendRAAcombo2->SetTextFont(42);
      legendRAAcombo2->SetMargin(0.17);
      legendRAAcombo2->SetTextSize(FontSize);
      legendRAAcombo2->SetHeader(collisionSystem2760GeV.Data());

      for(Int_t c=0; c<5; c++){
        DrawGammaSetMarkerTGraphAsym(graphRAASysPCM[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
        graphRAASysPCM[c]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRAAPCM[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRAAPCM[c]);
        graphRAAPCM[c]->Draw("p,same,e");
        legendRAAcombo2->AddEntry(graphRAASysPCM[c],cent[c],"fp");
      }
      legendRAAcombo2->Draw();

      DrawGammaLines(0.3, 30. , 1., 1.,0.5,   kGray);
      TLatex *labelSystRaa;
      if(meson.CompareTo("Pi0")==0){     labelSystRaa= new TLatex(0.6,0.83,"#pi^{0} #rightarrow #gamma#gamma");}
      else if(meson.CompareTo("Eta")==0){labelSystRaa= new TLatex(0.6,0.83,"#eta #rightarrow #gamma#gamma");}
      SetStyleTLatex( labelSystRaa,FontSize,4);
      labelSystRaa->Draw();
      boxErrorNorm0005->Draw();
      boxErrorNorm0510->Draw();
      boxErrorNorm0010->Draw();
      boxErrorNorm2040->Draw();
      boxErrorNorm2050->Draw();

//     canvasRAAMeasurements->SaveAs(Form("%s/%s_recalculatedRAA_DataOnly.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    canvasRAAMeasurements->cd();
    histo2DRAADummy->DrawCopy();

      for(Int_t c=0; c<5; c++){
        DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAASysPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
        grapharrayPCMRAASysPbPb2760GeV[c]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAAStatPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(grapharrayPCMRAAStatPbPb2760GeV[c]);
        grapharrayPCMRAAStatPbPb2760GeV[c]->Draw("p,same,e");
      }
      legendRAAcombo2->Draw();

      DrawGammaLines(0., 16. , 1., 1.,0.5,   kGray);
      labelSystRaa->Draw();
      boxErrorNorm0005->Draw();
      boxErrorNorm0510->Draw();
      boxErrorNorm0010->Draw();
      boxErrorNorm2040->Draw();
      boxErrorNorm2050->Draw();

    canvasRAAMeasurements->SaveAs(Form("%s/%s_RAA_DataOnly.%s",outputDir.Data(),meson.Data(),suffix.Data()));


    for(Int_t c=0; c<5; c++){
      canvasRAAMeasurements->cd();
      histo2DRAADummy->DrawCopy();


        TLegend* legendRAATheoryPbPb = new TLegend(0.6,0.75,0.91,0.8);
        legendRAATheoryPbPb->SetFillColor(0);
        legendRAATheoryPbPb->SetLineColor(0);
        legendRAATheoryPbPb->SetTextFont(42);
        legendRAATheoryPbPb->SetTextSize(FontSize);
//         legendRAATheoryPbPb->SetMargin(0.2);
//         legendRAATheoryPbPb->SetNColumns(2);
        legendRAATheoryPbPb->AddEntry(graphRAASysPCM[c],cent[c],"fp");

        if(meson.CompareTo("Pi0")==0){

          if(c==0){

            DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_0005, markerStyleCommmonSpectrum0005,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
            gWHDG_Raa_0005->SetFillStyle(fillStyleWHDG);
            gWHDG_Raa_0005->SetFillColor(colorWHDG0005);
            gWHDG_Raa_0005->Draw("3 same");
            DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_0005, markerStyleCommmonSpectrum0005,markerSizeSpectrum, colorVitevBas0005 ,colorVitevBas0005,2);
            Vitev_Bas_Raa_0005->SetFillStyle(fillStyleVitev);
            Vitev_Bas_Raa_0005->SetFillColor(colorVitevBas0005);
            Vitev_Bas_Raa_0005->Draw("3 same");
            Vitev_Bas_Raa_0005->Print();

            TLegend* legTheory1 = new TLegend(0.13,0.7,0.3,0.8);
            legTheory1->SetFillColor(0);
            legTheory1->SetLineColor(0);
            legTheory1->SetTextFont(42);
            legTheory1->SetTextSize(FontSize);
            legTheory1->AddEntry(gWHDG_Raa_0005,"WHDG","f");
            legTheory1->AddEntry(Vitev_Bas_Raa_0005,"Vitev","f");
            legTheory1->Draw();

          } else if(c==1){

            DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_0510, markerStyleCommmonSpectrum0510,markerSizeSpectrum, colorVitevBas0510 ,colorVitevBas0510,2);
            Vitev_Bas_Raa_0510->SetFillStyle(fillStyleVitev);
            Vitev_Bas_Raa_0510->SetFillColor(colorVitevBas0510);
            Vitev_Bas_Raa_0510->Draw("3 same");
            DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_0510, markerStyleCommmonSpectrum0510,markerSizeSpectrum, colorWHDG0510, colorWHDG0510,0.8);
            gWHDG_Raa_0510->SetFillStyle(fillStyleWHDG);
            gWHDG_Raa_0510->SetFillColor(colorWHDG0510);
            gWHDG_Raa_0510->Draw("3 same");

            TLegend* legTheory2 = new TLegend(0.13,0.7,0.3,0.8);
            legTheory2->SetFillColor(0);
            legTheory2->SetLineColor(0);
            legTheory2->SetTextFont(42);
            legTheory2->SetTextSize(FontSize);
            legTheory2->AddEntry(gWHDG_Raa_0510,"WHDG","f");
            legTheory2->AddEntry(Vitev_Bas_Raa_0510,"Vitev","f");
            legTheory2->Draw();

          } else if(c==2){

            DrawGammaSetMarkerTGraphAsym(graphPi0RAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
            graphPi0RAAJetQuenching_0010->Draw("3,same");
            DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
            gWHDG_Raa_0010->SetFillStyle(fillStyleWHDG);
            gWHDG_Raa_0010->SetFillColor(colorWHDG0005);
            gWHDG_Raa_0010->Draw("3 same");
            DrawGammaSetMarkerTGraphAsym(graphPi0Djordjevic_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorVitevBas0005, colorVitevBas0005,widthLinesBoxes, kTRUE);
            graphPi0Djordjevic_0010->SetFillStyle(fillStyleVitev);
            graphPi0Djordjevic_0010->SetFillColor(colorVitevBas0005);
            graphPi0Djordjevic_0010->Draw("3 same");

            TLegend* legTheory3 = new TLegend(0.13,0.65,0.3,0.8);
            legTheory3->SetFillColor(0);
            legTheory3->SetLineColor(0);
            legTheory3->SetTextFont(42);
            legTheory3->SetTextSize(FontSize);
            legTheory3->AddEntry(graphPi0RAAJetQuenching_0010,"NLO DCZW","f");
            legTheory3->AddEntry(gWHDG_Raa_0010,"WHDG","f");
            legTheory3->AddEntry(graphPi0Djordjevic_0010,"Djordjevic","f");
            legTheory3->Draw();

          } else if(c==3){

            DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_2040, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorVitevBas2040 ,colorVitevBas2040,2);
            Vitev_Bas_Raa_2040->SetFillStyle(fillStyleVitev);
            Vitev_Bas_Raa_2040->SetFillColor(colorVitevBas2040);
            Vitev_Bas_Raa_2040->Draw("3 same");
            DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_2040, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG2040, colorWHDG2040,0.8);
            gWHDG_Raa_2040->SetFillStyle(fillStyleWHDG);
            gWHDG_Raa_2040->SetFillColor(colorWHDG2040);
            gWHDG_Raa_2040->Draw("3 same");

            TLegend* legTheory4 = new TLegend(0.13,0.7,0.3,0.8);
            legTheory4->SetFillColor(0);
            legTheory4->SetLineColor(0);
            legTheory4->SetTextFont(42);
            legTheory4->SetTextSize(FontSize);
            legTheory4->AddEntry(gWHDG_Raa_2040,"WHDG","f");
            legTheory4->AddEntry(Vitev_Bas_Raa_2040,"Vitev","f");
            legTheory4->Draw();

          } else if(c==4){

            DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
            gWHDG_Raa_2050->SetFillStyle(fillStyleWHDG);
            gWHDG_Raa_2050->SetFillColor(colorWHDG4060);
            gWHDG_Raa_2050->Draw("3 same");
            DrawGammaSetMarkerTGraphAsym(graphPi0Djordjevic_2050, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorVitevBas4060,colorVitevBas4060,widthLinesBoxes, kTRUE);
            graphPi0Djordjevic_2050->SetFillStyle(fillStyleVitev);
            graphPi0Djordjevic_2050->SetFillColor(colorVitevBas4060);
            graphPi0Djordjevic_2050->Draw("3 same");

            TLegend* legTheory5 = new TLegend(0.13,0.7,0.3,0.8);
            legTheory5->SetFillColor(0);
            legTheory5->SetLineColor(0);
            legTheory5->SetTextFont(42);
            legTheory5->SetTextSize(FontSize);
            legTheory5->AddEntry(gWHDG_Raa_2050,"WHDG","f");
            legTheory5->AddEntry(graphPi0Djordjevic_2050,"Djordjevic","f");
            legTheory5->Draw();
          }
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAASysPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
          grapharrayPCMRAASysPbPb2760GeV[c]->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAAStatPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(grapharrayPCMRAAStatPbPb2760GeV[c]);
          grapharrayPCMRAAStatPbPb2760GeV[c]->Draw("p,same,e");


        } else if(meson.CompareTo("Eta")==0){

          if(c==2){
            DrawGammaSetMarkerTGraphAsym(graphEtaRAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
            graphEtaRAAJetQuenching_0010->Draw("3,same");

            DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
            gWHDG_Eta_Raa_0010->SetFillStyle(fillStyleWHDG);
            gWHDG_Eta_Raa_0010->SetFillColor(colorWHDG0005);
            gWHDG_Eta_Raa_0010->Draw("3 same");

              TLegend* legTheory6 = new TLegend(0.13,0.7,0.3,0.8);
              legTheory6->SetFillColor(0);
              legTheory6->SetLineColor(0);
              legTheory6->SetTextFont(42);
              legTheory6->SetTextSize(FontSize);
              legTheory6->AddEntry(graphEtaRAAJetQuenching_0010,"NLO DCZW","f");
              legTheory6->AddEntry(gWHDG_Eta_Raa_0010,"WHDG","f");
              legTheory6->Draw();

          } else if(c==4){

            DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
            gWHDG_Eta_Raa_2050->SetFillStyle(fillStyleWHDG);
            gWHDG_Eta_Raa_2050->SetFillColor(colorWHDG4060);
            gWHDG_Eta_Raa_2050->Draw("3 same");

              TLegend* legTheory7 = new TLegend(0.13,0.75,0.3,0.8);
              legTheory7->SetFillColor(0);
              legTheory7->SetLineColor(0);
              legTheory7->SetTextFont(42);
              legTheory7->SetTextSize(FontSize);
              legTheory7->AddEntry(gWHDG_Eta_Raa_2050,"WHDG","f");
              legTheory7->Draw();
          }
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAASysPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
          grapharrayPCMRAASysPbPb2760GeV[c]->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAAStatPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(grapharrayPCMRAAStatPbPb2760GeV[c]);
          grapharrayPCMRAAStatPbPb2760GeV[c]->Draw("p,same,e");

        }
        legendRAATheoryPbPb->Draw();
        labelSystRaa->Draw();
        DrawGammaLines(0., 16. , 1, 1 ,1,kGray, 2);

      canvasRAAMeasurements->SaveAs(Form("%s/%s_RAA_combinedWithTheoryModels_%d.%s",outputDir.Data(),meson.Data(),c,suffix.Data()));


      canvasRAAMeasurements->cd();
      histo2DRAADummy->DrawCopy();

        labelSystRaa->Draw();
        if(meson.CompareTo("Pi0")==0){

          DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA[c], 24,2, colorCharged,colorCharged, 0.1, kFALSE);
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedPionRAA[c]);
          graphChargedPionRAA[c]->Draw("p,same");
          DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys[c],24,2, colorCharged,colorCharged, 1, kTRUE);
          graphChargedPionRAASys[c]->Draw("2same");

          TLegend* legCharged1 = new TLegend(0.55,0.6,0.91,0.8);
          legCharged1->SetFillColor(0);
          legCharged1->SetLineColor(0);
          legCharged1->SetTextFont(42);
          legCharged1->SetMargin(0.17);
          legCharged1->SetTextSize(FontSize);
          legCharged1->SetHeader(Form("%s %s",cent[c].Data(),collisionSystem2760GeV.Data()));
          legCharged1->AddEntry(grapharrayPCMRAASysPbPb2760GeV[c],"#pi^{0}","fp");
          legCharged1->AddEntry(graphChargedPionRAA[c],"#pi^{#pm}","fp");
          legCharged1->AddEntry((TObject*)0,"PLB 736 (2014)","");
          legCharged1->Draw();

          DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAASysPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
          grapharrayPCMRAASysPbPb2760GeV[c]->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAAStatPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(grapharrayPCMRAAStatPbPb2760GeV[c]);
          grapharrayPCMRAAStatPbPb2760GeV[c]->Draw("p,same,e");

        } else if(meson.CompareTo("Eta")==0){

          DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA[c], 24,2, colorCharged,colorCharged, 0.1, kFALSE);
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAA[c]);
          graphChargedKaonRAA[c]->Draw("p,same");
          DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys[c],24,2, colorCharged,colorCharged, 1, kTRUE);
          graphChargedKaonRAASys[c]->Draw("2same");

          TLegend* legCharged2 = new TLegend(0.55,0.6,0.91,0.8);
          legCharged2->SetFillColor(0);
          legCharged2->SetLineColor(0);
          legCharged2->SetTextFont(42);
          legCharged2->SetMargin(0.17);
          legCharged2->SetTextSize(FontSize);
          legCharged2->SetHeader(Form("%s %s",cent[c].Data(),collisionSystem2760GeV.Data()));
          legCharged2->AddEntry(grapharrayPCMRAASysPbPb2760GeV[c],"#eta","fp");
          legCharged2->AddEntry(graphChargedKaonRAA[c],"K^{#pm}","fp");
          legCharged2->AddEntry((TObject*)0,"PLB 736 (2014)","");
          legCharged2->Draw();

          DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAASysPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
          grapharrayPCMRAASysPbPb2760GeV[c]->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMRAAStatPbPb2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(grapharrayPCMRAAStatPbPb2760GeV[c]);
          grapharrayPCMRAAStatPbPb2760GeV[c]->Draw("p,same,e");

        }
        labelSystRaa->Draw();
        DrawGammaLines(0., 16. , 1, 1 ,1,kGray, 2);

      canvasRAAMeasurements->SaveAs(Form("%s/%s_RAA_combinedWithCharged_%d.%s",outputDir.Data(),meson.Data(),c,suffix.Data()));


    }

    if(meson.CompareTo("Eta")==0){

      TFile* fEtatoPi0input = new TFile("EtaToPi0InputsForCombination.root");
      TGraphAsymmErrors *graphCombEtaToPi0RatioSysErrpp7TeV = (TGraphAsymmErrors*)fEtatoPi0input->Get("graphCombEtaToPi0RatioSysErrpp7TeV");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RatioSysErrpp7TeV, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);
      TGraphAsymmErrors *graphCombEtaToPi0Ratiopp7TeVNoXErrors = (TGraphAsymmErrors*)fEtatoPi0input->Get("graphCombEtaToPi0Ratiopp7TeVNoXErrors");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Ratiopp7TeVNoXErrors, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);

      //**********************************************************************************************************************//
      //********************************************   Eta to Pi0 combination    *********************************************//
      //**********************************************************************************************************************//
      //eta to pi0 ratio

      TCanvas* canvasEtatoPi0combo = new TCanvas("canvasEtatoPi0combo","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
      DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.09, 0.03, 0.035, 0.1);
      canvasEtatoPi0combo->SetLogx();

      TH2F * histo2DEtatoPi0combo = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",11000,0.23,25.,1000,0.01,1.2);
      SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}",0.035,0.04, 0.035,0.04, 1.2,1.);
      histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
      histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
      histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.01,maxYEtatoPi0);
      histo2DEtatoPi0combo->Draw("copy");

        TLegend* legendEtatoPi0combo_onlyPbPb = new TLegend(0.12,0.73,0.5,0.92);
        legendEtatoPi0combo_onlyPbPb->SetFillColor(0);
        legendEtatoPi0combo_onlyPbPb->SetLineColor(0);
        legendEtatoPi0combo_onlyPbPb->SetTextFont(42);
        legendEtatoPi0combo_onlyPbPb->SetTextSize(0.037);
        legendEtatoPi0combo_onlyPbPb->SetMargin(0.17);
        legendEtatoPi0combo_onlyPbPb->SetHeader(collisionSystem2760GeV.Data());
        for(Int_t c=0; c<5; c++){
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMEtatoPi0Sys2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]), widthLinesBoxes, kTRUE);
          grapharrayPCMEtatoPi0Sys2760GeV[c]->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(grapharrayPCMEtatoPi0Stat2760GeV[c], markerStyleDet[c] ,markerSizeDet[c]*0.5, GetColorDefaultColor(energy.Data(),"",cent[c]), GetColorDefaultColor(energy.Data(),"",cent[c]));
          if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(grapharrayPCMEtatoPi0Stat2760GeV[c]);
          grapharrayPCMEtatoPi0Stat2760GeV[c]->Draw("p,same,e");
          legendEtatoPi0combo_onlyPbPb->AddEntry(grapharrayPCMEtatoPi0Sys2760GeV[c],cent[c],"fp");
        }
        legendEtatoPi0combo_onlyPbPb->Draw();

      canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioDataOnly.%s",outputDir.Data(),suffix.Data()));
      canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_DataOnly.%s",paperPlots.Data(),suffix.Data()));


      canvasEtatoPi0combo->cd();
      histo2DEtatoPi0combo->Draw("copy");

        for(Int_t c=0; c<5; c++){
          grapharrayPCMEtatoPi0Sys2760GeV[c]->Draw("E2same");
          grapharrayPCMEtatoPi0Stat2760GeV[c]->Draw("p,same,e");
        }
        legendEtatoPi0combo_onlyPbPb->Draw();

        graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
        graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");

        TLegend* legendEtatoPi0combo_withPP = new TLegend(0.55,0.15,0.95,0.26);
        legendEtatoPi0combo_withPP->SetFillColor(0);
        legendEtatoPi0combo_withPP->SetLineColor(0);
        legendEtatoPi0combo_withPP->SetTextFont(42);
        legendEtatoPi0combo_withPP->SetTextSize(0.037);
        legendEtatoPi0combo_withPP->SetMargin(0.17);
        legendEtatoPi0combo_withPP->SetHeader(Form("#eta/#pi^{0} %s",collisionSystemPP7TeV.Data()));
        legendEtatoPi0combo_withPP->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,"PLB 717 (2012) 162","fp");//"Phys. Lett. B 717 (2012) 162-172","fp");
        legendEtatoPi0combo_withPP->Draw();

      canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_DataOnlyWithPP.%s",outputDir.Data(),suffix.Data()));
      canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_DataOnlyWithPP.%s",PubNotePlots.Data(),suffix.Data()));

      for(Int_t c=0; c<5; c++){
        canvasEtatoPi0combo->cd();
        histo2DEtatoPi0combo->Draw("copy");
        histo2DEtatoPi0combo->GetYaxis()->SetTitleOffset(1.1);
        histo2DEtatoPi0combo->GetYaxis()->SetTitle("Particle ratio");

          grapharrayPCMEtatoPi0Sys2760GeV[c]->Draw("E2same");
          grapharrayPCMEtatoPi0Stat2760GeV[c]->Draw("p,same,e");

          DrawGammaSetMarkerTGraphAsym(graphChargedKaonToPion[c], 24,markerSizeComb, colorCharged, colorCharged);
          graphChargedKaonToPion[c]->Draw("p,same");
          DrawGammaSetMarkerTGraphAsym(graphChargedKaonToPionSys[c],24,markerSizeComb, colorCharged, colorCharged, 1, kTRUE);
          graphChargedKaonToPionSys[c]->Draw("2same");

          TLegend* legendChargedRatio = new TLegend(0.12,0.76,0.5,0.92);
          legendChargedRatio->SetFillColor(0);
          legendChargedRatio->SetLineColor(0);
          legendChargedRatio->SetTextFont(42);
          legendChargedRatio->SetTextSize(0.037);
          legendChargedRatio->SetMargin(0.17);
          legendChargedRatio->SetHeader(collisionSystem2760GeV.Data());
//           legendChargedRatio->SetNColumns(2);
          legendChargedRatio->AddEntry(grapharrayPCMEtatoPi0Sys2760GeV[c],Form("#eta/#pi^{0},   %s",cent[c].Data()),"fp");
//           legendChargedRatio->AddEntry((TObject*)0,"","");
          legendChargedRatio->AddEntry(graphChargedKaonToPionSys[c],Form("K^{#pm}/#pi^{#pm}, %s",cent[c].Data()),"fp");
          legendChargedRatio->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
          legendChargedRatio->Draw();

        canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0Ratio_withKaonsToPions_%d.%s",outputDir.Data(),c,suffix.Data()));
      }
    }
//     //******************************************************************************//
//     //************************* Saving of final results ****************************//
//     //******************************************************************************//
//
//     TFile *fCombResults = new TFile(Form("%s/CombinedResultsPaperPbPb2760GeV_%s.root", outputDir.Data(),dateForOutput.Data()), "UPDATE");
//
//         //x shifted spectra
//         graphCombInvYieldTotPbPb2760GeV_0010->Write(Form("graphInvYield%sCombPbPb2760GeV_0010",meson.Data()));
//         graphCombInvYieldStatPbPb2760GeV_0010->Write(Form("graphInvYield%sCombPbPb2760GeVStatErr_0010",meson.Data()));
//         graphCombInvYieldSysPbPb2760GeV_0010->Write(Form("graphInvYield%sCombPbPb2760GeVSysErr_0010",meson.Data()));
//         graphCombInvYieldTotPbPb2760GeV_2050->Write(Form("graphInvYield%sCombPbPb2760GeV_2050",meson.Data()));
//         graphCombInvYieldStatPbPb2760GeV_2050->Write(Form("graphInvYield%sCombPbPb2760GeVStatErr_2050",meson.Data()));
//         graphCombInvYieldSysPbPb2760GeV_2050->Write(Form("graphInvYield%sCombPbPb2760GeVSysErr_2050",meson.Data()));
//
//         graphPCMPi0InvYieldStatPbPb2760GeV_0010->Write(Form("graphInvYield%sPCMPbPb2760GeVStatErr_0010",meson.Data()));
//         graphPCMPi0InvYieldSysPbPb2760GeV_0010->Write(Form("graphInvYield%sPCMPbPb2760GeVSysErr_0010",meson.Data()));
//         graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Write(Form("graphInvYield%sEMCalPbPb2760GeVStatErr_0010",meson.Data()));
//         graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Write(Form("graphInvYield%sEMCalPbPb2760GeVSysErr_0010",meson.Data()));
//         graphPHOSPi0InvYieldStatPbPb2760GeV_0010->Write(Form("graphInvYield%sPHOSPbPb2760GeVStatErr_0010",meson.Data()));
//         graphPHOSPi0InvYieldSysPbPb2760GeV_0010->Write(Form("graphInvYield%sPHOSPbPb2760GeVSysErr_0010",meson.Data()));
//         graphPCMPi0InvYieldStatPbPb2760GeV_2050->Write(Form("graphInvYield%sPCMPbPb2760GeVStatErr_2050",meson.Data()));
//         graphPCMPi0InvYieldSysPbPb2760GeV_2050->Write(Form("graphInvYield%sPCMPbPb2760GeVSysErr_2050",meson.Data()));
//         graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Write(Form("graphInvYield%sEMCalPbPb2760GeVStatErr_2050",meson.Data()));
//         graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Write(Form("graphInvYield%sEMCalPbPb2760GeVSysErr_2050",meson.Data()));
//
//         //y shifted spectra
//         graphCombInvYieldSysPbPb2760GeVYShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVSysErrYshifted_0010",meson.Data()));
//         graphCombInvYieldStatPbPb2760GeVYShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVStatErrYshifted_0010",meson.Data()));
//         graphCombInvYieldSysPbPb2760GeVYShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVSysErrYshifted_2050",meson.Data()));
//         graphCombInvYieldStatPbPb2760GeVYShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVStatErrYshifted_2050",meson.Data()));
//
//         fitBylinkinPbPb2760GeVPtLHC11h_0010->Write(Form("FitToYield%s_0010",meson.Data()));
//         fitBylinkinPbPb2760GeVPtLHC11h_2050->Write(Form("FitToYield%s_2050",meson.Data()));
//
//         graphCombRAASysPbPb2760GeV_0010->Write(Form("graphRAA%sCombPbPb2760GeVSysErr_0010",meson.Data()));
//         graphCombRAAStatPbPb2760GeV_0010->Write(Form("graphRAA%sCombPbPb2760GeVStatErr_0010",meson.Data()));
//         graphCombRAASysPbPb2760GeV_2050->Write(Form("graphRAA%sCombPbPb2760GeVSysErr_2050",meson.Data()));
//         graphCombRAAStatPbPb2760GeV_2050->Write(Form("graphRAA%sCombPbPb2760GeVStatErr_2050",meson.Data()));
//
//         if(graphCombEtatoPi0SysPbPb2760GeV_0010)graphCombEtatoPi0SysPbPb2760GeV_0010->Write("graphEtaToPi0RAtioCombPbPb2760GeVSysErr_0010");
//         if(graphCombEtatoPi0StatPbPb2760GeV_0010)graphCombEtatoPi0StatPbPb2760GeV_0010->Write("graphEtaToPi0RAtioCombPbPb2760GeVStatErr_0010");
//         if(graphCombEtatoPi0SysPbPb2760GeV_2050)graphCombEtatoPi0SysPbPb2760GeV_2050->Write("graphEtaToPi0RAtioCombPbPb2760GeVSysErr_2050");
//         if(graphCombEtatoPi0StatPbPb2760GeV_2050)graphCombEtatoPi0StatPbPb2760GeV_2050->Write("graphEtaToPi0RAtioCombPbPb2760GeVStatErr_2050");
//
//         histoEPOSpred_0010->Write(Form("EPOSprediction%s_0010",meson.Data()));
//         histoEPOSpred_2050->Write(Form("EPOSprediction%s_2050",meson.Data()));
//         TheoryCracowLowPt_0010->Write(Form("Cracowprediction%s_0010",meson.Data()));
//         TheoryCracowLowPt_2050->Write(Form("Cracowprediction%s_2050",meson.Data()));
//
//         //for comparison to charged pions:
//         graphPCMPi0InvYieldSysPbPb2760GeV_2040->Write("graphPCMPi0InvYieldSysPbPb2760GeV_2040");
//         graphPCMPi0InvYieldStatPbPb2760GeV_2040->Write("graphPCMPi0InvYieldStatPbPb2760GeV_2040");
//
//         graphPCMEtaInvYieldSysPbPb2760GeV_2040->Write("graphPCMEtaInvYieldSysPbPb2760GeV_2040");
//         graphPCMEtaInvYieldStatPbPb2760GeV_2040->Write("graphPCMEtaInvYieldStatPbPb2760GeV_2040");
//
//
//         fCombResults->mkdir("UnshiftedSpectra");
//         TDirectoryFile* directoryNoShift = (TDirectoryFile*)fCombResults->Get("UnshiftedSpectra");
//         fCombResults->cd("UnshiftedSpectra");
//             //unshifted spectra
//             graphCombInvYieldTotPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVNoShift_0010",meson.Data()));
//             graphCombInvYieldStatPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVStatErrNoShift_0010",meson.Data()));
//             graphCombInvYieldSysPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sCombPbPb2760GeVSysErrNoShift_0010",meson.Data()));
//             graphCombInvYieldTotPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVNoShift_2050",meson.Data()));
//             graphCombInvYieldStatPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVStatErrNoShift_2050",meson.Data()));
//             graphCombInvYieldSysPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sCombPbPb2760GeVSysErrNoShift_2050",meson.Data()));
//
//             if(graphPCMInvYieldStatPbPb2760GeVUnShifted_0010)graphPCMInvYieldStatPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sPCMPbPb2760GeVStatErrNoShift_0010",meson.Data()));
//             if(graphPCMInvYieldSysPbPb2760GeVUnShifted_0010)graphPCMInvYieldSysPbPb2760GeVUnShifted_0010->Write(Form("graphInvYield%sPCMPbPb2760GeVSysErrNoShift_0010",meson.Data()));
//             if(graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010)graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010->Write(Form("graphInvYield%sEMCalPbPb2760GeVStatErrNoShift_0010",meson.Data()));
//             if(graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010)graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010->Write(Form("graphInvYield%sEMCalPbPb2760GeVSysErrNoShift_0010",meson.Data()));
//             if(graphPHOSInvYieldStatPbPb2760GeVUnshifted_0010)graphPHOSInvYieldStatPbPb2760GeVUnshifted_0010->Write(Form("graphInvYield%sPHOSPbPb2760GeVStatErrNoShift_0010",meson.Data()));
//             if(graphPHOSInvYieldSysPbPb2760GeVUnshifted_0010)graphPHOSInvYieldSysPbPb2760GeVUnshifted_0010->Write(Form("graphInvYield%sPHOSPbPb2760GeVSysErrNoShift_0010",meson.Data()));
//             if(graphPCMInvYieldStatPbPb2760GeVUnShifted_2050)graphPCMInvYieldStatPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sPCMPbPb2760GeVStatErrNoShift_2050",meson.Data()));
//             if(graphPCMInvYieldSysPbPb2760GeVUnShifted_2050)graphPCMInvYieldSysPbPb2760GeVUnShifted_2050->Write(Form("graphInvYield%sPCMPbPb2760GeVSysErrNoShift_2050",meson.Data()));
//             if(graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050)graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050->Write(Form("graphInvYield%sEMCalPbPb2760GeVStatErrNoShift_2050",meson.Data()));
//             if(graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050)graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050->Write(Form("graphInvYield%sEMCalPbPb2760GeVSysErrNoShift_2050",meson.Data()));
//
//             if(graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_0010)graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_0010->Write(Form("UnshiftedSysforRAAPCM%s_0010",meson.Data()));
//             if(graphPHOSInvYieldSysPbPb2760GeVforRAAUnshifted_0010)graphPHOSInvYieldSysPbPb2760GeVforRAAUnshifted_0010->Write(Form("UnshiftedSysforRAAPHOS%s_0010",meson.Data()));
//             if(graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_0010)graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_0010->Write(Form("UnshiftedSysforRAAEMCal%s_0010",meson.Data()));
//             if(graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_2050)graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_2050->Write(Form("UnshiftedSysWOMatPCM%s_2050",meson.Data()));
//             if(graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_2050)graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_2050->Write(Form("UnshiftedSysforRAAEMCal%s_2050",meson.Data()));
//
//
//     fCombResults->Close();

}