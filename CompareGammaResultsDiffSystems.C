#include <Riostream.h>
#include <fstream>
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
#include "TGaxis.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/ExtractSignalBinning.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

Color_t GetDefaultDifferentSystemsAndCents( TString centrality,
                                            TString collisionsEnergy,
                                            TString experiment          = "ALICE"
                                          ){


    if (!collisionsEnergy.CompareTo("PbPb_2.76TeV") && !experiment.CompareTo("ALICE")){
        if (!centrality.CompareTo("0-20%")){
            return kRed+1;
        } else if (!centrality.CompareTo("20-40%")){
            return kRed-3;
        } else if (!centrality.CompareTo("40-80%")){
            return kRed-5;
        }
    } else if (!collisionsEnergy.CompareTo("pPb_5.023TeV") && !experiment.CompareTo("ALICE")){
        if (!centrality.CompareTo("0-20%")){
            return kRed+4;
        } else if (!centrality.CompareTo("20-40%")){
            return kRed+3;
        } else if (!centrality.CompareTo("40-60%")){
            return kRed-2;
        } else if (!centrality.CompareTo("60-100%")){
            return kRed-3;
        } else if (!centrality.CompareTo("0-100%")){
            return kRed+2;
        }
    } else if (!collisionsEnergy.CompareTo("pPb_8.16TeV") && !experiment.CompareTo("ATLAS")){
        if (!centrality.CompareTo("0-100%")){
            return kOrange+9;
        }
    } else if (!collisionsEnergy.CompareTo("PbPb_2.76TeV") && !experiment.CompareTo("CMS")){
        if (!centrality.CompareTo("0-10%")){
            return kRed+3;
        } else if (!centrality.CompareTo("10-30%")){
            return kRed+2;
        } else if (!centrality.CompareTo("30-100%")){
            return kRed-2;
        }
    } else if (!collisionsEnergy.CompareTo("AuAu_0.2TeV") && !experiment.CompareTo("PHENIX")){
        if (!centrality.CompareTo("0-20%")){
            return kGreen+2;
        } else if (!centrality.CompareTo("20-40%")){
            return kGreen-2;
        } else if (!centrality.CompareTo("40-60%")){
            return kGreen+3;
        } else if (!centrality.CompareTo("60-92%")){
            return kGreen-1;
        } else if (!centrality.CompareTo("0-92%")){
            return kGreen-5;
        }
    } else if (!collisionsEnergy.CompareTo("pAu_0.2TeV") && !experiment.CompareTo("PHENIX")){
        if (!centrality.CompareTo("0-20%")){
            return kGreen+2;
        } else if (!centrality.CompareTo("20-40%")){
            return kGreen-2;
        } else if (!centrality.CompareTo("40-60%")){
            return kGreen+3;
        } else if (!centrality.CompareTo("60-100%")){
            return kGreen-1;
        } else if (!centrality.CompareTo("0-100%")){
            return kGreen-5;
        }
    } else if (!collisionsEnergy.CompareTo("AuAu_0.2TeV") && !experiment.CompareTo("STAR")){
        if (!centrality.CompareTo("0-20%")){
            return kYellow+2;
        } else if (!centrality.CompareTo("20-40%")){
            return kYellow-2;
        } else if (!centrality.CompareTo("40-60%")){
            return kYellow+3;
        } else if (!centrality.CompareTo("60-80%")){
            return kYellow-1;
        } else if (!centrality.CompareTo("0-80%")){
            return kYellow-5;
        }
    } else if (!collisionsEnergy.CompareTo("CuCu_0.2TeV") && !experiment.CompareTo("PHENIX")){
        if (!centrality.CompareTo("0-40%")){
            return kTeal-6;
        } else if (!centrality.CompareTo("0-94%")){
            return kTeal-1;
        }
    } else if (!collisionsEnergy.CompareTo("dAu_0.2TeV") && !experiment.CompareTo("PHENIX")){
        if (!centrality.CompareTo("0-100%")){
            return kTeal-6;
        } else if (!centrality.CompareTo("0-94%")){
            return kTeal-1;
        }
    } else if (!collisionsEnergy.CompareTo("AuAu_62.4GeV") && !experiment.CompareTo("PHENIX")){
        if (!centrality.CompareTo("0-20%")){
            return kAzure+2;
        } else if (!centrality.CompareTo("20-40%")){
            return kAzure-6;
        } else if (!centrality.CompareTo("0-86%")){
            return kAzure+1;
        }
    } else if (!collisionsEnergy.CompareTo("AuAu_39GeV") && !experiment.CompareTo("PHENIX")){
        if (!centrality.CompareTo("0-86%")){
            return kViolet-5;
        }
    } else if (!collisionsEnergy.CompareTo("PbPb_17.2GeV") && !experiment.CompareTo("WA98")){
        if (!centrality.CompareTo("0-10%")){
            return kViolet+3;
        }
    }

}


void CompareGammaResultsDiffSystems(    TString inputFileNamePP2760GeV      = "",
                                        TString inputFileNamePP7TeV         = "",
                                        TString inputFileNamePP8TeV         = "",
                                        TString inputFileNamePPb5TeV        = "",
                                        TString inputFileNamePbPb2760GeV    = "",
                                        TString suffix                      = "eps",
                                        Bool_t enablePrelim                   = kFALSE
                                    ){



    //*******************************************************************************************************************************************
    //*********************************************************** Set main stylabelScalingDirGammale choices ********************************************************
    //*******************************************************************************************************************************************
    StyleSettingsThesis();
    SetPlotStyle();

    //*******************************************************************************************************************************************
    //********************************************* Create output directory and copy input files there ******************************************
    //*******************************************************************************************************************************************
    TString dateForOutput                                       = ReturnDateStringForOutput();
    TString outputDir                                           = Form("%s/%s/CompareDiffSystems",suffix.Data(),dateForOutput.Data());

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputGammapp2760GeV.root", inputFileNamePP2760GeV.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputGammapp8TeV.root", inputFileNamePP8TeV.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputGammapp7TeV.root", inputFileNamePP7TeV.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputGammapPb5TeV.root", inputFileNamePPb5TeV.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputGammaPbPb2760GeV.root", inputFileNamePbPb2760GeV.Data(), outputDir.Data()));

    TString fileNameTheoryPP                                    = "ExternalInput/Theory/TheoryCompilationPP.root";
    TString fileNameTheoryPPb                                   = "ExternalInputPPb/Theory/TheoryCompilationPPb.root";
    TString fileNameTheoryPbPb                                  = "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";
    TString fileNameExperimentPP                                = "ExternalInput/OtherExperiments/DataCompilationFromOtherEnergiesPP.root";
    TString fileNameExperimentPA                                = "ExternalInputpPb/OtherExperiments/DataCompilationFromOtherEnergiesPA.root";
    TString fileNameExperimentAA                                = "ExternalInputPbPb/OtherExperiments/DataCompilationFromOtherEnergiesAA.root";
    //*******************************************************************************************************************************************
    //******************************************************* set ranges for plotting ***********************************************************
    //*******************************************************************************************************************************************
    Double_t doubleRatio[2];
    Double_t doubleRatioX[2];
    doubleRatio[0]      = 0.75;     doubleRatio[1]      = 1.65;
    doubleRatioX[0]     = 0.23;     doubleRatioX[1]     = 25.5;

    Color_t colorCocktailPi0                        = kRed+2;
    Color_t colorCocktailEta                        = kBlue+1;
    Color_t colorCocktailEtaP                       = kOrange+1;
    Color_t colorCocktailOmega                      = kYellow+2;
    Color_t colorCocktailPhi                        = kViolet;
    Color_t colorCocktailRho0                       = kAzure-2;
    Color_t colorCocktailSigma0                     = kGray+1;

    Color_t  colorNLOWerner                         = kAzure+2;
    Color_t  colorNLOWernerBand                     = kAzure-9;
    Color_t  colorNLOMcGill                         = kGreen+2;
    Color_t  colorNLOMcGillBand                     = kGreen-6;
    Style_t  styleMarkerNLOWerner                   = 24;
    Style_t  styleLineNLOWerner                     = 5;
    Style_t  styleLineMcGill                        = 7;
    Width_t  widthLineNLO                           = 2.;

    // pPb settings
    TString namepPbCent[5]                         = {"0-20%", "20-40%", "40-60%", "60-100%", "0-100%"};
    TString namepPbCentOut[5]                      = {"0020", "2040", "4060", "60100", "0100"};
    TString namepPbCentForNColl[5]                      = {"0020_V0A", "2040_V0A", "4060_V0A", "60100_V0A", "0100_V0A"};
    TString labelXTALICEpPb[5]                         = {"ALICE prel. (0-20% p-Pb 5.02TeV)", "ALICE prel. (20-40% p-Pb 5.02TeV)", "ALICE prel. (40-60% p-Pb 5.02TeV)", "ALICE prel. (60-100% p-Pb 5.02TeV)", "ALICE prel. (0-100% p-Pb 5.02TeV)"};
    TString labelALICEpPb[5]                           = {"ALICE prel. (0-20%)", "ALICE prel. (20-40%)", "ALICE prel. (40-60%)", "ALICE prel. (60-100%)", "ALICE prel. (0-100%)"};
    Double_t xLegALICEXTpPb[5]                         = {0.43, 0.445, 0.46, 0.475, 0.49};
    Double_t yLegALICEXTpPb[5]                         = {0.88, 0.86, 0.84, 0.82, 0.80};
    TString nameLabelScalepPb[5]                   = {"x 10^{1}", "x 10^{0}", "x 10^{-2}", "x 10^{-3}", "x 10^{-4}"};
    Double_t scaleFacPlotApPb[5]                            = {10e8, 10e6, 10e4, 10e2, 1};


    Color_t colorCombpPb[5];
    Color_t colorCombpPbBox[5];
    Style_t markerStyleCombpPb[5];
    Style_t markerStyleCombpPbMC[5];
    Size_t markerSizeCombpPb[5];
    Double_t nCollpPb5[5];
    for (Int_t ncent = 0; ncent<5; ncent++){
        colorCombpPb[ncent]                            = GetDefaultDifferentSystemsAndCents(namepPbCent[ncent].Data(),"pPb_5.023TeV", "ALICE");
        colorCombpPbBox[ncent]                         = GetColorDefaultColor("pPb_5.023TeV", "", namepPbCent[ncent].Data(), kTRUE);
        markerStyleCombpPb[ncent]                      = GetDefaultMarkerStyle("pPb_5.023TeV", "", namepPbCent[ncent].Data());
        markerStyleCombpPbMC[ncent]                    = GetDefaultMarkerStyle("pPb_5.023TeV", "MC", namepPbCent[ncent].Data());
        markerSizeCombpPb[ncent]                       = GetDefaultMarkerSize("pPb_5.023TeV", "", namepPbCent[ncent].Data());
        nCollpPb5[ncent]                               = GetNCollFromName(namepPbCentForNColl[ncent].Data(),"pPb_5.023TeV");
    }

    TString nameCentspPb8TeVATLAS                         = "0-100%";
    TString nameCentspPb8TeVATLASOut                      = "0100";
    TString labelXTpPb8TeVATLAS                         = "ATLAS prel. (0-100% p-Pb 8.16TeV)";
    TString labelpPb8TeVATLAS                         = "ATLAS prel. (0-100%)";
    Double_t xLegpPb8TeVATLASXT                         = 0.65;
    Double_t yLegpPb8TeVATLASXT                         = 0.58;
    Color_t colorCombpPb8ATLAS   = GetDefaultDifferentSystemsAndCents(nameCentspPb8TeVATLAS.Data(),"pPb_8.16TeV", "ATLAS");
    Color_t colorCombpPb8ATLASBox= GetDefaultDifferentSystemsAndCents(nameCentspPb8TeVATLAS.Data(),"pPb_8.16TeV", "ATLAS");
    Style_t markerStyleCombpPb8ATLAS= GetDefaultMarkerStyle("pPb_5.023TeV", "", nameCentspPb8TeVATLAS.Data());
    Style_t markerStyleCombpPb8ATLASMC= GetDefaultMarkerStyle("pPb_5.023TeV", "MC", nameCentspPb8TeVATLAS.Data());
    Size_t markerSizeCombpPb8ATLAS =GetDefaultMarkerSize("pPb_5.023TeV", "", nameCentspPb8TeVATLAS.Data());
    Double_t nCollpPb8ATLAS = 7.09;

    // PbPb settings
    TString namePbPbCent[3]                         = {"0-20%", "20-40%", "40-80%"};
    TString namePbPbCentOut[3]                      = {"0020", "2040", "4080"};
    TString labelXTALICE[3]                         = {"ALICE (0-20% Pb-Pb 2.76TeV)", "ALICE (20-40% Pb-Pb 2.76TeV)", "ALICE (40-80% Pb-Pb 2.76TeV)"};
    TString labelALICE[3]                           = {"ALICE (0-20%)", "ALICE (20-40%)", "ALICE (40-80%)"};
    Double_t xLegALICEXT[3]                         = {0.43, 0.445, 0.46};
    Double_t yLegALICEXT[3]                         = {0.88, 0.86, 0.84};
    TString nameLabelScalePbPb[3]                   = {"x 10^{1}", "x 10^{0}", "x 10^{-2}"};
    Double_t scaleFacPlotAPbPb[3]                   = {10e1, 1, 10e-2};

    Color_t colorComb[3];
    Color_t colorCombPbALICE[3];
    Style_t markerStyleComb[3];
    Style_t markerStyleCombMC[3];
    Size_t markerSizeComb[3];
    Double_t nColl[3];
    Color_t colorCombAr[3]                          = {kRed-5, kYellow-5, kAzure-7};
    for (Int_t centPb = 0; centPb < 3; centPb++){
        colorComb[centPb]                           = GetColorDefaultColor("PbPb_2.76TeV", "", namePbPbCent[centPb].Data());
        colorCombPbALICE[centPb]                    = GetDefaultDifferentSystemsAndCents(namePbPbCent[centPb].Data(), "PbPb_2.76TeV", "ALICE");
        markerStyleComb[centPb]                     = GetDefaultMarkerStyle("PbPb_2.76TeV", "", namePbPbCent[centPb].Data());
        markerStyleCombMC[centPb]                   = GetDefaultMarkerStyle("PbPb_2.76TeV", "MC", namePbPbCent[centPb].Data());
        markerSizeComb[centPb]                      = GetDefaultMarkerSize("PbPb_2.76TeV", "", namePbPbCent[centPb].Data());
        nColl[centPb]                               = GetNCollFromName(namePbPbCentOut[centPb].Data());
    }
    TString namePbPbCentCMS[3]                      = {"0-10%", "10-30%", "30-100%"};
    TString namePbPbCentCMSOut[3]                   = {"0010", "1030", "30100"};
    TString labelXTCMS[3]                           = {"CMS (0-10% Pb-Pb 2.76TeV)", "CMS (10-30% Pb-Pb 2.76TeV)", "CMS (30-100% Pb-Pb 2.76TeV)"};
    TString labelCMS[3]                             = {"CMS (0-10%)", "CMS (10-30%)", "CMS (30-100%)"};
    Double_t xLegCMSXT[3]                           = {0.25, 0.265, 0.28};
    Double_t yLegCMSXT[3]                           = {0.56, 0.54, 0.52};

    Color_t colorCombCMS[3];
    Style_t markerStyleCombCMS[3];
    Style_t markerStyleCombMCCMS[3];
    Size_t markerSizeCombCMS[3];
    Double_t nCollCMS[3];
    for (Int_t centPb = 0; centPb < 3; centPb++){
        colorCombCMS[centPb]                           = GetDefaultDifferentSystemsAndCents(namePbPbCentCMS[centPb].Data(), "PbPb_2.76TeV", "CMS");
        markerStyleCombCMS[centPb]                     = GetDefaultMarkerStyle("PbPb_2.76TeV", "", namePbPbCentCMS[centPb].Data());
        markerStyleCombMCCMS[centPb]                   = GetDefaultMarkerStyle("PbPb_2.76TeV", "MC", namePbPbCentCMS[centPb].Data());
        markerSizeCombCMS[centPb]                      = GetDefaultMarkerSize("PbPb_2.76TeV", "", namePbPbCentCMS[centPb].Data());
        nCollCMS[centPb]                               = GetNCollFromName(namePbPbCentCMSOut[centPb].Data());
    }

    Color_t colorCombpAu200[2];
    Style_t markerStyleCombpAu200[2];
    Style_t markerStyleCombMCpAu200[2];
    Size_t markerSizeCombpAu200[2];
    Double_t nCollpAu200[2];
    TString nameCentspAu200Out[2]                  = {"0020", "0100"};
    TString nameCentspAu200[2]                     = {"0-20%", "0-100%"};
    TString labelXTPHENIXpAu200[2]                     = {"PHENIX prel. (0-20% p-Au 0.2TeV)", "PHENIX prel. (0-100% p-Au 0.2TeV)"};
    TString labelPHENIXpAu200[2]                       = {"PHENIX prel. (0-20%)", "PHENIX prel. (0-100%)"};
    Double_t xLegPHENIXpAu200XT[2]                     = {0.33, 0.345};
    Double_t yLegPHENIXpAu200XT[2]                     = {0.43, 0.41};
    for (Int_t centAu = 0; centAu < 2; centAu++){
        colorCombpAu200[centAu]                    = GetDefaultDifferentSystemsAndCents(nameCentspAu200[centAu].Data(), "pAu_0.2TeV", "PHENIX");
        markerStyleCombpAu200[centAu]              = GetDefaultMarkerStyle("pPb_5.023TeV", "", nameCentspAu200[centAu].Data());
        markerStyleCombMCpAu200[centAu]            = GetDefaultMarkerStyle("pPb_5.023TeV", "MC", nameCentspAu200[centAu].Data());
        markerSizeCombpAu200[centAu]               = GetDefaultMarkerSize("pPb_5.023TeV", "", nameCentspAu200[centAu].Data());
        nCollpAu200[centAu]                        = GetNCollFromName(nameCentspAu200Out[centAu].Data(), "pAu_0.2TeV");
    }


    Color_t colorCombdAu200;
    Style_t markerStyleCombdAu200;
    Style_t markerStyleCombMCdAu200;
    Size_t markerSizeCombdAu200;
    Double_t nColldAu200;
    TString nameCentsdAu200Out                  = "0100";
    TString nameCentsdAu200                     = "0-100%";
    TString labelXTPHENIXdAu200                     = "PHENIX (0-100% d-Au 0.2TeV)";
    TString labelPHENIXdAu200                       = "PHENIX (0-100%)";
    Double_t xLegPHENIXdAu200XT                     = 0.395;
    Double_t yLegPHENIXdAu200XT                     = 0.37;
        colorCombdAu200                    = GetDefaultDifferentSystemsAndCents(nameCentsdAu200.Data(), "dAu_0.2TeV", "PHENIX");
        markerStyleCombdAu200              = GetDefaultMarkerStyle("PbPb_2.76TeV", "", nameCentsdAu200.Data());
        markerStyleCombMCdAu200            = 28;
        markerSizeCombdAu200               = GetDefaultMarkerSize("PbPb_2.76TeV", "", nameCentsdAu200.Data());
        nColldAu200                        = GetNCollFromName(nameCentsdAu200Out.Data(), "dAu_0.2TeV");

    Color_t colorCombAuAu200[5];
    Style_t markerStyleCombAuAu200[5];
    Style_t markerStyleCombMCAuAu200[5];
    Size_t markerSizeCombAuAu200[5];
    Double_t nCollAuAu200[5];
    TString nameCentsAuAu200Out[5]                  = {"0020", "2040", "4060", "6092", "0092"};
    TString nameCentsAuAu200[5]                     = {"0-20%", "20-40%", "40-60%", "60-92%", "0-92%"};
    TString labelXTPHENIX200[5]                     = {"PHENIX (0-20% Au-Au 0.2TeV)", "PHENIX (20-40% Au-Au 0.2TeV)", "PHENIX (40-60% Au-Au 0.2TeV)", "PHENIX (60-92% Au-Au 0.2TeV)", "PHENIX (0-92% Au-Au 0.2TeV)"};
    TString labelPHENIX200[5]                       = {"PHENIX (0-20%)", "PHENIX (20-40%)", "PHENIX (40-60%)", "PHENIX (60-92%)", "PHENIX (0-92%)"};
    Double_t xLegPHENIX200XT[5]                     = {0.535, 0.555, 0.575, 0.595, 0.615};
    Double_t yLegPHENIX200XT[5]                     = {0.78, 0.76, 0.74, 0.72, 0.70};
    for (Int_t centAu = 0; centAu < 5; centAu++){
        colorCombAuAu200[centAu]                    = GetDefaultDifferentSystemsAndCents(nameCentsAuAu200[centAu].Data(), "AuAu_0.2TeV", "PHENIX");
        markerStyleCombAuAu200[centAu]              = GetDefaultMarkerStyle("PbPb_2.76TeV", "", nameCentsAuAu200[centAu].Data());
        markerStyleCombMCAuAu200[centAu]            = GetDefaultMarkerStyle("PbPb_2.76TeV", "MC", nameCentsAuAu200[centAu].Data());
        markerSizeCombAuAu200[centAu]               = GetDefaultMarkerSize("PbPb_2.76TeV", "", nameCentsAuAu200[centAu].Data());
        nCollAuAu200[centAu]                        = GetNCollFromName(nameCentsAuAu200Out[centAu].Data(), "AuAu_0.2TeV");
    }

    Color_t colorCombAuAu200STAR[5];
    Style_t markerStyleCombAuAu200STAR[5];
    Style_t markerStyleCombMCAuAu200STAR[5];
    Size_t markerSizeCombAuAu200STAR[5];
    Double_t nCollAuAu200STAR[5];
    TString nameCentsAuAu200STAROut[5]                  = {"0020", "2040", "4060", "6080", "0080"};
    TString nameCentsAuAu200STAR[5]                     = {"0-20%", "20-40%", "40-60%", "60-80%", "0-80%"};
    TString labelXTSTAR200[5]                           = {"STAR (0-20% Au-Au 0.2TeV)", "STAR (20-40% Au-Au 0.2TeV)", "STAR (40-60% Au-Au 0.2TeV)", "STAR (60-80% Au-Au 0.2TeV)", "STAR (0-80% Au-Au 0.2TeV)"};
    TString labelSTAR200[5]                             = {"STAR (0-20%)", "STAR (20-40%)", "STAR (40-60%)", "STAR (60-80%)", "STAR (0-80%)"};
    Double_t xLegSTAR200XT[5]                           = {0.33, 0.345, 0.36, 0.375, 0.39};
    Double_t yLegSTAR200XT[5]                           = {0.43, 0.41, 0.39, 0.37, 0.35};
    for (Int_t centAu = 0; centAu < 5; centAu++){
        colorCombAuAu200STAR[centAu]                    = GetDefaultDifferentSystemsAndCents(nameCentsAuAu200STAR[centAu].Data(), "AuAu_0.2TeV", "STAR");
        markerStyleCombAuAu200STAR[centAu]              = GetDefaultMarkerStyle("PbPb_2.76TeV", "", nameCentsAuAu200STAR[centAu].Data());
        markerStyleCombMCAuAu200STAR[centAu]            = GetDefaultMarkerStyle("PbPb_2.76TeV", "MC", nameCentsAuAu200STAR[centAu].Data());
        markerSizeCombAuAu200STAR[centAu]               = GetDefaultMarkerSize("PbPb_2.76TeV", "", nameCentsAuAu200STAR[centAu].Data());
        nCollAuAu200STAR[centAu]                        = GetNCollFromName(nameCentsAuAu200STAROut[centAu].Data(), "AuAu_0.2TeV_STAR");
    }


    Color_t colorCombAuAu62[3];
    Style_t markerStyleCombAuAu62[3];
    Style_t markerStyleCombMCAuAu62[3];
    Size_t markerSizeCombAuAu62[3];
    Double_t nCollAuAu62[3];
    TString nameCentsAuAu62Out[3]                  = {"0020", "2040", "0086"};
    TString nameCentsAuAu62[3]                     = {"0-20%", "20-40%", "0-86%"};
    TString labelXTPHENIX62[3]                     = {"PHENIX (0-20% Au-Au 62.4GeV)", "PHENIX (20-40% Au-Au 62.4GeV)", "PHENIX (0-86% Au-Au 62.4GeV)"};
    TString labelPHENIXAu62[3]                     = {"PHENIX (0-20%)", "PHENIX (20-40%)", "PHENIX (0-86%)"};
    Double_t xLegPHENIX62XT[3]                     = {0.665, 0.675, 0.685};
    Double_t yLegPHENIX62XT[3]                     = {0.62, 0.60, 0.58};
    for (Int_t centAu = 0; centAu < 3; centAu++){
        colorCombAuAu62[centAu]                    = GetDefaultDifferentSystemsAndCents(nameCentsAuAu62[centAu].Data(), "AuAu_62.4GeV", "PHENIX");
        markerStyleCombAuAu62[centAu]              = GetDefaultMarkerStyle("PbPb_2.76TeV", "", nameCentsAuAu62[centAu].Data());
        markerStyleCombMCAuAu62[centAu]            = GetDefaultMarkerStyle("PbPb_2.76TeV", "MC", nameCentsAuAu62[centAu].Data());
        markerSizeCombAuAu62[centAu]               = GetDefaultMarkerSize("PbPb_2.76TeV", "", nameCentsAuAu62[centAu].Data());
        nCollAuAu62[centAu]                        = GetNCollFromName(nameCentsAuAu62Out[centAu].Data(), "AuAu_62.4GeV");
    }


    TString nameCentsAuAu39                         = "0-86%";
    TString nameCentsAuAu39Out                      = "0086";
    TString labelXTPHENIX39                         = "PHENIX (0-86% Au-Au 39GeV)";
    TString labelPHENIXAu39                         = "PHENIX (0-86%)";
    Double_t xLegPHENIX39XT                         = 0.71;
    Double_t yLegPHENIX39XT                         = 0.55;
    Double_t nCollAuAu39                            = GetNCollFromName("0086", "AuAu_39GeV");
    Color_t colorCombAuAu39                         = GetDefaultDifferentSystemsAndCents(nameCentsAuAu39.Data(), "AuAu_39GeV", "PHENIX");
    Style_t markerStyleCombAuAu39                   = kStar;
    Style_t markerStyleCombMCAuAu39                 = 28;
    Size_t markerSizeCombAuAu39                     = 2.5;

    TString nameCentsPbPb17                         = "0-10%";
    TString nameCentsPbPb17Out                      = "0010";
    TString labelXTWA98                             = "WA98 (0-10% Pb-Pb 17.2GeV)";
    TString labelWA98                               = "WA98 (0-10%)";
    Double_t xLegWA98XT                             = 0.49;
    Double_t yLegWA98XT                             = 0.20;
    Double_t nCollPbPb17                            = GetNCollFromName("0010", "PbPb_17.2GeV");
    Color_t colorCombPbPb17                         = GetDefaultDifferentSystemsAndCents(nameCentsPbPb17.Data(), "PbPb_17.2GeV", "WA98");
    Style_t markerStyleCombPbPb17                   = kStar;
    Style_t markerStyleCombMCPbPb17                 = kOpenStar;
    Size_t markerSizeCombPbPb17                     = 2.5;

    Color_t colorCombCuCu200[2];
    Style_t markerStyleCombCuCu200[2];
    Style_t markerStyleCombMCCuCu200[2];
    Size_t markerSizeCombCuCu200[2];
    Double_t nCollCuCu200[2];
    TString nameCentsCuCu200Out[2]                  = {"0040", "0094"};
    TString nameCentsCuCu200[2]                     = {"0-40%", "0-94%"};
    TString labelXTPHENIXCuCu200[2]                 = {"PHENIX (0-40% Cu-Cu 0.2TeV)", "PHENIX (0-92% Cu-Cu 0.2TeV)"};
    TString labelPHENIXCu200[2]                     = {"PHENIX (0-40%)", "PHENIX (0-92%)"};
    Double_t xLegPHENIXCuCu200CXT[2]                = {0.635, 0.65};
    Double_t yLegPHENIXCuCu200XT[2]                 = {0.67, 0.65};
    for (Int_t centCu = 0; centCu < 2; centCu++){
        colorCombCuCu200[centCu]                    = GetDefaultDifferentSystemsAndCents(nameCentsCuCu200[centCu].Data(), "CuCu_0.2TeV", "PHENIX");
        markerStyleCombCuCu200[centCu]              = GetDefaultMarkerStyle("PbPb_2.76TeV", "", nameCentsCuCu200[centCu].Data());
        markerStyleCombMCCuCu200[centCu]            = GetDefaultMarkerStyle("PbPb_2.76TeV", "MC", nameCentsCuCu200[centCu].Data());
        markerSizeCombCuCu200[centCu]               = GetDefaultMarkerSize("PbPb_2.76TeV", "", nameCentsCuCu200[centCu].Data());
        nCollCuCu200[centCu]                        = GetNCollFromName(nameCentsCuCu200Out[centCu].Data(), "CuCu_0.2TeV");
    }

    TString collisionSystempp2760GeV                = "pp, #sqrt{#it{s}} = 2.76 TeV";
    TString collisionSystempp8TeV                   = "pp, #sqrt{#it{s}} = 8 TeV";
    TString collisionSystempp7TeV                   = "pp, #sqrt{#it{s}} = 7 TeV";
    TString collisionSystemPbPb2760GeV               = "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystempAu200GeV                = "p-Au, #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";
    TString collisionSystemdAu200GeV                = "d-Au, #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";
    TString collisionSystemAuAu200GeV               = "Au-Au, #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";
    TString collisionSystemAuAu39GeV                = "Au-Au, #sqrt{#it{s}_{_{NN}}} = 39 GeV";
    TString collisionSystemAuAu62GeV                = "Au-Au, #sqrt{#it{s}_{_{NN}}} = 62.4 GeV";
    TString collisionSystemCuCu200GeV               = "Cu-Cu, #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";
    TString collisionSystemPbPb17GeV                = "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 17.2 GeV";
    TString collisionSystempPb5TeV                  = "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    TString collisionSystempPb8TeV                  = "p-Pb, #sqrt{#it{s}_{_{NN}}} = 8.16 TeV";
    TString textALICE                               = "ALICE";

    Color_t colorCombpp2760GeV                      = kBlack;
    Style_t markerStyleCombpp2760GeV                = 20;
    Size_t markerSizeCombpp2760GeV                  = 1.8;
    Width_t widthLinesBoxes                         = 1.4;
    Width_t widthCommonFit                          = 2.4;

    Color_t colorEnergy[6]                          = {kRed+1, kGreen+2, kBlue+2, kRed+3, kYellow+2, kViolet-2 };

    Color_t colorEpp2760GeV                         = GetColorDefaultColor("2.76TeV", "", "");
    Style_t markerStyleEpp2760GeV                   = GetDefaultMarkerStyle("2.76TeV", "", "");
    Size_t markerSizeEpp2760GeV                     = GetDefaultMarkerSize("2.76TeV", "", "");

    Color_t colorEpp8TeV                            = GetColorDefaultColor("8TeV", "", "");
    Style_t markerStyleEpp8TeV                      = GetDefaultMarkerStyle("8TeV", "", "");
    Size_t markerSizeEpp8TeV                        = GetDefaultMarkerSize("8TeV", "", "");

    Color_t colorEpp7TeV                            = GetColorDefaultColor("7TeV", "", "");
    Style_t markerStyleEpp7TeV                      = GetDefaultMarkerStyle("7TeV", "", "");
    Size_t markerSizeEpp7TeV                        = GetDefaultMarkerSize("7TeV", "", "")*0.8;

    const Int_t nOtherEnergies                            = 23;
    TString nameOtherExp[nOtherEnergies]            = { "ATLAS","ATLAS", "ALICE prel.", "ATLAS", "ATLAS", "CMS", "CMS",
                                                        "CMS", "CDF", "D0", "UA1", "UA2",
                                                        "UA1", "PHENIX", "R807", "R807",
                                                        "E706", "E706", "NA24", "WA70",
                                                        "UA6", "UA6" , "E704" };
    TString nameOtherExpRead[nOtherEnergies]       = { "ATLAS","ATLAS", "ALICE", "ATLAS_1", "ATLAS", "CMS_1", "CMS",
                                                        "CMS", "CDF", "D0", "UA1", "UA2",
                                                        "UA1", "PHENIX", "R807", "R807",
                                                        "E706", "E706", "NA24", "WA70",
                                                        "UA6_0", "UA6_1", "E704" };
    TString nameOtherEnergyRead[nOtherEnergies]    = { "13TeV", "8TeV", "7TeV", "7TeV", "7TeV", "7TeV", "7TeV",
                                                        "2.76TeV",  "1.8TeV", "1.8TeV", "630GeV", "630GeV",
                                                        "546GeV", "200GeV", "63GeV", "62.4GeV",
                                                        "38.8GeV", "31.6GeV", "23.8GeV", "23GeV",
                                                        "24.3GeV", "24.3GeV", "19.4GeV" };
    TString nameOtherEnergy[nOtherEnergies]        = { "13 TeV","8 TeV", "7 TeV", "7 TeV", "7 TeV", "7 TeV", "7 TeV",
                                                        "2.76 TeV", "1.8 TeV", "1.8 TeV", "630 GeV", "630 GeV",
                                                        "546 GeV", "200 GeV", "63 GeV", "62.4 GeV",
                                                        "38.8 GeV", "31.6 GeV", "23.8 GeV", "23 GeV",
                                                        "24.3 GeV", "24.3 GeV", "19.4 GeV"};
    Color_t colorOtherExp[nOtherEnergies]          = {  kViolet-3, kGreen+2, kBlue+3, kBlue+1, kBlue+2, kBlue+2, kBlue+1,    kMagenta+1, kCyan+2, kCyan+2, kOrange+7, kOrange+7,
                                                        kOrange+7, kViolet+2, kGray, kGray, kGray+1,     kGray+1, kGray+2, kGray+2, kGray+2,
                                                        kGray+2, 1 };
    Style_t markerStyleOtherExp[nOtherEnergies]     = { 21, 27, 25, 24, 28, 30, 25, 24, 34, 33, 29, 25,
                                                        24, 20, 28, 27, 30, 25, 20, 21, 34, 33,
                                                        24  };
    Size_t markerSizeOtherExp[nOtherEnergies]       = { 1.3, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
                                                        1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
                                                        1.5, 1.5};
    Double_t legendOtherExtX[nOtherEnergies]        = { 0.60, 0.65, 0.30, 0.404, 0.382, 0.446, 0.426,   0.7, 0.74, 0.76, 0.5, 0.52,
                                                        0.54, 0.79, 0.57, 0.59, 0.82,       0.83, 0.62, 0.64, 0.66, 0.67,
                                                        0.84  };
    Double_t legendOtherExtY[nOtherEnergies]        = { 0.60, 0.55, 0.58, 0.52, 0.55, 0.46, 0.49,       0.50, 0.45, 0.42, 0.40, 0.37,
                                                        0.34, 0.37, 0.30, 0.27, 0.32,       0.29, 0.23, 0.20, 0.17, 0.14,
                                                        0.26 };
    //*************************************************************************************************************************************************
    //*************************************** Read in data ********************************************************************************************
    //*************************************************************************************************************************************************
    //--------------------------------------- pp 2.76TeV --------------------------------------------
    TFile* fileCombPP2760GeV                                = new TFile( inputFileNamePP2760GeV.Data());
    TGraphAsymmErrors* graphDRStatpp2760GeV                 = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphRGammaNonFitCombStatErr");
    TGraphAsymmErrors* graphDRSyspp2760GeV                  = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphRGammaNonFitCombSysErr");
    TGraphAsymmErrors* graphDRTotpp2760GeV                  = AddErrorsQuadraticallyTGraph(graphDRStatpp2760GeV, graphDRSyspp2760GeV);
    TGraphAsymmErrors* graphInvYieldDirGammaStatpp2760GeV   = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldDirGammaNonFitStatErr");
    TGraphAsymmErrors* graphInvYieldDirGammaSyspp2760GeV    = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldDirGammaNonFitSysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpp2760GeV  = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldDirGammaNonFitSumErrAr");
    TGraphAsymmErrors* graphXSecDirGammaStatpp2760GeV       = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvXSectionDirGammaNonFitStatErr");
    TGraphAsymmErrors* graphXSecDirGammaSyspp2760GeV        = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvXSectionDirGammaNonFitSysErr");
    TGraphAsymmErrors* graphXSecDirGammaTotArpp2760GeV      = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvXSectionDirGammaNonFitSumErrAr");
    TGraphAsymmErrors* graphXSecDirGammaStatpp2760GeV_xT    = xTScalePhoton(graphXSecDirGammaStatpp2760GeV, 2760.);
    TGraphAsymmErrors* graphXSecDirGammaSyspp2760GeV_xT     = xTScalePhoton(graphXSecDirGammaSyspp2760GeV, 2760.);
    TGraphAsymmErrors* graphXSecDirGammaTotArpp2760GeV_xT   = xTScalePhoton(graphXSecDirGammaTotArpp2760GeV, 2760.);
    TGraphAsymmErrors* graphInvYieldIncGammaStatpp2760GeV   = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldIncGammaStatErr");
    TGraphAsymmErrors* graphInvYieldIncGammaSyspp2760GeV    = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldIncGammaSysErr");
    TF1* fitInvYieldIncGammapp2760GeV                       = (TF1*) fileCombPP2760GeV->Get("Gamma2.76TeV/Fits/TwoComponentModelFitGamma");

    //--------------------------------------- pp 2.76TeV --------------------------------------------
    TFile* fileCombPP7TeV                                = new TFile( inputFileNamePP7TeV.Data());

    TH1D* histDRStatpp7TeV                               = (TH1D*)fileCombPP7TeV->Get("Gamma_7TeV_pp/DoubleRatioStatError");
    TGraphAsymmErrors* graphDRStatpp7TeV                 = new TGraphAsymmErrors(histDRStatpp7TeV);
    TGraphAsymmErrors* graphDRSyspp7TeV                  = (TGraphAsymmErrors*) fileCombPP7TeV->Get("Gamma_7TeV_pp/DoubleRatioSystError");
    while (graphDRStatpp7TeV->GetX()[0] < 0.51) graphDRStatpp7TeV->RemovePoint(0);
    while (graphDRSyspp7TeV->GetX()[0] < 0.51) graphDRSyspp7TeV->RemovePoint(0);
    while (graphDRStatpp7TeV->GetX()[graphDRStatpp7TeV->GetN()-1] > 10) graphDRStatpp7TeV->RemovePoint(graphDRStatpp7TeV->GetN()-1);
    while (graphDRSyspp7TeV->GetX()[graphDRSyspp7TeV->GetN()-1] > 10) graphDRSyspp7TeV->RemovePoint(graphDRSyspp7TeV->GetN()-1);

    TGraphAsymmErrors* graphDRTotpp7TeV                  = AddErrorsQuadraticallyTGraph(graphDRStatpp7TeV, graphDRSyspp7TeV);
    TGraphAsymmErrors* graphInvYieldDirGammaStatpp7TeV   = (TGraphAsymmErrors*) fileCombPP7TeV->Get("Gamma_7TeV_pp/graphDirGammaSpectrumStat");
    TGraphAsymmErrors* graphInvYieldDirGammaSyspp7TeV    = (TGraphAsymmErrors*) fileCombPP7TeV->Get("Gamma_7TeV_pp/graphDirGammaSpectrumSyst");
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpp7TeV  = (TGraphAsymmErrors*) fileCombPP7TeV->Get("Gamma_7TeV_pp/graphDirGammaSpectrumSummedAr");
    TGraphAsymmErrors* graphXSecDirGammaStatpp7TeV       = ScaleGraph(graphInvYieldDirGammaStatpp7TeV, xSection7TeV*recalcBarn);
    TGraphAsymmErrors* graphXSecDirGammaSyspp7TeV        = ScaleGraph(graphInvYieldDirGammaSyspp7TeV, xSection7TeV*recalcBarn);
    TGraphAsymmErrors* graphXSecDirGammaTotArpp7TeV      = ScaleGraph(graphInvYieldDirGammaTotArpp7TeV, xSection7TeV*recalcBarn);
    while (graphXSecDirGammaStatpp7TeV->GetX()[0] < 0.51) graphXSecDirGammaStatpp7TeV->RemovePoint(0);
    while (graphXSecDirGammaSyspp7TeV->GetX()[0] < 0.51) graphXSecDirGammaSyspp7TeV->RemovePoint(0);
    while (graphXSecDirGammaTotArpp7TeV->GetX()[0] < 0.51) graphXSecDirGammaTotArpp7TeV->RemovePoint(0);
    while (graphXSecDirGammaStatpp7TeV->GetX()[graphXSecDirGammaStatpp7TeV->GetN()-1] > 10) graphXSecDirGammaStatpp7TeV->RemovePoint(graphXSecDirGammaStatpp7TeV->GetN()-1);
    while (graphXSecDirGammaSyspp7TeV->GetX()[graphXSecDirGammaSyspp7TeV->GetN()-1] > 10) graphXSecDirGammaSyspp7TeV->RemovePoint(graphXSecDirGammaSyspp7TeV->GetN()-1);
    while (graphXSecDirGammaTotArpp7TeV->GetX()[graphXSecDirGammaTotArpp7TeV->GetN()-1] > 10) graphXSecDirGammaTotArpp7TeV->RemovePoint(graphXSecDirGammaTotArpp7TeV->GetN()-1);

    TGraphAsymmErrors* graphXSecDirGammaStatpp7TeV_xT    = xTScalePhoton(graphXSecDirGammaStatpp7TeV, 7000.);
    TGraphAsymmErrors* graphXSecDirGammaSyspp7TeV_xT     = xTScalePhoton(graphXSecDirGammaSyspp7TeV, 7000.);
    TGraphAsymmErrors* graphXSecDirGammaTotArpp7TeV_xT   = xTScalePhoton(graphXSecDirGammaTotArpp7TeV, 7000.);
    TGraphAsymmErrors* graphInvYieldIncGammaStatpp7TeV   = (TGraphAsymmErrors*) fileCombPP7TeV->Get("Gamma_7TeV_pp/graphInvYieldIncGammaStatErr");
    TGraphAsymmErrors* graphInvYieldIncGammaSyspp7TeV    = (TGraphAsymmErrors*) fileCombPP7TeV->Get("Gamma_7TeV_pp/graphInvYieldIncGammaSysErr");


    //--------------------------------------- pp 8TeV --------------------------------------------
    TFile* fileCombPP8TeV                                   = new TFile( inputFileNamePP8TeV.Data());
    TGraphAsymmErrors* graphDRStatpp8TeV                    = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphRGammaCombNonFitStatErr");
    TGraphAsymmErrors* graphDRSyspp8TeV                     = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphRGammaCombNonFitSysErr");
    TGraphAsymmErrors* graphDRTotpp8TeV                     = AddErrorsQuadraticallyTGraph(graphDRStatpp8TeV, graphDRSyspp8TeV);
    TGraphAsymmErrors* graphInvYieldDirGammaStatpp8TeV      = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldDirGammaNonFitStatErr");
    TGraphAsymmErrors* graphInvYieldDirGammaSyspp8TeV       = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldDirGammaNonFitSysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpp8TeV     = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldDirGammaNonFitSumErrAr");
    TGraphAsymmErrors* graphXSecDirGammaStatpp8TeV          = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvXSectionDirGammaNonFitStatErr");
    TGraphAsymmErrors* graphXSecDirGammaSyspp8TeV           = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvXSectionDirGammaNonFitSysErr");
    TGraphAsymmErrors* graphXSecDirGammaTotArpp8TeV         = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvXSectionDirGammaNonFitSumErrAr");
    TGraphAsymmErrors* graphXSecDirGammaStatpp8TeV_xT       = xTScalePhoton(graphXSecDirGammaStatpp8TeV, 8000.);
    TGraphAsymmErrors* graphXSecDirGammaSyspp8TeV_xT        = xTScalePhoton(graphXSecDirGammaSyspp8TeV, 8000.);
    TGraphAsymmErrors* graphXSecDirGammaTotArpp8TeV_xT      = xTScalePhoton(graphXSecDirGammaTotArpp8TeV, 8000.);
    TGraphAsymmErrors* graphInvYieldIncGammaStatpp8TeV      = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldIncGammaStatErr");
    TGraphAsymmErrors* graphInvYieldIncGammaSyspp8TeV       = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldIncGammaSysErr");
    TF1* fitInvYieldIncGammapp8TeV                          = (TF1*) fileCombPP8TeV->Get("Gamma8TeV/Fits/TwoComponentModelFitGamma");
    //--------------------------------------- pPb 5TeV --------------------------------------------
    TFile* fileCombPPb5TeV                                  = new TFile( inputFileNamePPb5TeV.Data());
    TGraphAsymmErrors* graphDRStatpPb5TeV[5]       = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphDRSyspPb5TeV[5]       = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaStatpPb5TeV[5]       = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaSyspPb5TeV[5]       = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpPb5TeV[5]       = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaStatpPb5TeVNCollScaled[5]       = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaSyspPb5TeVNCollScaled[5]        = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpPb5TeVNCollScaled[5]      = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotpPb5TeVNCollScaled[5]      = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpPb5TeVNCollScaledXT[5]    = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotpPb5TeVNCollScaledXT[5]      = {NULL, NULL, NULL, NULL, NULL};

    for (Int_t ncent = 0; ncent < 5; ncent++){
        graphDRStatpPb5TeV[ncent]                            = (TGraphAsymmErrors*) fileCombPPb5TeV->Get(Form("Gamma_pPb5TeV%s/graphRGammaCombStatErr",namepPbCent[ncent].Data()));
        graphDRSyspPb5TeV[ncent]                             = (TGraphAsymmErrors*) fileCombPPb5TeV->Get(Form("Gamma_pPb5TeV%s/graphRGammaCombSysErr",namepPbCent[ncent].Data()));
        graphInvYieldDirGammaStatpPb5TeV[ncent]              = (TGraphAsymmErrors*) fileCombPPb5TeV->Get(Form("Gamma_pPb5TeV%s/graphInvYieldDirGammaStatErr",namepPbCent[ncent].Data()));
        graphInvYieldDirGammaSyspPb5TeV[ncent]               = (TGraphAsymmErrors*) fileCombPPb5TeV->Get(Form("Gamma_pPb5TeV%s/graphInvYieldDirGammaSysErr",namepPbCent[ncent].Data()));
        graphInvYieldDirGammaTotArpPb5TeV[ncent]             = (TGraphAsymmErrors*) fileCombPPb5TeV->Get(Form("Gamma_pPb5TeV%s/graphInvYieldDirGammaSumErrAr",namepPbCent[ncent].Data()));
        if (graphInvYieldDirGammaStatpPb5TeV[ncent])
            graphInvYieldDirGammaStatpPb5TeVNCollScaled[ncent]    = ScaleGraph(graphInvYieldDirGammaStatpPb5TeV[ncent], 1./nCollpPb5[ncent]);
        if (graphInvYieldDirGammaSyspPb5TeV[ncent])
            graphInvYieldDirGammaSyspPb5TeVNCollScaled[ncent]     = ScaleGraph(graphInvYieldDirGammaSyspPb5TeV[ncent], 1./nCollpPb5[ncent]);
        if (graphInvYieldDirGammaTotArpPb5TeV[ncent])
            graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent]   = ScaleGraph(graphInvYieldDirGammaTotArpPb5TeV[ncent], 1./nCollpPb5[ncent]);

        if (graphInvYieldDirGammaStatpPb5TeVNCollScaled[ncent] && graphInvYieldDirGammaSyspPb5TeVNCollScaled[ncent] && graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent] ){
            for (Int_t pt = 0; pt < graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent]->GetN();pt++){
                Int_t currPt = 0;
                while (currPt < graphInvYieldDirGammaStatpPb5TeVNCollScaled[ncent]->GetN()){
                    if ( graphInvYieldDirGammaStatpPb5TeVNCollScaled[ncent]->GetX()[currPt] == graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent]->GetX()[pt] )
                        graphInvYieldDirGammaStatpPb5TeVNCollScaled[ncent]->RemovePoint(currPt);
                    else
                        currPt++;
                }
                currPt = 0;
                while (currPt < graphInvYieldDirGammaSyspPb5TeVNCollScaled[ncent]->GetN()){
                    if ( graphInvYieldDirGammaSyspPb5TeVNCollScaled[ncent]->GetX()[currPt] == graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent]->GetX()[pt] )
                        graphInvYieldDirGammaSyspPb5TeVNCollScaled[ncent]->RemovePoint(currPt);
                    else
                        currPt++;
                }
            }
            graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent] = AddErrorsQuadraticallyTGraph(graphInvYieldDirGammaStatpPb5TeVNCollScaled[ncent],graphInvYieldDirGammaSyspPb5TeVNCollScaled[ncent]);
            graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent]->Print();

            graphInvYieldDirGammaTotArpPb5TeVNCollScaledXT[ncent]     = xTScalePhoton(graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent], 5020.);
            graphInvYieldDirGammaTotpPb5TeVNCollScaledXT[ncent]       = xTScalePhoton(graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent], 5020.);

        } else if (graphInvYieldDirGammaStatpPb5TeVNCollScaled[ncent] && graphInvYieldDirGammaSyspPb5TeVNCollScaled[ncent]){
            graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent] = AddErrorsQuadraticallyTGraph(graphInvYieldDirGammaStatpPb5TeVNCollScaled[ncent],graphInvYieldDirGammaSyspPb5TeVNCollScaled[ncent]);
            graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent]->Print();
            graphInvYieldDirGammaTotpPb5TeVNCollScaledXT[ncent]       = xTScalePhoton(graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent], 5020.);
        }
    }
    //--------------------------------------- PbPb 2.76TeV --------------------------------------------
    TFile* fileCombPbPb2760GeV                              = new TFile( inputFileNamePbPb2760GeV.Data());
    TGraphAsymmErrors* graphDRStatPbPb[3]                           = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphDRSysPbPb[3]                            = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaStatPbPb[3]             = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaSysPbPb[3]              = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPb[3]            = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaStatPbPbNCollScaled[3]  = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaSysPbPbNCollScaled[3]   = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPbNCollScaled[3] = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotPbPbNCollScaled[3]   = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPbNCollScaledXT[3] = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotPbPbNCollScaledXT[3]   = {NULL, NULL, NULL};
    for (Int_t centPb   = 0; centPb < 3; centPb++){
        graphDRStatPbPb[centPb]                             = (TGraphAsymmErrors*) fileCombPbPb2760GeV->Get(Form("Gamma_PbPb_2.76TeV_%s/DR_comb_StatErr", namePbPbCent[centPb].Data()));
        graphDRSysPbPb[centPb]                              = (TGraphAsymmErrors*) fileCombPbPb2760GeV->Get(Form("Gamma_PbPb_2.76TeV_%s/DR_comb_SysErr", namePbPbCent[centPb].Data()));
        graphInvYieldDirGammaStatPbPb[centPb]               = (TGraphAsymmErrors*) fileCombPbPb2760GeV->Get(Form("Gamma_PbPb_2.76TeV_%s/DirGammaSpec_comb_StatErr", namePbPbCent[centPb].Data()));
        graphInvYieldDirGammaSysPbPb[centPb]                = (TGraphAsymmErrors*) fileCombPbPb2760GeV->Get(Form("Gamma_PbPb_2.76TeV_%s/DirGammaSpec_comb_SysErr", namePbPbCent[centPb].Data()));
        graphInvYieldDirGammaTotArPbPb[centPb]              = (TGraphAsymmErrors*) fileCombPbPb2760GeV->Get(Form("Gamma_PbPb_2.76TeV_%s/DirGammaSpec_comb_upperLimits", namePbPbCent[centPb].Data()));
        if (graphInvYieldDirGammaTotArPbPb[centPb]){
            if (!graphInvYieldDirGammaTotArPbPb[centPb]->GetN() > 0)
                graphInvYieldDirGammaTotArPbPb[centPb]      = NULL;
        }
        if (graphInvYieldDirGammaStatPbPb[centPb])
            graphInvYieldDirGammaStatPbPbNCollScaled[centPb]    = ScaleGraph(graphInvYieldDirGammaStatPbPb[centPb], 1./nColl[centPb]);
        if (graphInvYieldDirGammaSysPbPb[centPb])
            graphInvYieldDirGammaSysPbPbNCollScaled[centPb]     = ScaleGraph(graphInvYieldDirGammaSysPbPb[centPb], 1./nColl[centPb]);
        if (graphInvYieldDirGammaTotArPbPb[centPb])
            graphInvYieldDirGammaTotArPbPbNCollScaled[centPb]   = ScaleGraph(graphInvYieldDirGammaTotArPbPb[centPb], 1./nColl[centPb]);
        if (graphInvYieldDirGammaStatPbPbNCollScaled[centPb] && graphInvYieldDirGammaSysPbPbNCollScaled[centPb] && graphInvYieldDirGammaTotArPbPbNCollScaled[centPb] ){
            for (Int_t pt = 0; pt < graphInvYieldDirGammaTotArPbPbNCollScaled[centPb]->GetN();pt++){
                Int_t currPt = 0;
                while (currPt < graphInvYieldDirGammaStatPbPbNCollScaled[centPb]->GetN()){
                    if ( graphInvYieldDirGammaStatPbPbNCollScaled[centPb]->GetX()[currPt] == graphInvYieldDirGammaTotArPbPbNCollScaled[centPb]->GetX()[pt] )
                        graphInvYieldDirGammaStatPbPbNCollScaled[centPb]->RemovePoint(currPt);
                    else
                        currPt++;
                }
                currPt = 0;
                while (currPt < graphInvYieldDirGammaSysPbPbNCollScaled[centPb]->GetN()){
                    if ( graphInvYieldDirGammaSysPbPbNCollScaled[centPb]->GetX()[currPt] == graphInvYieldDirGammaTotArPbPbNCollScaled[centPb]->GetX()[pt] )
                        graphInvYieldDirGammaSysPbPbNCollScaled[centPb]->RemovePoint(currPt);
                    else
                        currPt++;
                }
            }
            graphInvYieldDirGammaTotPbPbNCollScaled[centPb] = AddErrorsQuadraticallyTGraph(graphInvYieldDirGammaStatPbPbNCollScaled[centPb],graphInvYieldDirGammaSysPbPbNCollScaled[centPb]);
            graphInvYieldDirGammaTotPbPbNCollScaled[centPb]->Print();

            graphInvYieldDirGammaTotArPbPbNCollScaledXT[centPb]     = xTScalePhoton(graphInvYieldDirGammaTotArPbPbNCollScaled[centPb], 2760.);
            graphInvYieldDirGammaTotPbPbNCollScaledXT[centPb]       = xTScalePhoton(graphInvYieldDirGammaTotPbPbNCollScaled[centPb], 2760.);

        } else if (graphInvYieldDirGammaStatPbPbNCollScaled[centPb] && graphInvYieldDirGammaSysPbPbNCollScaled[centPb]){
            graphInvYieldDirGammaTotPbPbNCollScaled[centPb] = AddErrorsQuadraticallyTGraph(graphInvYieldDirGammaStatPbPbNCollScaled[centPb],graphInvYieldDirGammaSysPbPbNCollScaled[centPb]);
            graphInvYieldDirGammaTotPbPbNCollScaled[centPb]->Print();
            graphInvYieldDirGammaTotPbPbNCollScaledXT[centPb]       = xTScalePhoton(graphInvYieldDirGammaTotPbPbNCollScaled[centPb], 2760.);
        }
    }

    TFile* fileTheoryPP                                     = new TFile( fileNameTheoryPP.Data());
    TGraph* graphTheoryNLODRpp2760GeVPaquettCenter          = (TGraph*) fileTheoryPP->Get("DirectPhoton/graphDRNLOPaquett_2760GeV_ALICECocktail");

    //*************************************************************************************************************************************************
    //*************************************** Read pp graphs for scaling plots from other experiments *************************************************
    //*************************************************************************************************************************************************
    TFile* fileOtherExperimentsPP                           = new TFile(fileNameExperimentPP.Data());
    TGraphAsymmErrors* xTOtherExperiments[nOtherEnergies];
    TGraphAsymmErrors* xSectionOtherExperiments[nOtherEnergies];
    for (Int_t i = 0; i < nOtherEnergies; i++){
        cout << "trying to find: " << Form("%s_tot_Gamma_%s_xT",nameOtherExpRead[i].Data(), nameOtherEnergyRead[i].Data()) << endl;
        xTOtherExperiments[i]                               = (TGraphAsymmErrors*)fileOtherExperimentsPP->Get(Form("%s_tot_Gamma_%s_xT",nameOtherExpRead[i].Data(), nameOtherEnergyRead[i].Data()));
        xSectionOtherExperiments[i]                         = (TGraphAsymmErrors*)fileOtherExperimentsPP->Get(Form("%s_tot_Gamma_%s",nameOtherExpRead[i].Data(), nameOtherEnergyRead[i].Data()));
    }
    TGraphAsymmErrors* graphRGammaLMEE7TeVStat              = (TGraphAsymmErrors*)fileOtherExperimentsPP->Get("RGamma_stat_LMEE_7TeV");
    TGraphAsymmErrors* graphRGammaLMEE7TeVSys               = (TGraphAsymmErrors*)fileOtherExperimentsPP->Get("RGamma_sys_LMEE_7TeV");
    TGraphAsymmErrors* graphRGammaLMEE7TeVTot               = (TGraphAsymmErrors*)fileOtherExperimentsPP->Get("RGamma_tot_LMEE_7TeV");

    //*************************************************************************************************************************************************
    //*************************************** Read p(d)-A graphs for scaling plots from other experiments ************************************************
    //*************************************************************************************************************************************************
    TFile* fileOtherExperimentsPA                           = new TFile(fileNameExperimentPA.Data());
    TGraphErrors* graphYieldDirGammapAu200GeVStat[2]                = {NULL, NULL};
    TGraphErrors* graphYieldDirGammapAu200GeVSys[2]                 = {NULL, NULL};
    TGraphErrors* graphYieldDirGammapAu200GeVTot[2]                 = {NULL, NULL};
    TGraphErrors* graphYieldDirGammapAu200GeVStatNCollScaled[2]     = {NULL, NULL};
    TGraphErrors* graphYieldDirGammapAu200GeVSysNCollScaled[2]      = {NULL, NULL};
    TGraphErrors* graphYieldDirGammapAu200GeVTotNCollScaled[2]      = {NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammapAu200GeVTotArNCollScaled[2]       = {NULL, NULL};
    TGraphErrors* graphYieldDirGammapAu200GeVTotNCollScaledXT[2]            = {NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammapAu200GeVTotArNCollScaledXT[2]     = {NULL, NULL};
    for (Int_t ncent = 0; ncent < 2; ncent++){
        graphYieldDirGammapAu200GeVStat[ncent]      = (TGraphErrors*)fileOtherExperimentsPA->Get(Form("Gamma/graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Stat_%s", nameCentspAu200Out[ncent].Data()));
        graphYieldDirGammapAu200GeVSys[ncent]       = (TGraphErrors*)fileOtherExperimentsPA->Get(Form("Gamma/graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Sys_%s", nameCentspAu200Out[ncent].Data()));
        if (graphYieldDirGammapAu200GeVStat[ncent] && graphYieldDirGammapAu200GeVSys[ncent])
            graphYieldDirGammapAu200GeVTot[ncent] = AddErrorsQuadraticallyTGraph(graphYieldDirGammapAu200GeVStat[ncent],graphYieldDirGammapAu200GeVSys[ncent]);
        // graphYieldDirGammapAu200GeVTot[ncent]       = (TGraphErrors*)fileOtherExperimentsPA->Get(Form("Gamma/graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Tot_%s", nameCentspAu200Out[ncent].Data()));
        if (graphYieldDirGammapAu200GeVStat[ncent])
            graphYieldDirGammapAu200GeVStatNCollScaled[ncent]   = ScaleGraph(graphYieldDirGammapAu200GeVStat[ncent], 1./(4520.0144*nCollpAu200[ncent]));
        if (graphYieldDirGammapAu200GeVSys[ncent])
            graphYieldDirGammapAu200GeVSysNCollScaled[ncent]    = ScaleGraph(graphYieldDirGammapAu200GeVSys[ncent], 1./(4520.0144*nCollpAu200[ncent]));
        if (graphYieldDirGammapAu200GeVTot[ncent])
            graphYieldDirGammapAu200GeVTotNCollScaled[ncent]    = ScaleGraph(graphYieldDirGammapAu200GeVTot[ncent], 1./(4520.0144*nCollpAu200[ncent]));

        if (graphYieldDirGammapAu200GeVTotNCollScaled[ncent]){
            Int_t nUpperLimits = 0;
            for (Int_t i=0; i< graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetN(); i++){
                if (graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetEY()[i]/graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetY()[i] > 0.8)
                    nUpperLimits++;
            }
            if (nUpperLimits>0){
                graphYieldDirGammapAu200GeVTotArNCollScaled[ncent] = new TGraphAsymmErrors(nUpperLimits);
                Int_t currPtBinNew = 0;
                for (Int_t i=0; i< graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetN(); i++){
                    if (graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetEY()[i]/graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetY()[i] > 0.8){
                        graphYieldDirGammapAu200GeVTotArNCollScaled[ncent]->SetPoint(currPtBinNew, graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetX()[i], graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetY()[i]+graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetEY()[i]);
                        graphYieldDirGammapAu200GeVTotArNCollScaled[ncent]->SetPointError(currPtBinNew, graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetEX()[i], graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetEX()[i], graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->GetY()[i]*1.6, 0);
                        graphYieldDirGammapAu200GeVTotNCollScaled[ncent]->RemovePoint(i);
                        i--;
                        currPtBinNew++;
                    }
                }
                graphYieldDirGammapAu200GeVTotArNCollScaledXT[ncent] = xTScalePhoton(graphYieldDirGammapAu200GeVTotArNCollScaled[ncent], 200.);
            }
            graphYieldDirGammapAu200GeVTotNCollScaledXT[ncent]  = xTScalePhoton(graphYieldDirGammapAu200GeVTotNCollScaled[ncent], 200.);
        }
    }
    TGraphErrors* graphYieldDirGammadAu200GeVStat                =NULL;
    TGraphErrors* graphYieldDirGammadAu200GeVSys                 =NULL;
    TGraphErrors* graphYieldDirGammadAu200GeVTot                 =NULL;
    TGraphErrors* graphYieldDirGammadAu200GeVStatNCollScaled     =NULL;
    TGraphErrors* graphYieldDirGammadAu200GeVSysNCollScaled      =NULL;
    TGraphErrors* graphYieldDirGammadAu200GeVTotNCollScaled      =NULL;
    TGraphAsymmErrors* graphYieldDirGammadAu200GeVTotArNCollScaled       =NULL;
    TGraphErrors* graphYieldDirGammadAu200GeVTotNCollScaledXT            =NULL;
    TGraphAsymmErrors* graphYieldDirGammadAu200GeVTotArNCollScaledXT     =NULL;
    for (Int_t ncent = 0; ncent < 2; ncent++){
        graphYieldDirGammadAu200GeVStat      = (TGraphErrors*)fileOtherExperimentsPA->Get(Form("Gamma/graph_InvYieldDirGamma_PHENIX_dAu_200GeV_Stat_%s", nameCentsdAu200Out.Data()));
        graphYieldDirGammadAu200GeVSys       = (TGraphErrors*)fileOtherExperimentsPA->Get(Form("Gamma/graph_InvYieldDirGamma_PHENIX_dAu_200GeV_Sys_%s", nameCentsdAu200Out.Data()));
        graphYieldDirGammadAu200GeVTot       = (TGraphErrors*)fileOtherExperimentsPA->Get(Form("Gamma/graph_InvYieldDirGamma_PHENIX_dAu_200GeV_Tot_%s", nameCentsdAu200Out.Data()));
        if (graphYieldDirGammadAu200GeVStat)
            graphYieldDirGammadAu200GeVStatNCollScaled   = ScaleGraph(graphYieldDirGammadAu200GeVStat, 1./(5624.336*nColldAu200));
        if (graphYieldDirGammadAu200GeVSys)
            graphYieldDirGammadAu200GeVSysNCollScaled    = ScaleGraph(graphYieldDirGammadAu200GeVSys, 1./(5624.336*nColldAu200));
        if (graphYieldDirGammadAu200GeVTot)
            graphYieldDirGammadAu200GeVTotNCollScaled    = ScaleGraph(graphYieldDirGammadAu200GeVTot, 1./(5624.336*nColldAu200));

        if (graphYieldDirGammadAu200GeVTotNCollScaled){
            Int_t nUpperLimits = 0;
            for (Int_t i=0; i< graphYieldDirGammadAu200GeVTotNCollScaled->GetN(); i++){
                if (graphYieldDirGammadAu200GeVTotNCollScaled->GetEY()[i]/graphYieldDirGammadAu200GeVTotNCollScaled->GetY()[i] > 0.8)
                    nUpperLimits++;
            }
            if (nUpperLimits>0){
                graphYieldDirGammadAu200GeVTotArNCollScaled = new TGraphAsymmErrors(nUpperLimits);
                Int_t currPtBinNew = 0;
                for (Int_t i=0; i< graphYieldDirGammadAu200GeVTotNCollScaled->GetN(); i++){
                    if (graphYieldDirGammadAu200GeVTotNCollScaled->GetEY()[i]/graphYieldDirGammadAu200GeVTotNCollScaled->GetY()[i] > 0.8){
                        graphYieldDirGammadAu200GeVTotArNCollScaled->SetPoint(currPtBinNew, graphYieldDirGammadAu200GeVTotNCollScaled->GetX()[i], graphYieldDirGammadAu200GeVTotNCollScaled->GetY()[i]+graphYieldDirGammadAu200GeVTotNCollScaled->GetEY()[i]);
                        graphYieldDirGammadAu200GeVTotArNCollScaled->SetPointError(currPtBinNew, graphYieldDirGammadAu200GeVTotNCollScaled->GetEX()[i], graphYieldDirGammadAu200GeVTotNCollScaled->GetEX()[i], graphYieldDirGammadAu200GeVTotNCollScaled->GetY()[i]*1.6, 0);
                        graphYieldDirGammadAu200GeVTotNCollScaled->RemovePoint(i);
                        i--;
                        currPtBinNew++;
                    }
                }
                graphYieldDirGammadAu200GeVTotArNCollScaledXT = xTScalePhoton(graphYieldDirGammadAu200GeVTotArNCollScaled, 200.);
            }
            graphYieldDirGammadAu200GeVTotNCollScaledXT  = xTScalePhoton(graphYieldDirGammadAu200GeVTotNCollScaled, 200.);
        }
    }

    Int_t npoints_ATLAS_dirgamma_pPb8TeV                  = 17;
    Double_t pT_ATLAS_dirgamma_pPb8TeV[17]                = { 30, 40, 50, 60, 70, 80, 95, 115, 137.5, 162.5, 187.5, 225, 275, 325, 375, 435, 510};
    // Double_t pTErr_ATLAS_dirgamma_pPb8TeV[17]                = { 5, 5, 5, 5, 5, 5, 10, 10, 12.5, 12.5, 12.5, 25, 25, 25, 25, 35, 40};
    Double_t pTErr_ATLAS_dirgamma_pPb8TeV[17]                = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Double_t xsection_ATLAS_dirgamma_pPb8TeV[17]      = { 157610, 44823.2, 17190.3, 7723.29, 3746.62, 2025.61, 957.837, 392.197, 158.167, 73.0265, 39.4647, 13.7071, 4.51454, 1.92182, 0.695431, 0.102295, 0.0921991};
    Double_t xsectionErr_ATLAS_dirgamma_pPb8TeV[17]       = { 1836.26, 457.206, 94.5415, 64.3389, 45.6204, 33.7393, 16.9768, 11.108, 6.29703, 4.32938, 3.22617, 1.35363, 0.797988, 0.564711, 0.528743, 0.358239, 0.353141};

    TGraphErrors* graphYieldDirGammapPb8TeV8000GeVATLASTot         = new TGraphErrors(npoints_ATLAS_dirgamma_pPb8TeV, pT_ATLAS_dirgamma_pPb8TeV, xsection_ATLAS_dirgamma_pPb8TeV,
                                                                            pTErr_ATLAS_dirgamma_pPb8TeV, xsectionErr_ATLAS_dirgamma_pPb8TeV);
    TGraphErrors* graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled      =NULL;
    TGraphErrors* graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaledXT            =NULL;
    TGraphAsymmErrors* graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled            =NULL;
    TGraphAsymmErrors* graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaledXT            =NULL;
        if (graphYieldDirGammapPb8TeV8000GeVATLASTot)
            graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled    = ScaleGraph(graphYieldDirGammapPb8TeV8000GeVATLASTot, 1./(70*1e12*7.09*2*TMath::Pi())); //70mb=179.7733 GeV/c2, Ncoll pPb8Tev = 7.09

        if (graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled){
            Int_t nUpperLimits = 0;
            for (Int_t i=0; i< graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetN(); i++){
                if (graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetEY()[i]/graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetY()[i] > 0.8)
                    nUpperLimits++;
            }
            if (nUpperLimits>0){
                graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled = new TGraphAsymmErrors(nUpperLimits);
                Int_t currPtBinNew = 0;
                for (Int_t i=0; i< graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetN(); i++){
                    if (graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetEY()[i]/graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetY()[i] > 0.8){
                        graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled->SetPoint(currPtBinNew, graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetX()[i], graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetY()[i]+graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetEY()[i]);
                        graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled->SetPointError(currPtBinNew, graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetEX()[i], graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetEX()[i], graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->GetY()[i]*3.6, 0);
                        graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->RemovePoint(i);
                        i--;
                        currPtBinNew++;
                    }
                }
                graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaledXT = xTScalePhoton(graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled, 8000.);
            }
            graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaledXT  = xTScalePhoton(graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled, 8000.);
        }

    //*************************************************************************************************************************************************
    //*************************************** Read A-A graphs for scaling plots from other experiments ************************************************
    //*************************************************************************************************************************************************
    TFile* fileOtherExperimentsAA                           = new TFile(fileNameExperimentAA.Data());
    TGraphErrors* graphYieldDirGammaAuAu200GeVStat[5][3]                = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
    TGraphErrors* graphYieldDirGammaAuAu200GeVSys[5][3]                 = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
    TGraphErrors* graphYieldDirGammaAuAu200GeVTot[5][3]                 = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
    TGraphErrors* graphYieldDirGammaAuAu200GeVStatNCollScaled[5][3]     = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
    TGraphErrors* graphYieldDirGammaAuAu200GeVSysNCollScaled[5][3]      = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
    TGraphErrors* graphYieldDirGammaAuAu200GeVTotNCollScaled[5][3]      = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVTotArNCollScaled[5][3]       = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
    TGraphErrors* graphYieldDirGammaAuAu200GeVTotNCollScaledXT[5][3]            = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVTotArNCollScaledXT[5][3]     = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
    for (Int_t centAu = 0; centAu < 5; centAu++){
        for (Int_t meth = 0; meth < 2; meth++){
            cout << Form("Gamma/graph_PHENIX_%d_AuAu200GeV_Stat_%s", meth, nameCentsAuAu200Out[centAu].Data()) << "\t" << nCollAuAu200[centAu] << endl;
            graphYieldDirGammaAuAu200GeVStat[centAu][meth]      = (TGraphErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_PHENIX_%d_AuAu200GeV_Stat_%s", meth, nameCentsAuAu200Out[centAu].Data()));
            graphYieldDirGammaAuAu200GeVSys[centAu][meth]       = (TGraphErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_PHENIX_%d_AuAu200GeV_Sys_%s", meth, nameCentsAuAu200Out[centAu].Data()));
            graphYieldDirGammaAuAu200GeVTot[centAu][meth]       = (TGraphErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_PHENIX_%d_AuAu200GeV_Tot_%s", meth, nameCentsAuAu200Out[centAu].Data()));
            if (graphYieldDirGammaAuAu200GeVStat[centAu][meth])
                graphYieldDirGammaAuAu200GeVStatNCollScaled[centAu][meth]   = ScaleGraph(graphYieldDirGammaAuAu200GeVStat[centAu][meth], 1./nCollAuAu200[centAu]);
            if (graphYieldDirGammaAuAu200GeVSys[centAu][meth])
                graphYieldDirGammaAuAu200GeVSysNCollScaled[centAu][meth]    = ScaleGraph(graphYieldDirGammaAuAu200GeVSys[centAu][meth], 1./nCollAuAu200[centAu]);
            if (graphYieldDirGammaAuAu200GeVTot[centAu][meth])
                graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]    = ScaleGraph(graphYieldDirGammaAuAu200GeVTot[centAu][meth], 1./nCollAuAu200[centAu]);

            if (graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]){
                Int_t nUpperLimits = 0;
                for (Int_t i=0; i< graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetN(); i++){
                    if (graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetEY()[i]/graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetY()[i] > 0.8)
                        nUpperLimits++;
                }
                if (nUpperLimits>0){
                    graphYieldDirGammaAuAu200GeVTotArNCollScaled[centAu][meth] = new TGraphAsymmErrors(nUpperLimits);
                    Int_t currPtBinNew = 0;
                    for (Int_t i=0; i< graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetN(); i++){
                        if (graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetEY()[i]/graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetY()[i] > 0.8){
                            graphYieldDirGammaAuAu200GeVTotArNCollScaled[centAu][meth]->SetPoint(currPtBinNew, graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetX()[i], graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetY()[i]+graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetEY()[i]);
                            graphYieldDirGammaAuAu200GeVTotArNCollScaled[centAu][meth]->SetPointError(currPtBinNew, graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetEX()[i], graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetEX()[i], graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->GetY()[i]*1.6, 0);
                            graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->RemovePoint(i);
                            i--;
                            currPtBinNew++;
                        }
                    }
                    graphYieldDirGammaAuAu200GeVTotArNCollScaledXT[centAu][meth] = xTScalePhoton(graphYieldDirGammaAuAu200GeVTotArNCollScaled[centAu][meth], 200.);
                }
                graphYieldDirGammaAuAu200GeVTotNCollScaledXT[centAu][meth]  = xTScalePhoton(graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth], 200.);
            }
        }
    }

    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVSTARStat[5]                = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVSTARSys[5]                 = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVSTARTot[5]                 = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVSTARStatNCollScaled[5]     = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVSTARSysNCollScaled[5]      = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[5]      = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[5]       = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVSTARTotNCollScaledXT[5]       = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu200GeVSTARTotArNCollScaledXT[5]     = {NULL, NULL, NULL, NULL, NULL};
    for (Int_t centAu = 0; centAu < 5; centAu++){
        cout << Form("Gamma/graph_STAR_AuAu200GeV_Stat_%s", nameCentsAuAu200STAROut[centAu].Data()) << "\t" << nCollAuAu200STAR[centAu] << endl;
        graphYieldDirGammaAuAu200GeVSTARStat[centAu]      = (TGraphAsymmErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_STAR_AuAu200GeV_Stat_%s", nameCentsAuAu200STAROut[centAu].Data()));
        graphYieldDirGammaAuAu200GeVSTARSys[centAu]       = (TGraphAsymmErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_STAR_AuAu200GeV_Sys_%s", nameCentsAuAu200STAROut[centAu].Data()));
        graphYieldDirGammaAuAu200GeVSTARTot[centAu]       = (TGraphAsymmErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_STAR_AuAu200GeV_Tot_%s", nameCentsAuAu200STAROut[centAu].Data()));
        if (graphYieldDirGammaAuAu200GeVSTARStat[centAu]){
            graphYieldDirGammaAuAu200GeVSTARStatNCollScaled[centAu]   = ScaleGraph(graphYieldDirGammaAuAu200GeVSTARStat[centAu], 1./nCollAuAu200STAR[centAu]);
        }
        if (graphYieldDirGammaAuAu200GeVSTARSys[centAu])
            graphYieldDirGammaAuAu200GeVSTARSysNCollScaled[centAu]    = ScaleGraph(graphYieldDirGammaAuAu200GeVSTARSys[centAu], 1./nCollAuAu200STAR[centAu]);
        if (graphYieldDirGammaAuAu200GeVSTARTot[centAu])
            graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]    = ScaleGraph(graphYieldDirGammaAuAu200GeVSTARTot[centAu], 1./nCollAuAu200STAR[centAu]);

        if (graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]){
            Int_t nUpperLimits = 0;
            for (Int_t i=0; i< graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetN(); i++){
                if (graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetEYlow()[i]/graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetY()[i] > 0.8)
                    nUpperLimits++;
            }
            if (nUpperLimits>0){
                graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu] = new TGraphAsymmErrors(nUpperLimits);
                Int_t currPtBinNew = 0;
                for (Int_t i=0; i< graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetN(); i++){
                    if (graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetEYlow()[i]/graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetY()[i] > 0.8){
                        graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu]->SetPoint(currPtBinNew, graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetX()[i], graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetY()[i]+graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetEYhigh()[i]);
                        graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu]->SetPointError(currPtBinNew, graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetEXlow()[i], graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetEXhigh()[i], graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->GetY()[i]*1.6, 0);
                        graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->RemovePoint(i);
                        i--;
                        currPtBinNew++;
                    }
                }
                graphYieldDirGammaAuAu200GeVSTARTotArNCollScaledXT[centAu] = xTScalePhoton(graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu], 200.);
            }
            graphYieldDirGammaAuAu200GeVSTARTotNCollScaledXT[centAu]  = xTScalePhoton(graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu], 200.);
        }
    }


    TGraphErrors* graphYieldDirGammaAuAu39GeVStat   = (TGraphErrors*)fileOtherExperimentsAA->Get("Gamma/graph_PHENIX_1_AuAu39GeV_Stat_0086");
    TGraphErrors* graphYieldDirGammaAuAu39GeVSys    = (TGraphErrors*)fileOtherExperimentsAA->Get("Gamma/graph_PHENIX_1_AuAu39GeV_Sys_0086");
    TGraphErrors* graphYieldDirGammaAuAu39GeVTot    = (TGraphErrors*)fileOtherExperimentsAA->Get("Gamma/graph_PHENIX_1_AuAu39GeV_Tot_0086");
    TGraphErrors* graphYieldDirGammaAuAu39GeVStatNCollScaled            = NULL;
    TGraphErrors* graphYieldDirGammaAuAu39GeVSysNCollScaled             = NULL;
    TGraphErrors* graphYieldDirGammaAuAu39GeVTotNCollScaled             = NULL;
    TGraphAsymmErrors* graphYieldDirGammaAuAu39GeVTotArNCollScaled      = NULL;
    TGraphErrors* graphYieldDirGammaAuAu39GeVTotNCollScaledXT           = NULL;
    TGraphAsymmErrors* graphYieldDirGammaAuAu39GeVTotArNCollScaledXT    = NULL;
    if (graphYieldDirGammaAuAu39GeVStat)
        graphYieldDirGammaAuAu39GeVStatNCollScaled   = ScaleGraph(graphYieldDirGammaAuAu39GeVStat, 1./nCollAuAu39);
    if (graphYieldDirGammaAuAu39GeVSys)
        graphYieldDirGammaAuAu39GeVSysNCollScaled    = ScaleGraph(graphYieldDirGammaAuAu39GeVSys, 1./nCollAuAu39);
    if (graphYieldDirGammaAuAu39GeVTot)
        graphYieldDirGammaAuAu39GeVTotNCollScaled    = ScaleGraph(graphYieldDirGammaAuAu39GeVTot, 1./nCollAuAu39);

    if (graphYieldDirGammaAuAu39GeVTotNCollScaled){
        Int_t nUpperLimits = 0;
        for (Int_t i=0; i< graphYieldDirGammaAuAu39GeVTotNCollScaled->GetN(); i++){
            if (graphYieldDirGammaAuAu39GeVTotNCollScaled->GetEY()[i]/graphYieldDirGammaAuAu39GeVTotNCollScaled->GetY()[i] > 0.8)
                nUpperLimits++;
        }
        if (nUpperLimits>0){
            graphYieldDirGammaAuAu39GeVTotArNCollScaled = new TGraphAsymmErrors(nUpperLimits);
            Int_t currPtBinNew = 0;
            for (Int_t i=0; i< graphYieldDirGammaAuAu39GeVTotNCollScaled->GetN(); i++){
                if (graphYieldDirGammaAuAu39GeVTotNCollScaled->GetEY()[i]/graphYieldDirGammaAuAu39GeVTotNCollScaled->GetY()[i] > 0.8){
                    graphYieldDirGammaAuAu39GeVTotArNCollScaled->SetPoint(currPtBinNew, graphYieldDirGammaAuAu39GeVTotNCollScaled->GetX()[i], graphYieldDirGammaAuAu39GeVTotNCollScaled->GetY()[i]+graphYieldDirGammaAuAu39GeVTotNCollScaled->GetEY()[i]);
                    graphYieldDirGammaAuAu39GeVTotArNCollScaled->SetPointError(currPtBinNew, graphYieldDirGammaAuAu39GeVTotNCollScaled->GetEX()[i], graphYieldDirGammaAuAu39GeVTotNCollScaled->GetEX()[i], graphYieldDirGammaAuAu39GeVTotNCollScaled->GetY()[i]*1.6, 0);
                    graphYieldDirGammaAuAu39GeVTotNCollScaled->RemovePoint(i);
                    i--;
                    currPtBinNew++;
                }
            }
            graphYieldDirGammaAuAu39GeVTotArNCollScaledXT = xTScalePhoton(graphYieldDirGammaAuAu39GeVTotArNCollScaled, 39.);
        }
        graphYieldDirGammaAuAu39GeVTotNCollScaledXT = xTScalePhoton(graphYieldDirGammaAuAu39GeVTotNCollScaled, 39.);
    } else {
      cout << "something went wrong" << endl;
      return;
    }

    TGraphErrors* graphYieldDirGammaAuAu62GeVStat[3]                   = {NULL, NULL, NULL};
    TGraphErrors* graphYieldDirGammaAuAu62GeVSys[3]                    = {NULL, NULL, NULL};
    TGraphErrors* graphYieldDirGammaAuAu62GeVTot[3]                    = {NULL, NULL, NULL};
    TGraphErrors* graphYieldDirGammaAuAu62GeVStatNCollScaled[3]        = {NULL, NULL, NULL};
    TGraphErrors* graphYieldDirGammaAuAu62GeVSysNCollScaled[3]         = {NULL, NULL, NULL};
    TGraphErrors* graphYieldDirGammaAuAu62GeVTotNCollScaled[3]         = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu62GeVTotArNCollScaled[3]  = {NULL, NULL, NULL};
    TGraphErrors* graphYieldDirGammaAuAu62GeVTotNCollScaledXT[3]         = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaAuAu62GeVTotArNCollScaledXT[3]  = {NULL, NULL, NULL};

    for (Int_t centAu = 0; centAu < 3; centAu++){
        cout << Form("Gamma/graph_PHENIX_2_AuAu62GeV_Stat_%s", nameCentsAuAu62Out[centAu].Data()) << "\t" << nCollAuAu62[centAu] << endl;
        graphYieldDirGammaAuAu62GeVStat[centAu]      = (TGraphErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_PHENIX_1_AuAu62GeV_Stat_%s", nameCentsAuAu62Out[centAu].Data()));
        graphYieldDirGammaAuAu62GeVSys[centAu]       = (TGraphErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_PHENIX_1_AuAu62GeV_Sys_%s", nameCentsAuAu62Out[centAu].Data()));
        graphYieldDirGammaAuAu62GeVTot[centAu]       = (TGraphErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_PHENIX_1_AuAu62GeV_Tot_%s", nameCentsAuAu62Out[centAu].Data()));
        if (graphYieldDirGammaAuAu62GeVStat[centAu])
            graphYieldDirGammaAuAu62GeVStatNCollScaled[centAu]   = ScaleGraph(graphYieldDirGammaAuAu62GeVStat[centAu], 1./nCollAuAu62[centAu]);
        if (graphYieldDirGammaAuAu62GeVSys[centAu])
            graphYieldDirGammaAuAu62GeVSysNCollScaled[centAu]    = ScaleGraph(graphYieldDirGammaAuAu62GeVSys[centAu], 1./nCollAuAu62[centAu]);
        if (graphYieldDirGammaAuAu62GeVTot[centAu])
            graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]    = ScaleGraph(graphYieldDirGammaAuAu62GeVTot[centAu], 1./nCollAuAu62[centAu]);

        if (graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]){
            Int_t nUpperLimits = 0;
            for (Int_t i=0; i< graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetN(); i++){
                if (graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetEY()[i]/graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetY()[i] > 0.8)
                    nUpperLimits++;
            }
            if (nUpperLimits>0){
                graphYieldDirGammaAuAu62GeVTotArNCollScaled[centAu] = new TGraphAsymmErrors(nUpperLimits);
                Int_t currPtBinNew = 0;
                for (Int_t i=0; i< graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetN(); i++){
                    if (graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetEY()[i]/graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetY()[i] > 0.8){
                        graphYieldDirGammaAuAu62GeVTotArNCollScaled[centAu]->SetPoint(currPtBinNew, graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetX()[i], graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetY()[i]+graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetEY()[i]);
                        graphYieldDirGammaAuAu62GeVTotArNCollScaled[centAu]->SetPointError(currPtBinNew, graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetEX()[i], graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetEX()[i], graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->GetY()[i]*1.6, 0);
                        graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->RemovePoint(i);
                        i--;
                        currPtBinNew++;
                    }
                }
                graphYieldDirGammaAuAu62GeVTotArNCollScaledXT[centAu]     = xTScalePhoton(graphYieldDirGammaAuAu62GeVTotArNCollScaled[centAu], 62.4);
            }
            graphYieldDirGammaAuAu62GeVTotNCollScaledXT[centAu]     = xTScalePhoton(graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu], 62.4);
        }
    }


    TGraphErrors* graphYieldDirGammaCuCu200GeVStat[2]                   = {NULL, NULL};
    TGraphErrors* graphYieldDirGammaCuCu200GeVSys[2]                    = {NULL, NULL};
    TGraphErrors* graphYieldDirGammaCuCu200GeVTot[2]                    = {NULL, NULL};
    TGraphErrors* graphYieldDirGammaCuCu200GeVStatNCollScaled[2]        = {NULL, NULL};
    TGraphErrors* graphYieldDirGammaCuCu200GeVSysNCollScaled[2]         = {NULL, NULL};
    TGraphErrors* graphYieldDirGammaCuCu200GeVTotNCollScaled[2]         = {NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCuCu200GeVTotArNCollScaled[2]  = {NULL, NULL};
    TGraphErrors* graphYieldDirGammaCuCu200GeVTotNCollScaledXT[2]         = {NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCuCu200GeVTotArNCollScaledXT[2]  = {NULL, NULL};

    for (Int_t centCu = 0; centCu < 2; centCu++){
        cout << Form("Gamma/graph_PHENIX_2_CuCu200GeV_Stat_%s", nameCentsCuCu200Out[centCu].Data()) << "\t" << nCollCuCu200[centCu] << endl;
        graphYieldDirGammaCuCu200GeVStat[centCu]      = (TGraphErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_PHENIX_2_CuCu200GeV_Stat_%s", nameCentsCuCu200Out[centCu].Data()));
        graphYieldDirGammaCuCu200GeVSys[centCu]       = (TGraphErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_PHENIX_2_CuCu200GeV_Sys_%s", nameCentsCuCu200Out[centCu].Data()));
        graphYieldDirGammaCuCu200GeVTot[centCu]       = (TGraphErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_PHENIX_2_CuCu200GeV_Tot_%s", nameCentsCuCu200Out[centCu].Data()));
        if (graphYieldDirGammaCuCu200GeVStat[centCu])
            graphYieldDirGammaCuCu200GeVStatNCollScaled[centCu]   = ScaleGraph(graphYieldDirGammaCuCu200GeVStat[centCu], 1./nCollCuCu200[centCu]);
        if (graphYieldDirGammaCuCu200GeVSys[centCu])
            graphYieldDirGammaCuCu200GeVSysNCollScaled[centCu]    = ScaleGraph(graphYieldDirGammaCuCu200GeVSys[centCu], 1./nCollCuCu200[centCu]);
        if (graphYieldDirGammaCuCu200GeVTot[centCu])
            graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]    = ScaleGraph(graphYieldDirGammaCuCu200GeVTot[centCu], 1./nCollCuCu200[centCu]);

        if (graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]){
            Int_t nUpperLimits = 0;
            for (Int_t i=0; i< graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetN(); i++){
                if (graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetEY()[i]/graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetY()[i] > 0.8)
                    nUpperLimits++;
            }
            if (nUpperLimits>0){
                graphYieldDirGammaCuCu200GeVTotArNCollScaled[centCu] = new TGraphAsymmErrors(nUpperLimits);
                Int_t currPtBinNew = 0;
                for (Int_t i=0; i< graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetN(); i++){
                    if (graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetEY()[i]/graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetY()[i] > 0.8){
                        graphYieldDirGammaCuCu200GeVTotArNCollScaled[centCu]->SetPoint(currPtBinNew, graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetX()[i], graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetY()[i]+graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetEY()[i]);
                        graphYieldDirGammaCuCu200GeVTotArNCollScaled[centCu]->SetPointError(currPtBinNew, graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetEX()[i], graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetEX()[i], graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->GetEY()[i]*1.6, 0);
                        graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->RemovePoint(i);
                        i--;
                        currPtBinNew++;
                    }
                }
                graphYieldDirGammaCuCu200GeVTotArNCollScaledXT[centCu]     = xTScalePhoton(graphYieldDirGammaCuCu200GeVTotArNCollScaled[centCu], 200.);
            }
            graphYieldDirGammaCuCu200GeVTotNCollScaledXT[centCu]     = xTScalePhoton(graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu], 200.);
        }
    }


    TGraphAsymmErrors* graphYieldDirGammaCMSPbPb2760GeVStat[3]                   = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCMSPbPb2760GeVSys[3]                    = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCMSPbPb2760GeVTot[3]                    = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCMSPbPb2760GeVStatNCollScaled[3]        = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCMSPbPb2760GeVSysNCollScaled[3]         = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[3]         = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[3]       = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCMSPbPb2760GeVTotNCollScaledXT[3]         = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaledXT[3]       = {NULL, NULL, NULL};
    for (Int_t centPb = 0; centPb < 3; centPb++){
        cout << Form("Gamma/graph_CMS_PbPb2760GeV_Stat_%s", namePbPbCentCMSOut[centPb].Data()) << "\t" << nCollCMS[centPb] << endl;
        graphYieldDirGammaCMSPbPb2760GeVStat[centPb]      = (TGraphAsymmErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_CMS_PbPb2760GeV_Stat_%s", namePbPbCentCMSOut[centPb].Data()));
        graphYieldDirGammaCMSPbPb2760GeVSys[centPb]       = (TGraphAsymmErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_CMS_PbPb2760GeV_Sys_%s", namePbPbCentCMSOut[centPb].Data()));
        graphYieldDirGammaCMSPbPb2760GeVTot[centPb]       = (TGraphAsymmErrors*)fileOtherExperimentsAA->Get(Form("Gamma/graph_CMS_PbPb2760GeV_Tot_%s", namePbPbCentCMSOut[centPb].Data()));
        if (graphYieldDirGammaCMSPbPb2760GeVStat[centPb])
            graphYieldDirGammaCMSPbPb2760GeVStatNCollScaled[centPb]   = ScaleGraph(graphYieldDirGammaCMSPbPb2760GeVStat[centPb], 1./nCollCMS[centPb]);
        if (graphYieldDirGammaCMSPbPb2760GeVSys[centPb])
            graphYieldDirGammaCMSPbPb2760GeVSysNCollScaled[centPb]    = ScaleGraph(graphYieldDirGammaCMSPbPb2760GeVSys[centPb], 1./nCollCMS[centPb]);
        if (graphYieldDirGammaCMSPbPb2760GeVTot[centPb])
            graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]    = ScaleGraph(graphYieldDirGammaCMSPbPb2760GeVTot[centPb], 1./nCollCMS[centPb]);

        if (graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]){
            Int_t nUpperLimits = 0;
            for (Int_t i=0; i< graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetN(); i++){
                if (graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetEYlow()[i]/graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetY()[i] > 0.8)
                    nUpperLimits++;
            }
            if (nUpperLimits>0){
                graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb] = new TGraphAsymmErrors(nUpperLimits);
                Int_t currPtBinNew = 0;
                for (Int_t i=0; i< graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetN(); i++){
                    if (graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetEYlow()[i]/graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetY()[i] > 0.8){
                        graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb]->SetPoint(currPtBinNew, graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetX()[i], graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetY()[i]+graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetEYhigh()[i]);
                        graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb]->SetPointError(currPtBinNew, graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetEXlow()[i], graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetEXhigh()[i], graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->GetY()[i]*1.6, 0);
                        graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->RemovePoint(i);
                        i--;
                        currPtBinNew++;
                    }
                }
                graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaledXT[centPb]     = xTScalePhoton(graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb], 2760.);
            }
            graphYieldDirGammaCMSPbPb2760GeVTotNCollScaledXT[centPb]  = xTScalePhoton(graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb], 2760.);
        }
    }

    TGraphAsymmErrors* graphYieldDirGammaWA98Pb17GeVTot                     = (TGraphAsymmErrors*)fileOtherExperimentsAA->Get("Gamma/graph_WA98_PbPb17.2GeV_Tot_0010");
    TGraphAsymmErrors* graphYieldDirGammaWA98Pb17GeVTotNCollScaled          = NULL;
    TGraphAsymmErrors* graphYieldDirGammaWA98Pb17GeVTotArNCollScaled        = NULL;
    TGraphAsymmErrors* graphYieldDirGammaWA98Pb17GeVTotNCollScaledXT        = NULL;
    TGraphAsymmErrors* graphYieldDirGammaWA98Pb17GeVTotArNCollScaledXT      = NULL;
    if (graphYieldDirGammaWA98Pb17GeVTot)
        graphYieldDirGammaWA98Pb17GeVTotNCollScaled    = ScaleGraph(graphYieldDirGammaWA98Pb17GeVTot, 1./nCollPbPb17);
    if (graphYieldDirGammaWA98Pb17GeVTotNCollScaled){
        Int_t nUpperLimits = 0;
        for (Int_t i=0; i< graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetN(); i++){
            if (graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetEYlow()[i]/graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetY()[i] > 0.8)
                nUpperLimits++;
        }
        if (nUpperLimits>0){
            graphYieldDirGammaWA98Pb17GeVTotArNCollScaled = new TGraphAsymmErrors(nUpperLimits);
            Int_t currPtBinNew = 0;
            for (Int_t i=0; i< graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetN(); i++){
                if (graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetEYlow()[i]/graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetY()[i] > 0.8){
                    graphYieldDirGammaWA98Pb17GeVTotArNCollScaled->SetPoint(currPtBinNew, graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetX()[i], graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetY()[i]+graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetEYhigh()[i]);
//                     if (graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetX()[i] == )
                    graphYieldDirGammaWA98Pb17GeVTotArNCollScaled->SetPointError(currPtBinNew, graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetEXlow()[i], graphYieldDirGammaWA98Pb17GeVTotNCollScaled->GetEXlow()[i], graphYieldDirGammaWA98Pb17GeVTotArNCollScaled->GetY()[currPtBinNew]*0.8, 0);
                    graphYieldDirGammaWA98Pb17GeVTotNCollScaled->RemovePoint(i);
                    i--;
                    currPtBinNew++;
                }
            }
            cout << "WA 98 upper:" << endl;
            graphYieldDirGammaWA98Pb17GeVTotArNCollScaled->Print();
            graphYieldDirGammaWA98Pb17GeVTotArNCollScaledXT = xTScalePhoton(graphYieldDirGammaWA98Pb17GeVTotArNCollScaled, 17.2);
//             graphYieldDirGammaWA98Pb17GeVTotArNCollScaledXT->Print();
        }
        cout << "WA 98 tot:" << endl;
        graphYieldDirGammaWA98Pb17GeVTotNCollScaled->Print();
        graphYieldDirGammaWA98Pb17GeVTotNCollScaledXT = xTScalePhoton(graphYieldDirGammaWA98Pb17GeVTotNCollScaled, 17.2);
    } else {
        cout << "something went wrong" << endl;
        return;
    }


    //*************************************************************************************************************************************************
    //*************************************** Prepare for plotting ************************************************************************************
    //*************************************************************************************************************************************************

    TGraphAsymmErrors* graphInvYieldDirGammaStatPbPbPlot[3]             = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaSysPbPbPlot[3]              = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPbPlot[3]            = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpp2760GeVScaled[3]     = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaSyspp2760GeVScaled[3]       = {NULL, NULL, NULL};
    TGraphAsymmErrors* graphInvYieldDirGammaStatpp2760GeVScaled[3]      = {NULL, NULL, NULL};

    for (Int_t centPb = 0; centPb < 3; centPb++){
        if (graphInvYieldDirGammaStatPbPb[centPb]) graphInvYieldDirGammaStatPbPbPlot[centPb]        = ScaleGraph(graphInvYieldDirGammaStatPbPb[centPb],scaleFacPlotAPbPb[centPb]);
        if (graphInvYieldDirGammaStatPbPbPlot[centPb]) ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPbPlot[centPb]);
        if (graphInvYieldDirGammaSysPbPb[centPb]) graphInvYieldDirGammaSysPbPbPlot[centPb]          = ScaleGraph(graphInvYieldDirGammaSysPbPb[centPb],scaleFacPlotAPbPb[centPb]);
        if (graphInvYieldDirGammaTotArPbPb[centPb]) graphInvYieldDirGammaTotArPbPbPlot[centPb]      = ScaleGraph(graphInvYieldDirGammaTotArPbPb[centPb],scaleFacPlotAPbPb[centPb]);
        graphInvYieldDirGammaStatpp2760GeVScaled[centPb]                                            = ScaleGraph(graphInvYieldDirGammaStatpp2760GeV,nColl[centPb]*scaleFacPlotAPbPb[centPb]);
        if (graphInvYieldDirGammaStatpp2760GeVScaled[centPb]) ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatpp2760GeVScaled[centPb]);
        graphInvYieldDirGammaSyspp2760GeVScaled[centPb]                                             = ScaleGraph(graphInvYieldDirGammaSyspp2760GeV,nColl[centPb]*scaleFacPlotAPbPb[centPb]);
        graphInvYieldDirGammaTotArpp2760GeVScaled[centPb]                                           = ScaleGraph(graphInvYieldDirGammaTotArpp2760GeV,nColl[centPb]*scaleFacPlotAPbPb[centPb]);
    }

    TCanvas *canvasDirGamma = new TCanvas("canvasDirGamma","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasDirGamma, 0.175, 0.01, 0.01, 0.073);
    canvasDirGamma->SetLogy();
    canvasDirGamma->SetLogx();

    Int_t textSizeLabelsPixelDirGam = 48;
    Double_t textsizeLabelsDirGamma = 0;
    if (canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) < canvasDirGamma->YtoPixel(canvasDirGamma->GetY1())){
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) ;
    } else {
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->YtoPixel(canvasDirGamma->GetY1());
    }

    TH1D* dummyDirGamma = new TH1D("dummyDirGamma", "dummyDirGamma", 1000, 0., 22.);
    SetStyleHistoTH1ForGraphs( dummyDirGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{inel}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.8);
    dummyDirGamma->GetYaxis()->SetRangeUser( 1.2e-9,1.5e1);
    dummyDirGamma->GetXaxis()->SetLabelOffset(-0.015);
    dummyDirGamma->GetXaxis()->SetTickLength(0.025);
    dummyDirGamma->GetYaxis()->SetTickLength(0.025);
    dummyDirGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDirGamma->DrawCopy();

        for (Int_t centPb = 0; centPb < 3; centPb++){
            if (graphInvYieldDirGammaSysPbPb[centPb]){
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSysPbPb[centPb], markerStyleComb[centPb], markerSizeComb[centPb], colorComb[centPb] , colorComb[centPb], widthLinesBoxes, kTRUE);
                graphInvYieldDirGammaSysPbPb[centPb]->Draw("E2same");
            }
            if (graphInvYieldDirGammaStatPbPb[centPb]){
                ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb[centPb]);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatPbPb[centPb], markerStyleComb[centPb], markerSizeComb[centPb], colorComb[centPb] , colorComb[centPb]);
                graphInvYieldDirGammaStatPbPb[centPb]->Draw("p,E1Z,same");
            }
            if (graphInvYieldDirGammaTotArPbPb[centPb]){
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArPbPb[centPb] , 1, 3, colorComb[centPb], colorComb[centPb], 1.8, kTRUE);
                graphInvYieldDirGammaTotArPbPb[centPb]->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPb[centPb]);
            }
        }

        if (graphInvYieldDirGammaSyspp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp2760GeV , markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV, colorCombpp2760GeV, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspp2760GeV->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatpp2760GeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatpp2760GeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV, colorCombpp2760GeV);
            graphInvYieldDirGammaStatpp2760GeV->Draw("p,E1Z,same");
        }
        if (graphInvYieldDirGammaTotArpp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp2760GeV , 1, 3, colorCombpp2760GeV, colorCombpp2760GeV, 1.8, kTRUE);
            graphInvYieldDirGammaTotArpp2760GeV->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeV);
        }
        if (graphInvYieldDirGammaSyspPb5TeV[4]){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspPb5TeV[4], markerStyleCombpPb[4], markerSizeCombpPb[4], colorCombpPb[4] , colorCombpPb[4], widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspPb5TeV[4]->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatpPb5TeV[4]){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatpPb5TeV[4]);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpPb5TeV[4], markerStyleCombpPb[4], markerSizeCombpPb[4], colorCombpPb[4] , colorCombpPb[4]);
            graphInvYieldDirGammaStatpPb5TeV[4]->Draw("p,E1Z,same");
        }
        if (graphInvYieldDirGammaTotArpPb5TeV[4]){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpPb5TeV[4] , 1, 3, colorCombpPb[4], colorCombpPb[4], 1.8, kTRUE);
            graphInvYieldDirGammaTotArpPb5TeV[4]->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpPb5TeV[4]);
        }

        TLatex *labelALICEDirGamma = new TLatex(0.95,0.94,textALICE);
        SetStyleTLatex( labelALICEDirGamma, 42, 4, 1, 43, kTRUE, 31);
        labelALICEDirGamma->Draw();

        TLegend* legendDirGammaPP = GetAndSetLegend2(0.21, 0.10+(6*textsizeLabelsDirGamma*0.85), 0.21+0.21, 0.10+(8*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1,
                                                     collisionSystempp2760GeV.Data(), 42, 0.25);
        legendDirGammaPP->AddEntry(graphInvYieldDirGammaSyspp2760GeV,"ALICE","pf");

        legendDirGammaPP->Draw();
        TLegend* legendDirGammaPPb = GetAndSetLegend2(0.21, 0.10+(4*textsizeLabelsDirGamma*0.85), 0.21+0.21, 0.10+(6*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1,
                                                      collisionSystempPb5TeV.Data(), 42, 0.25);
        legendDirGammaPPb->AddEntry(graphInvYieldDirGammaSyspPb5TeV[4],"ALICE","pf");
        legendDirGammaPPb->Draw();

        TLegend* legendDirGamma = GetAndSetLegend2(0.21, 0.10, 0.21+0.21, 0.10+(4*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, collisionSystemPbPb2760GeV.Data(), 42, 0.25);
        legendDirGamma->AddEntry(graphInvYieldDirGammaSysPbPb[0],"  0-20% ALICE","pf");
        legendDirGamma->AddEntry(graphInvYieldDirGammaSysPbPb[1],"20-40% ALICE","pf");
        legendDirGamma->AddEntry(graphInvYieldDirGammaSysPbPb[2],"40-80% ALICE","pf");
        legendDirGamma->Draw();


    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled.pdf",outputDir.Data()));
    dummyDirGamma->DrawCopy();

        for (Int_t centPb = 0; centPb < 3; centPb++){
            if (graphInvYieldDirGammaSysPbPb[centPb]){
                graphInvYieldDirGammaSysPbPb[centPb]->Draw("E2same");
            }
            if (graphInvYieldDirGammaStatPbPb[centPb]){
                graphInvYieldDirGammaStatPbPb[centPb]->Draw("p,E1Z,same");
            }
            if (graphInvYieldDirGammaTotArPbPb[centPb]){
                graphInvYieldDirGammaTotArPbPb[centPb]->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPb[centPb]);
            }
        }
        legendDirGamma->Draw();


    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled_OnlyPbPb.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled_OnlyPbPb.pdf",outputDir.Data()));

    dummyDirGamma->DrawCopy();

        for (Int_t centPb = 0; centPb < 3; centPb++){
            if (graphInvYieldDirGammaSysPbPb[centPb]){
                graphInvYieldDirGammaSysPbPb[centPb]->Draw("E2same");
            }
            if (graphInvYieldDirGammaStatPbPb[centPb]){
                graphInvYieldDirGammaStatPbPb[centPb]->Draw("p,E1Z,same");
            }
            if (graphInvYieldDirGammaTotArPbPb[centPb]){
                graphInvYieldDirGammaTotArPbPb[centPb]->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPb[centPb]);
            }
        }
        for (Int_t centAu = 0; centAu < 4; centAu++){
            for (Int_t meth = 0; meth < 2; meth++){
//                 cout << "trying to plot: " << centAu << "\t" << meth << endl;
                if (graphYieldDirGammaAuAu200GeVSys[centAu][meth]){
                    DrawGammaSetMarkerTGraphErr(graphYieldDirGammaAuAu200GeVSys[centAu][meth], markerStyleCombAuAu200[centAu], markerSizeCombAuAu200[centAu], colorCombAuAu200[centAu] , colorCombAuAu200[centAu], widthLinesBoxes, kTRUE);
                    graphYieldDirGammaAuAu200GeVSys[centAu][meth]->Draw("E2same");
//                     graphYieldDirGammaAuAu200GeVSys[centAu][meth]->Print();
                }
                if (graphYieldDirGammaAuAu200GeVStat[centAu][meth]){
//                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVStat[centAu][meth]);
                    DrawGammaSetMarkerTGraphErr(graphYieldDirGammaAuAu200GeVStat[centAu][meth], markerStyleCombAuAu200[centAu], markerSizeCombAuAu200[centAu], colorCombAuAu200[centAu] , colorCombAuAu200[centAu]);
                    graphYieldDirGammaAuAu200GeVStat[centAu][meth]->Draw("p,E1Z,same");
                }
            }
        }
        legendDirGamma->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled_AA.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled_AA.pdf",outputDir.Data()));



    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    TCanvas *canvasDirGammapPb = new TCanvas("canvasDirGammapPb","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasDirGammapPb, 0.18, 0.01, 0.01, 0.08);
    canvasDirGammapPb->SetLogy();
    canvasDirGammapPb->SetLogx();

    TH1D* dummyDirGammapPb1 = new TH1D("dummyDirGammapPb", "dummyDirGammapPb", 1000, 0.09, 820.);
    SetStyleHistoTH1ForGraphs( dummyDirGammapPb1, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{coll}} #frac{1}{2#pi #it{N}_{inel}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.78, 1.85);
    dummyDirGammapPb1->GetYaxis()->SetRangeUser( 1.2e-18,4e2);
    dummyDirGammapPb1->GetXaxis()->SetLabelOffset(-0.012);
    dummyDirGammapPb1->GetXaxis()->SetTickLength(0.025);
    dummyDirGammapPb1->GetYaxis()->SetTickLength(0.025);
//     dummyDirGammapPb1->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDirGammapPb1->DrawCopy();

        TLatex *labelXTscaledGammapPb = new TLatex(0.95,0.94,"p + A #rightarrow #gamma + X");
        SetStyleTLatex( labelXTscaledGammapPb, textSizeLabelsPixelDirGam, 4, 1, 43, kTRUE, 31);
        labelXTscaledGammapPb->Draw();
        TLatex *labelXTscaledGammapPb2 = new TLatex(0.95,0.90,"y #approx 0");
        SetStyleTLatex( labelXTscaledGammapPb2, textSizeLabelsPixelDirGam, 4, 1, 43, kTRUE, 31);
        labelXTscaledGammapPb2->Draw();

        Int_t nRowsGammapPb              = 2;
        Int_t nRowsGammapPb8000           = 2;
        Int_t nRowsGammadAu200           = 2;
        Int_t nRowsGammapAu200           = 3;
        Int_t nColumnsGammapPb           = 3;
        Double_t yMinLegpPb                = 0.1;
        Double_t xMinLegpPb                = 0.22;
        Double_t widthColumnpPb            = 0.23;
        TLegend* legendInvYieldGammapPb8000      = GetAndSetLegend2( xMinLegpPb, yMinLegpPb, xMinLegpPb+widthColumnpPb*2, yMinLegpPb+(nRowsGammapPb*textsizeLabelsDirGamma*0.5) ,0.5*textsizeLabelsDirGamma,
                                                                nColumnsGammapPb, Form("#scale[1.25]{%s}",collisionSystempPb8TeV.Data()), 42, 0.15);
        TLegend* legendInvYieldGammapPb   = GetAndSetLegend2( xMinLegpPb, yMinLegpPb+(nRowsGammapPb)*textsizeLabelsDirGamma*0.5+0.01, xMinLegpPb+widthColumnpPb*3,
                                                                yMinLegpPb+((nRowsGammapPb+nRowsGammapAu200)*textsizeLabelsDirGamma*0.5)+0.01 ,0.5*textsizeLabelsDirGamma,
                                                                nColumnsGammapPb, Form("#scale[1.25]{%s}",collisionSystempPb5TeV.Data()), 42, 0.15);
        TLegend* legendInvYieldGammapAu200   = GetAndSetLegend2( xMinLegpPb, yMinLegpPb+(nRowsGammapPb+nRowsGammapAu200)*textsizeLabelsDirGamma*0.5+0.02, xMinLegpPb+widthColumnpPb*2,
                                                                yMinLegpPb+((nRowsGammapPb+nRowsGammapAu200+nRowsGammadAu200)*textsizeLabelsDirGamma*0.5)+0.02 ,0.5*textsizeLabelsDirGamma,
                                                                nColumnsGammapPb, Form("#scale[1.25]{%s}",collisionSystempAu200GeV.Data()), 42, 0.15);
        TLegend* legendInvYieldGammadAu200    = GetAndSetLegend2( xMinLegpPb, yMinLegpPb+(nRowsGammapPb+nRowsGammapPb8000+nRowsGammapAu200)*textsizeLabelsDirGamma*0.5+0.03, xMinLegpPb+widthColumnpPb*nColumnsGammapPb,
                                                                yMinLegpPb+((nRowsGammapPb+nRowsGammapPb8000+nRowsGammadAu200+nRowsGammapAu200)*textsizeLabelsDirGamma*0.5)+0.03 ,0.5*textsizeLabelsDirGamma,
                                                                nColumnsGammapPb, Form("#scale[1.25]{%s}",collisionSystemdAu200GeV.Data()), 42, 0.15);

        for (Int_t ncent = 0; ncent < 5; ncent++){
            if (graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent]){
                ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent]);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent], markerStyleCombpPbMC[ncent], markerSizeCombpPb[ncent], colorCombpPb[ncent] , colorCombpPb[ncent]);
                graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent]->Draw("p,E1,same");
                legendInvYieldGammapPb->AddEntry(graphInvYieldDirGammaTotpPb5TeVNCollScaled[ncent], labelALICEpPb[ncent].Data(), "p");
            }
            if (graphInvYieldDirGammaTotArpPb5TeV[ncent]){
                cout << "Pb:    " << ncent << endl;
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent] , 1, 3, colorCombpPb[ncent], colorCombpPb[ncent], 1.8, kTRUE);
                graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent]->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent], 0.05);
                graphInvYieldDirGammaTotArpPb5TeVNCollScaled[ncent]->Print();
            }
        }

        // plotting PHENIX pAu 200
        for (Int_t centAu = 0; centAu < 2; centAu++){
            if (graphYieldDirGammapAu200GeVTotNCollScaled[centAu]){
                //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammapAu200GeVStat[centAu]);
                DrawGammaSetMarkerTGraphErr(graphYieldDirGammapAu200GeVTotNCollScaled[centAu], markerStyleCombMCpAu200[centAu], markerSizeCombpAu200[centAu], colorCombpAu200[centAu], colorCombpAu200[centAu] );
                graphYieldDirGammapAu200GeVTotNCollScaled[centAu]->Draw("p,E1,same");
                legendInvYieldGammapAu200->AddEntry(graphYieldDirGammapAu200GeVTotNCollScaled[centAu], labelPHENIXpAu200[centAu].Data(), "p");
            }
            if (graphYieldDirGammapAu200GeVTotArNCollScaled[centAu]){
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammapAu200GeVTotArNCollScaled[centAu] , 1, 3, colorCombpAu200[centAu], colorCombpAu200[centAu], 1.8, kTRUE);
                graphYieldDirGammapAu200GeVTotArNCollScaled[centAu]->Draw(">,same");
                graphYieldDirGammapAu200GeVTotArNCollScaled[centAu]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammapAu200GeVTotArNCollScaled[centAu], 0.1);
            }
        }

        // plotting PHENIX dAu 200
        if (graphYieldDirGammadAu200GeVTotNCollScaled){
            //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammadAu200GeVStat);
            DrawGammaSetMarkerTGraphErr(graphYieldDirGammadAu200GeVTotNCollScaled, markerStyleCombMCdAu200, markerSizeCombdAu200, colorCombdAu200, colorCombdAu200 );
            graphYieldDirGammadAu200GeVTotNCollScaled->Draw("p,E1,same");
            legendInvYieldGammadAu200->AddEntry(graphYieldDirGammadAu200GeVTotNCollScaled, labelPHENIXdAu200.Data(), "p");
        }
        if (graphYieldDirGammadAu200GeVTotArNCollScaled){
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammadAu200GeVTotArNCollScaled , 1, 3, colorCombdAu200, colorCombdAu200, 1.8, kTRUE);
            graphYieldDirGammadAu200GeVTotArNCollScaled->Draw(">,same");
            graphYieldDirGammadAu200GeVTotArNCollScaled->Print();
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammadAu200GeVTotArNCollScaled, 0.1);
        }

        // plotting ATLAS pPb 8000
        if (graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled){
            //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammadAu200GeVStat);
            DrawGammaSetMarkerTGraphErr(graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled, markerStyleCombpPb8ATLASMC, markerSizeCombpPb8ATLAS, colorCombpPb8ATLAS, colorCombpPb8ATLAS );
            graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled->Draw("p,E1,same");
            legendInvYieldGammapPb8000->AddEntry(graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaled, labelpPb8TeVATLAS.Data(), "p");
        }
        if (graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled){
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled , 1, 3, colorCombpPb8ATLAS, colorCombpPb8ATLAS, 1.8, kTRUE);
            graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled->Draw(">,same");
            graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled->Print();
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaled, 0.1);
        }

        legendInvYieldGammapPb->Draw();
        legendInvYieldGammapAu200->Draw();
        legendInvYieldGammadAu200->Draw();
        legendInvYieldGammapPb8000->Draw();

    canvasDirGammapPb->Print(Form("%s/DirGammaSpectra_SmallSystems_NCollScaled.%s",outputDir.Data(),suffix.Data()));
    canvasDirGammapPb->Print(Form("%s/DirGammaSpectra_SmallSystems_NCollScaled.pdf",outputDir.Data()));



    canvasDirGammapPb->cd();
    TH1D* dummyDirGammapPb2 = new TH1D("dummyDirGamma", "dummyDirGamma", 1000, 3e-5, 3);
    SetStyleHistoTH1ForGraphs( dummyDirGammapPb2, "#it{x}_{T} = 2#it{p}_{T}/#sqrt{#it{s}_{_{NN}}}", "#sqrt{(#it{s}_{_{NN}}/GeV)}^{n} #frac{1}{#it{N}_{coll}} #frac{1}{2#pi #it{N}_{inel}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.77, 1.8);
    dummyDirGammapPb2->GetYaxis()->SetRangeUser( 1.2e-5,2.3e17);
    dummyDirGammapPb2->GetXaxis()->SetLabelOffset(-0.012);
    dummyDirGammapPb2->GetXaxis()->SetTickLength(0.025);
    dummyDirGammapPb2->GetYaxis()->SetTickLength(0.025);
    dummyDirGammapPb2->DrawCopy();

        labelXTscaledGammapPb->Draw();
        labelXTscaledGammapPb2->Draw();
        TLatex *labelXTscaledGammaNpPb = new TLatex(0.95,0.86,"n = 4.5");
        SetStyleTLatex( labelXTscaledGammaNpPb, textSizeLabelsPixelDirGam, 4, 1, 43, kTRUE, 31);
        labelXTscaledGammaNpPb->Draw();

        // Plotting ALICE
        for (Int_t ncent = 0; ncent < 5; ncent++){
            if (graphInvYieldDirGammaTotpPb5TeVNCollScaledXT[ncent]){
                ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaTotpPb5TeVNCollScaledXT[ncent]);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotpPb5TeVNCollScaledXT[ncent], markerStyleCombpPbMC[ncent], markerSizeCombpPb[ncent], colorCombpPb[ncent] , colorCombpPb[ncent]);
                graphInvYieldDirGammaTotpPb5TeVNCollScaledXT[ncent]->Draw("p,E1,same");

                TMarker* markerCurrent             = CreateMarkerFromGraph(graphInvYieldDirGammaTotpPb5TeVNCollScaledXT[ncent], xLegALICEXTpPb[ncent]-0.015 , yLegALICEXTpPb[ncent],1);
                markerCurrent->SetNDC(kTRUE);
                markerCurrent->Draw("same,p");
                TLatex *labelCurrent               = new TLatex(xLegALICEXTpPb[ncent] , yLegALICEXTpPb[ncent], labelXTALICEpPb[ncent]);
                SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
                labelCurrent->Draw();

            }
            if (graphInvYieldDirGammaTotArpPb5TeV[ncent]){
                cout << "Pb:    " << ncent << endl;
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpPb5TeVNCollScaledXT[ncent] , 1, 3, colorCombpPb[ncent], colorCombpPb[ncent], 1.8, kTRUE);
                graphInvYieldDirGammaTotArpPb5TeVNCollScaledXT[ncent]->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpPb5TeVNCollScaledXT[ncent], 0.00001);
                graphInvYieldDirGammaTotArpPb5TeVNCollScaledXT[ncent]->Print();
            }
        }
             // plotting PHENIX pAu
        for (Int_t centAu = 0; centAu < 2; centAu++){
            if (graphYieldDirGammapAu200GeVTotNCollScaledXT[centAu]){
                DrawGammaSetMarkerTGraphErr(graphYieldDirGammapAu200GeVTotNCollScaledXT[centAu], markerStyleCombMCpAu200[centAu], markerSizeCombpAu200[centAu], colorCombpAu200[centAu], colorCombpAu200[centAu] );
                graphYieldDirGammapAu200GeVTotNCollScaledXT[centAu]->Draw("p,E1,same");

                TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammapAu200GeVTotNCollScaledXT[centAu], xLegPHENIXpAu200XT[centAu]-0.015 , yLegPHENIXpAu200XT[centAu],1);
                markerCurrent->SetNDC(kTRUE);
                markerCurrent->Draw("same,p");
                TLatex *labelCurrent               = new TLatex(xLegPHENIXpAu200XT[centAu] , yLegPHENIXpAu200XT[centAu], labelXTPHENIXpAu200[centAu]);
                SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
                labelCurrent->Draw();
            }
            if (graphYieldDirGammapAu200GeVTotArNCollScaledXT[centAu]){
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammapAu200GeVTotArNCollScaledXT[centAu] , 1, 3, colorCombpAu200[centAu], colorCombpAu200[centAu], 1.8, kTRUE);
                graphYieldDirGammapAu200GeVTotArNCollScaledXT[centAu]->Draw(">,same");
                graphYieldDirGammapAu200GeVTotArNCollScaledXT[centAu]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammapAu200GeVTotArNCollScaledXT[centAu], 0.001);
            }
        }
            // plotting PHENIX dAu
        if (graphYieldDirGammadAu200GeVTotNCollScaledXT){
            DrawGammaSetMarkerTGraphErr(graphYieldDirGammadAu200GeVTotNCollScaledXT, markerStyleCombMCdAu200, markerSizeCombdAu200, colorCombdAu200, colorCombdAu200 );
            graphYieldDirGammadAu200GeVTotNCollScaledXT->Draw("p,E1,same");

            TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammadAu200GeVTotNCollScaledXT, xLegPHENIXdAu200XT-0.015 , yLegPHENIXdAu200XT,1);
            markerCurrent->SetNDC(kTRUE);
            markerCurrent->Draw("same,p");
            TLatex *labelCurrent               = new TLatex(xLegPHENIXdAu200XT , yLegPHENIXdAu200XT, labelXTPHENIXdAu200);
            SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
            labelCurrent->Draw();
        }
        if (graphYieldDirGammadAu200GeVTotArNCollScaledXT){
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammadAu200GeVTotArNCollScaledXT , 1, 3, colorCombdAu200, colorCombdAu200, 1.8, kTRUE);
            graphYieldDirGammadAu200GeVTotArNCollScaledXT->Draw(">,same");
            graphYieldDirGammadAu200GeVTotArNCollScaledXT->Print();
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammadAu200GeVTotArNCollScaledXT, 0.001);
        }
            // plotting ATLAS pPb 8 TeV
        if (graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaledXT){
            DrawGammaSetMarkerTGraphErr(graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaledXT, markerStyleCombpPb8ATLASMC, markerSizeCombpPb8ATLAS, colorCombpPb8ATLAS, colorCombpPb8ATLAS );
            graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaledXT->Draw("p,E1,same");

            TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammapPb8TeV8000GeVATLASTotNCollScaledXT, xLegpPb8TeVATLASXT-0.015 , yLegpPb8TeVATLASXT,1);
            markerCurrent->SetNDC(kTRUE);
            markerCurrent->Draw("same,p");
            TLatex *labelCurrent               = new TLatex(xLegpPb8TeVATLASXT , yLegpPb8TeVATLASXT, labelXTpPb8TeVATLAS);
            SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
            labelCurrent->Draw();
        }
        if (graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaledXT){
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaledXT , 1, 3, colorCombpPb8ATLAS, colorCombpPb8ATLAS, 1.8, kTRUE);
            graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaledXT->Draw(">,same");
            graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaledXT->Print();
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammapPb8TeV8000GeVATLASTotArNCollScaledXT, 0.001);
        }
    canvasDirGammapPb->Print(Form("%s/DirGammaSpectra_SmallSystems_NCollScaled_xT.%s",outputDir.Data(),suffix.Data()));
    canvasDirGammapPb->Print(Form("%s/DirGammaSpectra_SmallSystems_NCollScaled_xT.pdf",outputDir.Data()));
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    
    TCanvas *canvasDirGammaPb = new TCanvas("canvasDirGammaPb","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasDirGammaPb, 0.18, 0.01, 0.01, 0.08);
    canvasDirGammaPb->SetLogy();
    canvasDirGammaPb->SetLogx();

    TH1D* dummyDirGamma2 = new TH1D("dummyDirGamma", "dummyDirGamma", 1000, 0.09, 120.);
    SetStyleHistoTH1ForGraphs( dummyDirGamma2, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{coll}} #frac{1}{2#pi #it{N}_{inel}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.78, 1.85);
    dummyDirGamma2->GetYaxis()->SetRangeUser( 1.2e-15,4e2);
    dummyDirGamma2->GetXaxis()->SetLabelOffset(-0.012);
    dummyDirGamma2->GetXaxis()->SetTickLength(0.025);
    dummyDirGamma2->GetYaxis()->SetTickLength(0.025);
//     dummyDirGamma2->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDirGamma2->DrawCopy();

        TLatex *labelXTscaledGammaPb = new TLatex(0.95,0.94,"A + A #rightarrow #gamma + X");
        SetStyleTLatex( labelXTscaledGammaPb, textSizeLabelsPixelDirGam, 4, 1, 43, kTRUE, 31);
        labelXTscaledGammaPb->Draw();
        TLatex *labelXTscaledGammaPb2 = new TLatex(0.95,0.90,"y #approx 0");
        SetStyleTLatex( labelXTscaledGammaPb2, textSizeLabelsPixelDirGam, 4, 1, 43, kTRUE, 31);
        labelXTscaledGammaPb2->Draw();

        Int_t nRowsGammaPb              = 3;
        Int_t nRowsGammaAu200           = 4;
        Int_t nRowsGammaCu200           = 2;
        Int_t nRowsGammaAu62            = 3;
        Int_t nRowsGammaAu39            = 2;
        Int_t nRowsGammaPb17            = 2;
        Int_t nColumnsGammaPb           = 3;
        Int_t nColumnsGammaCu           = 2;
        Int_t nColumnsGammaAu62         = 2;
        Int_t nColumnsGammaAu39         = 1;
        Double_t yMinLeg                = 0.1;
        Double_t xMinLeg                = 0.22;
        Double_t widthColumn            = 0.18;
        TLegend* legendInvYieldGammaPb      = GetAndSetLegend2( xMinLeg, yMinLeg, xMinLeg+widthColumn*nColumnsGammaPb, yMinLeg+(nRowsGammaPb*textsizeLabelsDirGamma*0.5) ,0.5*textsizeLabelsDirGamma,
                                                                nColumnsGammaPb, Form("#scale[1.25]{%s}",collisionSystemPbPb2760GeV.Data()), 42, 0.15);
        TLegend* legendInvYieldGammaAu200   = GetAndSetLegend2( xMinLeg, yMinLeg+(nRowsGammaPb)*textsizeLabelsDirGamma*0.5+0.01, xMinLeg+widthColumn*nColumnsGammaPb,
                                                                yMinLeg+((nRowsGammaPb+nRowsGammaAu200)*textsizeLabelsDirGamma*0.5)+0.01 ,0.5*textsizeLabelsDirGamma,
                                                                nColumnsGammaPb, Form("#scale[1.25]{%s}",collisionSystemAuAu200GeV.Data()), 42, 0.15);
        TLegend* legendInvYieldGammaCu200   = GetAndSetLegend2( xMinLeg, yMinLeg+(nRowsGammaPb+nRowsGammaAu200)*textsizeLabelsDirGamma*0.5+0.02, xMinLeg+widthColumn*nColumnsGammaCu,
                                                                yMinLeg+((nRowsGammaPb+nRowsGammaAu200+nRowsGammaCu200)*textsizeLabelsDirGamma*0.5)+0.02 ,0.5*textsizeLabelsDirGamma,
                                                                nColumnsGammaCu, Form("#scale[1.25]{%s}",collisionSystemCuCu200GeV.Data()), 42, 0.15);
        TLegend* legendInvYieldGammaAu62    = GetAndSetLegend2( xMinLeg, yMinLeg+(nRowsGammaPb+nRowsGammaAu200+nRowsGammaCu200)*textsizeLabelsDirGamma*0.5+0.03, xMinLeg+widthColumn*nColumnsGammaAu62,
                                                                yMinLeg+((nRowsGammaPb+nRowsGammaAu200+nRowsGammaCu200+nRowsGammaAu62)*textsizeLabelsDirGamma*0.5)+0.03 ,0.5*textsizeLabelsDirGamma,
                                                                nColumnsGammaAu62, Form("#scale[1.25]{%s}",collisionSystemAuAu62GeV.Data()), 42, 0.15);
        TLegend* legendInvYieldGammaAu39    = GetAndSetLegend2( xMinLeg, yMinLeg+(nRowsGammaPb+nRowsGammaAu200+nRowsGammaCu200+nRowsGammaAu62)*textsizeLabelsDirGamma*0.5+0.04, xMinLeg+widthColumn*nColumnsGammaAu39,
                                                                yMinLeg+((nRowsGammaPb+nRowsGammaAu200+nRowsGammaCu200+nRowsGammaAu62+nRowsGammaAu39)*textsizeLabelsDirGamma*0.5)+0.04 ,0.5*textsizeLabelsDirGamma,
                                                                nColumnsGammaAu39,  Form("#scale[1.25]{%s}",collisionSystemAuAu39GeV.Data()), 42, 0.15);
        TLegend* legendInvYieldGammaPb17    = GetAndSetLegend2( xMinLeg, yMinLeg+(nRowsGammaPb+nRowsGammaAu200+nRowsGammaCu200+nRowsGammaAu62+nRowsGammaAu39)*textsizeLabelsDirGamma*0.5+0.05,
                                                                xMinLeg+widthColumn*nColumnsGammaAu39,
                                                                yMinLeg+((nRowsGammaPb+nRowsGammaAu200+nRowsGammaCu200+nRowsGammaAu62+nRowsGammaAu39+nRowsGammaPb17)*textsizeLabelsDirGamma*0.5)+0.05 ,
                                                                0.5*textsizeLabelsDirGamma, nColumnsGammaAu39, Form("#scale[1.25]{%s}",collisionSystemPbPb17GeV.Data()), 42, 0.15);

        for (Int_t centPb = 0; centPb < 3; centPb++){
            if (graphInvYieldDirGammaTotPbPbNCollScaled[centPb]){
                ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaTotPbPbNCollScaled[centPb]);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotPbPbNCollScaled[centPb], markerStyleCombMC[centPb], markerSizeComb[centPb], colorCombPbALICE[centPb] , colorCombPbALICE[centPb]);
                graphInvYieldDirGammaTotPbPbNCollScaled[centPb]->Draw("p,E1,same");
                legendInvYieldGammaPb->AddEntry(graphInvYieldDirGammaTotPbPbNCollScaled[centPb], labelALICE[centPb].Data(), "p");
            }
            if (graphInvYieldDirGammaTotArPbPb[centPb]){
                cout << "Pb:    " << centPb << endl;
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArPbPbNCollScaled[centPb] , 1, 3, colorCombPbALICE[centPb], colorCombPbALICE[centPb], 1.8, kTRUE);
                graphInvYieldDirGammaTotArPbPbNCollScaled[centPb]->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPbNCollScaled[centPb], 0.1);
                graphInvYieldDirGammaTotArPbPbNCollScaled[centPb]->Print();
            }
        }
        // plotting CMS
        for (Int_t centPb = 0; centPb < 3; centPb++){
            if (graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]){
                ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]);
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb], markerStyleCombMCCMS[centPb], markerSizeCombCMS[centPb], colorCombCMS[centPb], colorCombCMS[centPb] );
                graphYieldDirGammaCMSPbPb2760GeVTotNCollScaled[centPb]->Draw("p,E1,same");
                legendInvYieldGammaPb->AddEntry(graphInvYieldDirGammaTotPbPbNCollScaled[centPb], labelCMS[centPb].Data(), "p");
            }
            if (graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb]){
                ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb]);
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb] , 1, 3, colorCombCMS[centPb], colorCombCMS[centPb], 1.8, kTRUE);
                graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb]->Draw(">,same");
                graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaled[centPb], 0.1);
            }

        }

        // plotting PHENIX AuAu 200
        for (Int_t centAu = 0; centAu < 5; centAu++){
            for (Int_t meth = 0; meth < 2; meth++){
//                 cout << "trying to plot: " << centAu << "\t" << meth << endl;
                if (graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]){
                    //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVStat[centAu][meth]);
                    DrawGammaSetMarkerTGraphErr(graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth], markerStyleCombMCAuAu200[centAu], markerSizeCombAuAu200[centAu], colorCombAuAu200[centAu], colorCombAuAu200[centAu] );
                    graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth]->Draw("p,E1,same");
                    if (meth == 0)legendInvYieldGammaAu200->AddEntry(graphYieldDirGammaAuAu200GeVTotNCollScaled[centAu][meth], labelPHENIX200[centAu].Data(), "p");
                }
                if (graphYieldDirGammaAuAu200GeVTotArNCollScaled[centAu][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu200GeVTotArNCollScaled[centAu][meth] , 1, 3, colorCombAuAu200[centAu], colorCombAuAu200[centAu], 1.8, kTRUE);
                    graphYieldDirGammaAuAu200GeVTotArNCollScaled[centAu][meth]->Draw(">,same");
                    graphYieldDirGammaAuAu200GeVTotArNCollScaled[centAu][meth]->Print();
                    PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaAuAu200GeVTotArNCollScaled[centAu][meth], 0.1);
                }
            }
        }
        // plotting STAR AuAu 200
        for (Int_t centAu = 0; centAu < 4; centAu++){
            //             cout << "trying to plot: " << centAu << "\t" << endl;
            if (graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]){
                ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]);
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu], markerStyleCombMCAuAu200STAR[centAu], markerSizeCombAuAu200STAR[centAu], colorCombAuAu200STAR[centAu], colorCombAuAu200STAR[centAu]);
                graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu]->Draw("p,E1,same");
                legendInvYieldGammaAu200->AddEntry(graphYieldDirGammaAuAu200GeVSTARTotNCollScaled[centAu], labelSTAR200[centAu].Data(), "p");
            }
            if (graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu]){
                ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu]);
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu] , 1, 3, colorCombAuAu200STAR[centAu], colorCombAuAu200STAR[centAu], 1.8, kTRUE);
                graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu]->Draw(">,same");
                graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaAuAu200GeVSTARTotArNCollScaled[centAu], 0.1);
            }
        }

        // plotting PHENIX CuCu 200
        for (Int_t centCu = 0; centCu < 2; centCu++){
//             cout << "trying to plot: " << centCu  << endl;
            if (graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]){
                //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaCuCu200GeVStat[centCu]);
                DrawGammaSetMarkerTGraphErr(graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu], markerStyleCombMCCuCu200[centCu], markerSizeCombCuCu200[centCu], colorCombCuCu200[centCu], colorCombCuCu200[centCu] );
                graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu]->Draw("p,E1,same");
                legendInvYieldGammaCu200->AddEntry(graphYieldDirGammaCuCu200GeVTotNCollScaled[centCu], labelPHENIXCu200[centCu].Data(), "p");
            }
            if (graphYieldDirGammaCuCu200GeVTotArNCollScaled[centCu]){
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaCuCu200GeVTotArNCollScaled[centCu] , 1, 3, colorCombCuCu200[centCu], colorCombCuCu200[centCu], 1.8, kTRUE);
                graphYieldDirGammaCuCu200GeVTotArNCollScaled[centCu]->Draw(">,same");
                graphYieldDirGammaCuCu200GeVTotArNCollScaled[centCu]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaCuCu200GeVTotArNCollScaled[centCu]);
            }
        }

        // plotting PHENIX AuAu 62
        for (Int_t centAu = 0; centAu < 3; centAu++){
//             cout << "trying to plot: " << centAu << "\t" << endl;
            if (graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]){
                //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVStat[centAu]);
                DrawGammaSetMarkerTGraphErr(graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu], markerStyleCombMCAuAu62[centAu], markerSizeCombAuAu62[centAu], colorCombAuAu62[centAu], colorCombAuAu62[centAu] );
                graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu]->Draw("p,E1,same");
                legendInvYieldGammaAu62->AddEntry(graphYieldDirGammaAuAu62GeVTotNCollScaled[centAu], labelPHENIXAu62[centAu].Data(), "p");
            }
            if (graphYieldDirGammaAuAu62GeVTotArNCollScaled[centAu]){
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu62GeVTotArNCollScaled[centAu] , 1, 3, colorCombAuAu62[centAu], colorCombAuAu62[centAu], 1.8, kTRUE);
                graphYieldDirGammaAuAu62GeVTotArNCollScaled[centAu]->Draw(">,same");
                graphYieldDirGammaAuAu62GeVTotArNCollScaled[centAu]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaAuAu62GeVTotArNCollScaled[centAu], 0.03);
            }
        }

        // plotting PHENIX 39
        if (graphYieldDirGammaAuAu39GeVTotNCollScaled){
            cout << "have 39GeV" << endl;
            graphYieldDirGammaAuAu39GeVTotNCollScaled->Print();
            //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVStat);
            DrawGammaSetMarkerTGraphErr(graphYieldDirGammaAuAu39GeVTotNCollScaled, markerStyleCombMCAuAu39, markerSizeCombAuAu39, colorCombAuAu39, colorCombAuAu39 );
            graphYieldDirGammaAuAu39GeVTotNCollScaled->Draw("p,E1,same");
            legendInvYieldGammaAu39->AddEntry(graphYieldDirGammaAuAu39GeVTotNCollScaled, labelPHENIXAu39.Data(), "p");
        }
        if (graphYieldDirGammaAuAu39GeVTotArNCollScaled){
            cout << "39 GeV" << endl;
            graphYieldDirGammaAuAu39GeVTotArNCollScaled->Print();
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu39GeVTotArNCollScaled , 1, 3, colorCombAuAu39, colorCombAuAu39, 1.8, kTRUE);
            graphYieldDirGammaAuAu39GeVTotArNCollScaled->Draw(">,same");
            graphYieldDirGammaAuAu39GeVTotArNCollScaled->Print();
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaAuAu39GeVTotArNCollScaled, 0.03);
        }
        // plotting WA98 17.2
        if (graphYieldDirGammaWA98Pb17GeVTotArNCollScaled){
            cout << "have 39GeV" << endl;
            graphYieldDirGammaWA98Pb17GeVTotNCollScaled->Print();
            ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaWA98Pb17GeVTotArNCollScaled);
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaWA98Pb17GeVTotNCollScaled, markerStyleCombMCPbPb17, markerSizeCombPbPb17, colorCombPbPb17, colorCombPbPb17 );
            graphYieldDirGammaWA98Pb17GeVTotNCollScaled->Draw("p,E1,same");
            legendInvYieldGammaPb17->AddEntry(graphYieldDirGammaWA98Pb17GeVTotNCollScaled, labelWA98.Data(), "p");
        }
        if (graphYieldDirGammaWA98Pb17GeVTotArNCollScaled){
            cout << "39 GeV" << endl;
            graphYieldDirGammaWA98Pb17GeVTotArNCollScaled->Print();
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaWA98Pb17GeVTotArNCollScaled , 1, 3, colorCombPbPb17, colorCombPbPb17, 1.8, kTRUE);
            graphYieldDirGammaWA98Pb17GeVTotArNCollScaled->Draw(">,same");
            graphYieldDirGammaWA98Pb17GeVTotArNCollScaled->Print();
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaWA98Pb17GeVTotArNCollScaled, 0.02);
        }

        legendInvYieldGammaPb->Draw();
        legendInvYieldGammaAu200->Draw();
        legendInvYieldGammaCu200->Draw();
        legendInvYieldGammaAu62->Draw();
        legendInvYieldGammaAu39->Draw();
        legendInvYieldGammaPb17->Draw();

    canvasDirGammaPb->Print(Form("%s/DirGammaSpectra_NCollScaled.%s",outputDir.Data(),suffix.Data()));
    canvasDirGammaPb->Print(Form("%s/DirGammaSpectra_NCollScaled.pdf",outputDir.Data()));

    canvasDirGamma->cd();
    TH1D* dummyDirGamma3 = new TH1D("dummyDirGamma", "dummyDirGamma", 1000, 0.0003, 3);
    SetStyleHistoTH1ForGraphs( dummyDirGamma3, "#it{x}_{T} = 2#it{p}_{T}/#sqrt{#it{s}_{_{NN}}}", "#sqrt{(#it{s}_{_{NN}}/GeV)}^{n} #frac{1}{#it{N}_{coll}} #frac{1}{2#pi #it{N}_{inel}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.77, 1.8);
    dummyDirGamma3->GetYaxis()->SetRangeUser( 1.2e-5,2.3e13);
    dummyDirGamma3->GetXaxis()->SetLabelOffset(-0.012);
    dummyDirGamma3->GetXaxis()->SetTickLength(0.025);
    dummyDirGamma3->GetYaxis()->SetTickLength(0.025);
    dummyDirGamma3->DrawCopy();

        labelXTscaledGammaPb->Draw();
        labelXTscaledGammaPb2->Draw();
        TLatex *labelXTscaledGammaN = new TLatex(0.95,0.86,"n = 4.5");
        SetStyleTLatex( labelXTscaledGammaN, textSizeLabelsPixelDirGam, 4, 1, 43, kTRUE, 31);
        labelXTscaledGammaN->Draw();

        // Plotting ALICE
        for (Int_t centPb = 0; centPb < 3; centPb++){
            if (graphInvYieldDirGammaTotPbPbNCollScaledXT[centPb]){
                ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaTotPbPbNCollScaledXT[centPb]);
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotPbPbNCollScaledXT[centPb], markerStyleCombMC[centPb], markerSizeComb[centPb], colorCombPbALICE[centPb] , colorCombPbALICE[centPb]);
                graphInvYieldDirGammaTotPbPbNCollScaledXT[centPb]->Draw("p,E1,same");

                TMarker* markerCurrent             = CreateMarkerFromGraph(graphInvYieldDirGammaTotPbPbNCollScaledXT[centPb], xLegALICEXT[centPb]-0.015 , yLegALICEXT[centPb],1);
                markerCurrent->SetNDC(kTRUE);
                markerCurrent->Draw("same,p");
                TLatex *labelCurrent               = new TLatex(xLegALICEXT[centPb] , yLegALICEXT[centPb], labelXTALICE[centPb]);
                SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
                labelCurrent->Draw();

            }
            if (graphInvYieldDirGammaTotArPbPb[centPb]){
                cout << "Pb:    " << centPb << endl;
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArPbPbNCollScaledXT[centPb] , 1, 3, colorCombPbALICE[centPb], colorCombPbALICE[centPb], 1.8, kTRUE);
                graphInvYieldDirGammaTotArPbPbNCollScaledXT[centPb]->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPbNCollScaledXT[centPb], 0.0001);
                graphInvYieldDirGammaTotArPbPbNCollScaledXT[centPb]->Print();
            }
        }
        // Plotting CMS
        for (Int_t centPb = 0; centPb < 3; centPb++){
            if (graphYieldDirGammaCMSPbPb2760GeVTotNCollScaledXT[centPb]){
                ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaCMSPbPb2760GeVTotNCollScaledXT[centPb]);
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaCMSPbPb2760GeVTotNCollScaledXT[centPb], markerStyleCombMCCMS[centPb], markerSizeCombCMS[centPb], colorCombCMS[centPb], colorCombCMS[centPb] );
                graphYieldDirGammaCMSPbPb2760GeVTotNCollScaledXT[centPb]->Draw("p,E1,same");
                TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammaCMSPbPb2760GeVTotNCollScaledXT[centPb], xLegCMSXT[centPb]-0.015 , yLegCMSXT[centPb],1);
                markerCurrent->SetNDC(kTRUE);
                markerCurrent->Draw("same,p");
                TLatex *labelCurrent               = new TLatex(xLegCMSXT[centPb] , yLegCMSXT[centPb], labelXTCMS[centPb]);
                SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
                labelCurrent->Draw();
            }
            if (graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaledXT[centPb]){
                ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaledXT[centPb]);
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaledXT[centPb] , 1, 3, colorCombCMS[centPb], colorCombCMS[centPb], 1.8, kTRUE);
                graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaledXT[centPb]->Draw(">,same");
                graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaledXT[centPb]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaCMSPbPb2760GeVTotArNCollScaledXT[centPb], 0.005);
            }

        }

        // plotting PHENIX
        for (Int_t centAu = 0; centAu < 5; centAu++){
            for (Int_t meth = 0; meth < 2; meth++){
//                 cout << "trying to plot: " << centAu << "\t" << meth << endl;
                if (graphYieldDirGammaAuAu200GeVTotNCollScaledXT[centAu][meth]){
                    //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVStat[centAu][meth]);
                    DrawGammaSetMarkerTGraphErr(graphYieldDirGammaAuAu200GeVTotNCollScaledXT[centAu][meth], markerStyleCombMCAuAu200[centAu], markerSizeCombAuAu200[centAu], colorCombAuAu200[centAu], colorCombAuAu200[centAu] );
                    graphYieldDirGammaAuAu200GeVTotNCollScaledXT[centAu][meth]->Draw("p,E1,same");

                    if (meth == 0){
                        TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammaAuAu200GeVTotNCollScaledXT[centAu][meth], xLegPHENIX200XT[centAu]-0.015 , yLegPHENIX200XT[centAu],1);
                        markerCurrent->SetNDC(kTRUE);
                        markerCurrent->Draw("same,p");
                        TLatex *labelCurrent               = new TLatex(xLegPHENIX200XT[centAu] , yLegPHENIX200XT[centAu], labelXTPHENIX200[centAu]);
                        SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
                        labelCurrent->Draw();
                    }
                }
                if (graphYieldDirGammaAuAu200GeVTotArNCollScaledXT[centAu][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu200GeVTotArNCollScaledXT[centAu][meth] , 1, 3, colorCombAuAu200[centAu], colorCombAuAu200[centAu], 1.8, kTRUE);
                    graphYieldDirGammaAuAu200GeVTotArNCollScaledXT[centAu][meth]->Draw(">,same");
                    graphYieldDirGammaAuAu200GeVTotArNCollScaledXT[centAu][meth]->Print();
                    PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaAuAu200GeVTotArNCollScaledXT[centAu][meth], 0.001);
                }
            }
        }
        // plotting STAR
        for (Int_t centAu = 0; centAu < 5; centAu++){
//             cout << "trying to plot: " << centAu << "\t" << endl;
            if (graphYieldDirGammaAuAu200GeVSTARTotNCollScaledXT[centAu]){
//                 graphYieldDirGammaAuAu200GeVSTARTotNCollScaledXT[centAu]->Print();
                ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVSTARTotNCollScaledXT[centAu]);
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu200GeVSTARTotNCollScaledXT[centAu], markerStyleCombMCAuAu200STAR[centAu], markerSizeCombAuAu200STAR[centAu], colorCombAuAu200STAR[centAu], colorCombAuAu200STAR[centAu]);
                graphYieldDirGammaAuAu200GeVSTARTotNCollScaledXT[centAu]->Draw("p,E1,same");
                TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammaAuAu200GeVSTARTotNCollScaledXT[centAu], xLegSTAR200XT[centAu]-0.015 , yLegSTAR200XT[centAu],1);
                markerCurrent->SetNDC(kTRUE);
                markerCurrent->Draw("same,p");
                TLatex *labelCurrent               = new TLatex(xLegSTAR200XT[centAu] , yLegSTAR200XT[centAu], labelXTSTAR200[centAu]);
                SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
                labelCurrent->Draw();
            }
            if (graphYieldDirGammaAuAu200GeVSTARTotArNCollScaledXT[centAu]){
                ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVSTARTotArNCollScaledXT[centAu]);
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu200GeVSTARTotArNCollScaledXT[centAu] , 1, 3, colorCombAuAu200STAR[centAu], colorCombAuAu200STAR[centAu], 1.8, kTRUE);
                graphYieldDirGammaAuAu200GeVSTARTotArNCollScaledXT[centAu]->Draw(">,same");
                graphYieldDirGammaAuAu200GeVSTARTotArNCollScaledXT[centAu]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaAuAu200GeVSTARTotArNCollScaledXT[centAu], 0.005);
            }
        }
        for (Int_t centCu = 0; centCu < 2; centCu++){
//             cout << "trying to plot: " << centCu  << endl;
            if (graphYieldDirGammaCuCu200GeVTotNCollScaledXT[centCu]){
                //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaCuCu200GeVStat[centCu]);
                DrawGammaSetMarkerTGraphErr(graphYieldDirGammaCuCu200GeVTotNCollScaledXT[centCu], markerStyleCombMCCuCu200[centCu], markerSizeCombCuCu200[centCu], colorCombCuCu200[centCu], colorCombCuCu200[centCu] );
                graphYieldDirGammaCuCu200GeVTotNCollScaledXT[centCu]->Draw("p,E1,same");

                TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammaCuCu200GeVTotNCollScaledXT[centCu], xLegPHENIXCuCu200CXT[centCu]-0.015 , yLegPHENIXCuCu200XT[centCu],1);
                markerCurrent->SetNDC(kTRUE);
                markerCurrent->Draw("same,p");
                TLatex *labelCurrent               = new TLatex(xLegPHENIXCuCu200CXT[centCu] , yLegPHENIXCuCu200XT[centCu], labelXTPHENIXCuCu200[centCu]);
                SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
                labelCurrent->Draw();

            }
            if (graphYieldDirGammaCuCu200GeVTotArNCollScaledXT[centCu]){
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaCuCu200GeVTotArNCollScaledXT[centCu] , 1, 3, colorCombCuCu200[centCu], colorCombCuCu200[centCu], 1.8, kTRUE);
                graphYieldDirGammaCuCu200GeVTotArNCollScaledXT[centCu]->Draw(">,same");
                graphYieldDirGammaCuCu200GeVTotArNCollScaledXT[centCu]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaCuCu200GeVTotArNCollScaledXT[centCu], 0.001);
            }
        }

        // plotting PHENIX
        for (Int_t centAu = 0; centAu < 3; centAu++){
//             cout << "trying to plot: " << centAu << "\t" << endl;
            if (graphYieldDirGammaAuAu62GeVTotNCollScaledXT[centAu]){
                //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVStat[centAu]);
                DrawGammaSetMarkerTGraphErr(graphYieldDirGammaAuAu62GeVTotNCollScaledXT[centAu], markerStyleCombMCAuAu62[centAu], markerSizeCombAuAu62[centAu], colorCombAuAu62[centAu], colorCombAuAu62[centAu] );
                graphYieldDirGammaAuAu62GeVTotNCollScaledXT[centAu]->Draw("p,E1,same");
                TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammaAuAu62GeVTotNCollScaledXT[centAu], xLegPHENIX62XT[centAu]-0.015 , yLegPHENIX62XT[centAu],1);
                markerCurrent->SetNDC(kTRUE);
                markerCurrent->Draw("same,p");
                TLatex *labelCurrent               = new TLatex(xLegPHENIX62XT[centAu] , yLegPHENIX62XT[centAu], labelXTPHENIX62[centAu]);
                SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
                labelCurrent->Draw();
            }
            if (graphYieldDirGammaAuAu62GeVTotArNCollScaledXT[centAu]){
                DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu62GeVTotArNCollScaledXT[centAu] , 1, 3, colorCombAuAu62[centAu], colorCombAuAu62[centAu], 1.8, kTRUE);
                graphYieldDirGammaAuAu62GeVTotArNCollScaledXT[centAu]->Draw(">,same");
//                 graphYieldDirGammaAuAu62GeVTotArNCollScaledXT[centAu]->Print();
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaAuAu62GeVTotArNCollScaledXT[centAu], 0.001);
            }
        }


        if (graphYieldDirGammaAuAu39GeVTotNCollScaledXT){
            //                     ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaAuAu200GeVStat);
            DrawGammaSetMarkerTGraphErr(graphYieldDirGammaAuAu39GeVTotNCollScaledXT, markerStyleCombMCAuAu39, markerSizeCombAuAu39, colorCombAuAu39, colorCombAuAu39 );
            graphYieldDirGammaAuAu39GeVTotNCollScaledXT->Draw("p,E1,same");
            TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammaAuAu39GeVTotNCollScaledXT, xLegPHENIX39XT-0.015 , yLegPHENIX39XT,1);
            markerCurrent->SetNDC(kTRUE);
            markerCurrent->Draw("same,p");
            TLatex *labelCurrent               = new TLatex(xLegPHENIX39XT , yLegPHENIX39XT, labelXTPHENIX39);
            SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
            labelCurrent->Draw();
        }
        if (graphYieldDirGammaAuAu39GeVTotArNCollScaledXT){
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaAuAu39GeVTotArNCollScaledXT , 1, 3, colorCombAuAu39, colorCombAuAu39, 1.8, kTRUE);
            graphYieldDirGammaAuAu39GeVTotArNCollScaledXT->Draw(">,same");
            graphYieldDirGammaAuAu39GeVTotArNCollScaledXT->Print();
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaAuAu39GeVTotArNCollScaledXT, 0.0001);
        }

        if (graphYieldDirGammaWA98Pb17GeVTotNCollScaledXT){
            ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaWA98Pb17GeVTotNCollScaledXT);
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaWA98Pb17GeVTotNCollScaledXT, markerStyleCombMCPbPb17, markerSizeCombPbPb17, colorCombPbPb17, colorCombPbPb17 );
            graphYieldDirGammaWA98Pb17GeVTotNCollScaledXT->Draw("p,E1,same");
            TMarker* markerCurrent             = CreateMarkerFromGraph(graphYieldDirGammaWA98Pb17GeVTotNCollScaledXT, xLegWA98XT-0.015 , yLegWA98XT,1);
            markerCurrent->SetNDC(kTRUE);
            markerCurrent->Draw("same,p");
            TLatex *labelCurrent               = new TLatex(xLegWA98XT , yLegWA98XT, labelXTWA98);
            SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
            labelCurrent->Draw();
        }
        if (graphYieldDirGammaWA98Pb17GeVTotArNCollScaledXT){
            ProduceGraphAsymmWithoutXErrors(graphYieldDirGammaWA98Pb17GeVTotArNCollScaledXT);
            DrawGammaSetMarkerTGraphAsym(graphYieldDirGammaWA98Pb17GeVTotArNCollScaledXT , 1, 3, colorCombPbPb17, colorCombPbPb17, 1.8, kTRUE);
            graphYieldDirGammaWA98Pb17GeVTotArNCollScaledXT->Draw(">,same");
            graphYieldDirGammaWA98Pb17GeVTotArNCollScaledXT->Print();
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphYieldDirGammaWA98Pb17GeVTotArNCollScaledXT, 0.001);
        }

    canvasDirGamma->Print(Form("%s/DirGammaSpectra_NCollScaled_xT.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectra_NCollScaled_xT.pdf",outputDir.Data()));



    dummyDirGamma->GetYaxis()->SetRangeUser( 1.2e-10,8.5);
    dummyDirGamma->DrawCopy();

        if (graphInvYieldDirGammaSyspp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp2760GeV , markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspp2760GeV->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatpp2760GeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatpp2760GeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp2760GeV, markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV);
            graphInvYieldDirGammaStatpp2760GeV->Draw("p,E1Z,same");
        }

        if (graphInvYieldDirGammaTotArpp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp2760GeV , 1, 3, colorEpp2760GeV, colorEpp2760GeV, 1.8, kTRUE);
            graphInvYieldDirGammaTotArpp2760GeV->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeV);
        }
        if (graphInvYieldDirGammaSyspp8TeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp8TeV , markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspp8TeV->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatpp8TeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatpp8TeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp8TeV, markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV);
            graphInvYieldDirGammaStatpp8TeV->Draw("p,E1Z,same");
        }

        if (graphInvYieldDirGammaTotArpp8TeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp8TeV , 1, 3, colorEpp8TeV, colorEpp8TeV, 1.8, kTRUE);
            graphInvYieldDirGammaTotArpp8TeV->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp8TeV);
        }

        TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.23, 0.11+0.04*4, "#gamma_{dir} ALICE");
        SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyInvYieldPaperAll->Draw();
        TLegend* legendDirGammaPPonly2 = GetAndSetLegend2(0.22, 0.11+(0*textsizeLabelsDirGamma*0.85), 0.21+0.25, 0.11+(4*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        legendDirGammaPPonly2->AddEntry(graphInvYieldDirGammaSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDirGammaPPonly2->AddEntry((TObject*)0,"norm. unc. 2.5%","");
        legendDirGammaPPonly2->AddEntry(graphInvYieldDirGammaSyspp8TeV,collisionSystempp8TeV.Data(),"pf");
        legendDirGammaPPonly2->AddEntry((TObject*)0,"norm. unc. 2.6%","");
        legendDirGammaPPonly2->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectraPP.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectraPP.pdf",outputDir.Data()));

    TH1D* dummyIncGamma = new TH1D("dummyIncGamma", "dummyIncGamma", 1000, 0., 22.);
    SetStyleHistoTH1ForGraphs( dummyIncGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{inel}} #frac{d^{2}#it{N}_{#gamma_{inc}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.8);
    dummyIncGamma->GetYaxis()->SetRangeUser( 1.2e-9,1.5e1);
    dummyIncGamma->GetXaxis()->SetLabelOffset(-0.015);
    dummyIncGamma->GetXaxis()->SetTickLength(0.025);
    dummyIncGamma->GetYaxis()->SetTickLength(0.025);
    dummyIncGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyIncGamma->DrawCopy();
        if (graphInvYieldIncGammaSyspp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldIncGammaSyspp2760GeV , markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV, widthLinesBoxes, kTRUE);
            graphInvYieldIncGammaSyspp2760GeV->Draw("E2same");
        }
        if (graphInvYieldIncGammaStatpp2760GeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldIncGammaStatpp2760GeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldIncGammaStatpp2760GeV, markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV);
            graphInvYieldIncGammaStatpp2760GeV->Draw("p,E1Z,same");
        }

        if (graphInvYieldIncGammaSyspp8TeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldIncGammaSyspp8TeV , markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV, widthLinesBoxes, kTRUE);
            graphInvYieldIncGammaSyspp8TeV->Draw("E2same");
        }
        if (graphInvYieldIncGammaStatpp8TeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldIncGammaStatpp8TeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldIncGammaStatpp8TeV, markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV);
            graphInvYieldIncGammaStatpp8TeV->Draw("p,E1Z,same");
        }

        TLatex *labelEnergyIncInvYieldPaperAll = new TLatex(0.23, 0.11+0.04*4, "#gamma_{inc} ALICE");
        SetStyleTLatex( labelEnergyIncInvYieldPaperAll, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyIncInvYieldPaperAll->Draw();
        legendDirGammaPPonly2->Draw();


    canvasDirGamma->Print(Form("%s/IncGammaSpectraPP.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/IncGammaSpectraPP.pdf",outputDir.Data()));

    TCanvas *canvasDirGamma2 = new TCanvas("canvasDirGamma2","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasDirGamma2, 0.157, 0.01, 0.01, 0.07);
    canvasDirGamma2->SetLogy();
    canvasDirGamma2->SetLogx();

    TH1D* dummyXSecGamma = new TH1D("dummyXSecGamma", "dummyXSecGamma", 1000, 0.1, 3000.);
    SetStyleHistoTH1ForGraphs( dummyXSecGamma, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.57);
    dummyXSecGamma->GetYaxis()->SetRangeUser( 1.2e-10,1.5e11);
    dummyXSecGamma->GetXaxis()->SetLabelOffset(-0.01);
    dummyXSecGamma->GetXaxis()->SetTickLength(0.025);
    dummyXSecGamma->GetYaxis()->SetTickLength(0.025);
//     dummyXSecGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyXSecGamma->DrawCopy();

    TLatex *labelXTscaledGamma = new TLatex(0.95,0.94,"p + p(#bar{p}) #rightarrow #gamma + X");
    SetStyleTLatex( labelXTscaledGamma, textSizeLabelsPixelDirGam, 4, 1, 43, kTRUE, 31);
    labelXTscaledGamma->Draw();
    TLatex *labelXTscaledGamma2 = new TLatex(0.95,0.90,"y #approx 0");
    SetStyleTLatex( labelXTscaledGamma2, textSizeLabelsPixelDirGam, 4, 1, 43, kTRUE, 31);
    labelXTscaledGamma2->Draw();

        Int_t nRowsGamma            = 8;
        if (enablePrelim)
            nRowsGamma++;
        Int_t nColumnsGamma         = 3;
        TLegend* legendXSecGamma    = GetAndSetLegend2(0.185, 0.095, 0.185+0.19*nColumnsGamma, 0.095+(nRowsGamma*textsizeLabelsDirGamma*0.5) ,0.5*textsizeLabelsDirGamma, nColumnsGamma, "", 42, 0.15);

        if ( enablePrelim ){
            if (graphXSecDirGammaSyspp7TeV){
                DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaSyspp7TeV , markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV, colorEpp7TeV, widthLinesBoxes, kTRUE);
                graphXSecDirGammaSyspp7TeV->Draw("E2same");
            }
            if (graphXSecDirGammaStatpp7TeV){
                ProduceGraphAsymmWithoutXErrors(graphXSecDirGammaStatpp7TeV);
                DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaStatpp7TeV, markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV, colorEpp7TeV);
                graphXSecDirGammaStatpp7TeV->Draw("p,E1Z,same");
            }

            if (graphXSecDirGammaTotArpp7TeV){
                DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaTotArpp7TeV , 1, 3, colorEpp7TeV, colorEpp7TeV, 1.8, kTRUE);
                graphXSecDirGammaTotArpp7TeV->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphXSecDirGammaTotArpp7TeV);
            }
        }
        if (graphXSecDirGammaSyspp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaSyspp2760GeV , markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV, widthLinesBoxes, kTRUE);
            graphXSecDirGammaSyspp2760GeV->Draw("E2same");
        }
        if (graphXSecDirGammaStatpp2760GeV){
            ProduceGraphAsymmWithoutXErrors(graphXSecDirGammaStatpp2760GeV);
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaStatpp2760GeV, markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV);
            graphXSecDirGammaStatpp2760GeV->Draw("p,E1Z,same");
        }

        if (graphXSecDirGammaTotArpp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaTotArpp2760GeV , 1, 3, colorEpp2760GeV, colorEpp2760GeV, 1.8, kTRUE);
            graphXSecDirGammaTotArpp2760GeV->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphXSecDirGammaTotArpp2760GeV);
        }
        if (graphXSecDirGammaSyspp8TeV){
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaSyspp8TeV , markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV, widthLinesBoxes, kTRUE);
            graphXSecDirGammaSyspp8TeV->Draw("E2same");
        }
        if (graphXSecDirGammaStatpp8TeV){
            ProduceGraphAsymmWithoutXErrors(graphXSecDirGammaStatpp8TeV);
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaStatpp8TeV, markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV);
            graphXSecDirGammaStatpp8TeV->Draw("p,E1Z,same");
        }

        if (graphXSecDirGammaTotArpp8TeV){
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaTotArpp8TeV , 1, 3, colorEpp8TeV, colorEpp8TeV, 1.8, kTRUE);
            graphXSecDirGammaTotArpp8TeV->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphXSecDirGammaTotArpp8TeV);
        }

        for (Int_t i = 0; i < nOtherEnergies; i++){
            if (i == 2 & !enablePrelim)
                continue;
            if (i ==1) legendXSecGamma->AddEntry(graphXSecDirGammaStatpp8TeV, Form("%s (%s)", "ALICE", "8 TeV" ), "p");
            if (i ==2 && enablePrelim) legendXSecGamma->AddEntry(graphXSecDirGammaStatpp7TeV, Form("%s (%s)", "ALICE prel.", "7 TeV" ), "p");
            if (i ==7) legendXSecGamma->AddEntry(graphXSecDirGammaStatpp2760GeV, Form("%s (%s)", "ALICE ", "2.76 TeV" ), "p");
            if (xSectionOtherExperiments[i]){
                DrawGammaSetMarkerTGraphAsym(xSectionOtherExperiments[i] , markerStyleOtherExp[i], markerSizeOtherExp[i], colorOtherExp[i], colorOtherExp[i], widthLinesBoxes);
                xSectionOtherExperiments[i]->Draw("p,E1Z,same");
                legendXSecGamma->AddEntry(xSectionOtherExperiments[i], Form("%s (%s)", nameOtherExp[i].Data(), nameOtherEnergy[i].Data() ), "p");
            }
        }
        legendXSecGamma->Draw();

    if (enablePrelim){
        canvasDirGamma2->Print(Form("%s/InvXSecSpectraPP_withPrelim.%s",outputDir.Data(),suffix.Data()));
        canvasDirGamma2->Print(Form("%s/InvXSecSpectraPP_withPrelim.pdf",outputDir.Data()));
    } else {
        canvasDirGamma2->Print(Form("%s/InvXSecSpectraPP.%s",outputDir.Data(),suffix.Data()));
        canvasDirGamma2->Print(Form("%s/InvXSecSpectraPP.pdf",outputDir.Data()));
    }

    TH1D* dummyXSecGamma_xT = new TH1D("dummyXSecGamma_xT", "dummyXSecGamma_xT", 5000, 0.00005, 5.);
    SetStyleHistoTH1ForGraphs( dummyXSecGamma_xT, "#it{x}_{T} = 2#it{p}_{T}/#sqrt{#it{s}}", "#sqrt{(#it{s}/GeV)}^{n} #frac{#it{E}d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.7);
    dummyXSecGamma_xT->GetYaxis()->SetRangeUser( 9e2,1.5e29);
    dummyXSecGamma_xT->GetXaxis()->SetLabelOffset(-0.01);
    dummyXSecGamma_xT->GetXaxis()->SetTickLength(0.025);
    dummyXSecGamma_xT->GetYaxis()->SetTickLength(0.025);
    dummyXSecGamma_xT->DrawCopy();

        labelXTscaledGamma->Draw();
        labelXTscaledGamma2->Draw();

        labelXTscaledGammaN->Draw();

        if ( enablePrelim ){
            if (graphXSecDirGammaSyspp7TeV_xT){
                DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaSyspp7TeV_xT , markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV, colorEpp7TeV, widthLinesBoxes, kTRUE);
                graphXSecDirGammaSyspp7TeV_xT->Draw("E2same");
            }
            if (graphXSecDirGammaStatpp7TeV_xT){
                ProduceGraphAsymmWithoutXErrors(graphXSecDirGammaStatpp7TeV_xT);
                DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaStatpp7TeV_xT, markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV, colorEpp7TeV);
                graphXSecDirGammaStatpp7TeV_xT->Draw("p,E1Z,same");
            }

            if (graphXSecDirGammaTotArpp7TeV_xT){
                DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaTotArpp7TeV_xT , 1, 3, colorEpp7TeV, colorEpp7TeV, 1.8, kTRUE);
                graphXSecDirGammaTotArpp7TeV_xT->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphXSecDirGammaTotArpp7TeV_xT,0.00001);
            }
            TMarker* markerCurrent             = CreateMarkerFromGraph(graphXSecDirGammaStatpp7TeV_xT, 0.225-0.012 ,0.69,1);
            markerCurrent->SetNDC(kTRUE);
            markerCurrent->Draw("same,p");
//             markerCurrent->Print();
            TLatex *labelCurrent               = new TLatex(0.225,0.69,Form("%s (%s)", "ALICE prel.", "7 TeV" ));
            SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
            labelCurrent->Draw();
        }

        if (graphXSecDirGammaSyspp2760GeV_xT){
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaSyspp2760GeV_xT , markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV, widthLinesBoxes, kTRUE);
            graphXSecDirGammaSyspp2760GeV_xT->Draw("E2same");
        }
        if (graphXSecDirGammaStatpp2760GeV_xT){
            ProduceGraphAsymmWithoutXErrors(graphXSecDirGammaStatpp2760GeV_xT);
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaStatpp2760GeV_xT, markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV);
            graphXSecDirGammaStatpp2760GeV_xT->Draw("p,E1Z,same");
        }

        if (graphXSecDirGammaTotArpp2760GeV_xT){
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaTotArpp2760GeV_xT , 1, 3, colorEpp2760GeV, colorEpp2760GeV, 1.8, kTRUE);
            graphXSecDirGammaTotArpp2760GeV_xT->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphXSecDirGammaTotArpp2760GeV_xT,0.00001);
        }
        TMarker* markerCurrent2760GeV             = CreateMarkerFromGraph(graphXSecDirGammaStatpp2760GeV_xT, 0.55-0.012 ,0.7,1);
        markerCurrent2760GeV->SetNDC(kTRUE);
        markerCurrent2760GeV->Draw("same,p");
//         markerCurrent2760GeV->Print();
        TLatex *labelCurrent2760GeV               = new TLatex(0.55,0.7,Form("%s (%s)", "ALICE", "2.76 TeV" ));
        SetStyleTLatex( labelCurrent2760GeV, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
        labelCurrent2760GeV->Draw();

        if (graphXSecDirGammaSyspp8TeV_xT){
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaSyspp8TeV_xT , markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV, widthLinesBoxes, kTRUE);
            graphXSecDirGammaSyspp8TeV_xT->Draw("E2same");
        }
        if (graphXSecDirGammaStatpp8TeV_xT){
            ProduceGraphAsymmWithoutXErrors(graphXSecDirGammaStatpp8TeV_xT);
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaStatpp8TeV_xT, markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV);
            graphXSecDirGammaStatpp8TeV_xT->Draw("p,E1Z,same");
        }

        if (graphXSecDirGammaTotArpp8TeV_xT){
            DrawGammaSetMarkerTGraphAsym(graphXSecDirGammaTotArpp8TeV_xT , 1, 3, colorEpp8TeV, colorEpp8TeV, 1.8, kTRUE);
            graphXSecDirGammaTotArpp8TeV_xT->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphXSecDirGammaTotArpp8TeV_xT,0.00001);
        }
        TMarker* markerCurrent8TeV             = CreateMarkerFromGraph(graphXSecDirGammaStatpp8TeV_xT, 0.39-0.012 ,0.85,1);
        markerCurrent8TeV->SetNDC(kTRUE);
        markerCurrent8TeV->Draw("same,p");
//         markerCurrent8TeV->Print();
        TLatex *labelCurrent8TeV               = new TLatex(0.39,0.85,Form("%s (%s)", "ALICE", "8 TeV" ));
        SetStyleTLatex( labelCurrent8TeV, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
        labelCurrent8TeV->Draw();


        for (Int_t i = 0; i < nOtherEnergies; i++){
            if (xTOtherExperiments[i]){
                if (i == 2 & !enablePrelim)
                    continue;
                DrawGammaSetMarkerTGraphAsym(xTOtherExperiments[i] , markerStyleOtherExp[i], markerSizeOtherExp[i], colorOtherExp[i], colorOtherExp[i], widthLinesBoxes);
                xTOtherExperiments[i]->Draw("p,E1Z,same");

                cout << Form("%s (%s)", nameOtherExp[i].Data(), nameOtherEnergy[i].Data() ) << endl;
//                 xTOtherExperiments[i]->Print();
                TMarker* markerCurrent             = CreateMarkerFromGraph(xTOtherExperiments[i],legendOtherExtX[i]-0.012 ,legendOtherExtY[i],1);
                markerCurrent->SetNDC(kTRUE);
                markerCurrent->Draw("same,p");
//                 markerCurrent->Print();
                TLatex *labelCurrent               = new TLatex(legendOtherExtX[i],legendOtherExtY[i],Form("%s (%s)", nameOtherExp[i].Data(), nameOtherEnergy[i].Data() ));
                SetStyleTLatex( labelCurrent, textSizeLabelsPixelDirGam*0.5, 4, 1, 43, kTRUE, 12);
                labelCurrent->Draw();
            }
        }

    if (enablePrelim){
        canvasDirGamma2->Print(Form("%s/XTSpectraPP_withPrelim.%s",outputDir.Data(),suffix.Data()));
        canvasDirGamma2->Print(Form("%s/XTSpectraPP_withPrelim.pdf",outputDir.Data()));
    } else {
        canvasDirGamma2->Print(Form("%s/XTSpectraPP.%s",outputDir.Data(),suffix.Data()));
        canvasDirGamma2->Print(Form("%s/XTSpectraPP.pdf",outputDir.Data()));
    }

    dummyIncGamma->DrawCopy();

        if (graphInvYieldIncGammaSyspp2760GeV){
            graphInvYieldIncGammaSyspp2760GeV->Draw("E2same");
        }
        if (graphInvYieldIncGammaStatpp2760GeV){
            graphInvYieldIncGammaStatpp2760GeV->Draw("p,E1Z,same");
        }

        if (graphInvYieldIncGammaSyspp8TeV){
            graphInvYieldIncGammaSyspp8TeV->Draw("E2same");
        }
        if (graphInvYieldIncGammaStatpp8TeV){
            graphInvYieldIncGammaStatpp8TeV->Draw("p,E1Z,same");
        }
        if(fitInvYieldIncGammapp2760GeV){
          DrawGammaSetMarkerTF1( fitInvYieldIncGammapp2760GeV, 3, 2, colorEpp2760GeV);
          fitInvYieldIncGammapp2760GeV->SetRange(0.45,9.7);
          fitInvYieldIncGammapp2760GeV->Draw("same");
        }
        if(fitInvYieldIncGammapp8TeV){
          DrawGammaSetMarkerTF1( fitInvYieldIncGammapp8TeV, 3, 2, colorEpp8TeV);
          fitInvYieldIncGammapp8TeV->SetRange(0.3,18.);
          fitInvYieldIncGammapp8TeV->Draw("same");
        }
        TF1* fitInvYieldIncGammaDummy = (TF1*)fitInvYieldIncGammapp2760GeV->Clone("fitInvYieldIncGammaDummy");
        DrawGammaSetMarkerTF1( fitInvYieldIncGammaDummy, 3, 2, kGray+2);
        TLatex *labelEnergyIncInvYieldPaperAllwFit = new TLatex(0.23, 0.11+0.04*5, "#gamma_{inc} ALICE");
        SetStyleTLatex( labelEnergyIncInvYieldPaperAllwFit, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyIncInvYieldPaperAllwFit->Draw();
        TLegend* legendDirGammaPPonlywFits = GetAndSetLegend2(0.22, 0.11+(0*textsizeLabelsDirGamma*0.85), 0.21+0.25, 0.11+(5*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        legendDirGammaPPonlywFits->AddEntry(graphInvYieldDirGammaSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDirGammaPPonlywFits->AddEntry((TObject*)0,"norm. unc. 2.5%","");
        legendDirGammaPPonlywFits->AddEntry(graphInvYieldDirGammaSyspp8TeV,collisionSystempp8TeV.Data(),"pf");
        legendDirGammaPPonlywFits->AddEntry((TObject*)0,"norm. unc. 2.6%","");
        legendDirGammaPPonlywFits->AddEntry(fitInvYieldIncGammaDummy,"TCM fit","l");
        legendDirGammaPPonlywFits->Draw();


    canvasDirGamma->Print(Form("%s/IncGammaSpectraPP_wFit.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/IncGammaSpectraPP_wFit.pdf",outputDir.Data()));

    dummyDirGamma->GetYaxis()->SetRangeUser( 1.2e-10,8.5e4);
    dummyDirGamma->DrawCopy();

        TLatex *labelScalingDirGamma0020 = new TLatex(12.2,1.2E-3,nameLabelScalePbPb[0]);
        SetStyleTLatex( labelScalingDirGamma0020, 0.85*textsizeLabelsDirGamma,4,colorComb[0],42,kFALSE);
        labelScalingDirGamma0020->Draw();
        TLatex *labelScalingDirGamma2040 = new TLatex(12.2,6.0E-6,nameLabelScalePbPb[1]);
        SetStyleTLatex( labelScalingDirGamma2040, 0.85*textsizeLabelsDirGamma,4,colorComb[1],42,kFALSE);
        labelScalingDirGamma2040->Draw();
        TLatex *labelScalingDirGamma4080 = new TLatex(12.2,7.5E-9,nameLabelScalePbPb[2]);
        SetStyleTLatex( labelScalingDirGamma4080, 0.85*textsizeLabelsDirGamma,4,colorComb[2],42,kFALSE);
        labelScalingDirGamma4080->Draw();

        for (Int_t centPb = 0; centPb < 3; centPb++){
            if (graphInvYieldDirGammaSysPbPbPlot[centPb]){
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSysPbPbPlot[centPb], markerStyleComb[centPb], markerSizeComb[centPb], colorComb[centPb] , colorComb[centPb], widthLinesBoxes, kTRUE);
                graphInvYieldDirGammaSysPbPbPlot[centPb]->Draw("E2same");
            }
            if (graphInvYieldDirGammaStatPbPbPlot[centPb]){
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatPbPbPlot[centPb], markerStyleComb[centPb], markerSizeComb[centPb], colorComb[centPb] , colorComb[centPb]);
                graphInvYieldDirGammaStatPbPbPlot[centPb]->Draw("p,E1Z,same");
            }
            if (graphInvYieldDirGammaTotArPbPbPlot[centPb]){
                if (graphInvYieldDirGammaTotArPbPbPlot[centPb]->GetN() > 0){
                    DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArPbPbPlot[centPb] , 1, 3, colorComb[centPb], colorComb[centPb], 1.8, kTRUE);
                    graphInvYieldDirGammaTotArPbPbPlot[centPb]->Draw(">,same");
                    PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPbPlot[centPb]);
                }
            }

            if (graphInvYieldDirGammaTotArpp2760GeVScaled[centPb]){
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp2760GeVScaled[centPb] , 1, 3, colorCombAr[centPb], colorCombAr[centPb], 2.0, kTRUE);
                graphInvYieldDirGammaTotArpp2760GeVScaled[centPb]->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeVScaled[centPb]);
            }
            if (graphInvYieldDirGammaSyspp2760GeVScaled[centPb]){
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp2760GeVScaled[centPb] , markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV, colorCombAr[centPb], colorCombAr[centPb], widthLinesBoxes, kTRUE);
                graphInvYieldDirGammaSyspp2760GeVScaled[centPb]->Draw("E2,same");
            }
            if (graphInvYieldDirGammaStatpp2760GeVScaled[centPb]){
                DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp2760GeVScaled[centPb] , markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV, colorCombAr[centPb], colorCombAr[centPb]);
                graphInvYieldDirGammaStatpp2760GeVScaled[centPb]->Draw("p,E1Z,same");
            }
        }

        TLegend* legendDirGammaPP2 = GetAndSetLegend2(0.21, 0.10+(4*textsizeLabelsDirGamma*0.85), 0.21+0.21, 0.10+(8*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1,
                                                     collisionSystempp2760GeV.Data(), 42, 0.25);
        legendDirGammaPP2->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled[0],Form("ALICE x %1.1f",nColl[0]),"pf");
        legendDirGammaPP2->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled[1],Form("ALICE x %1.1f",nColl[1]),"pf");
        legendDirGammaPP2->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled[2],Form("ALICE x %1.1f",nColl[2]),"pf");
        legendDirGammaPP2->Draw();
        legendDirGamma->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectra_WithScaledPP.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectra_WithScaledPP.pdf",outputDir.Data()));

    dummyDirGamma->GetYaxis()->SetRangeUser( 1.2e-10,8.5e5);
    dummyDirGamma->DrawCopy();

        labelScalingDirGamma0020->Draw();
        labelScalingDirGamma2040->Draw();
        labelScalingDirGamma4080->Draw();

        for (Int_t centPb = 0; centPb < 3; centPb++){
            if (graphInvYieldDirGammaSysPbPbPlot[centPb]){
                graphInvYieldDirGammaSysPbPbPlot[centPb]->Draw("E2same");
            }
            if (graphInvYieldDirGammaStatPbPbPlot[centPb]){
                graphInvYieldDirGammaStatPbPbPlot[centPb]->Draw("p,E1Z,same");
            }
            if (graphInvYieldDirGammaTotArPbPbPlot[centPb]){
                if (graphInvYieldDirGammaTotArPbPbPlot[centPb]->GetN() > 0){
                    graphInvYieldDirGammaTotArPbPbPlot[centPb]->Draw(">,same");
                    PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPbPlot[centPb]);
                }
            }
            if (graphInvYieldDirGammaTotArpp2760GeVScaled[centPb]){
                graphInvYieldDirGammaTotArpp2760GeVScaled[centPb]->Draw(">,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeVScaled[centPb]);
            }
            if (graphInvYieldDirGammaSyspp2760GeVScaled[centPb]){
                graphInvYieldDirGammaSyspp2760GeVScaled[centPb]->Draw("E2,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaSyspp2760GeVScaled[centPb]);
            }
            if (graphInvYieldDirGammaStatpp2760GeVScaled[centPb]){
                graphInvYieldDirGammaStatpp2760GeVScaled[centPb]->Draw("p,E1Z,same");
                PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaStatpp2760GeVScaled[centPb]);
            }
        }

        TLatex *labelEnergyScaledPP = new TLatex(0.22, 0.095+0.04*21, "ALICE");
        SetStyleTLatex( labelEnergyScaledPP, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyScaledPP->Draw();
        TLatex *labelEnergyScaledPP2 = new TLatex(0.22, 0.10+0.04*3, "pp #times <#it{N}_{coll}> , #sqrt{#it{s}} = 2.76 TeV");
        SetStyleTLatex( labelEnergyScaledPP2, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyScaledPP2->Draw();
        TLegend* legendDirGammaPP2Split = GetAndSetLegend2(0.21, 0.10, 0.21+0.22, 0.10+(3*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        legendDirGammaPP2Split->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled[0],Form("#gamma_{dir} #times %1.1f",nColl[0]),"pf");
        legendDirGammaPP2Split->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled[1],Form("#gamma_{dir} #times %1.1f",nColl[1]),"pf");
        legendDirGammaPP2Split->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled[2],Form("#gamma_{dir} #times %1.1f",nColl[2]),"pf");
        legendDirGammaPP2Split->Draw();

        // TLegend* legendDirGammaSplit = GetAndSetLegend2(0.21, 0.10, 0.21+0.21, 0.10+(4*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, collisionSystemPbPb2760GeV.Data(), 42, 0.25);
        TLatex *labelEnergyScaledPBPB2 = new TLatex(0.95, 0.095+0.04*21, Form("%s",collisionSystemPbPb2760GeV.Data()));
        SetStyleTLatex( labelEnergyScaledPBPB2, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 31);
        labelEnergyScaledPBPB2->Draw();
        TLegend* legendDirGammaSplit = GetAndSetLegend2(0.71, 0.10+(21*textsizeLabelsDirGamma*0.85), 0.71+0.21, 0.10+(24*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        legendDirGammaSplit->AddEntry(graphInvYieldDirGammaSysPbPb[0],"  0-20% #gamma_{dir}","pf");
        legendDirGammaSplit->AddEntry(graphInvYieldDirGammaSysPbPb[1],"20-40% #gamma_{dir}","pf");
        legendDirGammaSplit->AddEntry(graphInvYieldDirGammaSysPbPb[2],"40-80% #gamma_{dir}","pf");
        legendDirGammaSplit->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectra_WithScaledPP_Split.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectra_WithScaledPP_Split.pdf",outputDir.Data()));


    //*******************************************************************************************************************************************
    //******************************************* DR plot with pp and PbPb measurements *********************************************************
    //*******************************************************************************************************************************************
    TGraphAsymmErrors* graphDRStatPbPbPlot[3];
    for (Int_t centPb = 0; centPb < 3; centPb++){
        graphDRStatPbPbPlot[centPb] = (TGraphAsymmErrors*)graphDRStatPbPb[centPb]->Clone(Form("graphDRStatPbPbPlot%s",namePbPbCentOut[centPb].Data()));
        ProduceGraphAsymmWithoutXErrors(graphDRStatPbPbPlot[centPb]);
    }
    TGraphAsymmErrors* graphDRStatpp2760GeVPlot = (TGraphAsymmErrors*)graphDRStatpp2760GeV->Clone("graphDRStatpp2760GeVPlot");
    ProduceGraphAsymmWithoutXErrors(graphDRStatpp2760GeVPlot);
    TGraphAsymmErrors* graphDRTotpp2760GeVPlot  = (TGraphAsymmErrors*)graphDRTotpp2760GeV->Clone("graphDRTotpp2760GeVPlot");
    ProduceGraphAsymmWithoutXErrors(graphDRTotpp2760GeVPlot);

    Double_t arrayBoundsXIndMeasRatio[2];
    Double_t arrayBoundsYIndMeasRatio[4];
    Double_t relativeMarginsIndMeasRatioX[3];
    Double_t relativeMarginsIndMeasRatioY[3];
    ReturnCorrectValuesForCanvasScaling(1200, 1400, 1, 3, 0.09, 0.01, 0.01, 0.065, arrayBoundsXIndMeasRatio, arrayBoundsYIndMeasRatio, relativeMarginsIndMeasRatioX, relativeMarginsIndMeasRatioY);

    TCanvas * canvasRatioIndDR = new TCanvas("canvasRatioIndDR","",10,10,1200,1400);  // gives the page size
    canvasRatioIndDR->cd();

    TPad* padPartRatioInDR1 = new TPad("padPartRatioInDR1", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
    DrawGammaPadSettings( padPartRatioInDR1, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
    padPartRatioInDR1->Draw();
    TPad* padPartRatioInDR2 = new TPad("padPartRatioInDR2", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
    DrawGammaPadSettings( padPartRatioInDR2, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[1]);
    padPartRatioInDR2->Draw();
    TPad* padPartRatioInDR3 = new TPad("padPartRatioInDR3", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[3], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2],-1, -1, -2);
    DrawGammaPadSettings( padPartRatioInDR3, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
    padPartRatioInDR3->Draw();

    //_______________________________________________________________ define text sizes _________________________________________________________
    Int_t textSizeLabelsPixel = 43;
    Double_t margin = relativeMarginsIndMeasRatioX[0]*1200;
    Double_t textsizeLabelsPad1 = 0;
    Double_t textsizeFacPad1 = 0;
    if (padPartRatioInDR1->XtoPixel(padPartRatioInDR1->GetX2()) < padPartRatioInDR1->YtoPixel(padPartRatioInDR1->GetY1())){
        textsizeLabelsPad1 = (Double_t)textSizeLabelsPixel/padPartRatioInDR1->XtoPixel(padPartRatioInDR1->GetX2()) ;
        textsizeFacPad1 = (Double_t)1./padPartRatioInDR1->XtoPixel(padPartRatioInDR1->GetX2()) ;
    } else {
        textsizeLabelsPad1 = (Double_t)textSizeLabelsPixel/padPartRatioInDR1->YtoPixel(padPartRatioInDR1->GetY1());
        textsizeFacPad1 = (Double_t)1./padPartRatioInDR1->YtoPixel(padPartRatioInDR1->GetY1());
    }
    Double_t textsizeLabelsPad2 = 0;
    Double_t textsizeFacPad2 = 0;
    if (padPartRatioInDR2->XtoPixel(padPartRatioInDR2->GetX2()) <padPartRatioInDR2->YtoPixel(padPartRatioInDR2->GetY1()) ){
        textsizeLabelsPad2 = (Double_t)textSizeLabelsPixel/padPartRatioInDR2->XtoPixel(padPartRatioInDR2->GetX2()) ;
        textsizeFacPad2 = (Double_t)1./padPartRatioInDR2->XtoPixel(padPartRatioInDR2->GetX2()) ;
    } else {
        textsizeLabelsPad2 = (Double_t)textSizeLabelsPixel/padPartRatioInDR2->YtoPixel(padPartRatioInDR2->GetY1());
        textsizeFacPad2 = (Double_t)1./padPartRatioInDR2->YtoPixel(padPartRatioInDR2->GetY1());
    }
    Double_t textsizeLabelsPad3 = 0;
    Double_t textsizeFacPad3= 0;
    if (padPartRatioInDR3->XtoPixel(padPartRatioInDR3->GetX2()) <padPartRatioInDR3->YtoPixel(padPartRatioInDR3->GetY1()) ){
        textsizeLabelsPad3 = (Double_t)textSizeLabelsPixel/padPartRatioInDR3->XtoPixel(padPartRatioInDR3->GetX2()) ;
        textsizeFacPad3 = (Double_t)1./padPartRatioInDR3->XtoPixel(padPartRatioInDR3->GetX2()) ;
    } else {
        textsizeLabelsPad3 = (Double_t)textSizeLabelsPixel/padPartRatioInDR3->YtoPixel(padPartRatioInDR3->GetY1());
        textsizeFacPad3 = (Double_t)1./padPartRatioInDR3->YtoPixel(padPartRatioInDR3->GetY1());
    }

    //_______________________________________________________________ 0-20% dummy upper panel ___________________________________________________
    TH2D *dummyDR1 ;
    dummyDR1 = new TH2D("dummyDR1", "dummyDR1", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR1, "#it{p}_{T} (GeV/#it{c})", "",
                            0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.95,0.10/(textsizeFacPad1*margin), 510, 505);
    dummyDR1->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR1->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR1->GetXaxis()->SetTickLength(0.06);
    dummyDR1->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 20-40% dummy middle panel _________________________________________________
    TH2D *dummyDR2 ;
    dummyDR2 = new TH2D("dummyDR2", "dummyDR2", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR2, "#it{p}_{T} (GeV/#it{c})", "#it{R}_{#gamma}", // = (#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})
                            0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.95,0.10/(textsizeFacPad2*margin), 510, 505);
    dummyDR2->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR2->GetYaxis()->CenterTitle(kTRUE);
    dummyDR2->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR2->GetXaxis()->SetTickLength(0.06);
    dummyDR2->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 40-80% dummy lower panel __________________________________________________
    TH2D *dummyDR3 ;
    dummyDR3 = new TH2D("dummyDR3", "dummyDR3", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR3, "#it{p}_{T} (GeV/#it{c})", "",
                            0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.92,0.10/(textsizeFacPad3*margin), 510, 505);
    dummyDR3->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR3->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR3->GetXaxis()->SetTickLength(0.055);
    dummyDR3->GetYaxis()->SetTickLength(0.035);

    //_______________________________________________________________ 0-20% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphDRSysPbPb[0], markerStyleComb[0], markerSizeComb[0], colorComb[0] , colorComb[0]);
        graphDRSysPbPb[0]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRSyspp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRSyspp2760GeV->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatPbPbPlot[0], markerStyleComb[0], markerSizeComb[0], colorComb[0] , colorComb[0]);
        graphDRStatPbPbPlot[0]->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp2760GeVPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");

        TLegend* legendDR0020 = GetAndSetLegend2(0.12, 0.9-(3*textsizeLabelsPad1), 0.12+0.21, 0.9,textsizeLabelsPad1, 1,
                                                 textALICE, 42, 0.3);
        legendDR0020->AddEntry(graphDRSysPbPb[0],Form("%s %s", "0-20%", collisionSystemPbPb2760GeV.Data()),"pf");
        legendDR0020->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDR0020->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphDRSysPbPb[1], markerStyleComb[1], markerSizeComb[1], colorComb[1] , colorComb[1]);
        graphDRSysPbPb[1]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRSyspp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRSyspp2760GeV->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatPbPbPlot[1], markerStyleComb[1], markerSizeComb[1], colorComb[1] , colorComb[1]);
        graphDRStatPbPbPlot[1]->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp2760GeVPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");

        TLegend* legendDR2040 = GetAndSetLegend2(0.12, 0.94-(2*textsizeLabelsPad2), 0.12+0.21, 0.94,textsizeLabelsPad2, 1,
                                                 "", 42, 0.3);
        legendDR2040->AddEntry(graphDRSysPbPb[1],Form("%s %s", "20-40%", collisionSystemPbPb2760GeV.Data()),"pf");
        legendDR2040->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDR2040->Draw();

    //_______________________________________________________________ 40-80% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);

        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphDRSysPbPb[2], markerStyleComb[2], markerSizeComb[2], colorComb[2] , colorComb[2]);
        graphDRSysPbPb[2]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRSyspp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRSyspp2760GeV->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatPbPbPlot[2], markerStyleComb[2], markerSizeComb[2], colorComb[2] , colorComb[2]);
        graphDRStatPbPbPlot[2]->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp2760GeVPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");

        TLegend* legendDR4080 = GetAndSetLegend2(0.12, 0.94-(2*textsizeLabelsPad3),0.12+0.21, 0.94,textsizeLabelsPad3, 1,
                                                      "", 42, 0.3);
        legendDR4080->AddEntry(graphDRSysPbPb[2],Form("%s %s","40-80%",collisionSystemPbPb2760GeV.Data()),"pf");
        legendDR4080->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDR4080->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_PbPbWithRefPP.%s", outputDir.Data(), suffix.Data()));
    canvasRatioIndDR->SaveAs(Form("%s/DR_PbPbWithRefPP.pdf", outputDir.Data()));



    TGraphAsymmErrors* graphDRStatpp8TeVPlot = (TGraphAsymmErrors*)graphDRStatpp8TeV->Clone("graphDRStatpp8TeVPlot");
    ProduceGraphAsymmWithoutXErrors(graphDRStatpp8TeVPlot);
    TGraphAsymmErrors* graphDRTotpp8TeVPlot = (TGraphAsymmErrors*)graphDRTotpp8TeV->Clone("graphDRTotpp8TeVPlot");
    ProduceGraphAsymmWithoutXErrors(graphDRTotpp8TeVPlot);
    TGraphAsymmErrors* graphDRStatpp7TeVPlot = (TGraphAsymmErrors*)graphDRStatpp7TeV->Clone("graphDRStatpp7TeVPlot");
    ProduceGraphAsymmWithoutXErrors(graphDRStatpp7TeVPlot);
    TGraphAsymmErrors* graphDRTotpp7TeVPlot = (TGraphAsymmErrors*)graphDRTotpp7TeV->Clone("graphDRTotpp7TeVPlot");
    ProduceGraphAsymmWithoutXErrors(graphDRTotpp7TeVPlot);
    TGraphAsymmErrors* graphRGammaLMEE7TeVStatPlot = (TGraphAsymmErrors*)graphRGammaLMEE7TeVStat->Clone("graphRGammaLMEE7TeVStat");
    ProduceGraphAsymmWithoutXErrors(graphRGammaLMEE7TeVStatPlot);
    TGraphAsymmErrors* graphRGammaLMEE7TeVTotPlot = (TGraphAsymmErrors*)graphRGammaLMEE7TeVTot->Clone("graphRGammaLMEE7TeVTot");
    ProduceGraphAsymmWithoutXErrors(graphRGammaLMEE7TeVTotPlot);


    // double ratio combined
    TCanvas *canvasDoubleRatio = new TCanvas("canvasDoubleRatio","",0.095,0.09,1000,815);
    DrawGammaCanvasSettings( canvasDoubleRatio, 0.086, 0.01, 0.01, 0.105);
    canvasDoubleRatio->cd();
    canvasDoubleRatio->SetLogx();

    widthLinesBoxes                         = 1.4;
    Double_t textSizeSinglePad               = 0.05;
    TH2F * hist2DDRDummySingle       = new TH2F("hist2DDRDummySingle","hist2DDRDummySingle",1000,doubleRatioX[0], doubleRatioX[1],1000,0.72, 1.55);
    SetStyleHistoTH2ForGraphs(hist2DDRDummySingle, "#it{p}_{T} (GeV/#it{c})","#it{R}_{#gamma}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    hist2DDRDummySingle->GetXaxis()->SetNoExponent();
    hist2DDRDummySingle->GetXaxis()->SetMoreLogLabels(kTRUE);
    hist2DDRDummySingle->DrawCopy();

        Int_t nRowsDRLegend = 3;
        if (enablePrelim) nRowsDRLegend++;

        TLegend* legendDRSingle = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*(nRowsDRLegend+0.3),0.43,0.953, textSizeSinglePad, 1, "ALICE", 42, 0.3);
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        if (enablePrelim){
            DrawGammaSetMarkerTGraphAsym(graphDRSyspp7TeV, markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV , colorEpp7TeV,widthLinesBoxes, kTRUE);
            graphDRSyspp7TeV->Draw("E2same");
        }

        DrawGammaSetMarkerTGraphAsym(graphDRSyspp2760GeV, markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV , colorEpp2760GeV,widthLinesBoxes, kTRUE);
        graphDRSyspp2760GeV->Draw("E2same");
        legendDRSingle->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        if (enablePrelim){
            legendDRSingle->AddEntry(graphDRSyspp7TeV,Form("%s (prelim.)", collisionSystempp7TeV.Data()),"pf");
        }
        DrawGammaSetMarkerTGraphAsym(graphDRSyspp8TeV, markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV , colorEpp8TeV,widthLinesBoxes, kTRUE);
        graphDRSyspp8TeV->Draw("E2same");
        legendDRSingle->AddEntry(graphDRSyspp8TeV,collisionSystempp8TeV.Data(),"pf");

        if (enablePrelim){
            DrawGammaSetMarkerTGraphAsym(graphDRStatpp7TeVPlot,  markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV , colorEpp7TeV);
            graphDRStatpp7TeVPlot->Draw("pp,E1Z,same");
        }
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp2760GeVPlot,  markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV , colorEpp2760GeV);
        graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp8TeVPlot,  markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV , colorEpp8TeV);
        graphDRStatpp8TeVPlot->Draw("pp,E1Z,same");

        legendDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    if (enablePrelim){
        canvasDoubleRatio->Print(Form("%s/DR_PP_withPrelim.%s", outputDir.Data(), suffix.Data()));
        canvasDoubleRatio->Print(Form("%s/DR_PP_withPrelim.pdf", outputDir.Data()));
    } else {
        canvasDoubleRatio->Print(Form("%s/DR_PP.%s", outputDir.Data(), suffix.Data()));
        canvasDoubleRatio->Print(Form("%s/DR_PP.pdf", outputDir.Data()));
    }
    hist2DDRDummySingle->DrawCopy();

    TLegend* legendDRSingle3 = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*(nRowsDRLegend+0.3),0.43,0.953, textSizeSinglePad, 1, "ALICE", 42, 0.15);
    TLatex *labelDisplace = new TLatex(0.95,0.14,"data points displaced in #it{p}_{T} for visualization");
    SetStyleTLatex( labelDisplace, textSizeLabelsPixel*0.65, 4, 1, 43, kTRUE, 31);
    labelDisplace->Draw();
    TLatex *labelUnc = new TLatex(0.95,0.17,"error bars represent total uncertainties");
    SetStyleTLatex( labelUnc, textSizeLabelsPixel*0.65, 4, 1, 43, kTRUE, 31);
    labelUnc->Draw();

        if (enablePrelim){
            ProduceGraphAsymmDisplacedInX(graphDRTotpp7TeVPlot, -0.015);
            DrawGammaSetMarkerTGraphAsym(graphDRTotpp7TeVPlot,  markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV , colorEpp7TeV);
            graphDRTotpp7TeVPlot->Draw("pp,E1,same");
        }
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphDRTotpp2760GeVPlot,  markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV , colorEpp2760GeV);
        graphDRTotpp2760GeVPlot->Draw("p,E1,same");
        legendDRSingle3->AddEntry(graphDRTotpp2760GeVPlot,collisionSystempp2760GeV.Data(),"p");
        if (enablePrelim){
            legendDRSingle3->AddEntry(graphDRTotpp7TeVPlot,Form("%s (prelim.)", collisionSystempp7TeV.Data()),"p");
        }
        ProduceGraphAsymmDisplacedInX(graphDRTotpp8TeVPlot, +0.015);
        DrawGammaSetMarkerTGraphAsym(graphDRTotpp8TeVPlot,  markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV , colorEpp8TeV);
        graphDRTotpp8TeVPlot->Draw("pp,E1,same");
        legendDRSingle3->AddEntry(graphDRTotpp8TeVPlot,collisionSystempp8TeV.Data(),"p");

        legendDRSingle3->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    if (enablePrelim){
        canvasDoubleRatio->Print(Form("%s/DR_PP_tot_withPrelim.%s", outputDir.Data(), suffix.Data()));
        canvasDoubleRatio->Print(Form("%s/DR_PP_tot_withPrelim.pdf", outputDir.Data()));
    } else {
        canvasDoubleRatio->Print(Form("%s/DR_PP_tot.%s", outputDir.Data(), suffix.Data()));
        canvasDoubleRatio->Print(Form("%s/DR_PP_tot.pdf", outputDir.Data()));
    }

       hist2DDRDummySingle->DrawCopy();

        nRowsDRLegend = 4;
        if (enablePrelim) nRowsDRLegend++;

        TLegend* legendDRSingle5 = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*(nRowsDRLegend+0.3),0.43,0.953, textSizeSinglePad, 1, "ALICE", 42, 0.3);
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        if (enablePrelim){
            DrawGammaSetMarkerTGraphAsym(graphDRSyspp7TeV, markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV , colorEpp7TeV,widthLinesBoxes, kTRUE);
            graphDRSyspp7TeV->Draw("E2same");
        }

        DrawGammaSetMarkerTGraphAsym(graphDRSyspp2760GeV, markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV , colorEpp2760GeV,widthLinesBoxes, kTRUE);
        graphDRSyspp2760GeV->Draw("E2same");
        legendDRSingle5->AddEntry(graphDRSyspp2760GeV,Form("#gamma, %s", collisionSystempp2760GeV.Data()),"pf");
        if (enablePrelim){
            legendDRSingle5->AddEntry(graphDRSyspp7TeV,Form("#gamma, %s (prelim.)", collisionSystempp7TeV.Data()),"pf");
        }

        DrawGammaSetMarkerTGraphAsym(graphRGammaLMEE7TeVSys, markerStyleEpp7TeV+4, markerSizeEpp7TeV, colorEpp7TeV+1 , colorEpp7TeV+1, widthLinesBoxes, kTRUE);
        graphRGammaLMEE7TeVSys->Draw("E2same");
        legendDRSingle5->AddEntry(graphRGammaLMEE7TeVSys,Form("#gamma*, %s", collisionSystempp7TeV.Data()),"pf");

        DrawGammaSetMarkerTGraphAsym(graphDRSyspp8TeV, markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV , colorEpp8TeV,widthLinesBoxes, kTRUE);
        graphDRSyspp8TeV->Draw("E2same");
        legendDRSingle5->AddEntry(graphDRSyspp8TeV,Form("#gamma, %s", collisionSystempp8TeV.Data()),"pf");

        if (enablePrelim){
            DrawGammaSetMarkerTGraphAsym(graphDRStatpp7TeVPlot,  markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV , colorEpp7TeV);
            graphDRStatpp7TeVPlot->Draw("p,E1Z,same");
        }
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp2760GeVPlot,  markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV , colorEpp2760GeV);
        graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphRGammaLMEE7TeVStatPlot,  markerStyleEpp7TeV+4, markerSizeEpp7TeV, colorEpp7TeV+1 , colorEpp7TeV+1);
        graphRGammaLMEE7TeVStatPlot->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp8TeVPlot,  markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV , colorEpp8TeV);
        graphDRStatpp8TeVPlot->Draw("p,E1Z,same");
        graphRGammaLMEE7TeVStatPlot->Draw("p,E1Z,same");
        legendDRSingle5->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    if (enablePrelim){
        canvasDoubleRatio->Print(Form("%s/DRwLMEE_PP_withPrelim.%s", outputDir.Data(), suffix.Data()));
        canvasDoubleRatio->Print(Form("%s/DRwLMEE_PP_withPrelim.pdf", outputDir.Data()));
    } else {
        canvasDoubleRatio->Print(Form("%s/DRwLMEE_PP.%s", outputDir.Data(), suffix.Data()));
        canvasDoubleRatio->Print(Form("%s/DRwLMEE_PP.pdf", outputDir.Data()));
    }
    hist2DDRDummySingle->DrawCopy();

    TLegend* legendDRSingle6 = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*(nRowsDRLegend+0.3),0.43,0.953, textSizeSinglePad, 1, "ALICE", 42, 0.15);
        labelDisplace->Draw();
        labelUnc->Draw();
        if (enablePrelim){
            DrawGammaSetMarkerTGraphAsym(graphDRTotpp7TeVPlot,  markerStyleEpp7TeV, markerSizeEpp7TeV, colorEpp7TeV , colorEpp7TeV);
            graphDRTotpp7TeVPlot->Draw("pp,E1,same");
        }
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphDRTotpp2760GeVPlot,  markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV , colorEpp2760GeV);
        graphDRTotpp2760GeVPlot->Draw("p,E1,same");
        legendDRSingle6->AddEntry(graphDRTotpp2760GeVPlot,Form("#gamma, %s", collisionSystempp2760GeV.Data()),"p");
        if (enablePrelim){
            legendDRSingle6->AddEntry(graphDRTotpp7TeVPlot,Form("#gamma, %s (prelim.)", collisionSystempp7TeV.Data()),"p");
        }
        DrawGammaSetMarkerTGraphAsym(graphRGammaLMEE7TeVTotPlot,  markerStyleEpp7TeV+4, markerSizeEpp7TeV, colorEpp7TeV+1 , colorEpp7TeV+1);
        graphRGammaLMEE7TeVTotPlot->Draw("pp,E1Z,same");
        legendDRSingle6->AddEntry(graphRGammaLMEE7TeVTotPlot, Form("#gamma*, %s", collisionSystempp7TeV.Data()),"p");

        DrawGammaSetMarkerTGraphAsym(graphDRTotpp8TeVPlot,  markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV , colorEpp8TeV);
        graphDRTotpp8TeVPlot->Draw("pp,E1,same");
        legendDRSingle6->AddEntry(graphDRTotpp8TeVPlot,Form("#gamma, %s", collisionSystempp8TeV.Data()),"p");
        graphRGammaLMEE7TeVTotPlot->Draw("pp,E1Z,same");
        legendDRSingle6->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    if (enablePrelim){
        canvasDoubleRatio->Print(Form("%s/DRwLMEE_PP_tot_withPrelim.%s", outputDir.Data(), suffix.Data()));
        canvasDoubleRatio->Print(Form("%s/DRwLMEE_PP_tot_withPrelim.pdf", outputDir.Data()));
    } else {
        canvasDoubleRatio->Print(Form("%s/DRwLMEE_PP_tot.%s", outputDir.Data(), suffix.Data()));
        canvasDoubleRatio->Print(Form("%s/DRwLMEE_PP_tot.pdf", outputDir.Data()));
    }

    TFile* fileTheory                               = new TFile("ExternalInput/Theory/TheoryCompilationPP.root");
    TGraph* graphTheoryNLODRpp7TeVPromptThermalLiuWerner  = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaThermalAndPromptDirectPhotonLiuWerner_pp7TeV_ALICECocktail");
    while(graphTheoryNLODRpp7TeVPromptThermalLiuWerner->GetX()[graphTheoryNLODRpp7TeVPromptThermalLiuWerner->GetN()-1] > 20.) graphTheoryNLODRpp7TeVPromptThermalLiuWerner->RemovePoint(graphTheoryNLODRpp7TeVPromptThermalLiuWerner->GetN()-1);

    hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRSingle2 = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*5.3,0.43,0.953, textSizeSinglePad, 1, "ALICE", 42, 0.3);
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        graphDRSyspp2760GeV->Draw("E2same");
        legendDRSingle2->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");

        graphDRSyspp8TeV->Draw("E2same");
        legendDRSingle2->AddEntry(graphDRSyspp8TeV,collisionSystempp8TeV.Data(),"pf");

        graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");
        graphDRStatpp8TeVPlot->Draw("pp,E1Z,same");

        if (graphTheoryNLODRpp7TeVPromptThermalLiuWerner){
            DrawGammaNLOTGraph( graphTheoryNLODRpp7TeVPromptThermalLiuWerner, 2, 5, kPink+2 );
            graphTheoryNLODRpp7TeVPromptThermalLiuWerner->Draw("lc,same");
            legendDRSingle2->AddEntry(graphTheoryNLODRpp7TeVPromptThermalLiuWerner,"NLO pQCD, #scale[0.75]{Prompt + Thermal}","l");
            legendDRSingle2->AddEntry((TObject*)0,"for pp, #sqrt{#it{s}} = 7 TeV","");
        }

        legendDRSingle2->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_PP_wThermal.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_PP_wThermal.pdf", outputDir.Data()));
}
