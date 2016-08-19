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
#include "CombineMesonMeasurementsPbPbLHC11hV1.h"

extern TRandom* gRandom;
extern TBenchmark*  gBenchmark;
extern TSystem* gSystem;
extern TMinuit*     gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //  TString name;
};

void CombineMesonMeasurementsPbPbAllCentAndMeas(TString meson = "Eta",
                                                TString fileNamePCM = "",
                                                TString suffix = "pdf",
                                                Bool_t noXerrorBars = kTRUE,
                                                TString bWCorrection="X"){

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput           = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;

    Int_t nCent = 8;
    Int_t method = 3;
    TString meth[method];
    meth[0] = "PCM";
    meth[1] = "PHOS";
    meth[2] = "EMCal";
    TString cent[nCent];
    cent[0]     = "0005";
    cent[1]     = "0510";
    cent[2]     = "0010";
    cent[3]     = "1020";
    cent[4]     = "2040";
    cent[5]     = "2050";
    cent[6]     = "4060";
    cent[7]     = "6080";

    TString collisionSystPbPb[nCent];
    collisionSystPbPb[0]     = "0#font[122]{-}5% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
    collisionSystPbPb[1]     = "5#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
    collisionSystPbPb[2]     = "0#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
    collisionSystPbPb[3]     = "10#font[122]{-}20% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
    collisionSystPbPb[4]     = "20#font[122]{-}40% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
    collisionSystPbPb[5]     = "20#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
    collisionSystPbPb[6]     = "40#font[122]{-}60% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
    collisionSystPbPb[7]     = "60#font[122]{-}80% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";

    Color_t colorCent[nCent];
    colorCent[0] = kRed+1;
    colorCent[1] = kOrange+8;
    colorCent[2] = kOrange+1;
    colorCent[3] = kOrange+1;
    colorCent[4] = kGreen+2;
    colorCent[5] = kAzure+1;
    colorCent[6] = kMagenta+2;
    colorCent[7] = kViolet;

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

    //___________________________________ Declaration of files _____________________________________________
    TString fileNamePHOS                    = "ExternalInputPbPb/PHOS/LHC10h_PHOS_pi0_PbPb_06022014.root";
    TString fileNameEMCal                   = "ExternalInputPbPb/EMCAL/LHC11h_EMCal_pi0eta_PbPb_10092015.root";
    TString fileNamePCMPub                  = "LHC11hInputFiles/data_PCMResults_PbPb_2.76TeV_LHC10h.root";

    TString outputDir                       = Form("%s/%s/NeutralMesonMeasurementsPbPb2760GeVAllCent%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
    TString rootFiles                       = Form("%s/rootFiles",outputDir.Data());
    TString nameFinalResDat                 = Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+rootFiles);
    gSystem->Exec(Form("cp %s %s/InputPCM.root",        fileNamePCM.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root",       fileNamePHOS.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCal.root",      fileNameEMCal.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMPub.root",     fileNamePCMPub.Data(), rootFiles.Data()));

    TH1D *histoarrayPi0InvYieldPbPb2760GeV[nCent][method];
    TH1D* histoarrayPi0InvYieldPbPb2760GeVYshifted[nCent][method];
    TH1D *histoarrayEtaInvYieldPbPb2760GeV[nCent][method];
    TH1D* histoarrayEtaInvYieldPbPb2760GeVYshifted[nCent][method];
    TH1D* histoarrayEtatoPi0Stat2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayPi0InvYieldStatPbPb2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayPi0InvYieldSysPbPb2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayEtaInvYieldStatPbPb2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayEtaInvYieldSysPbPb2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayEtaToPi0Stat2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayEtaToPi0Sys2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayInvYieldTotPbPb2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayInvYieldStatPbPb2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayInvYieldSysPbPb2760GeV[nCent][method];
    for(Int_t c=0; c<nCent; c++){
      for(Int_t m=0; m<3; m++){
        histoarrayPi0InvYieldPbPb2760GeV[c][m] = NULL;
        grapharrayPi0InvYieldStatPbPb2760GeV[c][m] = NULL;
        grapharrayPi0InvYieldSysPbPb2760GeV[c][m] = NULL;

        histoarrayEtaInvYieldPbPb2760GeV[c][m] = NULL;
        grapharrayEtaInvYieldStatPbPb2760GeV[c][m] = NULL;
        grapharrayEtaInvYieldSysPbPb2760GeV[c][m] = NULL;
        histoarrayEtatoPi0Stat2760GeV[c][m] = NULL;
        grapharrayEtaToPi0Stat2760GeV[c][m] = NULL;
        grapharrayEtaToPi0Sys2760GeV[c][m] = NULL;

        grapharrayInvYieldTotPbPb2760GeV[c][m] = NULL;
        grapharrayInvYieldStatPbPb2760GeV[c][m] = NULL;
        grapharrayInvYieldSysPbPb2760GeV[c][m] = NULL;
      }
    }

    TFile* filePCM      = new TFile(fileNamePCM.Data());
    TFile* filePHOS     = new TFile(fileNamePHOS);
    TFile* fileEMCal    = new TFile(fileNameEMCal.Data());
    TFile* filePCMPub   = new TFile(fileNamePCMPub.Data());

    cout << "For the Pi0 in 0-5% " << endl;
    cout << "PCM" << endl;
    TDirectoryFile* directoryPCMPi0PbPb2760GeV_0005         = (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_0-5%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_0005                   = (TH1D*)directoryPCMPi0PbPb2760GeV_0005->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0005       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0005->Get("Pi0SystError");
    cout << "PHOS" << endl;
    TDirectory* directoryPHOSPi0PbPb0005 = (TDirectory*)filePHOS->Get("pi0_PbPb_2760_Centrality_0-5%");
    TH1D*histoPi0PHOSPbPb0005 =      (TH1D*)directoryPHOSPi0PbPb0005->Get("hPi0_PbPb_cen0_NoBW_Stat");
    TH1D*histoPi0PHOSSysPbPb0005 =   (TH1D*)directoryPHOSPi0PbPb0005->Get("hPi0_PbPb_cen0_NoBW_Syst");
    TGraphAsymmErrors*graphPHOSYieldPi0SysErrPbPb0005 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb0005);

    histoarrayPi0InvYieldPbPb2760GeV[0][0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_0005->Clone("histoPCMPi0InvYieldPbPb2760GeV_0005");
    grapharrayPi0InvYieldStatPbPb2760GeV[0][0] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[0][0]);
    grapharrayPi0InvYieldSysPbPb2760GeV[0][0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0005->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_0005");

    histoarrayPi0InvYieldPbPb2760GeV[0][1] = (TH1D*)histoPi0PHOSPbPb0005->Clone("histoPHOSPi0InvYieldPbPb2760GeV_0005");
    grapharrayPi0InvYieldStatPbPb2760GeV[0][1] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[0][1]);
    grapharrayPi0InvYieldSysPbPb2760GeV[0][1] = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb0005->Clone("graphPHOSPi0InvYieldSysPbPb2760GeV_0005");


    cout << "For the Pi0 in 5-10% " << endl;
    cout << "PCM" << endl;
    TDirectoryFile* directoryPCMPi0PbPb2760GeV_0510         = (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_5-10%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_0510                   = (TH1D*)directoryPCMPi0PbPb2760GeV_0510->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0510       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0510->Get("Pi0SystError");
    cout << "PHOS" << endl;
    TDirectory* directoryPHOSPi0PbPb0510 = (TDirectory*)filePHOS->Get("pi0_PbPb_2760_Centrality_5-10%");
    TH1D*histoPi0PHOSPbPb0510 =      (TH1D*)directoryPHOSPi0PbPb0510->Get("hPi0_PbPb_cen1_NoBW_Stat");
    TH1D*histoPi0PHOSSysPbPb0510 =   (TH1D*)directoryPHOSPi0PbPb0510->Get("hPi0_PbPb_cen1_NoBW_Syst");
    TGraphAsymmErrors*graphPHOSYieldPi0SysErrPbPb0510 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb0510);

    histoarrayPi0InvYieldPbPb2760GeV[1][0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_0510->Clone("histoPCMPi0InvYieldPbPb2760GeV_0510");
    grapharrayPi0InvYieldStatPbPb2760GeV[1][0] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[1][0]);
    grapharrayPi0InvYieldSysPbPb2760GeV[1][0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0510->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_0510");

    histoarrayPi0InvYieldPbPb2760GeV[1][1] = (TH1D*)histoPi0PHOSPbPb0510->Clone("histoPHOSPi0InvYieldPbPb2760GeV_0510");
    grapharrayPi0InvYieldStatPbPb2760GeV[1][1] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[1][1]);
    grapharrayPi0InvYieldSysPbPb2760GeV[1][1] = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb0510->Clone("graphPHOSPi0InvYieldSysPbPb2760GeV_0510");


    cout << "For the Pi0 in 0-10% " << endl;
    cout << "PCM" << endl;
    TDirectoryFile* directoryPCMPi0PbPb2760GeV_0010         = (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_0-10%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_0010                   = (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("CorrectedYieldPi0");
    TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_0010                   = (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("CorrectedYieldBinShiftedPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0010       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0SystError");
    cout << "PHOS" << endl;
    TDirectory* directoryPHOSPi0PbPb0010 = (TDirectory*)filePHOS->Get("pi0_PbPb_2760_Centrality_0-10%");
    TH1D*histoPi0PHOSPbPb0010 =      (TH1D*)directoryPHOSPi0PbPb0010->Get("hPi0_PbPb_cen7_NoBW_Stat");
    TH1D*histoPi0PHOSSysPbPb0010 =   (TH1D*)directoryPHOSPi0PbPb0010->Get("hPi0_PbPb_cen7_NoBW_Syst");
    TGraphAsymmErrors*graphPHOSYieldPi0SysErrPbPb0010 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb0010);
    cout << "EMCal" << endl;
    TDirectory* directoryEMCalPi0PbPb2760GeV                = (TDirectory*)fileEMCal->Get("Pi02.76TeV_PbPb");
    TH1D* histoEMCalPi0InvYieldPbPb2760GeV_0010         = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("LuciaBinningInvYieldPbPbStatErrPi_0010");
    TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_0010        = new TGraphAsymmErrors(histoEMCalPi0InvYieldPbPb2760GeV_0010);
    graphEMCalPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_0010->GetN()-1);
    graphEMCalPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_0010->GetN()-1);
    TH1D* histoEMCalPi0InvYieldSysPbPb2760GeV_0010  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("LuciaBinningInvYieldPbPbSysErrPi_0010");
    TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_0010     = new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_0010);
    graphEMCalPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);
    graphEMCalPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);

    histoarrayPi0InvYieldPbPb2760GeV[2][0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_0010->Clone("histoPCMPi0InvYieldPbPb2760GeV_0010");
    grapharrayPi0InvYieldStatPbPb2760GeV[2][0] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[2][0]);
    grapharrayPi0InvYieldSysPbPb2760GeV[2][0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0010->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_0010");

    histoarrayPi0InvYieldPbPb2760GeV[2][1] = (TH1D*)histoPi0PHOSPbPb0010->Clone("histoPHOSPi0InvYieldPbPb2760GeV_0010");
    grapharrayPi0InvYieldStatPbPb2760GeV[2][1] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[2][1]);
    grapharrayPi0InvYieldSysPbPb2760GeV[2][1] = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb0010->Clone("graphPHOSPi0InvYieldSysPbPb2760GeV_0010");

    histoarrayPi0InvYieldPbPb2760GeV[2][2] = (TH1D*)histoEMCalPi0InvYieldPbPb2760GeV_0010->Clone("histoEMCalPi0InvYieldPbPb2760GeV_0010");
    grapharrayPi0InvYieldStatPbPb2760GeV[2][2] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[2][2]);
    grapharrayPi0InvYieldStatPbPb2760GeV[2][2]->RemovePoint(grapharrayPi0InvYieldStatPbPb2760GeV[2][2]->GetN()-1);
    grapharrayPi0InvYieldStatPbPb2760GeV[2][2]->RemovePoint(grapharrayPi0InvYieldStatPbPb2760GeV[2][2]->GetN()-1);
    grapharrayPi0InvYieldSysPbPb2760GeV[2][2] = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Clone("graphEMCalPi0InvYieldSysPbPb2760GeV_0010");


    cout << "For the Pi0 in 10-20% " << endl;
    cout << "PCM" << endl;
    TDirectoryFile* directoryPCMPi0PbPb2760GeV_1020         = (TDirectoryFile*)filePCMPub->Get("Pi0_PbPb_2.76TeV_10-20%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_1020                   = (TH1D*)directoryPCMPi0PbPb2760GeV_1020->Get("CorrectedYieldPi0");
    TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_1020                   = (TH1D*)directoryPCMPi0PbPb2760GeV_1020->Get("CorrectedYieldBinShiftedPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_1020       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_1020->Get("Pi0SystError");
    cout << "PHOS" << endl;
    TDirectory* directoryPHOSPi0PbPb1020 = (TDirectory*)filePHOS->Get("pi0_PbPb_2760_Centrality_10-20%");
    TH1D*histoPi0PHOSPbPb1020 =      (TH1D*)directoryPHOSPi0PbPb1020->Get("hPi0_PbPb_cen2_NoBW_Stat");
    TH1D*histoPi0PHOSSysPbPb1020 =   (TH1D*)directoryPHOSPi0PbPb1020->Get("hPi0_PbPb_cen2_NoBW_Syst");
    TGraphAsymmErrors*graphPHOSYieldPi0SysErrPbPb1020 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb1020);

    histoarrayPi0InvYieldPbPb2760GeV[3][0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_1020->Clone("histoPCMPi0InvYieldPbPb2760GeV_1020");
    grapharrayPi0InvYieldStatPbPb2760GeV[3][0] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[3][0]);
    grapharrayPi0InvYieldSysPbPb2760GeV[3][0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_1020->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_1020");

    histoarrayPi0InvYieldPbPb2760GeV[3][1] = (TH1D*)histoPi0PHOSPbPb1020->Clone("histoPHOSPi0InvYieldPbPb2760GeV_1020");
    grapharrayPi0InvYieldStatPbPb2760GeV[3][1] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[3][1]);
    grapharrayPi0InvYieldSysPbPb2760GeV[3][1] = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb1020->Clone("graphPHOSPi0InvYieldSysPbPb2760GeV_1020");


    cout << "For the Pi0 in 20-40% " << endl;
    cout << "PCM" << endl;
    TDirectoryFile* directoryPCMPi0PbPb2760GeV_2040                 = (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_20-40%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_2040                   = (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("CorrectedYieldPi0");
    TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_2040                   = (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("CorrectedYieldBinShiftedPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_2040       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0SystError");
    cout << "PHOS" << endl;
    TDirectory* directoryPHOSPi0PbPb2040 = (TDirectory*)filePHOS->Get("pi0_PbPb_2760_Centrality_20-40%");
    TH1D*histoPi0PHOSPbPb2040 =      (TH1D*)directoryPHOSPi0PbPb2040->Get("hPi0_PbPb_cen3_NoBW_Stat");
    TH1D*histoPi0PHOSSysPbPb2040 =   (TH1D*)directoryPHOSPi0PbPb2040->Get("hPi0_PbPb_cen3_NoBW_Syst");
    TGraphAsymmErrors*graphPHOSYieldPi0SysErrPbPb2040 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb2040);

    histoarrayPi0InvYieldPbPb2760GeV[4][0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_2040->Clone("histoPCMPi0InvYieldPbPb2760GeV_2040");
    grapharrayPi0InvYieldStatPbPb2760GeV[4][0] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[4][0]);
    grapharrayPi0InvYieldSysPbPb2760GeV[4][0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2040->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_2040");

    histoarrayPi0InvYieldPbPb2760GeV[4][1] = (TH1D*)histoPi0PHOSPbPb2040->Clone("histoPHOSPi0InvYieldPbPb2760GeV_2040");
    grapharrayPi0InvYieldStatPbPb2760GeV[4][1] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[4][1]);
    grapharrayPi0InvYieldSysPbPb2760GeV[4][1] = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb2040->Clone("graphPHOSPi0InvYieldSysPbPb2760GeV_2040");


    cout << "For the Pi0 in 20-50% " << endl;
    cout << "PCM" << endl;
    TDirectory* directoryPCMPi0PbPb2760GeV_2050                     = (TDirectory*)filePCM->Get("Pi0_PbPb_2.76TeV_20-50%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_2050                   = (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("CorrectedYieldPi0");
    TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_2050                  = (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("CorrectedYieldBinShiftedPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_2050       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0SystError");
    cout << "EMCal" << endl;
    TH1D* histoEMCalPi0InvYieldPbPb2760GeV_2050         = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("LuciaBinningInvYieldPbPbStatErrPi_2050");
    TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_2050        = new TGraphAsymmErrors(histoEMCalPi0InvYieldPbPb2760GeV_2050);
    graphEMCalPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_2050->GetN()-1);
    graphEMCalPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_2050->GetN()-1);
    TH1D* histoEMCalPi0InvYieldSysPbPb2760GeV_2050  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("LuciaBinningInvYieldPbPbSysErrPi_2050");
    TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_2050     = new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_2050);
    graphEMCalPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_2050->GetN()-1);
    graphEMCalPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_2050->GetN()-1);

    histoarrayPi0InvYieldPbPb2760GeV[5][0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_2050->Clone("histoPCMPi0InvYieldPbPb2760GeV_2050");
    grapharrayPi0InvYieldStatPbPb2760GeV[5][0] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[5][0]);
    grapharrayPi0InvYieldSysPbPb2760GeV[5][0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2050->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_2050");

    grapharrayPi0InvYieldStatPbPb2760GeV[5][1] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[4][1]);
    grapharrayPi0InvYieldSysPbPb2760GeV[5][1] = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb2040->Clone("graphPHOSPi0InvYieldSysPbPb2760GeV_2050");
    grapharrayPi0InvYieldStatPbPb2760GeV[5][1]->Set(0);
    grapharrayPi0InvYieldSysPbPb2760GeV[5][1]->Set(0);

    histoarrayPi0InvYieldPbPb2760GeV[5][2] = (TH1D*)histoEMCalPi0InvYieldPbPb2760GeV_2050->Clone("histoEMCalPi0InvYieldPbPb2760GeV_2050");
    grapharrayPi0InvYieldStatPbPb2760GeV[5][2] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[5][2]);
    grapharrayPi0InvYieldStatPbPb2760GeV[5][2]->RemovePoint(grapharrayPi0InvYieldStatPbPb2760GeV[5][2]->GetN()-1);
    grapharrayPi0InvYieldStatPbPb2760GeV[5][2]->RemovePoint(grapharrayPi0InvYieldStatPbPb2760GeV[5][2]->GetN()-1);
    grapharrayPi0InvYieldSysPbPb2760GeV[5][2] = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Clone("graphEMCalPi0InvYieldSysPbPb2760GeV_2050");

    cout << "For the Pi0 in 40-60% " << endl;
    cout << "PCM" << endl;
    TDirectoryFile* directoryPCMPi0PbPb2760GeV_4060         = (TDirectoryFile*)filePCMPub->Get("Pi0_PbPb_2.76TeV_40-60%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_4060                   = (TH1D*)directoryPCMPi0PbPb2760GeV_4060->Get("CorrectedYieldPi0");
    TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_4060                   = (TH1D*)directoryPCMPi0PbPb2760GeV_4060->Get("CorrectedYieldBinShiftedPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_4060       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_4060->Get("Pi0SystError");
    cout << "PHOS" << endl;
    TDirectory* directoryPHOSPi0PbPb4060 = (TDirectory*)filePHOS->Get("pi0_PbPb_2760_Centrality_40-60%");
    TH1D*histoPi0PHOSPbPb4060 =      (TH1D*)directoryPHOSPi0PbPb4060->Get("hPi0_PbPb_cen4_NoBW_Stat");
    TH1D*histoPi0PHOSSysPbPb4060 =   (TH1D*)directoryPHOSPi0PbPb4060->Get("hPi0_PbPb_cen4_NoBW_Syst");
    TGraphAsymmErrors*graphPHOSYieldPi0SysErrPbPb4060 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb4060);

    histoarrayPi0InvYieldPbPb2760GeV[6][0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_4060->Clone("histoPCMPi0InvYieldPbPb2760GeV_4060");
    grapharrayPi0InvYieldStatPbPb2760GeV[6][0] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[6][0]);
    grapharrayPi0InvYieldSysPbPb2760GeV[6][0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_4060->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_4060");

    histoarrayPi0InvYieldPbPb2760GeV[6][1] = (TH1D*)histoPi0PHOSPbPb4060->Clone("histoPHOSPi0InvYieldPbPb2760GeV_4060");
    grapharrayPi0InvYieldStatPbPb2760GeV[6][1] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[6][1]);
    grapharrayPi0InvYieldSysPbPb2760GeV[6][1] = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb4060->Clone("graphPHOSPi0InvYieldSysPbPb2760GeV_4060");


    cout << "For the Pi0 in 60-80% " << endl;
    cout << "PCM" << endl;
    TDirectoryFile* directoryPCMPi0PbPb2760GeV_6080         = (TDirectoryFile*)filePCMPub->Get("Pi0_PbPb_2.76TeV_60-80%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_6080                   = (TH1D*)directoryPCMPi0PbPb2760GeV_6080->Get("CorrectedYieldPi0");
    TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_6080                   = (TH1D*)directoryPCMPi0PbPb2760GeV_6080->Get("CorrectedYieldBinShiftedPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_6080       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_6080->Get("Pi0SystError");
    cout << "PHOS" << endl;
    TDirectory* directoryPHOSPi0PbPb6080 = (TDirectory*)filePHOS->Get("pi0_PbPb_2760_Centrality_60-80%");
    TH1D*histoPi0PHOSPbPb6080 =      (TH1D*)directoryPHOSPi0PbPb6080->Get("hPi0_PbPb_cen5_NoBW_Stat");
    TH1D*histoPi0PHOSSysPbPb6080 =   (TH1D*)directoryPHOSPi0PbPb6080->Get("hPi0_PbPb_cen5_NoBW_Syst");
    TGraphAsymmErrors*graphPHOSYieldPi0SysErrPbPb6080 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb6080);

    histoarrayPi0InvYieldPbPb2760GeV[7][0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_6080->Clone("histoPCMPi0InvYieldPbPb2760GeV_6080");
    grapharrayPi0InvYieldStatPbPb2760GeV[7][0] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[7][0]);
    grapharrayPi0InvYieldSysPbPb2760GeV[7][0] = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_6080->Clone("graphPCMPi0InvYieldSysPbPb2760GeV_6080");

    histoarrayPi0InvYieldPbPb2760GeV[7][1] = (TH1D*)histoPi0PHOSPbPb6080->Clone("histoPHOSPi0InvYieldPbPb2760GeV_6080");
    grapharrayPi0InvYieldStatPbPb2760GeV[7][1] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[7][1]);
    grapharrayPi0InvYieldSysPbPb2760GeV[7][1] = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb6080->Clone("graphPHOSPi0InvYieldSysPbPb2760GeV_6080");


    for(Int_t c=0; c<nCent; c++){
      if(c==0 || c==1 || c==3 || c==4 || c==6 || c==7){
        grapharrayPi0InvYieldStatPbPb2760GeV[c][2] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[2][2]);
        grapharrayPi0InvYieldSysPbPb2760GeV[c][2] = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Clone("graphEMCalEtaInvYieldSysPbPb2760GeV_2040");
        grapharrayPi0InvYieldStatPbPb2760GeV[c][2]->Set(0);
        grapharrayPi0InvYieldSysPbPb2760GeV[c][2]->Set(0);
      }
    }




    cout << "For the Eta in 0-5% " << endl;
    cout << "PCM" << endl;
    TDirectory* directoryPCMEtaPbPb2760GeV_0005                     = (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_0-5%");
    TH1D* histoPCMEtaInvYieldPbPb2760GeV_0005                   = (TH1D*)directoryPCMEtaPbPb2760GeV_0005->Get("CorrectedYieldEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_0005       = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0005->Get("EtaSystError");
    TH1D* histoPCMEtatoPi0Stat2760GeV_0005                  = (TH1D*)directoryPCMEtaPbPb2760GeV_0005->Get("EtatoPi0Ratio");
    TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_0005  = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0005->Get("EtatoPi0RatioSys");

    histoarrayEtaInvYieldPbPb2760GeV[0][0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_0005->Clone("histoPCMEtaInvYieldPbPb2760GeV_0005");
    grapharrayEtaInvYieldStatPbPb2760GeV[0][0] = new TGraphAsymmErrors(histoarrayEtaInvYieldPbPb2760GeV[0][0]);
    grapharrayEtaInvYieldSysPbPb2760GeV[0][0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0005->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_0005");
    histoarrayEtatoPi0Stat2760GeV[0][0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_0005->Clone("histoPCMEtatoPi0Stat2760GeV_0005");
    grapharrayEtaToPi0Stat2760GeV[0][0] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[0][0]);
    grapharrayEtaToPi0Stat2760GeV[0][0]->RemovePoint(0);
    grapharrayEtaToPi0Stat2760GeV[0][0]->RemovePoint(0);
    grapharrayEtaToPi0Sys2760GeV[0][0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_0005->Clone("graphPCMEtatoPi0Sys2760GeV_0005");

    cout << "For the Eta in 5-10% " << endl;
    cout << "PCM" << endl;
    TDirectory* directoryPCMEtaPbPb2760GeV_0510                     = (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_5-10%");
    TH1D* histoPCMEtaInvYieldPbPb2760GeV_0510                   = (TH1D*)directoryPCMEtaPbPb2760GeV_0510->Get("CorrectedYieldEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_0510       = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0510->Get("EtaSystError");
    TH1D* histoPCMEtatoPi0Stat2760GeV_0510                  = (TH1D*)directoryPCMEtaPbPb2760GeV_0510->Get("EtatoPi0Ratio");
    TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_0510  = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0510->Get("EtatoPi0RatioSys");

    histoarrayEtaInvYieldPbPb2760GeV[1][0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_0510->Clone("histoPCMEtaInvYieldPbPb2760GeV_0510");
    grapharrayEtaInvYieldStatPbPb2760GeV[1][0] = new TGraphAsymmErrors(histoarrayEtaInvYieldPbPb2760GeV[1][0]);
    grapharrayEtaInvYieldSysPbPb2760GeV[1][0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0510->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_0510");
    histoarrayEtatoPi0Stat2760GeV[1][0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_0510->Clone("histoPCMEtatoPi0Stat2760GeV_0510");
    grapharrayEtaToPi0Stat2760GeV[1][0] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[1][0]);
    grapharrayEtaToPi0Stat2760GeV[1][0]->RemovePoint(0);
    grapharrayEtaToPi0Stat2760GeV[1][0]->RemovePoint(0);
    grapharrayEtaToPi0Sys2760GeV[1][0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_0510->Clone("graphPCMEtatoPi0Sys2760GeV_0510");


    cout << "For the Eta in 0-10% " << endl;
    cout << "PCM" << endl;
    TDirectory* directoryPCMEtaPbPb2760GeV_0010                     = (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_0-10%");
    TH1D* histoPCMEtaInvYieldPbPb2760GeV_0010                   = (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("CorrectedYieldEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_0010       = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaSystError");
    TH1D* histoPCMEtatoPi0Stat2760GeV_0010                  = (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("EtatoPi0Ratio");
    TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_0010  = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtatoPi0RatioSys");
    cout << "EMCal" << endl;
    TDirectory* directoryEMCalEtaPbPb2760GeV                = (TDirectory*)fileEMCal->Get("Eta2.76TeV_PbPb");
    TH1D* histoEMCalEtaInvYieldPbPb2760GeV_0010         = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbStatErrEta_0010");
    TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_0010        = new TGraphAsymmErrors(histoEMCalEtaInvYieldPbPb2760GeV_0010);
    graphEMCalEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_0010->GetN()-1);
    graphEMCalEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_0010->GetN()-1);
    TH1D*   histoEMCalEtaInvYieldSysPbPb2760GeV_0010    = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbSysErrEta_0010");
    TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_0010     = new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_0010);
    graphEMCalEtaInvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_0010->GetN()-1);
    graphEMCalEtaInvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_0010->GetN()-1);
    TH1D*   histoEMCalEtatoPi0StatPbPb2760GeV_0010  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("StatErrEtatoPi0Ratio_0010");
    TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_0010   = new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV_0010);
    graphEMCalEtatoPi0Stat2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_0010->GetN()-1);
    graphEMCalEtatoPi0Stat2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_0010->GetN()-1);
    TH1D*   histoEMCalEtatoPi0SysPbPb2760GeV_0010   = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("SysErrEtatoPi0Ratio_0010");
    TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_0010    = new TGraphAsymmErrors(histoEMCalEtatoPi0SysPbPb2760GeV_0010);
    graphEMCalEtatoPi0Sys2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_0010->GetN()-1);
    graphEMCalEtatoPi0Sys2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_0010->GetN()-1);

    histoarrayEtaInvYieldPbPb2760GeV[2][0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_0010->Clone("histoPCMEtaInvYieldPbPb2760GeV_0010");
    grapharrayEtaInvYieldStatPbPb2760GeV[2][0] = new TGraphAsymmErrors(histoarrayEtaInvYieldPbPb2760GeV[2][0]);
    grapharrayEtaInvYieldSysPbPb2760GeV[2][0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0010->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_0010");
    histoarrayEtatoPi0Stat2760GeV[2][0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_0010->Clone("histoPCMEtatoPi0Stat2760GeV_0010");
    grapharrayEtaToPi0Stat2760GeV[2][0] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[2][0]);
    grapharrayEtaToPi0Stat2760GeV[2][0]->RemovePoint(0);
    grapharrayEtaToPi0Stat2760GeV[2][0]->RemovePoint(0);
    grapharrayEtaToPi0Sys2760GeV[2][0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_0010->Clone("graphPCMEtatoPi0Sys2760GeV_0010");

    histoarrayEtaInvYieldPbPb2760GeV[2][2] = (TH1D*)histoEMCalEtaInvYieldPbPb2760GeV_0010->Clone("histoEMCalEtaInvYieldPbPb2760GeV_0010");
    grapharrayEtaInvYieldStatPbPb2760GeV[2][2] = new TGraphAsymmErrors(histoarrayEtaInvYieldPbPb2760GeV[2][2]);
    grapharrayEtaInvYieldStatPbPb2760GeV[2][2]->RemovePoint(grapharrayEtaInvYieldStatPbPb2760GeV[2][2]->GetN()-1);
    grapharrayEtaInvYieldStatPbPb2760GeV[2][2]->RemovePoint(grapharrayEtaInvYieldStatPbPb2760GeV[2][2]->GetN()-1);
    grapharrayEtaInvYieldSysPbPb2760GeV[2][2] = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Clone("graphEMCalEtaInvYieldSysPbPb2760GeV_0010");
    histoarrayEtatoPi0Stat2760GeV[2][2] = (TH1D*)histoEMCalEtatoPi0StatPbPb2760GeV_0010->Clone("histoEMCalEtatoPi0StatPbPb2760GeV_0010");
    grapharrayEtaToPi0Stat2760GeV[2][2] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[2][2]);
    grapharrayEtaToPi0Stat2760GeV[2][2]->RemovePoint(grapharrayEtaToPi0Stat2760GeV[2][2]->GetN()-1);
    grapharrayEtaToPi0Stat2760GeV[2][2]->RemovePoint(grapharrayEtaToPi0Stat2760GeV[2][2]->GetN()-1);
    grapharrayEtaToPi0Sys2760GeV[2][2] = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_0010->Clone("graphEMCalEtatoPi0Sys2760GeV_0010");

    cout << "For the Eta in 20-40% " << endl;
    cout << "PCM" << endl;
    TDirectory* directoryPCMEtaPbPb2760GeV_2040                    = (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_20-40%");
    TH1D* histoPCMEtaInvYieldPbPb2760GeV_2040                   = (TH1D*)directoryPCMEtaPbPb2760GeV_2040->Get("CorrectedYieldEta");
    TH1D* histoPCMEtaInvYieldPbPb2760GeVYshifted_2040                  = (TH1D*)directoryPCMEtaPbPb2760GeV_2040->Get("CorrectedYieldBinShiftedEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_2040       = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2040->Get("EtaSystError");
    TH1D* histoPCMEtatoPi0Stat2760GeV_2040                  = (TH1D*)directoryPCMEtaPbPb2760GeV_2040->Get("EtatoPi0Ratio");
    TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_2040  = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2040->Get("EtatoPi0RatioSys");

    histoarrayEtaInvYieldPbPb2760GeV[4][0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_2040->Clone("histoPCMEtaInvYieldPbPb2760GeV_2040");
    grapharrayEtaInvYieldStatPbPb2760GeV[4][0] = new TGraphAsymmErrors(histoarrayEtaInvYieldPbPb2760GeV[4][0]);
    grapharrayEtaInvYieldSysPbPb2760GeV[4][0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2040->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_2040");
    histoarrayEtatoPi0Stat2760GeV[4][0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_2040->Clone("histoPCMEtatoPi0Stat2760GeV_2040");
    grapharrayEtaToPi0Stat2760GeV[4][0] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[4][0]);
    grapharrayEtaToPi0Sys2760GeV[4][0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_2040->Clone("graphPCMEtatoPi0Sys2760GeV_2040");


    cout << "For the Eta in 20-50% " << endl;
    cout << "PCM" << endl;
    TDirectory* directoryPCMEtaPbPb2760GeV_2050                     = (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_20-50%");
    TH1D* histoPCMEtaInvYieldPbPb2760GeV_2050                   = (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("CorrectedYieldEta");
    TH1D* histoPCMEtaInvYieldPbPb2760GeVYshifted_2050                  = (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("CorrectedYieldBinShiftedEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_2050       = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaSystError");
    TH1D* histoPCMEtatoPi0Stat2760GeV_2050                  = (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("EtatoPi0Ratio");
    TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_2050  = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtatoPi0RatioSys");
    cout << "EMCal" << endl;
    TH1D* histoEMCalEtaInvYieldPbPb2760GeV_2050         = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbStatErrEta_2050");
    TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_2050        = new TGraphAsymmErrors(histoEMCalEtaInvYieldPbPb2760GeV_2050);
    graphEMCalEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_2050->GetN()-1);
    graphEMCalEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_2050->GetN()-1);
    TH1D*   histoEMCalEtaInvYieldSysPbPb2760GeV_2050    = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbSysErrEta_2050");
    TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_2050     = new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_2050);
    graphEMCalEtaInvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_2050->GetN()-1);
    graphEMCalEtaInvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_2050->GetN()-1);
    TH1D*   histoEMCalEtatoPi0StatPbPb2760GeV_2050  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("StatErrEtatoPi0Ratio_2050");
    TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_2050   = new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV_2050);
    graphEMCalEtatoPi0Stat2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_2050->GetN()-1);
    graphEMCalEtatoPi0Stat2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_2050->GetN()-1);
    TH1D*   histoEMCalEtatoPi0SysPbPb2760GeV_2050   = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("SysErrEtatoPi0Ratio_2050");
    TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_2050    = new TGraphAsymmErrors(histoEMCalEtatoPi0SysPbPb2760GeV_2050);
    graphEMCalEtatoPi0Sys2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_2050->GetN()-1);
    graphEMCalEtatoPi0Sys2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_2050->GetN()-1);

    histoarrayEtaInvYieldPbPb2760GeV[5][0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_2050->Clone("histoPCMEtaInvYieldPbPb2760GeV_2050");
    grapharrayEtaInvYieldStatPbPb2760GeV[5][0] = new TGraphAsymmErrors(histoarrayEtaInvYieldPbPb2760GeV[5][0]);
    grapharrayEtaInvYieldSysPbPb2760GeV[5][0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2050->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_2050");
    histoarrayEtatoPi0Stat2760GeV[5][0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_2050->Clone("histoPCMEtatoPi0Stat2760GeV_2050");
    grapharrayEtaToPi0Stat2760GeV[5][0] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[5][0]);
    grapharrayEtaToPi0Stat2760GeV[5][0]->RemovePoint(0);
    grapharrayEtaToPi0Stat2760GeV[5][0]->RemovePoint(0);
    grapharrayEtaToPi0Sys2760GeV[5][0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_2050->Clone("graphPCMEtatoPi0Sys2760GeV_2050");

    histoarrayEtaInvYieldPbPb2760GeV[5][2] = (TH1D*)histoEMCalEtaInvYieldPbPb2760GeV_2050->Clone("histoEMCalEtaInvYieldPbPb2760GeV_2050");
    grapharrayEtaInvYieldStatPbPb2760GeV[5][2] = new TGraphAsymmErrors(histoarrayEtaInvYieldPbPb2760GeV[5][2]);
    grapharrayEtaInvYieldStatPbPb2760GeV[5][2]->RemovePoint(grapharrayEtaInvYieldStatPbPb2760GeV[5][2]->GetN()-1);
    grapharrayEtaInvYieldStatPbPb2760GeV[5][2]->RemovePoint(grapharrayEtaInvYieldStatPbPb2760GeV[5][2]->GetN()-1);
    grapharrayEtaInvYieldSysPbPb2760GeV[5][2] = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Clone("graphEMCalEtaInvYieldSysPbPb2760GeV_2050");
    histoarrayEtatoPi0Stat2760GeV[5][2] = (TH1D*)histoEMCalEtatoPi0StatPbPb2760GeV_2050->Clone("histoEMCalEtatoPi0StatPbPb2760GeV_2050");
    grapharrayEtaToPi0Stat2760GeV[5][2] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[5][2]);
    grapharrayEtaToPi0Stat2760GeV[5][2]->RemovePoint(grapharrayEtaToPi0Stat2760GeV[5][2]->GetN()-1);
    grapharrayEtaToPi0Stat2760GeV[5][2]->RemovePoint(grapharrayEtaToPi0Stat2760GeV[5][2]->GetN()-1);
    grapharrayEtaToPi0Sys2760GeV[5][2] = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_2050->Clone("graphEMCalEtatoPi0Sys2760GeV_2050");

    for(Int_t c=0; c<nCent; c++){
      if(c==3 || c==6 || c==7){
        grapharrayEtaInvYieldStatPbPb2760GeV[c][0] = new TGraphAsymmErrors(histoarrayEtaInvYieldPbPb2760GeV[4][0]);
        grapharrayEtaInvYieldSysPbPb2760GeV[c][0] = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2040->Clone("graphPCMEtaInvYieldSysPbPb2760GeV_2040");
        grapharrayEtaInvYieldStatPbPb2760GeV[c][0]->Set(0);
        grapharrayEtaInvYieldSysPbPb2760GeV[c][0]->Set(0);
        grapharrayEtaToPi0Stat2760GeV[c][0] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[4][0]);
        grapharrayEtaToPi0Sys2760GeV[c][0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_2040->Clone("graphPCMEtatoPi0Sys2760GeV_2040");
        grapharrayEtaToPi0Stat2760GeV[c][0]->Set(0);
        grapharrayEtaToPi0Sys2760GeV[c][0]->Set(0);
      }

      if(c==0 || c==1 || c==3 || c==4 || c==6 || c==7){
        grapharrayEtaInvYieldStatPbPb2760GeV[c][2] = new TGraphAsymmErrors(histoarrayEtaInvYieldPbPb2760GeV[2][2]);
        grapharrayEtaInvYieldSysPbPb2760GeV[c][2] = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Clone("graphEMCalEtaInvYieldSysPbPb2760GeV_2040");
        grapharrayEtaInvYieldStatPbPb2760GeV[c][2]->Set(0);
        grapharrayEtaInvYieldSysPbPb2760GeV[c][2]->Set(0);
        grapharrayEtaToPi0Stat2760GeV[c][2] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[2][2]);
        grapharrayEtaToPi0Sys2760GeV[c][2] = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_0010->Clone("graphEMCalEtatoPi0Sys2760GeV_0010");
        grapharrayEtaToPi0Stat2760GeV[c][2]->Set(0);
        grapharrayEtaToPi0Sys2760GeV[c][2]->Set(0);
      }

      grapharrayEtaInvYieldStatPbPb2760GeV[c][1] = new TGraphAsymmErrors(histoarrayPi0InvYieldPbPb2760GeV[4][1]);
      grapharrayEtaInvYieldSysPbPb2760GeV[c][1] = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb2040->Clone("graphPHOSPi0InvYieldSysPbPb2760GeV_2050");
      grapharrayEtaInvYieldStatPbPb2760GeV[c][1]->Set(0);
      grapharrayEtaInvYieldSysPbPb2760GeV[c][1]->Set(0);
      grapharrayEtaToPi0Stat2760GeV[c][1] = new TGraphAsymmErrors(histoarrayEtatoPi0Stat2760GeV[2][2]);
      grapharrayEtaToPi0Sys2760GeV[c][1] = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_0010->Clone("graphEMCalEtatoPi0Sys2760GeV_0010");
      grapharrayEtaToPi0Stat2760GeV[c][1]->Set(0);
      grapharrayEtaToPi0Sys2760GeV[c][1]->Set(0);
    }


    //syst and stat relative errors
    TH1D *statErrorCollectionLHC11h[nCent][11];
    TH1D *statErrorCollectionEtatoPi0LHC11h[nCent][11];
    TGraphAsymmErrors *sysErrorCollectionLHC11h[nCent][11];
    TGraphAsymmErrors *sysErrorCollectionEtatoPi0LHC11h[nCent][11];

    for (Int_t i = 0; i< 11; i++){
      for (Int_t c = 0; c< nCent; c++){
        statErrorCollectionLHC11h[c][i] = NULL;
        sysErrorCollectionLHC11h[c][i] = NULL;
        statErrorCollectionEtatoPi0LHC11h[c][i] = NULL;
        sysErrorCollectionEtatoPi0LHC11h[c][i] = NULL;
      }
    }
    if(meson.CompareTo("Pi0")==0){
      statErrorCollectionLHC11h[0][0] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[0][0]->Clone("statErrPCMPi0_0005");
      statErrorCollectionLHC11h[1][0] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[1][0]->Clone("statErrPCMPi0_0510");
      statErrorCollectionLHC11h[2][0] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[2][0]->Clone("statErrPCMPi0_0010");
      statErrorCollectionLHC11h[3][0] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[3][0]->Clone("statErrPCMPi0_1020");
      statErrorCollectionLHC11h[4][0] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[4][0]->Clone("statErrPCMPi0_2040");
      statErrorCollectionLHC11h[5][0] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[5][0]->Clone("statErrPCMPi0_2050");
      statErrorCollectionLHC11h[6][0] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[6][0]->Clone("statErrPCMPi0_4060");
      statErrorCollectionLHC11h[7][0] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[7][0]->Clone("statErrPCMPi0_6080");

      sysErrorCollectionLHC11h[0][0] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[0][0]->Clone("sysErrPCMPi0_0005");
      sysErrorCollectionLHC11h[1][0] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[1][0]->Clone("sysErrPCMPi0_0510");
      sysErrorCollectionLHC11h[2][0] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[2][0]->Clone("sysErrPCMPi0_0010");
      sysErrorCollectionLHC11h[3][0] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[3][0]->Clone("sysErrPCMPi0_1020");
      sysErrorCollectionLHC11h[4][0] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[4][0]->Clone("sysErrPCMPi0_2040");
      sysErrorCollectionLHC11h[5][0] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[5][0]->Clone("sysErrPCMPi0_2050");
      sysErrorCollectionLHC11h[6][0] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[6][0]->Clone("sysErrPCMPi0_4060");
      sysErrorCollectionLHC11h[7][0] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[7][0]->Clone("sysErrPCMPi0_6080");

      statErrorCollectionLHC11h[0][1] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[0][1]->Clone("statErrPHOSPi0_0005");
      statErrorCollectionLHC11h[1][1] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[1][1]->Clone("statErrPHOSPi0_0510");
      statErrorCollectionLHC11h[2][1] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[2][1]->Clone("statErrPHOSPi0_0010");
      statErrorCollectionLHC11h[3][1] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[3][1]->Clone("statErrPHOSPi0_1020");
      statErrorCollectionLHC11h[4][1] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[4][1]->Clone("statErrPHOSPi0_2040");
      statErrorCollectionLHC11h[6][1] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[6][1]->Clone("statErrPHOSPi0_4060");
      statErrorCollectionLHC11h[7][1] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[7][1]->Clone("statErrPHOSPi0_6080");

      sysErrorCollectionLHC11h[0][1] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[0][1]->Clone("sysErrPHOSPi0_0005");
      sysErrorCollectionLHC11h[1][1] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[1][1]->Clone("sysErrPHOSPi0_0510");
      sysErrorCollectionLHC11h[2][1] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[2][1]->Clone("sysErrPHOSPi0_0010");
      sysErrorCollectionLHC11h[3][1] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[3][1]->Clone("sysErrPHOSPi0_1020");
      sysErrorCollectionLHC11h[4][1] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[4][1]->Clone("sysErrPHOSPi0_2040");
      sysErrorCollectionLHC11h[6][1] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[6][1]->Clone("sysErrPHOSPi0_4060");
      sysErrorCollectionLHC11h[7][1] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[7][1]->Clone("sysErrPHOSPi0_6080");

      statErrorCollectionLHC11h[2][2] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[2][2]->Clone("statErrEMCalPi0_0010");
      statErrorCollectionLHC11h[5][2] = (TH1D*)histoarrayPi0InvYieldPbPb2760GeV[5][2]->Clone("statErrEMCalPi0_2050");

      sysErrorCollectionLHC11h[2][2] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[2][2]->Clone("sysErrEMCalPi0_0010");
      sysErrorCollectionLHC11h[5][2] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[5][2]->Clone("sysErrEMCalPi0_2050");

    } else if(meson.CompareTo("Eta")==0) {

      statErrorCollectionLHC11h[0][0] = (TH1D*)histoarrayEtaInvYieldPbPb2760GeV[0][0]->Clone("statErrPCMEta_0005");
      statErrorCollectionLHC11h[1][0] = (TH1D*)histoarrayEtaInvYieldPbPb2760GeV[1][0]->Clone("statErrPCMEta_0510");
      statErrorCollectionLHC11h[2][0] = (TH1D*)histoarrayEtaInvYieldPbPb2760GeV[2][0]->Clone("statErrPCMEta_0010");
      statErrorCollectionLHC11h[4][0] = (TH1D*)histoarrayEtaInvYieldPbPb2760GeV[4][0]->Clone("statErrPCMEta_2040");
      statErrorCollectionLHC11h[5][0] = (TH1D*)histoarrayEtaInvYieldPbPb2760GeV[5][0]->Clone("statErrPCMEta_2050");

      sysErrorCollectionLHC11h[0][0] = (TGraphAsymmErrors*)grapharrayEtaInvYieldSysPbPb2760GeV[0][0]->Clone("sysErrPCMEta_0005");
      sysErrorCollectionLHC11h[1][0] = (TGraphAsymmErrors*)grapharrayEtaInvYieldSysPbPb2760GeV[1][0]->Clone("sysErrPCMEta_0510");
      sysErrorCollectionLHC11h[2][0] = (TGraphAsymmErrors*)grapharrayEtaInvYieldSysPbPb2760GeV[2][0]->Clone("sysErrPCMEta_0010");
      sysErrorCollectionLHC11h[4][0] = (TGraphAsymmErrors*)grapharrayEtaInvYieldSysPbPb2760GeV[4][0]->Clone("sysErrPCMEta_2040");
      sysErrorCollectionLHC11h[5][0] = (TGraphAsymmErrors*)grapharrayEtaInvYieldSysPbPb2760GeV[5][0]->Clone("sysErrPCMEta_2050");

      statErrorCollectionEtatoPi0LHC11h[0][0] = (TH1D*)histoarrayEtatoPi0Stat2760GeV[0][0]->Clone("statErrPCMEtatoPi0_0005");
      statErrorCollectionEtatoPi0LHC11h[1][0] = (TH1D*)histoarrayEtatoPi0Stat2760GeV[1][0]->Clone("statErrPCMEtatoPi0_0510");
      statErrorCollectionEtatoPi0LHC11h[2][0] = (TH1D*)histoarrayEtatoPi0Stat2760GeV[2][0]->Clone("statErrPCMEtatoPi0_0010");
      statErrorCollectionEtatoPi0LHC11h[4][0] = (TH1D*)histoarrayEtatoPi0Stat2760GeV[4][0]->Clone("statErrPCMEtatoPi0_2040");
      statErrorCollectionEtatoPi0LHC11h[5][0] = (TH1D*)histoarrayEtatoPi0Stat2760GeV[5][0]->Clone("statErrPCMEtatoPi0_2050");

      sysErrorCollectionEtatoPi0LHC11h[0][0] = (TGraphAsymmErrors*)grapharrayEtaToPi0Sys2760GeV[0][0]->Clone("sysErrPCMEtatoPi0_0005");
      sysErrorCollectionEtatoPi0LHC11h[1][0] = (TGraphAsymmErrors*)grapharrayEtaToPi0Sys2760GeV[1][0]->Clone("sysErrPCMEtatoPi0_0510");
      sysErrorCollectionEtatoPi0LHC11h[2][0] = (TGraphAsymmErrors*)grapharrayEtaToPi0Sys2760GeV[2][0]->Clone("sysErrPCMEtatoPi0_0010");
      sysErrorCollectionEtatoPi0LHC11h[4][0] = (TGraphAsymmErrors*)grapharrayEtaToPi0Sys2760GeV[4][0]->Clone("sysErrPCMEtatoPi0_2040");
      sysErrorCollectionEtatoPi0LHC11h[5][0] = (TGraphAsymmErrors*)grapharrayEtaToPi0Sys2760GeV[5][0]->Clone("sysErrPCMEtatoPi0_2050");

      statErrorCollectionLHC11h[2][2] = (TH1D*)histoarrayEtaInvYieldPbPb2760GeV[2][2]->Clone("statErrEMCalEta_0010");
      statErrorCollectionLHC11h[5][2] = (TH1D*)histoarrayEtaInvYieldPbPb2760GeV[5][2]->Clone("statErrEMCalEta_2050");

      sysErrorCollectionLHC11h[2][2] = (TGraphAsymmErrors*)grapharrayEtaInvYieldSysPbPb2760GeV[2][2]->Clone("sysErrEMCalEta_0010");
      sysErrorCollectionLHC11h[5][2] = (TGraphAsymmErrors*)grapharrayEtaInvYieldSysPbPb2760GeV[5][2]->Clone("sysErrEMCalEta_2050");

      statErrorCollectionEtatoPi0LHC11h[2][2] = (TH1D*)histoarrayEtatoPi0Stat2760GeV[2][2]->Clone("statErrEMCalEtatoPi0_0010");
      statErrorCollectionEtatoPi0LHC11h[5][2] = (TH1D*)histoarrayEtatoPi0Stat2760GeV[5][2]->Clone("statErrEMCalEtatoPi0_2050");

      sysErrorCollectionEtatoPi0LHC11h[2][2] = (TGraphAsymmErrors*)grapharrayEtaToPi0Sys2760GeV[2][2]->Clone("sysErrEMCalEtatoPi0_0010");
      sysErrorCollectionEtatoPi0LHC11h[5][2] = (TGraphAsymmErrors*)grapharrayEtaToPi0Sys2760GeV[5][2]->Clone("sysErrEMCalEtatoPi0_2050");

    }


    //Combining the spectra together
    TString fileNameOutputWeighting  = Form("%s/WeightingMethod%s.dat",outputDir.Data(),meson.Data());
    TGraphAsymmErrors* grapharrayInvYieldCombTotPbPb2760GeV[nCent];
    TGraphAsymmErrors* grapharrayInvYieldCombStatPbPb2760GeV[nCent];
    TGraphAsymmErrors* grapharrayInvYieldCombSysPbPb2760GeV[nCent];
    for(Int_t c=0; c<nCent; c++){
        grapharrayInvYieldCombTotPbPb2760GeV[c] = NULL;
        grapharrayInvYieldCombStatPbPb2760GeV[c] = NULL;
        grapharrayInvYieldCombSysPbPb2760GeV[c] = NULL;
    }
    Int_t offSetsCombPi0[11]    =   { 0, 2, 15, 0, 0,
                                  0, 0, 0, 0, 0, 0};
    Int_t offSetsCombPi0Sys[11]=    { 1, 2, 15, 0, 0,
                                  0, 0, 0, 0, 0, 1};
    Int_t offSetsCombEta[11]    =   { 0, 0, 6, 0, 0,
                                  0, 0, 0, 0, 0, 0};
    Int_t offSetsCombEtaSys[11]=    { 2, 0, 6, 0, 0,
                                  0, 0, 0, 0, 0, 0};

    if(meson.CompareTo("Pi0")==0){
        // Declaration & calculation of combined spectrum
      for(Int_t c=0; c<nCent; c++){
        cout << "***********************************************" << endl;
        cout << "**************** centrality " << c << " *****************" << endl;
        cout << "***********************************************" << endl;
        grapharrayInvYieldCombTotPbPb2760GeV[c] = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h[c], sysErrorCollectionLHC11h[c],
                                                                                     xPtLimitsPi0, 23, offSetsCombPi0, offSetsCombPi0Sys,
                                                                                     grapharrayInvYieldCombStatPbPb2760GeV[c], grapharrayInvYieldCombSysPbPb2760GeV[c],
                                                                                     fileNameOutputWeighting,1 );
      }

    } else  if(meson.CompareTo("Eta")==0){
        // Declaration & calculation of combined spectrum
      for(Int_t c=0; c<nCent; c++){
        cout << "***********************************************" << endl;
        cout << "**************** centrality " << c << " *****************" << endl;
        cout << "***********************************************" << endl;
        grapharrayInvYieldCombTotPbPb2760GeV[c] = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h[c], sysErrorCollectionLHC11h[c],
                                                                                     xPtLimitsEta, 13, offSetsCombEta, offSetsCombEtaSys,
                                                                                     grapharrayInvYieldCombStatPbPb2760GeV[c], grapharrayInvYieldCombSysPbPb2760GeV[c],
                                                                                     fileNameOutputWeighting,1 );

        grapharrayInvYieldCombStatPbPb2760GeV[c]->RemovePoint(0);
        grapharrayInvYieldCombStatPbPb2760GeV[c]->RemovePoint(0);
        grapharrayInvYieldCombSysPbPb2760GeV[c]->RemovePoint(0);
        grapharrayInvYieldCombSysPbPb2760GeV[c]->RemovePoint(0);
        grapharrayInvYieldCombTotPbPb2760GeV[c]->RemovePoint(0);
        grapharrayInvYieldCombTotPbPb2760GeV[c]->RemovePoint(0);
      }

    }


    //**********************************************************************************************************************//
    //************************************* Calculating bin shifted spectra & fitting **************************************//
    //**********************************************************************************************************************//
    //cloning spectra for shifting
    TGraphAsymmErrors* grapharrayInvYieldCombTotPbPb2760GeVUnshifted[nCent];
    TGraphAsymmErrors* grapharrayInvYieldCombStatPbPb2760GeVUnshifted[nCent];
    TGraphAsymmErrors* grapharrayInvYieldCombSysPbPb2760GeVUnshifted[nCent];
    TGraphAsymmErrors* grapharrayInvYieldTotPbPb2760GeVUnshifted[nCent][method];
    TGraphAsymmErrors* grapharrayInvYieldStatPbPb2760GeVUnshifted[nCent][method];
    TGraphAsymmErrors* grapharrayInvYieldSysPbPb2760GeVUnshifted[nCent][method];
    TGraphAsymmErrors* grapharrayRatioToFitStatPbPb2760GeV[nCent][method];
    TGraphAsymmErrors* grapharrayRatioToFitSysPbPb2760GeV[nCent][method];
    for(Int_t c=0; c<nCent; c++){
      grapharrayInvYieldCombTotPbPb2760GeVUnshifted[c] = NULL;
      grapharrayInvYieldCombStatPbPb2760GeVUnshifted[c] = NULL;
      grapharrayInvYieldCombSysPbPb2760GeVUnshifted[c] = NULL;
      for(Int_t m=0; m<3; m++){
        grapharrayInvYieldTotPbPb2760GeVUnshifted[c][m] = NULL;
        grapharrayInvYieldStatPbPb2760GeVUnshifted[c][m] = NULL;
        grapharrayInvYieldSysPbPb2760GeVUnshifted[c][m] = NULL;
        grapharrayRatioToFitStatPbPb2760GeV[c][m] = NULL;
        grapharrayRatioToFitSysPbPb2760GeV[c][m] = NULL;
      }
    }
    for(Int_t c = 0; c< nCent; c++){
      grapharrayInvYieldCombTotPbPb2760GeVUnshifted[c] = (TGraphAsymmErrors*)grapharrayInvYieldCombTotPbPb2760GeV[c]->Clone(Form("UnshiftedSpectrumCombTot_%d",c));
      grapharrayInvYieldCombStatPbPb2760GeVUnshifted[c] = (TGraphAsymmErrors*)grapharrayInvYieldCombStatPbPb2760GeV[c]->Clone(Form("UnshiftedSpectrumCombStat_%d",c));
      grapharrayInvYieldCombSysPbPb2760GeVUnshifted[c] = (TGraphAsymmErrors*)grapharrayInvYieldCombSysPbPb2760GeV[c]->Clone(Form("UnshiftedSpectrumCombSys_%d",c));

      for(Int_t m=0; m<3; m++){
        cout << "cent " << c << " method " << m << endl;
        if(meson.CompareTo("Pi0")==0){
          grapharrayInvYieldStatPbPb2760GeV[c][m] = (TGraphAsymmErrors*)grapharrayPi0InvYieldStatPbPb2760GeV[c][m]->Clone("grapharrayInvYieldStatPbPb2760GeV");
          grapharrayInvYieldSysPbPb2760GeV[c][m] = (TGraphAsymmErrors*)grapharrayPi0InvYieldSysPbPb2760GeV[c][m]->Clone("grapharrayInvYieldSysPbPb2760GeV");
        } else if(meson.CompareTo("Eta")==0){
          grapharrayInvYieldStatPbPb2760GeV[c][m] = (TGraphAsymmErrors*)grapharrayEtaInvYieldStatPbPb2760GeV[c][m]->Clone("grapharrayInvYieldStatPbPb2760GeV");
          grapharrayInvYieldSysPbPb2760GeV[c][m] = (TGraphAsymmErrors*)grapharrayEtaInvYieldSysPbPb2760GeV[c][m]->Clone("grapharrayInvYieldSysPbPb2760GeV");
        }
        grapharrayInvYieldStatPbPb2760GeVUnshifted[c][m] = (TGraphAsymmErrors*)grapharrayInvYieldStatPbPb2760GeV[c][m]->Clone(Form("UnshiftedSpectrumStat_%d_%d",c,m));
        grapharrayInvYieldSysPbPb2760GeVUnshifted[c][m] = (TGraphAsymmErrors*)grapharrayInvYieldSysPbPb2760GeV[c][m]->Clone(Form("UnshiftedSpectrumSys_%d_%d",c,m));
      }
    }

    if(meson.CompareTo("Pi0")==0){
        for(Int_t c=0; c<nCent; c++){
//           cout << "cent " << c << " combined" << endl;
          if(!(c==2 || c==5)){
            if(c==3 || c==6 || c==7){
              grapharrayInvYieldCombTotPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN()-1);
              grapharrayInvYieldCombStatPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombStatPbPb2760GeV[c]->GetN()-1);
              grapharrayInvYieldCombSysPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombSysPbPb2760GeV[c]->GetN()-1);
            }
            grapharrayInvYieldCombTotPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombTotPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombStatPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombStatPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombSysPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombSysPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombStatPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombStatPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombSysPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombSysPbPb2760GeV[c]->GetN()-1);
          }
//           grapharrayInvYieldCombStatPbPb2760GeV[c]->Print();
//           grapharrayInvYieldCombSysPbPb2760GeV[c]->Print();
          for(Int_t m=0; m<3; m++){
//             cout << "cent " << c << " method " << m << endl;
            if(m==0){
              grapharrayInvYieldStatPbPb2760GeV[c][m]->RemovePoint(0);
            } else if(m==1){
              grapharrayInvYieldStatPbPb2760GeV[c][m]->RemovePoint(grapharrayInvYieldStatPbPb2760GeV[c][m]->GetN()-1);
              grapharrayInvYieldStatPbPb2760GeV[c][m]->RemovePoint(grapharrayInvYieldStatPbPb2760GeV[c][m]->GetN()-1);
              grapharrayInvYieldSysPbPb2760GeV[c][m]->RemovePoint(grapharrayInvYieldSysPbPb2760GeV[c][m]->GetN()-1);
              grapharrayInvYieldSysPbPb2760GeV[c][m]->RemovePoint(grapharrayInvYieldSysPbPb2760GeV[c][m]->GetN()-1);
              grapharrayInvYieldStatPbPb2760GeV[c][m]->RemovePoint(0);
              grapharrayInvYieldStatPbPb2760GeV[c][m]->RemovePoint(0);
              grapharrayInvYieldSysPbPb2760GeV[c][m]->RemovePoint(0);
              grapharrayInvYieldSysPbPb2760GeV[c][m]->RemovePoint(0);
            }
//             grapharrayInvYieldStatPbPb2760GeV[c][m]->Print();
//             grapharrayInvYieldSysPbPb2760GeV[c][m]->Print();
          }
        }
    } else if(meson.CompareTo("Eta")==0){
        for(Int_t c=0; c<nCent; c++){
//           cout << "cent " << c << " combined" << endl;
          if(c==3 || c==6 || c==7){
            grapharrayInvYieldCombTotPbPb2760GeV[c]->Set(0);
            grapharrayInvYieldCombStatPbPb2760GeV[c]->Set(0);
            grapharrayInvYieldCombSysPbPb2760GeV[c]->Set(0);
          }
          if(!(c==2 || c==5)){
            grapharrayInvYieldCombTotPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombTotPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombTotPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombTotPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombStatPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombStatPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombStatPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombStatPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombStatPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombStatPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombStatPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombStatPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombSysPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombSysPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombSysPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombSysPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombSysPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombSysPbPb2760GeV[c]->GetN()-1);
            grapharrayInvYieldCombSysPbPb2760GeV[c]->RemovePoint(grapharrayInvYieldCombSysPbPb2760GeV[c]->GetN()-1);
          }
//           grapharrayInvYieldCombStatPbPb2760GeV[c]->Print();
//           grapharrayInvYieldCombSysPbPb2760GeV[c]->Print();
          for(Int_t m=0; m<3; m++){
//             cout << "cent " << c << " method " << m << endl;
            if(m==0){
              grapharrayInvYieldStatPbPb2760GeV[c][m]->RemovePoint(0);
              grapharrayInvYieldStatPbPb2760GeV[c][m]->RemovePoint(0);
            }
//             grapharrayInvYieldStatPbPb2760GeV[c][m]->Print();
//             grapharrayInvYieldSysPbPb2760GeV[c][m]->Print();
          }
        }
    }

    Double_t fitStart;
    Double_t fitStop;
    Int_t rangePi0LowComb[3] = {0,3,14};
    Int_t rangePi0HighComb[3] = {19,18,21};
    Int_t rangeEtaLowComb[3] = {0,0,4};
    Int_t rangeEtaHighComb[3] = {6,0,10};
    if(meson.CompareTo("Pi0")==0){
      fitStart = 0.4;
      fitStop = 14.;
    } else if(meson.CompareTo("Eta")==0){
      fitStart = 1.;
      fitStop = 10.;
    }



    for(Int_t c = 0; c< nCent; c++){

      //Calculating binshifts
      if(bWCorrection.CompareTo("X")==0 ){

        if(meson.CompareTo("Pi0")==0){

          TF1* fitFunctionShiftingX = FitObject("tcmpt","tcmptPi0","Pi0",grapharrayInvYieldCombTotPbPb2760GeV[c]);
          //combined
          grapharrayInvYieldCombTotPbPb2760GeV[c]         = ApplyXshift(grapharrayInvYieldCombTotPbPb2760GeV[c], fitFunctionShiftingX,"Pi0");
          grapharrayInvYieldCombStatPbPb2760GeV[c]        = ApplyXshiftIndividualSpectra (grapharrayInvYieldCombTotPbPb2760GeV[c],
                                                                                          grapharrayInvYieldCombStatPbPb2760GeV[c],
                                                                                          fitFunctionShiftingX,
                                                                                          0, grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN(),"Pi0");

          grapharrayInvYieldCombSysPbPb2760GeV[c]         = ApplyXshiftIndividualSpectra (grapharrayInvYieldCombTotPbPb2760GeV[c],
                                                                                          grapharrayInvYieldCombSysPbPb2760GeV[c],
                                                                                          fitFunctionShiftingX,
                                                                                          0, grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN(),"Pi0");

          for(Int_t m=0; m<3; m++){
              cout << "cent " << c << " method " << m << endl;
//               if(m==0 && (c==3 || c==6 ||c==7) ) rangePi0HighComb[m] = 18;
              if(m==1 && c==5) cout << "not for this" << endl;
              else if(m==2 && !(c==2 || c==5)) cout << "not for this" << endl;
              else {
                grapharrayInvYieldStatPbPb2760GeV[c][m]         = ApplyXshiftIndividualSpectra( grapharrayInvYieldCombTotPbPb2760GeV[c],
                                                                                                grapharrayInvYieldStatPbPb2760GeV[c][m],
                                                                                                fitFunctionShiftingX,
                                                                                                rangePi0LowComb[m], rangePi0HighComb[m],"Pi0");

                grapharrayInvYieldSysPbPb2760GeV[c][m]          = ApplyXshiftIndividualSpectra( grapharrayInvYieldCombTotPbPb2760GeV[c],
                                                                                                grapharrayInvYieldSysPbPb2760GeV[c][m],
                                                                                                fitFunctionShiftingX,
                                                                                                rangePi0LowComb[m], rangePi0HighComb[m],"Pi0");
              }
          }
        } else if(meson.CompareTo("Eta")==0){

          if(!(c==3 || c==6 || c==7)){
            TF1* fitFunctionShiftingX = FitObject("tcmpt","tcmptEta","Eta",grapharrayInvYieldCombTotPbPb2760GeV[c]);
            //combined
            grapharrayInvYieldCombTotPbPb2760GeV[c]         = ApplyXshift(grapharrayInvYieldCombTotPbPb2760GeV[c], fitFunctionShiftingX,"Eta");
            grapharrayInvYieldCombStatPbPb2760GeV[c]        = ApplyXshiftIndividualSpectra (grapharrayInvYieldCombTotPbPb2760GeV[c],
                                                                                            grapharrayInvYieldCombStatPbPb2760GeV[c],
                                                                                            fitFunctionShiftingX,
                                                                                            0, grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN(),"Eta");

            grapharrayInvYieldCombSysPbPb2760GeV[c]         = ApplyXshiftIndividualSpectra (grapharrayInvYieldCombTotPbPb2760GeV[c],
                                                                                            grapharrayInvYieldCombSysPbPb2760GeV[c],
                                                                                            fitFunctionShiftingX,
                                                                                            0, grapharrayInvYieldCombTotPbPb2760GeV[c]->GetN(),"Eta");

            for(Int_t m=0; m<3; m++){
              if(!m==1){
                  cout << "cent " << c << " method " << m << endl;
    //               if(m==0 && (c==3 || c==6 ||c==7) ) rangeEtaHighComb[m] = 18;
                  if(m==1 && c==5) cout << "not for this" << endl;
                  else if(m==2 && !(c==2 || c==5)) cout << "not for this" << endl;
                  else {
                    grapharrayInvYieldStatPbPb2760GeV[c][m]         = ApplyXshiftIndividualSpectra( grapharrayInvYieldCombTotPbPb2760GeV[c],
                                                                                                    grapharrayInvYieldStatPbPb2760GeV[c][m],
                                                                                                    fitFunctionShiftingX,
                                                                                                    rangeEtaLowComb[m], rangeEtaHighComb[m],"Eta");

                    grapharrayInvYieldSysPbPb2760GeV[c][m]          = ApplyXshiftIndividualSpectra( grapharrayInvYieldCombTotPbPb2760GeV[c],
                                                                                                    grapharrayInvYieldSysPbPb2760GeV[c][m],
                                                                                                    fitFunctionShiftingX,
                                                                                                    rangeEtaLowComb[m], rangeEtaHighComb[m],"Eta");
                  }
              }
            }
          }
        }
      }

      TCanvas* canvasDummy3 = new TCanvas("canvasDummy3","",200,10,1350,1350*1.15);  // gives the page size
      DrawGammaCanvasSettings( canvasDummy3, 0.16, 0.02, 0.02, 0.09);
      canvasDummy3->SetLogx();
      canvasDummy3->SetLogy();

      TH2F * histo2DDummy3 = new TH2F("histo2DDummy3","histo2DDummy3",11000,minPtRange-0.2,maxPtRange,1000,minYaxisYields,maxYaxisYields);
          SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
      histo2DDummy3->GetXaxis()->SetMoreLogLabels();
      histo2DDummy3->GetXaxis()->SetLabelOffset(-0.01);
      histo2DDummy3->Draw("copy");

        if(grapharrayInvYieldCombStatPbPb2760GeVUnshifted[c]){
          DrawGammaSetMarkerTGraphAsym(grapharrayInvYieldCombStatPbPb2760GeVUnshifted[c], 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
          grapharrayInvYieldCombStatPbPb2760GeVUnshifted[c]->Draw("pEsame");
          DrawGammaSetMarkerTGraphAsym(grapharrayInvYieldCombStatPbPb2760GeV[c], 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
          grapharrayInvYieldCombStatPbPb2760GeV[c]->Draw("pEsame");
        }

        TLegend* legendYdummy = new TLegend(0.55,0.8,0.85,0.95);
        legendYdummy->SetFillColor(0);
        legendYdummy->SetLineColor(0);
        legendYdummy->SetTextFont(42);
        legendYdummy->SetTextSize(FontSize);
        legendYdummy->AddEntry(grapharrayInvYieldCombStatPbPb2760GeVUnshifted[c],"combined unshifted","p");
        legendYdummy->AddEntry(grapharrayInvYieldCombStatPbPb2760GeV[c],"combined shifted","p");
        legendYdummy->Draw();

      canvasDummy3->Update();
      canvasDummy3->Print(Form("%s/%s_ComparisonShifted_%d.%s",outputDir.Data(),meson.Data(),c,suffix.Data()));


      for(Int_t m = 0; m<3; m++){

        TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1350,1350*1.15);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2, 0.16, 0.02, 0.02, 0.09);
        canvasDummy2->SetLogx();
        canvasDummy2->SetLogy();

        TH2F * histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",11000,minPtRange-0.2,maxPtRange,1000,minYaxisYields,maxYaxisYields);
            SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
        histo2DDummy2->GetXaxis()->SetMoreLogLabels();
        histo2DDummy2->GetXaxis()->SetLabelOffset(-0.01);
        histo2DDummy2->Draw("copy");

        TLegend* legendXdummy = new TLegend(0.55,0.85,0.85,0.95);
        legendXdummy->SetFillColor(0);
        legendXdummy->SetLineColor(0);
        legendXdummy->SetTextFont(42);
        legendXdummy->SetTextSize(FontSize);
        if(grapharrayInvYieldStatPbPb2760GeVUnshifted[c][m]){
          DrawGammaSetMarkerTGraphAsym(grapharrayInvYieldStatPbPb2760GeVUnshifted[c][m], 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
          grapharrayInvYieldStatPbPb2760GeVUnshifted[c][m]->Draw("pEsame");
          DrawGammaSetMarkerTGraphAsym(grapharrayInvYieldStatPbPb2760GeV[c][m], 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
          grapharrayInvYieldStatPbPb2760GeV[c][m]->Draw("pEsame");
          legendXdummy->AddEntry(grapharrayInvYieldStatPbPb2760GeVUnshifted[c][m],"combined unshifted","lp");
          legendXdummy->AddEntry(grapharrayInvYieldStatPbPb2760GeV[c][m],"combined shifted","lp");
        }
        legendXdummy->Draw();
        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/%s_ComparisonXShifted_%d_%d.%s",outputDir.Data(),meson.Data(),c,m,suffix.Data()));
      }
    }


    TGraphAsymmErrors* grapharrayEtaToPi0CombTotPbPb2760GeV[nCent];
    TGraphAsymmErrors* grapharrayEtaToPi0CombStatPbPb2760GeV[nCent];
    TGraphAsymmErrors* grapharrayEtaToPi0CombSysPbPb2760GeV[nCent];
    for(Int_t c=0; c<nCent; c++){
        grapharrayEtaToPi0CombTotPbPb2760GeV[c] = NULL;
        grapharrayEtaToPi0CombStatPbPb2760GeV[c] = NULL;
        grapharrayEtaToPi0CombSysPbPb2760GeV[c] = NULL;
    }

    if(meson.CompareTo("Eta")==0){
        // Declaration & calculation of combined spectrum
      for(Int_t c=0; c<nCent; c++){
        if(c==0 || c==1 || c==2 || c==4 || c==5){
          cout << "***********************************************" << endl;
          cout << "**************** centrality " << c << " *****************" << endl;
          cout << "***********************************************" << endl;
          grapharrayEtaToPi0CombTotPbPb2760GeV[c] = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEtatoPi0LHC11h[c], sysErrorCollectionEtatoPi0LHC11h[c],
                                                                                      xPtLimitsEta, 13, offSetsCombEta, offSetsCombEtaSys,
                                                                                      grapharrayEtaToPi0CombStatPbPb2760GeV[c], grapharrayEtaToPi0CombSysPbPb2760GeV[c],
                                                                                      fileNameOutputWeighting,1 );
          grapharrayEtaToPi0CombTotPbPb2760GeV[c]->RemovePoint(0);
          grapharrayEtaToPi0CombStatPbPb2760GeV[c]->RemovePoint(0);
          grapharrayEtaToPi0CombSysPbPb2760GeV[c]->RemovePoint(0);

          if(!(c==2 || c==5)){
              grapharrayEtaToPi0CombTotPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombTotPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombTotPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombTotPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombTotPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombTotPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombTotPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombTotPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombStatPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombStatPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombStatPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombStatPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombStatPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombStatPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombStatPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombStatPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombSysPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombSysPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombSysPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombSysPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombSysPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombSysPbPb2760GeV[c]->GetN()-1);
              grapharrayEtaToPi0CombSysPbPb2760GeV[c]->RemovePoint(grapharrayEtaToPi0CombSysPbPb2760GeV[c]->GetN()-1);
          }
        }
      }
    }

    TFile *fCombResults = new TFile(Form("%s/NeutralMesonInputPbPb2760GeV_%s.root", outputDir.Data(),dateForOutput.Data()), "UPDATE");
    fCombResults->cd();

    for(Int_t c = 0; c<nCent; c++){
      cout << "writing for centrality " << cent[c].Data() << endl;
      if(!(grapharrayInvYieldCombStatPbPb2760GeV[c]->GetY()==0x0))grapharrayInvYieldCombStatPbPb2760GeV[c]->Write(Form("graphInvYield%sCombPbPb2760GeVStatErr_%s",meson.Data(),cent[c].Data()));
      if(!(grapharrayInvYieldCombSysPbPb2760GeV[c]->GetY()==0x0))grapharrayInvYieldCombSysPbPb2760GeV[c]->Write(Form("graphInvYield%sCombPbPb2760GeVSysErr_%s",meson.Data(),cent[c].Data()));

      for(Int_t m = 0; m<3; m++){
          if(!(grapharrayInvYieldStatPbPb2760GeV[c][m]->GetY()==0x0))grapharrayInvYieldStatPbPb2760GeV[c][m]->Write(Form("graphInvYield%s%sPbPb2760GeVStatErr_%s",meson.Data(),meth[m].Data(),cent[c].Data()));
          if(!(grapharrayInvYieldSysPbPb2760GeV[c][m]->GetY()==0x0))grapharrayInvYieldSysPbPb2760GeV[c][m]->Write(Form("graphInvYield%s%sPbPb2760GeVSysErr_%s",meson.Data(),meth[m].Data(),cent[c].Data()));
      }

      if(meson.CompareTo("Eta")==0){
        if(!(c==3 || c==6 || c==7)){
          if(!(grapharrayEtaToPi0CombStatPbPb2760GeV[c]->GetY()==0x0))grapharrayEtaToPi0CombStatPbPb2760GeV[c]->Write(Form("graphEtaToPi0CombPbPb2760GeVStatErr_%s",cent[c].Data()));
          if(!(grapharrayEtaToPi0CombSysPbPb2760GeV[c]->GetY()==0x0))grapharrayEtaToPi0CombSysPbPb2760GeV[c]->Write(Form("graphEtaToPi0CombPbPb2760GeVSysErr_%s",cent[c].Data()));
          for(Int_t m = 0; m<3; m++){
            if(!(m==1)){
              cout << "writing for method " << meth[m].Data() << endl;
              if(!(grapharrayEtaToPi0Stat2760GeV[c][m]->GetY()==0x0))grapharrayEtaToPi0Stat2760GeV[c][m]->Write(Form("graphEtaToPi0%sPbPb2760GeVStatErr_%s",meth[m].Data(),cent[c].Data()));
              if(!(grapharrayEtaToPi0Sys2760GeV[c][m]->GetY()==0x0))grapharrayEtaToPi0Sys2760GeV[c][m]->Write(Form("graphEtaToPi0%sPbPb2760GeVSysErr_%s",meth[m].Data(),cent[c].Data()));
            }
          }
        }
      }
    }
    fCombResults->Close();

}

