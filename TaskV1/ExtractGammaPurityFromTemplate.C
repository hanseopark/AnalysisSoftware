// provided by Gamma Conversion Group, $ALICE_PHYSCIS/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

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
#include <string>
#include "TGaxis.h"
#include "TFractionFitter.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "THStack.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h"
#include "TEllipse.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/PlottingMeson.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "ExtractGammaPurityFromTemplate.h"

void ExtractGammaPurityFromTemplate( TString meson              = "",
                                     TString fileData           = "",
                                     TString fileMC             = "",
                                     TString cutSelection       = "",
                                     TString suffix             = "",
                                     TString isMC               = "",
                                     TString option             = "",
                                     TString directphotonPlots  = "",
                                     TString periodData         = "",
                                     TString periodMC           = "",
                                     Int_t   numberOfBins       = 30,
                                     Int_t   moreBG             = 2,
                                     Bool_t  addSig             = 0,
                                     Int_t   mode               = 12)
{
    if (fileMC.CompareTo("") == 0) { cout << "no MC file specified!" << endl; return; }

    //********************************* Catch modes which are not supported ****************************
    if (mode == 9) {
        cout << "ERROR: this mode is not supported anymore" << endl;
        return;
    } else if ( mode == 1 ){
        cout << "ERROR: you can't run the photon extraction in the Dalitz mode" << endl;
        return;
    } else if ( mode == 2 ||  mode == 3 ){
        cout << "WARNING: you are running in hybrid mode the software is still under construction for this one" << endl;
        fEnablePCM                              = 1;
        fEnableCalo                             = 1;
    } else if ( mode == 4 ||  mode == 5 ){
        cout << "WARNING: you are running in calo mode the software is still under construction for this one" << endl;
        fEnableCalo                             = 1;
    } else if ( mode == 0){
        fEnablePCM                              = 1;
    }

    //************************************ Set general style settings **********************************
    StyleSettingsThesis();
    SetPlotStyle();

    //************************************ Define Output directory *************************************
    fOutputDir = Form("%s/%s/%s/ExtractGammaPurityFromTemplate",cutSelection.Data(),option.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+fOutputDir);
    TString fAllPlots = Form("%s/%s/%s/ExtractGammaPurityFromTemplate/AllContributions",cutSelection.Data(),option.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+fAllPlots);
    TString fMainPlots = Form("%s/%s/%s/ExtractGammaPurityFromTemplate/MainContributions",cutSelection.Data(),option.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+fMainPlots);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    //************************************ Set global variables ****************************************
    fDate                                       = ReturnDateString();
    fDirectPhoton                               = directphotonPlots;
    fEnergyFlag                                 = option;
    fCollisionSystem                            = ReturnFullCollisionsSystem(fEnergyFlag);
    fPrefix                                     = meson;
    fSuffix                                     = suffix;
    fMeson                                      = meson;
    fMode                                       = mode;
    cout << "Pictures are saved as " << suffix.Data() << endl;

    //************************************ Detect correct folder name **********************************
    TString nameMainDir                         = "GammaConvV1_70_v2";

    //************************************ Separate cutstrings *****************************************
    fCutSelection                               = cutSelection;
    fCutSelectionRead                           = cutSelection;
    ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutNumber, fGammaCutNumber, fClusterCutNumber, fElectronCutNumber, fMesonCutNumber, fMode);
    fEventCutSelection                          = fEventCutNumber.Data();
    fGammaCutSelection                          = fGammaCutNumber.Data();
    TString centralityString                    = GetCentralityString(fEventCutSelection);
    TString InfoSystem                          = "";
    if (fCollisionSystem.Contains("pp"))
        InfoSystem                              = Form("%s",fCollisionSystem.Data());
    else
        InfoSystem                              = Form("%s %s",centralityString.Data(),fCollisionSystem.Data());

    //************************************ Specification Data/MC ***************************************
    if(isMC.CompareTo("kTRUE") == 0 || fileData.CompareTo("") == 0){
        fIsMC                                   = 1;
        fPrefix2                                = "MC";
    } else {
        fIsMC                                   = 0;
        fPrefix2                                = "data";
    }

    //************************************ Load binning for spectrum ***********************************
    Initialize(fMeson, fEnergyFlag, numberOfBins, fMode, addSig);

    Bool_t      kAllKappaPlots                  = 0;
    if(fIsMC)   kAllKappaPlots                  = 1;
    Int_t       nTemplate;
    if(moreBG==1)  nTemplate = 6;
    else        nTemplate = 4;

    Double_t value[nTemplate];
    Double_t error[nTemplate];
    Double_t valuesElEl[fNBinsPt];
    Double_t valuesRest[fNBinsPt];
    Double_t valuesPiPi[fNBinsPt];
    Double_t valuesElPi[fNBinsPt];
    Double_t valuesPartialRest[fNBinsPt];
    Double_t valuesElPK[fNBinsPt];
    Double_t valuesPiPK[fNBinsPt];
    Double_t ratiofitfractionElEl[fNBinsPt];
    Double_t ratiofitfractionRest[fNBinsPt];
    Double_t ratiofitfractionPiPi[fNBinsPt];
    Double_t ratiofitfractionElPi[fNBinsPt];
    Double_t ratiofitfractionPartialRest[fNBinsPt];
    Double_t ratiofitfractionElPK[fNBinsPt];
    Double_t ratiofitfractionPiPK[fNBinsPt];
    Double_t errorsElEl[fNBinsPt];
    Double_t errorsRest[fNBinsPt];
    Double_t errorsPiPi[fNBinsPt];
    Double_t errorsElPi[fNBinsPt];
    Double_t errorsPartialRest[fNBinsPt];
    Double_t errorsElPK[fNBinsPt];
    Double_t errorsPiPK[fNBinsPt];
    Double_t fracElEl[fNBinsPt];
    Double_t fracRest[fNBinsPt];
    Double_t fracPiPi[fNBinsPt];
    Double_t fracElPi[fNBinsPt];
    Double_t fracPartialRest[fNBinsPt];
    Double_t fracElPK[fNBinsPt];
    Double_t fracPiPK[fNBinsPt];

    //************************************ Read files **************************************************
    TFile fMC(fileMC.Data());
    TList *TopDirMC                             = (TList*)fMC.Get(nameMainDir.Data());
    if (!TopDirMC)                              { cout << "ERROR: TopDirMC not Found" << endl; return; }
    TList* HistosGammaConversionMC              = (TList*)TopDirMC->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if (!HistosGammaConversionMC)               { cout << "ERROR: folder with Cutnumber " <<   fCutSelectionRead.Data() << " not contained in file " << fileMC.Data() << endl; return; }
    TList* ESDContainerMC                       = (TList*)HistosGammaConversionMC->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));

    TFile fData(fileData.Data());
    TList *TopDirData                           = (TList*)fData.Get(nameMainDir.Data());
    if (!TopDirData)                            { cout << "ERROR: TopDirData not Found" << endl; return; }
    TList* HistosGammaConversionData            = (TList*)TopDirData->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if (!HistosGammaConversionData)             { cout << "ERROR: folder with Cutnumber " <<   fCutSelectionRead.Data() << " not contained in file " << fileData.Data() << endl; return; }
    TList* ESDContainerData                     = (TList*)HistosGammaConversionData->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));

    TH2F*       hKappaTPCPtAfterCut             = NULL;
    if(fIsMC)   hKappaTPCPtAfterCut             = (TH2F*)ESDContainerMC->FindObject("KappaTPC_Pt_after");
    else        hKappaTPCPtAfterCut             = (TH2F*)ESDContainerData->FindObject("KappaTPC_Pt_after");

    TH2F* hKappaTPCPtElEl                       = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp0_ee");
    TH2F* hKappaTPCPtPiPi                       = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp1_pipi");
    TH2F* hKappaTPCPtElPi                       = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp2_pie");
    TH2F* hKappaTPCPtPiK                        = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp3_piK");
    TH2F* hKappaTPCPtPiP                        = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp4_pip");
    TH2F* hKappaTPCPtElK                        = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp5_eK");
    TH2F* hKappaTPCPtElP                        = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp6_ep");
    TH2F* hKappaTPCPtKK                         = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp7_KK");
    TH2F* hKappaTPCPtHad                        = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp8_had");
    TH2F* hKappaTPCPtRest                       = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp9_rem4");
    TH2F* hKappaTPCPtLeftoverRest               = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp10_rem10");

    hKappaTPCAfterCut                           = new TH1D*[fNBinsPt];
    hKappaTPCElEl                               = new TH1D*[fNBinsPt];
    hKappaTPCPiPi                               = new TH1D*[fNBinsPt];
    hKappaTPCElPi                               = new TH1D*[fNBinsPt];
    hKappaTPCPiK                                = new TH1D*[fNBinsPt];
    hKappaTPCPiP                                = new TH1D*[fNBinsPt];
    hKappaTPCElK                                = new TH1D*[fNBinsPt];
    hKappaTPCElP                                = new TH1D*[fNBinsPt];
    hKappaTPCKK                                 = new TH1D*[fNBinsPt];
    hKappaTPCHad                                = new TH1D*[fNBinsPt];
    hKappaTPCRest                               = new TH1D*[fNBinsPt];
    hKappaTPCLeftoverRest                       = new TH1D*[fNBinsPt];
    hKappaTPCPartialRest                        = new TH1D*[fNBinsPt];
    hKappaTPCElPK                               = new TH1D*[fNBinsPt];
    hKappaTPCPiPK                               = new TH1D*[fNBinsPt];

    hKappaTPCSum                                = new TH1D*[fNBinsPt];
    hKappaTPCTotalSum                           = new TH1D*[fNBinsPt];
    hKappaTPCSumBkg                             = new TH1D*[fNBinsPt];

    hTemplateElEl                               = new TH1D*[fNBinsPt];
    hTemplateRest                               = new TH1D*[fNBinsPt];
    hTemplatePiPi                               = new TH1D*[fNBinsPt];
    hTemplateElPi                               = new TH1D*[fNBinsPt];
    hTemplateSum                                = new TH1D*[fNBinsPt];
    hTemplateElPK                               = new TH1D*[fNBinsPt];
    hTemplatePiPK                               = new TH1D*[fNBinsPt];
    hTemplatePartialRest                        = new TH1D*[fNBinsPt];

    Double_t kappaRangeElEl1[]                  = { -3.0,   5.0}; //narrow
    Double_t kappaRangeElEl2[]                  = { -5.0,  10.0}; //wide
    Double_t kappaRangeElEl3[]                  = { -3.0,  10.0}; //asymm
    Double_t kappaRangePiPi[]                   = {-20.0, -13.0};
    Double_t kappaRangeElPi[]                   = {-11.0,  -6.0};
    Double_t kappaRangeTail[]                   = { 11.0,  20.0};

    hValuesElEl                                 = new TH1D("hValuesElEl", "", fNBinsPt, fBinsPt);
    hValuesElPi                                 = new TH1D("hValuesElPi", "", fNBinsPt, fBinsPt);
    hValuesPiPi                                 = new TH1D("hValuesPiPi", "", fNBinsPt, fBinsPt);
    hValuesRest                                 = new TH1D("hValuesRest", "", fNBinsPt, fBinsPt);
    hValuesElPK                                 = new TH1D("hValuesElPK", "", fNBinsPt, fBinsPt);
    hValuesPiPK                                 = new TH1D("hValuesPiPK", "", fNBinsPt, fBinsPt);
    hValuesPartialRest                          = new TH1D("hValuesPartialRest", "", fNBinsPt, fBinsPt);
    hFractionElEl                               = new TH1D("hFractionElEl", "", fNBinsPt, fBinsPt);
    hFractionElPi                               = new TH1D("hFractionElPi", "", fNBinsPt, fBinsPt);
    hFractionPiPi                               = new TH1D("hFractionPiPi", "", fNBinsPt, fBinsPt);
    hFractionRest                               = new TH1D("hFractionRest", "", fNBinsPt, fBinsPt);
    hFractionElPK                               = new TH1D("hFractionElPK", "", fNBinsPt, fBinsPt);
    hFractionPiPK                               = new TH1D("hFractionPiPK", "", fNBinsPt, fBinsPt);
    hFractionPartialRest                        = new TH1D("hFractionPartialRest", "", fNBinsPt, fBinsPt);

    //************************************ Purity histogram ********************************************
    hSignalPurity1                              = new TH1D("hSignalPurity1", "", fNBinsPt, fBinsPt);    // narrow
    hSignalPurity2                              = new TH1D("hSignalPurity2", "", fNBinsPt, fBinsPt);    // wide
    hSignalPurity3                              = new TH1D("hSignalPurity3", "", fNBinsPt, fBinsPt);    // asymm
    hBackgroundPiPi1                            = new TH1D("hBackgroundPiPi1", "", fNBinsPt, fBinsPt);
    hBackgroundPiPi2                            = new TH1D("hBackgroundPiPi2", "", fNBinsPt, fBinsPt);
    hBackgroundPiPi3                            = new TH1D("hBackgroundPiPi3", "", fNBinsPt, fBinsPt);
    hBackgroundElPi1                            = new TH1D("hBackgroundElPi1", "", fNBinsPt, fBinsPt);
    hBackgroundElPi2                            = new TH1D("hBackgroundElPi2", "", fNBinsPt, fBinsPt);
    hBackgroundElPi3                            = new TH1D("hBackgroundElPi3", "", fNBinsPt, fBinsPt);
    hBackgroundRest1                            = new TH1D("hBackgroundRest1", "", fNBinsPt, fBinsPt);
    hBackgroundRest2                            = new TH1D("hBackgroundRest2", "", fNBinsPt, fBinsPt);
    hBackgroundRest3                            = new TH1D("hBackgroundRest3", "", fNBinsPt, fBinsPt);
    hBackgroundElPK1                            = new TH1D("hBackgroundElPK1", "", fNBinsPt, fBinsPt);
    hBackgroundElPK2                            = new TH1D("hBackgroundElPK2", "", fNBinsPt, fBinsPt);
    hBackgroundElPK3                            = new TH1D("hBackgroundElPK3", "", fNBinsPt, fBinsPt);
    hBackgroundPiPK1                            = new TH1D("hBackgroundPiPK1", "", fNBinsPt, fBinsPt);
    hBackgroundPiPK2                            = new TH1D("hBackgroundPiPK2", "", fNBinsPt, fBinsPt);
    hBackgroundPiPK3                            = new TH1D("hBackgroundPiPK3", "", fNBinsPt, fBinsPt);
    hBackgroundPartialRest1                     = new TH1D("hBackgroundPartialRest1", "", fNBinsPt, fBinsPt);
    hBackgroundPartialRest2                     = new TH1D("hBackgroundPartialRest2", "", fNBinsPt, fBinsPt);
    hBackgroundPartialRest3                     = new TH1D("hBackgroundPartialRest3", "", fNBinsPt, fBinsPt);

    Int_t errorCounter                          = 0;
    Int_t errorBin[fNBinsPt];
    Int_t errorStatus[fNBinsPt];
    Bool_t      kUseFitConstraints              = 1;

    Double_t ptconstraintHigh[fNBinsPt]; //{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
////                                0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20
    Double_t ptconstraintHigh0[] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10};
    Double_t ptconstraintHigh1[] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10,15,20,20,20,30,40,40};
    Double_t ptconstraintHigh2[] = {10,10,10,10,10,10,10,10,10,10,10,10,15,20,20,20,20,20,20,20,20};
    Double_t ptconstraintHigh3[] = {10,10,10,10,10,10,10,10,10,10,10,10,15,30,15,20, 5, 5,10,10,10};

    for(Int_t i = fStartPtBin; i<fNBinsPt; i++){
        errorBin[i]                             = 0;
        errorStatus[i]                          = 0;

        Int_t startBin                          = hKappaTPCPtAfterCut->GetYaxis()->FindBin(fBinsPt[i]);
        Int_t endBin                            = hKappaTPCPtAfterCut->GetYaxis()->FindBin(fBinsPt[i+1]);

        cout << "Projecting histograms" << endl;
        hKappaTPCAfterCut[i]                    = (TH1D*)hKappaTPCPtAfterCut->ProjectionX(Form("KappaTPC_AfterCut_Pt_%.2d",i),startBin,endBin,"e");

        hKappaTPCElEl[i]                        = (TH1D*)hKappaTPCPtElEl->ProjectionX(Form("KappaTPC_ElEl_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCElPi[i]                        = (TH1D*)hKappaTPCPtElPi->ProjectionX(Form("KappaTPC_ElPi_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCPiPi[i]                        = (TH1D*)hKappaTPCPtPiPi->ProjectionX(Form("KappaTPC_PiPi_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCRest[i]                        = (TH1D*)hKappaTPCPtRest->ProjectionX(Form("KappaTPC_Rest_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCElK[i]                         = (TH1D*)hKappaTPCPtElK->ProjectionX(Form("KappaTPC_ElK_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCElP[i]                         = (TH1D*)hKappaTPCPtElP->ProjectionX(Form("KappaTPC_ElP_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCPiP[i]                         = (TH1D*)hKappaTPCPtPiP->ProjectionX(Form("KappaTPC_PiP_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCPiK[i]                         = (TH1D*)hKappaTPCPtPiK->ProjectionX(Form("KappaTPC_PiK_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCKK[i]                          = (TH1D*)hKappaTPCPtKK->ProjectionX(Form("KappaTPC_KK_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCHad[i]                         = (TH1D*)hKappaTPCPtHad->ProjectionX(Form("KappaTPC_Had_Pt_%.2d",i),startBin,endBin,"e");
        hKappaTPCLeftoverRest[i]                = (TH1D*)hKappaTPCPtLeftoverRest->ProjectionX(Form("KappaTPC_LeftoverRest_Pt_%.2d",i),startBin,endBin,"e");

        hKappaTPCSum[i]                         = (TH1D*)hKappaTPCElEl[i]->Clone(Form("hKappaTPCSum%d",i));
        hKappaTPCSum[i]->Sumw2();
        if(moreBG==1){

          hKappaTPCElPK[i]               = (TH1D*)hKappaTPCElP[i]->Clone(Form("hKappaTPCElPK%d",i));
          hKappaTPCElPK[i]->Sumw2();
          hKappaTPCElPK[i]->Add(hKappaTPCElK[i]);
          hKappaTPCPiPK[i]               = (TH1D*)hKappaTPCPiP[i]->Clone(Form("hKappaTPCPiPK%d",i));
          hKappaTPCPiPK[i]->Sumw2();
          hKappaTPCPiPK[i]->Add(hKappaTPCPiK[i]);

          hKappaTPCPartialRest[i]               = (TH1D*)hKappaTPCKK[i]->Clone(Form("hKappaTPCPartialRest%d",i));
          hKappaTPCPartialRest[i]->Sumw2();
          hKappaTPCPartialRest[i]->Add(hKappaTPCHad[i]);
          hKappaTPCPartialRest[i]->Add(hKappaTPCLeftoverRest[i]);

          hKappaTPCSum[i]->Add(hKappaTPCPiPi[i]);
          hKappaTPCSum[i]->Add(hKappaTPCElPi[i]);
          hKappaTPCSum[i]->Add(hKappaTPCElPK[i]);
          hKappaTPCSum[i]->Add(hKappaTPCPiPK[i]);
          hKappaTPCSum[i]->Add(hKappaTPCPartialRest[i]);

        } else if(moreBG==2){

          hKappaTPCElPK[i]               = (TH1D*)hKappaTPCElPi[i]->Clone(Form("hKappaTPCElPK%d",i));
          hKappaTPCElPK[i]->Sumw2();
          hKappaTPCElPK[i]->Add(hKappaTPCElK[i]);
          hKappaTPCElPK[i]->Add(hKappaTPCElP[i]);
          hKappaTPCPiPK[i]               = (TH1D*)hKappaTPCPiPi[i]->Clone(Form("hKappaTPCPiPK%d",i));
          hKappaTPCPiPK[i]->Sumw2();
          hKappaTPCPiPK[i]->Add(hKappaTPCPiK[i]);
          hKappaTPCPiPK[i]->Add(hKappaTPCPiP[i]);

          hKappaTPCPartialRest[i]               = (TH1D*)hKappaTPCKK[i]->Clone(Form("hKappaTPCPartialRest%d",i));
          hKappaTPCPartialRest[i]->Sumw2();
          hKappaTPCPartialRest[i]->Add(hKappaTPCHad[i]);
          hKappaTPCPartialRest[i]->Add(hKappaTPCLeftoverRest[i]);

          hKappaTPCSum[i]                         = (TH1D*)hKappaTPCElEl[i]->Clone(Form("hKappaTPCSum%d",i));
          hKappaTPCSum[i]->Sumw2();
          hKappaTPCSum[i]->Add(hKappaTPCElPK[i]);
          hKappaTPCSum[i]->Add(hKappaTPCPiPK[i]);
          hKappaTPCSum[i]->Add(hKappaTPCPartialRest[i]);

        } else {
          hKappaTPCSum[i]                         = (TH1D*)hKappaTPCElEl[i]->Clone(Form("hKappaTPCSum%d",i));
          hKappaTPCSum[i]->Sumw2();
          hKappaTPCSum[i]->Add(hKappaTPCPiPi[i]);
          hKappaTPCSum[i]->Add(hKappaTPCElPi[i]);
          hKappaTPCSum[i]->Add(hKappaTPCRest[i]);
        }

        hKappaTPCTotalSum[i]                    = (TH1D*)hKappaTPCElEl[i]->Clone(Form("hKappaTPCTotalSum%d",i));
        hKappaTPCTotalSum[i]->Sumw2();
        hKappaTPCTotalSum[i]->Add(hKappaTPCPiPi[i]);
        hKappaTPCTotalSum[i]->Add(hKappaTPCElPi[i]);
        hKappaTPCTotalSum[i]->Add(hKappaTPCPiK[i]);
        hKappaTPCTotalSum[i]->Add(hKappaTPCElK[i]);
        hKappaTPCTotalSum[i]->Add(hKappaTPCPiP[i]);
        hKappaTPCTotalSum[i]->Add(hKappaTPCElP[i]);
        hKappaTPCTotalSum[i]->Add(hKappaTPCKK[i]);
        hKappaTPCTotalSum[i]->Add(hKappaTPCHad[i]);
        hKappaTPCTotalSum[i]->Add(hKappaTPCLeftoverRest[i]);

        Double_t NTotalMC;
        if(moreBG==1) NTotalMC = hKappaTPCElEl[i]->GetEntries() + hKappaTPCPiPi[i]->GetEntries() + hKappaTPCElPi[i]->GetEntries() + hKappaTPCElPK[i]->GetEntries() + hKappaTPCPiPK[i]->GetEntries() + hKappaTPCPartialRest[i]->GetEntries();
        else if(moreBG==2) NTotalMC = hKappaTPCElEl[i]->GetEntries() + hKappaTPCElPK[i]->GetEntries() + hKappaTPCPiPK[i]->GetEntries() + hKappaTPCPartialRest[i]->GetEntries();
        else NTotalMC = hKappaTPCElEl[i]->GetEntries() + hKappaTPCPiPi[i]->GetEntries() + hKappaTPCElPi[i]->GetEntries() + hKappaTPCRest[i]->GetEntries();

        Double_t NEntriesAfterCut               = hKappaTPCAfterCut[i]->GetEntries();

        if(fIsMC){
            if(moreBG==1){
              fracElEl[i]                         = hKappaTPCElEl[i]->GetEntries()/NEntriesAfterCut;
              fracPiPi[i]                         = hKappaTPCPiPi[i]->GetEntries()/NEntriesAfterCut;
              fracElPi[i]                         = hKappaTPCElPi[i]->GetEntries()/NEntriesAfterCut;
              fracElPK[i]                       = hKappaTPCElPK[i]->GetEntries()/NEntriesAfterCut;
              fracPiPK[i]                       = hKappaTPCPiPK[i]->GetEntries()/NEntriesAfterCut;
              fracPartialRest[i]                = hKappaTPCPartialRest[i]->GetEntries()/NEntriesAfterCut;
            } else if(moreBG==2){
              fracElEl[i]                         = hKappaTPCElEl[i]->GetEntries()/NEntriesAfterCut;
              fracElPK[i]                       = hKappaTPCElPK[i]->GetEntries()/NEntriesAfterCut;
              fracPiPK[i]                       = hKappaTPCPiPK[i]->GetEntries()/NEntriesAfterCut;
              fracPartialRest[i]                = hKappaTPCPartialRest[i]->GetEntries()/NEntriesAfterCut;
            } else {
              fracElEl[i]                         = hKappaTPCElEl[i]->GetEntries()/NEntriesAfterCut;
              fracPiPi[i]                         = hKappaTPCPiPi[i]->GetEntries()/NEntriesAfterCut;
              fracElPi[i]                         = hKappaTPCElPi[i]->GetEntries()/NEntriesAfterCut;
              fracRest[i]                         = hKappaTPCRest[i]->GetEntries()/NEntriesAfterCut;
            }
        } else {
            if(moreBG==1){
              fracElEl[i]                         = hKappaTPCElEl[i]->GetEntries()/NTotalMC;
              fracPiPi[i]                         = hKappaTPCPiPi[i]->GetEntries()/NTotalMC;
              fracElPi[i]                         = hKappaTPCElPi[i]->GetEntries()/NTotalMC;
              fracElPK[i]                       = hKappaTPCElPK[i]->GetEntries()/NTotalMC;
              fracPiPK[i]                       = hKappaTPCPiPK[i]->GetEntries()/NTotalMC;
              fracPartialRest[i]                = hKappaTPCPartialRest[i]->GetEntries()/NTotalMC;
            } else if(moreBG==2){
              fracElEl[i]                       = hKappaTPCElEl[i]->GetEntries()/NTotalMC;
              fracElPK[i]                       = hKappaTPCElPK[i]->GetEntries()/NTotalMC;
              fracPiPK[i]                       = hKappaTPCPiPK[i]->GetEntries()/NTotalMC;
              fracPartialRest[i]                = hKappaTPCPartialRest[i]->GetEntries()/NTotalMC;
            } else {
              fracElEl[i]                         = hKappaTPCElEl[i]->GetEntries()/NTotalMC;
              fracPiPi[i]                         = hKappaTPCPiPi[i]->GetEntries()/NTotalMC;
              fracElPi[i]                         = hKappaTPCElPi[i]->GetEntries()/NTotalMC;
              fracRest[i]                         = hKappaTPCRest[i]->GetEntries()/NTotalMC;
            }
        }

        //==================================================================================//
        //                       TEMPLATE FITTING and PLOTTING                              //
        //==================================================================================//


        TObjArray *kappaArray                   = NULL;
        TFractionFitter* fit                    = NULL;
        cout << "TFractionFitter initialization, then plotting" << endl;
        cout << "fitting " << Form("%2.1f < #it{p}_{T} < %2.1f (GeV/#it{c})",fBinsPt[i],fBinsPt[i+1]) << endl;
        if(moreBG==1){
            //array of all MC that build up the data
            kappaArray                   = new TObjArray(nTemplate);
            kappaArray->Add(hKappaTPCElEl[i]); //primary + secondary
            kappaArray->Add(hKappaTPCPiPi[i]); //pion+pion
            kappaArray->Add(hKappaTPCElPi[i]); //pion+electron
            kappaArray->Add(hKappaTPCElPK[i]); //electron+kaon+proton
            kappaArray->Add(hKappaTPCPiPK[i]); //pion+proton+kaon
            kappaArray->Add(hKappaTPCPartialRest[i]); //remaining

            //configure the TFractionFitter
            fit                    = new TFractionFitter(hKappaTPCAfterCut[i], kappaArray);
            if(kUseFitConstraints){
                fit->Constrain(1,fracElEl[i]*constrainLow,fracElEl[i]*constrainHigh);
                fit->Constrain(2,fracPiPi[i]*constrainLow,fracPiPi[i]*constrainHigh);
                fit->Constrain(3,fracElPi[i]*constrainLow,fracElPi[i]*constrainHigh);
                fit->Constrain(4,fracElPK[i]*constrainLow,fracElPK[i]*constrainHigh);
                fit->Constrain(5,fracPiPK[i]*constrainLow,fracPiPK[i]*constrainHigh);
                fit->Constrain(6,fracPartialRest[i]*constrainLow,fracPartialRest[i]*constrainHigh);
            }
        } else if(moreBG==2){
            //array of all MC that build up the data
            kappaArray                   = new TObjArray(nTemplate);
            kappaArray->Add(hKappaTPCElEl[i]); //primary + secondary
            kappaArray->Add(hKappaTPCElPK[i]); //electron+kaon+proton
            kappaArray->Add(hKappaTPCPiPK[i]); //pion+proton+kaon
            kappaArray->Add(hKappaTPCPartialRest[i]); //remaining

            //configure the TFractionFitter
            fit                    = new TFractionFitter(hKappaTPCAfterCut[i], kappaArray);
            if(kUseFitConstraints){
                fit->Constrain(1,fracElEl[i]*constrainLow,fracElEl[i]*ptconstraintHigh0[i]);
                fit->Constrain(2,fracElPK[i]*constrainLow,fracElPK[i]*ptconstraintHigh1[i]);
                fit->Constrain(3,fracPiPK[i]*constrainLow,fracPiPK[i]*ptconstraintHigh2[i]);
                fit->Constrain(4,fracPartialRest[i]*constrainLow,fracPartialRest[i]*ptconstraintHigh3[i]);
            }

        } else {
            //array of all MC that build up the data
            kappaArray                   = new TObjArray(nTemplate);
            kappaArray->Add(hKappaTPCElEl[i]); //primary + secondary
            kappaArray->Add(hKappaTPCPiPi[i]); //pion+pion
            kappaArray->Add(hKappaTPCElPi[i]); //pion+electron
            kappaArray->Add(hKappaTPCRest[i]); //remaining

            //configure the TFractionFitter
            fit                    = new TFractionFitter(hKappaTPCAfterCut[i], kappaArray);
            if(kUseFitConstraints){
                fit->Constrain(1,fracElEl[i]*constrainLow,fracElEl[i]*constrainHigh);
                fit->Constrain(2,fracPiPi[i]*constrainLow,fracPiPi[i]*constrainHigh);
                fit->Constrain(3,fracElPi[i]*constrainLow,fracElPi[i]*constrainHigh);
                fit->Constrain(4,fracRest[i]*constrainLow,fracRest[i]*constrainHigh);
            }
        }

        //Fit the templates
        Int_t status                            = fit->Fit();
        //cout << "fit status: " << status << endl;

        //plot 6 pad projection
        TCanvas* canvasPartKappaProj           = new TCanvas("canvasPartKappaProj","",0,0,1200,800);
        DrawGammaCanvasSettings( canvasPartKappaProj,  0.15, 0.05, 0.2, 0.15);
        if(moreBG==0 || moreBG==2) canvasPartKappaProj->Divide(3,2,0.001,0.001);
        else if(moreBG==1) canvasPartKappaProj->Divide(4,2,0.001,0.001);

        TPaveText * pave                        = new TPaveText(0.15,0.76,0.45,0.93,"NDC");
        SetStylePave(pave,42,11,0.035);
        pave->InsertText(InfoSystem.Data());
//         if(fIsMC) pave->InsertText(periodMC.Data());
//         else pave->InsertText(periodData.Data());
        pave->InsertText(SigmaStarForm.Data());
        pave->InsertText(Form("%2.1f < #it{p}_{T} < %2.1f (GeV/#it{c})",fBinsPt[i],fBinsPt[i+1]));

        TPaveText * paveElPi                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
        SetStylePave(paveElPi,62,11,0.04);
        TPaveText * pavePiPi                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
        SetStylePave(pavePiPi,62,11,0.04);

        canvasPartKappaProj->cd(1);
          TPaveText * paveSignal                  = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
          SetStylePave(paveSignal,62,11,0.04);
          paveSignal->InsertText("Signal");
          paveSignal->InsertText(Form("%f",hKappaTPCElEl[i]->GetEntries()/NTotalMC));
          HistoPlotSettings(hKappaTPCElEl[i],"K","Counts",0.00001,hKappaTPCElEl[i]->GetMaximum()*1.2,-20,20,12,12,3004);
          hKappaTPCElEl[i]->Draw("histo");
          pave->Draw();
          paveSignal->Draw();

        canvasPartKappaProj->cd(2);
          if(moreBG==0 || moreBG ==1){
              paveElPi->InsertText("#pi^{#pm} + e^{#mp}");
              paveElPi->InsertText(Form("%f",hKappaTPCElPi[i]->GetEntries()/NTotalMC));
              HistoPlotSettings(hKappaTPCElPi[i],"K","Counts",0.00001,hKappaTPCElPi[i]->GetMaximum()*1.4,-20,20,kMagenta-2,kMagenta-2,3004);
              hKappaTPCElPi[i]->Draw("histo");
              paveElPi->Draw();
          } else {
              TPaveText * paveElectrons                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
              SetStylePave(paveElectrons,62,11,0.04);
              paveElectrons->InsertText("e^{#pm} + #pi^{#mp}/K^{#mp}/p");
              paveElectrons->InsertText(Form("%f",hKappaTPCElPK[i]->GetEntries()/NTotalMC));
              HistoPlotSettings(hKappaTPCElPK[i],"K","Counts",0.00001,hKappaTPCElPK[i]->GetMaximum()*1.4,-20,20,kMagenta-2,kMagenta-2,3004);
              hKappaTPCElPK[i]->Draw("histo");
              paveElectrons->Draw();
          }

        canvasPartKappaProj->cd(3);
          if(moreBG==0 || moreBG==1){
              pavePiPi->InsertText("#pi^{#pm} + #pi^{#mp}");
              pavePiPi->InsertText(Form("%f",hKappaTPCPiPi[i]->GetEntries()/NTotalMC));
              HistoPlotSettings(hKappaTPCPiPi[i],"K","Counts",0.00001,hKappaTPCPiPi[i]->GetMaximum()*1.2,-20,20,kBlue-7,kBlue-7,3005);
              hKappaTPCPiPi[i]->Draw("histo");
              pavePiPi->Draw();
          } else {
              TPaveText * pavePions                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
              SetStylePave(pavePions,62,11,0.04);
              pavePions->InsertText("#pi^{#pm} + #pi^{#mp}/K^{#mp}/p");
              pavePions->InsertText(Form("%f",hKappaTPCPiPK[i]->GetEntries()/NTotalMC));
              HistoPlotSettings(hKappaTPCPiPK[i],"K","Counts",0.00001,hKappaTPCPiPK[i]->GetMaximum()*1.2,-20,20,kBlue-7,kBlue-7,3005);
              hKappaTPCPiPK[i]->Draw("histo");
              pavePions->Draw();
          }

        canvasPartKappaProj->cd(4);
          if(moreBG==0 || moreBG==2){
              TPaveText * paveRest                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
              SetStylePave(paveRest,62,11,0.04);
              paveRest->InsertText("Remaining");
              if(moreBG==0){
                  paveRest->InsertText(Form("%f",hKappaTPCRest[i]->GetEntries()/NTotalMC));
                  HistoPlotSettings(hKappaTPCRest[i],"K","Counts",0.00001,hKappaTPCRest[i]->GetMaximum()*1.2,-20,20,kRed+1,kRed+1,3005);
                  hKappaTPCRest[i]->Draw("histo");
              } else if(moreBG==2){
                  paveRest->InsertText(Form("%f",hKappaTPCPartialRest[i]->GetEntries()/NTotalMC));
                  HistoPlotSettings(hKappaTPCPartialRest[i],"K","Counts",0.00001,hKappaTPCPartialRest[i]->GetMaximum()*1.2,-20,20,kRed+1,kRed+1,3005);
                  hKappaTPCPartialRest[i]->Draw("histo");
              }
              paveRest->Draw();
          } else {
              TPaveText * paveElPK                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
              SetStylePave(paveElPK,62,11,0.04);
              paveElPK->InsertText("e^{#pm} + K^{#mp} + p");
              paveElPK->InsertText(Form("%f",hKappaTPCElPK[i]->GetEntries()/NTotalMC));
              HistoPlotSettings(hKappaTPCElPK[i],"K","Counts",0.00001,hKappaTPCElPK[i]->GetMaximum()*1.4,-20,20,kGreen-5,kGreen-5,3305);
              hKappaTPCElPK[i]->Draw("histo");
              paveElPK->Draw();
          }

        canvasPartKappaProj->cd(5);
          if(moreBG==0 || moreBG==2){
            canvasPartKappaProj->cd(5)->SetLogy();

              TLegend* leg                            = new TLegend(0.7,0.7,0.9,0.93);
              leg->SetTextSize(0.04);
              leg->SetBorderSize(0);

              HistoPlotSettings(hKappaTPCSum[i],"K","Counts",1,hKappaTPCSum[i]->GetMaximum()*10,-20,20,kBlack,kBlack,0);
              hKappaTPCSum[i]->SetLineWidth(2.0);
              hKappaTPCSum[i]->Draw("histo");
              leg->AddEntry(hKappaTPCSum[i],"total","lp");
              hKappaTPCElEl[i]->Draw("histo,same");
              leg->AddEntry(hKappaTPCElEl[i],"real e^{+}e^{-}","f");
              if(moreBG==0){
                hKappaTPCPiPi[i]->Draw("histo,same");
                leg->AddEntry(hKappaTPCPiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
                hKappaTPCElPi[i]->Draw("histo,same");
                leg->AddEntry(hKappaTPCElPi[i],"#pi^{#pm} + e^{#mp}","f");
                hKappaTPCRest[i]->Draw("histo,same");
                leg->AddEntry(hKappaTPCRest[i],"remaining","f");
              } else if(moreBG==2){
                hKappaTPCElPK[i]->Draw("histo,same");
                leg->AddEntry(hKappaTPCElPK[i],"e^{#pm} + #pi^{#mp}/K^{#mp}/p","f");
                hKappaTPCPiPK[i]->Draw("histo,same");
                leg->AddEntry(hKappaTPCPiPK[i],"#pi^{#pm} + #pi^{#mp}/K^{#mp}/p","f");
                hKappaTPCPartialRest[i]->Draw("histo,same");
                leg->AddEntry(hKappaTPCPartialRest[i],"remaining","f");
              }
              leg->Draw();

          } else {

              TPaveText * pavePiPK                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
              SetStylePave(pavePiPK,62,11,0.04);
              pavePiPK->InsertText("#pi^{#pm} + K^{#mp} + p");
              pavePiPK->InsertText(Form("%f",hKappaTPCPiPK[i]->GetEntries()/NTotalMC));
              HistoPlotSettings(hKappaTPCPiPK[i],"K","Counts",0.00001,hKappaTPCPiPK[i]->GetMaximum()*1.2,-20,20,kOrange+1,kOrange+1,3007);
              hKappaTPCPiPK[i]->Draw("histo");
              pavePiPK->Draw();

          }

        canvasPartKappaProj->cd(6);
          if(moreBG==0 || moreBG==2){
            canvasPartKappaProj->cd(6)->SetLogy();

              HistoPlotSettings(hKappaTPCAfterCut[i],"K","Counts",/*0.0000*/1,hKappaTPCAfterCut[i]->GetMaximum()*10/*1.2*/,-20,20,kBlue,kBlue,0);
              hKappaTPCAfterCut[i]->Draw("histo");

              TLegend* leg2                           = new TLegend(0.7,0.8,0.9,0.9);
              leg2->SetBorderSize(0);
              leg2->SetTextSize(0.04);
              leg2->AddEntry(hKappaTPCAfterCut[i],"data","lp");
              if (status == 0) {
                  TH1D *result                        = (TH1D*)fit->GetPlot();
                  result->SetLineColor(kRed);
                  result->Draw("hist,same");
                  leg2->AddEntry(result,"fit","lp");
              }
              leg2->Draw();

          } else {

              TPaveText * paveRest                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
              SetStylePave(paveRest,62,11,0.04);
              paveRest->InsertText("Remaining");
              HistoPlotSettings(hKappaTPCPartialRest[i],"K","Counts",0.00001,hKappaTPCPartialRest[i]->GetMaximum()*1.2,-20,20,kRed+1,kRed+1,3005);
              hKappaTPCPartialRest[i]->Draw("histo");
              paveRest->Draw();

            canvasPartKappaProj->cd(7);
            canvasPartKappaProj->cd(7)->SetLogy();

              TLegend* leg8                            = new TLegend(0.7,0.68,0.9,0.93);
              leg8->SetTextSize(0.04);
              leg8->SetBorderSize(0);

              HistoPlotSettings(hKappaTPCSum[i],"K","Counts",1,hKappaTPCSum[i]->GetMaximum()*20,-20,20,kBlack,kBlack,0);
              hKappaTPCSum[i]->SetLineWidth(2.0);
              hKappaTPCSum[i]->Draw("histo");
              hKappaTPCElEl[i]->Draw("histo,same");
              hKappaTPCPiPi[i]->Draw("histo,same");
              hKappaTPCElPi[i]->Draw("histo,same");
              hKappaTPCElPK[i]->Draw("histo,same");
              hKappaTPCPiPK[i]->Draw("histo,same");
              hKappaTPCPartialRest[i]->Draw("histo,same");

              leg8->AddEntry(hKappaTPCSum[i],"total","lp");
              leg8->AddEntry(hKappaTPCElEl[i],"real e^{+}e^{-}","f");
              leg8->AddEntry(hKappaTPCPartialRest[i],"remaining","f");
              leg8->AddEntry(hKappaTPCPiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
              leg8->AddEntry(hKappaTPCElPi[i],"#pi^{#pm} + e^{#mp}","f");
              leg8->AddEntry(hKappaTPCElPK[i],"e^{#pm} + K^{#mp}(p)","f");
              leg8->AddEntry(hKappaTPCPiPK[i],"#pi^{#pm} + K^{#mp}(p)","f");
              leg8->Draw();

            canvasPartKappaProj->cd(8);
            canvasPartKappaProj->cd(8)->SetLogy();

              HistoPlotSettings(hKappaTPCAfterCut[i],"K","Counts",/*0.0000*/1,hKappaTPCAfterCut[i]->GetMaximum()*10/*1.2*/,-20,20,kBlue,kBlue,0);
              hKappaTPCAfterCut[i]->Draw("histo");
              TLegend* leg2                           = new TLegend(0.7,0.8,0.9,0.9);
              leg2->SetBorderSize(0);
              leg2->SetTextSize(0.04);
              leg2->AddEntry(hKappaTPCAfterCut[i],"data","lp");
              if (status == 0) {
                  TH1D *result                        = (TH1D*)fit->GetPlot();
                  result->SetLineColor(kRed);
                  result->Draw("hist,same");
                  leg2->AddEntry(result,"fit","lp");
              }
              leg2->Draw();
          }

        canvasPartKappaProj->SaveAs(Form("%s/KappaProj_%.2d.%s",fMainPlots.Data(),i,suffix.Data()));

//         if(!moreBG){
//             TCanvas* canvasRatioFitDataProj = new TCanvas("canvasRatioFitDataProj","",800,800);
//             DrawGammaCanvasSettings( canvasRatioFitDataProj,  0.12, 0.03, 0.03, 0.1);
//
//               TH1D* hRatioDataToFitFinerBG = (TH1D*)hKappaTPCAfterCut[i]->Clone("hRatioDataToFitFinerBG");
//               if (status == 0) {
//                   TH1D *result                        = (TH1D*)fit->GetPlot();
//                   result->SetLineColor(kRed);
//                   hRatioDataToFitFinerBG->Divide(hKappaTPCAfterCut[i],result, 1., 1., "B");
//                   hRatioDataToFitFinerBG->Draw("hist");
//               }
//
//             pave->Draw();
//             canvasRatioFitDataProj->SaveAs(Form("%s/RatioFitToData_%.2d.%s",fMainPlots.Data(),i,suffix.Data()));
//         }
//
//         if(moreBG){
//             //plot 6 pad projection
//             TCanvas* canvas8PartKappaProj           = new TCanvas("canvas8PartKappaProj","",0,0,1200,800);
//             DrawGammaCanvasSettings( canvas8PartKappaProj,  0.15, 0.05, 0.2, 0.15);
//             canvas8PartKappaProj->Divide(4,2,0.001,0.001);
//
//             canvas8PartKappaProj->cd(1);
//
//               HistoPlotSettings(hKappaTPCElEl[i],"K","Counts",0.00001,hKappaTPCElEl[i]->GetMaximum()*1.2,-20,20,12,12,3004);
//               hKappaTPCElEl[i]->Draw("histo");
//               pave->Draw();
//               paveSignal->Draw();
//
//             canvas8PartKappaProj->cd(2);
//
//               HistoPlotSettings(hKappaTPCElPi[i],"K","Counts",0.00001,hKappaTPCElPi[i]->GetMaximum()*1.4,-20,20,kMagenta-2,kMagenta-2,3004);
//               hKappaTPCElPi[i]->Draw("histo");
//               paveElPi->Draw();
//
//             canvas8PartKappaProj->cd(3);
//
//               HistoPlotSettings(hKappaTPCPiPi[i],"K","Counts",0.00001,hKappaTPCPiPi[i]->GetMaximum()*1.2,-20,20,kBlue-7,kBlue-7,3005);
//               hKappaTPCPiPi[i]->Draw("histo");
//               pavePiPi->Draw();
//
//             canvas8PartKappaProj->cd(4);
//
//               TPaveText * paveElPK                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
//               SetStylePave(paveElPK,62,11,0.04);
//               paveElPK->InsertText("#pi^{#pm} + #pi^{#mp}");
//               paveElPK->InsertText(Form("%f",hKappaTPCElPK[i]->GetEntries()/NTotalMC));
//
//               HistoPlotSettings(hKappaTPCElPK[i],"K","Counts",0.00001,hKappaTPCElPK[i]->GetMaximum()*1.4,-20,20,kGreen-5,kGreen-5,3305);
//               hKappaTPCElPK[i]->Draw("histo");
//               paveElPK->Draw();
//
//             canvas8PartKappaProj->cd(5);
//
//               TPaveText * pavePiPK                    = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
//               SetStylePave(pavePiPK,62,11,0.04);
//               pavePiPK->InsertText("#pi^{#pm} + #pi^{#mp}");
//               pavePiPK->InsertText(Form("%f",hKappaTPCPiPK[i]->GetEntries()/NTotalMC));
//
//               HistoPlotSettings(hKappaTPCPiPK[i],"K","Counts",0.00001,hKappaTPCPiPK[i]->GetMaximum()*1.2,-20,20,kOrange+1,kOrange+1,3007);
//               hKappaTPCPiPK[i]->Draw("histo");
//               pavePiPK->Draw();
//
//             canvas8PartKappaProj->cd(6);
//
//               HistoPlotSettings(hKappaTPCPartialRest[i],"K","Counts",0.00001,hKappaTPCPartialRest[i]->GetMaximum()*1.2,-20,20,kRed+1,kRed+1,3005);
//               hKappaTPCPartialRest[i]->Draw("histo");
//               paveRest->Draw();
//
//             canvas8PartKappaProj->cd(7);
//             canvas8PartKappaProj->cd(7)->SetLogy();
//
//               TLegend* leg8                            = new TLegend(0.7,0.68,0.9,0.93);
//               leg8->SetTextSize(0.04);
//               leg8->SetBorderSize(0);
//
//                 HistoPlotSettings(hKappaTPCSum[i],"K","Counts",1,hKappaTPCSum[i]->GetMaximum()*10,-20,20,kBlack,kBlack,0);
//                 hKappaTPCSum[i]->SetLineWidth(2.0);
//
//                 hKappaTPCSum[i]->Draw("histo");
//                 hKappaTPCElEl[i]->Draw("histo,same");
//                 hKappaTPCPiPi[i]->Draw("histo,same");
//                 hKappaTPCElPi[i]->Draw("histo,same");
//                 hKappaTPCElPK[i]->Draw("histo,same");
//                 hKappaTPCPiPK[i]->Draw("histo,same");
//                 hKappaTPCPartialRest[i]->Draw("histo,same");
//
//                 leg8->AddEntry(hKappaTPCSum[i],"total","lp");
//                 leg8->AddEntry(hKappaTPCElEl[i],"real e^{+}e^{-}","f");
//                 leg8->AddEntry(hKappaTPCPartialRest[i],"remaining","f");
//                 leg8->AddEntry(hKappaTPCPiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
//                 leg8->AddEntry(hKappaTPCElPi[i],"#pi^{#pm} + e^{#mp}","f");
//                 leg8->AddEntry(hKappaTPCElPK[i],"e^{#pm} + K^{#mp}(p(#bar{p}))","f");
//                 leg8->AddEntry(hKappaTPCPiPK[i],"#pi^{#pm} + K^{#mp}(p(#bar{p}))","f");
//
//               leg8->Draw();
//
//             canvas8PartKappaProj->cd(8);
//             canvas8PartKappaProj->cd(8)->SetLogy();
//
//               HistoPlotSettings(hKappaTPCAfterCut[i],"K","Counts",/*0.0000*/1,hKappaTPCAfterCut[i]->GetMaximum()*10/*1.2*/,-20,20,kBlue,kBlue,0);
//               hKappaTPCAfterCut[i]->Draw("histo");
//
//               if (status == 0) {
//                   TH1D *result                        = (TH1D*)fit->GetPlot();
//                   result->SetLineColor(kRed);
//                   result->Draw("hist,same");
//               }
//               leg2->Draw();
//
//             canvas8PartKappaProj->SaveAs(Form("%s/KappaProj8Pads_%.2d.%s",fMainPlots.Data(),i,suffix.Data()));
//
//             TCanvas* canvasRatioFitDataProj8Pads = new TCanvas("canvasRatioFitDataProj8Pads","",800,800);
//             DrawGammaCanvasSettings( canvasRatioFitDataProj8Pads,  0.12, 0.03, 0.03, 0.1);
//
//               TH1D* hRatioDataToFitFinerBG = (TH1D*)hKappaTPCAfterCut[i]->Clone("hRatioDataToFitFinerBG");
//               if (status == 0) {
//                   TH1D *result                        = (TH1D*)fit->GetPlot();
//                   result->SetLineColor(kRed);
//                   hRatioDataToFitFinerBG->Divide(hKappaTPCAfterCut[i],result, 1., 1., "B");
//                   hRatioDataToFitFinerBG->Draw("hist");
//               }
//
//             pave->Draw();
//             canvasRatioFitDataProj8Pads->SaveAs(Form("%s/RatioFitToDataFinerBG_%.2d.%s",fMainPlots.Data(),i,suffix.Data()));
//
//         }

        if(kAllKappaPlots){

            cout << "Plotting all the background contributions" << endl;
            // ==================================================================================
            TCanvas* canvasAllKappaProj         = new TCanvas("canvasAllKappaProj","",0,0,1200,800);
            DrawGammaCanvasSettings( canvasAllKappaProj,  0.15, 0.05, 0.2, 0.15);
            canvasAllKappaProj->Divide(4,3,0.001,0.001);

            TPaveText * paveAll                 = new TPaveText(0.15,0.7,0.45,0.93,"NDC");
            SetStylePave(paveAll,42,11,0.035);
            paveAll->InsertText(InfoSystem.Data());
            if(isMC) paveAll->InsertText(periodMC.Data());
            else paveAll->InsertText(periodData.Data());
            paveAll->InsertText(SigmaStarForm.Data());
            paveAll->InsertText(Form("%2.1f < #it{p}_{T} < %2.1f (GeV/#it{c})",fBinsPt[i],fBinsPt[i+1]));

            canvasAllKappaProj->cd(1);
            HistoPlotSettings(hKappaTPCElEl[i],"K","Counts",0.00001,hKappaTPCElEl[i]->GetMaximum()*1.2,-20,20,12,12,3004);
            hKappaTPCElEl[i]->Draw("histo");
            paveAll->Draw();
            paveSignal->Draw();

            canvasAllKappaProj->cd(2);
            HistoPlotSettings(hKappaTPCElPi[i],"K","Counts",0.00001,hKappaTPCElPi[i]->GetMaximum()*1.4,-20,20,kMagenta-2,kMagenta-2,3004);
            hKappaTPCElPi[i]->Draw("histo");
            paveElPi->InsertText("#pi^{#pm} + e^{#mp}");
            paveElPi->InsertText(Form("%f",hKappaTPCElPi[i]->GetEntries()/NTotalMC));
            paveElPi->Draw();

            canvasAllKappaProj->cd(3);
            HistoPlotSettings(hKappaTPCPiPi[i],"K","Counts",0.00001,hKappaTPCPiPi[i]->GetMaximum()*1.2,-20,20,kBlue-7,kBlue-7,3005);
            hKappaTPCPiPi[i]->Draw("histo");
            pavePiPi->InsertText("#pi^{#pm} + #pi^{#mp}");
            pavePiPi->InsertText(Form("%f",hKappaTPCPiPi[i]->GetEntries()/NTotalMC));
            pavePiPi->Draw();

            canvasAllKappaProj->cd(4);
            TPaveText * pavePiK                 = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(pavePiK,62,11,0.04);
            pavePiK->InsertText("#pi^{#pm} + K^{#mp}");
            pavePiK->InsertText(Form("%f",hKappaTPCPiK[i]->GetEntries()/NTotalMC));

            HistoPlotSettings(hKappaTPCPiK[i],"K","Counts",0.00001,hKappaTPCPiK[i]->GetMaximum()*1.2,-20,20,kBlue+2,kBlue+2,3006);
            hKappaTPCPiK[i]->Draw("histo");
            pavePiK->Draw();

            canvasAllKappaProj->cd(5);
            TPaveText * pavePiP                 = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(pavePiP,62,11,0.04);
            pavePiP->InsertText("#pi^{#pm} + p(#bar{p})");
            pavePiP->InsertText(Form("%f",hKappaTPCPiP[i]->GetEntries()/NTotalMC));

            HistoPlotSettings(hKappaTPCPiP[i],"K","Counts",0.00001,hKappaTPCPiP[i]->GetMaximum()*1.2,-20,20,kOrange+1,kOrange+1,3007);
            hKappaTPCPiP[i]->Draw("histo");
            pavePiP->Draw();

            canvasAllKappaProj->cd(6);
            TPaveText * paveElK                 = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveElK,62,11,0.04);
            paveElK->InsertText("e^{#pm} + K^{#mp}");
            paveElK->InsertText(Form("%f",hKappaTPCElK[i]->GetEntries()/NTotalMC));

            HistoPlotSettings(hKappaTPCElK[i],"K","Counts",0.00001,hKappaTPCElK[i]->GetMaximum()*1.4,-20,20,kGreen-5,kGreen-5,3305);
            hKappaTPCElK[i]->Draw("histo");
            paveElK->Draw();

            canvasAllKappaProj->cd(7);
            TPaveText * paveElP                 = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveElP,62,11,0.04);
            paveElP->InsertText("e^{#pm} + p(#bar{p})");
            paveElP->InsertText(Form("%f",hKappaTPCElP[i]->GetEntries()/NTotalMC));

            HistoPlotSettings(hKappaTPCElP[i],"K","Counts",0.00001,hKappaTPCElP[i]->GetMaximum()*1.4,-20,20,kRed-4,kRed-4,3395);
            hKappaTPCElP[i]->Draw("histo");
            paveElP->Draw();

            canvasAllKappaProj->cd(8);
            TPaveText * paveKK                  = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveKK,62,11,0.04);
            paveKK->InsertText("K^{#pm} + K^{#mp}");
            paveKK->InsertText(Form("%f",hKappaTPCKK[i]->GetEntries()/NTotalMC));

            HistoPlotSettings(hKappaTPCKK[i],"K","Counts",0.00001,hKappaTPCKK[i]->GetMaximum()*1.4,-20,20,kCyan-2,kCyan-2,0);
            hKappaTPCKK[i]->Draw("histo");
            paveKK->Draw();

            canvasAllKappaProj->cd(9);
            TPaveText * paveHad                 = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveHad,62,11,0.04);
            paveHad->InsertText("Hadronic");
            paveHad->InsertText(Form("%f",hKappaTPCHad[i]->GetEntries()/NTotalMC));

            HistoPlotSettings(hKappaTPCHad[i],"K","Counts",0.00001,hKappaTPCHad[i]->GetMaximum()*1.2,-20,20,kGreen+2,kGreen+2,3005);
            hKappaTPCHad[i]->Draw("histo");
            paveHad->Draw();

            canvasAllKappaProj->cd(10);
            TPaveText * paveAllRest             = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveAllRest,62,11,0.04);
            paveAllRest->InsertText("Remaining");
            paveAllRest->InsertText(Form("%f",hKappaTPCLeftoverRest[i]->GetEntries()/NTotalMC));

            HistoPlotSettings(hKappaTPCLeftoverRest[i],"K","Counts",0.00001,hKappaTPCLeftoverRest[i]->GetMaximum()*1.2,-20,20,kRed+1,kRed+1,0);
            hKappaTPCLeftoverRest[i]->Draw("histo");
            paveAllRest->Draw();

            canvasAllKappaProj->cd(11)->SetLogy();;

            HistoPlotSettings(hKappaTPCTotalSum[i],"K","Counts",1,hKappaTPCTotalSum[i]->GetMaximum()*5,-20,20,kBlack,kBlack,0);
            hKappaTPCTotalSum[i]->SetLineWidth(2.0);
            hKappaTPCTotalSum[i]->Draw("histo");
            hKappaTPCElEl[i]->Draw("histo,same");
            hKappaTPCPiPi[i]->Draw("histo,same");
            hKappaTPCElPi[i]->Draw("histo,same");
            hKappaTPCPiK[i]->Draw("histo,same");
            hKappaTPCElK[i]->Draw("histo,same");
            hKappaTPCPiP[i]->Draw("histo,same");
            hKappaTPCElP[i]->Draw("histo,same");
            hKappaTPCKK[i]->Draw("histo,same");
            hKappaTPCHad[i]->Draw("histo,same");
            hKappaTPCLeftoverRest[i]->Draw("histo,same");

            canvasAllKappaProj->cd(12);

            TLegend* legAll                     = new TLegend(0.2,0.2,0.9,0.9);
            legAll->SetTextSize(0.05);
            legAll->AddEntry(hKappaTPCTotalSum[i],"total","lp");
            legAll->AddEntry(hKappaTPCElEl[i],"real e^{+}e^{-}","f");
            legAll->AddEntry(hKappaTPCPiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
            legAll->AddEntry(hKappaTPCElPi[i],"#pi^{#pm} + e^{#mp}","f");
            legAll->AddEntry(hKappaTPCPiK[i],"#pi^{#pm} + K^{#mp}","f");
            legAll->AddEntry(hKappaTPCElK[i],"e^{#pm} + K^{#mp}","f");
            legAll->AddEntry(hKappaTPCPiP[i],"#pi^{#pm} + p(#bar{p})","f");
            legAll->AddEntry(hKappaTPCElP[i],"e^{#pm} + p(#bar{p})","f");
            legAll->AddEntry(hKappaTPCKK[i],"K^{#pm} + K^{#mp}","f");
            legAll->AddEntry(hKappaTPCHad[i],"hadronic","f");
            legAll->AddEntry(hKappaTPCLeftoverRest[i],"remaining","f");

            legAll->SetBorderSize(0);
            legAll->Draw();

            canvasAllKappaProj->SaveAs(Form("%s/KappaAllProj_%.2d.%s",fAllPlots.Data(),i,suffix.Data()));

            // ==================================================================================
            TCanvas* canvasAllProjSinglePlot = new TCanvas("canvasAllProjSinglePlot","",1000,800);
            DrawGammaCanvasSettings( canvasAllProjSinglePlot,  0.12, 0.03, 0.03, 0.1);
            canvasAllProjSinglePlot->SetLogy();

            hKappaTPCTotalSum[i]->GetYaxis()->SetRangeUser(1,1e6);
            hKappaTPCTotalSum[i]->Draw("histo");
            hKappaTPCElEl[i]->Draw("histo,same");
            hKappaTPCPiPi[i]->Draw("histo,same");
            hKappaTPCElPi[i]->Draw("histo,same");
            hKappaTPCPiK[i]->Draw("histo,same");
            hKappaTPCElK[i]->Draw("histo,same");
            hKappaTPCPiP[i]->Draw("histo,same");
            hKappaTPCElP[i]->Draw("histo,same");
            hKappaTPCKK[i]->Draw("histo,same");
            hKappaTPCHad[i]->Draw("histo,same");
            hKappaTPCLeftoverRest[i]->Draw("histo,same");

            TLegend* legAllSinglePlot = new TLegend(0.15,0.85,0.95,0.95);
            legAllSinglePlot->SetTextSize(0.035);
            legAllSinglePlot->SetFillStyle(0);
            legAllSinglePlot->SetNColumns(6);
            legAllSinglePlot->AddEntry(hKappaTPCTotalSum[i],"total","lp");
            legAllSinglePlot->AddEntry(hKappaTPCElEl[i],"real e^{+}e^{-}","f");
            legAllSinglePlot->AddEntry(hKappaTPCPiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
            legAllSinglePlot->AddEntry(hKappaTPCElPi[i],"#pi^{#pm} + e^{#mp}","f");
            legAllSinglePlot->AddEntry(hKappaTPCPiK[i],"#pi^{#pm} + K^{#mp}","f");
            legAllSinglePlot->AddEntry(hKappaTPCElK[i],"e^{#pm} + K^{#mp}","f");
            legAllSinglePlot->AddEntry(hKappaTPCPiP[i],"#pi^{#pm} + p(#bar{p})","f");
            legAllSinglePlot->AddEntry(hKappaTPCElP[i],"e^{#pm} + p(#bar{p})","f");
            legAllSinglePlot->AddEntry(hKappaTPCKK[i],"K^{#pm} + K^{#mp}","f");
            legAllSinglePlot->AddEntry(hKappaTPCHad[i],"hadronic","f");
            legAllSinglePlot->AddEntry(hKappaTPCLeftoverRest[i],"rest","f");
            legAllSinglePlot->SetBorderSize(0);
            legAllSinglePlot->Draw();

            canvasAllProjSinglePlot->SaveAs(Form("%s/KappaAllProjSinglePlot_%.2d.%s",fAllPlots.Data(),i,suffix.Data()));
        }

        TCanvas* canvasProj = new TCanvas("canvasProj","",800,800);
        DrawGammaCanvasSettings( canvasProj,  0.12, 0.03, 0.03, 0.1);
        canvasProj->SetLogy();

        TLegend* legSingleProj = new TLegend(0.7,0.7,0.9,0.93);
        legSingleProj->SetTextSize(0.04);
        legSingleProj->SetBorderSize(0);
        hKappaTPCSum[i]->Draw("histo");
        legSingleProj->AddEntry(hKappaTPCSum[i],"total","lp");
        hKappaTPCElEl[i]->Draw("histo,same");
        legSingleProj->AddEntry(hKappaTPCElEl[i],"real e^{+}e^{-}","f");
        if(moreBG==1){
          hKappaTPCPiPi[i]->Draw("histo,same");
          hKappaTPCElPi[i]->Draw("histo,same");
          hKappaTPCElPK[i]->Draw("histo,same");
          hKappaTPCPiPK[i]->Draw("histo,same");
          hKappaTPCPartialRest[i]->Draw("histo,same");
          legSingleProj->AddEntry(hKappaTPCPiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
          legSingleProj->AddEntry(hKappaTPCElPi[i],"#pi^{#pm} + e^{#mp}","f");
          legSingleProj->AddEntry(hKappaTPCElPK[i],"e^{#pm} + #pi^{#mp}, K^{#mp}, p","f");
          legSingleProj->AddEntry(hKappaTPCPiPK[i],"#pi^{#pm} + #pi^{#mp}, K^{#mp}, p","f");
          legSingleProj->AddEntry(hKappaTPCPartialRest[i],"remaining","f");
        } else if(moreBG==2){
          hKappaTPCElPK[i]->Draw("histo,same");
          hKappaTPCPiPK[i]->Draw("histo,same");
          hKappaTPCPartialRest[i]->Draw("histo,same");
          legSingleProj->AddEntry(hKappaTPCElPK[i],"e^{#pm} + #pi^{#mp}, K^{#mp}, p","f");
          legSingleProj->AddEntry(hKappaTPCPiPK[i],"#pi^{#pm} + #pi^{#mp}, K^{#mp}, p","f");
          legSingleProj->AddEntry(hKappaTPCPartialRest[i],"remaining","f");
        } else{
          hKappaTPCPiPi[i]->Draw("histo,same");
          hKappaTPCElPi[i]->Draw("histo,same");
          hKappaTPCRest[i]->Draw("histo,same");
          legSingleProj->AddEntry(hKappaTPCPiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
          legSingleProj->AddEntry(hKappaTPCElPi[i],"#pi^{#pm} + e^{#mp}","f");
          legSingleProj->AddEntry(hKappaTPCRest[i],"remaining","f");
        }
        legSingleProj->Draw();
        pave->Draw();
        canvasProj->SaveAs(Form("%s/KappaProjSinglePlot_%.2d.%s",fMainPlots.Data(),i,suffix.Data()));

        if (status == 0) {

            cout << "Bin range " << startBin << " - " << endBin << endl;
            cout << "pT range " << fBinsPt[i] << " - " << fBinsPt[i+1] << endl;

            for (int j=0; j<nTemplate; j++) fit->GetResult(j,value[j],error[j]);

            valuesElEl[i] = value[0];
            errorsElEl[i] = error[0];
            cout << " valuesElEl[" << i << "] = " << valuesElEl[i] << endl;
            cout << "   fracElEl[" << i << "] = " << fracElEl[i] << endl;
            ratiofitfractionElEl[i] = valuesElEl[i]/fracElEl[i];
            hValuesElEl->SetBinContent(i+1, valuesElEl[i]);
            hValuesElEl->SetBinError(i+1, errorsElEl[i]);
            hFractionElEl->SetBinContent(i+1, ratiofitfractionElEl[i]);
            hFractionElEl->SetBinError(i+1, 0);
            hTemplateElEl[i] = (TH1D*)hKappaTPCElEl[i]->Clone();
            hTemplateElEl[i]->Scale(valuesElEl[i]/fracElEl[i]*(NEntriesAfterCut/NTotalMC));

            if(moreBG==0){

              valuesPiPi[i] = value[1];
              errorsPiPi[i] = error[1];
              cout << " valuesPiPi[" << i << "] = " << valuesPiPi[i] << endl;
              cout << "   fracPiPi[" << i << "] = " << fracPiPi[i] << endl;
              ratiofitfractionPiPi[i] = valuesPiPi[i]/fracPiPi[i];
              hValuesPiPi->SetBinContent(i+1, valuesPiPi[i]);
              hValuesPiPi->SetBinError(i+1, errorsPiPi[i]);
              hFractionPiPi->SetBinContent(i+1, ratiofitfractionPiPi[i]);
              hFractionPiPi->SetBinError(i+1, 0);
              hTemplatePiPi[i] = (TH1D*)hKappaTPCPiPi[i]->Clone();
              hTemplatePiPi[i]->Scale(valuesPiPi[i]/fracPiPi[i]*(NEntriesAfterCut/NTotalMC));

              valuesElPi[i] = value[2];
              errorsElPi[i] = error[2];
              cout << " valuesElPi[" << i << "] = " << valuesElPi[i] << endl;
              cout << "   fracElPi[" << i << "] = " << fracElPi[i] << endl;
              ratiofitfractionElPi[i] = valuesElPi[i]/fracElPi[i];
              hValuesElPi->SetBinContent(i+1, valuesElPi[i]);
              hValuesElPi->SetBinError(i+1, errorsElPi[i]);
              hFractionElPi->SetBinContent(i+1, ratiofitfractionElPi[i]);
              hFractionElPi->SetBinError(i+1, 0);
              hTemplateElPi[i] = (TH1D*)hKappaTPCElPi[i]->Clone();
              hTemplateElPi[i]->Scale(valuesElPi[i]/fracElPi[i]*(NEntriesAfterCut/NTotalMC));

              valuesRest[i] = value[3];
              errorsRest[i] = error[3];
              cout << " valuesRest[" << i << "] = " << valuesRest[i] << endl;
              cout << "   fracRest[" << i << "] = " << fracRest[i] << endl;
              ratiofitfractionRest[i] = valuesRest[i]/fracRest[i];
              hValuesRest->SetBinContent(i+1, valuesRest[i]);
              hValuesRest->SetBinError(i+1, errorsRest[i]);
              hFractionRest->SetBinContent(i+1, ratiofitfractionRest[i]);
              hFractionRest->SetBinError(i+1, 0);
              hTemplateRest[i] = (TH1D*)hKappaTPCRest[i]->Clone();
              hTemplateRest[i]->Scale(valuesRest[i]/fracRest[i]*(NEntriesAfterCut/NTotalMC));

            } else if(moreBG==1){

              valuesPiPi[i] = value[1];
              errorsPiPi[i] = error[1];
              cout << " valuesPiPi[" << i << "] = " << valuesPiPi[i] << endl;
              cout << "   fracPiPi[" << i << "] = " << fracPiPi[i] << endl;
              ratiofitfractionPiPi[i] = valuesPiPi[i]/fracPiPi[i];
              hValuesPiPi->SetBinContent(i+1, valuesPiPi[i]);
              hValuesPiPi->SetBinError(i+1, errorsPiPi[i]);
              hFractionPiPi->SetBinContent(i+1, ratiofitfractionPiPi[i]);
              hFractionPiPi->SetBinError(i+1, 0);
              hTemplatePiPi[i] = (TH1D*)hKappaTPCPiPi[i]->Clone();
              hTemplatePiPi[i]->Scale(valuesPiPi[i]/fracPiPi[i]*(NEntriesAfterCut/NTotalMC));

              valuesElPi[i] = value[2];
              errorsElPi[i] = error[2];
              cout << " valuesElPi[" << i << "] = " << valuesElPi[i] << endl;
              cout << "   fracElPi[" << i << "] = " << fracElPi[i] << endl;
              ratiofitfractionElPi[i] = valuesElPi[i]/fracElPi[i];
              hValuesElPi->SetBinContent(i+1, valuesElPi[i]);
              hValuesElPi->SetBinError(i+1, errorsElPi[i]);
              hFractionElPi->SetBinContent(i+1, ratiofitfractionElPi[i]);
              hFractionElPi->SetBinError(i+1, 0);
              hTemplateElPi[i] = (TH1D*)hKappaTPCElPi[i]->Clone();
              hTemplateElPi[i]->Scale(valuesElPi[i]/fracElPi[i]*(NEntriesAfterCut/NTotalMC));

              valuesElPK[i] = value[3];
              errorsElPK[i] = error[3];
              cout << " valuesElPK[" << i << "] = " << valuesElPK[i] << endl;
              cout << "   fracElPK[" << i << "] = " << fracElPK[i] << endl;
              ratiofitfractionElPK[i] = valuesElPK[i]/fracElPK[i];
              hValuesElPK->SetBinContent(i+1, valuesElPK[i]);
              hValuesElPK->SetBinError(i+1, errorsElPK[i]);
              hFractionElPK->SetBinContent(i+1, ratiofitfractionElPK[i]);
              hFractionElPK->SetBinError(i+1, 0);
              hTemplateElPK[i] = (TH1D*)hKappaTPCElPK[i]->Clone();
              hTemplateElPK[i]->Scale(valuesElPK[i]/fracElPK[i]*(NEntriesAfterCut/NTotalMC));

              valuesPiPK[i] = value[4];
              errorsPiPK[i] = error[4];
              cout << " valuesPiPK[" << i << "] = " << valuesPiPK[i] << endl;
              cout << "   fracPiPK[" << i << "] = " << fracPiPK[i] << endl;
              ratiofitfractionPiPK[i] = valuesPiPK[i]/fracPiPK[i];
              hValuesPiPK->SetBinContent(i+1, valuesPiPK[i]);
              hValuesPiPK->SetBinError(i+1, errorsPiPK[i]);
              hFractionPiPK->SetBinContent(i+1, ratiofitfractionPiPK[i]);
              hFractionPiPK->SetBinError(i+1, 0);
              hTemplatePiPK[i] = (TH1D*)hKappaTPCPiPK[i]->Clone();
              hTemplatePiPK[i]->Scale(valuesPiPK[i]/fracPiPK[i]*(NEntriesAfterCut/NTotalMC));

              valuesPartialRest[i] = value[5];
              errorsPartialRest[i] = error[5];
              cout << " valuesPartialRest[" << i << "] = " << valuesPartialRest[i] << endl;
              cout << "   fracRest[" << i << "] = " << fracPartialRest[i] << endl;
              ratiofitfractionPartialRest[i] = valuesPartialRest[i]/fracPartialRest[i];
              hValuesPartialRest->SetBinContent(i+1, valuesPartialRest[i]);
              hValuesPartialRest->SetBinError(i+1, errorsPartialRest[i]);
              hFractionPartialRest->SetBinContent(i+1, ratiofitfractionPartialRest[i]);
              hFractionPartialRest->SetBinError(i+1, 0);
              hTemplatePartialRest[i] = (TH1D*)hKappaTPCPartialRest[i]->Clone();
              hTemplatePartialRest[i]->Scale(valuesPartialRest[i]/fracPartialRest[i]*(NEntriesAfterCut/NTotalMC));

            } else if(moreBG==2){

              valuesElPK[i] = value[1];
              errorsElPK[i] = error[1];
              cout << " valuesElPK[" << i << "] = " << valuesElPK[i] << endl;
              cout << "   fracElPK[" << i << "] = " << fracElPK[i] << endl;
              ratiofitfractionElPK[i] = valuesElPK[i]/fracElPK[i];
              hValuesElPK->SetBinContent(i+1, valuesElPK[i]);
              hValuesElPK->SetBinError(i+1, errorsElPK[i]);
              hFractionElPK->SetBinContent(i+1, ratiofitfractionElPK[i]);
              hFractionElPK->SetBinError(i+1, 0);
              hTemplateElPK[i] = (TH1D*)hKappaTPCElPK[i]->Clone();
              hTemplateElPK[i]->Scale(valuesElPK[i]/fracElPK[i]*(NEntriesAfterCut/NTotalMC));

              valuesPiPK[i] = value[2];
              errorsPiPK[i] = error[2];
              cout << " valuesPiPK[" << i << "] = " << valuesPiPK[i] << endl;
              cout << "   fracPiPK[" << i << "] = " << fracPiPK[i] << endl;
              ratiofitfractionPiPK[i] = valuesPiPK[i]/fracPiPK[i];
              hValuesPiPK->SetBinContent(i+1, valuesPiPK[i]);
              hValuesPiPK->SetBinError(i+1, errorsPiPK[i]);
              hFractionPiPK->SetBinContent(i+1, ratiofitfractionPiPK[i]);
              hFractionPiPK->SetBinError(i+1, 0);
              hTemplatePiPK[i] = (TH1D*)hKappaTPCPiPK[i]->Clone();
              hTemplatePiPK[i]->Scale(valuesPiPK[i]/fracPiPK[i]*(NEntriesAfterCut/NTotalMC));

              valuesPartialRest[i] = value[3];
              errorsPartialRest[i] = error[3];
              cout << " valuesPartialRest[" << i << "] = " << valuesPartialRest[i] << endl;
              cout << "   fracRest[" << i << "] = " << fracPartialRest[i] << endl;
              ratiofitfractionPartialRest[i] = valuesPartialRest[i]/fracPartialRest[i];
              hValuesPartialRest->SetBinContent(i+1, valuesPartialRest[i]);
              hValuesPartialRest->SetBinError(i+1, errorsPartialRest[i]);
              hFractionPartialRest->SetBinContent(i+1, ratiofitfractionPartialRest[i]);
              hFractionPartialRest->SetBinError(i+1, 0);
              hTemplatePartialRest[i] = (TH1D*)hKappaTPCPartialRest[i]->Clone();
              hTemplatePartialRest[i]->Scale(valuesPartialRest[i]/fracPartialRest[i]*(NEntriesAfterCut/NTotalMC));
            }

            if(fIsMC) cout << "NEntriesAfterCut from MC histo = " << NEntriesAfterCut << endl;
            else  cout << "NEntriesAfterCut from data histo = " << NEntriesAfterCut << endl;
            cout << "NTotalMC = " << NTotalMC << endl;

            hTemplateSum[i] = (TH1D*)hTemplateElEl[i]->Clone();
            hTemplateSum[i]->GetXaxis()->SetTitleOffset(0.8);
            hTemplateSum[i]->SetFillColor(0);
            hTemplateSum[i]->SetLineColor(kBlack);
            hTemplateSum[i]->SetLineWidth(2.0);
            if(moreBG==0){
              hTemplateSum[i]->Add(hTemplatePiPi[i]);
              hTemplateSum[i]->Add(hTemplateElPi[i]);
              hTemplateSum[i]->Add(hTemplateRest[i]);
            } else if(moreBG==1){
              hTemplateSum[i]->Add(hTemplatePiPi[i]);
              hTemplateSum[i]->Add(hTemplateElPi[i]);
              hTemplateSum[i]->Add(hTemplateElPK[i]);
              hTemplateSum[i]->Add(hTemplatePiPK[i]);
              hTemplateSum[i]->Add(hTemplatePartialRest[i]);
            } else if(moreBG==2){
              hTemplateSum[i]->Add(hTemplateElPK[i]);
              hTemplateSum[i]->Add(hTemplatePiPK[i]);
              hTemplateSum[i]->Add(hTemplatePartialRest[i]);
            }

            TCanvas* canvasTemplates = new TCanvas("canvasTemplates","",800,800);
            DrawGammaCanvasSettings( canvasTemplates,  0.12, 0.03, 0.04, 0.1);
            canvasTemplates->SetLogy();

            TLegend* legTemplate = new TLegend(0.7,0.66,0.9,0.93);
            legTemplate->SetTextSize(0.04);
            legTemplate->SetBorderSize(0);

            HistoPlotSettings(hTemplateSum[i],"K","Counts",1,hTemplateSum[i]->GetMaximum()*30,-20,20,kBlack,kBlack,0);
            hTemplateSum[i]->Draw("histo");
            legTemplate->AddEntry(hTemplateSum[i],"total","lp");
            hTemplateElEl[i]->Draw("histo,same");
            legTemplate->AddEntry(hTemplateElEl[i],"real e^{+}e^{-}","f");
            if(moreBG==0){
              hTemplatePiPi[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplatePiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
              hTemplateElPi[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplateElPi[i],"#pi^{#pm} + e^{#mp}","f");
              hTemplateRest[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplateRest[i],"remaining","f");
            } else if(moreBG==1){
              hTemplatePiPi[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplatePiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
              hTemplateElPi[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplateElPi[i],"#pi^{#pm} + e^{#mp}","f");
              hTemplateElPK[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplateElPK[i],"e^{#pm} + K^{#mp} + p","f");
              hTemplatePiPK[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplatePiPK[i],"#pi^{#pm} + K^{#mp} + p","f");
              hTemplatePartialRest[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplatePartialRest[i],"remaining","f");
            } else if(moreBG==2){
              hTemplateElPK[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplateElPK[i],"e^{#pm} + #pi^{#mp}, K^{#mp}, p","f");
              hTemplatePiPK[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplatePiPK[i],"#pi^{#pm} + #pi^{#mp}, K^{#mp}, p","f");
              hTemplatePartialRest[i]->Draw("histo,same");
              legTemplate->AddEntry(hTemplatePartialRest[i],"remaining","f");
            }

            if (status == 0) {
                TH1D *result                        = (TH1D*)fit->GetPlot();
                result->SetLineColor(kRed);
                result->Draw("hist,same");
                legTemplate->AddEntry(result,"fit","lp");
            }
            hKappaTPCAfterCut[i]->Draw("histo,same");
            legTemplate->AddEntry(hKappaTPCAfterCut[i],"data","lp");
            legTemplate->Draw();
            pave->Draw();
            canvasTemplates->SaveAs(Form("%s/Templates_%.2d.%s",fMainPlots.Data(),i,suffix.Data()));

            // calculate signal purity
            signalPurityNum1    = hTemplateElEl[i]->Integral(hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
            if(moreBG==0){
              backgroundPiPiNum1  = hTemplatePiPi[i]->Integral(hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
              backgroundElPiNum1  = hTemplateElPi[i]->Integral(hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
              backgroundRestNum1  = hTemplateRest[i]->Integral(hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
            } else if(moreBG==1){
              backgroundPiPiNum1  = hTemplatePiPi[i]->Integral(hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
              backgroundElPiNum1  = hTemplateElPi[i]->Integral(hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
              backgroundElPKNum1  = hTemplateElPK[i]->Integral(hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
              backgroundPiPKNum1  = hTemplatePiPK[i]->Integral(hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
              backgroundRestNum1  = hTemplatePartialRest[i]->Integral(hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
            } else if(moreBG==2){
              backgroundElPKNum1  = hTemplateElPK[i]->Integral(hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
              backgroundPiPKNum1  = hTemplatePiPK[i]->Integral(hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
              backgroundRestNum1  = hTemplatePartialRest[i]->Integral(hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
            }

            signalPurityNum2    = hTemplateElEl[i]->Integral(hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
            if(moreBG==0){
              backgroundPiPiNum2  = hTemplatePiPi[i]->Integral(hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
              backgroundElPiNum2  = hTemplateElPi[i]->Integral(hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
              backgroundRestNum2  = hTemplateRest[i]->Integral(hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
            } else if(moreBG==1){
              backgroundPiPiNum2  = hTemplatePiPi[i]->Integral(hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
              backgroundElPiNum2  = hTemplateElPi[i]->Integral(hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
              backgroundElPKNum2  = hTemplateElPK[i]->Integral(hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
              backgroundPiPKNum2  = hTemplatePiPK[i]->Integral(hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
              backgroundRestNum2  = hTemplatePartialRest[i]->Integral(hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
            } else if(moreBG==2){
              backgroundElPKNum2  = hTemplateElPK[i]->Integral(hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
              backgroundPiPKNum2  = hTemplatePiPK[i]->Integral(hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
              backgroundRestNum2  = hTemplatePartialRest[i]->Integral(hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
            }

            signalPurityNum3    = hTemplateElEl[i]->Integral(hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
            if(moreBG==0){
              backgroundPiPiNum3  = hTemplatePiPi[i]->Integral(hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
              backgroundElPiNum3  = hTemplateElPi[i]->Integral(hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
              backgroundRestNum3  = hTemplateRest[i]->Integral(hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
            } else if(moreBG==1){
              backgroundPiPiNum3  = hTemplatePiPi[i]->Integral(hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
              backgroundElPiNum3  = hTemplateElPi[i]->Integral(hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
              backgroundElPKNum3  = hTemplateElPK[i]->Integral(hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
              backgroundPiPKNum3  = hTemplatePiPK[i]->Integral(hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
              backgroundRestNum3  = hTemplatePartialRest[i]->Integral(hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
            } else if(moreBG==2){
              backgroundElPKNum3  = hTemplateElPK[i]->Integral(hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateElPK[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
              backgroundPiPKNum3  = hTemplatePiPK[i]->Integral(hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplatePiPK[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
              backgroundRestNum3  = hTemplatePartialRest[i]->Integral(hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplatePartialRest[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
            }

            signalPurityDenom1 = hTemplateSum[i]->Integral(hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
            signalPurityDenom2 = hTemplateSum[i]->Integral(hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
            signalPurityDenom3 = hTemplateSum[i]->Integral(hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));

            signalPurity1                    = signalPurityNum1 / signalPurityDenom1;
            if(moreBG==0){
              backgroundPiPi1                  = backgroundPiPiNum1 / signalPurityDenom1;
              backgroundElPi1                  = backgroundElPiNum1 / signalPurityDenom1;
              backgroundRest1                  = backgroundRestNum1 / signalPurityDenom1;
            } else if(moreBG==1){
              backgroundPiPi1                  = backgroundPiPiNum1 / signalPurityDenom1;
              backgroundElPi1                  = backgroundElPiNum1 / signalPurityDenom1;
              backgroundElPK1                  = backgroundElPKNum1 / signalPurityDenom1;
              backgroundPiPK1                  = backgroundPiPKNum1 / signalPurityDenom1;
              backgroundRest1                  = backgroundRestNum1 / signalPurityDenom1;
            } else if(moreBG==2){
              backgroundElPK1                  = backgroundElPKNum1 / signalPurityDenom1;
              backgroundPiPK1                  = backgroundPiPKNum1 / signalPurityDenom1;
              backgroundRest1                  = backgroundRestNum1 / signalPurityDenom1;
            }
//             cout << "------------------------------------------------------------" << endl;
//             cout << "Purity for bin " << i << " in narrow range is " << signalPurity1 << endl;
//             cout << "BG PiPi for bin " << i << " in narrow range is " << backgroundPiPi1 << endl;
//             cout << "BG ElPi for bin " << i << " in narrow range is " << backgroundElPi1 << endl;
//             cout << "BG rest for bin " << i << " in narrow range is " << backgroundRest1 << endl;
//             cout << "------------------------------------------------------------" << endl;

            signalPurity2                    = signalPurityNum2 / signalPurityDenom2;
            signalPurity2                    = signalPurityNum2 / signalPurityDenom2;
            if(moreBG==0){
              backgroundPiPi2                  = backgroundPiPiNum2 / signalPurityDenom2;
              backgroundElPi2                  = backgroundElPiNum2 / signalPurityDenom2;
              backgroundRest2                  = backgroundRestNum2 / signalPurityDenom2;
            } else if(moreBG==1){
              backgroundPiPi2                  = backgroundPiPiNum2 / signalPurityDenom2;
              backgroundElPi2                  = backgroundElPiNum2 / signalPurityDenom2;
              backgroundElPK2                  = backgroundElPKNum2 / signalPurityDenom2;
              backgroundPiPK2                  = backgroundPiPKNum2 / signalPurityDenom2;
              backgroundRest2                  = backgroundRestNum2 / signalPurityDenom2;
            } else if(moreBG==2){
              backgroundElPK2                  = backgroundElPKNum2 / signalPurityDenom2;
              backgroundPiPK2                  = backgroundPiPKNum2 / signalPurityDenom2;
              backgroundRest2                  = backgroundRestNum2 / signalPurityDenom2;
            }
//             cout << "------------------------------------------------------------" << endl;
//             cout << "Purity for bin " << i << " in wide range is " << signalPurity2 << endl;
//             cout << "BG PiPi for bin " << i << " in wide range is " << backgroundPiPi2 << endl;
//             cout << "BG ElPi for bin " << i << " in wide range is " << backgroundElPi2 << endl;
//             cout << "BG rest for bin " << i << " in wide range is " << backgroundRest2 << endl;
//             cout << "------------------------------------------------------------" << endl;

            signalPurity3                    = signalPurityNum3 / signalPurityDenom3;
            if(moreBG==0){
              backgroundPiPi3                  = backgroundPiPiNum3 / signalPurityDenom3;
              backgroundElPi3                  = backgroundElPiNum3 / signalPurityDenom3;
              backgroundRest3                  = backgroundRestNum3 / signalPurityDenom3;
            } else if(moreBG==1){
              backgroundPiPi3                  = backgroundPiPiNum3 / signalPurityDenom3;
              backgroundElPi3                  = backgroundElPiNum3 / signalPurityDenom3;
              backgroundElPK3                  = backgroundElPKNum3 / signalPurityDenom3;
              backgroundPiPK3                  = backgroundPiPKNum3 / signalPurityDenom3;
              backgroundRest3                  = backgroundRestNum3 / signalPurityDenom3;
            } else if(moreBG==2){
              backgroundElPK3                  = backgroundElPKNum3 / signalPurityDenom3;
              backgroundPiPK3                  = backgroundPiPKNum3 / signalPurityDenom3;
              backgroundRest3                  = backgroundRestNum3 / signalPurityDenom3;
            }
//             cout << "------------------------------------------------------------" << endl;
//             cout << "Purity for bin " << i << " in asymm range is " << signalPurity3 << endl;
//             cout << "BG PiPi for bin " << i << " in asymm range is " << backgroundPiPi3 << endl;
//             cout << "BG ElPi for bin " << i << " in asymm range is " << backgroundElPi3 << endl;
//             cout << "BG rest for bin " << i << " in asymm range is " << backgroundRest3 << endl;
//             cout << "------------------------------------------------------------" << endl;

            //signalPurityErr                 = TMath::Sqrt(signalPurityNumErr*signalPurityNumErr + signalPurityDenomErr*signalPurityDenomErr);

            hSignalPurity1->SetBinContent(i+1, signalPurity1);
            hSignalPurity1->SetBinError(i+1, 0);
            hSignalPurity2->SetBinContent(i+1, signalPurity2);
            hSignalPurity2->SetBinError(i+1, 0);
            hSignalPurity3->SetBinContent(i+1, signalPurity3);
            hSignalPurity3->SetBinError(i+1, 0);

            if(moreBG==0){
              hBackgroundPiPi1->SetBinContent(i+1, backgroundPiPi1);
              hBackgroundPiPi1->SetBinError(i+1, 0);
              hBackgroundPiPi2->SetBinContent(i+1, backgroundPiPi2);
              hBackgroundPiPi2->SetBinError(i+1, 0);
              hBackgroundPiPi3->SetBinContent(i+1, backgroundPiPi3);
              hBackgroundPiPi3->SetBinError(i+1, 0);

              hBackgroundElPi1->SetBinContent(i+1, backgroundElPi1);
              hBackgroundElPi1->SetBinError(i+1, 0);
              hBackgroundElPi2->SetBinContent(i+1, backgroundElPi2);
              hBackgroundElPi2->SetBinError(i+1, 0);
              hBackgroundElPi3->SetBinContent(i+1, backgroundElPi3);
              hBackgroundElPi3->SetBinError(i+1, 0);

            } else if(moreBG==1){
              hBackgroundPiPi1->SetBinContent(i+1, backgroundPiPi1);
              hBackgroundPiPi1->SetBinError(i+1, 0);
              hBackgroundPiPi2->SetBinContent(i+1, backgroundPiPi2);
              hBackgroundPiPi2->SetBinError(i+1, 0);
              hBackgroundPiPi3->SetBinContent(i+1, backgroundPiPi3);
              hBackgroundPiPi3->SetBinError(i+1, 0);

              hBackgroundElPi1->SetBinContent(i+1, backgroundElPi1);
              hBackgroundElPi1->SetBinError(i+1, 0);
              hBackgroundElPi2->SetBinContent(i+1, backgroundElPi2);
              hBackgroundElPi2->SetBinError(i+1, 0);
              hBackgroundElPi3->SetBinContent(i+1, backgroundElPi3);
              hBackgroundElPi3->SetBinError(i+1, 0);

              hBackgroundElPK1->SetBinContent(i+1, backgroundElPK1);
              hBackgroundElPK1->SetBinError(i+1, 0);
              hBackgroundElPK2->SetBinContent(i+1, backgroundElPK2);
              hBackgroundElPK2->SetBinError(i+1, 0);
              hBackgroundElPK3->SetBinContent(i+1, backgroundElPK3);
              hBackgroundElPK3->SetBinError(i+1, 0);

              hBackgroundPiPK1->SetBinContent(i+1, backgroundPiPK1);
              hBackgroundPiPK1->SetBinError(i+1, 0);
              hBackgroundPiPK2->SetBinContent(i+1, backgroundPiPK2);
              hBackgroundPiPK2->SetBinError(i+1, 0);
              hBackgroundPiPK3->SetBinContent(i+1, backgroundPiPK3);
              hBackgroundPiPK3->SetBinError(i+1, 0);

              hBackgroundRest1->SetBinContent(i+1, backgroundRest1);
              hBackgroundRest1->SetBinError(i+1, 0);
              hBackgroundRest2->SetBinContent(i+1, backgroundRest2);
              hBackgroundRest2->SetBinError(i+1, 0);
              hBackgroundRest3->SetBinContent(i+1, backgroundRest3);
              hBackgroundRest3->SetBinError(i+1, 0);

            } else if(moreBG==2){
              hBackgroundElPK1->SetBinContent(i+1, backgroundElPK1);
              hBackgroundElPK1->SetBinError(i+1, 0);
              hBackgroundElPK2->SetBinContent(i+1, backgroundElPK2);
              hBackgroundElPK2->SetBinError(i+1, 0);
              hBackgroundElPK3->SetBinContent(i+1, backgroundElPK3);
              hBackgroundElPK3->SetBinError(i+1, 0);

              hBackgroundPiPK1->SetBinContent(i+1, backgroundPiPK1);
              hBackgroundPiPK1->SetBinError(i+1, 0);
              hBackgroundPiPK2->SetBinContent(i+1, backgroundPiPK2);
              hBackgroundPiPK2->SetBinError(i+1, 0);
              hBackgroundPiPK3->SetBinContent(i+1, backgroundPiPK3);
              hBackgroundPiPK3->SetBinError(i+1, 0);

              hBackgroundRest1->SetBinContent(i+1, backgroundRest1);
              hBackgroundRest1->SetBinError(i+1, 0);
              hBackgroundRest2->SetBinContent(i+1, backgroundRest2);
              hBackgroundRest2->SetBinError(i+1, 0);
              hBackgroundRest3->SetBinContent(i+1, backgroundRest3);
              hBackgroundRest3->SetBinError(i+1, 0);

            }

        } //end status == 0 if
        else {
            errorCounter++;
            errorBin[i] = 1;
            errorStatus[i] = status;
        }

        delete kappaArray;
        delete fit;
    } //end of pt loop

    TCanvas* canvasPurity   = new TCanvas("canvasPurity","",0,0,1000,800);  // gives the page size
    DrawGammaCanvasSettings( canvasPurity,  0.1, 0.02, 0.02, 0.09);
    TH2D *histo2DPurity     = new TH2D("histo2DPurity", "histo2DPurity", 1000,0.,fDeltaPt->GetXaxis()->GetXmax(),1000,0.7,1.01);
    SetStyleHistoTH2ForGraphs(histo2DPurity, "#it{p}_{T} (GeV/#it{c})","Purity", 0.035,0.04, 0.035,0.04, 1.,1.05);
    histo2DPurity->DrawCopy();

    DrawGammaSetMarker(hSignalPurity1,20,1.5,kBlack,kBlack);
    hSignalPurity1->Draw("same,p");
    DrawGammaSetMarker(hSignalPurity2,20,1.5,kGreen+2,kGreen+2);
    hSignalPurity2->Draw("same,p");
    DrawGammaSetMarker(hSignalPurity3,20,1.5,kBlue+1,kBlue+1);
    hSignalPurity3->Draw("same,p");

    TLegend* legPurity = new TLegend(0.2,0.2,0.5,0.4);
    legPurity->SetTextSize(0.04);
    legPurity->SetHeader("Purity");
    legPurity->AddEntry(hSignalPurity1,Form("%.2f < K < %.2f",kappaRangeElEl1[0],kappaRangeElEl1[1]),"p");
    legPurity->AddEntry(hSignalPurity2,Form("%.2f < K < %.2f",kappaRangeElEl2[0],kappaRangeElEl2[1]),"p");
    legPurity->AddEntry(hSignalPurity3,Form("%.2f < K < %.2f",kappaRangeElEl3[0],kappaRangeElEl3[1]),"p");
    legPurity->SetBorderSize(0);
    legPurity->Draw("same");

    canvasPurity->SaveAs(Form("%s/Purity.%s",fOutputDir.Data(),suffix.Data()));

    TCanvas* canvasBackground       = new TCanvas("canvasBackground","",0,0,1500,600);  // gives the page size
    DrawGammaCanvasSettings( canvasBackground,  0.07, 0.0, 0.02, 0.09);
    canvasBackground->Divide(3,1,0.000001,0.00001);
    TH2D *histo2DBackground         = new TH2D("histo2DBackground", "histo2DBackground", 1000,0.,fDeltaPt->GetXaxis()->GetXmax(),10,-.001,.1);
    SetStyleHistoTH2ForGraphs(histo2DBackground, "#it{p}_{T} (GeV/#it{c})","Background", 0.035,0.04, 0.035,0.04, 1.,1.4);

    DrawGammaSetMarker(hBackgroundElPi1,25,1.,kMagenta-2,kMagenta-2);
    DrawGammaSetMarker(hBackgroundElPi2,25,1.,kMagenta-2,kMagenta-2);
    DrawGammaSetMarker(hBackgroundElPi3,25,1.,kMagenta-2,kMagenta-2);
    DrawGammaSetMarker(hBackgroundPiPi1,24,1.,kBlue-7,kBlue-7);
    DrawGammaSetMarker(hBackgroundPiPi2,24,1.,kBlue-7,kBlue-7);
    DrawGammaSetMarker(hBackgroundPiPi3,24,1.,kBlue-7,kBlue-7);
    DrawGammaSetMarker(hBackgroundElPK1,21,1.,kMagenta-2,kMagenta-2);
    DrawGammaSetMarker(hBackgroundElPK2,21,1.,kMagenta-2,kMagenta-2);
    DrawGammaSetMarker(hBackgroundElPK3,21,1.,kMagenta-2,kMagenta-2);
    DrawGammaSetMarker(hBackgroundPiPK1,20,1.,kBlue-7,kBlue-7);
    DrawGammaSetMarker(hBackgroundPiPK2,20,1.,kBlue-7,kBlue-7);
    DrawGammaSetMarker(hBackgroundPiPK3,20,1.,kBlue-7,kBlue-7);
    DrawGammaSetMarker(hBackgroundRest1,27,1.,kRed,kRed);
    DrawGammaSetMarker(hBackgroundRest2,27,1.,kRed,kRed);
    DrawGammaSetMarker(hBackgroundRest3,27,1.,kRed,kRed);

    canvasBackground->cd(1);
    histo2DBackground->DrawCopy();

    TLegend* legBackground1 = new TLegend(0.2,0.7,0.5,0.9);
    legBackground1->SetTextSize(0.04);
    legBackground1->SetHeader(Form("Contamination in %.2f < K < %.2f",kappaRangeElEl1[0],kappaRangeElEl1[1]));

    if(moreBG==0){
      hBackgroundPiPi1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundPiPi1,"#pi^{#pm}#pi^{#mp}","p");
      hBackgroundElPi1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundElPi1,"e^{#pm}#pi^{#mp}","p");
      hBackgroundRest1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundRest1,"rest","p");
    } else if(moreBG==1){
      hBackgroundPiPi1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundPiPi1,"#pi^{#pm}#pi^{#mp}","p");
      hBackgroundElPi1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundElPi1,"e^{#pm}#pi^{#mp}","p");
      hBackgroundElPK1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundElPK1,"e^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundPiPK1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundPiPK1,"#pi^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundRest1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundRest1,"rest","p");
    } else if(moreBG==2){
      hBackgroundElPK1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundElPK1,"e^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundPiPK1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundPiPK1,"#pi^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundRest1->Draw("same,p");
      legBackground1->AddEntry(hBackgroundRest1,"rest","p");
    }

    legBackground1->SetBorderSize(0);
    legBackground1->Draw();

    canvasBackground->cd(2);
    histo2DBackground->DrawCopy();

    TLegend* legBackground2 = new TLegend(0.2,0.7,0.5,0.9);
    legBackground2->SetTextSize(0.04);
    legBackground2->SetHeader(Form("Contamination in %.2f < K < %.2f",kappaRangeElEl2[0],kappaRangeElEl2[2]));

    if(moreBG==0){
      hBackgroundPiPi2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundPiPi2,"#pi^{#pm}#pi^{#mp}","p");
      hBackgroundElPi2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundElPi2,"e^{#pm}#pi^{#mp}","p");
      hBackgroundRest2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundRest2,"rest","p");
    } else if(moreBG==2){
      hBackgroundPiPi2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundPiPi2,"#pi^{#pm}#pi^{#mp}","p");
      hBackgroundElPi2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundElPi2,"e^{#pm}#pi^{#mp}","p");
      hBackgroundElPK2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundElPK2,"e^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundPiPK2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundPiPK2,"#pi^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundRest2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundRest2,"rest","p");
    } else if(moreBG==2){
      hBackgroundElPK2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundElPK2,"e^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundPiPK2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundPiPK2,"#pi^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundRest2->Draw("same,p");
      legBackground2->AddEntry(hBackgroundRest2,"rest","p");
    }

    legBackground2->SetBorderSize(0);
    legBackground2->Draw();

    canvasBackground->cd(3);
    histo2DBackground->DrawCopy();

    TLegend* legBackground3 = new TLegend(0.2,0.7,0.5,0.9);
    legBackground3->SetTextSize(0.04);
    legBackground3->SetHeader(Form("Contamination in %.2f < K < %.2f",kappaRangeElEl3[0],kappaRangeElEl3[3]));

    if(moreBG==0){
      hBackgroundPiPi3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundPiPi3,"#pi^{#pm}#pi^{#mp}","p");
      hBackgroundElPi3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundElPi3,"e^{#pm}#pi^{#mp}","p");
      hBackgroundRest3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundRest3,"rest","p");
    } else if(moreBG==3){
      hBackgroundPiPi3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundPiPi3,"#pi^{#pm}#pi^{#mp}","p");
      hBackgroundElPi3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundElPi3,"e^{#pm}#pi^{#mp}","p");
      hBackgroundElPK3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundElPK3,"e^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundPiPK3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundPiPK3,"#pi^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundRest3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundRest3,"rest","p");
    } else if(moreBG==2){
      hBackgroundElPK3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundElPK3,"e^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundPiPK3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundPiPK3,"#pi^{#pm} + #pi^{#mp},K^{#mp},p","p");
      hBackgroundRest3->Draw("same,p");
      legBackground3->AddEntry(hBackgroundRest3,"rest","p");
    }

    legBackground3->SetBorderSize(0);
    legBackground3->Draw();

    canvasBackground->SaveAs(Form("%s/Background.%s",fOutputDir.Data(),suffix.Data()));

    HistoMarkerPlotSettings(hValuesElEl,"#it{p}_{T} (GeV/#it{c})","values",0.,1.1,0.,14.,12,12,20);
    HistoMarkerPlotSettings(hValuesElPi,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kMagenta-2,kMagenta-2,20);
    HistoMarkerPlotSettings(hValuesPiPi,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kBlue-7,kBlue-7,20);
    HistoMarkerPlotSettings(hValuesElPK,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kMagenta-2,kMagenta-2,21);
    HistoMarkerPlotSettings(hValuesPiPK,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kBlue-7,kBlue-7,21);
    HistoMarkerPlotSettings(hValuesRest,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kRed+1,kRed+1,20);

    HistoMarkerPlotSettings(hFractionElEl,"#it{p}_{T} (GeV/#it{c})","fractions",0.5,1.5,0.,14.,12,12,20);
    HistoMarkerPlotSettings(hFractionElPi,"#it{p}_{T} (GeV/#it{c})","",0.5,1.5,0.,14.,kMagenta-2,kMagenta-2,20);
    HistoMarkerPlotSettings(hFractionPiPi,"#it{p}_{T} (GeV/#it{c})","",0.5,1.5,0.,14.,kBlue-7,kBlue-7,20);
    HistoMarkerPlotSettings(hFractionElPK,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kMagenta-2,kMagenta-2,21);
    HistoMarkerPlotSettings(hFractionPiPK,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kBlue-7,kBlue-7,21);
    HistoMarkerPlotSettings(hFractionRest,"#it{p}_{T} (GeV/#it{c})","",0.5,1.5,0.,14.,kRed+1,kRed+1,20);

    TCanvas* canvasFractions     = new TCanvas("canvasFractions","",0,0,1000,500);  // gives the page size
    DrawGammaCanvasSettings( canvasFractions,  0.07, 0.02, 0.02, 0.09);
    canvasFractions->Divide(2,1,0.0001,0.0001);

    canvasFractions->cd(1);

    hValuesElEl->Draw("p");
    TLegend* legFracTemplate = new TLegend(0.5,0.35,0.9,0.55);
    legFracTemplate->SetTextSize(0.04);
    legFracTemplate->AddEntry(hValuesElEl,"e^{#pm}e^{#mp}","p");

    if(hValuesPiPi){
      hValuesPiPi->Draw("p,same");
      legFracTemplate->AddEntry(hValuesPiPi,"#pi^{#pm}#pi^{#mp}","p");
    }
    if(hValuesElPi){
      hValuesElPi->Draw("p,same");
      legFracTemplate->AddEntry(hValuesElPi,"e^{#pm}#pi^{#mp}","p");

    }
    if(hValuesElPK){
      hValuesElPK->Draw("p,same");
      legFracTemplate->AddEntry(hValuesElPK,"e^{#pm} + #pi^{#mp}, K^{#mp}, p","p");
    }
    if(hValuesPiPK){
      hValuesPiPK->Draw("p,same");
      legFracTemplate->AddEntry(hValuesPiPK,"#pi^{#pm} + #pi^{#mp}, K^{#mp}, p","p");
    }
    if(hValuesRest){
      hValuesRest->Draw("p,same");
      legFracTemplate->AddEntry(hValuesRest,"rest","p");
    }

    legFracTemplate->SetBorderSize(0);
    legFracTemplate->Draw();
    TPaveText * pave1 = new TPaveText(0.5,0.56,0.9,0.7,"NDC");
    SetStylePave(pave1,42,11,0.035);
    pave1->InsertText(InfoSystem.Data());
    if(isMC) pave1->InsertText(periodMC.Data());
    else pave1->InsertText(periodData.Data());
    pave1->Draw();

    canvasFractions->cd(2);

    hFractionElEl->Draw("p");
    if(hFractionElPi->GetEntries()!=0)hFractionElPi->Draw("p,same");
    if(hFractionPiPi->GetEntries()!=0)hFractionPiPi->Draw("p,same");
    if(hFractionElPK->GetEntries()!=0)hFractionElPK->Draw("p,same");
    if(hFractionPiPK->GetEntries()!=0)hFractionPiPK->Draw("p,same");
    if(hFractionRest->GetEntries()!=0)hFractionRest->Draw("p,same");

    canvasFractions->SaveAs(Form("%s/Fractions.%s",fOutputDir.Data(),suffix.Data()));

    // status != 0
    cout << "=====================================================================================" << endl;
    cout << "fit failed " << errorCounter << " times!" << endl;
    for (Int_t i=fStartPtBin; i<fNBinsPt; i++) {
        if (errorBin[i]) {
            cout << "bin " << i << ":\t" << fBinsPt[i] << "\tto\t" << fBinsPt[i+1] << "\twith status " << errorStatus[i] << endl;
        }
    }
    cout << "=====================================================================================" << endl;

    // save histograms
//     SaveHistos(cutSelection, moreBG);
}


// *****************************************************************************************************
// *********************** Initialize histograms and binning *******************************************
// *****************************************************************************************************
void Initialize(TString setPi0, TString energy , Int_t numberOfBins, Int_t mode, Bool_t addSig){

    InitializeBinning(setPi0, numberOfBins, energy, fDirectPhoton, fMode, fEventCutNumber, fClusterCutNumber);

    fDeltaPt                                            = new TH1D("deltaPt","",fNBinsPt,fBinsPt);
    for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){
        fDeltaPt->SetBinContent(    iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
        fDeltaPt->SetBinError(      iPt,0);
    }
}

//**************************************************************************************************
//******************************** Function to rebin spectra histos ********************************
//**************************************************************************************************
void RebinSpectrum(TH1D *Spectrum, TString NewName){
    if(NewName.CompareTo("")) NewName = Spectrum->GetName();

    if (fUseAnalysisBinning) {
        *Spectrum = *((TH1D*)Spectrum->Rebin(fNBinsPt,NewName,fBinsPt));
        Spectrum->Divide(fDeltaPt);
    } else {
        *Spectrum = *((TH1D*)Spectrum->Rebin(fNBinsPt,NewName,fBinsPt));
        Spectrum->Divide(fDeltaPt);
    }
}

// *****************************************************************************************************
// *********************** Save histograms to file *****************************************************
// *****************************************************************************************************
void SaveHistos(TString cutString, Bool_t moreBG){
    const char* nameOutput  = Form("%s/%s/PurityFromTemplate_%s.root",cutString.Data(),fEnergyFlag.Data(), cutString.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *file             = new TFile(nameOutput,"UPDATE");

    // write kappa projections in pt bins
    for (Int_t i=0; i<fNBinsPt; i++) {
        if (hKappaTPCAfterCut[i])   hKappaTPCAfterCut[i]->Write(hKappaTPCAfterCut[i]->GetName(),TObject::kOverwrite);
        if (hKappaTPCElEl[i])       hKappaTPCElEl[i]->Write(hKappaTPCElEl[i]->GetName(),TObject::kOverwrite);
        if (hKappaTPCElPi[i])       hKappaTPCElPi[i]->Write(hKappaTPCElPi[i]->GetName(),TObject::kOverwrite);
        if (hKappaTPCPiPi[i])       hKappaTPCPiPi[i]->Write(hKappaTPCPiPi[i]->GetName(),TObject::kOverwrite);
        if(moreBG){
          if (hKappaTPCElPK[i])       hKappaTPCElPK[i]->Write(hKappaTPCElPK[i]->GetName(),TObject::kOverwrite);
          if (hKappaTPCPiPK[i])       hKappaTPCPiPK[i]->Write(hKappaTPCPiPK[i]->GetName(),TObject::kOverwrite);
          if (hKappaTPCPartialRest[i])hKappaTPCPartialRest[i]->Write(hKappaTPCPartialRest[i]->GetName(),TObject::kOverwrite);
        } else
          if (hKappaTPCRest[i])       hKappaTPCRest[i]->Write(hKappaTPCRest[i]->GetName(),TObject::kOverwrite);
    }

    // write templates
    for (Int_t i=0; i<fNBinsPt; i++) {
        if (hTemplateElEl[i])   hTemplateElEl[i]->Write(hTemplateElEl[i]->GetName(),TObject::kOverwrite);
        if (hTemplateElPi[i])   hTemplateElPi[i]->Write(hTemplateElPi[i]->GetName(),TObject::kOverwrite);
        if (hTemplatePiPi[i])   hTemplatePiPi[i]->Write(hTemplatePiPi[i]->GetName(),TObject::kOverwrite);
        if(moreBG){
          if (hTemplateElPK[i])       hTemplateElPK[i]->Write(hTemplateElPK[i]->GetName(),TObject::kOverwrite);
          if (hTemplatePiPK[i])       hTemplatePiPK[i]->Write(hTemplatePiPK[i]->GetName(),TObject::kOverwrite);
          if (hTemplatePartialRest[i])hTemplatePartialRest[i]->Write(hTemplatePartialRest[i]->GetName(),TObject::kOverwrite);
        } else
          if (hTemplateRest[i])   hTemplateRest[i]->Write(hTemplateRest[i]->GetName(),TObject::kOverwrite);
    }

    // write purity
    if (hSignalPurity1)          hSignalPurity1->Write(hSignalPurity1->GetName(),TObject::kOverwrite);
    if (hSignalPurity2)          hSignalPurity2->Write(hSignalPurity2->GetName(),TObject::kOverwrite);
    if (hSignalPurity3)          hSignalPurity3->Write(hSignalPurity3->GetName(),TObject::kOverwrite);

}

void HistoPlotSettings(TH1* histo1,
                       TString XTitle,
                       TString YTitle,
                       Double_t YMin,
                       Double_t YMax,
                       Double_t XMin,
                       Double_t XMax,
                       Color_t lineColor,
                       Color_t fillColor,
                       Style_t fillStyle)
{

    histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    histo1->GetXaxis()->SetRangeUser(XMin, XMax);

    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }

    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
    histo1->GetYaxis()->SetLabelSize(0.035);
    histo1->GetYaxis()->SetTitleSize(0.043);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetXaxis()->SetTitleSize(0.043);
    histo1->GetXaxis()->SetLabelSize(0.035);
    histo1->GetXaxis()->SetTitleOffset(0.8);
    histo1->GetYaxis()->SetTitleOffset(1.2);

    histo1->SetFillColor(fillColor);
    histo1->SetLineColor(lineColor);
    histo1->SetFillStyle(fillStyle);
}

void HistoMarkerPlotSettings(TH1* histo1,
                             TString XTitle,
                             TString YTitle,
                             Double_t YMin,
                             Double_t YMax,
                             Double_t XMin,
                             Double_t XMax,
                             Color_t lineColor,
                             Color_t markerColor,
                             Style_t markerStyle)
{

    histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    histo1->GetXaxis()->SetRangeUser(XMin, XMax);

    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }

    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
    histo1->GetYaxis()->SetLabelSize(0.035);
    histo1->GetYaxis()->SetTitleSize(0.043);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetXaxis()->SetTitleSize(0.043);
    histo1->GetXaxis()->SetLabelSize(0.035);
    histo1->GetXaxis()->SetTitleOffset(0.8);
    histo1->GetYaxis()->SetTitleOffset(1.2);

    histo1->SetMarkerColor(markerColor);
    histo1->SetLineColor(lineColor);
    histo1->SetMarkerStyle(markerStyle);
}

void SetStylePave(TPaveText *pave,
                  Style_t textFontStyle,
                  Style_t textAlign,
                  Size_t textSize)
{
    pave->SetTextFont(textFontStyle);
    pave->SetTextSize(textSize);
    pave->SetBorderSize(0);
    pave->SetFillColor(kWhite);
    pave->SetFillStyle(0);
    pave->SetTextAlign(textAlign);
}