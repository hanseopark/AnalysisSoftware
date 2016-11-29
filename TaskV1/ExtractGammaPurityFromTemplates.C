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

void ExtractGammaPurityFromTemplates(Int_t Trainconfig    = 51,
                                     TString fileData     = "/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20161116-1/GA_PbPb_MC-326_20161117-1504/merge_runlist_1/GammaConvFlow_51.root",
                                     TString fileMC       = "/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20161116-1/GA_PbPb_MC-326_20161117-1504/merge_runlist_1/GammaConvFlow_51.root",
                                     TString cutSelection = "50200013_00200009007000008250400000",
                                     Bool_t isMC          = 1,
                                     TString suffix       = "pdf",
                                     TString option       = "PbPb_2.76TeV",
                                     TString periodData   = "LHC11h",
                                     TString periodMC     = "LHC14a1a",
                                     Int_t   numberOfBins = 21,
                                     Int_t   fMode        = 12)
{

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
    fEnergyFlag                                 = option;
    fCollisionSystem    = ReturnFullCollisionsSystem(fEnergyFlag);

    fSuffix                                     = suffix;
    cout << "Pictures are saved as " << suffix.Data() << endl;
    //************************************ Load binning for spectrum ***********************************
    //meson and mode are hardcoded
    Initialize("Pi0", fEnergyFlag, numberOfBins, 0, kFALSE);
    //************************************ Detect correct folder name **********************************
    TString nameMainDir = Form("GammaConvV1_%i_v2",Trainconfig);

    //************************************ Separate cutstrings *****************************************
    fCutSelection                               = cutSelection;
    fCutSelectionRead                           = cutSelection;
    ReturnSeparatedCutNumberAdvanced( fCutSelection,fEventCutNumber, fGammaCutNumber, fClusterCutNumber, fElectronCutNumber, fMesonCutNumber, fMode);

    fEventCutSelection                        = fEventCutNumber.Data();
    fGammaCutSelection                        = fGammaCutNumber.Data();

    TString centralityString        = GetCentralityString(fEventCutSelection);
    TString InfoSystem              = Form("%s %s",centralityString.Data(),fCollisionSystem.Data());

    Bool_t kAllKappaPlots = 0;
    Bool_t kUseFitConstraints = 1;
    if(isMC){
      kAllKappaPlots = 1;
    }

    Double_t value[4];
    Double_t error[4];
    Double_t valuesElEl[fNBinsPt];
    Double_t valuesRest[fNBinsPt];
    Double_t valuesPiPi[fNBinsPt];
    Double_t valuesElPi[fNBinsPt];
    Double_t ratiofitfractionElEl[fNBinsPt];
    Double_t ratiofitfractionRest[fNBinsPt];
    Double_t ratiofitfractionPiPi[fNBinsPt];
    Double_t ratiofitfractionElPi[fNBinsPt];
    Double_t errorsElEl[fNBinsPt];
    Double_t errorsRest[fNBinsPt];
    Double_t errorsPiPi[fNBinsPt];
    Double_t errorsElPi[fNBinsPt];
    Double_t fracElEl[fNBinsPt];
    Double_t fracRest[fNBinsPt];
    Double_t fracPiPi[fNBinsPt];
    Double_t fracElPi[fNBinsPt];

    //************************************ Read files **************************************************
    TFile fMC(fileMC.Data());
    TList *TopDirMC                             = (TList*)fMC.Get(nameMainDir.Data());
    if (!TopDirMC) {
        cout << "ERROR: TopDirMC not Found" << endl;
        return;
    }
    TList* HistosGammaConversionMC              = (TList*)TopDirMC->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if (!HistosGammaConversionMC) {
        cout << "ERROR: folder with Cutnumber " <<   fCutSelectionRead.Data() << " not contained in file " << fileMC.Data() << endl;
        return;
    }
    TList* ESDContainerMC                       = (TList*)HistosGammaConversionMC->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));

    TFile fData(fileData.Data());
    TList *TopDirData                             = (TList*)fData.Get(nameMainDir.Data());
    if (!TopDirData) {
        cout << "ERROR: TopDirData not Found" << endl;
        return;
    }
    TList* HistosGammaConversionData              = (TList*)TopDirData->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if (!HistosGammaConversionData) {
        cout << "ERROR: folder with Cutnumber " <<   fCutSelectionRead.Data() << " not contained in file " << fileData.Data() << endl;
        return;
    }
    TList* ESDContainerData                       = (TList*)HistosGammaConversionData->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));

    TH2F* hKappaTPCPtAfterCut;
    if(isMC) hKappaTPCPtAfterCut = (TH2F*)ESDContainerMC->FindObject("KappaTPC_Pt_after");
    else     hKappaTPCPtAfterCut = (TH2F*)ESDContainerData->FindObject("KappaTPC_Pt_after");

    TH2F* hKappaTPCPtElEl     = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp0_ee");
    TH2F* hKappaTPCPtPiPi     = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp1_pipi");
    TH2F* hKappaTPCPtElPi     = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp2_pie");
    TH2F* hKappaTPCPtPiK      = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp3_piK");
    TH2F* hKappaTPCPtPiP      = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp4_pip");
    TH2F* hKappaTPCPtElK      = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp5_eK");
    TH2F* hKappaTPCPtElP      = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp6_ep");
    TH2F* hKappaTPCPtKK       = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp7_KK");
    TH2F* hKappaTPCPtHad      = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp8_had");
    TH2F* hKappaTPCPtRest     = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp9_rem4");
    TH2F* hKappaTPCPtLeftoverRest = (TH2F*)ESDContainerMC->FindObject("hKappaTPC_Temp10_rem10");

    hKappaTPCAfterCut                       = new TH1D*[fNBinsPt];
    hKappaTPCElEl                           = new TH1D*[fNBinsPt];
    hKappaTPCPiPi                           = new TH1D*[fNBinsPt];
    hKappaTPCElPi                           = new TH1D*[fNBinsPt];
    hKappaTPCPiK                            = new TH1D*[fNBinsPt];
    hKappaTPCPiP                            = new TH1D*[fNBinsPt];
    hKappaTPCElK                            = new TH1D*[fNBinsPt];
    hKappaTPCElP                            = new TH1D*[fNBinsPt];
    hKappaTPCKK                             = new TH1D*[fNBinsPt];
    hKappaTPCHad                            = new TH1D*[fNBinsPt];
    hKappaTPCRest                           = new TH1D*[fNBinsPt];
    hKappaTPCLeftoverRest                   = new TH1D*[fNBinsPt];

    hKappaTPCSum                            = new TH1D*[fNBinsPt];
    hKappaTPCTotalSum                       = new TH1D*[fNBinsPt];
    hKappaTPCSumBkg                         = new TH1D*[fNBinsPt];

    hTemplateElEl                           = new TH1D*[fNBinsPt];
    hTemplateRest                           = new TH1D*[fNBinsPt];
    hTemplatePiPi                           = new TH1D*[fNBinsPt];
    hTemplateElPi                           = new TH1D*[fNBinsPt];
    hTemplateSum                            = new TH1D*[fNBinsPt];

    Double_t kappaRangeElEl1[] = { -3.0,   5.0}; //narrow
    Double_t kappaRangeElEl2[] = { -5.0,  10.0}; //wide
    Double_t kappaRangeElEl3[] = { -3.0,  10.0}; //asymm
    Double_t kappaRangePiPi[]  = {-20.0, -13.0};
    Double_t kappaRangeElPi[]  = {-11.0,  -6.0};
    Double_t kappaRangeTail[]  = { 11.0,  20.0};


    hValuesElEl                           = new TH1D("hValuesElEl", "", fNBinsPt, fBinsPt);
    hValuesElPi                           = new TH1D("hValuesElPi", "", fNBinsPt, fBinsPt);
    hValuesPiPi                           = new TH1D("hValuesPiPi", "", fNBinsPt, fBinsPt);
    hValuesRest                           = new TH1D("hValuesRest", "", fNBinsPt, fBinsPt);
    hFractionElEl                         = new TH1D("hFractionElEl", "", fNBinsPt, fBinsPt);
    hFractionElPi                         = new TH1D("hFractionElPi", "", fNBinsPt, fBinsPt);
    hFractionPiPi                         = new TH1D("hFractionPiPi", "", fNBinsPt, fBinsPt);
    hFractionRest                         = new TH1D("hFractionRest", "", fNBinsPt, fBinsPt);

    //************************************ Purity histogram ********************************************
    hSignalPurity1                           = new TH1D("hSignalPurity1", "", fNBinsPt, fBinsPt);
    hSignalPurity2                           = new TH1D("hSignalPurity2", "", fNBinsPt, fBinsPt);
    hSignalPurity3                           = new TH1D("hSignalPurity3", "", fNBinsPt, fBinsPt);
    hBackgroundPiPi1                         = new TH1D("hBackgroundPiPi1", "", fNBinsPt, fBinsPt);
    hBackgroundPiPi2                         = new TH1D("hBackgroundPiPi2", "", fNBinsPt, fBinsPt);
    hBackgroundPiPi3                         = new TH1D("hBackgroundPiPi3", "", fNBinsPt, fBinsPt);
    hBackgroundElPi1                         = new TH1D("hBackgroundElPi1", "", fNBinsPt, fBinsPt);
    hBackgroundElPi2                         = new TH1D("hBackgroundElPi2", "", fNBinsPt, fBinsPt);
    hBackgroundElPi3                         = new TH1D("hBackgroundElPi3", "", fNBinsPt, fBinsPt);
    hBackgroundRest1                         = new TH1D("hBackgroundRest1", "", fNBinsPt, fBinsPt);
    hBackgroundRest2                         = new TH1D("hBackgroundRest2", "", fNBinsPt, fBinsPt);
    hBackgroundRest3                         = new TH1D("hBackgroundRest3", "", fNBinsPt, fBinsPt);

  for(Int_t i = fStartPtBin; i<fNBinsPt; i++){

      Int_t startBin                      = hKappaTPCPtAfterCut->GetYaxis()->FindBin(fBinsPt[i]);
      Int_t endBin                        = hKappaTPCPtAfterCut->GetYaxis()->FindBin(fBinsPt[i+1]);
      cout << "projecting for bin range " << startBin << " - " << endBin << endl;
      cout << "pt range " << fBinsPt[i] << " - " << fBinsPt[i+1] << endl;
      hKappaTPCAfterCut[i] = (TH1D*)hKappaTPCPtAfterCut->ProjectionX(Form("KappaTPC_AfterCut_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCElEl[i]     = (TH1D*)hKappaTPCPtElEl->ProjectionX(Form("KappaTPC_ElEl_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCElPi[i]     = (TH1D*)hKappaTPCPtElPi->ProjectionX(Form("KappaTPC_ElPi_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCPiPi[i]     = (TH1D*)hKappaTPCPtPiPi->ProjectionX(Form("KappaTPC_PiPi_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCRest[i]     = (TH1D*)hKappaTPCPtRest->ProjectionX(Form("KappaTPC_Rest_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCElK[i]      = (TH1D*)hKappaTPCPtElK->ProjectionX(Form("KappaTPC_ElK_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCElP[i]      = (TH1D*)hKappaTPCPtElP->ProjectionX(Form("KappaTPC_ElP_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCPiP[i]      = (TH1D*)hKappaTPCPtPiP->ProjectionX(Form("KappaTPC_PiP_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCPiK[i]      = (TH1D*)hKappaTPCPtPiK->ProjectionX(Form("KappaTPC_PiK_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCKK[i]       = (TH1D*)hKappaTPCPtKK->ProjectionX(Form("KappaTPC_KK_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCHad[i]      = (TH1D*)hKappaTPCPtHad->ProjectionX(Form("KappaTPC_Had_Pt_%.2d",i),startBin,endBin,"e");
      hKappaTPCLeftoverRest[i] = (TH1D*)hKappaTPCPtLeftoverRest->ProjectionX(Form("KappaTPC_LeftoverRest_Pt_%.2d",i),startBin,endBin,"e");

      hKappaTPCSum[i]      = (TH1D*)hKappaTPCElEl[i]->Clone();
      hKappaTPCSum[i]->Add(hKappaTPCRest[i]);
      hKappaTPCSum[i]->Add(hKappaTPCPiPi[i]);
      hKappaTPCSum[i]->Add(hKappaTPCElPi[i]);

      hKappaTPCTotalSum[i] = (TH1D*)hKappaTPCElEl[i]->Clone();
      hKappaTPCTotalSum[i]->Add(hKappaTPCPiPi[i]);
      hKappaTPCTotalSum[i]->Add(hKappaTPCElPi[i]);
      hKappaTPCTotalSum[i]->Add(hKappaTPCPiK[i]);
      hKappaTPCTotalSum[i]->Add(hKappaTPCElK[i]);
      hKappaTPCTotalSum[i]->Add(hKappaTPCPiP[i]);
      hKappaTPCTotalSum[i]->Add(hKappaTPCElP[i]);
      hKappaTPCTotalSum[i]->Add(hKappaTPCKK[i]);
      hKappaTPCTotalSum[i]->Add(hKappaTPCHad[i]);
      hKappaTPCTotalSum[i]->Add(hKappaTPCLeftoverRest[i]);

      Double_t NTotal = hKappaTPCElEl[i]->GetEntries() + hKappaTPCRest[i]->GetEntries() + hKappaTPCPiPi[i]->GetEntries() + hKappaTPCElPi[i]->GetEntries();
      Double_t NEntriesAfterCut = hKappaTPCAfterCut[i]->GetEntries();

      if(isMC){
        fracElEl[i] = hKappaTPCElEl[i]->GetEntries()/NEntriesAfterCut;
        fracRest[i] = hKappaTPCRest[i]->GetEntries()/NEntriesAfterCut;
        fracPiPi[i] = hKappaTPCPiPi[i]->GetEntries()/NEntriesAfterCut;
        fracElPi[i] = hKappaTPCElPi[i]->GetEntries()/NEntriesAfterCut;

      } else {

        fracElEl[i] = hKappaTPCElEl[i]->GetEntries()/NTotal;
        fracRest[i] = hKappaTPCRest[i]->GetEntries()/NTotal;
        fracPiPi[i] = hKappaTPCPiPi[i]->GetEntries()/NTotal;
        fracElPi[i] = hKappaTPCElPi[i]->GetEntries()/NTotal;
      }

      //==================================================================================//
      //                       TEMPLATE FITTING and PLOTTING                              //
      //==================================================================================//

      //array of all MC that build up the data
      TObjArray *kappaArray = new TObjArray(4);
      kappaArray->Add(hKappaTPCElEl[i]); //primary + secondary
      kappaArray->Add(hKappaTPCRest[i]); //remaining
      kappaArray->Add(hKappaTPCPiPi[i]); //pion+pion
      kappaArray->Add(hKappaTPCElPi[i]); //pion+electron

      //configure the TFractionFitter
      TFractionFitter* fit = new TFractionFitter(hKappaTPCAfterCut[i], kappaArray);
      if(kUseFitConstraints){
        fit->Constrain(1,fracElEl[i]*constrainLow,fracElEl[i]*constrainHigh);
        fit->Constrain(2,fracRest[i]*constrainLow,fracRest[i]*constrainHigh);
        fit->Constrain(3,fracPiPi[i]*constrainLow,fracPiPi[i]*constrainHigh);
        fit->Constrain(4,fracElPi[i]*constrainLow,fracElPi[i]*constrainHigh);
      }

      //Fit the templates
      Int_t status = fit->Fit();
      //cout << "fit status: " << status << endl;

      //plot 6 pad projection
      TCanvas* canvas6PartKappaProj     = new TCanvas("canvas6PartKappaProj","",0,0,1200,800);
      DrawGammaCanvasSettings( canvas6PartKappaProj,  0.15, 0.05, 0.2, 0.15);
      canvas6PartKappaProj->Divide(3,2,0.001,0.001);

      TPaveText * pave = new TPaveText(0.15,0.7,0.45,0.93,"NDC");
      SetStylePave(pave,42,11,0.035);
      pave->InsertText(InfoSystem.Data());
      if(isMC) pave->InsertText(periodMC.Data());
      else pave->InsertText(periodData.Data());
      pave->InsertText(SigmaStarForm.Data());
      pave->InsertText(Form("%2.1f < #it{p}_{T} < %2.1f (GeV/#it{c})",fBinsPt[i],fBinsPt[i+1]));

      canvas6PartKappaProj->cd(1);
        TPaveText * paveSignal = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
        SetStylePave(paveSignal,62,11,0.04);
        paveSignal->InsertText("Signal");
        paveSignal->InsertText(Form("%f",hKappaTPCElEl[i]->GetEntries()/NTotal));

        HistoPlotSettings(hKappaTPCElEl[i],"K","Counts",0.00001,hKappaTPCElEl[i]->GetMaximum()*1.2,-20,20,12,12,3004);
        hKappaTPCElEl[i]->Draw("histo");
        pave->Draw();
        paveSignal->Draw();

      canvas6PartKappaProj->cd(2);
        TPaveText * paveElPi = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
        SetStylePave(paveElPi,62,11,0.04);
        paveElPi->InsertText("#pi^{#pm} + e^{#mp}");
        paveElPi->InsertText(Form("%f",hKappaTPCElPi[i]->GetEntries()/NTotal));

        HistoPlotSettings(hKappaTPCElPi[i],"K","Counts",0.00001,hKappaTPCElPi[i]->GetMaximum()*1.4,-20,20,kMagenta-2,kMagenta-2,3004);
        hKappaTPCElPi[i]->Draw("histo");
        paveElPi->Draw();

      canvas6PartKappaProj->cd(3);
        TPaveText * pavePiPi = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
        SetStylePave(pavePiPi,62,11,0.04);
        pavePiPi->InsertText("#pi^{#pm} + #pi^{#mp}");
        pavePiPi->InsertText(Form("%f",hKappaTPCPiPi[i]->GetEntries()/NTotal));

        HistoPlotSettings(hKappaTPCPiPi[i],"K","Counts",0.00001,hKappaTPCPiPi[i]->GetMaximum()*1.2,-20,20,kBlue-7,kBlue-7,3005);
        hKappaTPCPiPi[i]->Draw("histo");
        pavePiPi->Draw();

      canvas6PartKappaProj->cd(4);
        TPaveText * paveRest = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
        SetStylePave(paveRest,62,11,0.04);
        paveRest->InsertText("Remaining");
        paveRest->InsertText(Form("%f",hKappaTPCRest[i]->GetEntries()/NTotal));

        HistoPlotSettings(hKappaTPCRest[i],"K","Counts",0.00001,hKappaTPCRest[i]->GetMaximum()*1.2,-20,20,kRed+1,kRed+1,3005);
        hKappaTPCRest[i]->Draw("histo");
        paveRest->Draw();

      canvas6PartKappaProj->cd(5);
        canvas6PartKappaProj->cd(5)->SetLogy();

        HistoPlotSettings(hKappaTPCSum[i],"K","Counts",1,hKappaTPCSum[i]->GetMaximum()*10,-20,20,kBlack,kBlack,0);
        hKappaTPCSum[i]->SetLineWidth(2.0);
        hKappaTPCSum[i]->Draw("histo");
        hKappaTPCElEl[i]->Draw("histo,same");
        hKappaTPCPiPi[i]->Draw("histo,same");
        hKappaTPCElPi[i]->Draw("histo,same");
        hKappaTPCRest[i]->Draw("histo,same");

        TLegend* leg = new TLegend(0.7,0.7,0.9,0.93);
        leg->SetTextSize(0.04);
        leg->AddEntry(hKappaTPCSum[i],"total","lp");
        leg->AddEntry(hKappaTPCElEl[i],"real e^{+}e^{-}","f");
        leg->AddEntry(hKappaTPCRest[i],"remaining","f");
        leg->AddEntry(hKappaTPCPiPi[i],"#pi^{#pm} + #pi^{#mp}","f");
        leg->AddEntry(hKappaTPCElPi[i],"#pi^{#pm} + e^{#mp}","f");
        leg->SetBorderSize(0);
        leg->Draw();

      canvas6PartKappaProj->cd(6);
        HistoPlotSettings(hKappaTPCAfterCut[i],"K","Counts",0.00001,hKappaTPCAfterCut[i]->GetMaximum()*1.2,-20,20,kBlue,kBlue,0);
        hKappaTPCAfterCut[i]->Draw("histo");

        TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
        leg2->AddEntry(hKappaTPCAfterCut[i],"KappaTPC after cut","lp");

        if (status == 0) {
            TH1D *result = (TH1D*)fit->GetPlot();
            result->SetLineColor(kRed);
            result->Draw("hist,same");
            leg2->AddEntry(result,"fit","lp");
        }

        leg2->SetBorderSize(0);
        leg2->SetTextSize(0.04);
        leg2->Draw();

      canvas6PartKappaProj->SaveAs(Form("%s/KappaProj_%.2d.%s",fMainPlots.Data(),i,suffix.Data()));

      if(kAllKappaPlots){
          TCanvas* canvasAllKappaProj     = new TCanvas("canvasAllKappaProj","",0,0,1200,800);
          DrawGammaCanvasSettings( canvasAllKappaProj,  0.15, 0.05, 0.2, 0.15);
          canvasAllKappaProj->Divide(4,3,0.001,0.001);

          TPaveText * paveAll = new TPaveText(0.15,0.7,0.45,0.93,"NDC");
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
            paveElPi->Draw();

          canvasAllKappaProj->cd(3);
            HistoPlotSettings(hKappaTPCPiPi[i],"K","Counts",0.00001,hKappaTPCPiPi[i]->GetMaximum()*1.2,-20,20,kBlue-7,kBlue-7,3005);
            hKappaTPCPiPi[i]->Draw("histo");
            pavePiPi->Draw();

          canvasAllKappaProj->cd(4);
            TPaveText * pavePiK = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(pavePiK,62,11,0.04);
            pavePiK->InsertText("#pi^{#pm} + K^{#mp}");
            pavePiK->InsertText(Form("%f",hKappaTPCPiK[i]->GetEntries()/NTotal));

            HistoPlotSettings(hKappaTPCPiK[i],"K","Counts",0.00001,hKappaTPCPiK[i]->GetMaximum()*1.2,-20,20,kBlue+2,kBlue+2,3006);
            hKappaTPCPiK[i]->Draw("histo");
            pavePiK->Draw();

          canvasAllKappaProj->cd(5);
            TPaveText * pavePiP = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(pavePiP,62,11,0.04);
            pavePiP->InsertText("#pi^{#pm} + p(#bar{p})");
            pavePiP->InsertText(Form("%f",hKappaTPCPiP[i]->GetEntries()/NTotal));

            HistoPlotSettings(hKappaTPCPiP[i],"K","Counts",0.00001,hKappaTPCPiP[i]->GetMaximum()*1.2,-20,20,kOrange+1,kOrange+1,3007);
            hKappaTPCPiP[i]->Draw("histo");
            pavePiP->Draw();

          canvasAllKappaProj->cd(6);
            TPaveText * paveElK = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveElK,62,11,0.04);
            paveElK->InsertText("e^{#pm} + K^{#mp}");
            paveElK->InsertText(Form("%f",hKappaTPCElK[i]->GetEntries()/NTotal));

            HistoPlotSettings(hKappaTPCElK[i],"K","Counts",0.00001,hKappaTPCElK[i]->GetMaximum()*1.4,-20,20,kGreen-5,kGreen-5,3305);
            hKappaTPCElK[i]->Draw("histo");
            paveElK->Draw();

          canvasAllKappaProj->cd(7);
            TPaveText * paveElP = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveElP,62,11,0.04);
            paveElP->InsertText("e^{#pm} + p(#bar{p})");
            paveElP->InsertText(Form("%f",hKappaTPCElP[i]->GetEntries()/NTotal));

            HistoPlotSettings(hKappaTPCElP[i],"K","Counts",0.00001,hKappaTPCElP[i]->GetMaximum()*1.4,-20,20,kRed-4,kRed-4,3395);
            hKappaTPCElP[i]->Draw("histo");
            paveElP->Draw();

          canvasAllKappaProj->cd(8);
            TPaveText * paveKK = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveKK,62,11,0.04);
            paveKK->InsertText("K^{#pm} + K^{#mp}");
            paveKK->InsertText(Form("%f",hKappaTPCKK[i]->GetEntries()/NTotal));

            HistoPlotSettings(hKappaTPCKK[i],"K","Counts",0.00001,hKappaTPCKK[i]->GetMaximum()*1.4,-20,20,kCyan-2,kCyan-2,0);
            hKappaTPCKK[i]->Draw("histo");
            paveKK->Draw();

          canvasAllKappaProj->cd(9);
            TPaveText * paveHad = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveHad,62,11,0.04);
            paveHad->InsertText("Hadronic");
            paveHad->InsertText(Form("%f",hKappaTPCHad[i]->GetEntries()/NTotal));

            HistoPlotSettings(hKappaTPCHad[i],"K","Counts",0.00001,hKappaTPCHad[i]->GetMaximum()*1.2,-20,20,kGreen+2,kGreen+2,3005);
            hKappaTPCHad[i]->Draw("histo");
            paveHad->Draw();

          canvasAllKappaProj->cd(10);
            TPaveText * paveAllRest = new TPaveText(0.7,0.8,0.89,0.93,"NDC");
            SetStylePave(paveAllRest,62,11,0.04);
            paveAllRest->InsertText("Remaining");
            paveAllRest->InsertText(Form("%f",hKappaTPCLeftoverRest[i]->GetEntries()/NTotal));

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

            TLegend* legAll = new TLegend(0.2,0.2,0.9,0.9);
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
        hKappaTPCSum[i]->GetYaxis()->SetRangeUser(1.,hKappaTPCSum[i]->GetMaximum()*1e2);
        hKappaTPCSum[i]->Draw("histo");
        hKappaTPCElEl[i]->Draw("histo,same");
        hKappaTPCRest[i]->Draw("histo,same");
        hKappaTPCPiPi[i]->Draw("histo,same");
        hKappaTPCElPi[i]->Draw("histo,same");
        leg->Draw();
        pave->Draw();

      canvasProj->SaveAs(Form("%s/KappaProjSinglePlot_%.2d.%s",fMainPlots.Data(),i,suffix.Data()));

      if (status == 0) {
          for(int j=0; j<4; j++){
              fit->GetResult(j,value[j],error[j]);
          }

          valuesElEl[i] = value[0];
          errorsElEl[i] = error[0];
          valuesRest[i] = value[1];
          errorsRest[i] = error[1];
          valuesPiPi[i] = value[2];
          errorsPiPi[i] = error[2];
          valuesElPi[i] = value[3];
          errorsElPi[i] = error[3];

          cout << " value[0] = " << value[0] << endl;
          cout << " value[1] = " << value[1] << endl;
          cout << " value[2] = " << value[2] << endl;
          cout << " value[3] = " << value[3] << endl;

          cout << " fracElEl" << i << "] = " << fracElEl[i] << endl;
          cout << " fracRest[" << i << "] = " << fracRest[i] << endl;
          cout << " fracPiPi[" << i << "] = " << fracPiPi[i] << endl;
          cout << " fracElPi[" << i << "] = " << fracElPi[i] << endl;

          ratiofitfractionElEl[i]=valuesElEl[i]/fracElEl[i];
          ratiofitfractionRest[i]=valuesRest[i]/fracRest[i];
          ratiofitfractionPiPi[i]=valuesPiPi[i]/fracPiPi[i];
          ratiofitfractionElPi[i]=valuesElPi[i]/fracElPi[i];

          cout << "NEntriesAfterCut = " << NEntriesAfterCut << endl;
          cout << "NTotal = " << NTotal << endl;

          hValuesElEl->SetBinContent(i+1, valuesElEl[i]);
          hValuesElEl->SetBinError(i+1, errorsElEl[i]);
          hFractionElEl->SetBinContent(i+1, ratiofitfractionElEl[i]);
          hFractionElEl->SetBinError(i+1, 0);
          hValuesPiPi->SetBinContent(i+1, valuesPiPi[i]);
          hValuesPiPi->SetBinError(i+1, errorsPiPi[i]);
          hFractionPiPi->SetBinContent(i+1, ratiofitfractionPiPi[i]);
          hFractionPiPi->SetBinError(i+1, 0);
          hValuesElPi->SetBinContent(i+1, valuesElPi[i]);
          hValuesElPi->SetBinError(i+1, errorsElPi[i]);
          hFractionElPi->SetBinContent(i+1, ratiofitfractionElPi[i]);
          hFractionElPi->SetBinError(i+1, 0);
          hValuesRest->SetBinContent(i+1, valuesRest[i]);
          hValuesRest->SetBinError(i+1, errorsRest[i]);
          hFractionRest->SetBinContent(i+1, ratiofitfractionRest[i]);
          hFractionRest->SetBinError(i+1, 0);

          hTemplateElEl[i] = (TH1D*)hKappaTPCElEl[i]->Clone();
          hTemplateElEl[i]->Scale(value[0]/fracElEl[i]*(NEntriesAfterCut/NTotal));
          hTemplateRest[i] = (TH1D*)hKappaTPCRest[i]->Clone();
          hTemplateRest[i]->Scale(value[1]/fracRest[i]*(NEntriesAfterCut/NTotal));
          hTemplatePiPi[i] = (TH1D*)hKappaTPCPiPi[i]->Clone();
          hTemplatePiPi[i]->Scale(value[2]/fracPiPi[i]*(NEntriesAfterCut/NTotal));
          hTemplateElPi[i] = (TH1D*)hKappaTPCElPi[i]->Clone();
          hTemplateElPi[i]->Scale(value[3]/fracElPi[i]*(NEntriesAfterCut/NTotal));

          hTemplateSum[i] = (TH1D*)hTemplateElEl[i]->Clone();
          hTemplateSum[i]->GetXaxis()->SetTitleOffset(0.8);
          hTemplateSum[i]->SetFillColor(0);
          hTemplateSum[i]->SetLineColor(kBlack);
          hTemplateSum[i]->SetLineWidth(2.0);
          hTemplateSum[i]->Add(hTemplateRest[i]);
          hTemplateSum[i]->Add(hTemplatePiPi[i]);
          hTemplateSum[i]->Add(hTemplateElPi[i]);

          TCanvas* canvasTemplates = new TCanvas("canvasTemplates","",800,800);
          DrawGammaCanvasSettings( canvasTemplates,  0.12, 0.03, 0.03, 0.1);
            hTemplateSum[i]->Draw("histo");
            hTemplateElEl[i]->Draw("histo,same");
            hTemplateRest[i]->Draw("histo,same");
            hTemplatePiPi[i]->Draw("histo,same");
            hTemplateElPi[i]->Draw("histo,same");
            leg->Draw();
            pave->Draw();
          canvasTemplates->SaveAs(Form("%s/Template_%.2d.%s",fMainPlots.Data(),i,suffix.Data()));


          // calculate signal purity
          signalPurityNum1 = hTemplateElEl[i]->Integral(hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
          backgroundPiPiNum1 = hTemplatePiPi[i]->Integral(hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
          backgroundElPiNum1 = hTemplateElPi[i]->Integral(hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
          backgroundRestNum1 = hTemplateRest[i]->Integral(hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));

          signalPurityNum2 = hTemplateElEl[i]->Integral(hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
          backgroundPiPiNum2 = hTemplatePiPi[i]->Integral(hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
          backgroundElPiNum2 = hTemplateElPi[i]->Integral(hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
          backgroundRestNum2 = hTemplateRest[i]->Integral(hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));

          signalPurityNum3 = hTemplateElEl[i]->Integral(hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateElEl[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
          backgroundPiPiNum3 = hTemplatePiPi[i]->Integral(hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplatePiPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
          backgroundElPiNum3 = hTemplateElPi[i]->Integral(hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateElPi[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
          backgroundRestNum3 = hTemplateRest[i]->Integral(hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateRest[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));

          signalPurityDenom1 = hTemplateSum[i]->Integral(hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl1[0]),hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl1[1]));
          signalPurityDenom2 = hTemplateSum[i]->Integral(hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl2[0]),hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl2[1]));
          signalPurityDenom3 = hTemplateSum[i]->Integral(hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl3[0]),hTemplateSum[i]->GetXaxis()->FindBin(kappaRangeElEl3[1]));
          signalPurity1                    = signalPurityNum1 / signalPurityDenom1;
          backgroundPiPi1                  = backgroundPiPiNum1 / signalPurityDenom1;
          backgroundElPi1                  = backgroundElPiNum1 / signalPurityDenom1;
          backgroundRest1                  = backgroundRestNum1 / signalPurityDenom1;
          cout << "Purity for bin " << i << " in narrow range is " << signalPurity1 << endl;
          cout << "BG PiPi for bin " << i << " in narrow range is " << backgroundPiPi1 << endl;
          cout << "BG ElPi for bin " << i << " in narrow range is " << backgroundElPi1 << endl;
          cout << "BG rest for bin " << i << " in narrow range is " << backgroundRest1 << endl;

          signalPurity2                    = signalPurityNum2 / signalPurityDenom2;
          backgroundPiPi2                  = backgroundPiPiNum2 / signalPurityDenom2;
          backgroundElPi2                  = backgroundElPiNum2 / signalPurityDenom2;
          backgroundRest2                  = backgroundRestNum2 / signalPurityDenom2;
          cout << "Purity for bin " << i << " in wide range is " << signalPurity2 << endl;
          cout << "BG PiPi for bin " << i << " in wide range is " << backgroundPiPi2 << endl;
          cout << "BG ElPi for bin " << i << " in wide range is " << backgroundElPi2 << endl;
          cout << "BG rest for bin " << i << " in wide range is " << backgroundRest2 << endl;

          signalPurity3                    = signalPurityNum3 / signalPurityDenom3;
          backgroundPiPi3                  = backgroundPiPiNum3 / signalPurityDenom3;
          backgroundElPi3                  = backgroundElPiNum3 / signalPurityDenom3;
          backgroundRest3                  = backgroundRestNum3 / signalPurityDenom3;
          cout << "Purity for bin " << i << " in asymm range is " << signalPurity3 << endl;
          cout << "BG PiPi for bin " << i << " in asymm range is " << backgroundPiPi3 << endl;
          cout << "BG ElPi for bin " << i << " in asymm range is " << backgroundElPi3 << endl;
          cout << "BG rest for bin " << i << " in asymm range is " << backgroundRest3 << endl;

//             signalPurityErr                 = TMath::Sqrt(signalPurityNumErr*signalPurityNumErr + signalPurityDenomErr*signalPurityDenomErr);

          hSignalPurity1->SetBinContent(i+1, signalPurity1);
          hSignalPurity1->SetBinError(i+1, 0);
          hSignalPurity2->SetBinContent(i+1, signalPurity2);
          hSignalPurity2->SetBinError(i+1, 0);
          hSignalPurity3->SetBinContent(i+1, signalPurity3);
          hSignalPurity3->SetBinError(i+1, 0);

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

          hBackgroundRest1->SetBinContent(i+1, backgroundRest1);
          hBackgroundRest1->SetBinError(i+1, 0);
          hBackgroundRest2->SetBinContent(i+1, backgroundRest2);
          hBackgroundRest2->SetBinError(i+1, 0);
          hBackgroundRest3->SetBinContent(i+1, backgroundRest3);
          hBackgroundRest3->SetBinError(i+1, 0);

      } //end status == 0 if


      delete kappaArray;
      delete fit;

  } //end of pt loop

    TCanvas* canvasPurity     = new TCanvas("canvasPurity","",0,0,1000,800);  // gives the page size
    DrawGammaCanvasSettings( canvasPurity,  0.07, 0.02, 0.02, 0.09);
    TH2D *histo2DPurity = new TH2D("histo2DPurity", "histo2DPurity", 1000,0.,14.,1000,0.7,1.01);
    SetStyleHistoTH2ForGraphs(histo2DPurity, "#it{p}_{T} (GeV/#it{c})","Purity", 0.035,0.04, 0.035,0.04, 1.,.8);
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

    TCanvas* canvasBackground     = new TCanvas("canvasBackground","",0,0,1500,600);  // gives the page size
    DrawGammaCanvasSettings( canvasBackground,  0.07, 0.0, 0.02, 0.09);
    canvasBackground->Divide(3,1,0.000001,0.00001);
    TH2D *histo2DBackground = new TH2D("histo2DBackground", "histo2DBackground", 1000,0.,14.,10,-.001,.1);
    SetStyleHistoTH2ForGraphs(histo2DBackground, "#it{p}_{T} (GeV/#it{c})","Background", 0.035,0.04, 0.035,0.04, 1.,1.4);

      DrawGammaSetMarker(hBackgroundElPi1,25,1.,kMagenta-2,kMagenta-2);
      DrawGammaSetMarker(hBackgroundPiPi1,24,1.,kBlue-7,kBlue-7);
      DrawGammaSetMarker(hBackgroundPiPi2,24,1.,kBlue-7,kBlue-7);
      DrawGammaSetMarker(hBackgroundElPi2,25,1.,kMagenta-2,kMagenta-2);
      DrawGammaSetMarker(hBackgroundElPi3,25,1.,kMagenta-2,kMagenta-2);
      DrawGammaSetMarker(hBackgroundPiPi3,24,1.,kBlue-7,kBlue-7);
      DrawGammaSetMarker(hBackgroundRest1,27,1.,kRed,kRed);
      DrawGammaSetMarker(hBackgroundRest2,27,1.,kRed,kRed);
      DrawGammaSetMarker(hBackgroundRest3,27,1.,kRed,kRed);

    canvasBackground->cd(1);
    histo2DBackground->DrawCopy();

      hBackgroundPiPi1->Draw("same,p");
      hBackgroundElPi1->Draw("same,p");
      hBackgroundRest1->Draw("same,p");

      TLegend* legBackground1 = new TLegend(0.2,0.7,0.5,0.9);
      legBackground1->SetTextSize(0.04);
      legBackground1->SetHeader(Form("Contamination in %.2f < K < %.2f",kappaRangeElEl1[0],kappaRangeElEl1[1]));
      legBackground1->AddEntry(hBackgroundPiPi1,"#pi^{#pm}#pi^{#mp}","p");
      legBackground1->AddEntry(hBackgroundElPi1,"e^{#pm}#pi^{#mp}","p");
      legBackground1->AddEntry(hBackgroundRest1,"rest","p");
      legBackground1->SetBorderSize(0);
      legBackground1->Draw();

    canvasBackground->cd(2);
    histo2DBackground->DrawCopy();

      hBackgroundElPi2->Draw("same,p");
      hBackgroundPiPi2->Draw("same,p");
      hBackgroundRest2->Draw("same,p");

      TLegend* legBackground2 = new TLegend(0.2,0.7,0.5,0.9);
      legBackground2->SetTextSize(0.04);
      legBackground2->SetHeader(Form("Contamination in %.2f < K < %.2f",kappaRangeElEl2[0],kappaRangeElEl2[1]));
      legBackground2->AddEntry(hBackgroundPiPi2,"#pi^{#pm}#pi^{#mp}","p");
      legBackground2->AddEntry(hBackgroundElPi2,"e^{#pm}#pi^{#mp}","p");
      legBackground2->AddEntry(hBackgroundRest2,"rest","p");
      legBackground2->SetBorderSize(0);
      legBackground2->Draw();

    canvasBackground->cd(3);
    histo2DBackground->DrawCopy();

      hBackgroundPiPi3->Draw("same,p");
      hBackgroundElPi3->Draw("same,p");
      hBackgroundRest3->Draw("same,p");

      TLegend* legBackground3 = new TLegend(0.2,0.7,0.5,0.9);
      legBackground3->SetTextSize(0.04);
      legBackground3->SetHeader(Form("Contamination in %.2f < K < %.2f",kappaRangeElEl3[0],kappaRangeElEl3[1]));
      legBackground3->AddEntry(hBackgroundPiPi3,"#pi^{#pm}#pi^{#mp}","p");
      legBackground3->AddEntry(hBackgroundElPi3,"e^{#pm}#pi^{#mp}","p");
      legBackground3->AddEntry(hBackgroundRest3,"rest","p");
      legBackground3->SetBorderSize(0);
      legBackground3->Draw();

    canvasBackground->SaveAs(Form("%s/Background.%s",fOutputDir.Data(),suffix.Data()));


    TCanvas* canvasFractions     = new TCanvas("canvasFractions","",0,0,1000,500);  // gives the page size
    DrawGammaCanvasSettings( canvasFractions,  0.07, 0.02, 0.02, 0.09);
    canvasFractions->Divide(2,1,0.0001,0.0001);

    canvasFractions->cd(1);
        HistoMarkerPlotSettings(hValuesElEl,"#it{p}_{T} (GeV/#it{c})","Template fraction",0.,1.1,0.,14.,12,12,20);
        HistoMarkerPlotSettings(hValuesElPi,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kMagenta-2,kMagenta-2,20);
        HistoMarkerPlotSettings(hValuesPiPi,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kBlue-7,kBlue-7,20);
        HistoMarkerPlotSettings(hValuesRest,"#it{p}_{T} (GeV/#it{c})","",0.,1.1,0.,14.,kRed+1,kRed+1,20);
        hValuesElEl->Draw("p");
        hValuesPiPi->Draw("p,same");
        hValuesElPi->Draw("p,same");
        hValuesRest->Draw("p,same");

      TLegend* legFracTemplate = new TLegend(0.5,0.35,0.9,0.55);
      legFracTemplate->SetTextSize(0.04);
      legFracTemplate->AddEntry(hValuesElEl,"e^{#pm}e^{#mp}","p");
      legFracTemplate->AddEntry(hValuesPiPi,"#pi^{#pm}#pi^{#mp}","p");
      legFracTemplate->AddEntry(hValuesElPi,"e^{#pm}#pi^{#mp}","p");
      legFracTemplate->AddEntry(hValuesRest,"rest","p");
      legFracTemplate->SetBorderSize(0);
      legFracTemplate->Draw();
      TPaveText * pave1 = new TPaveText(0.5,0.56,0.9,0.7,"NDC");
      SetStylePave(pave1,42,11,0.035);
      pave1->InsertText(InfoSystem.Data());
      if(isMC) pave1->InsertText(periodMC.Data());
      else pave1->InsertText(periodData.Data());
      pave1->Draw();

    canvasFractions->cd(2);
        HistoMarkerPlotSettings(hFractionElEl,"#it{p}_{T} (GeV/#it{c})","fit/histo",0.5,1.5,0.,14.,12,12,20);
        HistoMarkerPlotSettings(hFractionElPi,"#it{p}_{T} (GeV/#it{c})","",0.5,1.5,0.,14.,kMagenta-2,kMagenta-2,20);
        HistoMarkerPlotSettings(hFractionPiPi,"#it{p}_{T} (GeV/#it{c})","",0.5,1.5,0.,14.,kBlue-7,kBlue-7,20);
        HistoMarkerPlotSettings(hFractionRest,"#it{p}_{T} (GeV/#it{c})","",0.5,1.5,0.,14.,kRed+1,kRed+1,20);
        hFractionElEl->Draw("p");
        hFractionElPi->Draw("p,same");
        hFractionPiPi->Draw("p,same");
        hFractionRest->Draw("p,same");


    canvasFractions->SaveAs(Form("%s/Fractions.%s",fOutputDir.Data(),suffix.Data()));



//     SaveHistos(cutSelection);

}



// *****************************************************************************************************
// *********************** Initialize histograms and binning *******************************************
// *****************************************************************************************************
void Initialize(TString setPi0, TString energy , Int_t numberOfBins, Int_t mode, Bool_t addSig){

  InitializeBinning(setPi0, numberOfBins, energy, fDirectPhoton, 0, fEventCutNumber, fClusterCutNumber);

}


// *****************************************************************************************************
// *********************** Save histograms to file *****************************************************
// *****************************************************************************************************
void SaveHistos(TString cutString){
    const char* nameOutput  = Form("%s/%s/PurityFromTemplate_%s.root",cutString.Data(),fEnergyFlag.Data(), cutString.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *file             = new TFile(nameOutput,"UPDATE");

    // write kappa projections in pt bins
    for (Int_t i=0; i<fNBinsPt; i++) {
        if (hKappaTPCAfterCut[i])   hKappaTPCAfterCut[i]->Write(hKappaTPCAfterCut[i]->GetName(),TObject::kOverwrite);
        if (hKappaTPCElEl[i])       hKappaTPCElEl[i]->Write(hKappaTPCElEl[i]->GetName(),TObject::kOverwrite);
        if (hKappaTPCElPi[i])       hKappaTPCElPi[i]->Write(hKappaTPCElPi[i]->GetName(),TObject::kOverwrite);
        if (hKappaTPCPiPi[i])       hKappaTPCPiPi[i]->Write(hKappaTPCPiPi[i]->GetName(),TObject::kOverwrite);
        if (hKappaTPCRest[i])       hKappaTPCRest[i]->Write(hKappaTPCRest[i]->GetName(),TObject::kOverwrite);
    }

    // write templates
    for (Int_t i=0; i<fNBinsPt; i++) {
        if (hTemplateElEl[i])   hTemplateElEl[i]->Write(hTemplateElEl[i]->GetName(),TObject::kOverwrite);
        if (hTemplateElPi[i])   hTemplateElPi[i]->Write(hTemplateElPi[i]->GetName(),TObject::kOverwrite);
        if (hTemplatePiPi[i])   hTemplatePiPi[i]->Write(hTemplatePiPi[i]->GetName(),TObject::kOverwrite);
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
                       Style_t fillStyle
//                        Width_t lineWidth
//                        Style_t lineStyle = 1,
                       )
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

//     histo1->DrawCopy("hist");
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
                       Style_t markerStyle
//                        Width_t lineWidth
//                        Style_t lineStyle = 1,
                       )
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

//     histo1->DrawCopy("hist");
}

void SetStylePave(TPaveText *pave,
                 Style_t textFontStyle,
                 Style_t textAlign,
                 Size_t textSize
                 )
{
        pave->SetTextFont(textFontStyle);
        pave->SetTextSize(textSize);
        pave->SetBorderSize(0);
        pave->SetFillColor(kWhite);
        pave->SetFillStyle(0);
        pave->SetTextAlign(textAlign);


}