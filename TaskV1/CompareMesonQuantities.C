//  **********************************************************************************
//  ******     provided by Gamma Conversion Group, PWGGA,                        *****
//  ******     Friederike Bock, friederike.bock@cern.ch                          *****
//  **********************************************************************************

#include <stdlib.h>
#include <iostream>
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
#include "TPad.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TMarker.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
// #include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/PlottingMeson.h"

using namespace std;
Bool_t debugMesonQual = kFALSE;

//**********************************************************************************
//******************* return minimum for 1D histo  *********************************
//**********************************************************************************
Double_t FindSmallestEntryIn1D(TH1* histo){
    Double_t minimum = 1;
    if (histo){
        for (Int_t i = 1; i<histo->GetNbinsX(); i++){
            if (histo->GetBinContent(i) < minimum ){
                minimum = histo->GetBinContent(i);
            }
        }
    }
    return minimum;
}

void Delete(){
    if (fBinsPt)                                                delete[] fBinsPt;
    if (fNRebin)                                                delete fNRebin;
}

TF1* FitSubtractedInvMassLinRemBG(TH1D* histoMappingSignalInvMassPtBinSingle,
                                  Double_t* mesonMassFitRange,
                                  Double_t* mesonMassPlotRange,
                                  Double_t& chi2ndf,
                                  Int_t ptBin,
                                  Bool_t debug = kFALSE)
{

    if(debug) cout << "FitSubtractedInvMassLinRemBG for pT bin " << ptBin << endl;

    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(mesonMassPlotRange[0],mesonMassPlotRange[1]);

    TF1 * fFitRemBck = NULL;
    fFitRemBck = new TF1("Linear","[0]+[1]*x",mesonMassFitRange[0],mesonMassFitRange[1]);

    histoMappingSignalInvMassPtBinSingle->Fit(fFitRemBck,"QRME0");
    Int_t fitStatus = -1;
    if(debug) fitStatus = histoMappingSignalInvMassPtBinSingle->Fit(fFitRemBck,"RME0");
    else fitStatus = histoMappingSignalInvMassPtBinSingle->Fit(fFitRemBck,"QRME0");

    if(debug){
        cout << "Fit status:  \t" << fitStatus << endl;
        cout << "Parameter 0: \t" << fFitRemBck->GetParameter(0) <<"+-" << fFitRemBck->GetParError(0) << endl;
        cout << "Parameter 1: \t" << fFitRemBck->GetParameter(1) <<"+-" << fFitRemBck->GetParError(1) << endl;
        cout << "Chi2/ndf   : \t" << fFitRemBck->GetChisquare()/fFitRemBck->GetNDF() << endl;
    }

    if(fitStatus>-1){ // fit OK
        chi2ndf = fFitRemBck->GetChisquare()/fFitRemBck->GetNDF();
        fFitRemBck->SetLineColor(kCyan+2);
        fFitRemBck->SetLineWidth(1);
        fFitRemBck->SetLineStyle(2);
        fFitRemBck->SetNpx(10000);
    } else {
        fFitRemBck->SetParameter(0,0);
        fFitRemBck->SetParameter(1,0);
        chi2ndf = 0.0;
    }

    return fFitRemBck;

}

TF1* FitSubtractedInvMassPol2RemBG(TH1D* histoMappingSignalInvMassPtBinSingle,
                                   Double_t* mesonMassFitRange,
                                   Double_t* mesonMassPlotRange,
                                   Double_t& chi2ndf,
                                   Int_t ptBin,
                                   Bool_t debug = kFALSE)
{

    if(debug) cout << "FitSubtractedInvMassPol2RemBG for pT bin " << ptBin << endl;

    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(mesonMassPlotRange[0],mesonMassPlotRange[1]);

    TF1* fFitRemBck   = NULL;
    fFitRemBck   = new TF1("BGfitPol2","[0]+[1]*x+[2]*x*x",mesonMassFitRange[0],mesonMassFitRange[1]);

    histoMappingSignalInvMassPtBinSingle->Fit(fFitRemBck,"QRME0");
    Int_t fitStatus = -1;
    if(debug) fitStatus = histoMappingSignalInvMassPtBinSingle->Fit(fFitRemBck,"RME0");
    else fitStatus = histoMappingSignalInvMassPtBinSingle->Fit(fFitRemBck,"QRME0");

    if(debug){
        cout << "Fit status:  \t" << fitStatus << endl;
        cout << "Parameter 0: \t" << fFitRemBck->GetParameter(0) <<"+-" << fFitRemBck->GetParError(0) << endl;
        cout << "Parameter 1: \t" << fFitRemBck->GetParameter(1) <<"+-" << fFitRemBck->GetParError(1) << endl;
        cout << "Parameter 2: \t" << fFitRemBck->GetParameter(2) <<"+-" << fFitRemBck->GetParError(2) << endl;
        cout << "Chi2/ndf   : \t" << fFitRemBck->GetChisquare()/fFitRemBck->GetNDF() << endl;
    }

    if(fitStatus>-1){ // fit OK
        chi2ndf = fFitRemBck->GetChisquare()/fFitRemBck->GetNDF();
        fFitRemBck->SetLineColor(kRed+1);
        fFitRemBck->SetLineWidth(1);
        fFitRemBck->SetLineStyle(2);
        fFitRemBck->SetNpx(10000);
    } else {
        fFitRemBck->SetParameter(0,0);
        fFitRemBck->SetParameter(1,0);
        fFitRemBck->SetParameter(2,0);
        chi2ndf = 0.0;
    }

    return fFitRemBck;

}

TF1* FitSubtractedInvMassPol3RemBG(TH1D* histoMappingSignalInvMassPtBinSingle,
                                   Double_t* mesonMassFitRange,
                                   Double_t* mesonMassPlotRange,
                                   Double_t& chi2ndf,
                                   Int_t ptBin,
                                   Bool_t debug = kFALSE)
{

    if(debug) cout << "FitSubtractedInvMassPol3RemBG for pT bin " << ptBin << endl;

    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(mesonMassPlotRange[0],mesonMassPlotRange[1]);

    TF1* fFitRemBck   = NULL;
    fFitRemBck   = new TF1("BGfitPol3","[0]+[1]*x+[2]*x*x+[3]*x*x*x",mesonMassFitRange[0],mesonMassFitRange[1]);

    histoMappingSignalInvMassPtBinSingle->Fit(fFitRemBck,"QRME0");
    Int_t fitStatus = -1;
    if(debug) fitStatus = histoMappingSignalInvMassPtBinSingle->Fit(fFitRemBck,"RME0");
    else fitStatus = histoMappingSignalInvMassPtBinSingle->Fit(fFitRemBck,"QRME0");

    if(debug){
        cout << "Fit status:  \t" << fitStatus << endl;
        cout << "Parameter 0: \t" << fFitRemBck->GetParameter(0) <<"+-" << fFitRemBck->GetParError(0) << endl;
        cout << "Parameter 1: \t" << fFitRemBck->GetParameter(1) <<"+-" << fFitRemBck->GetParError(1) << endl;
        cout << "Parameter 2: \t" << fFitRemBck->GetParameter(2) <<"+-" << fFitRemBck->GetParError(2) << endl;
        cout << "Parameter 2: \t" << fFitRemBck->GetParameter(3) <<"+-" << fFitRemBck->GetParError(3) << endl;
        cout << "Chi2/ndf   : \t" << fFitRemBck->GetChisquare()/fFitRemBck->GetNDF() << endl;
    }

    if(fitStatus>-1){ // fit OK
        chi2ndf = fFitRemBck->GetChisquare()/fFitRemBck->GetNDF();
        fFitRemBck->SetLineColor(kBlack);
        fFitRemBck->SetLineWidth(1);
        fFitRemBck->SetLineStyle(2);
        fFitRemBck->SetNpx(10000);
    } else {
        fFitRemBck->SetParameter(0,0);
        fFitRemBck->SetParameter(1,0);
        fFitRemBck->SetParameter(2,0);
        fFitRemBck->SetParameter(3,0);
        chi2ndf = 0.0;
    }

    return fFitRemBck;

}

void CompareMesonQuantities(    const char *dataFilename        = "rawSignalData",
                                const char *mcFilename          = "rawSignalMC",
                                TString fCutSelection           = "",
                                TString mesonType               = "Pi0",
                                TString fSuffix                 = "",
                                TString energyFlag              = "" ,
                                TString directPhoton            = "",
                                Int_t numberOfBins              = 25,
                                Int_t mode                      = 0
                            )
{
    gROOT->Reset();
    // mode:    0 // new output PCM-PCM
    //          1 // new output PCM dalitz
    //          2 // new output PCM-Calo
    //          3 // new output Calo-Calo
    //          4 // new output EMCAL-EMCAL
    //          5 // new output PHOS-PHOS
    //          9 // old output PCM-PCM


    StyleSettingsThesis(fSuffix);
    SetPlotStyle();

    TString DetectionChannel    = ReturnFullTextReconstructionProcess(mode);
    TString date                = ReturnDateString();
    TString textAlice           = "ALICE performance";
    TString textProcess         = ReturnMesonString (mesonType);
    TString decayChannel        = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
    TString energyText          = ReturnFullCollisionsSystem(energyFlag);
    if (energyText.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    if(textProcess.CompareTo("") == 0 ){
        cout << "Meson unknown" << endl;
        return ;
    }

    TFile* fileRawSignalData = new TFile(dataFilename, "READ");
    TFile* fileRawSignalMC = new TFile(mcFilename, "READ");


    cout << "dataFilename: " << dataFilename << endl;
    cout << "mcFilename: " << mcFilename << endl;
    cout << "fCutSelection: " << fCutSelection.Data() << endl;
    cout << "mesonType: " << mesonType.Data() << endl;
    cout << "fSuffix: " << fSuffix.Data()<< endl;
    cout << "energyFlag: " << energyFlag.Data() << endl;
    cout << "numberOfBins: " << numberOfBins << endl;

    TH1D* histoChi2Data                     = (TH1D*) fileRawSignalData->Get("histoChi2_0");
    TH1D* histoChi2_Pol2_Data               = (TH1D*) fileRawSignalData->Get("histoChi2_1");
    TH1D* histoChi2_Exp1_Data               = (TH1D*) fileRawSignalData->Get("histoChi2_2");
    TH1D* histoChi2_Exp2_Data               = (TH1D*) fileRawSignalData->Get("histoChi2_3");
    TH1D* histoConstResBGData               = (TH1D*) fileRawSignalData->Get("histoResidualBGcon");
    TH1D* histoLinResBGData                 = (TH1D*) fileRawSignalData->Get("histoResidualBGlin");
    TH1D* histoResBGYieldVsTotBGData        = (TH1D*) fileRawSignalData->Get("histoRatioResBGYield");
    TH1D* histoResBGYieldVsResBGPlSigData   = (TH1D*) fileRawSignalData->Get("histoRatioResBGYieldToSPlusResBG");
    TH1D* histoLambdaTailData               = (TH1D*) fileRawSignalData->Get("histoLambdaTail");
    TH1D* histoChi2MC                       = (TH1D*) fileRawSignalMC->Get("histoChi2_0");
    TH1D* histoChi2_Pol2_MC                 = (TH1D*) fileRawSignalMC->Get("histoChi2_1");
    TH1D* histoChi2_Exp1_MC                 = (TH1D*) fileRawSignalMC->Get("histoChi2_2");
    TH1D* histoChi2_Exp2_MC                 = (TH1D*) fileRawSignalMC->Get("histoChi2_3");
    TH1D* histoConstResBGMC                 = (TH1D*) fileRawSignalMC->Get("histoResidualBGcon");
    TH1D* histoLinResBGMC                   = (TH1D*) fileRawSignalMC->Get("histoResidualBGlin");
    TH1D* histoResBGYieldVsTotBGMC          = (TH1D*) fileRawSignalMC->Get("histoRatioResBGYield");
    TH1D* histoResBGYieldVsResBGPlSigMC     = (TH1D*) fileRawSignalMC->Get("histoRatioResBGYieldToSPlusResBG");
    TH1D* histoLambdaTailMC                 = (TH1D*) fileRawSignalMC->Get("histoLambdaTail");
    TH1D* histoRatioYieldLowMassDataDivMC   = (TH1D*) histoChi2Data->Clone("histoRatioYieldLowMassDataDivMC");

    TString outputDir                       = Form("%s/%s/%s/ExtractSignal",fCutSelection.Data(),energyFlag.Data(),fSuffix.Data());
    TString nameLineShapePlot               = Form("%s/%s_MesonLineShapeCompared_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data());
    TString nameLineShapePlotLeft           = Form("%s/%s_MesonLineShapeComparedLeft_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data());
    TString nameLineShapePlotFull           = Form("%s/%s_MesonCompLineShapeWithBG_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data());

    TString fEventCutSelection              = "";
    TString fGammaCutSelection              = "";
    TString fClusterCutSelection            = "";
    TString fElectronCutSelection           = "";
    TString fMesonCutSelection              = "";
    ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
    InitializeBinning(mesonType, numberOfBins, energyFlag, directPhoton, mode, fEventCutSelection, fClusterCutSelection, -1, kFALSE, "", "", fGammaCutSelection, kFALSE);

    Double_t fMesonRange[2]    = {0., 0.79};
    Double_t peakRange[2]      = {0.10,0.15};
    Double_t fMesonFitRange[2] = {0.02, 0.3};
    TF1* fitLinRemBG[200]      = {NULL};
    TF1* fitPol2RemBG[200]     = {NULL};
    TF1* fitPol3RemBG[200]     = {NULL};
    Double_t chi2ndfLin[200];
    Double_t chi2ndfPol2[200];
    Double_t chi2ndfPol3[200];
    TH1D* histoChi2RemLinBGFit  = NULL;
    TH1D* histoChi2RemPol2BGFit = NULL;
    TH1D* histoChi2RemPol3BGFit = NULL;

    if (mesonType.CompareTo("Pi0") == 0 || mesonType.CompareTo("Pi0EtaBinning") == 0){
        fMesonRange[0]      = 0.;
        fMesonRange[1]      = 0.3;
    } else if (mesonType.CompareTo("Eta") == 0){
        fMesonRange[0]      = 0.35;
        fMesonRange[1]      = 0.79;
        peakRange[0]        = 0.50;
        peakRange[0]        = 0.60;
        fMesonFitRange[0]   = fMesonRange[0];
        fMesonFitRange[1]   = fMesonRange[1];
    } else if (mesonType.CompareTo("EtaPrime") == 0){
        fMesonRange[0]      = 0.9;
        fMesonRange[1]      = 1.0;
        peakRange[0]        = 0.90;
        peakRange[0]        = 1.00;
    }

    cout << "Start bin for " << mesonType << " (mode " << mode << "): " << fStartPtBin << endl;

    //******************************* Reading histograms **************************************************************
    if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<<"initializing hists" << endl;
    TH1D*  histoSignalDataInvMassPtBin[200] = {NULL};
    TH1D*  histoSignalMCInvMassPtBin[200] = {NULL};
    TH1D*  histoTrueMCInvMassPtBin[200] = {NULL};
    TH1D*  histoMCRemBGInvMassPtBin[200] = {NULL};
    if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "initializing arrays" << endl;
    Double_t intLowMassData[200]                    = {0.};
    Double_t intLowMassMC[200]                      = {0.};
    Double_t intErrLowMassData[200]                 = {0.};
    Double_t intErrLowMassMC[200]                   = {0.};
    Double_t ratioLowMass[200]                      = {1.};
    Double_t ratioErrLowMass[200]                   = {1.};
    for(Int_t j=0;j<3;j++){
        if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "Beginning of Loop j = " << j << " fNBinsPt = " << fNBinsPt << endl;
        TString histonameSignal;
        TString histonameMCTruth;
        for(Int_t iPt=fStartPtBin; iPt<fNBinsPt && iPt < 200; iPt++){

            histoSignalDataInvMassPtBin[iPt]    = NULL;
            histoSignalMCInvMassPtBin[iPt]      = NULL;
            histoTrueMCInvMassPtBin[iPt]        = NULL;
            histoMCRemBGInvMassPtBin[iPt]       = NULL;

            Double_t startPt    = fBinsPt[iPt];
            Double_t endPt      = fBinsPt[iPt+1];
            Double_t ptValue    = (fBinsPt[iPt] + fBinsPt[iPt+1])/2.;
            if(j==0){
                histonameSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", iPt);  // normal
            } else if(j==1) {
                histonameSignal = Form("fHistoMappingSignalInvMassLeft_in_Pt_Bin%02d", iPt);  // left normalization
            } else if(j==2) {
                histonameSignal = Form("Mapping_GG_InvMass_in_Pt_Bin%02d", iPt); //with combinatorial BG
            }
            if (!(fileRawSignalData && fileRawSignalMC)){
                cout << "lost file info" << endl;
                continue;
            }
            if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "loading data " << fileRawSignalData << endl;
            histoSignalDataInvMassPtBin[iPt]    = (TH1D*)fileRawSignalData->Get(histonameSignal);
            if(!histoSignalDataInvMassPtBin[iPt] ) continue;
            if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "loading MC " << fileRawSignalMC << endl;
            histoSignalMCInvMassPtBin[iPt]      = (TH1D*)fileRawSignalMC->Get(histonameSignal);
            if(!histoSignalMCInvMassPtBin[iPt])   continue;
            histoMCRemBGInvMassPtBin[iPt]       = (TH1D*)histoSignalMCInvMassPtBin[iPt]->Clone();

            Double_t firstBinIntData = histoSignalDataInvMassPtBin[iPt]->FindBin(fMesonRange[0]+0.0001);
            Double_t lastBinIntData = histoSignalDataInvMassPtBin[iPt]->FindBin(fMesonRange[1]-0.0001);
            Double_t firstBinIntMC = histoSignalMCInvMassPtBin[iPt]->FindBin(fMesonRange[0]+0.0001);
            Double_t lastBinIntMC = histoSignalMCInvMassPtBin[iPt]->FindBin(fMesonRange[1]-0.0001);
            if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "iPt = " << iPt << " firstBinIntData= " << firstBinIntData << "\t lastBinIntData= " << lastBinIntData << endl;
            if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "iPt = " << iPt << " firstBinIntMC= " << firstBinIntMC << "\t lastBinIntMC= " << lastBinIntMC << endl;
            Double_t integralData = histoSignalDataInvMassPtBin[iPt]->Integral(firstBinIntData,lastBinIntData);
            Double_t integralMC   = histoSignalMCInvMassPtBin[iPt]->Integral(firstBinIntMC,lastBinIntMC);
            if (j == 2){
                integralData = histoSignalDataInvMassPtBin[iPt]->Integral(histoSignalDataInvMassPtBin[iPt]->FindBin(peakRange[0]+0.0001),histoSignalDataInvMassPtBin[iPt]->FindBin(peakRange[1]-0.0001) );
                integralMC   = histoSignalMCInvMassPtBin[iPt]->Integral(histoSignalMCInvMassPtBin[iPt]->FindBin(peakRange[0]+0.0001),histoSignalMCInvMassPtBin[iPt]->FindBin(peakRange[1]-0.0001) );
            }

            if (integralData < 0 || integralMC < 0 || j == 1){  // normalize with maximum
                integralData                    = histoSignalDataInvMassPtBin[iPt]->GetMaximum();
                integralMC                      = histoSignalMCInvMassPtBin[iPt]->GetMaximum();
            }
            if (j == 2){
                intLowMassData[iPt]             = histoSignalDataInvMassPtBin[iPt]->Integral( 1, histoSignalDataInvMassPtBin[iPt]->FindBin(0.05));
                intErrLowMassData[iPt]          = TMath::Sqrt(intLowMassData[iPt]);
                intLowMassMC[iPt]               = histoSignalMCInvMassPtBin[iPt]->Integral( 1, histoSignalMCInvMassPtBin[iPt]->FindBin(0.05));
                intErrLowMassMC[iPt]            = TMath::Sqrt(intLowMassMC[iPt]);
                if (intLowMassMC[iPt] != 0){
                    ratioLowMass[iPt]           = (intLowMassData[iPt]/integralData)/(intLowMassMC[iPt]/integralMC);
                    ratioErrLowMass[iPt]        = ratioLowMass[iPt]* TMath::Sqrt(TMath::Power(intErrLowMassData[iPt]/intLowMassData[iPt],2) + TMath::Power(intErrLowMassMC[iPt]/intLowMassMC[iPt],2));
                } else {
                    ratioLowMass[iPt]           = 1;
                    ratioErrLowMass[iPt]        = 0.1;
                }
            }

            if (j < 2){
                //histonameMCTruth = Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);  MC truth > MC does not make sense
                histonameMCTruth = Form("Mapping_TrueMeson_InvMassUnweighted_in_Pt_Bin%02d", iPt);
                histoTrueMCInvMassPtBin[iPt]=(TH1D*)fileRawSignalMC->Get(histonameMCTruth);
                if((histoSignalMCInvMassPtBin[iPt]->GetNbinsX()) == (histoTrueMCInvMassPtBin[iPt]->GetNbinsX())){
                    for(Int_t m=0; m<histoSignalMCInvMassPtBin[iPt]->GetNbinsX(); m++){
                        histoMCRemBGInvMassPtBin[iPt]->SetBinContent(m,histoSignalMCInvMassPtBin[iPt]->GetBinContent(m)-histoTrueMCInvMassPtBin[iPt]->GetBinContent(m));
                    }
                }

                if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "integrals: " << integralData << "\t" << integralMC << endl;
                if (integralData != 0) histoSignalDataInvMassPtBin[iPt]->Scale(1./integralData);
                if (integralMC != 0){
                    histoSignalMCInvMassPtBin[iPt]->Scale(1./integralMC);
                    if(j==0){ // why not plot MC truth for histo with left side normalization?
                        if(histoTrueMCInvMassPtBin[iPt])histoTrueMCInvMassPtBin[iPt]->Scale(1./integralMC);
                        if(histoMCRemBGInvMassPtBin[iPt])histoMCRemBGInvMassPtBin[iPt]->Scale(1./integralMC);
                    }
                }
            }
            if (j == 0){
                histoLinResBGData->SetBinContent(histoLinResBGData->FindBin(ptValue),histoLinResBGData->GetBinContent(histoLinResBGData->FindBin(ptValue))/integralData) ;
                histoLinResBGData->SetBinError(histoLinResBGData->FindBin(ptValue),histoLinResBGData->GetBinError(histoLinResBGData->FindBin(ptValue))/integralData) ;
                histoConstResBGData->SetBinContent(histoConstResBGData->FindBin(ptValue),histoConstResBGData->GetBinContent(histoConstResBGData->FindBin(ptValue))/integralData) ;
                histoConstResBGData->SetBinError(histoConstResBGData->FindBin(ptValue),histoConstResBGData->GetBinError(histoConstResBGData->FindBin(ptValue))/integralData) ;
                histoLinResBGMC->SetBinContent(histoLinResBGMC->FindBin(ptValue),histoLinResBGMC->GetBinContent(histoLinResBGMC->FindBin(ptValue))/integralMC) ;
                histoLinResBGMC->SetBinError(histoLinResBGMC->FindBin(ptValue),histoLinResBGMC->GetBinError(histoLinResBGMC->FindBin(ptValue))/integralMC) ;
                histoConstResBGMC->SetBinContent(histoConstResBGMC->FindBin(ptValue),histoConstResBGMC->GetBinContent(histoConstResBGMC->FindBin(ptValue))/integralMC) ;
                histoConstResBGMC->SetBinError(histoConstResBGMC->FindBin(ptValue),histoConstResBGMC->GetBinError(histoConstResBGMC->FindBin(ptValue))/integralMC) ;
                // fit remaining background
                fitLinRemBG[iPt]  = FitSubtractedInvMassLinRemBG(histoMCRemBGInvMassPtBin[iPt], fMesonFitRange, fMesonRange, chi2ndfLin[iPt], iPt, debugMesonQual);
                fitPol2RemBG[iPt] = FitSubtractedInvMassPol2RemBG(histoMCRemBGInvMassPtBin[iPt], fMesonFitRange, fMesonRange, chi2ndfPol2[iPt], iPt, debugMesonQual);
                fitPol3RemBG[iPt] = FitSubtractedInvMassPol3RemBG(histoMCRemBGInvMassPtBin[iPt], fMesonFitRange, fMesonRange, chi2ndfPol3[iPt], iPt, debugMesonQual);
            }
            if (j == 2){
                histoRatioYieldLowMassDataDivMC->SetBinContent(histoRatioYieldLowMassDataDivMC->FindBin(ptValue), ratioLowMass[iPt] );
                histoRatioYieldLowMassDataDivMC->SetBinError(histoRatioYieldLowMassDataDivMC->FindBin(ptValue), ratioErrLowMass[iPt]);
            }

        } // end of loop over pT bins

        if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "Loop over Pt Bins done" << endl;

        TCanvas * canvasLineShape = new TCanvas("CanvasLineShape","",2800,1800);  // gives the page size
        canvasLineShape->SetTopMargin(0.00);
        canvasLineShape->SetBottomMargin(0.00);
        canvasLineShape->SetRightMargin(0.00);
        canvasLineShape->SetLeftMargin(0.00);

        TPad * padLineShape = new TPad("PadLineShape","",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
        padLineShape->SetFillColor(0);
        padLineShape->GetFrame()->SetFillColor(0);
        padLineShape->SetBorderMode(0);
        padLineShape->SetLogy(0);
        padLineShape->Divide(fColumn,fRow,0.0,0.0);
        padLineShape->Draw();

        // cout<<"fColumn: "<<fColumn<<" fRow: "<<fRow<<endl;

        Double_t relWidthLogo;
        if (mesonType.CompareTo("Pi0") == 0){
            relWidthLogo=0.5;
        } else {
            relWidthLogo=0.3;
        }
        Double_t padXWidth = 1400/fColumn;
        Double_t padYWidth = 900/fRow;


        Int_t place = 0;
        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt && iPt < 200;iPt++){
            if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "Pt: "<<iPt<<" of "<<fNBinsPt<<endl;
            Double_t startPt = fBinsPt[iPt];
            Double_t endPt = fBinsPt[iPt+1];

            // cout << startPt << "\t" << endPt << endl;

            place = place + 1; //give the right place in the page
            if(place == fColumn) {

                iPt--;
                padLineShape->cd(place);

                Double_t nPixels = 13;
                Double_t textHeight = 0.07;

                Double_t startTextX = 0.10;
                Double_t startTextY = 0.8;
                Double_t differenceText = textHeight*1.25;

                PlotLabelsInvMassInPtPlots ( startTextX, startTextY, textHeight, differenceText, textAlice, date, energyText, decayChannel, DetectionChannel);

                TLegend* legendLineShape = new TLegend(startTextX,startTextY-8.5*differenceText,1,startTextY-5*differenceText);
                legendLineShape->SetTextSize(textHeight);
                legendLineShape->SetTextFont(62);
                legendLineShape->SetFillColor(0);
                legendLineShape->SetFillStyle(0);
                legendLineShape->SetLineWidth(0);
                legendLineShape->SetLineColor(0);
                legendLineShape->SetMargin(0.15);
                Size_t markersize = histoSignalDataInvMassPtBin[fStartPtBin]->GetMarkerSize();
                histoSignalDataInvMassPtBin[fStartPtBin]->SetMarkerSize(2*markersize);
                legendLineShape->AddEntry(histoSignalDataInvMassPtBin[fStartPtBin],"Data","ep");
                Size_t markersize2 = histoSignalMCInvMassPtBin[fStartPtBin]->GetMarkerSize();
                histoSignalMCInvMassPtBin[fStartPtBin]->SetMarkerSize(2*markersize2);
                legendLineShape->AddEntry(histoSignalMCInvMassPtBin[fStartPtBin],"MC reconstructed","ep");
                if (j == 0 && histoTrueMCInvMassPtBin[fStartPtBin]){
                    Size_t linesize = histoTrueMCInvMassPtBin[fStartPtBin]->GetLineWidth();
                    histoTrueMCInvMassPtBin[fStartPtBin]->SetLineWidth(linesize);
                    legendLineShape->AddEntry(histoTrueMCInvMassPtBin[fStartPtBin],"MC truth" ,"l");
                }
                if( j == 0 && histoMCRemBGInvMassPtBin[fStartPtBin]){
                    Size_t linesize = histoMCRemBGInvMassPtBin[fStartPtBin]->GetLineWidth();
                    histoMCRemBGInvMassPtBin[fStartPtBin]->SetLineWidth(linesize);
                    legendLineShape->AddEntry(histoMCRemBGInvMassPtBin[fStartPtBin],"MC rec - MC truth" ,"l");
                }
                legendLineShape->Draw();

            } else {

                padLineShape->cd(place);
                padLineShape->cd(place)->SetTopMargin(0.12);
                padLineShape->cd(place)->SetBottomMargin(0.15);
                padLineShape->cd(place)->SetRightMargin(0.05);
                padLineShape->cd(place)->SetLeftMargin(0.15);

                if(histoSignalDataInvMassPtBin[iPt] == NULL || histoSignalMCInvMassPtBin[iPt]==NULL) continue;
                Double_t maxY   = 0;
                if (histoTrueMCInvMassPtBin[iPt])
                    maxY        = histoTrueMCInvMassPtBin[iPt]->GetMaximum();
                if (maxY < histoSignalDataInvMassPtBin[iPt]->GetMaximum())
                    maxY        = histoSignalDataInvMassPtBin[iPt]->GetMaximum();
                if (maxY < histoSignalMCInvMassPtBin[iPt]->GetMaximum())
                    maxY        = histoSignalMCInvMassPtBin[iPt]->GetMaximum();
                maxY            = maxY*1.4;

                Double_t minY   = FindSmallestEntryIn1D(histoTrueMCInvMassPtBin[iPt]);
                if (minY > FindSmallestEntryIn1D(histoSignalDataInvMassPtBin[iPt]))
                    minY        = FindSmallestEntryIn1D(histoSignalDataInvMassPtBin[iPt]);
                if (minY > FindSmallestEntryIn1D(histoSignalMCInvMassPtBin[iPt]))
                    minY        = FindSmallestEntryIn1D(histoSignalMCInvMassPtBin[iPt]);

                if (j == 0){
                    // Comparison between Data, MC and validated MC
                    if (histoTrueMCInvMassPtBin[iPt]) {
                        histoTrueMCInvMassPtBin[iPt]->GetYaxis()->SetRangeUser(minY,maxY);
                        DrawGammaHistoColored( histoTrueMCInvMassPtBin[iPt],
                                Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
                                "M_{#gamma#gamma} (GeV/c^{2})", "",
                                fMesonRange[0],fMesonRange[1],1,634,-1);
                        histoTrueMCInvMassPtBin[iPt]->GetYaxis()->SetRangeUser(minY,maxY);
                        histoTrueMCInvMassPtBin[iPt]->Draw("hist");

                    }
                    if (histoSignalDataInvMassPtBin[iPt]) {
                        DrawGammaHistoColored( histoSignalDataInvMassPtBin[iPt],
                                Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
                                "M_{#gamma#gamma} (GeV/c^{2})", "",
                                fMesonRange[0],fMesonRange[1],0,1,20,0.8);
                    }
                    if (histoSignalMCInvMassPtBin[iPt]){
                        DrawGammaHistoColored( histoSignalMCInvMassPtBin[iPt],
                                Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
                                "M_{#gamma#gamma} (GeV/c^{2})", "",
                                fMesonRange[0],fMesonRange[1],0,860,1,0.8);
                    }
                    if(histoMCRemBGInvMassPtBin[iPt]){
                        DrawGammaHistoColored( histoMCRemBGInvMassPtBin[iPt],
                                               Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
                                               "M_{#gamma#gamma} (GeV/c^{2})", "",
                                               fMesonRange[0],fMesonRange[1],0,kGreen+2,-1);
                    }
                    if(fitLinRemBG[iPt]){
                        fitLinRemBG[iPt]->SetLineColor(kCyan+2);
                        fitLinRemBG[iPt]->Draw("same");
                    }
                    if(fitPol2RemBG[iPt]){
                        fitPol2RemBG[iPt]->SetLineColor(kRed+1);
                        fitPol2RemBG[iPt]->Draw("same");
                    }
                    if(fitPol3RemBG[iPt]){
                        fitPol3RemBG[iPt]->SetLineColor(kGray+1);
                        fitPol3RemBG[iPt]->Draw("same");
                    }
                } else { // Comparison between data and MC
                    if (histoSignalDataInvMassPtBin[iPt]) {
                        histoSignalDataInvMassPtBin[iPt]->GetYaxis()->SetRangeUser(minY,maxY);
                        DrawGammaHistoColored( histoSignalDataInvMassPtBin[iPt],
                                Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
                                "M_{#gamma#gamma} (GeV/c^{2})", "",
                                fMesonRange[0],fMesonRange[1],1,1,20,0.8);
                    }
                    if (histoSignalMCInvMassPtBin[iPt]){
                        DrawGammaHistoColored( histoSignalMCInvMassPtBin[iPt],
                                Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
                                "M_{#gamma#gamma} (GeV/c^{2})", "",
                                fMesonRange[0],fMesonRange[1],0,860,1,0.8);
                    }
                    if (j == 2){
                        TLatex *ratio = 		new TLatex(0.7, 0.75, Form("%2.2f",ratioLowMass[iPt]));
                        ratio->SetNDC();
                        ratio->SetTextColor(1);
                        ratio->SetTextSize(0.08);
                        ratio->Draw();
                        cout << ratioLowMass[iPt] << "\t plotted" << endl;
                        delete ratio;
                    }
                }
            }
        } // end of loop over pT bins
        if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; " << "saving" << endl;
        cout << "nameLineShapePlot: " << nameLineShapePlot.Data() << endl;
        if(j==0) {
            canvasLineShape->SaveAs(nameLineShapePlot.Data());
        } else if (j==1) {
            canvasLineShape->SaveAs(nameLineShapePlotLeft.Data());
        } else {
            canvasLineShape->SaveAs(nameLineShapePlotFull.Data());
        }

        if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; " << "deleting" << endl;
        if (padLineShape) delete padLineShape;
        if (canvasLineShape) delete canvasLineShape;
        if (debugMesonQual) cout << "Debug Line: " << __LINE__ <<"; "<< "Ending of Loop j = " << j << " fNBinsPt = " << fNBinsPt << endl;
    } // end of loop over j

    // create histograms with Chi2 of remaining background fits
    histoChi2RemLinBGFit  = new TH1D("histoChi2RemLinBGFit","", fNBinsPt, fBinsPt);
    histoChi2RemPol2BGFit = new TH1D("histoChi2RemPol2BGFit","", fNBinsPt, fBinsPt);
    histoChi2RemPol3BGFit = new TH1D("histoChi2RemPol3BGFit","", fNBinsPt, fBinsPt);
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt && iPt < 200;iPt++){
        histoChi2RemLinBGFit->SetBinContent(iPt,chi2ndfLin[iPt]);
        histoChi2RemLinBGFit->SetBinError(iPt,0);
        histoChi2RemPol2BGFit->SetBinContent(iPt,chi2ndfPol2[iPt]);
        histoChi2RemPol2BGFit->SetBinError(iPt,0);
        histoChi2RemPol3BGFit->SetBinContent(iPt,chi2ndfPol3[iPt]);
        histoChi2RemPol3BGFit->SetBinError(iPt,0);
    }
    // **************************************************************************************************************
    // ************************ Chi2/ndf compared MC vs Data ********************************************************
    // **************************************************************************************************************
    TCanvas* canvasChi2 = new TCanvas("canvasChi2","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasChi2, 0.092, 0.01, 0.02, 0.082);

        Double_t maxChi2    = histoChi2Data->GetMaximum();
        if (maxChi2 < histoChi2MC->GetMaximum())
            maxChi2         = histoChi2MC->GetMaximum();
        maxChi2             = maxChi2*1.2;

        histoChi2Data->GetYaxis()->SetRangeUser(0, maxChi2);
        DrawAutoGammaMesonHistos( histoChi2Data,
                                    "", "#it{p}_{T} (GeV/#it{c})", "#it{#chi}^{2}/ndf",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoChi2Data, 20, 2, kBlack, kBlack);
        histoChi2Data->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoChi2MC, 24, 2, kRed+1, kRed+1);
        histoChi2MC->DrawCopy("same,e1,p");

        TLegend* legendChi2 = GetAndSetLegend2(0.85, 0.13, 0.95, 0.13+(0.035*2), 0.035, 1, "", 42, 0.25);
        legendChi2->AddEntry(histoChi2Data,"Data");
        legendChi2->AddEntry(histoChi2MC,"MC");
        legendChi2->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());

    canvasChi2->Update();
    canvasChi2->SaveAs(Form("%s/%s_Chi2_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    canvasChi2->cd();
        maxChi2    = histoChi2_Pol2_Data->GetMaximum();
        if (maxChi2 < histoChi2_Pol2_MC->GetMaximum())
            maxChi2         = histoChi2_Pol2_MC->GetMaximum();
        maxChi2             = maxChi2*1.2;

        histoChi2_Pol2_Data->GetYaxis()->SetRangeUser(0, maxChi2);
        DrawAutoGammaMesonHistos( histoChi2_Pol2_Data,
                                    "", "#it{p}_{T} (GeV/#it{c})", "#it{#chi}^{2}/ndf",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoChi2_Pol2_Data, 20, 2, kBlack, kBlack);
        histoChi2_Pol2_Data->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoChi2_Pol2_MC, 24, 2, kRed+1, kRed+1);
        histoChi2_Pol2_MC->DrawCopy("same,e1,p");

        legendChi2->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());

    canvasChi2->Update();
    canvasChi2->SaveAs(Form("%s/%s_Chi2_WithPol2BG_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    canvasChi2->cd();
        maxChi2    = histoChi2_Exp1_Data->GetMaximum();
        if (maxChi2 < histoChi2_Exp1_MC->GetMaximum())
            maxChi2         = histoChi2_Exp1_MC->GetMaximum();
        maxChi2             = maxChi2*1.2;

        histoChi2_Exp1_Data->GetYaxis()->SetRangeUser(0, maxChi2);
        DrawAutoGammaMesonHistos( histoChi2_Exp1_Data,
                                    "", "#it{p}_{T} (GeV/#it{c})", "#it{#chi}^{2}/ndf",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoChi2_Exp1_Data, 20, 2, kBlack, kBlack);
        histoChi2_Exp1_Data->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoChi2_Exp1_MC, 24, 2, kRed+1, kRed+1);
        histoChi2_Exp1_MC->DrawCopy("same,e1,p");

        legendChi2->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());

    canvasChi2->Update();
    canvasChi2->SaveAs(Form("%s/%s_Chi2_WithExp1BG_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    canvasChi2->cd();
        maxChi2    = histoChi2_Exp2_Data->GetMaximum();
        if (maxChi2 < histoChi2_Exp2_MC->GetMaximum())
            maxChi2         = histoChi2_Exp2_MC->GetMaximum();
        maxChi2             = maxChi2*1.2;

        histoChi2_Exp2_Data->GetYaxis()->SetRangeUser(0, maxChi2);
        DrawAutoGammaMesonHistos( histoChi2_Exp2_Data,
                                    "", "#it{p}_{T} (GeV/#it{c})", "#it{#chi}^{2}/ndf",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoChi2_Exp2_Data, 20, 2, kBlack, kBlack);
        histoChi2_Exp2_Data->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoChi2_Exp2_MC, 24, 2, kRed+1, kRed+1);
        histoChi2_Exp2_MC->DrawCopy("same,e1,p");

        legendChi2->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());

    canvasChi2->Update();
    canvasChi2->SaveAs(Form("%s/%s_Chi2_WithExp2BG_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    // Plot remaining BG fit Chi2
    canvasChi2->cd();
        maxChi2    = histoChi2RemLinBGFit->GetMaximum();
        if (maxChi2 < histoChi2RemPol2BGFit->GetMaximum())
            maxChi2         = histoChi2RemPol2BGFit->GetMaximum();
        if (maxChi2 < histoChi2RemPol3BGFit->GetMaximum())
            maxChi2         = histoChi2RemPol3BGFit->GetMaximum();
        maxChi2             = maxChi2*1.2;

        histoChi2RemLinBGFit->GetYaxis()->SetRangeUser(0, maxChi2);
        DrawAutoGammaMesonHistos( histoChi2RemLinBGFit,
                                    "", "#it{p}_{T} (GeV/#it{c})", "#it{#chi}^{2}/ndf",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoChi2RemLinBGFit, 20, 2, kCyan+2, kCyan+2);
        histoChi2RemLinBGFit->DrawCopy("same,p,L");
        DrawGammaSetMarker(histoChi2RemPol2BGFit, 34, 2, kRed+1, kRed+1);
        histoChi2RemPol2BGFit->DrawCopy("same,p,L");
        DrawGammaSetMarker(histoChi2RemPol3BGFit, 24, 2, kGray+2, kGray+2);
        histoChi2RemPol3BGFit->DrawCopy("same,p,L");

        PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 0.035, energyText.Data(), decayChannel.Data(), Form("Int. %.2f-%.2f",fMesonFitRange[0],fMesonFitRange[1]));

        TLegend* legendChi2DiffFits = GetAndSetLegend2(0.85, 0.8, 0.95, 0.8+(0.035*3), 0.035, 1, "", 42, 0.25);  // xleft, ydown, xright, yup
        legendChi2DiffFits->AddEntry(histoChi2RemLinBGFit,"Linear");
        legendChi2DiffFits->AddEntry(histoChi2RemPol2BGFit,"Pol2");
        legendChi2DiffFits->AddEntry(histoChi2RemPol3BGFit,"Pol3");
        legendChi2DiffFits->Draw();

    canvasChi2->Update();
    canvasChi2->SaveAs(Form("%s/%s_Chi2_LinAndOtherRemBG_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    if (legendChi2) delete legendChi2;
    if (canvasChi2) delete canvasChi2;

    // **************************************************************************************************************
    // ************************ Res BG yield/ tot BG yield compared MC vs Data **************************************
    // **************************************************************************************************************
    TCanvas* canvasResBGYieldDivTotBG = new TCanvas("canvasResBGYieldDivTotBG","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasResBGYieldDivTotBG, 0.092, 0.01, 0.02, 0.082);

        Double_t maxResBGYieldDivTotBG    = histoResBGYieldVsTotBGData->GetMaximum();
        if (maxResBGYieldDivTotBG < histoResBGYieldVsTotBGMC->GetMaximum())
            maxResBGYieldDivTotBG         = histoResBGYieldVsTotBGMC->GetMaximum();
        maxResBGYieldDivTotBG             = maxResBGYieldDivTotBG*1.2;

        Double_t minResBGYieldDivTotBG    = histoResBGYieldVsTotBGData->GetMinimum();
        if (minResBGYieldDivTotBG > histoResBGYieldVsTotBGMC->GetMinimum())
            minResBGYieldDivTotBG         = histoResBGYieldVsTotBGMC->GetMinimum();
        if (minResBGYieldDivTotBG < 0)
            minResBGYieldDivTotBG         = minResBGYieldDivTotBG*1.4;
        else
            minResBGYieldDivTotBG         = minResBGYieldDivTotBG*0.6;

        histoResBGYieldVsTotBGData->GetYaxis()->SetRangeUser(minResBGYieldDivTotBG, maxResBGYieldDivTotBG);
        DrawAutoGammaMesonHistos( histoResBGYieldVsTotBGData,
                                    "", "#it{p}_{T} (GeV/#it{c})", "Res BG/ Tot BG",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoResBGYieldVsTotBGData, 20, 2, kBlack, kBlack);
        histoResBGYieldVsTotBGData->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoResBGYieldVsTotBGMC, 24, 2, kRed+1, kRed+1);
        histoResBGYieldVsTotBGMC->DrawCopy("same,e1,p");

        TLegend* legendResBGYieldDivTotBG = GetAndSetLegend2(0.85, 0.13, 0.95, 0.13+(0.035*2), 0.035, 1, "", 42, 0.25);
        legendResBGYieldDivTotBG->AddEntry(histoResBGYieldVsTotBGData,"Data");
        legendResBGYieldDivTotBG->AddEntry(histoResBGYieldVsTotBGMC,"MC");
        legendResBGYieldDivTotBG->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());

    canvasResBGYieldDivTotBG->Update();
    canvasResBGYieldDivTotBG->SaveAs(Form("%s/%s_ResBGYieldDivTotBG_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
    if (legendResBGYieldDivTotBG) delete legendResBGYieldDivTotBG;
    if (canvasResBGYieldDivTotBG) delete canvasResBGYieldDivTotBG;

    // **************************************************************************************************************
    // ************************ Res BG yield/ Res BG + Signal yield compared MC vs Data **************************************
    // **************************************************************************************************************
    TCanvas* canvasResBGYieldDivResBGPlSig = new TCanvas("canvasResBGYieldDivResBGPlSig","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasResBGYieldDivResBGPlSig, 0.092, 0.01, 0.02, 0.082);

        Double_t maxResBGYieldDivResBGPlSig    = histoResBGYieldVsResBGPlSigData->GetMaximum();
        if (maxResBGYieldDivResBGPlSig < histoResBGYieldVsResBGPlSigMC->GetMaximum())
            maxResBGYieldDivResBGPlSig         = histoResBGYieldVsResBGPlSigMC->GetMaximum();
        maxResBGYieldDivResBGPlSig             = maxResBGYieldDivResBGPlSig*1.2;

        Double_t minResBGYieldDivResBGPlSig    = histoResBGYieldVsResBGPlSigData->GetMinimum();
        if (minResBGYieldDivResBGPlSig > histoResBGYieldVsResBGPlSigMC->GetMinimum())
            minResBGYieldDivResBGPlSig         = histoResBGYieldVsResBGPlSigMC->GetMinimum();
        if (minResBGYieldDivResBGPlSig < 0)
            minResBGYieldDivResBGPlSig         = minResBGYieldDivResBGPlSig*1.4;
        else
            minResBGYieldDivResBGPlSig         = minResBGYieldDivResBGPlSig*0.6;


        histoResBGYieldVsResBGPlSigData->GetYaxis()->SetRangeUser(minResBGYieldDivResBGPlSig, maxResBGYieldDivResBGPlSig);
        DrawAutoGammaMesonHistos( histoResBGYieldVsResBGPlSigData,
                                    "", "#it{p}_{T} (GeV/#it{c})", "Res BG/ (Res BG + Signal)",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoResBGYieldVsResBGPlSigData, 20, 2, kBlack, kBlack);
        histoResBGYieldVsResBGPlSigData->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoResBGYieldVsResBGPlSigMC, 24, 2, kRed+1, kRed+1);
        histoResBGYieldVsResBGPlSigMC->DrawCopy("same,e1,p");

        TLegend* legendResBGYieldDivResBGPlSig = GetAndSetLegend2(0.85, 0.13, 0.95, 0.13+(0.035*2), 0.035, 1, "", 42, 0.25);
        legendResBGYieldDivResBGPlSig->AddEntry(histoResBGYieldVsResBGPlSigData,"Data");
        legendResBGYieldDivResBGPlSig->AddEntry(histoResBGYieldVsResBGPlSigMC,"MC");
        legendResBGYieldDivResBGPlSig->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());

    canvasResBGYieldDivResBGPlSig->Update();
    canvasResBGYieldDivResBGPlSig->SaveAs(Form("%s/%s_ResBGYieldDivResBGPlSig_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
    if (legendResBGYieldDivResBGPlSig) delete legendResBGYieldDivResBGPlSig;
    if (canvasResBGYieldDivResBGPlSig) delete canvasResBGYieldDivResBGPlSig;

    // **************************************************************************************************************
    // ************************ Res BG slope compared MC vs Data ****************************************************
    // **************************************************************************************************************
    TCanvas* canvasResBGSlope = new TCanvas("canvasResBGSlope","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasResBGSlope, 0.092, 0.01, 0.035, 0.082);

        Double_t maxResBGSlope    = histoLinResBGData->GetMaximum();
        if (maxResBGSlope < histoLinResBGMC->GetMaximum())
            maxResBGSlope         = histoLinResBGMC->GetMaximum();
        maxResBGSlope             = maxResBGSlope*1.2;

        Double_t minResBGSlope    = histoLinResBGData->GetMinimum();
        if (minResBGSlope > histoLinResBGMC->GetMinimum())
            minResBGSlope         = histoLinResBGMC->GetMinimum();
        if (minResBGSlope < 0)
            minResBGSlope         = minResBGSlope*1.4;
        else
            minResBGSlope         = minResBGSlope*0.6;


        histoLinResBGData->GetYaxis()->SetRangeUser(minResBGSlope, maxResBGSlope);
        DrawAutoGammaMesonHistos( histoLinResBGData,
                                    "", "#it{p}_{T} (GeV/#it{c})", "Res BG slope #it{b}",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoLinResBGData, 20, 2, kBlack, kBlack);
        histoLinResBGData->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoLinResBGMC, 24, 2, kRed+1, kRed+1);
        histoLinResBGMC->DrawCopy("same,e1,p");

        TLegend* legendResBGSlope = GetAndSetLegend2(0.85, 0.13, 0.95, 0.13+(0.035*2), 0.035, 1, "", 42, 0.25);
        legendResBGSlope->AddEntry(histoLinResBGData,"Data");
        legendResBGSlope->AddEntry(histoLinResBGMC,"MC");
        legendResBGSlope->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.95, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());

    canvasResBGSlope->Update();
    canvasResBGSlope->SaveAs(Form("%s/%s_ResBGSlope_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
    if (legendResBGSlope) delete legendResBGSlope;
    if (canvasResBGSlope) delete canvasResBGSlope;

    // **************************************************************************************************************
    // ************************ Res BG const compared MC vs Data ****************************************************
    // **************************************************************************************************************
    TCanvas* canvasResBGConst = new TCanvas("canvasResBGConst","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasResBGConst, 0.092, 0.01, 0.035, 0.082);

        Double_t maxResBGConst    = histoConstResBGData->GetMaximum();
        if (maxResBGConst < histoConstResBGMC->GetMaximum())
            maxResBGConst         = histoConstResBGMC->GetMaximum();
        maxResBGConst             = maxResBGConst*1.2;

        Double_t minResBGConst    = histoConstResBGData->GetMinimum();
        if (minResBGConst > histoConstResBGMC->GetMinimum())
            minResBGConst         = histoConstResBGMC->GetMinimum();
        if (minResBGConst < 0)
            minResBGConst         = minResBGConst*1.4;
        else
            minResBGConst         = minResBGConst*0.6;

        histoConstResBGData->GetYaxis()->SetRangeUser(minResBGConst, maxResBGConst);
        DrawAutoGammaMesonHistos( histoConstResBGData,
                                    "", "#it{p}_{T} (GeV/#it{c})", "Res BG const #it{a}",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoConstResBGData, 20, 2, kBlack, kBlack);
        histoConstResBGData->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoConstResBGMC, 24, 2, kRed+1, kRed+1);
        histoConstResBGMC->DrawCopy("same,e1,p");

        TLegend* legendResBGConst = GetAndSetLegend2(0.85, 0.13, 0.95, 0.13+(0.035*2), 0.035, 1, "", 42, 0.25);
        legendResBGConst->AddEntry(histoConstResBGData,"Data");
        legendResBGConst->AddEntry(histoConstResBGMC,"MC");
        legendResBGConst->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.95, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());

    canvasResBGConst->Update();
    canvasResBGConst->SaveAs(Form("%s/%s_ResBGConst_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
    if (legendResBGConst) delete legendResBGConst;
    if (canvasResBGConst) delete canvasResBGConst;


    // **************************************************************************************************************
    // ************************ Lambda tail compared MC vs Data ****************************************************
    // **************************************************************************************************************
    TCanvas* canvasLambda = new TCanvas("canvasLambda","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasLambda, 0.092, 0.01, 0.035, 0.082);

        Double_t maxLambda    = histoLambdaTailData->GetMaximum();
        if (maxLambda < histoLambdaTailMC->GetMaximum())
            maxLambda         = histoLambdaTailMC->GetMaximum();
        maxLambda             = maxLambda*1.2;

        Double_t minLambda    = histoLambdaTailData->GetMinimum();
        if (minLambda > histoLambdaTailMC->GetMinimum())
            minLambda         = histoLambdaTailMC->GetMinimum();
        if (minLambda < 0)
            minLambda         = minLambda*1.4;
        else
            minLambda         = minLambda*0.6;

        histoLambdaTailData->GetYaxis()->SetRangeUser(minLambda, maxLambda);
        DrawAutoGammaMesonHistos( histoLambdaTailData,
                                    "", "#it{p}_{T} (GeV/#it{c})", "#it{#lambda}",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoLambdaTailData, 20, 2, kBlack, kBlack);
        histoLambdaTailData->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoLambdaTailMC, 24, 2, kRed+1, kRed+1);
        histoLambdaTailMC->DrawCopy("same,e1,p");

        TLegend* legendLambda = GetAndSetLegend2(0.85, 0.13, 0.95, 0.13+(0.035*2), 0.035, 1, "", 42, 0.25);
        legendLambda->AddEntry(histoLambdaTailData,"Data");
        legendLambda->AddEntry(histoLambdaTailMC,"MC");
        legendLambda->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.95, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());

    canvasLambda->Update();
    canvasLambda->SaveAs(Form("%s/%s_LambdaTail_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
    if (legendLambda) delete legendLambda;
    if (canvasLambda) delete canvasLambda;

    // **************************************************************************************************************
    // ************************ Comparison low mass yield ***********************************************************
    // **************************************************************************************************************
    TCanvas* canvasLowMassYield = new TCanvas("canvasLowMassYield","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasLowMassYield, 0.092, 0.01, 0.035, 0.082);

        Double_t maxLowMassYield    = histoRatioYieldLowMassDataDivMC->GetMaximum();
        maxLowMassYield             = maxLowMassYield*1.2;

        Double_t minLowMassYield    = histoRatioYieldLowMassDataDivMC->GetMinimum();
        minLowMassYield             = minLowMassYield*0.6;

        histoRatioYieldLowMassDataDivMC->GetYaxis()->SetRangeUser(minLowMassYield, maxLowMassYield);
        DrawAutoGammaMesonHistos(   histoRatioYieldLowMassDataDivMC,
                                    "", "#it{p}_{T} (GeV/#it{c})", "access yield #it{M}_{#gamma#gamma} < 0.05 (GeV/#it{c}^{2})",
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7,
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoRatioYieldLowMassDataDivMC, 20, 2, kBlack, kBlack);
        histoRatioYieldLowMassDataDivMC->DrawCopy("same,e1,p");

        PutProcessLabelAndEnergyOnPlot(0.15, 0.95, 0.035, energyText.Data(), decayChannel.Data(), DetectionChannel.Data());
        DrawGammaLines( histoRatioYieldLowMassDataDivMC->GetBinCenter(1) - (histoRatioYieldLowMassDataDivMC->GetBinWidth(1)/2.),
                        histoRatioYieldLowMassDataDivMC->GetBinCenter(histoRatioYieldLowMassDataDivMC->GetNbinsX()) +
                        (histoRatioYieldLowMassDataDivMC->GetBinWidth(histoRatioYieldLowMassDataDivMC->GetNbinsX())/2.), 1.0, 1.0,2.0, kGray+2 ,7);

    canvasLowMassYield->Update();
    canvasLowMassYield->SaveAs(Form("%s/%s_LowMassYield_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    if (canvasLowMassYield) delete canvasLowMassYield;

    if (histoRatioYieldLowMassDataDivMC){
        // write separate output file with ratio off excess yield at low masses
        TFile* outputFile           = new TFile(Form("%s/%s/ExcessYieldAtLowInvMass.root",fCutSelection.Data(),energyFlag.Data()),"UPDATE");
        if (histoRatioYieldLowMassDataDivMC)histoRatioYieldLowMassDataDivMC->Write(Form("%s_ExcessYieldAtLowMass",mesonType.Data()), TObject::kOverwrite);
        outputFile->Close();
        delete outputFile;
    }

    Delete();
}



