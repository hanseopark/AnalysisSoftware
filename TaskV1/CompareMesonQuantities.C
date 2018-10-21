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

//**********************************************************************************
//******************* return minimum for 1D histo  *********************************
//**********************************************************************************
Double_t FindSmallestEntryIn1D(TH1* histo){
    Double_t minimum = 1;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        if (histo->GetBinContent(i) < minimum ){
            minimum = histo->GetBinContent(i);
        }
    }
    return minimum;
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

    TFile fileRawSignalData(dataFilename);
    TFile fileRawSignalMC(mcFilename);


    cout << dataFilename << endl;
    cout << mcFilename << endl;
    cout << fCutSelection.Data() << endl;
    cout << mesonType.Data() << endl;
    cout << fSuffix.Data()<< endl;
    cout << energyFlag.Data() << endl;
    cout << numberOfBins << endl;

    TH1D* histoChi2Data                     = (TH1D*) fileRawSignalData.Get("histoChi2_0");
    TH1D* histoChi2_Pol2_Data               = (TH1D*) fileRawSignalData.Get("histoChi2_1");
    TH1D* histoChi2_Exp1_Data               = (TH1D*) fileRawSignalData.Get("histoChi2_2");
    TH1D* histoChi2_Exp2_Data               = (TH1D*) fileRawSignalData.Get("histoChi2_3");
    TH1D* histoConstResBGData               = (TH1D*) fileRawSignalData.Get("histoResidualBGcon");
    TH1D* histoLinResBGData                 = (TH1D*) fileRawSignalData.Get("histoResidualBGlin");
    TH1D* histoResBGYieldVsTotBGData        = (TH1D*) fileRawSignalData.Get("histoRatioResBGYield");
    TH1D* histoResBGYieldVsResBGPlSigData   = (TH1D*) fileRawSignalData.Get("histoRatioResBGYieldToSPlusResBG");
    TH1D* histoLambdaTailData               = (TH1D*) fileRawSignalData.Get("histoLambdaTail");
    TH1D* histoChi2MC                       = (TH1D*) fileRawSignalMC.Get("histoChi2_0");
    TH1D* histoChi2_Pol2_MC                 = (TH1D*) fileRawSignalMC.Get("histoChi2_1");
    TH1D* histoChi2_Exp1_MC                 = (TH1D*) fileRawSignalMC.Get("histoChi2_2");
    TH1D* histoChi2_Exp2_MC                 = (TH1D*) fileRawSignalMC.Get("histoChi2_3");
    TH1D* histoConstResBGMC                 = (TH1D*) fileRawSignalMC.Get("histoResidualBGcon");
    TH1D* histoLinResBGMC                   = (TH1D*) fileRawSignalMC.Get("histoResidualBGlin");
    TH1D* histoResBGYieldVsTotBGMC          = (TH1D*) fileRawSignalMC.Get("histoRatioResBGYield");
    TH1D* histoResBGYieldVsResBGPlSigMC     = (TH1D*) fileRawSignalMC.Get("histoRatioResBGYieldToSPlusResBG");
    TH1D* histoLambdaTailMC                 = (TH1D*) fileRawSignalMC.Get("histoLambdaTail");
    TH1D* histoRatioYieldLowMassDataDivMC   = (TH1D*) histoChi2Data->Clone("histoRatioYieldLowMassDataDivMC");

    Double_t* fMesonRange                   = NULL;
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
    InitializeBinning(mesonType, numberOfBins, energyFlag, directPhoton, mode, fEventCutSelection, fClusterCutSelection, -1, kFALSE, "", "");

    Double_t peakRange[2]   = {0.10,0.15};
    if (mesonType.CompareTo("Pi0") == 0 || mesonType.CompareTo("Pi0EtaBinning") == 0){
        fMesonRange         = new Double_t[2];
        fMesonRange[0]      = 0.;
        fMesonRange[1]      = 0.3;
    } else if (mesonType.CompareTo("Eta") == 0){
        fMesonRange         = new Double_t[2];
        fMesonRange[0]      = 0.35;
        fMesonRange[1]      = 0.79;
        peakRange[0]        = 0.50;
        peakRange[0]        = 0.60;
    } else if (mesonType.CompareTo("EtaPrime") == 0){
        fMesonRange         = new Double_t[2];
        fMesonRange[0]      = 0.9;
        fMesonRange[1]      = 1.0;
        peakRange[0]        = 0.90;
        peakRange[0]        = 1.00;
    }

    cout << "Start bin for " << mesonType << " (mode " << mode << "): " << fStartPtBin << endl;

    //******************************* Reading histograms **************************************************************
    TH1D *  histoSignalDataInvMassPtBin[100];
    TH1D *  histoSignalMCInvMassPtBin[100];
    TH1D *  histoTrueMCInvMassPtBin[100];
    Double_t intLowMassData[100]                    = {0.};
    Double_t intLowMassMC[100]                      = {0.};
    Double_t intErrLowMassData[100]                 = {0.};
    Double_t intErrLowMassMC[100]                   = {0.};
    Double_t ratioLowMass[100]                      = {1.};
    Double_t ratioErrLowMass[100]                   = {1.};
    for(Int_t j=0;j<3;j++){
        TCanvas * canvasDummy = new TCanvas("canvasDummy","",2800,1800);  // gives the page size
        canvasDummy->SetTopMargin(0.02);
        canvasDummy->SetBottomMargin(0.02);
        canvasDummy->SetRightMargin(0.02);
        canvasDummy->SetLeftMargin(0.02);
        TString histonameSignal;
        TString histonameMCTruth;
        for(Int_t iPt=fStartPtBin; iPt<fNBinsPt; iPt++){
            Double_t startPt    = fBinsPt[iPt];
            Double_t endPt      = fBinsPt[iPt+1];
            Double_t ptValue    = (fBinsPt[iPt] + fBinsPt[iPt+1])/2.;
            if(j==0){
                histonameSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", iPt);
            } else if(j==1) {
                histonameSignal = Form("fHistoMappingSignalInvMassLeft_in_Pt_Bin%02d", iPt);
            } else if(j==2) {
                histonameSignal = Form("Mapping_GG_InvMass_in_Pt_Bin%02d", iPt);
            }
            histoSignalDataInvMassPtBin[iPt]    = (TH1D*)fileRawSignalData.Get(histonameSignal);
            histoSignalMCInvMassPtBin[iPt]      = (TH1D*)fileRawSignalMC.Get(histonameSignal);

            if (j == 2){
                histoSignalDataInvMassPtBin[iPt]->Rebin(4);
                histoSignalMCInvMassPtBin[iPt]->Rebin(4);
            }
            Double_t integralData = histoSignalDataInvMassPtBin[iPt]->Integral(histoSignalDataInvMassPtBin[iPt]->FindBin(fMesonRange[0]+0.0001),histoSignalDataInvMassPtBin[iPt]->FindBin(fMesonRange[1]-0.0001) );
            Double_t integralMC   = histoSignalMCInvMassPtBin[iPt]->Integral(histoSignalMCInvMassPtBin[iPt]->FindBin(fMesonRange[0]+0.0001),histoSignalMCInvMassPtBin[iPt]->FindBin(fMesonRange[1]-0.0001) );
            if (j == 2){
                integralData = histoSignalDataInvMassPtBin[iPt]->Integral(histoSignalDataInvMassPtBin[iPt]->FindBin(peakRange[0]+0.0001),histoSignalDataInvMassPtBin[iPt]->FindBin(peakRange[1]-0.0001) );
                integralMC   = histoSignalMCInvMassPtBin[iPt]->Integral(histoSignalMCInvMassPtBin[iPt]->FindBin(peakRange[0]+0.0001),histoSignalMCInvMassPtBin[iPt]->FindBin(peakRange[1]-0.0001) );
            }

            if (integralData < 0 || integralMC < 0 || j == 1){
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
                    ratioErrLowMass[iPt]        = ratioLowMass[iPt]* TMath::Sqrt(pow(intErrLowMassData[iPt]/intLowMassData[iPt],2) + pow(intErrLowMassMC[iPt]/intLowMassMC[iPt],2));
                } else {
                    ratioLowMass[iPt]           = 1;
                    ratioErrLowMass[iPt]        = 0.1;
                }
            }

            histoSignalDataInvMassPtBin[iPt]->Scale(1./integralData);
            histoSignalMCInvMassPtBin[iPt]->Scale(1./integralMC);

            if (j < 2){
                histonameMCTruth = Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);
                histoTrueMCInvMassPtBin[iPt]=(TH1D*)fileRawSignalMC.Get(histonameMCTruth);
                if(j==0)histoTrueMCInvMassPtBin[iPt]->Scale(1./integralMC);
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
            }
            if (j == 2){
                histoRatioYieldLowMassDataDivMC->SetBinContent(histoRatioYieldLowMassDataDivMC->FindBin(ptValue), ratioLowMass[iPt] );
                histoRatioYieldLowMassDataDivMC->SetBinError(histoRatioYieldLowMassDataDivMC->FindBin(ptValue), ratioErrLowMass[iPt]);
            }

        }

        delete canvasDummy;

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
        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
            // cout<<"Pt: "<<iPt<<" of "<<fNBinsPt<<endl;
            Double_t startPt = fBinsPt[iPt];
            Double_t endPt = fBinsPt[iPt+1];

            // cout << startPt << "\t" << endPt << endl;

            place = place + 1; //give the right place in the page
            if(place == fColumn) {

                iPt--;
                padLineShape->cd(place);

                Double_t nPixels = 13;
                Double_t textHeight = 0.08;

                Double_t startTextX = 0.10;
                Double_t startTextY = 0.8;
                Double_t differenceText = textHeight*1.25;

                PlotLabelsInvMassInPtPlots ( startTextX, startTextY, textHeight, differenceText, textAlice, date, energyText, decayChannel, DetectionChannel);

                TLegend* legendLineShape = new TLegend(startTextX,startTextY-4.75*differenceText,1,startTextY-(4.75+2.)*differenceText);
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
                if (j == 0){
                    Size_t linesize = histoTrueMCInvMassPtBin[fStartPtBin]->GetLineWidth();
                    histoTrueMCInvMassPtBin[fStartPtBin]->SetLineWidth(linesize);
                    legendLineShape->AddEntry(histoTrueMCInvMassPtBin[fStartPtBin],"MC truth" ,"l");
                }
                legendLineShape->Draw();

            } else {

                padLineShape->cd(place);
                padLineShape->cd(place)->SetTopMargin(0.12);
                padLineShape->cd(place)->SetBottomMargin(0.15);
                padLineShape->cd(place)->SetRightMargin(0.05);
                padLineShape->cd(place)->SetLeftMargin(0.15);

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
                } else {
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
                    }
                }
            }
        }
        // cout << "saving" << endl;
        cout << nameLineShapePlot.Data() << endl;
        if(j==0) {
            canvasLineShape->SaveAs(nameLineShapePlot.Data());
        } else if (j==1) {
            canvasLineShape->SaveAs(nameLineShapePlotLeft.Data());
        } else {
            canvasLineShape->SaveAs(nameLineShapePlotFull.Data());
        }

        // cout << "deleting" << endl;
        delete padLineShape;
        delete canvasLineShape;
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


    // **************************************************************************************************************
    // ************************ Res BG const compared MC vs Data ****************************************************
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

    // write separate output file with ratio off excess yield at low masses
    TFile* outputFile           = new TFile(Form("%s/%s/ExcessYieldAtLowInvMass.root",fCutSelection.Data(),energyFlag.Data()),"UPDATE");
        histoRatioYieldLowMassDataDivMC->Write(Form("%s_ExcessYieldAtLowMass",mesonType.Data()), TObject::kOverwrite);

    outputFile->Close();
    delete outputFile;

}



